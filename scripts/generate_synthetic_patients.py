#!/usr/bin/env python3
"""
HaploStats — Synthetic Patient Generator
Phase 6A: Generate 1,000 synthetic patients for concordance benchmarking.

Method
------
  1. Query all haplotypes + Global_freq from haplostats.db.
  2. For each patient:
     a. Sample 2 haplotypes (weighted by frequency, HWE).
     b. The chosen pair is the "Ground Truth" phased typing.
     c. Scramble into an unphased genotype dict.
     d. Randomly mask 1–3 loci (set to null) to simulate missing data.
  3. Export as CSV with Patient_ID, Observed_Unphased_Genotype (JSON),
     and True_Phased_Haplotypes.
"""

import csv
import json
import math
import random
import sys
import os
from pathlib import Path
from collections import Counter

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH  = BASE_DIR / "db" / "haplostats.db"
OUTPUT   = BASE_DIR / "data" / "clean" / "synthetic_truth_set.csv"

# Seed for reproducibility
RANDOM_SEED = 42

# 9 HLA loci
LOCI = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
        'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

LOCUS_LABELS = {
    'hla_a': 'HLA-A', 'hla_c': 'HLA-C', 'hla_b': 'HLA-B',
    'hla_drb345': 'HLA-DRB345', 'hla_drb1': 'HLA-DRB1',
    'hla_dqa1': 'HLA-DQA1', 'hla_dqb1': 'HLA-DQB1',
    'hla_dpa1': 'HLA-DPA1', 'hla_dpb1': 'HLA-DPB1',
}

N_PATIENTS = 1000
N_MASK_RANGE = (1, 3)  # randomly mask 1 to 3 loci per patient


def connect_db(db_path):
    """Load haplotypes and their global frequencies from SQLite."""
    import sqlite3
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cols = ', '.join(LOCI + ['Global_freq'])
    rows = cur.execute(
        f"SELECT {cols} FROM haplotypes "
        f"WHERE Global_freq IS NOT NULL"
    ).fetchall()
    conn.close()
    haplotypes = [dict(r) for r in rows]
    freqs = [h['Global_freq'] for h in haplotypes]

    # Normalise frequencies to sum to 1
    total = sum(freqs)
    if total > 0:
        freqs = [f / total for f in freqs]

    print(f"  [Generator] Loaded {len(haplotypes)} haplotypes, "
          f"total frequency = {total:.6f}")

    return haplotypes, freqs


def sample_patient(haplotypes, freqs, rng):
    """
    Generate one synthetic patient.
    Samples 2 haplotypes weighted by frequency (HWE: 2pq for heterozygotes).
    Returns (hap1, hap2, unphased_genotype, mask_summary).
    """
    # Sample two indices (with replacement allowed = homozygous possible)
    # Weighted random choice
    idx1 = rng.choices(range(len(haplotypes)), weights=freqs, k=1)[0]
    idx2 = rng.choices(range(len(haplotypes)), weights=freqs, k=1)[0]

    h1 = haplotypes[idx1]
    h2 = haplotypes[idx2]

    # Build true phased pair
    true_pair = {
        'haplotype_1': {loc: h1[loc] for loc in LOCI},
        'haplotype_2': {loc: h2[loc] for loc in LOCI},
        'h1_global_freq': h1['Global_freq'],
        'h2_global_freq': h2['Global_freq'],
        'h1_idx': idx1,
        'h2_idx': idx2,
        'is_homozygous': (idx1 == idx2),
    }

    # Build unphased genotype (scrambled: randomly assign each haplotype's
    # alleles to "allele1" or "allele2" per locus)
    unphased = {}
    for loc in LOCI:
        # Randomly decide order of H1/H2 for this locus
        if rng.random() < 0.5:
            a1, a2 = h1[loc], h2[loc]
        else:
            a1, a2 = h2[loc], h1[loc]

        # If both alleles are identical, store only one
        alleles = [a1, a2]
        # Deduplicate while preserving order
        seen = set()
        unique = []
        for a in alleles:
            if a not in seen:
                seen.add(a)
                unique.append(a)

        unphased[loc] = unique if len(unique) > 0 else None

    # Randomly mask 1 to 3 loci
    n_mask = rng.randint(N_MASK_RANGE[0], N_MASK_RANGE[1])
    loci_to_mask = rng.sample(LOCI, n_mask)

    mask_summary = []
    for loc in loci_to_mask:
        old_val = unphased[loc]
        unphased[loc] = None
        label = LOCUS_LABELS.get(loc, loc)
        if old_val and len(old_val) >= 2 and old_val[0] != old_val[1]:
            mask_summary.append(f"{label}(het:{'/'.join(old_val)})")
        elif old_val:
            mask_summary.append(f"{label}(hom:{old_val[0]})")
        else:
            mask_summary.append(f"{label}(already_missing)")

    return true_pair, unphased, mask_summary


def make_haplotype_label(hap_dict: dict) -> str:
    """Format a haplotype dict into a compact string."""
    parts = []
    for loc in LOCI:
        val = hap_dict.get(loc, '')
        # Compact JSON arrays
        if val and val.startswith('['):
            try:
                arr = json.loads(val)
                if len(arr) <= 2:
                    val = '/'.join(arr)
                else:
                    val = f"[{len(arr)} alleles]"
            except Exception:
                pass
        label = LOCUS_LABELS.get(loc, loc)
        parts.append(f"{label}={val}")
    return ' | '.join(parts)


def make_genotype_summary(unphased: dict) -> str:
    """Format unphased genotype into a human-readable summary string."""
    parts = []
    for loc in LOCI:
        alleles = unphased.get(loc)
        label = LOCUS_LABELS.get(loc, loc)
        if alleles is None:
            parts.append(f"{label}=<MISSING>")
        elif len(alleles) == 1:
            parts.append(f"{label}={alleles[0]} (hom)")
        else:
            parts.append(f"{label}={alleles[0]}/{alleles[1]}")
    return ' | '.join(parts)


def generate(
    n_patients: int = N_PATIENTS,
    output_path: str = None,
    seed: int = RANDOM_SEED,
) -> list:
    """Generate synthetic patients and write CSV."""
    out_path = Path(output_path or OUTPUT)
    rng = random.Random(seed)

    print(f"  Seed: {seed}")
    print(f"  Patients to generate: {n_patients}")
    print(f"  Mask range: {N_MASK_RANGE[0]}–{N_MASK_RANGE[1]} loci per patient")

    # Load reference
    haplotypes, freqs = connect_db(str(DB_PATH))

    patients = []

    for i in range(n_patients):
        pid = f"SYNTH_{i+1:04d}"

        true_pair, unphased, mask_summary = sample_patient(haplotypes, freqs, rng)

        # Build the "observed" JSON — serialise the unphased genotype
        # We store as a JSON string for the CSV
        observed_json = json.dumps(unphased, default=str)

        # Build the "true" haplotype JSON — both haplotypes with their alleles
        true_h1 = true_pair['haplotype_1']
        true_h2 = true_pair['haplotype_2']
        true_phased = {
            'haplotype_1': {loc: true_h1[loc] for loc in LOCI},
            'haplotype_2': {loc: true_h2[loc] for loc in LOCI},
        }
        true_json = json.dumps(true_phased, default=str)

        # Also create flat columns for easier parsing
        flat_row = {
            'patient_id': pid,
            'observed_json': observed_json,
            'true_phased_json': true_json,
        }

        # Add flat genotype columns (unphased observed)
        for loc in LOCI:
            alleles = unphased.get(loc)
            if alleles and len(alleles) >= 1:
                flat_row[f"{LOCUS_LABELS[loc]}_1"] = alleles[0]
                flat_row[f"{LOCUS_LABELS[loc]}_2"] = alleles[1] if len(alleles) > 1 else alleles[0]
            else:
                flat_row[f"{LOCUS_LABELS[loc]}_1"] = '.'
                flat_row[f"{LOCUS_LABELS[loc]}_2"] = '.'

        # Add flat truth columns (phased)
        for loc in LOCI:
            flat_row[f"TRUE_{LOCUS_LABELS[loc]}_1"] = true_h1.get(loc, '.')
            flat_row[f"TRUE_{LOCUS_LABELS[loc]}_2"] = true_h2.get(loc, '.')

        flat_row['num_masked'] = len(mask_summary)
        flat_row['mask_detail'] = '; '.join(mask_summary)

        patients.append({
            'pid': pid,
            'true_pair': true_pair,
            'unphased': unphased,
            'mask_summary': mask_summary,
            'flat_row': flat_row,
        })

    # ── Write CSV ──────────────────────────────────────────────────
    os.makedirs(out_path.parent, exist_ok=True)

    # Determine all possible column names from first patient
    first = patients[0]['flat_row']
    fieldnames = list(first.keys())

    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for p in patients:
            writer.writerow(p['flat_row'])

    file_size = out_path.stat().st_size

    print(f"\n  ✅ Generated {len(patients)} patients")
    print(f"  ✅ CSV written: {out_path} ({file_size:,} bytes)")

    return patients


def print_sample(patients: list, idx: int = 0):
    """Pretty-print one synthetic patient for verification."""
    p = patients[idx]
    print(f"\n{'='*72}")
    print(f"  SAMPLE PATIENT: {p['pid']}")
    print(f"{'='*72}")

    tp = p['true_pair']
    unphased = p['unphased']

    print(f"\n  🧬 Ground Truth (Phased):")
    h1_lbl = make_haplotype_label(tp['haplotype_1'])
    h2_lbl = make_haplotype_label(tp['haplotype_2'])
    print(f"    H1: {h1_lbl}")
    print(f"    H2: {h2_lbl}")
    print(f"    H1 freq (Global): {tp['h1_global_freq']:.6f}")
    print(f"    H2 freq (Global): {tp['h2_global_freq']:.6f}")
    print(f"    H1 != H2: {not tp['is_homozygous']}")

    print(f"\n  🎭 Observed Genotype (After Masking):")
    for loc in LOCI:
        alleles = unphased.get(loc)
        label = LOCUS_LABELS.get(loc, loc)
        if alleles is None:
            print(f"    {label:15}: <MISSING> (masked)")
        elif len(alleles) == 1:
            print(f"    {label:15}: {alleles[0]} (homozygous)")
        else:
            print(f"    {label:15}: {alleles[0]} / {alleles[1]}")

    print(f"\n  🎭 Loci masked: {len(p['mask_summary'])}")
    for m in p['mask_summary']:
        print(f"    • {m}")

    # HWE check: if H1 != H2, the expected product in the reference should be ~2pq
    if not tp['is_homozygous']:
        h1f = tp['h1_global_freq']
        h2f = tp['h2_global_freq']
        hwe_prob = 2 * h1f * h2f
        print(f"\n  📐 HWE check: 2×{h1f:.6f}×{h2f:.6f} = {hwe_prob:.8f}")

    # Show JSON payload
    print(f"\n  📦 Observed JSON:")
    observed_json = json.dumps(unphased, default=str, indent=2)
    print(f"    {observed_json}")

    print(f"\n  📦 Truth JSON (abbreviated):")
    true_summary = {
        'haplotype_1': {loc: tp['haplotype_1'][loc][:30] + ('...' if len(tp['haplotype_1'][loc]) > 30 else '')
                        for loc in LOCI},
        'haplotype_2': {loc: tp['haplotype_2'][loc][:30] + ('...' if len(tp['haplotype_2'][loc]) > 30 else '')
                        for loc in LOCI},
    }
    print(f"    {json.dumps(true_summary, indent=2)}")


def print_stats(patients: list):
    """Print aggregate statistics about the generated dataset."""
    n = len(patients)
    homozygotes = sum(1 for p in patients if p['true_pair']['is_homozygous'])
    het = n - homozygotes

    mask_counts = Counter(len(p['mask_summary']) for p in patients)

    # Per-locus missing rates
    locus_missing = Counter()
    for p in patients:
        for loc in LOCI:
            if p['unphased'].get(loc) is None:
                locus_missing[loc] += 1

    print(f"\n{'='*72}")
    print(f"  GENERATION STATISTICS")
    print(f"{'='*72}")
    print(f"  Total patients:                {n}")
    print(f"  Homozygous pairs:              {homozygotes} ({100*homozygotes/n:.1f}%)")
    print(f"  Heterozygous pairs:            {het} ({100*het/n:.1f}%)")
    print(f"\n  Masked loci distribution:")
    for k in sorted(mask_counts.keys()):
        print(f"    {k} locus masked: {mask_counts[k]} ({100*mask_counts[k]/n:.1f}%)")
    print(f"\n  Per-locus missing rate:")
    total_patients = n
    for loc in LOCI:
        label = LOCUS_LABELS.get(loc, loc)
        rate = locus_missing[loc] / total_patients * 100
        bar = '█' * int(rate / 5) + '░' * max(0, 20 - int(rate / 5))
        print(f"    {label:15} {locus_missing[loc]:4}/{total_patients} ({rate:5.1f}%) {bar}")


# ── Main ───────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 72)
    print("  HaploStats — Synthetic Patient Generator (Phase 6A)")
    print("=" * 72)

    patients = generate(
        n_patients=N_PATIENTS,
        output_path=str(OUTPUT),
        seed=RANDOM_SEED,
    )

    # Print one sample patient
    print_sample(patients, idx=0)

    # Print a second sample with different masking
    print_sample(patients, idx=len(patients)//3)

    # Print aggregate stats
    print_stats(patients)

    print(f"\n{'='*72}")
    print(f"  Generation complete. CSV ready at:")
    print(f"    {OUTPUT}")
    print(f"{'='*72}")
