#!/usr/bin/env python3
"""
HaploStats — Batch Benchmarker & Concordance Analyser
Phase 6: Asynchronous mass validation of the HaploStats imputation engine.

Pipeline
--------
  1. Read a large CSV of patient genotypes (unphased) + hidden truth set
  2. Fire async HTTP POST requests to the local /impute endpoint
  3. Compare each inferred top pair against the true phased high-resolution typing
  4. Write detailed per-patient results + aggregate concordance report

Expected input CSV columns:
  patient_id, ancestry,
  A1, A2, C1, C2, B1, B2, DRB345_1, DRB345_2,
  DRB1_1, DRB1_2, DQA1_1, DQA1_2, DQB1_1, DQB1_2,
  DPA1_1, DPA1_2, DPB1_1, DPB1_2,
  TRUE_A1, TRUE_A2, TRUE_C1, TRUE_C2, ... (truth: full 9-locus phased)
"""

import csv
import json
import sys
import os
import math
import time
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from typing import Optional

try:
    import aiohttp
    import asyncio
except ImportError:
    print("❌ aiohttp not installed. Run: pip install aiohttp")
    sys.exit(1)

# ── Config ─────────────────────────────────────────────────────────

BASE_DIR = Path(__file__).resolve().parent.parent
DEFAULT_INPUT  = BASE_DIR / "data" / "raw" / "Concordance_Truth_Set.csv"
DEFAULT_OUTPUT = BASE_DIR / "data" / "results" / "concordance_results.csv"
DEFAULT_REPORT = BASE_DIR / "data" / "results" / "concordance_report.txt"
API_URL_TEMPLATE = "http://localhost:{port}/impute?population={population}"

API_PORT = 8000
CONCURRENCY = 10  # max simultaneous API calls


# ── Locus definitions ──────────────────────────────────────────────

# 9 core HLA loci
LOCI = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
         'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

LOCUS_LABELS = {
    'hla_a': 'A', 'hla_c': 'C', 'hla_b': 'B',
    'hla_drb345': 'DRB345', 'hla_drb1': 'DRB1',
    'hla_dqa1': 'DQA1', 'hla_dqb1': 'DQB1',
    'hla_dpa1': 'DPA1', 'hla_dpb1': 'DPB1',
}

# CSV column names for the input genotype (unphased)
GENO_COLS = {
    'hla_a':      ('A1', 'A2'),
    'hla_c':      ('C1', 'C2'),
    'hla_b':      ('B1', 'B2'),
    'hla_drb345': ('DRB345_1', 'DRB345_2'),
    'hla_drb1':   ('DRB1_1', 'DRB1_2'),
    'hla_dqa1':   ('DQA1_1', 'DQA1_2'),
    'hla_dqb1':   ('DQB1_1', 'DQB1_2'),
    'hla_dpa1':   ('DPA1_1', 'DPA1_2'),
    'hla_dpb1':   ('DPB1_1', 'DPB1_2'),
}

# CSV columns for the truth set (phased high-resolution)
TRUTH_COLS = {
    'hla_a':      ('TRUE_A1', 'TRUE_A2'),
    'hla_c':      ('TRUE_C1', 'TRUE_C2'),
    'hla_b':      ('TRUE_B1', 'TRUE_B2'),
    'hla_drb345': ('TRUE_DRB345_1', 'TRUE_DRB345_2'),
    'hla_drb1':   ('TRUE_DRB1_1', 'TRUE_DRB1_2'),
    'hla_dqa1':   ('TRUE_DQA1_1', 'TRUE_DQA1_2'),
    'hla_dqb1':   ('TRUE_DQB1_1', 'TRUE_DQB1_2'),
    'hla_dpa1':   ('TRUE_DPA1_1', 'TRUE_DPA1_2'),
    'hla_dpb1':   ('TRUE_DPB1_1', 'TRUE_DPB1_2'),
}


# ── Data classes ───────────────────────────────────────────────────

@dataclass
class PatientRecord:
    """One row from the test CSV."""
    patient_id: str
    ancestry: str = "Global"
    genotype: dict = field(default_factory=dict)   # locus → [allele1, allele2]
    truth: dict = field(default_factory=dict)      # locus → [true_a1, true_a2]

    # Results from API
    top_haplotype_1: str = ""
    top_haplotype_2: str = ""
    top_posterior: float = 0.0
    total_pairs: int = 0
    entropy: float = 0.0
    api_status: str = "pending"
    api_error: str = ""
    response_time_ms: float = 0.0

    # Concordance
    locus_concordance: dict = field(default_factory=dict)
    full_concordance: bool = False
    phase_concordant: bool = False


@dataclass
class ConcordanceReport:
    """Aggregate statistics."""
    total_patients: int = 0
    api_success: int = 0
    api_failed: int = 0
    full_concordant: int = 0
    phase_concordant: int = 0
    no_match_found: int = 0
    locus_full_concordance: dict = field(default_factory=dict)
    locus_allele_concordance: dict = field(default_factory=dict)
    mean_response_time: float = 0.0
    p50_posterior: float = 0.0
    total_response_times: list = field(default_factory=list)
    posteriors: list = field(default_factory=list)
    per_patient: list = field(default_factory=list)


# ── Allele matching (same logic as em_algorithm.py) ────────────────

def _norm(a: str) -> str:
    return a.strip().replace(' ', '')

def _allele_match(pat: str, ref: str) -> bool:
    """Match patient allele against reference, handling JSON arrays and resolution."""
    pa = _norm(pat)
    ra = _norm(ref)
    if pa == ra:
        return True
    if ra.startswith('['):
        try:
            return any(_allele_match(pa, c) for c in json.loads(ra))
        except Exception:
            return False
    if pa.startswith('['):
        try:
            return any(_allele_match(c, ra) for c in json.loads(pa))
        except Exception:
            return False
    pp = pa.split(':')
    rp = ra.split(':')
    for a, b in zip(pp, rp):
        if a != b:
            return False
    return True


# ── CSV Parser ─────────────────────────────────────────────────────

def load_patients(csv_path: str) -> list[PatientRecord]:
    """Parse the input CSV into PatientRecord objects."""
    patients = []

    with open(csv_path, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = set(reader.fieldnames or [])

        for row_num, row in enumerate(reader, start=2):
            patient_id = row.get('patient_id', '').strip()
            if not patient_id:
                patient_id = f"PAT_{row_num}"

            ancestry = row.get('ancestry', 'Global').strip() or 'Global'

            # Build unphased genotype
            genotype = {}
            for loc in LOCI:
                c1, c2 = GENO_COLS[loc]
                a1 = row.get(c1, '').strip()
                a2 = row.get(c2, '').strip()
                alleles = []
                if a1 and a1 != '.' and a1 != '-':
                    alleles.append(a1)
                if a2 and a2 != '.' and a2 != '-':
                    alleles.append(a2)
                genotype[loc] = alleles if alleles else None

            # Build truth
            truth = {}
            has_truth = False
            for loc in LOCI:
                t1, t2 = TRUTH_COLS[loc]
                ta1 = row.get(t1, '').strip()
                ta2 = row.get(t2, '').strip()
                truth_all = []
                if ta1 and ta1 != '.' and ta1 != '-':
                    truth_all.append(ta1)
                if ta2 and ta2 != '.' and ta2 != '-':
                    truth_all.append(ta2)
                truth[loc] = truth_all if truth_all else None
                if truth[loc] is not None:
                    has_truth = True

            if not has_truth:
                truth = None  # No truth available — skip concordance

            patients.append(PatientRecord(
                patient_id=patient_id,
                ancestry=ancestry,
                genotype=genotype,
                truth=truth,
            ))

    return patients


# ── API Payload Builder ────────────────────────────────────────────

def build_payload(patient: PatientRecord) -> dict:
    """Convert a PatientRecord to the /impute JSON payload."""
    payload = {}
    for loc in LOCI:
        alleles = patient.genotype.get(loc)
        if alleles and len(alleles) > 0:
            payload[loc] = alleles
        else:
            payload[loc] = None
    return payload


# ── Concordance Logic ──────────────────────────────────────────────

def _parse_haplotype_label(label: str) -> dict:
    """Parse a formatted haplotype label back into locus→allele dict.
    e.g. 'HLA-A=A*01:01 | HLA-C=C*07:01 | ...' → {'hla_a': 'A*01:01', ...}
    """
    result = {}
    parts = label.split('|')
    for part in parts:
        part = part.strip()
        if '=' not in part:
            continue
        key, val = part.split('=', 1)
        key = key.strip().lower().replace('-', '_')
        # Map display keys back to internal locus names
        display_to_internal = {
            'hla_a': 'hla_a', 'hla_c': 'hla_c', 'hla_b': 'hla_b',
            'drb345': 'hla_drb345', 'drb1': 'hla_drb1',
            'dqa1': 'hla_dqa1', 'dqb1': 'hla_dqb1',
            'dpa1': 'hla_dpa1', 'dpb1': 'hla_dpb1',
        }
        internal_key = display_to_internal.get(key)
        if internal_key:
            result[internal_key] = val.strip()
    return result


def compute_concordance(patient: PatientRecord) -> PatientRecord:
    """
    Compare the API's top inferred haplotype pair against the truth set.

    For each locus, checks:
    - The inferred haplotypes' alleles match the true alleles (order-agnostic)
    - Phase concordance: H1 carries the correct chr6 copy

    Returns the PatientRecord with concordance fields populated.
    """
    if not patient.truth:
        patient.full_concordance = False
        patient.phase_concordant = False
        patient.locus_concordance = {}
        return patient

    if not patient.top_haplotype_1 or not patient.top_haplotype_2:
        patient.full_concordance = False
        patient.phase_concordant = False
        return patient

    # Parse the inferred haplotypes
    h1_alleles = _parse_haplotype_label(patient.top_haplotype_1)
    h2_alleles = _parse_haplotype_label(patient.top_haplotype_2)

    locus_conc = {}
    all_match = True          # all loci covered (order-agnostic)
    specific_phase = True     # H1→True_H1 AND H2→True_H2
    any_phase_possible = True  # either orientation works

    for loc in LOCI:
        true_alleles = patient.truth.get(loc)
        if not true_alleles or len(true_alleles) < 2:
            continue  # no truth for this locus, skip

        true_a1 = _norm(true_alleles[0])
        true_a2 = _norm(true_alleles[1])

        inf_a1 = _norm(h1_alleles.get(loc, ''))
        inf_a2 = _norm(h2_alleles.get(loc, ''))

        if not inf_a1 or not inf_a2:
            locus_conc[loc] = 0.0
            all_match = False
            specific_phase = False
            any_phase_possible = False
            continue

        # Order-agnostic: does H1 ∪ H2 cover {true_a1, true_a2}?
        cov_a1 = _allele_match(true_a1, inf_a1) or _allele_match(true_a1, inf_a2)
        cov_a2 = _allele_match(true_a2, inf_a1) or _allele_match(true_a2, inf_a2)

        if cov_a1 and cov_a2:
            locus_conc[loc] = 1.0
        else:
            locus_conc[loc] = 0.0
            all_match = False
            specific_phase = False
            any_phase_possible = False
            continue

        # Specific phase check: does TRUE_H1 = H1 AND TRUE_H2 = H2?
        h1_matches_true_a1 = _allele_match(true_a1, inf_a1)
        h2_matches_true_a2 = _allele_match(true_a2, inf_a2)
        if not (h1_matches_true_a1 and h2_matches_true_a2):
            specific_phase = False

    patient.locus_concordance = locus_conc
    patient.full_concordance = all_match
    patient.phase_concordant = specific_phase

    return patient


# ── Async API Client ───────────────────────────────────────────────

async def send_impute_request(
    session: aiohttp.ClientSession,
    patient: PatientRecord,
    semaphore: asyncio.Semaphore,
    port: int = API_PORT,
) -> PatientRecord:
    """Send a single patient to /impute and parse the response."""
    url = API_URL_TEMPLATE.format(port=port, population=patient.ancestry)
    payload = build_payload(patient)

    async with semaphore:
        t_start = time.monotonic()
        try:
            async with session.post(url, json=payload, timeout=aiohttp.ClientTimeout(total=120)) as resp:
                t_elapsed = (time.monotonic() - t_start) * 1000  # ms
                patient.response_time_ms = round(t_elapsed, 1)

                if resp.status == 200:
                    data = await resp.json()
                    patient.api_status = "success"
                    top_pairs = data.get('top_pairs', [])
                    if top_pairs:
                        top = top_pairs[0]
                        # The API response has haplotype_1 and haplotype_2 as strings
                        # Check the structure
                        ht1 = top.get('haplotype_1', '')
                        ht2 = top.get('haplotype_2', '')
                        if not ht1:
                            # maybe returned as model with haplotype_1
                            ht1 = top.get('haplotype_1', '')
                        patient.top_haplotype_1 = ht1
                        patient.top_haplotype_2 = ht2
                        patient.top_posterior = top.get('posterior', 0.0)
                    else:
                        patient.top_posterior = 0.0
                    patient.total_pairs = data.get('total_possible_pairs', 0)
                    patient.entropy = data.get('entropy', 0.0)
                else:
                    patient.api_status = "error"
                    patient.api_error = f"HTTP {resp.status}: {await resp.text()}"
        except asyncio.TimeoutError:
            patient.api_status = "error"
            patient.api_error = "timeout"
        except aiohttp.ClientError as e:
            patient.api_status = "error"
            patient.api_error = str(e)
        except Exception as e:
            patient.api_status = "error"
            patient.api_error = f"unexpected: {e}"
            t_elapsed = (time.monotonic() - t_start) * 1000
            patient.response_time_ms = round(t_elapsed, 1)

    return patient


async def run_batch(
    patients: list[PatientRecord],
    port: int = API_PORT,
    concurrency: int = CONCURRENCY,
) -> list[PatientRecord]:
    """Fire all patients asynchronously against the API."""
    semaphore = asyncio.Semaphore(concurrency)
    connector = aiohttp.TCPConnector(limit=concurrency, limit_per_host=concurrency)

    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = [
            send_impute_request(session, p, semaphore, port)
            for p in patients
        ]
        results = await asyncio.gather(*tasks)

    return list(results)


# ── Results Writer ─────────────────────────────────────────────────

def write_results_csv(patients: list[PatientRecord], output_path: str):
    """Write per-patient results to CSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fieldnames = [
        'patient_id', 'ancestry',
        'api_status', 'api_error', 'response_time_ms',
        'top_posterior', 'total_pairs', 'entropy',
        'top_haplotype_1', 'top_haplotype_2',
        'full_concordance', 'phase_concordant',
    ]
    # Add per-locus concordance
    for loc in LOCI:
        fieldnames.append(f'conc_{LOCUS_LABELS[loc]}')

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()

        for p in patients:
            row = asdict(p)
            # Flatten locus concordance
            for loc in LOCI:
                row[f'conc_{LOCUS_LABELS[loc]}'] = p.locus_concordance.get(loc, '')
            writer.writerow(row)

    print(f"  ✅ Results written to {output_path}")


def generate_report(patients: list[PatientRecord]) -> ConcordanceReport:
    """Compute aggregate concordance statistics."""
    report = ConcordanceReport()
    report.total_patients = len(patients)

    times = []
    posteriors = []
    matched_patients = [p for p in patients if p.truth is not None]
    report.total_response_times = []
    report.posteriors = []

    for p in patients:
        if p.api_status == "success":
            report.api_success += 1
            times.append(p.response_time_ms)
            posteriors.append(p.top_posterior)
            report.total_response_times.append(p.response_time_ms)
            report.posteriors.append(p.top_posterior)

            if p.truth:
                if p.full_concordance:
                    report.full_concordant += 1
                if p.phase_concordant:
                    report.phase_concordant += 1
                if p.top_posterior == 0.0 and not p.top_haplotype_1:
                    report.no_match_found += 1
        else:
            report.api_failed += 1

        report.per_patient.append(p)

    # Response time stats
    if times:
        times.sort()
        report.mean_response_time = round(sum(times) / len(times), 1)

    # Posterior stats
    if posteriors:
        posteriors.sort()
        report.p50_posterior = posteriors[len(posteriors) // 2]

    # Per-locus concordance
    for loc in LOCI:
        conc_values = [
            p.locus_concordance.get(loc, 0)
            for p in matched_patients
            if p.locus_concordance and loc in p.locus_concordance
        ]
        if conc_values:
            report.locus_full_concordance[loc] = round(
                sum(1 for v in conc_values if v >= 0.999) / len(conc_values), 4
            )
            report.locus_allele_concordance[loc] = round(
                sum(conc_values) / len(conc_values), 4
            )

    return report


def print_report(report: ConcordanceReport):
    """Pretty-print the concordance report to console."""
    lines = []
    lines.append("=" * 72)
    lines.append("  HaploStats — Concordance Benchmark Report")
    lines.append("=" * 72)
    lines.append(f"")
    lines.append(f"  Patients processed:       {report.total_patients}")
    lines.append(f"  API success:              {report.api_success}")
    lines.append(f"  API failed:               {report.api_failed}")
    lines.append(f"  Full 9-locus concordant:  {report.full_concordant} "
                 f"({round(100*report.full_concordant/report.api_success, 2) if report.api_success else 0}%)")
    lines.append(f"  Phase concordant:         {report.phase_concordant} "
                 f"({round(100*report.phase_concordant/report.api_success, 2) if report.api_success else 0}%)")
    lines.append(f"  No match found:           {report.no_match_found}")
    lines.append(f"")
    lines.append(f"  Mean response time:       {report.mean_response_time} ms")
    lines.append(f"  Median posterior:         {report.p50_posterior:.6f}")
    lines.append(f"")
    lines.append(f"  ── Per-Locus Allele Concordance ──")
    for loc in LOCI:
        label = LOCUS_LABELS[loc]
        full = report.locus_full_concordance.get(loc, 0)
        allele = report.locus_allele_concordance.get(loc, 0)
        lines.append(f"    {label:10}  full={full:.2%}  allele-wise={allele:.2%}")
    lines.append(f"")
    lines.append("=" * 72)

    report_str = '\n'.join(lines)
    print(report_str)

    # Also write to file
    os.makedirs(os.path.dirname(DEFAULT_REPORT), exist_ok=True)
    with open(DEFAULT_REPORT, 'w') as f:
        f.write(report_str)
    print(f"  Report saved to {DEFAULT_REPORT}")


# ── Main Entrypoint ────────────────────────────────────────────────

def main(input_csv: str = None, output_csv: str = None,
         port: int = API_PORT, concurrency: int = CONCURRENCY,
         run_async: bool = True):
    """
    Main benchmark orchestrator.

    Parameters
    ----------
    input_csv  : path to the concordance truth set CSV
    output_csv : path for per-patient results
    port       : the local API server port
    concurrency: max simultaneous async requests
    run_async  : if False, only parse and print stats without calling API
    """
    csv_path = input_csv or str(DEFAULT_INPUT)
    out_path = output_csv or str(DEFAULT_OUTPUT)

    # ── Parse CSV ──────────────────────────────────────────────────
    print(f"\n  📂 Loading patients from: {csv_path}")
    if not os.path.exists(csv_path):
        print(f"  ❌ File not found: {csv_path}")
        print(f"  ⏳ Waiting for Concordance_Truth_Set.csv upload...")
        return

    patients = load_patients(csv_path)
    print(f"  ✅ Loaded {len(patients)} patient records")

    # Stats on typing resolution
    locus_stats = defaultdict(lambda: {'typed': 0, 'untyped': 0, 'with_truth': 0})
    for p in patients:
        for loc in LOCI:
            if p.genotype.get(loc):
                locus_stats[loc]['typed'] += 1
            else:
                locus_stats[loc]['untyped'] += 1
            if p.truth and p.truth.get(loc):
                locus_stats[loc]['with_truth'] += 1

    print(f"\n  ── Locus Typing Coverage ──")
    for loc in LOCI:
        s = locus_stats[loc]
        total = len(patients)
        print(f"    {LOCUS_LABELS[loc]:10} typed={s['typed']:3}/{total}  "
              f"truth={s['with_truth']:3}/{total}")

    has_truth = sum(1 for p in patients if p.truth)
    print(f"  Patients with truth data: {has_truth}/{len(patients)}")

    if not run_async:
        print(f"\n  📋 Script ready. Pass run_async=True to execute.")
        return

    # ── Fire async API calls ───────────────────────────────────────
    print(f"\n  🚀 Firing {len(patients)} patients at "
          f"localhost:{port} (concurrency={concurrency})...")
    t0 = time.monotonic()

    results = asyncio.run(run_batch(patients, port, concurrency))

    elapsed = time.monotonic() - t0
    print(f"  ⏱  Elapsed: {elapsed:.1f}s ({(elapsed/len(results)):.2f}s/patient avg)")

    # ── Compute concordance ────────────────────────────────────────
    print(f"  🔬 Computing concordance...")
    for p in results:
        compute_concordance(p)

    # ── Write results ──────────────────────────────────────────────
    write_results_csv(results, out_path)

    # ── Generate report ────────────────────────────────────────────
    report = generate_report(results)
    print_report(report)

    return report


# ── CLI Entry ──────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="HaploStats Batch Benchmarker & Concordance Analyser"
    )
    parser.add_argument("--input", "-i", help="Input CSV path",
                        default=str(DEFAULT_INPUT))
    parser.add_argument("--output", "-o", help="Output CSV path",
                        default=str(DEFAULT_OUTPUT))
    parser.add_argument("--port", "-p", type=int, default=API_PORT,
                        help=f"API server port (default: {API_PORT})")
    parser.add_argument("--concurrency", "-c", type=int, default=CONCURRENCY,
                        help=f"Max concurrent requests (default: {CONCURRENCY})")
    parser.add_argument("--dry-run", "-n", action="store_true",
                        help="Parse CSV only, don't call API")

    args = parser.parse_args()

    main(
        input_csv=args.input,
        output_csv=args.output,
        port=args.port,
        concurrency=args.concurrency,
        run_async=not args.dry_run,
    )
