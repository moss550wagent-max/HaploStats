#!/usr/bin/env python3
"""
HaploStats — DPB1-Specific Statistical Model (Phase 7)
Addresses the DQ-DP recombination hotspot by decoupling DPB1 imputation
from the strict 9-locus chain and using DPA1-conditional probabilities.

Key changes to HaploEM:
  1. Load DPB1 marginal frequencies + DPA1→DPB1 conditional lookup table
  2. In block 3 (hla_dpa1, hla_dpb1), decouple DP from the core chain
  3. For masked DPB1: impute using P(DPB1 | DPA1) or marginal P(DPB1)
  4. Score pairs using core posterior × DP contribution
"""

import sqlite3
import json
import math
import sys
import os
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH = BASE_DIR / "db" / "haplostats.db"

# ── Loci ───────────────────────────────────────────────────────────

ALL_LOCI = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
            'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

LOCUS_LABEL = {
    'hla_a': 'HLA-A', 'hla_c': 'HLA-C', 'hla_b': 'HLA-B',
    'hla_drb345': 'DRB345', 'hla_drb1': 'DRB1',
    'hla_dqa1': 'DQA1', 'hla_dqb1': 'DQB1',
    'hla_dpa1': 'DPA1', 'hla_dpb1': 'DPB1',
}

# Blocks: block 3 now splits DPA1/DPB1 into a separate DP step
CORE_BLOCKS = [
    ('block_1_core',     ['hla_a', 'hla_b', 'hla_c']),
    ('block_2_drdq',     ['hla_drb1', 'hla_dqb1', 'hla_drb345', 'hla_dqa1']),
]

DP_BLOCK = ('block_3_dp', ['hla_dpa1', 'hla_dpb1'])


# ── Allele matching ────────────────────────────────────────────────

def _norm(a):
    return a.strip().replace(' ', '') if a else ''

def _allele_match(pat, ref):
    pa, ra = _norm(pat), _norm(ref)
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
    pp, rp = pa.split(':'), ra.split(':')
    for a, b in zip(pp, rp):
        if a != b:
            return False
    return True


# ════════════════════════════════════════════════════════════════════
# HaploEM — DP-Enhanced Version
# ════════════════════════════════════════════════════════════════════

class HaploEM:
    """
    EM-based haplotype imputation with decoupled DPB1 model.
    """

    def __init__(self, db_path=None, population='Global'):
        self.db_path = str(db_path or DB_PATH)
        self.population = population
        self.freq_col = f'{population}_freq'
        self.conn = None
        self.full_ref = None

        # DPB1 special tables
        self.dpb1_marginal = {}       # dpb1_allele → frequency
        self.dpa1_to_dpb1 = {}        # dpa1_allele → {dpb1_allele: freq}
        self.dpa1_marginal = {}       # dpa1_allele → frequency
        self.dpb1_total_entries = 0   # number of unique DPB1 across ref
        self.dpa1_total_entries = 0   # number of unique DPA1 across ref

    def connect(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._load_reference()
        self._build_dp_tables()
        return self

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    # ── Reference loading ──────────────────────────────────────────

    def _load_reference(self):
        cur = self.conn.cursor()
        cols = ', '.join(ALL_LOCI + [self.freq_col])
        rows = cur.execute(
            f"SELECT {cols} FROM haplotypes "
            f"WHERE {self.freq_col} IS NOT NULL"
        ).fetchall()
        self.full_ref = [dict(r) for r in rows]
        print(f"  [HaploEM-DP] Loaded {len(self.full_ref)} haplotypes "
              f"(pop={self.population})")

    def _build_dp_tables(self):
        """
        Build DPA1→DPB1 conditional probability tables from the reference.

        For each DPA1 allele in the reference, count the DPB1 alleles it
        appears with and weight by haplotype frequency.
        """
        dpa1_counts = defaultdict(float)      # DPA1 → total frequency
        dpa1_dpb1_counts = defaultdict(lambda: defaultdict(float))  # DPA1 → DPB1 → freq
        dpb1_counts = defaultdict(float)      # DPB1 → total frequency

        for h in self.full_ref:
            freq = h.get(self.freq_col, 0.0) or 0.0
            if freq <= 0:
                continue

            dpa1 = h.get('hla_dpa1', '') or ''
            dpb1 = h.get('hla_dpb1', '') or ''

            # Handle ambiguous alleles (JSON arrays) — each possibility
            # gets a fraction of the frequency
            dpa1_alleles = self._expand_ambiguity(dpa1)
            dpb1_alleles = self._expand_ambiguity(dpb1)

            for da in dpa1_alleles:
                weight = freq / len(dpa1_alleles)
                dpa1_counts[da] += weight
                for db in dpb1_alleles:
                    sub_weight = weight / len(dpb1_alleles)
                    dpa1_dpb1_counts[da][db] += sub_weight

            for db in dpb1_alleles:
                dpb1_counts[db] += freq / len(dpb1_alleles)

        # Normalise to conditional probabilities
        # P(DPB1 = X | DPA1 = Y) = count(Y→X) / count(Y)
        self.dpa1_to_dpb1 = {}
        for dpa1_allele, dpb1_map in dpa1_dpb1_counts.items():
            total = dpa1_counts.get(dpa1_allele, 0.0)
            if total > 0:
                cond = {db: cnt / total for db, cnt in dpb1_map.items()}
                # Sort by probability descending
                self.dpa1_to_dpb1[dpa1_allele] = dict(
                    sorted(cond.items(), key=lambda x: -x[1])
                )

        # Normalise marginal DPB1 frequencies
        total_dpb1 = sum(dpb1_counts.values())
        if total_dpb1 > 0:
            self.dpb1_marginal = {k: v / total_dpb1 for k, v in dpb1_counts.items()}
            self.dpb1_marginal = dict(
                sorted(self.dpb1_marginal.items(), key=lambda x: -x[1])
            )

        # DPA1 marginal
        total_dpa1 = sum(dpa1_counts.values())
        if total_dpa1 > 0:
            self.dpa1_marginal = {k: v / total_dpa1 for k, v in dpa1_counts.items()}

        self.dpb1_total_entries = len(self.dpb1_marginal)
        self.dpa1_total_entries = len(self.dpa1_marginal)

        print(f"  [HaploEM-DP] DPA1 alleles: {self.dpa1_total_entries}")
        print(f"  [HaploEM-DP] DPB1 alleles: {self.dpb1_total_entries}")
        print(f"  [HaploEM-DP] DPA1→DPB1 conditional tables: "
              f"{len(self.dpa1_to_dpb1)} entries")

    @staticmethod
    def _expand_ambiguity(allele_str):
        """If the allele is a JSON array, return all elements.
        Otherwise return [allele_str]."""
        if not allele_str:
            return ['']
        if allele_str.startswith('['):
            try:
                return json.loads(allele_str)
            except Exception:
                return [allele_str]
        return [allele_str]

    @staticmethod
    def _hap_key(h):
        return tuple(h.get(l, '') for l in ALL_LOCI)

    def _make_label(self, idx):
        h = self.full_ref[idx]
        parts = []
        for loc in ALL_LOCI:
            val = h.get(loc, '')
            if val.startswith('['):
                try:
                    arr = json.loads(val)
                    val = '/'.join(arr[:2])
                    if len(arr) > 2:
                        val += f'/...({len(arr)})'
                except Exception:
                    pass
            parts.append(f"{LOCUS_LABEL[loc]}={val}")
        return ' | '.join(parts)

    # ── DPB1 imputation ────────────────────────────────────────────

    def impute_dpb1_for_pair(self, dpa1_allele_h1, dpa1_allele_h2,
                             patient_dpa1=None, patient_dpb1=None):
        """
        Given DPA1 alleles from two candidate haplotypes and patient
        typing, return the best DPB1 alleles and a confidence weight.

        Returns: (dpb1_h1, dpb1_h2, dpb1_weight, weight_detail)
          dpb1_weight: 0-1 multiplier to apply to the pair's posterior
          weight_detail: explanation string
        """
        # Normalise
        dpa1_h1 = _norm(dpa1_allele_h1) if dpa1_allele_h1 else ''
        dpa1_h2 = _norm(dpa1_allele_h2) if dpa1_allele_h2 else ''

        # Case 1: Patient has DPB1 typed — use directly
        if patient_dpb1:
            # Check if the patient's DPB1 matches either haplotype's DPA1
            # Actually the patient's DPB1 is the ground truth here
            # We accept whatever the patient has
            return (patient_dpb1, patient_dpb1, 1.0,
                    "patient typed at DPB1")

        # Case 2: Patient has DPA1 typed — use P(DPB1 | DPA1)
        if patient_dpa1:
            best_dpb1_h1 = None
            best_dpb1_h2 = None
            best_weight_h1 = 0.0
            best_weight_h2 = 0.0

            for pat_dpa1 in patient_dpa1:
                # Try matching H1's DPA1
                if _allele_match(pat_dpa1, dpa1_h1):
                    cond = self.dpa1_to_dpb1.get(dpa1_h1, self.dpb1_marginal)
                    if cond:
                        top = list(cond.items())[0]
                        if top[1] > best_weight_h1:
                            best_dpb1_h1 = top[0]
                            best_weight_h1 = top[1]

                # Try matching H2's DPA1
                if _allele_match(pat_dpa1, dpa1_h2) and dpa1_h2 != dpa1_h1:
                    cond = self.dpa1_to_dpb1.get(dpa1_h2, self.dpb1_marginal)
                    if cond:
                        top = list(cond.items())[0]
                        if top[1] > best_weight_h2:
                            best_dpb1_h2 = top[0]
                            best_weight_h2 = top[1]

            if best_dpb1_h1 and best_dpb1_h2:
                avg_weight = (best_weight_h1 + best_weight_h2) / 2
                return (best_dpb1_h1, best_dpb1_h2, avg_weight,
                        f"DPA1-conditional: {avg_weight:.3f}")

            # Fallback: use marginal for unmatched DPA1
            if not best_dpb1_h1:
                if self.dpb1_marginal:
                    top = list(self.dpb1_marginal.items())[0]
                    best_dpb1_h1 = top[0]
                    best_weight_h1 = top[1]
            if not best_dpb1_h2:
                if self.dpb1_marginal:
                    top = list(self.dpb1_marginal.items())[0]
                    best_dpb1_h2 = top[0]
                    best_weight_h2 = top[1]

            if best_dpb1_h1 and best_dpb1_h2:
                avg_weight = (best_weight_h1 + best_weight_h2) / 2
                return (best_dpb1_h1, best_dpb1_h2, max(avg_weight, 0.01),
                        f"partial DPA1-cond: {avg_weight:.3f} (marginal fallback)")

        # Case 3: Neither DPB1 nor DPA1 typed — use marginal frequencies
        if self.dpb1_marginal:
            top = list(self.dpb1_marginal.items())[0]
            top2 = list(self.dpb1_marginal.items())[1] if len(self.dpb1_marginal) > 1 else top
            marginal_weight = top[1]
            return (top[0], top2[0], max(marginal_weight, 0.01),
                    f"marginal DPB1: {marginal_weight:.3f}")

        return (None, None, 0.0, "no DP data available")

    # ── Progressive insertion (modified block 3) ───────────────────

    def progressive_em(self, patient_genotype, trim_threshold=1e-4,
                       max_em_iterations=50):
        """Run EM with decoupled DPB1 at block 3."""
        gt = {}
        for loc in ALL_LOCI:
            val = patient_genotype.get(loc)
            if val is None or val == '' or (isinstance(val, list) and len(val) == 0):
                gt[loc] = None
            elif isinstance(val, str):
                gt[loc] = [val]
            else:
                gt[loc] = list(val)

        known = sum(1 for v in gt.values() if v is not None)
        print(f"\n  Patient typed loci: {known} / 9")
        for loc in ALL_LOCI:
            v = gt.get(loc)
            if v:
                print(f"    {LOCUS_LABEL[loc]:10} {v[0]:30} / {v[1] if len(v) > 1 else v[0]}")
            else:
                print(f"    {LOCUS_LABEL[loc]:10} <missing>")

        # Initial frequencies
        self.current_freqs = {}
        for h in self.full_ref:
            key = self._hap_key(h)
            self.current_freqs[key] = h.get(self.freq_col, 0.0)

        surviving_indices = set(range(len(self.full_ref)))
        block_log = []

        # ── Blocks 1-2: Core A/B/C/DR/DQ (unchanged) ──────────────
        for block_name, block_loci in CORE_BLOCKS:
            print(f"\n  ── {block_name.upper()} ({', '.join(block_loci)}) ──")
            known_in_block = [l for l in block_loci if gt.get(l) is not None]
            # Constrained = all loci from this block AND all earlier blocks
            constrained_so_far = []
            for bn, bl in CORE_BLOCKS:
                constrained_so_far.extend(bl)
                if bn == block_name:
                    break

            new_surviving = set()
            for idx in surviving_indices:
                h = self.full_ref[idx]
                ok = True
                for loc in known_in_block:
                    pat = gt[loc]
                    hap = h.get(loc, '')
                    if not hap:
                        ok = False
                        break
                    if not any(_allele_match(a, hap) for a in pat):
                        ok = False
                        break
                if ok:
                    new_surviving.add(idx)
            surviving_indices = new_surviving
            print(f"    Haplotypes matching: {len(surviving_indices)}")
            if len(surviving_indices) < 2:
                break

            pairs = self._build_pairs(list(surviving_indices), gt, constrained_so_far)
            print(f"    Pairs formed: {len(pairs)}")
            if len(pairs) < 1:
                break

            pairs = self._em_loop(pairs, max_iterations=max_em_iterations)
            pre = len(pairs)
            pairs = [p for p in pairs if p['posterior'] >= trim_threshold]
            print(f"    EM done. Pairs: {pre} → {len(pairs)} after trim")

            used = set()
            for p in pairs:
                used.add(p['h1_idx']); used.add(p['h2_idx'])
            surviving_indices = used
            block_log.append({
                'block': block_name, 'haplotypes': len(surviving_indices),
                'pairs_before_em': pre, 'pairs_after_trim': len(pairs),
                'converged_iterations': pairs[0].get('em_iterations', 0) if pairs else 0,
            })

        # ── Block 3: DP decoupled from core ────────────────────────
        print(f"\n  ── {DP_BLOCK[0].upper()} ({', '.join(DP_BLOCK[1])}) [DECOUPLED] ──")

        patient_dpa1 = gt.get('hla_dpa1')
        patient_dpb1 = gt.get('hla_dpb1')

        # Filter survivors to match DPA1 if known
        if patient_dpa1:
            new_surviving = set()
            for idx in surviving_indices:
                h = self.full_ref[idx]
                hap_dpa1 = h.get('hla_dpa1', '')
                if not hap_dpa1:
                    continue
                if any(_allele_match(a, hap_dpa1) for a in patient_dpa1):
                    new_surviving.add(idx)
            surviving_indices = new_surviving
            print(f"    Haplotypes matching DPA1: {len(surviving_indices)}")

        # Build pairs from survivors (constrained to core + DPA1)
        dp_constrained = []
        for bn, bl in CORE_BLOCKS:
            dp_constrained.extend(bl)
        dp_constrained.append('hla_dpa1')

        pairs = self._build_pairs(list(surviving_indices), gt, dp_constrained)
        print(f"    Core+DPA1 pairs formed: {len(pairs)}")

        if len(pairs) > 0:
            # Run EM on core+DPA1
            pairs = self._em_loop(pairs, max_iterations=max_em_iterations)

            # Now impute DPB1 for each pair using the decoupled model
            for p in pairs:
                h1 = self.full_ref[p['h1_idx']]
                h2 = self.full_ref[p['h2_idx']]
                dpa1_h1 = h1.get('hla_dpa1', '')
                dpa1_h2 = h2.get('hla_dpa1', '')
                dpb1_h1_truth = h1.get('hla_dpb1', '')
                dpb1_h2_truth = h2.get('hla_dpb1', '')

                # Impute DPB1 from DPA1
                imp_dpb1_h1, imp_dpb1_h2, dp_weight, dp_detail = \
                    self.impute_dpb1_for_pair(
                        dpa1_h1, dpa1_h2,
                        patient_dpa1, patient_dpb1
                    )

                # Apply DP weight to posterior
                p['dp_weight'] = dp_weight
                p['dp_detail'] = dp_detail
                p['imputed_dpb1_h1'] = imp_dpb1_h1
                p['imputed_dpb1_h2'] = imp_dpb1_h2
                p['true_dpb1_h1'] = dpb1_h1_truth
                p['true_dpb1_h2'] = dpb1_h2_truth
                p['posterior'] *= dp_weight

            # Re-normalise after DP weighting
            total = sum(p['posterior'] for p in pairs)
            if total > 0:
                for p in pairs:
                    p['posterior'] /= total

            pre = len(pairs)
            pairs = [p for p in pairs if p['posterior'] >= trim_threshold]
            print(f"    DP imputation done. Pairs: {pre} → {len(pairs)}")
        else:
            print(f"    ⚠️ No pairs survived DP filtering")

        block_log.append({
            'block': DP_BLOCK[0], 'haplotypes': len(surviving_indices),
            'pairs_before_em': len(pairs), 'pairs_after_trim': len(pairs),
            'converged_iterations': pairs[0].get('em_iterations', 0) if pairs else 0,
        })

        # ── Final sort ─────────────────────────────────────────────
        pairs.sort(key=lambda p: p['posterior'], reverse=True)

        entropy = 0.0
        for p in pairs:
            if p['posterior'] > 0:
                entropy -= p['posterior'] * math.log2(p['posterior'])

        top3 = []
        cum = 0.0
        for rank, p in enumerate(pairs[:3], 1):
            cum += p['posterior']
            top3.append({
                'rank': rank,
                'haplotype_1': self._make_label(p['h1_idx']),
                'haplotype_2': self._make_label(p['h2_idx']),
                'posterior': round(p['posterior'], 8),
                'cumulative': round(cum, 8),
                'h1_frequency': round(p['h1_freq'], 6),
                'h2_frequency': round(p['h2_freq'], 6),
                'dp_weight': round(p.get('dp_weight', 1.0), 4),
                'dp_detail': p.get('dp_detail', ''),
                'imputed_dpb1_h1': p.get('imputed_dpb1_h1', ''),
                'imputed_dpb1_h2': p.get('imputed_dpb1_h2', ''),
                'true_dpb1_h1': p.get('true_dpb1_h1', ''),
                'true_dpb1_h2': p.get('true_dpb1_h2', ''),
            })

        print(f"\n  {'='*68}")
        print(f"  FINAL: {len(pairs)} pairs | Entropy: {entropy:.4f} bits")

        return {
            'population': self.population,
            'patient_genotype': {k: v for k, v in gt.items()},
            'blocks': block_log,
            'total_pairs_final': len(pairs),
            'entropy': round(entropy, 4),
            'top_3': top3,
            'pairs': pairs,
        }

    # ── Pair builder ───────────────────────────────────────────────

    def _build_pairs(self, indices, gt, constrained_loci):
        pairs = []
        n = len(indices)
        for i in range(n):
            h1 = self.full_ref[indices[i]]
            for j in range(i, n):
                h2 = self.full_ref[indices[j]]
                valid = True
                for loc in constrained_loci:
                    pat = gt.get(loc)
                    if pat is None:
                        continue
                    a1, a2 = pat[0], pat[1] if len(pat) >= 2 else pat[0]
                    h1v, h2v = h1.get(loc, ''), h2.get(loc, '')
                    c1 = _allele_match(a1, h1v) or _allele_match(a1, h2v)
                    c2 = _allele_match(a2, h1v) or _allele_match(a2, h2v)
                    if not (c1 and c2):
                        valid = False
                        break
                if not valid:
                    continue
                f1 = h1.get(self.freq_col, 0.0) or 0.0
                f2 = h2.get(self.freq_col, 0.0) or 0.0
                if f1 <= 0 or f2 <= 0:
                    continue
                jp = f1 * f1 if i == j else 2.0 * f1 * f2
                pairs.append({
                    'h1_idx': indices[i], 'h2_idx': indices[j],
                    'h1_freq': f1, 'h2_freq': f2,
                    'joint_prob': jp, 'k_factor': 1.0 if i == j else 2.0,
                    'posterior': 0.0,
                })
        return pairs

    # ── EM loop ────────────────────────────────────────────────────

    def _em_loop(self, pairs, max_iterations=50, epsilon=1e-8):
        if not pairs:
            return pairs
        prev_ll = -float('inf')
        for iteration in range(max_iterations):
            total = sum(p['joint_prob'] for p in pairs)
            if total == 0:
                break
            for p in pairs:
                p['posterior'] = p['joint_prob'] / total
            ll = math.log(total)

            # M-step
            sums = {}
            for p in pairs:
                post = p['posterior']
                k1 = self._hap_key(self.full_ref[p['h1_idx']])
                k2 = self._hap_key(self.full_ref[p['h2_idx']])
                sums[k1] = sums.get(k1, 0.0) + post
                sums[k2] = sums.get(k2, 0.0) + post
            tw = sum(sums.values())
            if tw > 0:
                for k in sums:
                    sums[k] /= tw
            for p in pairs:
                k1 = self._hap_key(self.full_ref[p['h1_idx']])
                k2 = self._hap_key(self.full_ref[p['h2_idx']])
                f1 = sums.get(k1, 0.0)
                f2 = sums.get(k2, 0.0)
                p['h1_freq'], p['h2_freq'] = f1, f2
                p['joint_prob'] = f1 * f1 if p['h1_idx'] == p['h2_idx'] else 2.0 * f1 * f2
            delta = ll - prev_ll
            if abs(delta) < epsilon:
                for p in pairs:
                    p['em_iterations'] = iteration + 1
                    p['log_likelihood'] = ll
                t = sum(p['joint_prob'] for p in pairs)
                if t > 0:
                    for p in pairs:
                        p['posterior'] = p['joint_prob'] / t
                print(f"    EM converged in {iteration+1} iterations")
                return pairs
            prev_ll = ll
        if pairs:
            pairs[0]['em_iterations'] = max_iterations
            pairs[0]['log_likelihood'] = prev_ll
        t = sum(p['joint_prob'] for p in pairs)
        if t > 0:
            for p in pairs:
                p['posterior'] = p['joint_prob'] / t
        print(f"    EM reached max iterations ({max_iterations})")
        return pairs


# ════════════════════════════════════════════════════════════════════
# Mini-Benchmark (50 patients with DPB1 masked)
# ════════════════════════════════════════════════════════════════════

def run_mini_benchmark(n_patients=50):
    """Run 50 patients from synthetic set with DPB1 masked."""
    import csv

    syn_path = BASE_DIR / "data" / "clean" / "synthetic_truth_set.csv"
    engine = HaploEM(population='Global')
    engine.connect()

    # Load patients
    all_patients = []
    with open(syn_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse the observed JSON
            obs_json = row.get('observed_json', '{}')
            try:
                obs = json.loads(obs_json)
            except json.JSONDecodeError:
                continue
            # Keep only patients that have DPB1 in their observed data
            # (we want to test the case where DPB1 is masked)
            dpb1_is_null = obs.get('hla_dpb1') is None
            dpa1_is_null = obs.get('hla_dpa1') is None
            if dpb1_is_null:
                all_patients.append((row, obs))
            if len(all_patients) >= 200:  # collect more, we'll filter
                break

    # Ensure we get good variety: some with DPA1 known, some without
    dpb1_only = [p for p in all_patients if p[1].get('hla_dpb1') is None]
    filtered = dpb1_only[:n_patients]

    if len(filtered) < n_patients:
        print(f"⚠️ Only found {len(filtered)} patients with DPB1 masked (wanted {n_patients})")
    else:
        filtered = filtered[:n_patients]

    print(f"\n{'='*72}")
    print(f"  MINI-BENCHMARK: {len(filtered)} patients with DPB1 masked")
    print(f"{'='*72}")

    results = []
    for idx, (csv_row, obs) in enumerate(filtered):
        pid = csv_row.get('patient_id', f"PAT_{idx}")
        true_json = csv_row.get('true_phased_json', '{}')
        try:
            truth = json.loads(true_json)
        except json.JSONDecodeError:
            continue

        # Run DP-enhanced imputation
        result = engine.progressive_em(obs)

        # Check DPB1 accuracy
        dpb1_ok = False
        dp_conf = 0.0
        top_label = ""
        if result.get('top_3') and len(result['top_3']) > 0:
            t3 = result['top_3'][0]
            top_label = f"H1={t3.get('imputed_dpb1_h1','?')} H2={t3.get('imputed_dpb1_h2','?')}"
            dp_weight = t3.get('dp_weight', 0.0)
            dp_detail = t3.get('dp_detail', '')
            true_h1 = truth.get('haplotype_1', {})
            true_h2 = truth.get('haplotype_2', {})
            t_dpb1_1 = true_h1.get('hla_dpb1', '')
            t_dpb1_2 = true_h2.get('hla_dpb1', '')
            i_dpb1_1 = t3.get('imputed_dpb1_h1', '')
            i_dpb1_2 = t3.get('imputed_dpb1_h2', '')
            # Check: do imputed DPB1 alleles match truth (order-agnostic)?
            match_1 = _allele_match(t_dpb1_1, i_dpb1_1) or _allele_match(t_dpb1_1, i_dpb1_2)
            match_2 = _allele_match(t_dpb1_2, i_dpb1_1) or _allele_match(t_dpb1_2, i_dpb1_2)
            dpb1_ok = match_1 and match_2
            dp_conf = dp_weight

        results.append({
            'patient_id': pid,
            'dpb1_accurate': dpb1_ok,
            'dp_conf': dp_conf,
            'dp_detail': t3.get('dp_detail', '') if result.get('top_3') else '',
            'top_posterior': t3.get('posterior', 0.0) if result.get('top_3') else 0.0,
            'top_label': top_label,
            'n_pairs': result.get('total_pairs_final', 0),
        })

    # ── Results summary ────────────────────────────────────────────
    total = len(results)
    accurate = sum(1 for r in results if r['dpb1_accurate'])
    not_accurate = total - accurate
    rate = accurate / total * 100 if total > 0 else 0

    # Group by whether DPA1 was known or masked
    dpa1_known = []
    dpa1_masked = []
    for r_idx, r in enumerate(results):
        pid = r['patient_id']
        # Find this patient's obs to check DPA1
        for cr, ob in filtered:
            if cr.get('patient_id') == pid:
                if ob.get('hla_dpa1') is not None:
                    dpa1_known.append(r)
                else:
                    dpa1_masked.append(r)
                break

    print(f"\n  ── DPB1 Accuracy ──")
    print(f"  Total: {total}/{total} ({rate:.1f}%)")
    print(f"  Accurate: {accurate}")
    print(f"  Inaccurate: {not_accurate}")
    print()

    if dpa1_known:
        ak = sum(1 for r in dpa1_known if r['dpb1_accurate'])
        print(f"  When DPA1 known ({len(dpa1_known)}): "
              f"{ak}/{len(dpa1_known)} ({100*ak/len(dpa1_known):.1f}%)")
    if dpa1_masked:
        am = sum(1 for r in dpa1_masked if r['dpb1_accurate'])
        print(f"  When DPA1 masked ({len(dpa1_masked)}): "
              f"{am}/{len(dpa1_masked)} ({100*am/len(dpa1_masked):.1f}%)")

    print(f"\n  ── Sample Results (first 5) ──")
    print(f"  {'Patient':15} {'DPB1 OK':>8} {'Conf':>6} {'Posterior':>10} {'DP Detail'}")
    print(f"  {'-'*60}")
    for r in results[:5]:
        print(f"  {r['patient_id']:15} {str(r['dpb1_accurate']):>8} "
              f"{r['dp_conf']:>6.3f} {r['top_posterior']:>10.6f} "
              f"{r['dp_detail'][:30]}")

    # Compare against old results from benchmark
    print(f"\n  ── Comparison vs Original Engine ──")
    print(f"  Original DPB1 accuracy (Phase 6): 11.5%")
    print(f"  DP-enhanced DPB1 accuracy:         {rate:.1f}%")
    print(f"  Improvement:                        {rate - 11.5:.1f} pp")

    engine.close()
    return results


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 72)
    print("  HaploStats — DPB1 Statistical Model (Phase 7)")
    print("=" * 72)

    results = run_mini_benchmark(n_patients=50)

    print(f"\n{'='*72}")
    print(f"  Phase 7 complete. DPB1 decoupled model ready.")
    print(f"{'='*72}")
