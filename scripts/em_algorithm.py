#!/usr/bin/env python3
"""
HaploStats — Expectation-Maximization (EM) Engine with Progressive Insertion
Phase 4: Resolve highly ambiguous / missing genotype data without O(n²) blowup.

Algorithm
---------
  1. Define locus blocks (core → fine → full resolution)
  2. For each block:
     a. E-step:    Calculate posterior for every valid haplotype pair using
                   current frequency estimates (Bayes' theorem + HWE).
     b. M-step:    Update individual haplotype frequencies by marginalising
                   pair posteriors: f(H) = ½ Σ P({H, Hⱼ} | G)
     c. Trim:      Discard pairs with posterior < threshold (default 1e-4).
     d. Insert:    Expand surviving haplotypes to the next loci and repeat.
  3. Final EM loop at full resolution → sorted results.
"""

import sqlite3
import json
import math
import sys
from pathlib import Path
from copy import deepcopy

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH  = BASE_DIR / "db" / "haplostats.db"


# ── Allele matching utilities (shared with bayesian_calc.py) ───────

def _norm(a: str) -> str:
    return a.strip().replace(' ', '')

def _allele_match(pat: str, ref: str) -> bool:
    pa = _norm(pat)
    ra = _norm(ref)
    if pa == ra:
        return True
    # JSON array on ref side
    if ra.startswith('['):
        try:
            return any(_allele_match(pa, c) for c in json.loads(ra))
        except Exception:
            return False
    # JSON array on patient side
    if pa.startswith('['):
        try:
            return any(_allele_match(c, ra) for c in json.loads(pa))
        except Exception:
            return False
    # Prefix / resolution-level matching (e.g. 2-field vs 4-field)
    pp = pa.split(':')
    rp = ra.split(':')
    for a, b in zip(pp, rp):
        if a != b:
            return False
    return True


# ── Locus block definitions for progressive insertion ──────────────

# Ordered: core Class I → Class II DR/DQ → remaining Class II
LOCUS_BLOCKS = [
    ('block_1_core',     ['hla_a', 'hla_b', 'hla_c']),
    ('block_2_drdq',     ['hla_drb1', 'hla_dqb1']),
    ('block_3_fine',     ['hla_drb345', 'hla_dqa1', 'hla_dpa1', 'hla_dpb1']),
]

ALL_LOCI = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
            'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

LOCUS_LABEL = {
    'hla_a': 'HLA-A', 'hla_c': 'HLA-C', 'hla_b': 'HLA-B',
    'hla_drb345': 'DRB345', 'hla_drb1': 'DRB1',
    'hla_dqa1': 'DQA1', 'hla_dqb1': 'DQB1',
    'hla_dpa1': 'DPA1', 'hla_dpb1': 'DPB1',
}


# ════════════════════════════════════════════════════════════════════
# HaploEM Class
# ════════════════════════════════════════════════════════════════════

class HaploEM:
    """
    Expectation-Maximisation haplotype inference with progressive
    locus-block insertion to control combinatorial explosion.
    """

    def __init__(self, db_path: str = None, population: str = 'Global'):
        self.db_path    = str(db_path or DB_PATH)
        self.population = population
        self.freq_col   = f'{population}_freq'
        self.conn       = None
        self.full_ref   = None   # list[dict] — all 577 haplotypes
        self.current_freqs = None  # dict[tuple] → float during EM

    # ── public API ──────────────────────────────────────────────────

    def connect(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._load_reference()
        return self

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    def progressive_em(self, patient_genotype: dict,
                       trim_threshold: float = 1e-4,
                       max_em_iterations: int = 50) -> dict:
        """
        Run EM with progressive locus insertion.

        Parameters
        ----------
        patient_genotype : dict
            Locus → [allele1, allele2] or None for missing.
            e.g. {'hla_a': ['A*02:01', 'A*01:01'],
                  'hla_c': None, ...}

        Returns
        -------
        dict with keys: population, convergence, blocks, final_pairs, ...
        """
        # Normalise input
        gt = {}
        for loc in ALL_LOCI:
            val = patient_genotype.get(loc)
            if val is None or val == '':
                gt[loc] = None
            elif isinstance(val, str):
                gt[loc] = [val]
            elif isinstance(val, list) and len(val) == 0:
                gt[loc] = None
            else:
                gt[loc] = list(val)

        print("\n" + "=" * 72)
        print("  HaploEM — Progressive Insertion Engine")
        print("=" * 72)
        print(f"  Population:         {self.population}")
        print(f"  Freq column:        {self.freq_col}")
        known = sum(1 for v in gt.values() if v is not None)
        print(f"  Patient typed loci: {known} / 9")
        for loc in ALL_LOCI:
            v = gt.get(loc)
            if v:
                print(f"    {LOCUS_LABEL[loc]:10} {v[0]:30} / {v[1] if len(v)>1 else v[0]}")
            else:
                print(f"    {LOCUS_LABEL[loc]:10} <missing>")

        # Initial frequency set = reference frequencies
        self.current_freqs = {}
        for h in self.full_ref:
            key = self._hap_key(h)
            self.current_freqs[key] = h.get(self.freq_col, 0.0)

        # Track surviving index set across blocks
        # We work with indices into full_ref
        surviving_indices = set(range(len(self.full_ref)))

        block_log = []
        total_pairs_before_trim = 0

        for block_name, block_loci in LOCUS_BLOCKS:
            print(f"\n  ── {block_name.upper()} ({', '.join(block_loci)}) ──")

            # Insert: filter surviving haplotypes to those matching
            # known loci in this block
            known_in_block = [l for l in block_loci if gt.get(l) is not None]
            missing_in_block = [l for l in block_loci if gt.get(l) is None]

            new_surviving = set()
            for idx in surviving_indices:
                h = self.full_ref[idx]
                ok = True
                for loc in known_in_block:
                    pat_alleles = gt[loc]
                    hap_allele = h.get(loc, '')
                    if not hap_allele:
                        ok = False
                        break
                    # Haplotype must carry at least one of the patient's alleles
                    if not any(_allele_match(a, hap_allele) for a in pat_alleles):
                        ok = False
                        break
                if ok:
                    new_surviving.add(idx)

            surviving_indices = new_surviving
            print(f"    Haplotypes matching: {len(surviving_indices)}")

            if len(surviving_indices) < 2:
                print("    ⚠️  Too few haplotypes — cannot form pairs.")
                break

            # Build all valid pairs from survivors
            # Only check loci constrained SO FAR (all blocks up to this one)
            constrained_so_far = []
            for bn, bl in LOCUS_BLOCKS:
                constrained_so_far.extend(bl)
                if bn == block_name:
                    break

            pairs = self._build_pairs(list(surviving_indices), gt, constrained_so_far)
            total_pairs_before_trim += len(pairs)
            print(f"    Pairs formed:       {len(pairs)}")

            if len(pairs) < 1:
                print("    ⚠️  Zero valid pairs — aborting.")
                break

            # ── EM loop ─────────────────────────────────────────────
            pairs = self._em_loop(pairs, max_iterations=max_em_iterations)
            final_ll = pairs[0].get('log_likelihood', pairs[-1].get('log_likelihood', 0))

            # ── Trim ────────────────────────────────────────────────
            pre_trim = len(pairs)
            pairs = [p for p in pairs if p['posterior'] >= trim_threshold]
            post_trim = len(pairs)
            print(f"    EM converged. Pairs before trim: {pre_trim}, after: {post_trim}")

            # Update surviving indices to those in surviving pairs
            used_indices = set()
            for p in pairs:
                used_indices.add(p['h1_idx'])
                used_indices.add(p['h2_idx'])
            surviving_indices = used_indices

            block_log.append({
                'block': block_name,
                'haplotypes': len(surviving_indices),
                'pairs_before_em': pre_trim,
                'pairs_after_trim': post_trim,
                'converged_iterations': pairs[0].get('em_iterations', 0) if pairs else 0,
            })

        # ── Final sort and build response ──────────────────────────
        pairs.sort(key=lambda p: p['posterior'], reverse=True)

        # Normalise posteriors to sum to 1
        total_post = sum(p['posterior'] for p in pairs)
        if total_post > 0:
            for p in pairs:
                p['posterior'] /= total_post

        # Calculate entropy
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
            })

        print(f"\n  {'='*68}")
        print(f"  FINAL: {len(pairs)} pairs after progressive insertion + EM")
        print(f"  Entropy: {entropy:.4f} bits")

        return {
            'population': self.population,
            'patient_genotype': {k: v for k, v in gt.items()},
            'blocks': block_log,
            'total_pairs_before_trim': total_pairs_before_trim,
            'total_pairs_final': len(pairs),
            'entropy': round(entropy, 4),
            'top_3': top3,
            'pairs': [
                {
                    'rank': i+1,
                    'posterior': round(p['posterior'], 8),
                    'cumulative': round(sum(pairs[j]['posterior'] for j in range(i+1)), 8),
                    'h1_frequency': round(p['h1_freq'], 6),
                    'h2_frequency': round(p['h2_freq'], 6),
                    'k_factor': p['k_factor'],
                }
                for i, p in enumerate(pairs[:10])
            ],
        }

    # ── Internal helpers ────────────────────────────────────────────

    def _load_reference(self):
        cur = self.conn.cursor()
        cols = ', '.join(ALL_LOCI + [self.freq_col])
        rows = cur.execute(
            f"SELECT {cols} FROM haplotypes "
            f"WHERE {self.freq_col} IS NOT NULL "
            f"ORDER BY {self.freq_col} DESC"
        ).fetchall()
        self.full_ref = [dict(r) for r in rows]
        print(f"  [HaploEM] Loaded {len(self.full_ref)} reference haplotypes "
              f"(pop={self.population})")

    @staticmethod
    def _hap_key(h: dict) -> tuple:
        """Immutable tuple key for a full 9-locus haplotype."""
        return tuple(h.get(l, '') for l in ALL_LOCI)

    @staticmethod
    def _hap_key_short(h: dict, loci: list) -> tuple:
        """Immutable tuple key for a subset of loci."""
        return tuple(h.get(l, '') for l in loci)

    def _make_label(self, idx: int) -> str:
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

    def _build_pairs(self, indices: list, gt: dict,
                     constrained_loci: list) -> list:
        """
        Build all valid haplotype pairs, checking only constrained loci.
        Pairs where H1 == H2 are allowed (homozygous).
        """
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
                        continue  # missing → any allele works
                    a1, a2 = (pat[0], pat[1]) if len(pat) >= 2 else (pat[0], pat[0])
                    h1v = h1.get(loc, '')
                    h2v = h2.get(loc, '')
                    # Both patient alleles must be covered by H1 ∪ H2
                    c1 = _allele_match(a1, h1v) or _allele_match(a1, h2v)
                    c2 = _allele_match(a2, h1v) or _allele_match(a2, h2v)
                    if not (c1 and c2):
                        valid = False
                        break
                if not valid:
                    continue

                f1 = h1.get(self.freq_col, 0.0)
                f2 = h2.get(self.freq_col, 0.0)
                if f1 <= 0.0 or f2 <= 0.0:
                    continue

                if i == j:
                    joint_prob = f1 * f1
                    k_factor = 1.0
                else:
                    joint_prob = 2.0 * f1 * f2
                    k_factor = 2.0

                pairs.append({
                    'h1_idx': indices[i],
                    'h2_idx': indices[j],
                    'h1_freq': f1,
                    'h2_freq': f2,
                    'joint_prob': joint_prob,
                    'k_factor': k_factor,
                    'posterior': 0.0,  # set in EM
                })

        return pairs

    def _em_loop(self, pairs: list, max_iterations: int = 50,
                  epsilon: float = 1e-8) -> list:
        """
        Run EM on a fixed set of candidate pairs.

        E-step:  P(H_pair | G) = joint_prob / sum(all_joint_probs)
        M-step:  f(H) = (1/2) Σ_{pairs containing H} P(H, Hⱼ | G)
                 Then re-compute joint_probs using updated frequencies.
        """
        if not pairs:
            return pairs

        prev_ll = -float('inf')

        for iteration in range(max_iterations):
            # ── E-step: normalise posteriors ────────────────────────
            total_joint = sum(p['joint_prob'] for p in pairs)

            if total_joint == 0.0:
                break

            for p in pairs:
                p['posterior'] = p['joint_prob'] / total_joint

            # Compute log-likelihood for convergence check
            log_likelihood = math.log(total_joint)

            # ── M-step: update haplotype frequencies ────────────────
            # Accumulate weighted sums over pairs
            freq_sums = {}
            for p in pairs:
                post = p['posterior']
                h1_key = self._hap_key(self.full_ref[p['h1_idx']])
                h2_key = self._hap_key(self.full_ref[p['h2_idx']])
                freq_sums[h1_key] = freq_sums.get(h1_key, 0.0) + post
                freq_sums[h2_key] = freq_sums.get(h2_key, 0.0) + post

            # Normalise: each patient contributes 2 haplotype copies
            total_weight = sum(freq_sums.values())
            if total_weight > 0:
                for key in freq_sums:
                    freq_sums[key] /= total_weight

            # ── Re-compute joint probabilities with new freqs ───────
            for p in pairs:
                h1_key = self._hap_key(self.full_ref[p['h1_idx']])
                h2_key = self._hap_key(self.full_ref[p['h2_idx']])
                f1 = freq_sums.get(h1_key, 0.0)
                f2 = freq_sums.get(h2_key, 0.0)
                p['h1_freq'] = f1
                p['h2_freq'] = f2
                if p['h1_idx'] == p['h2_idx']:
                    p['joint_prob'] = f1 * f1
                else:
                    p['joint_prob'] = 2.0 * f1 * f2

            # ── Check convergence ───────────────────────────────────
            delta = log_likelihood - prev_ll
            if abs(delta) < epsilon:
                for p in pairs:
                    p['em_iterations'] = iteration + 1
                    p['log_likelihood'] = log_likelihood
                # Final E-step
                total = sum(p['joint_prob'] for p in pairs)
                if total > 0:
                    for p in pairs:
                        p['posterior'] = p['joint_prob'] / total
                print(f"    EM converged in {iteration+1} iterations "
                      f"(ΔLL={delta:.2e})")
                return pairs

            prev_ll = log_likelihood

        # Hit max iterations
        if pairs:
            pairs[0]['em_iterations'] = max_iterations
            pairs[0]['log_likelihood'] = prev_ll
        total = sum(p['joint_prob'] for p in pairs)
        if total > 0:
            for p in pairs:
                p['posterior'] = p['joint_prob'] / total
        print(f"    EM reached max iterations ({max_iterations})")
        return pairs


# ════════════════════════════════════════════════════════════════════
# Mock Patient: highly ambiguous, only typed at A, B, DRB1
# ════════════════════════════════════════════════════════════════════

def build_ambiguous_mock():
    """
    Create a mock patient typed ONLY at HLA-A, HLA-B, and HLA-DRB1.
    All other loci (C, DRB345, DQA1, DQB1, DPA1, DPB1) are missing.
    """
    return {
        'hla_a':      ['A*02:01:01:01', 'A*01:01:01:01'],
        'hla_c':      None,   # missing
        'hla_b':      ['B*44:02:01:01', 'B*08:01:01:01'],
        'hla_drb345': None,   # missing
        'hla_drb1':   ['DRB1*04:01:01:01SG', 'DRB1*03:01:01:01SG'],
        'hla_dqa1':   None,   # missing
        'hla_dqb1':   None,   # missing
        'hla_dpa1':   None,   # missing
        'hla_dpb1':   None,   # missing
    }


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import pprint

    patient = build_ambiguous_mock()

    engine = HaploEM(population='Global')
    engine.connect()

    result = engine.progressive_em(patient, trim_threshold=1e-4)

    print("\n\n📊 PROGRESSIVE INSERTION RESULTS")
    print("=" * 72)
    print(f"  Population:                     {result['population']}")
    print(f"  Total pairs before trim:        {result['total_pairs_before_trim']}")
    print(f"  Total pairs after progressive:  {result['total_pairs_final']}")
    print(f"  Entropy:                        {result['entropy']} bits")
    print()

    print("  ── Block Progression ──")
    for b in result['blocks']:
        print(f"  {b['block']:20}  haplotypes={b['haplotypes']:4}  "
              f"pairs before EM={b['pairs_before_em']:5}  "
              f"after trim={b['pairs_after_trim']:4}")

    print()
    print("  ── Top 3 Inferred 9-Locus Haplotype Pairs ──")
    print()
    for t3 in result['top_3']:
        print(f"  Rank {t3['rank']}  (posterior = {t3['posterior']:.6f}, "
              f"cumulative = {t3['cumulative']:.6f})")
        print(f"    H1: {t3['haplotype_1']}")
        print(f"    H2: {t3['haplotype_2']}")
        print(f"    H1 freq={t3['h1_frequency']:.6f}, H2 freq={t3['h2_frequency']:.6f}")
        print()

    print("\n  ── Top 10 Pair Summary (JSON payload) ──")
    payload = {
        'population': result['population'],
        'blocks': result['blocks'],
        'total_pairs_before_trim': result['total_pairs_before_trim'],
        'total_pairs_final': result['total_pairs_final'],
        'entropy': result['entropy'],
        'top_3': result['top_3'],
    }
    print(pprint.pformat(payload, indent=2, width=72))

    engine.close()
    print("\n✅ EM engine complete.")
