#!/usr/bin/env python3
"""
HaploStats — Bayesian Inference Module
Phase 3: Posterior probability calculation for unphased patient genotypes.

Given an unphased 9-locus HLA genotype, find all possible phased haplotype
pairs from the population reference and rank them by posterior probability.

Bayes: P(H_pair | G) ∝ P(G | H_pair) × P(H_pair)
  - P(H_pair) = k × P(H1) × P(H2), where k = 2 for heterozygotes, 1 for homozygotes
  - P(G | H_pair) = 1 if H1+H2 produce G, 0 otherwise
"""

import sqlite3
import json
import sys
from pathlib import Path
from collections import defaultdict
import math

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH = BASE_DIR / "db" / "haplostats.db"


# ────────────────────────────────────────────────────────────────────
# Utility: Allele Matching
# ────────────────────────────────────────────────────────────────────

def _normalise(allele: str) -> str:
    """Strip whitespace and standardise formatting."""
    return allele.strip().replace(' ', '')


def _allele_matches(patient_allele: str, ref_allele: str) -> bool:
    """
    Check if a patient allele matches a reference allele.
    Handles exact matches, JSON-encoded ambiguous arrays, and
    resolution-level prefix matching (e.g., 2-field vs 4-field).
    """
    pa = _normalise(patient_allele)
    ra = _normalise(ref_allele)

    # 1. Exact string match
    if pa == ra:
        return True

    # 2. Reference is a JSON array of possible alleles
    if ra.startswith('['):
        try:
            candidates = json.loads(ra)
        except (json.JSONDecodeError, TypeError):
            return False
        for candidate in candidates:
            if _allele_matches(pa, candidate):
                return True
        return False

    # 3. Patient is a JSON array (unlikely but handle symmetrically)
    if pa.startswith('['):
        try:
            candidates = json.loads(pa)
        except (json.JSONDecodeError, TypeError):
            return False
        for candidate in candidates:
            if _allele_matches(candidate, ra):
                return True
        return False

    # 4. Prefix / resolution-level matching
    #    e.g. patient says A*01:01, reference has A*01:01:01:01
    #    Match any level if one is a prefix of the other (separator-aware)
    #    Allele format: LOCUS*NN:NN[:NN[:NN]][SUFFIX]
    parts_p = pa.split(':')
    parts_r = ra.split(':')
    # Match up to the shorter length
    min_len = min(len(parts_p), len(parts_r))
    for i in range(min_len):
        if parts_p[i] != parts_r[i]:
            return False
    return True


def _haplotype_covers_genotype(haplotype: dict, patient_allele1: str,
                                patient_allele2: str, locus_col: str) -> bool:
    """Check if a single haplotype carries one of the two patient alleles
    at the given locus."""
    h_allele = haplotype.get(locus_col, '')
    if not h_allele or h_allele == '':
        return False
    return (_allele_matches(patient_allele1, h_allele) or
            _allele_matches(patient_allele2, h_allele))


# ────────────────────────────────────────────────────────────────────
# HaploMath Engine
# ────────────────────────────────────────────────────────────────────

class HaploMath:
    """
    Bayesian haplotype inference engine.

    Connects to haplostats.db and provides posterior probability
    calculations for unphased patient genotypes.
    """

    LOCI = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
            'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

    LOCUS_LABELS = ['HLA-A', 'HLA-C', 'HLA-B', 'HLA-DRB345', 'HLA-DRB1',
                    'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1']

    def __init__(self, db_path: str = None, population: str = 'Global'):
        self.db_path = str(db_path or DB_PATH)
        self.population = population
        self.freq_col = f'{population}_freq'
        self.conn = None
        self._all_haplotypes = None  # cached list

    def connect(self):
        """Open database connection and cache all haplotypes."""
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._load_haplotypes()
        return self

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    def _load_haplotypes(self):
        """Cache all haplotypes from the reference database."""
        cur = self.conn.cursor()
        cols = ', '.join(self.LOCI + [self.freq_col])
        rows = cur.execute(
            f"SELECT {cols} FROM haplotypes "
            f"WHERE {self.freq_col} IS NOT NULL "
            f"ORDER BY {self.freq_col} DESC"
        ).fetchall()
        self._all_haplotypes = [dict(r) for r in rows]
        print(f"  [HaploMath] Loaded {len(self._all_haplotypes)} haplotypes "
              f"(population={self.population}, freq_col={self.freq_col})")

    def calculate_posterior(self, genotype: dict) -> dict:
        """
        Main entry point.

        Parameters
        ----------
        genotype : dict
            Maps locus_col -> [allele1, allele2]
            e.g. {'hla_a': ['A*01:01', 'A*03:01'],
                  'hla_b': ['B*08:01', 'B*07:02'], ...}
            Missing loci are treated as 'any match' (not constrained).

        Returns
        -------
        dict with keys:
            'patient_genotype'   : input genotype (normalised)
            'population'         : population used
            'total_possible_pairs': int
            'pairs'              : list of dicts sorted by posterior desc
            'top_pair'           : dict (highest posterior)
            'entropy'            : float (uncertainty measure)
        """
        # ── Step 1: Filter haplotypes to those matching at least one
        #            patient allele at every constrained locus ────────
        constrained_loci = [loc for loc in self.LOCI if loc in genotype]

        candidates = []
        for hap in self._all_haplotypes:
            valid = True
            for loc in constrained_loci:
                a1, a2 = genotype[loc]
                if not _haplotype_covers_genotype(hap, a1, a2, loc):
                    valid = False
                    break
            if valid:
                candidates.append(hap)

        # ── Step 2: Build all valid haplotype pairs ─────────────────
        pairs = []
        n_cand = len(candidates)

        for i in range(n_cand):
            h1 = candidates[i]
            for j in range(i, n_cand):
                h2 = candidates[j]
                # Check that together (H1, H2) cover both patient alleles
                # at every constrained locus
                valid_pair = True
                for loc in constrained_loci:
                    a1, a2 = genotype[loc]
                    h1_val = h1.get(loc, '')
                    h2_val = h2.get(loc, '')
                    set_h = set()
                    for v in (h1_val, h2_val):
                        if v.startswith('['):
                            try:
                                for x in json.loads(v):
                                    set_h.add(x)
                            except Exception:
                                set_h.add(v)
                        else:
                            set_h.add(v)
                    # Both patient alleles must be covered
                    if not (_allele_matches(a1, h1_val) or _allele_matches(a1, h2_val)):
                        valid_pair = False
                        break
                    if not (_allele_matches(a2, h1_val) or _allele_matches(a2, h2_val)):
                        valid_pair = False
                        break
                if not valid_pair:
                    continue

                # ── Step 3: Calculate joint probability ─────────────
                f1 = h1.get(self.freq_col) or 0.0
                f2 = h2.get(self.freq_col) or 0.0
                if f1 <= 0.0 or f2 <= 0.0:
                    continue

                if i == j:
                    # Homozygous pair: P(H, H) = P(H)²
                    joint_prob = f1 * f1
                    k_factor = 1.0
                else:
                    # Heterozygous pair: P(H1, H2) = 2 × P(H1) × P(H2)
                    joint_prob = 2.0 * f1 * f2
                    k_factor = 2.0

                h1_label = self._make_haplotype_label(h1)
                h2_label = self._make_haplotype_label(h2)

                pairs.append({
                    'haplotype_1': h1_label,
                    'haplotype_2': h2_label,
                    'h1_freq': f1,
                    'h2_freq': f2,
                    'joint_prob': joint_prob,
                    'k_factor': k_factor,
                    'h1_id': i,
                    'h2_id': j,
                })

        # ── Step 4: Normalise to posterior probabilities ────────────
        total_prob = sum(p['joint_prob'] for p in pairs)

        if total_prob == 0.0:
            return {
                'patient_genotype': genotype,
                'population': self.population,
                'total_possible_pairs': 0,
                'pairs': [],
                'top_pair': None,
                'entropy': 0.0,
                'error': 'No matching haplotype pairs found.'
            }

        for p in pairs:
            p['posterior'] = p['joint_prob'] / total_prob

        # Sort descending by posterior
        pairs.sort(key=lambda x: x['posterior'], reverse=True)

        # ── Step 5: Rank and build response ─────────────────────────
        ranked_pairs = []
        cumulative = 0.0
        for rank_idx, p in enumerate(pairs, start=1):
            cumulative += p['posterior']
            ranked_pairs.append({
                'rank': rank_idx,
                'haplotype_1': p['haplotype_1'],
                'haplotype_2': p['haplotype_2'],
                'h1_frequency': round(p['h1_freq'], 6),
                'h2_frequency': round(p['h2_freq'], 6),
                'k_factor': p['k_factor'],
                'posterior': round(p['posterior'], 8),
                'cumulative': round(cumulative, 8),
            })

        # ── Shannon entropy of the posterior distribution ───────────
        entropy = 0.0
        for p in pairs:
            if p['posterior'] > 0:
                entropy -= p['posterior'] * math.log2(p['posterior'])

        top = ranked_pairs[0] if ranked_pairs else None

        return {
            'patient_genotype': {k: list(v) for k, v in genotype.items()},
            'population': self.population,
            'total_candidate_haplotypes': n_cand,
            'total_possible_pairs': len(ranked_pairs),
            'pairs': ranked_pairs,
            'top_pair': top,
            'entropy': round(entropy, 4),
        }

    @staticmethod
    def _make_haplotype_label(hap: dict) -> str:
        """Build a compact human-readable haplotype label."""
        parts = []
        labels = ['HLA-A', 'HLA-C', 'HLA-B', 'DRB345', 'DRB1',
                  'DQA1', 'DQB1', 'DPA1', 'DPB1']
        locus_keys = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
                      'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']
        for lbl, lk in zip(labels, locus_keys):
            val = hap.get(lk, '')
            if val:
                # Truncate JSON arrays for readability
                if val.startswith('['):
                    try:
                        arr = json.loads(val)
                        if len(arr) <= 2:
                            val = '/'.join(arr)
                        else:
                            val = f"{arr[0]}/...({len(arr)} alleles)"
                    except Exception:
                        pass
                parts.append(f"{lbl}={val}")
            else:
                parts.append(f"{lbl}=?")
        return ' | '.join(parts)


# ────────────────────────────────────────────────────────────────────
# Mock Patient Genotype (Hardcoded Test)
# ────────────────────────────────────────────────────────────────────

def build_mock_patient():
    """
    Create a hardcoded unphased 9-locus patient genotype.

    This patient has the two most common European haplotypes:
      H1: A*01:01:01:01 C*07:01:01:01 B*08:01:01:01 DRB3*01:01:02:01
          DRB1*03:01:01:01SG DQA1*05:01:01:02 DQB1*02:01:01
          DPA1*01:03:01:02 DPB1*04:01:01:01/DPB1*04:01:01:02  (f=0.021118)

      H2: A*03:01:01:01 C*07:02:01:03 B*07:02:01 DRB5*01:01:01
          DRB1*15:01:01:01SG DQA1*01:02:01:01SG DQB1*06:02:01
          DPA1*01:03:01:02 DPB1*04:01:01:01/DPB1*04:01:01:02  (f=0.017016)

    Each locus shows both patient alleles (unphased).
    """
    genotype = {
        'hla_a':      ['A*01:01:01:01', 'A*03:01:01:01'],
        'hla_c':      ['C*07:01:01:01', 'C*07:02:01:03'],
        'hla_b':      ['B*08:01:01:01', 'B*07:02:01'],
        'hla_drb345': ['DRB3*01:01:02:01', 'DRB5*01:01:01'],
        'hla_drb1':   ['DRB1*03:01:01:01SG', 'DRB1*15:01:01:01SG'],
        'hla_dqa1':   ['DQA1*05:01:01:02', 'DQA1*01:02:01:01SG'],
        'hla_dqb1':   ['DQB1*02:01:01', 'DQB1*06:02:01'],
        'hla_dpa1':   ['DPA1*01:03:01:02', 'DPA1*01:03:01:02'],  # homozygous
        'hla_dpb1':   ['DPB1*04:01:01:01', 'DPB1*04:01:01:01'],  # homozygous
    }
    return genotype


# ────────────────────────────────────────────────────────────────────
# Main: Run Bayesian Calculation on Mock Patient
# ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import pprint

    print("=" * 72)
    print("  HaploStats — Bayesian Inference Engine (Phase 3)")
    print("=" * 72)

    # Initialise engine
    engine = HaploMath(population='Global')
    engine.connect()

    # Build mock patient
    patient = build_mock_patient()
    print("\n📋 Mock Patient Genotype (Unphased):")
    print("-" * 72)
    for loc, alleles in patient.items():
        if alleles[0] == alleles[1]:
            print(f"  {loc:15} {alleles[0]:30} (homozygous)")
        else:
            print(f"  {loc:15} {alleles[0]:30} / {alleles[1]}")
    print()

    # Run Bayesian inference
    result = engine.calculate_posterior(patient)

    # ── Print Results ───────────────────────────────────────────────
    print(f"\n📊 Bayesian Posterior Results")
    print(f"   Population:          {result['population']}")
    print(f"   Total candidate haplotypes: {result['total_candidate_haplotypes']}")
    print(f"   Total possible pairs:       {result['total_possible_pairs']}")
    print(f"   Entropy (uncertainty):      {result['entropy']} bits")
    print()

    if result['pairs']:
        print(f"{'Rank':>5} {'Posterior':>12} {'Cumulative':>12} {'k':>3}  Haplotype Pair")
        print(f"   {'-'*68}")
        for p in result['pairs'][:10]:  # Top 10
            rank = p['rank']
            post = p['posterior']
            cum = p['cumulative']
            k = p['k_factor']
            h1 = p['haplotype_1']
            h2 = p['haplotype_2']
            # Truncate long labels for display
            if len(h1) > 60:
                h1 = h1[:57] + '...'
            if len(h2) > 60:
                h2 = h2[:57] + '...'
            print(f"  {rank:>4}  {post:>10.6f}  {cum:>10.6f}  {k:>1.0f}   H1: {h1}")
            print(f"  {'':>5} {'':>12} {'':>12} {'':>3}   H2: {h2}")
            print()
    else:
        print("  ❌ No matching haplotype pairs found.")
        print("  Check that patient alleles exist in the reference database.")

    # Print full JSON payload for verification
    print("\n📦 Full JSON Payload (Top 5 pairs):")
    print("-" * 72)
    payload = {
        'population': result['population'],
        'total_possible_pairs': result['total_possible_pairs'],
        'entropy': result['entropy'],
        'pairs': [
            {k: p[k] for k in ['rank', 'posterior', 'cumulative',
                               'h1_frequency', 'h2_frequency', 'k_factor']}
            for p in result['pairs'][:5]
        ],
        'top_pair': {
            'haplotype_1': result['pairs'][0]['haplotype_1'][:80] + '...',
            'haplotype_2': result['pairs'][0]['haplotype_2'][:80] + '...',
            'posterior': result['pairs'][0]['posterior'],
        } if result['pairs'] else None,
    }
    print(pprint.pformat(payload, indent=2, width=72))

    engine.close()
    print("\n✅ Bayesian calculation complete.")
