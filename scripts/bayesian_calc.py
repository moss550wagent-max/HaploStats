#!/usr/bin/env python3
"""
HaploStats — Bayesian Inference Module (Multi-Population)
Phase 6+: Posterior probability calculation with all-population diplotype
frequencies (Hardy-Weinberg 2pq / p²).

Given an unphased 9-locus HLA genotype, find all possible phased haplotype
pairs from the population reference and rank them by posterior probability.

For every matched pair, this module returns the diplotype frequency for ALL
populations in the reference database, not just one selected population.

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
import pprint

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH = BASE_DIR / "db" / "haplostats.db"

# ── Multi-Population Frequency Column Map ──────────────────────────
# These are the population frequency columns available in the reference DB.
# Each maps a human-readable population code to its SQL column.
# Uses standard transplant-registry codes where applicable.

POPULATION_COLUMNS = {
    'Global':   'Global_freq',
    'AFA':      'AfAm_freq',       # African American
    'API':      'API_freq',        # Asian / Pacific Islander
    'CAU':      'EuAm_freq',       # European American (Caucasian)
    'HIS':      'USA_Hispanic_freq',  # US Hispanic
    'NAM':      None,              # Native American — not in current reference
    'European': 'European_freq',
    'Spanish':  'Spanish_freq',
    'Mexican':  'Mexican_freq',
    'Arab':     'Arab_freq',
}

# Shorthand list for iteration and serialisation order
POPULATION_ORDER = ['Global', 'AFA', 'API', 'CAU', 'HIS', 'NAM',
                    'European', 'Spanish', 'Mexican', 'Arab']


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
    parts_p = pa.split(':')
    parts_r = ra.split(':')
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
# Hardy-Weinberg Diplotype Frequency
# ────────────────────────────────────────────────────────────────────

def _diplotype_frequency(h1_freq: float, h2_freq: float, is_homozygous: bool) -> float:
    """
    Compute the diplotype frequency under Hardy-Weinberg equilibrium.

    - Heterozygous pair: 2 × p × q  (k = 2)
    - Homozygous pair:   p²          (k = 1)
    """
    if h1_freq <= 0.0 or h2_freq <= 0.0:
        return 0.0
    if is_homozygous:
        return h1_freq * h1_freq
    return 2.0 * h1_freq * h2_freq


# ────────────────────────────────────────────────────────────────────
# HaploMath Engine
# ────────────────────────────────────────────────────────────────────

class HaploMath:
    """
    Bayesian haplotype inference engine (multi-population).

    Connects to haplostats.db and provides posterior probability
    calculations for unphased patient genotypes.  Unlike the original
    single-population engine, this version loads ALL population frequency
    columns and returns per-population diplotype frequencies (2pq / p²)
    for every matched haplotype pair.
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
        self._all_haplotypes = None  # cached list (all freq cols loaded)

    def connect(self):
        """Open database connection and cache all haplotypes
        with all population frequency columns."""
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._load_haplotypes()
        return self

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    def _load_haplotypes(self):
        """
        Cache all haplotypes from the reference database, loading EVERY
        population frequency column so we can compute per-population
        diplotype frequencies later.
        """
        cur = self.conn.cursor()

        # Build column list: loci + ALL frequency columns
        freq_cols = [col for col in POPULATION_COLUMNS.values() if col is not None]
        all_cols = ', '.join(self.LOCI + freq_cols)

        rows = cur.execute(
            f"SELECT {all_cols} FROM haplotypes "
            f"WHERE {self.freq_col} IS NOT NULL "
            f"ORDER BY {self.freq_col} DESC"
        ).fetchall()

        self._all_haplotypes = [dict(r) for r in rows]
        print(f"  [HaploMath] Loaded {len(self._all_haplotypes)} haplotypes "
              f"(all populations, {len(freq_cols)} freq columns)")

    # ── Public API ─────────────────────────────────────────────────

    def calculate_posterior(self, genotype: dict) -> dict:
        """
        Main entry point — multi-population version.

        Parameters
        ----------
        genotype : dict
            Maps locus_col -> [allele1, allele2]
            Missing loci are treated as unconstrained (any match).

        Returns
        -------
        dict with keys:
            'patient_genotype'    : input genotype (normalised)
            'population'          : selected population for ranking
            'total_possible_pairs': int
            'pairs'               : list of dicts sorted by primary-pop posterior
            'top_pair'            : dict (highest posterior)
            'entropy'             : float (Shannon entropy of posterior)
            'populations_available' : list of population codes returned
        """
        constrained_loci = [loc for loc in self.LOCI if loc in genotype]

        # ── Step 1: Filter haplotypes to those matching at least one
        #            patient allele at every constrained locus ────────
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
        pairs_raw = []
        n_cand = len(candidates)

        for i in range(n_cand):
            h1 = candidates[i]
            for j in range(i, n_cand):
                h2 = candidates[j]

                # Verify that (H1, H2) together cover both alleles at
                # every constrained locus
                valid_pair = True
                for loc in constrained_loci:
                    a1, a2 = genotype[loc]
                    h1_val = h1.get(loc, '')
                    h2_val = h2.get(loc, '')

                    # Collect all allele strings from both haplotypes
                    set_h = set()
                    for val in (h1_val, h2_val):
                        if val.startswith('['):
                            try:
                                for x in json.loads(val):
                                    set_h.add(_normalise(x))
                            except Exception:
                                set_h.add(_normalise(val))
                        else:
                            set_h.add(_normalise(val))

                    # Both patient alleles must appear in the combined set
                    if not (_allele_matches(a1, h1_val) or _allele_matches(a1, h2_val)):
                        valid_pair = False
                        break
                    if not (_allele_matches(a2, h1_val) or _allele_matches(a2, h2_val)):
                        valid_pair = False
                        break

                if not valid_pair:
                    continue

                # ── Step 3: Compute per-population diplotype frequencies ──
                is_homozygous = (i == j)
                pop_freqs = {}

                for pop_name, col_name in POPULATION_COLUMNS.items():
                    if col_name is None:
                        pop_freqs[pop_name] = 0.0
                    else:
                        f1 = h1.get(col_name) or 0.0
                        f2 = h2.get(col_name) or 0.0
                        freq = _diplotype_frequency(f1, f2, is_homozygous)
                        pop_freqs[pop_name] = round(freq, 8)

                # Use the primary population (self.population) for ranking
                # Fall back to Global if selected population has no column
                raw_col = POPULATION_COLUMNS.get(self.population)
                primary_col = raw_col if raw_col else 'Global_freq'
                f1_primary = h1.get(primary_col) or 0.0
                f2_primary = h2.get(primary_col) or 0.0
                primary_joint = _diplotype_frequency(f1_primary, f2_primary, is_homozygous)

                if primary_joint <= 0.0:
                    continue

                h1_label = self._make_haplotype_label(h1)
                h2_label = self._make_haplotype_label(h2)

                pairs_raw.append({
                    'haplotype_1': h1_label,
                    'haplotype_2': h2_label,
                    'joint_prob': primary_joint,
                    'is_homozygous': is_homozygous,
                    'population_frequencies': pop_freqs,
                    'h1_id': i,
                    'h2_id': j,
                })

        # ── Step 4: Normalise → posterior probabilities (primary pop) ──
        total_primary = sum(p['joint_prob'] for p in pairs_raw)

        if total_primary == 0.0:
            return {
                'patient_genotype': genotype,
                'population': self.population,
                'total_possible_pairs': 0,
                'pairs': [],
                'top_pair': None,
                'entropy': 0.0,
                'populations_available': POPULATION_ORDER,
                'error': 'No matching haplotype pairs found in any population.'
            }

        for p in pairs_raw:
            p['posterior'] = round(p['joint_prob'] / total_primary, 8)

        # Sort descending by primary-pop posterior
        pairs_raw.sort(key=lambda x: x['posterior'], reverse=True)

        # ── Step 5: Rank and build response ─────────────────────────
        ranked_pairs = []
        cumulative = 0.0

        for rank_idx, p in enumerate(pairs_raw, start=1):
            cumulative += p['posterior']

            ranked_pairs.append({
                'rank': rank_idx,
                'haplotype_1': p['haplotype_1'],
                'haplotype_2': p['haplotype_2'],
                'posterior': p['posterior'],
                'cumulative': round(cumulative, 8),
                'is_homozygous': p['is_homozygous'],
                'population_frequencies': p['population_frequencies'],
            })

        # ── Shannon entropy of the primary-population posterior ──────
        entropy = 0.0
        for p in pairs_raw:
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
            'populations_available': POPULATION_ORDER,
        }

    # ── Helpers ─────────────────────────────────────────────────────

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
      H1: A*01:01 C*07:01 B*08:01 DRB3*01:01 DRB1*03:01 DQA1*05:01
          DQB1*02:01 DPA1*01:03 DPB1*04:01
      H2: A*03:01 C*07:02 B*07:02 DRB5*01:01 DRB1*15:01 DQA1*01:02
          DQB1*06:02 DPA1*01:03 DPB1*04:01
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
    print("=" * 72)
    print("  HaploStats — Bayesian Inference Engine (Multi-Population)")
    print("=" * 72)

    engine = HaploMath(population='Global')
    engine.connect()

    patient = build_mock_patient()
    print("\n📋 Mock Patient Genotype (Unphased):")
    print("-" * 72)
    for loc, alleles in patient.items():
        if alleles[0] == alleles[1]:
            print(f"  {loc:15} {alleles[0]:30} (homozygous)")
        else:
            print(f"  {loc:15} {alleles[0]:30} / {alleles[1]}")
    print()

    result = engine.calculate_posterior(patient)

    # ── Print Results ───────────────────────────────────────────────
    print(f"\n📊 Bayesian Posterior Results (ranking by {result['population']})")
    print(f"   Total candidate haplotypes: {result['total_candidate_haplotypes']}")
    print(f"   Total possible pairs:       {result['total_possible_pairs']}")
    print(f"   Entropy (uncertainty):      {result['entropy']} bits")
    print(f"   Populations available:      {', '.join(result['populations_available'])}")
    print()

    if result['pairs']:
        print(f"{'Rank':>5} {'Posterior':>10} {'Cumulative':>10}  Pair + Population Frequencies")
        print(f"   {'='*72}")
        for p in result['pairs'][:5]:
            rank = p['rank']
            post = p['posterior']
            cum = p['cumulative']
            h1 = p['haplotype_1']
            h2 = p['haplotype_2']

            # Format population frequencies
            pop_freqs = p['population_frequencies']
            freqs_str = ' | '.join(
                f"{k}={v:.6f}" for k, v in pop_freqs.items() if v > 0
            )

            # Truncate long labels
            if len(h1) > 60:
                h1 = h1[:57] + '...'
            if len(h2) > 60:
                h2 = h2[:57] + '...'

            print(f"  {rank:>4}  {post:>8.6f}  {cum:>8.6f}   H1: {h1}")
            print(f"  {'':>5} {'':>10} {'':>10}   H2: {h2}")
            print(f"  {'':>5} {'':>10} {'':>10}   💠 {freqs_str}")
            print()
    else:
        print("  ❌ No matching haplotype pairs found.")
        print("  Check that patient alleles exist in the reference database.")

    # Print full JSON payload (top 3 pairs)
    print("\n📦 Full JSON Payload (Top 3 pairs) for API:")
    print("-" * 72)
    payload = {
        'population': result['population'],
        'total_possible_pairs': result['total_possible_pairs'],
        'entropy': result['entropy'],
        'populations_available': result['populations_available'],
        'pairs': [
            {k: p[k] for k in ['rank', 'posterior', 'cumulative',
                                'is_homozygous', 'population_frequencies']}
            for p in result['pairs'][:3]
        ],
    }
    print(pprint.pformat(payload, indent=2, width=100))

    engine.close()
    print("\n✅ Multi-population Bayesian calculation complete.")
