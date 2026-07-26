#!/usr/bin/env python3
"""
HaploStats — Bayesian Inference Module (Normalized DB)
Phase Final: Posterior probability calculation against the normalized 2-field
reference database (haplostats_normalized.db).

Given an unphased 9-locus HLA genotype, find all possible phased haplotype
pairs and rank them by posterior probability (Hardy-Weinberg 2pq / p²).

Bayes: P(H_pair | G) ∝ P(G | H_pair) × P(H_pair)
  - P(H_pair) = k × P(H1) × P(H2), k = 2 (het), 1 (hom)
  - P(G | H_pair) = 1 if H1+H2 produce G, 0 otherwise
"""

import sqlite3
import math
import sys
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH = BASE_DIR / "db" / "haplostats_normalized.db"

# ── Locus Definitions ──────────────────────────────────────────────

LOCI = [
    "hla_a", "hla_c", "hla_b", "hla_drb345",
    "hla_drb1", "hla_dqa1", "hla_dqb1", "hla_dpa1", "hla_dpb1",
]

LOCUS_LABELS = [
    "HLA-A", "HLA-C", "HLA-B", "HLA-DRB345",
    "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
]

LOCUS_LABEL_MAP = dict(zip(LOCI, LOCUS_LABELS))

# ── Population Frequency Column Map (new normalized DB) ──────────
# Maps display name → SQL column (hap_freq) → display label
POPULATION_COLUMNS = {
    "Global": "Global_hap_freq",
    "AFA":    "AFA_hap_freq",
    "ASI":    "ASI_hap_freq",
    "EUR":    "EUR_hap_freq",
    "HIS":    "HIS_hap_freq",
}

POPULATION_ORDER = ["Global", "AFA", "ASI", "EUR", "HIS"]

# ── Input Sanitizer ────────────────────────────────────────────────

GENE_PREFIX_MAP = {
    "hla_a":      "A",
    "hla_c":      "C",
    "hla_b":      "B",
    "hla_drb345": None,  # special handling below
    "hla_drb1":   "DRB1",
    "hla_dqa1":   "DQA1",
    "hla_dqb1":   "DQB1",
    "hla_dpa1":   "DPA1",
    "hla_dpb1":   "DPB1",
}


def sanitize_allele(raw: str, locus: str) -> str:
    """
    Sanitize a user-typed allele by prepending the correct gene prefix
    if missing. Returns a string matching the normalized DB format.

    For HLA-DRB345, accepts DRB3/DRB4/DRB5 (or 3/4/5) prefixes.
    For all other loci, prepends the standard prefix (e.g. A*, C*, B*).

    Examples:
      "02:01" + hla_a       → "A*02:01"
      "A*02:01" + hla_a     → "A*02:01"
      "DRB3*01:01" + drb345 → "DRB3*01:01"
      "3*01:01" + drb345    → "DRB3*01:01"
      "01:01" + hla_drb1    → "DRB1*01:01"
    """
    val = raw.strip()
    if not val:
        return ""

    # Already has a gene prefix with asterisk
    if "*" in val:
        gene_part = val.split("*", 1)[0]

        # DRB345: normalize shorthand (3→DRB3, 4→DRB4, 5→DRB5)
        if locus == "hla_drb345":
            if gene_part in ("3", "4", "5"):
                return f"DRB{gene_part}*{val.split('*', 1)[1]}"
            # Already has DRB3/4/5 prefix — return as-is
            if gene_part in ("DRB3", "DRB4", "DRB5"):
                return val
            # DRB1 also valid in this column (appears in reference)
            if gene_part == "DRB1":
                return val

        # Any other valid prefix — return as-is (should match DB)
        return val

    # No asterisk — prepend the correct prefix
    if locus == "hla_drb345":
        # DRB345 without a prefix is ambiguous — require prefix from user
        return val  # will fail to match, handled downstream
    else:
        prefix = GENE_PREFIX_MAP.get(locus)
        if prefix:
            return f"{prefix}*{val}"
        return val


# ── Hardy-Weinberg ─────────────────────────────────────────────────

def diplotype_frequency(h1_freq: float, h2_freq: float,
                         is_homozygous: bool) -> float:
    """Compute 2pq (het) or p² (hom) under HWE."""
    if h1_freq <= 0.0 or h2_freq <= 0.0:
        return 0.0
    if is_homozygous:
        return h1_freq * h1_freq
    return 2.0 * h1_freq * h2_freq


# ── HaploMath Engine ────────────────────────────────────────────────

class HaploMath:
    """Bayesian haplotype inference engine (normalized 2-field DB)."""

    def __init__(self, db_path: str = None, population: str = "Global"):
        self.db_path = str(db_path or DB_PATH)
        self.population = population if population in POPULATION_ORDER else "Global"
        self.conn = None
        self._all_haplotypes = None

    def connect(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._load_haplotypes()
        return self

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    def _load_haplotypes(self):
        """Cache all haplotypes with ALL population frequency columns."""
        cur = self.conn.cursor()
        freq_cols = [col for col in POPULATION_COLUMNS.values() if col]
        all_cols = ", ".join(LOCI + freq_cols)
        rows = cur.execute(
            f"SELECT {all_cols} FROM haplotypes "
            f"WHERE Global_hap_freq IS NOT NULL "
            f"ORDER BY Global_hap_freq DESC"
        ).fetchall()
        self._all_haplotypes = [dict(r) for r in rows]
        sys.stderr.write(
            f"  [HaploMath] Loaded {len(self._all_haplotypes)} haplotypes "
            f"({len(freq_cols)} pop freq cols)\n"
        )

    def _haplotype_matches(self, hap: dict, patient_alleles: list,
                            locus: str) -> bool:
        """Check if haplotype carries at least one patient allele at locus."""
        hap_val = hap.get(locus, "")
        if not hap_val:
            return False
        for pa in patient_alleles:
            if pa and pa == hap_val:
                return True
        return False

    def calculate_posterior(self, genotype: dict) -> dict:
        """
        Main entry point: compute all valid phased pairs and rank by posterior.

        Parameters
        ----------
        genotype : dict
            Maps locus → [allele1, allele2]  (already sanitized)

        Returns
        -------
        dict with keys: patient_genotype, population, total_possible_pairs,
                        entropy, populations_available, pairs[]
        """
        constrained_loci = [loc for loc in LOCI if loc in genotype]

        # Step 1: Filter haplotypes to those matching every constrained locus
        candidates = []
        for hap in self._all_haplotypes:
            ok = True
            for loc in constrained_loci:
                a1 = genotype[loc][0]
                a2 = genotype[loc][1] if len(genotype[loc]) > 1 else a1
                if a1 and not self._haplotype_matches(hap, [a1, a2], loc):
                    ok = False
                    break
            if ok:
                candidates.append(hap)

        # Step 2 & 3: Build valid pairs with per-population diplotype freqs
        pairs_raw = []
        n_cand = len(candidates)
        primary_col = POPULATION_COLUMNS.get(self.population, "Global_hap_freq")

        for i in range(n_cand):
            h1 = candidates[i]
            for j in range(i, n_cand):
                h2 = candidates[j]

                # Verify the pair covers both alleles at each constrained locus
                valid = True
                for loc in constrained_loci:
                    a1, a2 = genotype[loc][0], genotype[loc][1] if len(genotype[loc]) > 1 else genotype[loc][0]
                    h1v = h1.get(loc, "")
                    h2v = h2.get(loc, "")
                    # Both patient alleles must be covered by at least one haplotype
                    c1 = (a1 and (a1 == h1v or a1 == h2v))
                    c2 = (a2 and (a2 == h1v or a2 == h2v))
                    if not (c1 and c2):
                        valid = False
                        break
                if not valid:
                    continue

                # Per-population diplotype frequencies
                is_hom = (i == j)
                pop_freqs = {}
                for pop_name, col_name in POPULATION_COLUMNS.items():
                    if col_name is None:
                        pop_freqs[pop_name] = 0.0
                    else:
                        f1 = h1.get(col_name) or 0.0
                        f2 = h2.get(col_name) or 0.0
                        pop_freqs[pop_name] = round(
                            diplotype_frequency(f1, f2, is_hom), 10
                        )

                # Primary population joint probability (for ranking)
                f1p = h1.get(primary_col) or 0.0
                f2p = h2.get(primary_col) or 0.0
                primary_joint = diplotype_frequency(f1p, f2p, is_hom)

                if primary_joint <= 0.0:
                    continue

                h1_label = self._make_label(h1)
                h2_label = self._make_label(h2)

                pairs_raw.append({
                    "haplotype_1": h1_label,
                    "haplotype_2": h2_label,
                    "joint_prob": primary_joint,
                    "is_homozygous": is_hom,
                    "population_frequencies": pop_freqs,
                    "h1_id": i,
                    "h2_id": j,
                })

        # Step 4: Normalize → posterior probabilities
        total = sum(p["joint_prob"] for p in pairs_raw)
        if total == 0.0:
            return {
                "patient_genotype": genotype,
                "population": self.population,
                "total_possible_pairs": 0,
                "entropy": 0.0,
                "populations_available": POPULATION_ORDER,
                "error": "No matching haplotype pairs found.",
                "pairs": [],
            }

        for p in pairs_raw:
            p["posterior"] = round(p["joint_prob"] / total, 10)

        pairs_raw.sort(key=lambda x: x["posterior"], reverse=True)

        # Step 5: Rank
        ranked = []
        cumulative = 0.0
        for rank_idx, p in enumerate(pairs_raw, 1):
            cumulative += p["posterior"]
            ranked.append({
                "rank": rank_idx,
                "haplotype_1": p["haplotype_1"],
                "haplotype_2": p["haplotype_2"],
                "posterior": p["posterior"],
                "cumulative": round(cumulative, 10),
                "is_homozygous": p["is_homozygous"],
                "population_frequencies": p["population_frequencies"],
            })

        # Shannon entropy
        entropy = 0.0
        for p in pairs_raw:
            if p["posterior"] > 0:
                entropy -= p["posterior"] * math.log2(p["posterior"])

        return {
            "patient_genotype": genotype,
            "population": self.population,
            "total_candidate_haplotypes": n_cand,
            "total_possible_pairs": len(ranked),
            "pairs": ranked,
            "entropy": round(entropy, 4),
            "populations_available": POPULATION_ORDER,
        }

    def _make_label(self, hap: dict) -> str:
        """Build a concise haplotype label."""
        parts = []
        for loc, lbl in zip(LOCI, LOCUS_LABELS):
            val = hap.get(loc, "")
            if val:
                parts.append(f"{lbl}={val}")
            else:
                parts.append(f"{lbl}=?")
        return " | ".join(parts)


# ── Main (test) ─────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("  HaploStats — Bayesian Inference (Normalized DB)")
    print("=" * 72)

    engine = HaploMath(population="Global")
    engine.connect()

    # Test with a simple genotype using 2-field alleles
    patient = {
        "hla_a":      ["A*01:01", "A*03:01"],
        "hla_c":      ["C*07:01", "C*07:02"],
        "hla_b":      ["B*08:01", "B*07:02"],
        "hla_drb345": ["DRB3*01:01", "DRB5*01:01"],
        "hla_drb1":   ["DRB1*03:01", "DRB1*15:01"],
        "hla_dqa1":   ["DQA1*05:01", "DQA1*01:02"],
        "hla_dqb1":   ["DQB1*02:01", "DQB1*06:02"],
        "hla_dpa1":   ["DPA1*01:03", "DPA1*01:03"],
        "hla_dpb1":   ["DPB1*04:01", "DPB1*04:01"],
    }

    print("\n📋 Patient Genotype:")
    for loc, alleles in patient.items():
        tag = " (hom)" if alleles[0] == alleles[1] else ""
        print(f"  {loc:15} {alleles[0]:15} / {alleles[1]}{tag}")

    result = engine.calculate_posterior(patient)

    print(f"\n📊 Results (ranked by {result['population']}):")
    print(f"   Candidate haplotypes: {result.get('total_candidate_haplotypes', 'N/A')}")
    print(f"   Possible pairs:       {result['total_possible_pairs']}")
    print(f"   Entropy:              {result.get('entropy', 0)} bits")
    print()

    for p in result["pairs"][:5]:
        print(f"  #{p['rank']} posterior={p['posterior']:.6f} cum={p['cumulative']:.6f}")
        print(f"     H1: {p['haplotype_1'][:80]}")
        print(f"     H2: {p['haplotype_2'][:80]}")
        pf = p["population_frequencies"]
        parts = " | ".join(f"{k}={v:.2e}" for k, v in pf.items() if v > 0)
        print(f"     💠 {parts}")
        print()

    engine.close()
    print("✅ Bayesian engine ready.")
