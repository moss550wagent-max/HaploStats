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

# Tilde-format display order (genomic order: DRB1 before DRB345)
LOCI_TILDE = [
    "hla_a", "hla_c", "hla_b",
    "hla_drb1",
    "hla_drb345",
    "hla_dqa1", "hla_dqb1", "hla_dpa1", "hla_dpb1",
]

# ── Fuzzy Match / Tolerance Scoring ──────────────────────────────
#
# The reference DB (2026 rows) is too small to guarantee a 100% exact
# 9-locus match for every new patient. Instead of discarding rows when a
# flexible locus misses, we apply heuristic tolerance scoring:
#
#   GATEKEEPER loci (A, C, B, DRB345): exact match REQUIRED.
#     Rows failing any gatekeeper are discarded outright.
#   FLEXIBLE loci (DRB1, DQA1, DQB1, DPA1, DPB1): scored 0-20 each
#     on a 100-point scale and summed into a `match_percentage`.
#
GATEKEEPER_LOCI = ["hla_a", "hla_c", "hla_b", "hla_drb345"]
FLEXIBLE_LOCI = ["hla_drb1", "hla_dqa1", "hla_dqb1", "hla_dpa1", "hla_dpb1"]

POINTS_PER_LOCUS = 20      # exact match  (05:01 == 05:01)
ALLELE_GROUP_POINTS = 16   # group match  (05:01 vs 05:02)

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


# ── Fuzzy Match Scoring Utilities ──────────────────────────────────

def _match_score(patient_allele: str, db_allele: str) -> int:
    """
    Score one patient allele against one database allele (2-field strings).

    Exact match        (05:01 == 05:01): 20 points
    Allele group match (05:01 vs 05:02): 16 points  (field 1 equal)
    Total mismatch     (field 1 differs, or either side blank): 0
    """
    pa = (patient_allele or "").strip().upper()
    da = (db_allele or "").strip().upper()
    if not pa or not da:
        return 0
    if pa == da:
        return POINTS_PER_LOCUS
    pf = [f for f in pa.split(":") if f]
    df = [f for f in da.split(":") if f]
    if pf and df and pf[0] == df[0]:
        return ALLELE_GROUP_POINTS
    return 0


def _locus_pair_score(patient_alleles: list, h1_allele: str,
                      h2_allele: str) -> float:
    """
    Score one flexible locus (0-20) for a candidate haplotype pair.

    Each patient allele receives its best score against either haplotype
    allele (H1 or H2). Blank patient alleles are untyped wildcards and
    score full points — an untyped allele must not penalize the pair.
    """
    if not patient_alleles:
        return 0.0
    scores = []
    for pa in patient_alleles:
        if not pa:
            scores.append(float(POINTS_PER_LOCUS))
        else:
            scores.append(float(max(
                _match_score(pa, h1_allele),
                _match_score(pa, h2_allele),
            )))
    return sum(scores) / len(scores)


def _compute_match_percentage(genotype: dict, h1: dict, h2: dict) -> int:
    """
    Heuristic tolerance score (0-100) over the 5 flexible loci
    (DRB1, DQA1, DQB1, DPA1, DPB1) for a candidate haplotype pair.

    Untyped flexible loci are not scored — they get full credit so a
    partially-typed patient can still reach 100% on what they provided.
    """
    total = 0.0
    for loc in FLEXIBLE_LOCI:
        if loc not in genotype:
            total += POINTS_PER_LOCUS
        else:
            total += _locus_pair_score(
                genotype[loc], h1.get(loc, ""), h2.get(loc, "")
            )
    return int(round(total))


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
        """
        Check if haplotype carries at least one patient allele at locus.

        Blank patient alleles ("") are treated as wildcards:
          - If a haplotype also has a blank at this locus (null/absent),
            it matches the blank patient allele.
          - If a haplotype has a real allele and the patient has a blank,
            it does NOT match via the blank — needs a different real match.
        """
        hap_val = hap.get(locus, "")
        for pa in patient_alleles:
            if not pa:
                # Blank patient allele — compatible with blank haplotype
                if not hap_val:
                    return True
                continue  # blank allele doesn't match a non-blank haplotype
            if pa == hap_val:
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
        # Gatekeeper loci (A, C, B, DRB345) — exact match REQUIRED.
        # Flexible loci (DRB1, DQA1, DQB1, DPA1, DPB1) are scored later
        # via tolerance scoring and do NOT gate candidate selection.
        gatekeeper_loci = [loc for loc in GATEKEEPER_LOCI if loc in genotype]

        # Step 1: Filter haplotypes to those matching every gatekeeper locus
        candidates = []
        for hap in self._all_haplotypes:
            ok = True
            for loc in gatekeeper_loci:
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

                # Verify the pair covers both alleles at each GATEKEEPER locus.
                # Flexible loci are NOT verified — they are tolerance-scored.
                valid = True
                for loc in gatekeeper_loci:
                    a1, a2 = genotype[loc][0], genotype[loc][1] if len(genotype[loc]) > 1 else genotype[loc][0]
                    h1v = h1.get(loc, "")
                    h2v = h2.get(loc, "")
                    # Both patient alleles must be covered by at least one haplotype
                    # Blank/null alleles are always satisfied (hemizygous/null scenario)
                    c1 = (not a1) or (a1 == h1v or a1 == h2v)
                    c2 = (not a2) or (a2 == h1v or a2 == h2v)
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

                # Heuristic tolerance score over the 5 flexible loci (0-100)
                match_pct = _compute_match_percentage(genotype, h1, h2)

                pairs_raw.append({
                    "haplotype_1": h1_label,
                    "haplotype_2": h2_label,
                    "joint_prob": primary_joint,
                    "match_percentage": match_pct,
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

        # Sort: match_percentage (descending) FIRST, then population
        # frequency (descending) as the tie-breaker.
        pairs_raw.sort(key=lambda x: (x["match_percentage"], x["joint_prob"]),
                       reverse=True)

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
                "match_percentage": p["match_percentage"],
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
        """
        Build a single-line haplotype label in tilde-separated format.
        Example:  A*01:01 ~ C*07:01 ~ B*08:01 ~ DRB1*03:01 ~ DRB3*01:01
                         ~ DQA1*05:01 ~ DQB1*02:01 ~ DPA1*01:03 ~ DPB1*04:01
        Empty loci are omitted from the string entirely.
        """
        parts = []
        for loc in LOCI_TILDE:
            val = hap.get(loc, "")
            if val:
                parts.append(val)
        return " ~ ".join(parts)


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
        print(f"  #{p['rank']} posterior={p['posterior']:.6f} cum={p['cumulative']:.6f} "
              f"match={p['match_percentage']}%")
        print(f"     H1: {p['haplotype_1']}")
        print(f"     H2: {p['haplotype_2']}")
        pf = p["population_frequencies"]
        parts = " | ".join(f"{k}={v:.2e}" for k, v in pf.items() if v > 0)
        print(f"     💠 {parts}")
        print()

    engine.close()

    # ── Test 2: Hemizygous DRB345 ────────────────────────────────────
    print("\n" + "=" * 72)
    print(f"  ZYGOSITY TEST — Hemizygous DRB345 [DRB3*01:01, '']")
    print("=" * 72)

    engine2 = HaploMath(population="Global")
    engine2.connect()

    # Patient: one haplotype side has DRB3*01:01 (A*01:01/C*07:01/B*08:01/DRB1*03:01)
    # Other side: A*03:01/C*04:01/B*35:01/no DRB345/DRB1*01:01
    # These two haplotype patterns actually exist in the reference DB
    hemizygous_patient = {
        "hla_a":      ["A*01:01", "A*03:01"],
        "hla_c":      ["C*07:01", "C*04:01"],
        "hla_b":      ["B*08:01", "B*35:01"],
        "hla_drb345": ["DRB3*01:01", ""],  # one allele blank!
        "hla_drb1":   ["DRB1*03:01", "DRB1*01:01"],
        "hla_dqa1":   ["DQA1*05:01", "DQA1*01:01"],
        "hla_dqb1":   ["DQB1*02:01", "DQB1*05:01"],
        "hla_dpa1":   ["DPA1*02:01", "DPA1*01:03"],
        "hla_dpb1":   ["DPB1*01:01", "DPB1*04:02"],
    }

    print("\n📋 Hemizygous Patient Genotype:")
    for loc, alleles in hemizygous_patient.items():
        tag = " (empty)" if not alleles[1] else (" (hom)" if alleles[0] == alleles[1] else "")
        print(f"  {loc:15} {alleles[0]:15} / {alleles[1]}{tag}")

    result2 = engine2.calculate_posterior(hemizygous_patient)

    print(f"\n📊 Results (ranked by {result2['population']}):")
    print(f"   Candidate haplotypes: {result2.get('total_candidate_haplotypes', 'N/A')}")
    print(f"   Possible pairs:       {result2['total_possible_pairs']}")
    print(f"   Entropy:              {result2.get('entropy', 0)} bits")
    print()

    # Verify we got both hom and hemizygous scenarios
    hom_count = sum(1 for p in result2['pairs'][:10] if p['is_homozygous'])
    het_count = sum(1 for p in result2['pairs'][:10] if not p['is_homozygous'])
    print(f"  Homozygous DRB345 pairs in top 10:  {hom_count}")
    print(f"  Hemizygous DRB345 pairs in top 10: {het_count}")
    print()

    for p in result2['pairs'][:5]:
        drb345_h1 = " ⍰" if not any(a in p['haplotype_1'] for a in ['DRB3*','DRB4*','DRB5*']) else "  "
        drb345_h2 = " ⍰" if not any(a in p['haplotype_2'] for a in ['DRB3*','DRB4*','DRB5*']) else "  "
        print(f"  #{p['rank']} posterior={p['posterior']:.6f} hom={p['is_homozygous']} "
              f"match={p['match_percentage']}%")
        print(f"     H1{drb345_h1} {p['haplotype_1']}")
        print(f"     H2{drb345_h2} {p['haplotype_2']}")
        print()

    engine2.close()

    # ── Test 3: Fuzzy Tolerance Scoring ──────────────────────────────
    # Same patient as Test 2 but DQB1 allele 2 typed as 05:02 instead of
    # 05:01 — an ALLELE GROUP mismatch (field 1 matches, field 2 differs).
    # The engine must still return pairs, ranked by match % (100% first,
    # then 96% group matches, then lower).
    print("\n" + "=" * 72)
    print("  FUZZY TEST — Allele-group mismatch at DQB1 (05:02 vs 05:01)")
    print("=" * 72)

    engine3 = HaploMath(population="Global")
    engine3.connect()

    fuzzy_patient = dict(hemizygous_patient)
    fuzzy_patient["hla_dqb1"] = ["DQB1*02:01", "DQB1*05:02"]  # group mismatch

    print("\n📋 Fuzzy Patient Genotype (DQB1 allele 2 = 05:02):")
    for loc, alleles in fuzzy_patient.items():
        tag = " (empty)" if not alleles[1] else (" (hom)" if alleles[0] == alleles[1] else "")
        print(f"  {loc:15} {alleles[0]:15} / {alleles[1]}{tag}")

    result3 = engine3.calculate_posterior(fuzzy_patient)

    print(f"\n📊 Results (ranked by match % → {result3['population']} freq):")
    print(f"   Candidate haplotypes: {result3.get('total_candidate_haplotypes', 'N/A')}")
    print(f"   Possible pairs:       {result3['total_possible_pairs']}")
    print(f"   Entropy:              {result3.get('entropy', 0)} bits")
    print()

    pcts = sorted({p['match_percentage'] for p in result3['pairs']}, reverse=True)
    print(f"  Match % values present: {pcts}")

    # Verify strict sort: match % non-increasing down the ranking
    mp_list = [p['match_percentage'] for p in result3['pairs']]
    assert all(mp_list[i] >= mp_list[i+1] for i in range(len(mp_list)-1)), \
        "Sort violation: match_percentage not descending"
    print("  ✅ Sort verified: match_percentage strictly non-increasing")

    for p in result3['pairs'][:6]:
        print(f"  #{p['rank']} posterior={p['posterior']:.6f} match={p['match_percentage']}%")
        print(f"     H1: {p['haplotype_1']}")
        print(f"     H2: {p['haplotype_2']}")
        print()

    engine3.close()
    print("✅ Bayesian engine ready (exact + hemizygous + fuzzy tests passed).")
