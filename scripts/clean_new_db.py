#!/usr/bin/env python3
"""
clean_new_db.py — HaploStats Data Ingestion & Normalization Reboot

Reads the 2018 Global Full Haplotype Summary (79-column CSV exported from xlsx),
normalises EVERY allele across all loci down to strict 2-field resolution,
and writes a clean SQLite database.

Normalisation rules (applied in order):
  1. Split ambiguous strings on '/', keep first allele.
  2. Strip 'HLA-' prefix.
  3. DRB345 exception: preserve DRB3/DRB4/DRB5 gene identifier (including DRB1
     that appear in this column).
  4. General truncation: keep only the first two colon-separated fields
     (e.g., A*01:01:01:01 → A*01:01; DRB3*01:01:02:01 → DRB3*01:01).
  5. Handle null-allele suffixes (N, L, S etc.) — dropped by field truncation.
"""

import csv
import re
import sqlite3
import sys
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_CSV    = PROJECT_ROOT / "data" / "raw" / "Global_Full_Haplotype_Summary_2018.csv"
OUTPUT_DB    = PROJECT_ROOT / "db" / "haplostats_normalized.db"

# The 9 locus columns in the CSV
LOCUS_COLS = [
    "hla_a", "hla_c", "hla_b", "hla_drb345",
    "hla_drb1", "hla_dqa1", "hla_dqb1", "hla_dpa1", "hla_dpb1",
]

# Population column mapping: each pop has 5 sub-columns
# (hap_count, hap_freq, hap_rank, family_count, sample_count)
POP_GROUPS = [
    ("Global",       "Global_hap_count",  "Global_hap_freq",  "Global_hap_rank",
                      "Global_family_count", "Global_sample_count"),
    ("AUSTRIA",      "AUSTRIA_hap_count", "AUSTRIA_hap_freq", "AUSTRIA_hap_rank",
                      "AUSTRIA_family_count", "AUSTRIA_sample_count"),
    ("KUWAIT",       "KUWAIT_hap_count",  "KUWAIT_hap_freq",  "KUWAIT_hap_rank",
                      "KUWAIT_family_count", "KUWAIT_sample_count"),
    ("AFA",          "AFA_hap_count",     "AFA_hap_freq",     "AFA_hap_rank",
                      "AFA_family_count",    "AFA_sample_count"),
    ("ASI",          "ASI_hap_count",     "ASI_hap_freq",     "ASI_hap_rank",
                      "ASI_family_count",    "ASI_sample_count"),
    ("EUR",          "EUR_hap_count",     "EUR_hap_freq",     "EUR_hap_rank",
                      "EUR_family_count",    "EUR_sample_count"),
    ("HIS",          "HIS_hap_count",     "HIS_hap_freq",     "HIS_hap_rank",
                      "HIS_family_count",    "HIS_sample_count"),
    ("OTHER",        "OTHER_hap_count",   "OTHER_hap_freq",   "OTHER_hap_rank",
                      "OTHER_family_count",  "OTHER_sample_count"),
    ("ARGENTINA",    "ARGENTINA_hap_count", "ARGENTINA_hap_freq", "ARGENTINA_hap_rank",
                      "ARGENTINA_family_count", "ARGENTINA_sample_count"),
    ("EGYPT",        "EGYPT_hap_count",   "EGYPT_hap_freq",   "EGYPT_hap_rank",
                      "EGYPT_family_count",  "EGYPT_sample_count"),
    ("GERMANY",      "GERMANY_hap_count", "GERMANY_hap_freq", "GERMANY_hap_rank",
                      "GERMANY_family_count", "GERMANY_sample_count"),
    ("GREECE",       "GREECE_hap_count",  "GREECE_hap_freq",  "GREECE_hap_rank",
                      "GREECE_family_count", "GREECE_sample_count"),
    ("SWITZERLAND",  "SWITZERLAND_hap_count", "SWITZERLAND_hap_freq",
                      "SWITZERLAND_hap_rank", "SWITZERLAND_family_count",
                      "SWITZERLAND_sample_count"),
    ("CZECH",        "CZECH_hap_count",   "CZECH_hap_freq",   "CZECH_hap_rank",
                      "CZECH_family_count",  "CZECH_sample_count"),
]


# ── Allele Normalisation ───────────────────────────────────────────────────

def normalize_allele(raw: str, locus: str = "") -> str:
    """
    Normalise an allele string to strict 2-field resolution.

    Edge cases handled in order:
      1. Empty / null / whitespace-only → ""
      2. Slash-separated ambiguous → take first part, recurse
      3. Strip 'HLA-' prefix
      4. DRB345 exception: preserve DRB3/DRB4/DRB5 gene identifier
      5. General truncation: keep only first 2 colon-separated fields
      6. Preserve '-' designations like '-' or 'Abs'
    """
    if not raw or not str(raw).strip():
        return ""

    val = str(raw).strip()

    # ── Handle slashes (ambiguous alleles) ───────────────────────────────
    # e.g. "HLA-DRB1*03:01:01:01/HLA-DRB1*03:01:01:02"
    #      "HLA-DQA1*01:02:01:01/HLA-DQA1*01:02:01:03/HLA-DQA1*01:02:01:05"
    if "/" in val:
        first = val.split("/")[0].strip()
        return normalize_allele(first, locus)

    # ── Strip HLA- prefix ────────────────────────────────────────────────
    val = re.sub(r"^HLA-", "", val, flags=re.IGNORECASE)

    # ── Preserve non-allele sentinels ─────────────────────────────────────
    if val in ("-", "Abs", "abs", "NULL", ""):
        return val if val in ("-", "Abs", "abs") else ""

    if "*" not in val:
        return val  # no asterisk — return as-is

    gene_part, allele_part = val.split("*", 1)

    # ── DRB345 exception: preserve the specific gene identifier ──────────
    if locus == "hla_drb345":
        # Normalise gene prefix: "HLA-DRB3" → "DRB3"
        gene_part = re.sub(r"^HLA-", "", gene_part, flags=re.IGNORECASE)

        # Accept DRB1 (sometimes appears in this column)
        if gene_part in ("DRB3", "DRB4", "DRB5", "DRB1", "3", "4", "5"):
            if gene_part in ("3", "4", "5"):
                gene_part = f"DRB{gene_part}"
        # For anything else, keep the gene prefix as-is

        fields = allele_part.split(":")
        truncated = ":".join(fields[:2])
        # Strip trailing non-digit suffixes (N, L, S, etc.) from last field
        truncated = re.sub(r"(\d+)[A-Za-z]+$", r"\1", truncated)
        return f"{gene_part}*{truncated}"

    # ── All other loci: normalise gene prefix + truncate to 2 fields ─────
    gene_part = re.sub(r"^HLA-", "", gene_part, flags=re.IGNORECASE)
    fields = allele_part.split(":")
    truncated = ":".join(fields[:2])
    # Strip trailing non-digit suffixes (N, L, S, etc.) from last field
    truncated = re.sub(r"(\d+)[A-Za-z]+$", r"\1", truncated)
    return f"{gene_part}*{truncated}"


# ── Database Builder ───────────────────────────────────────────────────────

def build_column_defs():
    """Build SQL column definitions for the normalised table."""
    cols = [f'"{c}" TEXT' for c in LOCUS_COLS]
    for grp_name, *sub_cols in POP_GROUPS:
        for sc in sub_cols:
            if sc.endswith("_freq"):
                cols.append(f'"{sc}" REAL')
            elif sc.endswith("_rank") or sc.endswith("_count"):
                cols.append(f'"{sc}" INTEGER')
            else:
                cols.append(f'"{sc}" REAL')
    return ",\n  ".join(cols)


def build_insert_sql():
    """Build parameterised INSERT statement."""
    all_cols = list(LOCUS_COLS)
    for grp_name, *sub_cols in POP_GROUPS:
        all_cols.extend(sub_cols)
    placeholders = ", ".join(["?"] * len(all_cols))
    quoted = ", ".join(f'"{c}"' for c in all_cols)
    return f"INSERT INTO haplotypes ({quoted}) VALUES ({placeholders});"


SAFE_FLOAT_RX = re.compile(r"^-?\d+(\.\d+)?([eE][+-]?\d+)?$")


def safe_float(val):
    """Convert to float, returning 0.0 on failure."""
    if val is None or str(val).strip() in ("", "-"):
        return 0.0
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


def safe_int(val):
    """Convert to int, returning 0 on failure."""
    if val is None or str(val).strip() in ("", "-"):
        return 0
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return 0


# ── Test harness ───────────────────────────────────────────────────────────

def run_tests():
    """Run normalisation tests and print results."""
    tests = [
        # (raw, locus, expected)
        # 4-field → 2-field
        ("HLA-A*01:01:01:01",            "hla_a",       "A*01:01"),
        ("HLA-C*07:01:01:01",            "hla_c",       "C*07:01"),
        ("HLA-B*08:01:01:01",            "hla_b",       "B*08:01"),

        # 3-field → 2-field
        ("HLA-DPB1*04:01:01",            "hla_dpb1",    "DPB1*04:01"),

        # 2-field → 2-field (unchanged)
        ("HLA-DRB1*01:03",               "hla_drb1",    "DRB1*01:03"),

        # Slash-separated ambiguous → take first
        ("HLA-DRB1*03:01:01:01/HLA-DRB1*03:01:01:02",
                                          "hla_drb1",    "DRB1*03:01"),
        ("HLA-DQA1*01:02:01:01/HLA-DQA1*01:02:01:03/HLA-DQA1*01:02:01:05",
                                          "hla_dqa1",    "DQA1*01:02"),
        ("HLA-DPB1*04:01:01:01/HLA-DPB1*04:01:01:02",
                                          "hla_dpb1",    "DPB1*04:01"),

        # DRB345 → preserve DRB3
        ("HLA-DRB3*01:01:02:01",         "hla_drb345",  "DRB3*01:01"),
        ("HLA-DRB3*02:02:01:01",         "hla_drb345",  "DRB3*02:02"),
        ("HLA-DRB3*03:01:01",            "hla_drb345",  "DRB3*03:01"),

        # DRB345 → preserve DRB4
        ("HLA-DRB4*01:01:01:01",         "hla_drb345",  "DRB4*01:01"),
        ("HLA-DRB4*01:03:01:01",         "hla_drb345",  "DRB4*01:03"),
        ("HLA-DRB4*01:03:01:02N",        "hla_drb345",  "DRB4*01:03"),

        # DRB345 → preserve DRB5
        ("HLA-DRB5*01:01:01",            "hla_drb345",  "DRB5*01:01"),
        ("HLA-DRB5*01:02",               "hla_drb345",  "DRB5*01:02"),
        ("HLA-DRB5*01:08N",              "hla_drb345",  "DRB5*01:08"),

        # DRB345 → DRB1 (also appears in this column)
        ("HLA-DRB1*01:01:01",            "hla_drb345",  "DRB1*01:01"),

        # Slash in DRB345
        ("HLA-DRB3*01:01:02:01/HLA-DRB3*01:01:02:02",
                                          "hla_drb345",  "DRB3*01:01"),

        # Empty/null
        ("",                              "hla_a",       ""),
        ("  ",                            "hla_a",       ""),
    ]

    print("=" * 72)
    print("UNIT TESTS — Allele Normalisation")
    print("=" * 72)
    all_pass = True
    for raw, locus, expected in tests:
        result = normalize_allele(raw, locus)
        status = "✅" if result == expected else "❌"
        if result != expected:
            all_pass = False
        print(f"  {status} normalize({raw!r:55s}, {locus:15s}) → {result!r:20s}  (expected {expected!r})")

    print(f"\n  → {'ALL PASS' if all_pass else 'SOME FAILED'}")
    print()
    return all_pass


# ── Main pipeline ──────────────────────────────────────────────────────────

def run():
    # ── Run unit tests first ──────────────────────────────────────────────
    tests_pass = run_tests()
    if not tests_pass:
        print("❌ Unit tests failed — aborting.", file=sys.stderr)
        sys.exit(1)

    if not INPUT_CSV.exists():
        print(f"❌ Input CSV not found: {INPUT_CSV}", file=sys.stderr)
        sys.exit(1)

    OUTPUT_DB.parent.mkdir(parents=True, exist_ok=True)

    # ── Connect / create schema ───────────────────────────────────────────
    conn = sqlite3.connect(str(OUTPUT_DB))
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS haplotypes;")
    ddl = f"""
    CREATE TABLE haplotypes (
      {build_column_defs()}
    );
    """
    cur.execute(ddl)
    conn.commit()

    # ── Read & transform CSV ──────────────────────────────────────────────
    with open(INPUT_CSV, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    print(f"\n📄 Read {len(rows)} rows from {INPUT_CSV.name}")
    print()

    # ── Collect stats ─────────────────────────────────────────────────────
    stats = {
        "total": len(rows),
        "slash_resolved": 0,
        "drb3_gene": 0,
        "drb4_gene": 0,
        "drb5_gene": 0,
        "drb1_in_drb345": 0,
        "drb345_empty": 0,
    }

    # Build list of all population column names
    all_pop_cols = []
    for grp_name, *sub_cols in POP_GROUPS:
        all_pop_cols.extend(sub_cols)

    insert_sql = build_insert_sql()
    batch = []

    for row in rows:
        norm_loci = []
        for locus in LOCUS_COLS:
            raw = row.get(locus, "")
            norm = normalize_allele(raw, locus)
            norm_loci.append(norm)

            # Track stats (first allele only)
            if locus == "hla_drb345":
                if "/" in raw:
                    stats["slash_resolved"] += 1
                first_allele = raw.split("/")[0].strip() if "/" in raw else raw.strip()
                first_allele = re.sub(r"^HLA-", "", first_allele)
                if first_allele.startswith("DRB3"):
                    stats["drb3_gene"] += 1
                elif first_allele.startswith("DRB4"):
                    stats["drb4_gene"] += 1
                elif first_allele.startswith("DRB5"):
                    stats["drb5_gene"] += 1
                elif first_allele.startswith("DRB1"):
                    stats["drb1_in_drb345"] += 1
                elif not first_allele or first_allele in ("-", "Abs"):
                    stats["drb345_empty"] += 1
            elif locus == "hla_dpb1" and "/" in raw:
                stats["slash_resolved"] += 1

        # Population data
        pop_vals = []
        for col in all_pop_cols:
            raw_val = row.get(col, "")
            if col.endswith("_freq"):
                pop_vals.append(safe_float(raw_val))
            else:  # _count or _rank
                pop_vals.append(safe_int(raw_val))

        batch.append(tuple(norm_loci + pop_vals))

    # ── Bulk insert ────────────────────────────────────────────────────────
    cur.executemany(insert_sql, batch)
    conn.commit()

    # ── Create indexes ────────────────────────────────────────────────────
    print("📊 Creating indexes…")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_global_freq ON haplotypes(Global_hap_freq DESC);")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_hla_a ON haplotypes(hla_a);")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_hla_b ON haplotypes(hla_b);")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_hla_drb1 ON haplotypes(hla_drb1);")
    conn.commit()

    # ── Summary ───────────────────────────────────────────────────────────
    row_count = cur.execute("SELECT COUNT(*) FROM haplotypes").fetchone()[0]

    print(f"\n{'=' * 72}")
    print(f"✅ NORMALISATION COMPLETE")
    print(f"{'=' * 72}")
    print(f"  Database:   {OUTPUT_DB}")
    print(f"  Rows:       {row_count}")
    print(f"  Columns:    {len(LOCUS_COLS) + len(all_pop_cols)}")
    print()
    print("Table columns:")
    print(f"  Locus columns: {', '.join(LOCUS_COLS)}")
    print(f"  Population groups: {len(POP_GROUPS)} populations × 5 metrics each")
    for grp_name, *sub_cols in POP_GROUPS:
        print(f"    {grp_name:15s} → {', '.join(sub_cols)}")
    print()
    print("Normalisation stats:")
    print(f"  Slash-ambiguous resolved:  {stats['slash_resolved']}")
    print(f"  DRB3 gene entries:         {stats['drb3_gene']}")
    print(f"  DRB4 gene entries:         {stats['drb4_gene']}")
    print(f"  DRB5 gene entries:         {stats['drb5_gene']}")
    print(f"  DRB1 in DRB345 column:     {stats['drb1_in_drb345']}")
    if stats["drb345_empty"] > 0:
        print(f"  DRB345 empty/absent:       {stats['drb345_empty']}")
    print()

    # ── Verify sample rows ────────────────────────────────────────────────
    print("Sample of 3 normalised rows:")
    cur.execute("""
        SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
               hla_dqa1, hla_dqb1, hla_dpa1, hla_dpb1,
               Global_hap_freq, Global_hap_count
        FROM haplotypes
        LIMIT 3
    """)
    for r in cur.fetchall():
        print(f"  A={r[0]:8s} C={r[1]:8s} B={r[2]:8s} DRB345={r[3]:10s} "
              f"DRB1={r[4]:10s} DQA1={r[5]:8s} DQB1={r[6]:8s} "
              f"DPA1={r[7]:8s} DPB1={r[8]:8s} | "
              f"Global_freq={r[9]:.6e} Global_n={r[10]}")

    # ── Integrity checks ──────────────────────────────────────────────────
    print("\nIntegrity checks:")
    issues = 0
    for col in LOCUS_COLS:
        cur.execute(f'SELECT COUNT(*) FROM haplotypes WHERE "{col}" LIKE "%:%:%"')
        bad = cur.fetchone()[0]
        if bad > 0:
            print(f"  ⚠️  {col}: {bad} rows with >2 fields!")
            issues += 1
        else:
            print(f"  ✅ {col}: all 2-field")

    cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_a LIKE 'HLA-%' OR hla_c LIKE 'HLA-%' OR hla_b LIKE 'HLA-%'")
    if cur.fetchone()[0] == 0:
        print("  ✅ No stray HLA- prefixes")

    cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_drb345 LIKE '%/%' OR hla_dpb1 LIKE '%/%'")
    if cur.fetchone()[0] == 0:
        print("  ✅ No stray '/' in normalised data")

    if issues == 0:
        print("\n  🎉 All loci clean at 2-field resolution.")

    conn.close()


if __name__ == "__main__":
    run()
