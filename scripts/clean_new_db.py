#!/usr/bin/env python3
"""
clean_new_db.py — HaploStats Data Ingestion (XLSX → Normalized DB)

Reads the original raw spreadsheet (Global_Full_Haplotype_Summary_2018.xlsx)
with STRICT COLUMN ANCHORING by absolute Excel column index.

This spreadsheet has a multi-row header:
  Row 0: Parent Count (per population)
  Row 1: Subject Count (per population)
  Row 2: Parent HapCount (per population)
  Row 3: Population group names  (14 groups × 5 cols = 70 pop cols)
  Row 4: Sub-metric labels        (Hap, Hap, Hap, Family, Sample)
  Row 5: Column sub-headers       (Count, Frequency, Rank, Count, Count)
  Row 6+: Data rows

Features:
  • Column-anchored parsing — no frameshift possible.
  • Auto-correction of CSV-export frameshift (detects DRB1* in DRB345 col).
  • NaN DRB345 → stored as empty/NULL properly.
  • 2-field allele normalisation with DRB3/4/5 gene-prefix preservation.
  • 14 population groups with 5 metrics each.
  • Metadata rows (Parent/Subject/HapCount) preserved in a separate table.
"""

import re
import sqlite3
import sys
from pathlib import Path

import openpyxl

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_XLSX   = PROJECT_ROOT / "data" / "raw" / "Global_Full_Haplotype_Summary_2018.xlsx"
OUTPUT_DB    = PROJECT_ROOT / "db" / "haplostats_normalized.db"

# ── Sheet layout ───────────────────────────────────────────────────────────
HEADER_ROWS = 6          # rows 0-5 are header/metadata
DATA_START  = 6          # data starts at row 6 (0-indexed)

# ── HLA Gene Columns (absolute indices 0-8) ───────────────────────────────
HLA_COL_MAP = [
    (0, "hla_a"),
    (1, "hla_c"),
    (2, "hla_b"),
    (3, "hla_drb345"),
    (4, "hla_drb1"),
    (5, "hla_dqa1"),
    (6, "hla_dqb1"),
    (7, "hla_dpa1"),
    (8, "hla_dpb1"),
]
HLA_DB_COLS = [name for _, name in HLA_COL_MAP]

# ── Population Column Blocks ──────────────────────────────────────────────
# Each block = 5 consecutive columns: (hap_count, hap_freq, hap_rank, family_count, sample_count)
# start_col is the first column index for this population
POP_BLOCKS = [
    (  9, "Global",      "Global"),
    ( 14, "AUSTRIA",     "AUSTRIA"),
    ( 19, "KUWAIT",      "KUWAIT"),
    ( 24, "AFA",         "AFA"),
    ( 29, "ASI",         "ASI"),
    ( 34, "EUR",         "EUR"),
    ( 39, "HIS",         "HIS"),
    ( 44, "OTHER",       "OTHER"),
    ( 49, "ARGENTINA",   "ARGENTINA"),
    ( 54, "EGYPT",       "EGYPT"),
    ( 59, "GERMANY",     "GERMANY"),
    ( 64, "GREECE",      "GREECE"),
    ( 69, "SWITZERLAND", "SWITZERLAND"),
    ( 74, "CZECH",       "CZECH"),
]

SUB_METRICS = ["_hap_count", "_hap_freq", "_hap_rank", "_family_count", "_sample_count"]

# ── Build full list of population database column names ───────────────────
POP_DB_COLS = []
for start_col, code, display in POP_BLOCKS:
    for suffix in SUB_METRICS:
        POP_DB_COLS.append(f"{code}{suffix}")

# ── Allele Normalisation ───────────────────────────────────────────────────

def normalize_allele(raw, locus=""):
    """
    Normalise an allele string to strict 2-field resolution.

    Rules:
      1. Empty / null / whitespace / "Abs" / "-" → ""
      2. Slash-separated ambiguous → take first part
      3. Strip "HLA-" prefix
      4. DRB345: preserve DRB3/DRB4/DRB5 gene prefix (also DRB1 if frameshift)
      5. All other loci: standard prefix + 2-field truncation
      6. Strip non-digit suffixes (N, L, S, SG) from last field
    """
    if raw is None:
        return ""

    val = str(raw).strip()
    if not val or val in ("-", "Abs", "abs", "NULL", "N/A"):
        return ""

    # ── Slash handling ──
    if "/" in val:
        first = val.split("/")[0].strip()
        return normalize_allele(first, locus)

    # ── Strip HLA- ──
    val = re.sub(r"^HLA-", "", val, flags=re.IGNORECASE)

    if "*" not in val:
        return val  # preserve raw sentinel values

    gene_part, allele_part = val.split("*", 1)

    # Strip SG suffix etc from fields
    fields = allele_part.split(":")
    truncated = ":".join(fields[:2])
    # Strip trailing non-digit alpha suffixes from last field
    truncated = re.sub(r"(\d+)[A-Za-z]+$", r"\1", truncated)

    if locus == "hla_drb345":
        # Accept any DRB-prefixed gene
        if gene_part in ("DRB3", "DRB4", "DRB5", "DRB1"):
            pass  # preserve as-is
        elif gene_part in ("3", "4", "5"):
            gene_part = f"DRB{gene_part}"
        return f"{gene_part}*{truncated}"

    # All other loci — keep standard gene prefix
    return f"{gene_part}*{truncated}"


# ── Frameshift Auto-Correction ─────────────────────────────────────────────

def is_drb1_allele(val):
    """Check if a normalized value looks like a DRB1* allele."""
    return bool(re.match(r"^DRB1\*", str(val)))

def auto_correct_row(row_values):
    """
    Detect and fix a frameshift in the 9 HLA columns.

    A frameshift is detected when the DRB345 column contains a DRB1* allele
    (meaning an empty DRB345 in the source caused leftward column shift).

    Fix: Move the DRB1* allele from DRB345 → DRB1,
         and shift remaining locus values rightwards.
    Returns (corrected_row, was_corrected).
    """
    hla_values = row_values[:9]
    corrected = list(hla_values)

    drb345_raw = hla_values[3]
    drb1_raw   = hla_values[4]

    # Check if DRB345 contains a DRB1* allele
    drb345_norm = normalize_allele(drb345_raw, "hla_drb345")
    if is_drb1_allele(drb345_norm):
        # Frameshift detected!
        corrected[3] = ""        # DRB345 = empty
        corrected[4] = drb345_raw  # DRB1 = what was in DRB345
        for i in range(5, 9):
            corrected[i] = hla_values[i - 1]  # shift left
        return (corrected + row_values[9:], True)

    return (row_values, False)


# ── Safe Converters ────────────────────────────────────────────────────────

def safe_float(val):
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in ("-", "", "nan", "NaN"):
        return None
    try:
        return float(s)
    except (ValueError, TypeError):
        return None

def safe_int(val):
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in ("-", "", "nan", "NaN"):
        return None
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return None


# ── SQL Builders ───────────────────────────────────────────────────────────

def build_create_table_sql():
    """Build CREATE TABLE DDL for haplotypes."""
    cols = []
    for _, name in HLA_COL_MAP:
        cols.append(f'  "{name}" TEXT')
    for col_name in POP_DB_COLS:
        if col_name.endswith("_freq") or col_name.endswith("_sample_count") or col_name.endswith("_family_count"):
            cols.append(f'  "{col_name}" REAL')
        else:
            cols.append(f'  "{col_name}" INTEGER')
    return "CREATE TABLE haplotypes (\n" + ",\n".join(cols) + "\n);"

def build_metadata_table_sql():
    return """
    CREATE TABLE IF NOT EXISTS metadata (
      key   TEXT PRIMARY KEY,
      value TEXT
    );
    """

def build_insert_sql():
    all_cols = HLA_DB_COLS + POP_DB_COLS
    quoted = ", ".join(f'"{c}"' for c in all_cols)
    placeholders = ", ".join("?" * len(all_cols))
    return f"INSERT INTO haplotypes ({quoted}) VALUES ({placeholders});"


# ── Unit Tests ─────────────────────────────────────────────────────────────

def run_tests():
    tests = [
        # Basic 4-field → 2-field
        ("A*01:01:01:01",             "hla_a",       "A*01:01"),
        ("C*07:01:01:01",             "hla_c",       "C*07:01"),
        ("B*08:01:01:01",             "hla_b",       "B*08:01"),
        # With HLA- prefix
        ("HLA-A*01:01:01:01",         "hla_a",       "A*01:01"),
        ("HLA-C*07:01:01:01",         "hla_c",       "C*07:01"),
        # 3-field → 2-field
        ("DPB1*04:01:01",             "hla_dpb1",    "DPB1*04:01"),
        # DRB345 - preserve gene
        ("DRB3*01:01:02:01",          "hla_drb345",  "DRB3*01:01"),
        ("DRB4*01:01:01:01",          "hla_drb345",  "DRB4*01:01"),
        ("DRB5*01:01:01",             "hla_drb345",  "DRB5*01:01"),
        ("HLA-DRB3*01:01:02:01",      "hla_drb345",  "DRB3*01:01"),
        # Strip N/L/S/SG suffix
        ("DRB4*01:03:01:02N",         "hla_drb345",  "DRB4*01:03"),
        ("DRB1*03:01:01:01SG",        "hla_drb1",    "DRB1*03:01"),
        # Slash ambiguous
        ("DPB1*04:01:01:01/DPB1*04:01:01:02",
                                     "hla_dpb1",    "DPB1*04:01"),
        ("HLA-DQA1*01:01:01:02/HLA-DQA1*01:01:01:03",
                                     "hla_dqa1",    "DQA1*01:01"),
        # Empty / Abs
        ("",                           "hla_drb345",  ""),
        ("Abs",                        "hla_drb345",  ""),
        ("-",                          "hla_drb345",  ""),
        (None,                         "hla_drb345",  ""),
        # DRB345 → DRB1 (frameshift case)
        ("DRB1*01:01:01",             "hla_drb345",  "DRB1*01:01"),
        ("HLA-DRB1*01:01:01",         "hla_drb345",  "DRB1*01:01"),
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
        print(f"  {status} normalize({str(raw)!r:55s}, {locus:20s}) → {result!r:20s}  (expected {expected!r})")
    print(f"\n  → {'ALL PASS' if all_pass else 'SOME FAILED'}")
    print()
    return all_pass


# ── Main Pipeline ──────────────────────────────────────────────────────────

def run():
    # ── Run unit tests ────────────────────────────────────────────────────
    if not run_tests():
        print("❌ Unit tests failed — aborting.", file=sys.stderr)
        sys.exit(1)

    if not INPUT_XLSX.exists():
        print(f"❌ Input XLSX not found: {INPUT_XLSX}", file=sys.stderr)
        sys.exit(1)

    OUTPUT_DB.parent.mkdir(parents=True, exist_ok=True)

    # ── Read xlsx ─────────────────────────────────────────────────────────
    wb = openpyxl.load_workbook(str(INPUT_XLSX), data_only=True)
    ws = wb.active
    all_rows = list(ws.iter_rows(values_only=True))

    # Parse multi-row header
    parent_counts   = all_rows[0]   # Row 0: Parent Count
    subject_counts  = all_rows[1]   # Row 1: Subject Count
    parent_hapcount = all_rows[2]   # Row 2: Parent HapCount
    pop_names_row   = all_rows[3]   # Row 3: population group names
    sub_metrics_row = all_rows[4]   # Row 4: sub-metric labels
    col_headers_row = all_rows[5]   # Row 5: column sub-headers
    raw_data_rows   = all_rows[6:]  # Row 6+: data

    print(f"\n📄 File: {INPUT_XLSX.name}")
    print(f"   Metadata: Parent={parent_counts[0]}, Subject={subject_counts[0]}, HapCount={parent_hapcount[0]}")
    print(f"   Populations: {len(POP_BLOCKS)} groups")
    print(f"   Data rows: {len(raw_data_rows)}")
    print()

    # ── Stats counters ────────────────────────────────────────────────────
    stats = {
        "total_raw": len(raw_data_rows),
        "corrected": 0,
        "dropped": 0,
        "drb345_nan": 0,       # rows with NaN in DRB345
        "drb345_abs": 0,       # rows with "Abs" in DRB345
        "slash_resolved": 0,
        "drb3_gene": 0,
        "drb4_gene": 0,
        "drb5_gene": 0,
        "drb1_in_drb345": 0,
    }

    # ── Build DB ──────────────────────────────────────────────────────────
    conn = sqlite3.connect(str(OUTPUT_DB))
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS haplotypes;")
    cur.execute("DROP TABLE IF EXISTS metadata;")
    cur.execute(build_create_table_sql())
    cur.execute(build_metadata_table_sql())

    # Store metadata
    meta_key_labels = {
        0: ("parent_count", parent_counts),
        1: ("subject_count", subject_counts),
        2: ("parent_hapcount", parent_hapcount),
    }
    for row_idx, (key_prefix, data_row) in meta_key_labels.items():
        for start_col, code, display in POP_BLOCKS:
            val = data_row[start_col]  # first column of each block
            if val is not None:
                cur.execute(
                    "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
                    (f"{code}_{key_prefix}", str(val))
                )
    conn.commit()

    insert_sql = build_insert_sql()
    batch = []
    corrected_examples = []

    for row_idx, row in enumerate(raw_data_rows):
        # Ensure we have 79 values
        row_values = list(row)[:79]
        while len(row_values) < 79:
            row_values.append(None)

        # ── Step 1: Auto-correct frameshift ──────────────────────────────
        corrected_row, was_corrected = auto_correct_row(row_values)
        if was_corrected:
            stats["corrected"] += 1
            if len(corrected_examples) < 3:
                corrected_examples.append({
                    "raw": list(row_values[:9]),
                    "fixed": list(corrected_row[:9]),
                })
            row_values = corrected_row

        # ── Step 2: Normalise HLA columns ────────────────────────────────
        hla_values = []
        for col_idx, db_name in HLA_COL_MAP:
            raw_val = row_values[col_idx]
            norm = normalize_allele(raw_val, db_name)
            hla_values.append(norm)

        # ── Step 3: Validate — drop hopeless rows ────────────────────────
        valid_loci = sum(1 for v in hla_values[:3] if v and "*" in v)
        if valid_loci == 0:
            stats["dropped"] += 1
            continue

        # ── Step 4: Extract population data ──────────────────────────────
        pop_values = []
        for start_col, code, display in POP_BLOCKS:
            for idx_in_block, suffix in enumerate(SUB_METRICS):
                col_idx = start_col + idx_in_block
                raw_val = row_values[col_idx]
                if suffix == "_hap_freq":
                    pop_values.append(safe_float(raw_val))
                elif suffix in ("_family_count", "_sample_count"):
                    pop_values.append(safe_float(raw_val))  # can be float in xlsx
                else:  # _hap_count, _hap_rank
                    pop_values.append(safe_int(raw_val))

        # ── Step 5: Collect stats ────────────────────────────────────────
        raw_drb345 = row_values[3]
        if raw_drb345 is None or (isinstance(raw_drb345, float) and str(raw_drb345) == "nan"):
            stats["drb345_nan"] += 1
        raw_drb345_str = str(raw_drb345).strip() if raw_drb345 else ""
        if raw_drb345_str in ("Abs", "abs", "-"):
            stats["drb345_abs"] += 1

        # Slash detection across all 9 HLA cols
        for col_idx, _ in HLA_COL_MAP:
            raw_val = str(row_values[col_idx] or "")
            if "/" in raw_val:
                stats["slash_resolved"] += 1

        # DRB subtyping
        drb345_norm = hla_values[3]
        if drb345_norm.startswith("DRB3"):
            stats["drb3_gene"] += 1
        elif drb345_norm.startswith("DRB4"):
            stats["drb4_gene"] += 1
        elif drb345_norm.startswith("DRB5"):
            stats["drb5_gene"] += 1
        elif drb345_norm.startswith("DRB1"):
            stats["drb1_in_drb345"] += 1

        batch.append(tuple(hla_values + pop_values))

    # ── Bulk insert ────────────────────────────────────────────────────────
    print(f"📊 Inserting {len(batch)} rows…")
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
    print(f"  Database:    {OUTPUT_DB}")
    print(f"  Rows:        {row_count}")
    print(f"  Columns:     {len(HLA_DB_COLS) + len(POP_DB_COLS)} ({len(HLA_DB_COLS)} HLA + {len(POP_DB_COLS)} pop)")
    print()
    print("HLA columns:")
    for _, name in HLA_COL_MAP:
        print(f"    → {name}")
    print()
    print("Population columns (14 groups × 5 metrics):")
    for start_col, code, display in POP_BLOCKS:
        print(f"    {display:15s} → {code}_hap_count, {code}_hap_freq, {code}_hap_rank, {code}_family_count, {code}_sample_count")
    print()
    print("Normalisation stats:")
    print(f"  Rows read from xlsx:       {stats['total_raw']}")
    print(f"  Auto-corrected (frameshift): {stats['corrected']}")
    print(f"  Dropped (hopeless):        {stats['dropped']}")
    print(f"  DRB345 NaN (blank):        {stats['drb345_nan']}")
    print(f"  DRB345 = 'Abs':            {stats['drb345_abs']}")
    print(f"  Slash-ambiguous resolved:  {stats['slash_resolved']}")
    print(f"  DRB3 gene entries:         {stats['drb3_gene']}")
    print(f"  DRB4 gene entries:         {stats['drb4_gene']}")
    print(f"  DRB5 gene entries:         {stats['drb5_gene']}")
    print(f"  DRB1 in DRB345 column:     {stats['drb1_in_drb345']}")

    # ── Show auto-corrected examples ──────────────────────────────────────
    if corrected_examples:
        print(f"\n📌 Auto-corrected rows (showing {len(corrected_examples)} examples):")
        for idx, ex in enumerate(corrected_examples):
            raw_vals = ex["raw"]
            fixed_vals = ex["fixed"]
            print(f"\n  Example {idx + 1} (BEFORE → AFTER):")
            for i, (_, name) in enumerate(HLA_COL_MAP):
                arrow = " → " if str(raw_vals[i]) != str(fixed_vals[i]) else "   "
                print(f"    {name:15s} {str(raw_vals[i]):40s}{arrow}{str(fixed_vals[i]):40s}")
    else:
        print("\n  (No auto-corrections needed — all rows clean.)")

    # ── Verify sample rows ────────────────────────────────────────────────
    print("\n📌 Sample of 3 normalised rows (with DRB345 = empty):")
    cur.execute("""
        SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
               hla_dqa1, hla_dqb1, hla_dpa1, hla_dpb1,
               Global_hap_freq, Global_hap_count
        FROM haplotypes
        WHERE hla_drb345 = '' OR hla_drb345 IS NULL
        LIMIT 3
    """)
    rows = cur.fetchall()
    if rows:
        for r in rows:
            print(f"  A={r[0]:8s} C={r[1]:8s} B={r[2]:8s} | DRB345={r[3]:10s} | DRB1={r[4]:10s} | DQA1={r[5]:8s} | DQB1={r[6]:8s} | DPA1={r[7]:8s} | DPB1={r[8]:8s} | Global_freq={r[9]} Global_n={r[10]}")
    else:
        # Fallback: show any 3 rows
        cur.execute("""
            SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
                   hla_dqa1, hla_dqb1, hla_dpa1, hla_dpb1,
                   Global_hap_freq, Global_hap_count
            FROM haplotypes LIMIT 3
        """)
        for r in cur.fetchall():
            print(f"  A={r[0]:8s} C={r[1]:8s} B={r[2]:8s} DRB345={r[3]:8s} DRB1={r[4]:8s} DQA1={r[5]:8s} DQB1={r[6]:8s} DPA1={r[7]:8s} DPB1={r[8]:8s} | Global_freq={r[9]} Global_n={r[10]}")

    # ── Integrity checks ──────────────────────────────────────────────────
    print("\n📌 Integrity checks:")
    issues = 0

    for _, name in HLA_COL_MAP:
        cur.execute(f'SELECT COUNT(*) FROM haplotypes WHERE "{name}" LIKE "%:%:%"')
        bad = cur.fetchone()[0]
        if bad > 0:
            print(f"  ⚠️  {name}: {bad} rows with >2 fields!")
            issues += 1
        else:
            print(f"  ✅ {name}: all 2-field")

    cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_a LIKE 'HLA-%' OR hla_c LIKE 'HLA-%' OR hla_b LIKE 'HLA-%'")
    if cur.fetchone()[0] == 0:
        print("  ✅ No stray HLA- prefixes")

    cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_drb345 LIKE '%/%' OR hla_dpb1 LIKE '%/%'")
    if cur.fetchone()[0] == 0:
        print("  ✅ No stray '/' in normalised data")

    # Verify no DRB1 in DRB345
    cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_drb345 LIKE 'DRB1%'")
    drb1_leak = cur.fetchone()[0]
    if drb1_leak == 0:
        print("  ✅ No DRB1 leaking into DRB345 column")
    else:
        print(f"  ⚠️  {drb1_leak} DRB1 entries still in DRB345!")

    if issues == 0:
        print("\n  🎉 All loci clean at 2-field resolution.")

    conn.close()
    print()


if __name__ == "__main__":
    run()
