#!/usr/bin/env python3
"""
clean_new_db.py — HaploStats Data Ingestion (XLSX → Normalized DB)

Reads the original raw spreadsheet (Table14_ACBDRB345DRB1DQDP_0221.xlsx) with
STRICT COLUMN ANCHORING by absolute Excel column index.  No dynamic array
packing — every column is mapped by its fixed position.

Features:
  • Column-anchored parsing — no frameshift possible.
  • Auto-correction of legacy frameshift rows (e.g. DRB1 allele in DRB345 col).
  • Proper "Abs" / empty handling for DRB345.
  • 2-field allele normalisation with DRB345 gene-prefix preservation.
  • Population mapping: 9 xlsx groups → DB schema compatible with API.
"""

import re
import sqlite3
import sys
from pathlib import Path

import openpyxl

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_XLSX   = PROJECT_ROOT / "data" / "raw" / "Table14_ACBDRB345DRB1DQDP_0221.xlsx"
OUTPUT_DB    = PROJECT_ROOT / "db" / "haplostats_normalized.db"

# ── Column Definitions (by absolute 0-based index in the xlsx) ────────────
# HLA gene columns: indices 0–8
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

# Population columns: (start_col, pop_code, xlsx_label_fq, xlsx_label_n, xlsx_label_rank)
# Each group is 3 columns: Fq, n, rank
POP_COL_MAP = [
    # (xlsx_fq_col, xlsx_n_col, xlsx_rank_col, db_prefix, display_name)
    ( 9, 10, 11, "Global", "Global"),
    (12, 13, 14, "EuAm",   "EuAm"),
    (15, 16, 17, "AFA",    "AFA"),       # AfAm → AFA
    (18, 19, 20, "HIS",    "HIS"),       # USA_Hispanic → HIS
    (21, 22, 23, "EUR",    "EUR"),       # European → EUR
    (24, 25, 26, "Spanish", "Spanish"),
    (27, 28, 29, "Mexican", "Mexican"),
    (30, 31, 32, "ASI",    "ASI"),       # API → ASI
    (33, 34, 35, "Arab",   "Arab"),
]

# Build the list of DB output column names
POP_DB_COLS = []
for fq_col, n_col, rk_col, db_prefix, _ in POP_COL_MAP:
    POP_DB_COLS.append(f"{db_prefix}_hap_freq")
    POP_DB_COLS.append(f"{db_prefix}_hap_count")
    POP_DB_COLS.append(f"{db_prefix}_hap_rank")


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
        else:
            # Unknown prefix — still prefix with DRB
            pass
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
    (meaning an empty DRB345 caused leftward column shift).

    Fix: Move the DRB1* allele from DRB345 → DRB1,
         and shift remaining locus values rightwards.
    Returns (corrected_row, was_corrected).
    """
    hla_values = row_values[:9]  # 9 HLA cols
    corrected = list(hla_values)

    drb345_raw = hla_values[3]  # index 3 = hla_drb345
    drb1_raw   = hla_values[4]  # index 4 = hla_drb1

    # Check if DRB345 contains a DRB1* allele
    drb345_norm = normalize_allele(drb345_raw, "hla_drb345")
    if is_drb1_allele(drb345_norm):
        # Frameshift detected!
        # DRB345 → clear it (set to "")
        # DRB1 → take the DRB1* value that was in DRB345
        # DQA1, DQB1, DPA1, DPB1 → shift left (each takes the value of the next)
        # This means: the whole block DRB1..DPB1 was shifted left by 1
        corrected[3] = ""  # DRB345 = empty
        corrected[4] = drb345_raw  # DRB1 = what was in DRB345
        # Shift remaining: what was in DRB1 goes to DQA1, etc.
        for i in range(5, 9):
            corrected[i] = hla_values[i - 1]  # shift left

        return (corrected + row_values[9:], True)

    return (row_values, False)


# ── Database Builder ───────────────────────────────────────────────────────

def safe_float(val):
    """Convert to float, returning None on empty/bad."""
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
    """Convert to int, returning None on empty/bad."""
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in ("-", "", "nan", "NaN"):
        return None
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return None


def build_create_table_sql():
    """Build CREATE TABLE DDL."""
    cols = []
    for _, name in HLA_COL_MAP:
        cols.append(f'  "{name}" TEXT')
    for col_name in POP_DB_COLS:
        if col_name.endswith("_freq"):
            cols.append(f'  "{col_name}" REAL')
        else:
            cols.append(f'  "{col_name}" INTEGER')
    return "CREATE TABLE haplotypes (\n" + ",\n".join(cols) + "\n);"


def build_insert_sql():
    """Build parameterised INSERT statement."""
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
        # Strip N/L/S/SG suffix
        ("DRB4*01:03:01:02N",         "hla_drb345",  "DRB4*01:03"),
        ("DRB1*03:01:01:01SG",        "hla_drb1",    "DRB1*03:01"),
        # Slash ambiguous
        ("DPB1*04:01:01:01/DPB1*04:01:01:02",
                                     "hla_dpb1",    "DPB1*04:01"),
        # Empty / Abs
        ("",                           "hla_drb345",  ""),
        ("Abs",                        "hla_drb345",  ""),
        ("-",                          "hla_drb345",  ""),
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
        print(f"  {status} normalize({raw!r:55s}, {locus:20s}) → {result!r:20s}  (expected {expected!r})")
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

    # ── Read xlsx with strict column anchoring ────────────────────────────
    wb = openpyxl.load_workbook(str(INPUT_XLSX), data_only=True)
    ws = wb.active
    all_rows = list(ws.iter_rows(values_only=True))

    # Parse header: row 0 is the total haplotypes count, row 1 is column headers
    # Data starts at row index 2
    header_total_haplotypes = all_rows[0][0] if len(all_rows) > 0 else "?"
    header_row = all_rows[1] if len(all_rows) > 1 else []
    raw_data_rows = all_rows[2:]  # data rows

    print(f"\n📄 File: {INPUT_XLSX.name}")
    print(f"   Total haplotypes in header: {header_total_haplotypes}")
    print(f"   Columns (from header row 2): {[str(h) for h in header_row]}")
    print(f"   Data rows: {len(raw_data_rows)}")
    print()

    # ── Stats counters ────────────────────────────────────────────────────
    stats = {
        "total_raw": len(raw_data_rows),
        "corrected": 0,       # auto-corrected frameshift rows
        "dropped": 0,         # hopeless rows dropped
        "drb345_abs": 0,      # DRB345 = "Abs"
        "slash_resolved": 0,  # rows with / in any HLA col
        "drb3_gene": 0,
        "drb4_gene": 0,
        "drb5_gene": 0,
        "drb1_in_drb345": 0, # DRB1 appearing in DRB345 (frameshift indicator)
    }

    # ── Connect / create schema ───────────────────────────────────────────
    conn = sqlite3.connect(str(OUTPUT_DB))
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS haplotypes;")
    cur.execute(build_create_table_sql())
    conn.commit()

    insert_sql = build_insert_sql()
    batch = []
    corrected_examples = []  # store first few corrected rows for display

    for row_idx, row in enumerate(raw_data_rows):
        # Ensure we have 36 values (pad with None)
        row_values = list(row)[:36]
        while len(row_values) < 36:
            row_values.append(None)

        # ── Step 1: Auto-correct frameshift ──────────────────────────────
        corrected_row, was_corrected = auto_correct_row(row_values)
        if was_corrected:
            stats["corrected"] += 1
            if len(corrected_examples) < 3:
                corrected_examples.append({
                    "raw": row_values[:9],
                    "fixed": corrected_row[:9],
                })
            row_values = corrected_row

        # ── Step 2: Normalise HLA columns ────────────────────────────────
        hla_values = []
        for col_idx, db_name in HLA_COL_MAP:
            raw_val = row_values[col_idx]
            norm = normalize_allele(raw_val, db_name)
            hla_values.append(norm)

        # ── Step 3: Validate — drop hopeless rows ────────────────────────
        # A row is hopeless if:
        #   - hla_a is empty/invalid
        #   - None of the first 3 loci (A, C, B) have a valid allele
        valid_loci = sum(1 for v in hla_values[:3] if v and "*" in v)
        if valid_loci == 0:
            stats["dropped"] += 1
            continue

        # ── Step 4: Extract population data ──────────────────────────────
        pop_values = []
        for fq_col, n_col, rk_col, db_prefix, _ in POP_COL_MAP:
            fq = safe_float(row_values[fq_col])
            n  = safe_int(row_values[n_col])
            rk = safe_int(row_values[rk_col])
            pop_values.extend([fq, n, rk])

        # ── Step 5: Collect stats ────────────────────────────────────────
        # Count "Abs" in DRB345
        raw_drb345 = str(row_values[3]).strip() if row_values[3] else ""
        if raw_drb345 in ("Abs", "abs", "-"):
            stats["drb345_abs"] += 1

        # Slash detection across all 9 HLA cols
        for col_idx, _ in HLA_COL_MAP:
            raw = str(row_values[col_idx] or "")
            if "/" in raw:
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

        # ── Assemble row ─────────────────────────────────────────────────
        batch.append(tuple(hla_values + pop_values))

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
    print(f"  Database:    {OUTPUT_DB}")
    print(f"  Rows:        {row_count}")
    print(f"  Columns:     {len(HLA_DB_COLS) + len(POP_DB_COLS)} ({len(HLA_DB_COLS)} HLA + {len(POP_DB_COLS)} pop)")
    print()
    print("HLA columns:")
    for _, name in HLA_COL_MAP:
        print(f"    → {name}")
    print()
    print("Population columns:")
    for _, _, _, db_prefix, display in POP_COL_MAP:
        print(f"    {display:15s} → {db_prefix}_hap_freq, {db_prefix}_hap_count, {db_prefix}_hap_rank")
    print()
    print("Normalisation stats:")
    print(f"  Rows read from xlsx:       {stats['total_raw']}")
    print(f"  Auto-corrected (frameshift): {stats['corrected']}")
    print(f"  Dropped (hopeless):        {stats['dropped']}")
    print(f"  Slash-ambiguous resolved:  {stats['slash_resolved']}")
    print(f"  DRB3 gene entries:         {stats['drb3_gene']}")
    print(f"  DRB4 gene entries:         {stats['drb4_gene']}")
    print(f"  DRB5 gene entries:         {stats['drb5_gene']}")
    print(f"  DRB1 in DRB345 column:     {stats['drb1_in_drb345']}")
    print(f"  DRB345 = 'Abs':            {stats['drb345_abs']}")

    # ── Show auto-corrected example rows ─────────────────────────────────
    if corrected_examples:
        print(f"\n📌 Auto-corrected rows (showing {len(corrected_examples)} examples):")
        for idx, ex in enumerate(corrected_examples):
            raw_vals = ex["raw"]
            fixed_vals = ex["fixed"]
            print(f"\n  Example {idx + 1} (BEFORE → AFTER):")
            for i, (_, name) in enumerate(HLA_COL_MAP):
                arrow = " → " if raw_vals[i] != fixed_vals[i] else "   "
                print(f"    {name:15s} {str(raw_vals[i]):30s}{arrow}{str(fixed_vals[i]):30s}")
    else:
        print("\n  (No auto-corrections needed — all rows clean.)")

    # ── Verify sample rows ────────────────────────────────────────────────
    print("\n📌 Sample of 3 normalised rows:")
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
              f"Global_freq={r[9]} Global_n={r[10]}")

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

    # Verify DRB345 is clean
    cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_drb345 = '' OR hla_drb345 IS NULL")
    empty_count = cur.fetchone()[0]
    print(f"  ✅ DRB345 empty/null count: {empty_count} (expected: {stats['drb345_abs']})")

    if issues == 0:
        print("\n  🎉 All loci clean at 2-field resolution.")

    conn.close()
    print()


if __name__ == "__main__":
    run()
