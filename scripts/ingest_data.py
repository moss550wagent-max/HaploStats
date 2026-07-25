#!/usr/bin/env python3
"""
HaploStats Data Ingestion Script
Phase 2: Clean the Table14 9-locus XLSX and ingest into SQLite.

Rules enforced:
  1. Dual-header extraction — pull N values from header row metadata
  2. 'Abs' biological rule — never convert Abs to NaN/NULL
  3. Ambiguity splitting — /-separated alleles → JSON arrays for SQLite
"""

import pandas as pd
import numpy as np
import json
import sqlite3
import re
import os
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).resolve().parent.parent
RAW_PATH   = BASE_DIR / "data" / "raw" / "Table14_ACBDRB345DRB1DQDP_0221.xlsx"
DB_PATH    = BASE_DIR / "db" / "haplostats.db"
CLEAN_CSV  = BASE_DIR / "data" / "clean" / "haplotypes_clean.csv"
SUMMARY_TXT = BASE_DIR / "data" / "clean" / "ingestion_summary.txt"

os.makedirs(BASE_DIR / "data" / "clean", exist_ok=True)
os.makedirs(BASE_DIR / "db", exist_ok=True)

print("=" * 60)
print("HaploStats — Data Ingestion Pipeline")
print("=" * 60)

# ── 1. Read raw Excel ─────────────────────────────────────────────────────
print("\n[1/6] Reading raw Excel…")
df_raw = pd.read_excel(RAW_PATH, sheet_name="UPHD_9_loci", header=None)
print(f"  Raw shape: {df_raw.shape}")

# ── 2. Dual-Header Extraction ─────────────────────────────────────────────
print("\n[2/6] Extracting dual-header metadata…")

# Row 0 (index 0) = junk header with embedded N values
# Row 1 (index 1) = the actual HLA locus and frequency column names
header_row = df_raw.iloc[0].tolist()
name_row   = df_raw.iloc[1].tolist()

# The first 9 columns are the HLA loci (A, C, B, DRB345, DRB1, DQA1, DQB1, DPA1, DPB1)
# After that, every group of 3 columns is: Population_Fq, Population_n, Population_rank
# The N values (6582, 4384, 756, 96, 424, 530, 208, 80, 104) are scattered in the
# header_row at positions [9], [12], [15], [18], [21], [24], [27], [30], [33]

hla_loci = name_row[:9]  # First 9 = HLA names

# Population names derived from name_row positions 9 onward (every 3 cols => Fq, n, rank)
pop_names_raw = []
for i in range(9, len(name_row), 3):
    col_name = str(name_row[i])
    if '_Fq' in col_name:
        pop_names_raw.append(col_name.replace('_Fq', ''))

print(f"  HLA loci detected: {hla_loci}")
print(f"  Populations detected: {pop_names_raw}")

# Extract N values from the header row
# They sit at positions where numeric values appear after the label columns
n_values = {}
n_positions = [9, 12, 15, 18, 21, 24, 27, 30, 33]
for idx, pos in enumerate(n_positions):
    if pos < len(pop_names_raw) * 3 + 9:
        pop_name = pop_names_raw[idx] if idx < len(pop_names_raw) else f"Pop{idx}"
        val = header_row[pos]
        try:
            n_values[pop_name] = int(val)
        except (ValueError, TypeError):
            n_values[pop_name] = None

print(f"  Population sample sizes (N): {n_values}")

# ── 3. Build Clean Column Names ───────────────────────────────────────────
print("\n[3/6] Building clean column schema…")

# Locus columns (TEXT)
locus_cols = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
              'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

# Frequency/Count/Rank columns per population
stat_cols = []
for pop in pop_names_raw:
    stat_cols.append(f"{pop}_freq")
    stat_cols.append(f"{pop}_n")
    stat_cols.append(f"{pop}_rank")

all_cols = locus_cols + stat_cols

print(f"  Total columns: {len(all_cols)}")
print(f"  Locus columns: {len(locus_cols)}")
print(f"  Stat columns:  {len(stat_cols)}")

# ── 4. Extract Clean Data ─────────────────────────────────────────────────
print("\n[4/6] Extracting clean data rows…")

# Data starts at row 2 (0-indexed), i.e., skip first 2 header rows
df_data = df_raw.iloc[2:].copy()
df_data.columns = all_cols[:len(df_data.columns)]
df_data = df_data.reset_index(drop=True)

print(f"  Data rows extracted: {len(df_data)}")

# ── 5. Apply Biological Rules ─────────────────────────────────────────────
print("\n[5/6] Applying biological data rules…")

# 5a. Force 'Abs' to stay as literal "Abs" (never NaN)
# Check the DRB345 column
abs_before = (df_data['hla_drb345'] == 'Abs').sum()
print(f"  'Abs' entries in hla_drb345 before: {abs_before}")

# Replace any NaN/None in string columns with empty string temporarily
# but preserve "Abs"
# First, make sure Abs stays
df_data['hla_drb345'] = df_data['hla_drb345'].fillna('__MISSING__')
df_data['hla_drb345'] = df_data['hla_drb345'].replace('__MISSING__', '')  # real missing stays
# "Abs" is already "Abs" from the XLSX

abs_after = (df_data['hla_drb345'] == 'Abs').sum()
print(f"  'Abs' entries in hla_drb345 after:  {abs_after}")

# Verify Abs wasn't lost
if abs_after < abs_before:
    print("  ⚠️  WARNING: Some 'Abs' entries were lost!")
else:
    print("  ✅ 'Abs' entries preserved successfully")

# 5b. Ambiguity splitting — convert /-separated alleles to JSON arrays
locus_cols_for_split = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
                        'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

ambiguity_stats = {}
for col in locus_cols_for_split:
    # Find entries with /
    has_slash = df_data[col].str.contains('/', na=False)
    count_slash = has_slash.sum()
    ambiguity_stats[col] = int(count_slash)
    
    if count_slash > 0:
        # Replace with JSON array
        df_data[col] = df_data[col].apply(
            lambda x: json.dumps([allele.strip() for allele in str(x).split('/')])
            if pd.notna(x) and '/' in str(x) else x
        )

print(f"  Ambiguity splitting stats:")
for col, cnt in ambiguity_stats.items():
    if cnt > 0:
        print(f"    {col}: {cnt} ambiguous entries → JSON arrays")

total_ambiguous = sum(ambiguity_stats.values())
print(f"  Total ambiguous entries split: {total_ambiguous}")

# 5c. Convert frequency columns to float, n and rank to int
for col in df_data.columns:
    if col.endswith('_freq'):
        df_data[col] = pd.to_numeric(df_data[col], errors='coerce')
    elif col.endswith('_n') or col.endswith('_rank'):
        df_data[col] = pd.to_numeric(df_data[col], errors='coerce').astype(pd.Int64Dtype())

# Convert NaN frequencies to None
for col in df_data.columns:
    if col.endswith('_freq'):
        df_data[col] = df_data[col].where(df_data[col].notna(), None)

print("  ✅ Numeric conversions complete")

# ── 6. Write Outputs ──────────────────────────────────────────────────────
print("\n[6/6] Writing outputs…")

# 6a. Write clean CSV
df_data.to_csv(CLEAN_CSV, index=False)
csv_size = CLEAN_CSV.stat().st_size
print(f"  ✅ Clean CSV written: {CLEAN_CSV.name} ({csv_size:,} bytes)")

# 6b. Define SQLite schema
schema_sql = """
CREATE TABLE IF NOT EXISTS haplotypes (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    hla_a       TEXT,
    hla_c       TEXT,
    hla_b       TEXT,
    hla_drb345  TEXT,
    hla_drb1    TEXT,
    hla_dqa1    TEXT,
    hla_dqb1    TEXT,
    hla_dpa1    TEXT,
    hla_dpb1    TEXT,
    global_freq   REAL,
    global_n      INTEGER,
    global_rank   INTEGER,
    euam_freq     REAL,
    euam_n        INTEGER,
    euam_rank     INTEGER,
    afam_freq     REAL,
    afam_n        INTEGER,
    afam_rank     INTEGER,
    usa_hispanic_freq   REAL,
    usa_hispanic_n      INTEGER,
    usa_hispanic_rank   INTEGER,
    european_freq   REAL,
    european_n      INTEGER,
    european_rank   INTEGER,
    spanish_freq    REAL,
    spanish_n       INTEGER,
    spanish_rank    INTEGER,
    mexican_freq    REAL,
    mexican_n       INTEGER,
    mexican_rank    INTEGER,
    api_freq        REAL,
    api_n           INTEGER,
    api_rank        INTEGER,
    arab_freq       REAL,
    arab_n          INTEGER,
    arab_rank       INTEGER
);

CREATE INDEX IF NOT EXISTS idx_hla_a ON haplotypes(hla_a);
CREATE INDEX IF NOT EXISTS idx_hla_b ON haplotypes(hla_b);
CREATE INDEX IF NOT EXISTS idx_hla_drb1 ON haplotypes(hla_drb1);
CREATE INDEX IF NOT EXISTS idx_global_freq ON haplotypes(global_freq DESC);
"""

# 6c. Create SQLite database and write data
conn = sqlite3.connect(str(DB_PATH))
cursor = conn.cursor()

# Create tables
cursor.executescript(schema_sql)

# Write data
# Drop existing table and write data
cursor.execute("DROP TABLE IF EXISTS haplotypes")
df_data.to_sql('haplotypes', conn, if_exists='append', index=False)

# Recreate indexes (pandas to_sql drops them)
cursor.execute("CREATE INDEX IF NOT EXISTS idx_hla_a ON haplotypes(hla_a)")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_hla_c ON haplotypes(hla_c)")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_hla_b ON haplotypes(hla_b)")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_hla_drb1 ON haplotypes(hla_drb1)")
cursor.execute('CREATE INDEX IF NOT EXISTS idx_global_freq ON haplotypes("Global_freq" DESC)')
cursor.execute('CREATE INDEX IF NOT EXISTS idx_afam_freq ON haplotypes("AfAm_freq" DESC)')
cursor.execute('CREATE INDEX IF NOT EXISTS idx_euam_freq ON haplotypes("EuAm_freq" DESC)')

# Verify
row_count = cursor.execute("SELECT COUNT(*) FROM haplotypes").fetchone()[0]
col_count = len(cursor.execute("PRAGMA table_info(haplotypes)").fetchall())

# Create views for fast population lookup
cursor.execute("""
    CREATE VIEW IF NOT EXISTS global_summary AS
    SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
           hla_dqa1, hla_dqb1, hla_dpa1, hla_dpb1,
           global_freq, global_rank
    FROM haplotypes
    WHERE global_freq IS NOT NULL
    ORDER BY global_freq DESC
""")

cursor.execute("""
    CREATE VIEW IF NOT EXISTS afam_summary AS
    SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
           hla_dqa1, hla_dqb1, hla_dpa1, hla_dpb1,
           afam_freq, afam_rank
    FROM haplotypes
    WHERE afam_freq IS NOT NULL
    ORDER BY afam_freq DESC
""")

cursor.execute("""
    CREATE VIEW IF NOT EXISTS euam_summary AS
    SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
           hla_dqa1, hla_dqb1, hla_dpa1, hla_dpb1,
           euam_freq, euam_rank
    FROM haplotypes
    WHERE euam_freq IS NOT NULL
    ORDER BY euam_freq DESC
""")

conn.commit()

db_size = DB_PATH.stat().st_size
print(f"  ✅ SQLite database written: {DB_PATH.name} ({db_size:,} bytes)")
print(f"  ✅ Rows ingested: {row_count}")
print(f"  ✅ Columns: {col_count}")
print(f"  ✅ Views created: global_summary, afam_summary, euam_summary")

# ── 7. Write Summary ──────────────────────────────────────────────────────
summary_lines = [
    "=" * 60,
    "HaploStats — Ingestion Summary",
    "=" * 60,
    f"Data source: Table14_ACBDRB345DRB1DQDP_0221.xlsx",
    f"Sheet: UPHD_9_loci",
    f"",
    f"--- Structural ---",
    f"Total haplotypes ingested: {row_count}",
    f"Total database columns:    {col_count}",
    f"SQLite database path:      {DB_PATH}",
    f"Clean CSV path:            {CLEAN_CSV}",
    f"",
    f"--- HLA Loci ---",
    f"{', '.join(locus_cols)}",
    f"",
    f"--- Populations ---",
    f"DB field (prefix)  |  N (sample size)",
    f"-------------------+---------------",
]

for pop in pop_names_raw:
    n_val = n_values.get(pop, '?')
    summary_lines.append(f"{pop:<18} | {n_val}")

summary_lines.append("")
summary_lines.append("--- Biological Rule Enforcement ---")
summary_lines.append(f"'Abs' entries preserved in hla_drb345: {abs_after}")
summary_lines.append(f"Ambiguous (/ ) alleles split into JSON: {total_ambiguous}")

total_alleles_split_by_locus = {col: cnt for col, cnt in ambiguity_stats.items() if cnt > 0}
if total_alleles_split_by_locus:
    summary_lines.append("Per locus:")
    for col, cnt in total_alleles_split_by_locus.items():
        summary_lines.append(f"  {col}: {cnt}")

summary_lines.append("")
summary_lines.append("--- Status ---")
summary_lines.append("✅ Ingestion complete. Database ready for query.")

with open(SUMMARY_TXT, 'w') as f:
    f.write('\n'.join(summary_lines))

print(f"\n{'=' * 60}")
print("INGESTION COMPLETE")
print(f"{'=' * 60}")
with open(SUMMARY_TXT, 'r') as f:
    print(f.read())

conn.close()
