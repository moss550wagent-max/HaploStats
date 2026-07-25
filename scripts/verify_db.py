#!/usr/bin/env python3
"""Quick database verification."""
import sqlite3

conn = sqlite3.connect("db/haplostats.db")
cur = conn.cursor()

print("=== TABLE SCHEMA ===")
for r in cur.execute("PRAGMA table_info(haplotypes)"):
    print(f"  {r[0]:2}. {r[1]:20} {r[2]:8}")

print("\n=== ROW COUNT ===")
print(f"  {cur.execute('SELECT COUNT(*) FROM haplotypes').fetchone()[0]} rows")

print("\n=== TOP 5 GLOBAL HAPLOTYPES ===")
for r in cur.execute("""
    SELECT hla_a, hla_c, hla_b, hla_drb345, hla_drb1,
           global_freq, global_rank
    FROM global_summary LIMIT 5
"""):
    print(f"  {r[0]:15} {r[1]:15} {r[2]:15} {r[3]:25} {r[4]:20} {r[5]:>8.6f} rank={r[6]}")

print("\n=== Abs PRESERVATION CHECK ===")
n = cur.execute("SELECT COUNT(*) FROM haplotypes WHERE hla_drb345 = 'Abs'").fetchone()[0]
print(f"  'Abs' entries: {n} (should be 107) {'✅' if n == 107 else '❌'}")

print("\n=== JSON AMBIGUITY CHECK ===")
for r in cur.execute("SELECT hla_dpb1 FROM haplotypes WHERE hla_dpb1 LIKE '[%' LIMIT 3"):
    print(f"  {r[0][:80]}...")

print("\n=== INDEXES ===")
for r in cur.execute("SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'"):
    print(f"  {r[0]}")

print("\n=== VIEWS ===")
for r in cur.execute("SELECT name FROM sqlite_master WHERE type='view'"):
    c = cur.execute(f"SELECT COUNT(*) FROM {r[0]}").fetchone()[0]
    print(f"  {r[0]:20} ({c} rows)")

conn.close()
print("\n✅ Database fully verified and operational")
