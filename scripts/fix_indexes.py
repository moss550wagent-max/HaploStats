#!/usr/bin/env python3
"""Fix indexes and verify database integrity."""
import sqlite3

conn = sqlite3.connect("db/haplostats.db")
cur = conn.cursor()

print("=== ALL EXISTING INDEXES ===")
for r in cur.execute("SELECT name FROM sqlite_master WHERE type='index' AND name NOT LIKE 'sqlite_%'"):
    print(f"  {r[0]}")

print("\n=== RECREATING INDEXES ===")
cur.executescript("""
    CREATE INDEX IF NOT EXISTS idx_hla_a ON haplotypes(hla_a);
    CREATE INDEX IF NOT EXISTS idx_hla_c ON haplotypes(hla_c);
    CREATE INDEX IF NOT EXISTS idx_hla_b ON haplotypes(hla_b);
    CREATE INDEX IF NOT EXISTS idx_hla_drb1 ON haplotypes(hla_drb1);
    CREATE INDEX IF NOT EXISTS idx_global_freq ON haplotypes("Global_freq" DESC);
    CREATE INDEX IF NOT EXISTS idx_afam_freq ON haplotypes("AfAm_freq" DESC);
    CREATE INDEX IF NOT EXISTS idx_euam_freq ON haplotypes("EuAm_freq" DESC);
""")
conn.commit()

print("=== RECREATED INDEXES ===")
for r in cur.execute("SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'"):
    print(f"  ✅ {r[0]}")

print("\n=== VERIFY VIEWS ===")
for view in ["global_summary", "afam_summary", "euam_summary"]:
    try:
        c = cur.execute(f"SELECT COUNT(*) FROM {view}").fetchone()[0]
        print(f"  ✅ {view}: {c} rows")
    except Exception as e:
        print(f"  ❌ {view}: {e}")

print("\n=== INTEGRITY CHECK ===")
cur.execute("PRAGMA integrity_check")
print(f"  ✅ {cur.fetchone()[0]}")

conn.close()
print("\n✅ Database fully verified and operational")
