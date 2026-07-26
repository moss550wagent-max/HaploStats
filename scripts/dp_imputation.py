#!/usr/bin/env python3
"""
HaploStats — DPB1-Decoupled Imputation (Normalized DB)
Phase Final: DPB1 imputation using DPA1→DPB1 conditional probabilities
from the normalized 2-field reference database.

The DQ-DP recombination hotspot means DPB1 can be decoupled from the core
A-C-B-DRB345-DRB1-DQA1-DQB1 chain. This module provides a clean DPB1
imputation model: P(DPB1 | DPA1) when DPA1 is typed, or marginal P(DPB1)
otherwise.
"""

import sqlite3
import math
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH = BASE_DIR / "db" / "haplostats_normalized.db"

LOCI = ["hla_a", "hla_c", "hla_b", "hla_drb345",
        "hla_drb1", "hla_dqa1", "hla_dqb1", "hla_dpa1", "hla_dpb1"]


class HaploEM:
    """
    EM-based haplotype imputation with decoupled DPB1 model.
    Connects to the normalized 2-field DB.
    """

    def __init__(self, db_path=None, population="Global"):
        self.db_path = str(db_path or DB_PATH)
        self.population = population
        self.freq_col = f"{population}_hap_freq"
        self.conn = None
        self.full_ref = None

        # DPB1 tables
        self.dpb1_marginal = {}
        self.dpa1_to_dpb1 = {}
        self.dpa1_marginal = {}

    def connect(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._load_reference()
        self._build_dp_tables()
        return self

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    def _load_reference(self):
        cur = self.conn.cursor()
        cols = ", ".join(LOCI + [self.freq_col])
        rows = cur.execute(
            f"SELECT {cols} FROM haplotypes "
            f"WHERE {self.freq_col} IS NOT NULL"
        ).fetchall()
        self.full_ref = [dict(r) for r in rows]
        print(f"  [HaploEM] Loaded {len(self.full_ref)} haplotypes (pop={self.population})")

    def _build_dp_tables(self):
        dpa1_counts = defaultdict(float)
        dpa1_dpb1 = defaultdict(lambda: defaultdict(float))
        dpb1_counts = defaultdict(float)

        for h in self.full_ref:
            freq = h.get(self.freq_col, 0.0) or 0.0
            if freq <= 0:
                continue
            dpa1 = h.get("hla_dpa1", "") or ""
            dpb1 = h.get("hla_dpb1", "") or ""

            dpa1_counts[dpa1] += freq
            dpa1_dpb1[dpa1][dpb1] += freq
            dpb1_counts[dpb1] += freq

        # Normalize DPA1→DPB1 conditional
        self.dpa1_to_dpb1 = {}
        for da, db_map in dpa1_dpb1.items():
            total = dpa1_counts.get(da, 0.0)
            if total > 0:
                self.dpa1_to_dpb1[da] = dict(
                    sorted(
                        ((db, cnt / total) for db, cnt in db_map.items()),
                        key=lambda x: -x[1],
                    )
                )

        # Normalize marginal DPB1
        total_dpb1 = sum(dpb1_counts.values())
        if total_dpb1 > 0:
            self.dpb1_marginal = dict(
                sorted(
                    ((k, v / total_dpb1) for k, v in dpb1_counts.items()),
                    key=lambda x: -x[1],
                )
            )

        # DPA1 marginal
        total_dpa1 = sum(dpa1_counts.values())
        if total_dpa1 > 0:
            self.dpa1_marginal = {
                k: v / total_dpa1 for k, v in dpa1_counts.items()
            }

    def impute_dpb1(self, dpa1_h1: str, dpa1_h2: str,
                     patient_dpa1=None, patient_dpb1=None):
        """Return the best DPB1 assignment and confidence weight."""
        # Case 1: Patient has DPB1 typed
        if patient_dpb1:
            return (patient_dpb1, patient_dpb1, 1.0)

        # Case 2: Use P(DPB1 | DPA1) when DPA1 is typed
        if patient_dpa1:
            best_h1, best_h2 = None, None
            w1, w2 = 0.0, 0.0
            for pat_dpa in patient_dpa1:
                if pat_dpa == dpa1_h1:
                    cond = self.dpa1_to_dpb1.get(dpa1_h1, self.dpb1_marginal)
                    if cond:
                        top = list(cond.items())[0]
                        if top[1] > w1:
                            best_h1, w1 = top[0], top[1]
                if pat_dpa == dpa1_h2 and dpa1_h2 != dpa1_h1:
                    cond = self.dpa1_to_dpb1.get(dpa1_h2, self.dpb1_marginal)
                    if cond:
                        top = list(cond.items())[0]
                        if top[1] > w2:
                            best_h2, w2 = top[0], top[1]
            if best_h1 and best_h2:
                return (best_h1, best_h2, (w1 + w2) / 2)
            if self.dpb1_marginal:
                top = list(self.dpb1_marginal.items())[0]
                if not best_h1:
                    best_h1, w1 = top[0], top[1]
                if not best_h2:
                    best_h2, w2 = top[0], top[1]
                if best_h1 and best_h2:
                    return (best_h1, best_h2, max((w1 + w2) / 2, 0.01))

        # Case 3: Marginal only
        if self.dpb1_marginal:
            items = list(self.dpb1_marginal.items())
            top = items[0]
            top2 = items[1] if len(items) > 1 else top
            return (top[0], top2[0], max(top[1], 0.01))

        return (None, None, 0.0)


if __name__ == "__main__":
    engine = HaploEM(population="Global")
    engine.connect()
    print(f"  DPA1→DPB1 tables: {len(engine.dpa1_to_dpb1)} entries")
    print(f"  DPB1 marginal:    {len(engine.dpb1_marginal)} alleles")
    print(f"  Top DPB1:         {list(engine.dpb1_marginal.items())[:3]}")
    engine.close()
