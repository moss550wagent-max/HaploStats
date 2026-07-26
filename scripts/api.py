#!/usr/bin/env python3
"""
HaploStats — FastAPI Web Service (Normalized DB)
Phase Final: Population-aware REST API returning multi-population
diplotype frequencies (2pq / p²) for every imputed haplotype pair.

Connects to db/haplostats_normalized.db (strict 2-field resolution).

Endpoints:
  POST /impute  — submit unphased patient genotype → ranked phased haplotypes
  GET  /health  — service health check
  GET  /        — population explorer web dashboard
"""

import sys
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Optional
import uvicorn

from scripts.bayesian_calc import HaploMath, POPULATION_ORDER, sanitize_allele, LOCUS_LABEL_MAP

# ── App Initialisation ─────────────────────────────────────────────

app = FastAPI(
    title="HaploStats — HLA Haplotype Population Explorer",
    description="Clinical-grade Bayesian haplotype imputation engine "
                "returning per-population diplotype frequencies (2pq / p²).",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

STATIC_DIR = HERE.parent / "static"
if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.get("/")
def serve_dashboard():
    idx = STATIC_DIR / "index.html"
    if idx.exists():
        return FileResponse(str(idx))
    return {"status": "ok", "message": "HaploStats API running."}


# ── Locus Constants ────────────────────────────────────────────────

ALL_LOCI = [
    "hla_a", "hla_c", "hla_b", "hla_drb345",
    "hla_drb1", "hla_dqa1", "hla_dqb1", "hla_dpa1", "hla_dpb1",
]


# ── Global Engine ─────────────────────────────────────────────────

_engine: Optional[HaploMath] = None


def get_engine(population: str = "Global") -> HaploMath:
    global _engine
    if _engine is None:
        _engine = HaploMath(population=population)
        _engine.connect()
    return _engine


# ── Input Sanitizer ────────────────────────────────────────────────

def sanitize_genotype(genotype: dict) -> dict:
    """
    Intercept user input and prepend correct gene identifiers.

    Rules:
      - If allele has no '*' prefix, add the correct gene prefix
        (e.g., "02:01" for hla_a → "A*02:01")
      - For hla_drb345, require DRB3/DRB4/DRB5 (or 3/4/5) prefix
      - Preserve already-correctly-prefixed alleles

    Returns a new dict with sanitized allele lists.
    """
    sanitized = {}
    for locus, alleles in genotype.items():
        if alleles is None or len(alleles) == 0:
            continue
        sanitized[locus] = [sanitize_allele(a, locus) for a in alleles]
    return sanitized


# ── Pydantic Schemas ───────────────────────────────────────────────


class ImputeRequest(BaseModel):
    """Patient unphased genotype. Each field is optional."""
    hla_a:      Optional[list] = None
    hla_c:      Optional[list] = None
    hla_b:      Optional[list] = None
    hla_drb345: Optional[list] = None
    hla_drb1:   Optional[list] = None
    hla_dqa1:   Optional[list] = None
    hla_dqb1:   Optional[list] = None
    hla_dpa1:   Optional[list] = None
    hla_dpb1:   Optional[list] = None


class PopulationFrequencies(BaseModel):
    Global: float = 0.0
    AFA:    float = 0.0
    ASI:    float = 0.0
    EUR:    float = 0.0
    HIS:    float = 0.0


class ImputedPair(BaseModel):
    rank: int
    haplotype_1: str = ""
    haplotype_2: str = ""
    posterior: float = 0.0
    cumulative: float = 0.0
    is_homozygous: bool = False
    population_frequencies: PopulationFrequencies = PopulationFrequencies()


class ImputeResponse(BaseModel):
    status: str
    population: str
    patient_genotype: dict
    total_possible_pairs: int
    entropy: float
    populations_available: list = []
    imputed_pairs: list[ImputedPair] = []


# ── Endpoints ──────────────────────────────────────────────────────


@app.get("/health")
def health_check():
    eng = get_engine()
    return {
        "status": "ok",
        "service": "HaploStats",
        "version": "1.0.0",
        "population": eng.population,
        "haplotypes_in_reference": len(eng._all_haplotypes) if eng._all_haplotypes else 0,
        "populations_available": POPULATION_ORDER,
        "database": "haplostats_normalized.db",
    }


@app.post("/impute", response_model=ImputeResponse)
def impute(request: ImputeRequest, population: str = "Global"):
    """
    Impute phased haplotype pairs from an unphased patient genotype.
    Returns top-ranked pairs with per-population diplotype frequencies.
    """
    if population not in POPULATION_ORDER:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown population '{population}'. "
                   f"Available: {', '.join(POPULATION_ORDER)}"
        )

    # Build raw genotype dict
    raw_genotype = {}
    for loc in ALL_LOCI:
        val = getattr(request, loc, None)
        if val is not None and len(val) > 0:
            raw_genotype[loc] = [str(a).strip() for a in val if a is not None]

    # Validate: at least 1 typed locus
    if sum(1 for v in raw_genotype.values() if v) == 0:
        raise HTTPException(
            status_code=400,
            detail="No typed loci provided. At least one locus is required."
        )

    # ── SANITIZE INPUT ──────────────────────────────────────────────
    clean_genotype = sanitize_genotype(raw_genotype)

    # Log sanitization differences for debugging
    for loc in clean_genotype:
        for i, (raw, clean) in enumerate(zip(
            raw_genotype.get(loc, []),
            clean_genotype[loc]
        )):
            if raw != clean:
                sys.stderr.write(f"  [Sanitize] {loc}[{i}]: {raw} → {clean}\n")

    # Run Bayesian engine
    engine = get_engine(population)
    result = engine.calculate_posterior(clean_genotype)

    if "error" in result or not result.get("pairs"):
        return ImputeResponse(
            status="error" if "error" in result else "success",
            population=population,
            patient_genotype={
                LOCUS_LABEL_MAP.get(k, k): v
                for k, v in clean_genotype.items()
            },
            total_possible_pairs=0,
            entropy=0.0,
            populations_available=POPULATION_ORDER,
            imputed_pairs=[],
        )

    # Format imputed pairs
    imputed_pairs = []
    for p in result["pairs"]:
        pf = p["population_frequencies"]
        imputed_pairs.append(ImputedPair(
            rank=p["rank"],
            haplotype_1=p["haplotype_1"],
            haplotype_2=p["haplotype_2"],
            posterior=p["posterior"],
            cumulative=p["cumulative"],
            is_homozygous=p["is_homozygous"],
            population_frequencies=PopulationFrequencies(
                Global=pf.get("Global", 0.0),
                AFA=pf.get("AFA", 0.0),
                ASI=pf.get("ASI", 0.0),
                EUR=pf.get("EUR", 0.0),
                HIS=pf.get("HIS", 0.0),
            ),
        ))

    return ImputeResponse(
        status="success",
        population=population,
        patient_genotype={
            LOCUS_LABEL_MAP.get(k, k): v
            for k, v in clean_genotype.items()
        },
        total_possible_pairs=result["total_possible_pairs"],
        entropy=result.get("entropy", 0.0),
        populations_available=POPULATION_ORDER,
        imputed_pairs=imputed_pairs,
    )


# ── Main ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    port = int(os.environ.get("HAPLOSTATS_PORT", 8000))
    host = os.environ.get("HAPLOSTATS_HOST", "0.0.0.0")
    print(f"🦞 HaploStats v1.0.0 on {host}:{port}")
    uvicorn.run(app, host=host, port=port, log_level="info")
