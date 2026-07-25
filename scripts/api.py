#!/usr/bin/env python3
"""
HaploStats — FastAPI Web Service (Population Explorer)
Phase 11+: Population-aware REST API returning multi-population
diplotype frequencies (2pq / p²) for every imputed haplotype pair.

Endpoints:
  POST /impute  — submit unphased patient genotype → ranked phased haplotypes
                   with per-population diplotype frequencies (AFA, API, CAU,
                   HIS, Global, and more)
  GET  /health  — service health check
  GET  /        — population explorer web dashboard
"""

import sys
import os
from pathlib import Path

# Ensure HaploStats project root is on the path
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

from scripts.bayesian_calc import HaploMath, POPULATION_ORDER

# ── App Initialisation ─────────────────────────────────────────────

app = FastAPI(
    title="HaploStats — HLA Haplotype Population Explorer",
    description="Clinical-grade Bayesian haplotype imputation engine "
                "returning per-population diplotype frequencies (2pq / p²) "
                "for all reference populations (AFA, API, CAU, HIS, Global).",
    version="0.2.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files (web dashboard)
STATIC_DIR = HERE.parent / "static"
if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.get("/")
def serve_dashboard():
    """Serve the population explorer web dashboard at the root URL."""
    idx = STATIC_DIR / "index.html"
    if idx.exists():
        return FileResponse(str(idx))
    return {
        "status": "ok",
        "message": "HaploStats Population Explorer API running. "
                   "Dashboard static files not found."
    }


# ── Locus constants ────────────────────────────────────────────────

LOCUS_LABEL = {
    'hla_a': 'HLA-A', 'hla_c': 'HLA-C', 'hla_b': 'HLA-B',
    'hla_drb345': 'HLA-DRB345', 'hla_drb1': 'HLA-DRB1',
    'hla_dqa1': 'HLA-DQA1', 'hla_dqb1': 'HLA-DQB1',
    'hla_dpa1': 'HLA-DPA1', 'hla_dpb1': 'HLA-DPB1',
}

ALL_LOCI = list(LOCUS_LABEL.keys())


# ── Global engine instance (lazy-init on first request) ────────────

_engine: Optional[HaploMath] = None


def get_engine(population: str = "Global") -> HaploMath:
    global _engine
    if _engine is None:
        _engine = HaploMath(population=population)
        _engine.connect()
    return _engine


# ── Pydantic Schemas ───────────────────────────────────────────────


class ImputeRequest(BaseModel):
    """
    Patient unphased genotype. Each field is optional.
    Missing loci are treated as unknown.

    Examples:
      {"hla_a": ["01:01", "02:01"],
       "hla_b": ["08:01", "44:02"],
       "hla_drb1": ["03:01", "04:01"]}
    """
    hla_a: Optional[list] = None
    hla_c: Optional[list] = None
    hla_b: Optional[list] = None
    hla_drb345: Optional[list] = None
    hla_drb1: Optional[list] = None
    hla_dqa1: Optional[list] = None
    hla_dqb1: Optional[list] = None
    hla_dpa1: Optional[list] = None
    hla_dpb1: Optional[list] = None


class PopulationFrequencies(BaseModel):
    """
    Diplotype frequencies (Hardy-Weinberg 2pq / p²) across all
    reference populations for a single haplotype pair.
    """
    Global: float = 0.0
    AFA: float = 0.0
    API: float = 0.0
    CAU: float = 0.0
    HIS: float = 0.0
    European: float = 0.0
    Spanish: float = 0.0
    Mexican: float = 0.0
    Arab: float = 0.0


class ImputedPair(BaseModel):
    """A single phased haplotype pair with all-population frequencies."""
    rank: int
    haplotype_1: str = ""
    haplotype_2: str = ""
    posterior: float = 0.0
    cumulative: float = 0.0
    is_homozygous: bool = False
    population_frequencies: PopulationFrequencies = PopulationFrequencies()


class ImputeResponse(BaseModel):
    """Top-level API response with multi-population frequency data."""
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
    """Return service status and reference database summary."""
    eng = get_engine()
    return {
        "status": "ok",
        "service": "HaploStats",
        "version": "0.2.0",
        "population": eng.population,
        "haplotypes_in_reference": len(eng._all_haplotypes) if eng._all_haplotypes else 0,
        "populations_available": POPULATION_ORDER,
    }


@app.post("/impute", response_model=ImputeResponse)
def impute(request: ImputeRequest, population: str = "Global"):
    """
    Impute phased haplotype pairs from an unphased patient genotype.

    Returns the top-ranked pairs with per-population diplotype
    frequencies (2pq / p²) for all reference populations (AFA, API,
    CAU, HIS, Global, European, Spanish, Mexican, Arab).

    - Missing loci are treated as unconstrained
    - Population parameter selects which population's posterior
      is used for ranking
    """
    # Validate population parameter
    if population not in POPULATION_ORDER:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown population '{population}'. "
                   f"Available: {', '.join(POPULATION_ORDER)}"
        )

    # Build genotype dict from the pydantic model
    patient_genotype = {}
    for loc in ALL_LOCI:
        val = getattr(request, loc, None)
        if val is not None and len(val) > 0:
            patient_genotype[loc] = [str(a).strip() for a in val if a is not None]
        else:
            patient_genotype[loc] = None

    # Validate: at least 1 typed locus
    typed_count = sum(1 for v in patient_genotype.values() if v is not None)
    if typed_count == 0:
        raise HTTPException(
            status_code=400,
            detail="No typed loci provided. At least one locus is required."
        )

    # Remove None entries before passing to engine
    clean_genotype = {
        loc: alleles
        for loc, alleles in patient_genotype.items()
        if alleles is not None
    }

    # Run Bayesian engine
    engine = get_engine(population)
    result = engine.calculate_posterior(clean_genotype)

    # Handle empty result
    if 'error' in result:
        return ImputeResponse(
            status="error",
            population=population,
            patient_genotype={
                LOCUS_LABEL.get(k, k): v
                for k, v in clean_genotype.items()
            },
            total_possible_pairs=0,
            entropy=0.0,
            populations_available=POPULATION_ORDER,
            imputed_pairs=[],
        )

    # Format imputed pairs with full population frequency dictionaries
    imputed_pairs = []
    for p in result.get('pairs', []):
        # Build population frequencies dict (all values, not just non-zero)
        pop_freqs = {}
        for pop_name in POPULATION_ORDER:
            pop_freqs[pop_name] = p['population_frequencies'].get(pop_name, 0.0)

        imputed_pairs.append(ImputedPair(
            rank=p['rank'],
            haplotype_1=p['haplotype_1'],
            haplotype_2=p['haplotype_2'],
            posterior=p['posterior'],
            cumulative=p['cumulative'],
            is_homozygous=p.get('is_homozygous', False),
            population_frequencies=PopulationFrequencies(**pop_freqs),
        ))

    return ImputeResponse(
        status="success",
        population=population,
        patient_genotype={
            LOCUS_LABEL.get(k, k): v
            for k, v in clean_genotype.items()
        },
        total_possible_pairs=result.get('total_possible_pairs', 0),
        entropy=result.get('entropy', 0.0),
        populations_available=POPULATION_ORDER,
        imputed_pairs=imputed_pairs,
    )


# ── Main (for direct execution) ────────────────────────────────────

if __name__ == "__main__":
    port = int(os.environ.get("HAPLOSTATS_PORT", 8000))
    host = os.environ.get("HAPLOSTATS_HOST", "0.0.0.0")
    print(f"🦞 HaploStats Population Explorer on {host}:{port}")
    uvicorn.run(app, host=host, port=port, log_level="info")
