#!/usr/bin/env python3
"""
HaploStats — FastAPI Web Service
Phase 5: Clinical-grade API for haplotype imputation queries.

Endpoints:
  POST /impute  — submit unphased patient genotype → ranked phased haplotypes
  GET  /health  — service health check
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

from scripts.dp_imputation import HaploEM, ALL_LOCI, LOCUS_LABEL
from scripts.match_grader import MatchGrader, CORE_MATCH_LOCL, EXTENDED_MATCH_LOCL, FULL_MATCH_LOCL

# ── App Initialisation ─────────────────────────────────────────────

app = FastAPI(
    title="HaploStats — HLA Haplotype Imputation Engine",
    description="Clinical-grade Bayesian/EM engine for resolving unphased "
                "HLA genotypes into phased high-resolution haplotype pairs.",
    version="0.1.0",
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
    """Serve the web dashboard at the root URL."""
    idx = STATIC_DIR / "index.html"
    if idx.exists():
        return FileResponse(str(idx))
    return {"status": "ok", "message": "HaploStats API running. Dashboard static files not found."}

# Global engine instance (lazy-init on first request)
_engine: Optional[HaploEM] = None


def get_engine(population: str = "Global") -> HaploEM:
    global _engine
    if _engine is None:
        _engine = HaploEM(population=population)
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


class TopPair(BaseModel):
    rank: int
    haplotype_1: str = ""
    haplotype_2: str = ""
    posterior: float = 0.0
    cumulative: float = 0.0
    h1_frequency: float = 0.0
    h2_frequency: float = 0.0


class BlockInfo(BaseModel):
    block: str
    haplotypes: int = 0
    pairs_before_em: int = 0
    pairs_after_trim: int = 0
    converged_iterations: int = 0


class ImputeResponse(BaseModel):
    """Top-level API response."""
    status: str
    population: str
    patient_genotype: dict
    total_possible_pairs: int
    entropy: float
    blocks: list[BlockInfo] = []
    top_pairs: list[TopPair] = []


class CompareRequest(BaseModel):
    """Donor-recipient comparison request."""
    patient: ImputeRequest
    donor: ImputeRequest
    population: str = "Global"
    grades: str = "core"  # "core" (A,B,C,DRB1,DQB1), "extended" (+DPB1), "full" (all 9)


class LocusGrade(BaseModel):
    grade: str
    patient_alleles: list = []
    donor_alleles: list = []
    patient_posterior: float = 0.0
    explanation: str = ""


class MatchOverall(BaseModel):
    overall_grade: str
    compatibility: str
    total_loci_scored: int = 0
    allele_matches: int = 0
    potential_matches: int = 0
    mismatches: int = 0


class CompareResponse(BaseModel):
    status: str
    population: str
    match_overall: MatchOverall
    locus_grades: dict = {}
    patient_pairs: list[TopPair] = []
    donor_pairs: list[TopPair] = []


# ── Init match grader (singleton) ──────────────────────────────────

_match_grader: Optional[MatchGrader] = None


def get_grader() -> MatchGrader:
    global _match_grader
    if _match_grader is None:
        _match_grader = MatchGrader()
    return _match_grader


# ── Endpoints ──────────────────────────────────────────────────────

@app.get("/health")
def health_check():
    eng = get_engine()
    return {
        "status": "ok",
        "service": "HaploStats",
        "version": "0.1.0",
        "haplotypes_in_reference": len(eng.full_ref) if eng.full_ref else 577,
    }


@app.post("/impute", response_model=ImputeResponse)
def impute(request: ImputeRequest, population: str = "Global"):
    """
    Impute phased haplotype pairs from an unphased patient genotype.

    - Missing loci → treated as unknown, imputed with progressive EM
    - Population parameter selects reference frequency column
    """
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

    # Run engine
    engine = get_engine(population)
    result = engine.progressive_em(patient_genotype)

    # Format blocks
    blocks = []
    for b in result.get('blocks', []):
        blocks.append(BlockInfo(
            block=b.get('block', ''),
            haplotypes=b.get('haplotypes', 0),
            pairs_before_em=b.get('pairs_before_em', 0),
            pairs_after_trim=b.get('pairs_after_trim', 0),
            converged_iterations=b.get('converged_iterations', 0),
        ))

    # Format top pairs
    top_pairs = []
    for t3 in result.get('top_3', []):
        top_pairs.append(TopPair(
            rank=t3.get('rank', 0),
            haplotype_1=t3.get('haplotype_1', ''),
            haplotype_2=t3.get('haplotype_2', ''),
            posterior=t3.get('posterior', 0.0),
            cumulative=t3.get('cumulative', 0.0),
            h1_frequency=t3.get('h1_frequency', 0.0),
            h2_frequency=t3.get('h2_frequency', 0.0),
        ))

    return ImputeResponse(
        status="success",
        population=result.get('population', population),
        patient_genotype={
            LOCUS_LABEL.get(k, k): v
            for k, v in patient_genotype.items()
            if v is not None
        },
        total_possible_pairs=result.get('total_pairs_final', 0),
        entropy=result.get('entropy', 0.0),
        blocks=blocks,
        top_pairs=top_pairs,
    )


@app.post("/compare", response_model=CompareResponse)
def compare(request: CompareRequest):
    """
    Compare a patient and a donor for HLA matching.

    Both genotypes are imputed independently, then the top-ranked
    haplotype pairs are compared locus-by-locus and assigned
    clinical match grades (A = allele match, P = potential, M = mismatch).
    """
    population = request.population or "Global"
    engine = get_engine(population)
    grader = get_grader()

    # Select loci set
    if request.grades == "extended":
        match_loci = EXTENDED_MATCH_LOCL
    elif request.grades == "full":
        match_loci = FULL_MATCH_LOCL
    else:
        match_loci = CORE_MATCH_LOCL

    # Helper: build genotype from ImputeRequest
    def build_genotype(req: ImputeRequest) -> dict:
        gt = {}
        for loc in ALL_LOCI:
            val = getattr(req, loc, None)
            if val is not None and len(val) > 0:
                gt[loc] = [str(a).strip() for a in val if a is not None]
            else:
                gt[loc] = None
        return gt

    # Impute patient
    p_gt = build_genotype(request.patient)
    typed_p = sum(1 for v in p_gt.values() if v is not None)
    if typed_p == 0:
        raise HTTPException(400, detail="Patient has no typed loci")
    p_result = engine.progressive_em(p_gt)

    # Impute donor
    d_gt = build_genotype(request.donor)
    typed_d = sum(1 for v in d_gt.values() if v is not None)
    if typed_d == 0:
        raise HTTPException(400, detail="Donor has no typed loci")
    d_result = engine.progressive_em(d_gt)

    # Format top pairs
    def fmt_pairs(result: dict) -> list[TopPair]:
        pairs = []
        for t3 in result.get('top_3', []):
            pairs.append(TopPair(
                rank=t3.get('rank', 0),
                haplotype_1=t3.get('haplotype_1', ''),
                haplotype_2=t3.get('haplotype_2', ''),
                posterior=t3.get('posterior', 0.0),
                cumulative=t3.get('cumulative', 0.0),
                h1_frequency=t3.get('h1_frequency', 0.0),
                h2_frequency=t3.get('h2_frequency', 0.0),
            ))
        return pairs

    p_top = fmt_pairs(p_result)
    d_top = fmt_pairs(d_result)

    # Grade (handle empty results gracefully)
    if not p_top:
        patient_unmatchable = True
    else:
        patient_unmatchable = False

    if not d_top:
        donor_unmatchable = True
    else:
        donor_unmatchable = False

    if patient_unmatchable or donor_unmatchable:
        locus_grades = {}
        for loc in match_loci:
            label = LOCUS_LABEL.get(loc, loc)
            explanation = "Patient has no matching haplotypes in reference" if patient_unmatchable else ""
            if donor_unmatchable:
                explanation = "Donor has no matching haplotypes in reference" if not explanation else "Both unmatchable"
            locus_grades[label] = LocusGrade(
                grade='?',
                patient_alleles=[],
                donor_alleles=[],
                explanation=explanation,
            )
        return CompareResponse(
            status="success",
            population=population,
            match_overall=MatchOverall(
                overall_grade="UNMATCHABLE",
                compatibility="No match possible — haplotype(s) outside reference",
                total_loci_scored=0,
                allele_matches=0,
                potential_matches=0,
                mismatches=0,
            ),
            locus_grades=locus_grades,
            patient_pairs=p_top,
            donor_pairs=d_top,
        )

    match_result = grader.compare_profiles(
        p_top[0].model_dump() if hasattr(p_top[0], 'model_dump') else p_top[0],
        d_top[0].model_dump() if hasattr(d_top[0], 'model_dump') else d_top[0],
        loci=match_loci,
    )

    # Build locus grades
    locus_grades = {}
    for label, lr in match_result['loci'].items():
        locus_grades[label] = LocusGrade(
            grade=lr['grade'],
            patient_alleles=lr['patient_alleles'],
            donor_alleles=lr['donor_alleles'],
            patient_posterior=lr.get('patient_posterior', 0),
            explanation=lr['explanation'],
        )

    overall = match_result['overall']

    return CompareResponse(
        status="success",
        population=population,
        match_overall=MatchOverall(
            overall_grade=overall['overall_grade'],
            compatibility=overall['compatibility'],
            total_loci_scored=overall['total_loci_scored'],
            allele_matches=overall['allele_matches'],
            potential_matches=overall['potential_matches'],
            mismatches=overall['mismatches'],
        ),
        locus_grades=locus_grades,
        patient_pairs=p_top,
        donor_pairs=d_top,
    )


# ── Main (for direct execution) ────────────────────────────────────

if __name__ == "__main__":
    port = int(os.environ.get("HAPLOSTATS_PORT", 8000))
    host = os.environ.get("HAPLOSTATS_HOST", "0.0.0.0")
    print(f"🦞 HaploStats API starting on {host}:{port}")
    uvicorn.run(app, host=host, port=port, log_level="info")
