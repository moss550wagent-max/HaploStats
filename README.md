# 🧬 HaploStats — Clinical HLA Haplotype Imputation Engine

**HaploStats** is a clinical-grade web service that resolves unphased HLA genotypes into phased, high-resolution haplotype pairs. The live service runs a **Bayesian posterior calculator** (Hardy-Weinberg 2pq / p² diplotype frequencies) against a **normalized 2,026-row reference database** covering five population panels, ranked by a **heuristic tolerance score** (`match_percentage`, 0–100).

Because a 2,026-row reference cannot guarantee an exact 9-locus match for every new patient, the engine ships with **two search modes**:

- **Original Algorithm (strict)** — gatekeeper loci must match exactly; flexible loci are tolerance-scored.
- **Testing Mode (lenient)** — a one-click toggle for exploratory/QC work that drops all gatekeepers and returns any row matching **≥ 2 typed loci**, with the match percentage recomputed as *covered loci ÷ typed loci*.

Built for transplant immunogenetics, registry operations, and research applications where rapid, reference-driven haplotype imputation is required.

---

## Architecture

```
                    ┌─────────────────────────────────┐
                    │     Web Dashboard (static/)      │
                    │   (GET /)  single-panel form     │
                    │   Testing Mode toggle + green    │
                    │   matched-allele highlighting    │
                    └──────────┬──────────────────────┘
                               │
                    ┌──────────▼──────────────────────┐
                    │      FastAPI REST Service        │
                    │   POST /impute   GET /health     │
                    │   testing_mode flag supported    │
                    └──────────┬──────────────────────┘
                               │
                    ┌──────────▼──────────────────────┐
                    │   Bayesian Calculator (live)     │
                    │   scripts/bayesian_calc.py       │
                    │   • strict gatekeeper search     │
                    │   • lenient Testing Mode (≥2     │
                    │     matched loci)                │
                    │   • tolerance scoring 0–100      │
                    └──────────┬──────────────────────┘
                               │
                    ┌──────────▼──────────────────────┐
                    │  SQLite Reference DB             │
                    │  db/haplostats_normalized.db     │
                    │  2,026 haplotypes · 5 panels     │
                    │  Global / AFA / ASI / EUR / HIS  │
                    └─────────────────────────────────┘
```

> **Research modules** (kept in the repo, not part of the live web API):
> `scripts/em_algorithm.py` (Expectation-Maximisation solver), `scripts/dp_imputation.py` (Dirichlet Process + DPB1-decoupled model), and `scripts/match_grader.py` (NMDP HapLogic-style A/P/M grading).

### Core Technologies

| Component | Technology | Role |
|-----------|------------|------|
| **API Layer** | [FastAPI](https://fastapi.tiangolo.com/) | Web framework with auto-generated OpenAPI docs |
| **Server** | [Uvicorn](https://www.uvicorn.org/) | ASGI server, production-ready |
| **Database** | SQLite 3 | Normalized reference haplotype frequency DB (2,026 rows) |
| **Imputation** | Bayesian posterior (2pq / p²) | Phased pair ranking from unphased genotypes under HWE |
| **Scoring** | Heuristic tolerance (0–100) | Exact = 20 pts, allele-group = 16 pts per flexible locus |
| **Search Modes** | Strict gatekeeper / Lenient Testing Mode | ≥ 2 matched loci minimum in Testing Mode |
| **Frontend** | Static HTML + CSS + JS | Imputation form, ranked results table, green match highlighting |

---

## Search Modes & Match Scoring

### Original Algorithm (`testing_mode: false` — default)

The strict clinical path, matching the tolerance-scoring design:

- **Gatekeeper loci** — `HLA-A`, `HLA-C`, `HLA-B`, `HLA-DRB345` — must **match exactly**. Rows failing any gatekeeper are discarded outright.
- **Flexible loci** — `DRB1`, `DQA1`, `DQB1`, `DPA1`, `DPB1` — contribute 0–20 points each to a 100-point `match_percentage`:
  - **20 pts** exact match (`05:01` == `05:01`)
  - **16 pts** allele-group match (`05:01` vs `05:02` — field 1 identical)
  - **0 pts** mismatch
- Untyped flexible loci are not penalized (full credit), so a partially-typed patient can still reach 100 % on what they provided.
- Pairs are ranked by `match_percentage` (descending), then by population frequency (descending) as tie-breaker.

### Testing Mode (`testing_mode: true` — lenient)

A deliberate relaxation for exploratory / QC work on the small reference DB, where strict gatekeeper matching often returns **0 results**:

- **All gatekeeper requirements are removed.** A database row (haplotype) qualifies when **ANY two or more typed locus categories** carry a patient allele (e.g. HLA-A + HLA-B, or DRB1 + DQB1).
- **Dynamic scoring:** `match_percentage` = *(loci covered by the pair) ÷ (loci the patient typed)* × 100. Example: patient types 4 loci, pair covers 2 → **50 %**.
- **Performance guardrail:** the lenient scan can admit ~1,200 rows (~718k raw pairs) on the 2,026-row DB. Testing Mode therefore bounds the search to the **top-250 highest-frequency qualifying rows** and returns the **top-300 ranked pairs** (posteriors are still normalized over the full bounded pair space, so totals reported in the UI remain exact). Measured end-to-end: **~0.4 s**.
- In the dashboard, matched alleles (exact **or** allele-group) are highlighted in bold green (`#28a745`); the copy button always copies the clean raw haplotype string via `innerText`, ignoring highlight markup.

---

## Clinical Match Grades

The repository ships an NMDP HapLogic-style grading module (`scripts/match_grader.py`, standalone — not exposed by the live web API):

| Grade | Meaning | Criterion |
|-------|---------|-----------|
| **A** | Allele Match | Both patient alleles match both donor alleles at high-resolution |
| **P** | Potential Match | One allele matches, or posterior confidence < 0.5 for one pair |
| **M** | Mismatch | Neither allele matches |

**Scoring panels:**
- **Core (5/5):** HLA-A, -C, -B, DRB1, DQB1
- **Extended (6/6):** Core + DPB1
- **Full (9/9):** All loci including DRB345, DQA1, DPA1

The overall grade is reported as `n/N` (e.g., 8/10) for transparent clinical decision support.

---

## API Endpoints

### `GET /`
Serves the clinical web dashboard — a single-panel form for entering a patient's 9-locus HLA genotype, with the **Testing Mode (Lenient Match)** toggle next to the Calculate button and a ranked results table featuring match-percentage coloring and green matched-allele highlighting.

### `GET /health`
Returns service status, engine version, and reference database size.

```json
{
  "status": "ok",
  "service": "HaploStats",
  "version": "1.0.0",
  "population": "Global",
  "haplotypes_in_reference": 2026,
  "populations_available": ["Global", "AFA", "ASI", "EUR", "HIS"],
  "database": "haplostats_normalized.db"
}
```

### `POST /impute`
Resolves an unphased patient genotype into phased haplotype pairs, ranked by match percentage then population frequency.

**Request:** all nine loci are optional; each value is an `[allele1, allele2]` list (single-allele entries may be sent as `[allele, ""]` for hemizygous/null handling, e.g. DRB345). Alleles may be typed with or without the gene prefix (`"01:01"` → `"A*01:01"`).

```json
{
  "hla_a":      ["A*01:01", "A*03:01"],
  "hla_c":      ["C*07:01", "C*07:02"],
  "hla_b":      ["B*08:01", "B*07:02"],
  "hla_drb345": ["DRB3*01:01", "DRB5*01:01"],
  "hla_drb1":   ["DRB1*03:01", "DRB1*15:01"],
  "hla_dqa1":   ["DQA1*05:01", "DQA1*01:02"],
  "hla_dqb1":   ["DQB1*02:01", "DQB1*06:02"],
  "hla_dpa1":   ["DPA1*01:03", "DPA1*01:03"],
  "hla_dpb1":   ["DPB1*04:01", "DPB1*04:01"],
  "testing_mode": false
}
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `hla_a` … `hla_dpb1` | `list[str] \| null` | `null` | Patient alleles at each locus (optional) |
| `testing_mode` | `bool` | `false` | `false` → strict gatekeeper algorithm; `true` → lenient ≥ 2-locus search |
| `population` | query param | `Global` | Ranking population panel (`Global`, `AFA`, `ASI`, `EUR`, `HIS`) |

**Response:** ranked phased pairs with per-locus match scoring and per-population diplotype frequencies.

```json
{
  "status": "success",
  "population": "Global",
  "testing_mode": false,
  "total_possible_pairs": 302,
  "entropy": 3.214,
  "populations_available": ["Global", "AFA", "ASI", "EUR", "HIS"],
  "imputed_pairs": [
    {
      "rank": 1,
      "haplotype_1": "A*01:01 ~ C*07:01 ~ B*08:01 ~ DRB1*03:01 ~ DRB3*01:01 ~ DQA1*05:01 ~ DQB1*02:01 ~ DPA1*02:01 ~ DPB1*01:01",
      "haplotype_2": "A*03:01 ~ C*04:01 ~ B*35:01 ~ DRB1*01:01 ~ DQA1*01:01 ~ DQB1*05:01 ~ DPA1*01:03 ~ DPB1*04:02",
      "posterior": 0.1842,
      "cumulative": 0.1842,
      "match_percentage": 100,
      "is_homozygous": false,
      "population_frequencies": { "Global": 0.00017, "AFA": 0.00031, "ASI": 0.0, "EUR": 0.00012, "HIS": 0.00005 }
    }
  ]
}
```

When `testing_mode: true`, `match_percentage` reflects the *covered ÷ typed* locus ratio and the returned list is capped at the top 300 pairs (the `total_possible_pairs` field always reports the full bounded count).

### `POST /compare` *(legacy — not exposed by the live API)*
A patient-vs-donor match comparison exists in the research modules (`match_grader.py`, `em_algorithm.py`) but is **not wired into the current FastAPI app**. The live service is the single-patient `/impute` pipeline.

---

## Quick Start — Local Development

### Prerequisites
- Python 3.9+
- pip

### Setup
```bash
# Clone the repository
git clone https://github.com/moss550wagent-max/HaploStats.git
cd HaploStats

# Create a virtual environment
python -m venv venv
source venv/bin/activate   # Linux/Mac
venv\Scripts\activate      # Windows

# Install dependencies
pip install -r requirements.txt

# Start the server
uvicorn scripts.api:app --host 127.0.0.1 --port 8000 --reload
```

Open **http://localhost:8000** in your browser — the dashboard loads automatically.

API documentation (auto-generated by FastAPI) is available at **http://localhost:8000/docs**.

---

## Deployment to Render.com

HaploStats ships with a production-ready Dockerfile for cloud deployment.

### One-Click Deploy

1. Push the repository to GitHub
2. Log in to [Render.com](https://render.com)
3. Click **New +** → **Blueprint** or **Web Service**
4. Connect your GitHub repo
5. Render auto-detects the `Dockerfile` — no additional configuration needed
6. Set the **Service Name** (e.g., `haplostats`) and click **Create Web Service**

Your app will be live at `https://haplostats.onrender.com`.

### Environment Variables (Optional)

| Variable | Default | Purpose |
|----------|---------|---------|
| `HAPLOSTATS_PORT` | `8000` | HTTP port override |
| `HAPLOSTATS_HOST` | `0.0.0.0` | Bind address override |

---

## Project Structure

```
HaploStats/
├── Dockerfile              # Production container image
├── requirements.txt        # Python dependency manifest
├── .dockerignore           # Build context exclusions
├── README.md               # This file
├── db/
│   ├── haplostats_normalized.db  # Live reference DB (2,026 haplotypes, 5 panels)
│   └── haplostats.db             # Legacy 577-haplotype DB
├── scripts/
│   ├── api.py              # FastAPI web service (entry point)
│   ├── bayesian_calc.py    # Live Bayesian engine: strict + Testing Mode
│   ├── em_algorithm.py     # Core Expectation-Maximisation solver (research)
│   ├── dp_imputation.py    # Dirichlet Process + DPB1-decoupled model (research)
│   ├── match_grader.py     # Clinical match grade assignment A/P/M (standalone)
│   ├── clean_new_db.py     # Normalized DB build/cleanup pipeline
│   ├── ingest_data.py      # Reference data ingestion pipeline
│   ├── batch_benchmarker.py# Concordance and performance benchmarking
│   ├── generate_synthetic_patients.py  # Synthetic test data generator
│   ├── fix_indexes.py      # Database index maintenance
│   ├── verify_db.py        # Reference database integrity checks
│   ├── test_api.py         # API smoke tests
│   └── test_compare.py     # Compare endpoint test suite
└── static/
    └── index.html          # Clinical web dashboard (toggle, highlighting, copy)
```

---

## Reference Database

The live engine reads `db/haplostats_normalized.db` — a normalized 2-field (high-resolution) reference built from the original population tables. All **2,026 haplotypes** carry frequency data for every panel:

| Population | Code | Rows |
|------------|------|------|
| **Global (meta)** | `Global` | 2,026 |
| **African American** | `AFA` | 2,026 |
| **Asian** | `ASI` | 2,026 |
| **European American** | `EUR` | 2,026 |
| **Hispanic** | `HIS` | 2,026 |

The database is read-only at runtime. Additional population tables can be ingested via `scripts/ingest_data.py` and re-normalized with `scripts/clean_new_db.py`.

---

## Benchmarking

A built-in benchmarking suite validates the imputation engine against a synthetic truth set:

```bash
python scripts/batch_benchmarker.py
```

Reports are generated in `data/results/` and include per-locus concordance rates, block-level convergence statistics, and posterior probability distributions.

---

## License

*To be determined.* Clinical deployment should ensure compliance with applicable regulations and data use agreements for the underlying haplotype reference data.

---

## Acknowledgements

- Inspired by the NMDP HapLogic matching algorithm
- Reference frequencies derived from the 11th International Histocompatibility Workshop population data
- Built with [FastAPI](https://fastapi.tiangolo.com/) and [SQLite](https://www.sqlite.org/)

---

<p align="center"><strong>HaploStats</strong> — <em>bringing haplotype resolution to the bedside.</em></p>
