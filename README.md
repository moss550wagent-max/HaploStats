# 🧬 HaploStats — Clinical HLA Haplotype Imputation Engine

**HaploStats** is a clinical-grade web service that resolves unphased HLA genotypes into phased, high-resolution haplotype pairs using a Dirichlet Process-enhanced Expectation-Maximisation (EM) algorithm with a decoupled DPB1 conditional probability model. It assigns NMDP HapLogic-style match grades (A / P / M) for donor-recipient compatibility assessment.

Built for transplant immunogenetics, registry operations, and research applications where rapid, reference-driven haplotype imputation is required.

---

## Architecture

```
                    ┌─────────────────────────────────┐
                    │     Web Dashboard (static/)      │
                    │   (GET /)                        │
                    └──────────┬──────────────────────┘
                               │
                    ┌──────────▼──────────────────────┐
                    │      FastAPI REST Service        │
                    │   POST /impute  POST /compare    │
                    │   GET  /health                   │
                    └──────────┬──────────────────────┘
                               │
          ┌────────────────────┼────────────────────┐
          ▼                    ▼                    ▼
   ┌─────────────┐    ┌──────────────┐    ┌────────────────┐
   │  EM Engine   │    │  Bayesian    │    │  Match Grader   │
   │  (Phase 4)   │    │  Calculator  │    │  (Phase 9)     │
   │  + DPB1      │    │  (Phase 6)   │    │  A/P/M per     │
   │  decoupled   │    │              │    │  locus +       │
   │  model       │    │              │    │  overall 10/10 │
   └──────┬───────┘    └──────┬───────┘    └────────────────┘
          └──────────┬────────┘
                     ▼
           ┌─────────────────┐
           │  SQLite Ref DB   │
           │  577 haplotypes  │
           │  AfAm / EuAm /   │
           │  Global panels   │
           └─────────────────┘
```

### Core Technologies

| Component | Technology | Role |
|-----------|------------|------|
| **API Layer** | [FastAPI](https://fastapi.tiangolo.com/) | High-performance async web framework with auto-generated OpenAPI docs |
| **Server** | [Uvicorn](https://www.uvicorn.org/) | ASGI server, production-ready |
| **Database** | SQLite 3 | Embedded reference haplotype frequency database |
| **Imputation** | Expectation-Maximisation + Dirichlet Process | Iterative phase resolution from unphased genotypes |
| **DPB1 Model** | Decoupled conditional probability | Addresses the DQ–DP recombination hotspot using P(DPB1 \| DPA1) marginals |
| **Match Grading** | NMDP HapLogic standards | Per-locus A / P / M + x/10 overall score |
| **Frontend** | Static HTML + CSS + JS | Two-panel clinical dashboard (patient vs. donor) |

### Imputation Pipeline

1. **Block Partitioning** — The 9-locus HLA space is split into three blocks:
   - **Block 1 — Core** (HLA-A, -C, -B): the classical Class I chain
   - **Block 2 — DR/DQ** (DRB1, DQB1, DRB345, DQA1): Class II core
   - **Block 3 — DP** (DPA1, DPB1): decoupled to handle the DQ–DP recombination hotspot independently
2. **Progressive EM** — Each block is solved sequentially; posterior probabilities propagate forward as priors for subsequent blocks
3. **Haplotype Scoring** — Full 9-locus haplotype pairs are ranked by posterior probability with cumulative coverage reporting
4. **Match Grading** — Two independent imputation runs (patient, donor) are compared locus-by-locus against HapLogic-style match criteria

---

## Clinical Match Grades

HaploStats assigns grades per locus following the NMDP HapLogic convention:

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
Serves the clinical web dashboard — a two-panel interface for side-by-side patient vs. donor analysis.

### `GET /health`
Returns service status and reference database size.

```json
{
  "status": "ok",
  "service": "HaploStats",
  "version": "0.1.0",
  "haplotypes_in_reference": 577
}
```

### `POST /impute`
Resolves an unphased patient genotype into phased haplotype pairs.

**Request:**
```json
{
  "hla_a": ["01:01", "02:01"],
  "hla_b": ["08:01", "44:02"],
  "hla_drb1": ["03:01", "04:01"]
}
```

**Response:** Top 3 phased haplotype pairs ranked by posterior probability, with block-level convergence diagnostics and entropy.

### `POST /compare`
Compares a patient and a donor for HLA matching. Both genotypes are independently imputed, then the top-ranked pairs are compared per-locus.

**Request:**
```json
{
  "patient": { "hla_a": ["02:01", "03:01"], "hla_b": ["07:02", "44:02"] },
  "donor": { "hla_a": ["02:01", "24:02"], "hla_b": ["07:02", "44:02"] },
  "population": "Global",
  "grades": "core"
}
```

**Response:** Per-locus A/P/M grades + overall compatibility score.

---

## Quick Start — Local Development

### Prerequisites
- Python 3.9+
- pip

### Setup
```bash
# Clone the repository
git clone https://github.com/your-org/haplostats.git
cd haplostats

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
haplostats/
├── Dockerfile              # Production container image
├── requirements.txt        # Python dependency manifest
├── .dockerignore           # Build context exclusions
├── README.md               # This file
├── db/
│   └── haplostats.db       # SQLite reference database (577 haplotypes)
├── scripts/
│   ├── api.py              # FastAPI web service (entry point)
│   ├── dp_imputation.py    # Dirichlet Process + EM engine
│   ├── em_algorithm.py     # Core Expectation-Maximisation solver
│   ├── bayesian_calc.py    # Bayesian posterior calculator
│   ├── match_grader.py     # Clinical match grade assignment (A/P/M)
│   ├── ingest_data.py      # Reference data ingestion pipeline
│   ├── batch_benchmarker.py# Concordance and performance benchmarking
│   ├── generate_synthetic_patients.py  # Synthetic test data generator
│   ├── fix_indexes.py      # Database index maintenance
│   ├── verify_db.py        # Reference database integrity checks
│   ├── test_api.py         # API smoke tests
│   └── test_compare.py     # Compare endpoint test suite
└── static/
    └── index.html          # Clinical web dashboard
```

---

## Reference Database

The built-in reference database (`db/haplostats.db`) contains allele and haplotype frequency data derived from publicly available population tables. Populations currently available:

| Population | Haplotypes | Source Panel |
|------------|-----------|--------------|
| **African American** | ~190 | NMDP / 11th IHW |
| **European American** | ~190 | NMDP / 11th IHW |
| **Global (meta)** | ~577 | Combined |

The database is read-only at runtime. Additional population tables can be ingested via `scripts/ingest_data.py`.

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
