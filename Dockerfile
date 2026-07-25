# ── HaploStats — Production Dockerfile ─────────────────────────────
# Phase 11: Cloud deployment (Render.com, Railway, Fly.io, etc.)

FROM python:3.9-slim

LABEL org.opencontainers.image.title="HaploStats"
LABEL org.opencontainers.image.description="Clinical HLA Haplotype Imputation Engine"
LABEL org.opencontainers.image.version="0.1.0"

# Prevent Python from writing .pyc files & buffer stdout
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Copy dependency manifest first (leveraging Docker layer cache)
COPY requirements.txt .

# Install system dependencies + Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY db/        db/
COPY scripts/   scripts/
COPY static/    static/

# Expose the API port
EXPOSE 8000

# Run with Uvicorn (binding to 0.0.0.0 so Render/Railway ingress works)
CMD ["uvicorn", "scripts.api:app", "--host", "0.0.0.0", "--port", "8000"]
