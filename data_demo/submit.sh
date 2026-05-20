#!/usr/bin/env bash
#SBATCH --job-name=bulkannex-demo
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=data_demo/logs/bulkannex_demo_%j.out
#SBATCH --error=data_demo/logs/bulkannex_demo_%j.err
#
# SLURM submission script for the bulkAnnex demo dataset.
# Run from the project root:
#   sbatch data_demo/submit.sh
#
# Prerequisites (run once from an interactive node or login node):
#   python data_demo/generate_demo.py

set -euo pipefail

# SLURM copies the script to a temp location before executing it, so
# BASH_SOURCE[0] resolves to /var/spool/slurmd/. Use SLURM_SUBMIT_DIR
# instead, which is set to the directory where sbatch was called.
# This requires submitting from the project root (see usage above).
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

WORKDIR="/gpfs/scratch/${USER}/bulkannex_demo_work"
COUNTS_FILE="${PROJECT_DIR}/data_demo/salmon.merged.gene_counts.tsv"

# ── Environment ───────────────────────────────────────────────────────────────
module unload openjdk 2>/dev/null || true
module load nextflow

# Run Nextflow from WORKDIR so .nextflow/ state files land on scratch (writable)
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# ── Run ───────────────────────────────────────────────────────────────────────
nextflow run "${PROJECT_DIR}" \
    -profile singularity,apocrita \
    -w "${WORKDIR}" \
    --input      "${PROJECT_DIR}/data_demo/samplesheet.csv" \
    --counts     "${COUNTS_FILE}" \
    --contrasts  "${PROJECT_DIR}/data_demo/contrasts.csv" \
    --organism   human \
    --outdir     "${PROJECT_DIR}/results_demo" \
    --run_gsea   true \
    -resume
