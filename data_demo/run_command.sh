#!/usr/bin/env bash
#
# Run the bulkAnnex demo dataset locally.
# Run from the project root:
#   bash data_demo/run_command.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
COUNTS_FILE="${SCRIPT_DIR}/salmon.merged.gene_counts.tsv"

# Generate counts matrix if not present
if [ ! -f "${COUNTS_FILE}" ]; then
    echo "Generating demo counts matrix..."
    python "${SCRIPT_DIR}/generate_demo.py"
fi

nextflow run "${PROJECT_DIR}" \
    -profile conda \
    --input      "${SCRIPT_DIR}/samplesheet.csv" \
    --counts     "${COUNTS_FILE}" \
    --contrasts  "${SCRIPT_DIR}/contrasts.csv" \
    --organism   human \
    --outdir     "${PROJECT_DIR}/results_demo" \
    --run_gsea   true \
    --max_memory '8.GB' \
    -resume
