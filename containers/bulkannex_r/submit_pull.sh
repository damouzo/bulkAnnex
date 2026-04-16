#!/usr/bin/env bash
#SBATCH --job-name=pull_bulkannex_sif
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --output=pull_bulkannex_sif_%j.log

# submit_pull.sh — Pull the bulkAnnex container SIF on a compute node.
#
# Run once before your first pipeline execution:
#
#   cd /data/BCI-KRP/projects/bulkAnnex
#   sbatch containers/bulkannex_r/submit_pull.sh
#
# When the job finishes, check the log and then run the pipeline normally.
# The singularity profile in nextflow.config auto-detects the SIF and uses it
# instead of pulling from DockerHub on every run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF="${SCRIPT_DIR}/bulkannex_r_1.0.0.sif"
IMAGE="docker://damouzo/bulkannex_r:1.0.0"

# Use a scratch tmp dir — /tmp is too small for mksquashfs
SCRATCH_TMP="/gpfs/scratch/${USER}/singularity_tmp"
mkdir -p "${SCRATCH_TMP}"
export SINGULARITY_TMPDIR="${SCRATCH_TMP}"
export APPTAINER_TMPDIR="${SCRATCH_TMP}"

echo "========================================"
echo "  bulkAnnex — Singularity SIF pull"
echo "  Node   : $(hostname)"
echo "  Image  : ${IMAGE}"
echo "  Target : ${SIF}"
echo "  TMP    : ${SCRATCH_TMP}"
echo "========================================"
echo ""

singularity pull --force "${SIF}" "${IMAGE}"

echo ""
echo "Done. SIF written to: ${SIF}"
echo "You can now run the pipeline with -profile apocrita,singularity"
