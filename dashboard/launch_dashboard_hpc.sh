#!/usr/bin/env bash
# launch_dashboard_hpc.sh — launch bulkAnnex dashboard on Apocrita (QMUL SLURM HPC)
#
# Starts the Shiny app inside the bulkAnnex Singularity container via srun,
# then prints SSH tunnel instructions for BOTH VSCode Remote SSH and local terminal.
#
# Usage:
#   bash launch_dashboard_hpc.sh [results_dir] [port] [time] [mem] [sif]
#
#   results_dir  Path to Nextflow outdir  (default: parent of this script)
#   port         Port for the app; 0 = auto-detect a free port (default: 0)
#   time         SLURM time limit         (default: 04:00:00)
#   mem          SLURM memory request     (default: 8G)
#   sif          Path to Singularity SIF  (default: $BULKANNEX_SIF env var,
#                then auto-detect, then docker://damouzo/bulkannex_r:1.0.0)
#
# Example:
#   BULKANNEX_SIF=/path/to/bulkannex_r_1.0.0.sif \
#     bash launch_dashboard_hpc.sh demo_results

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-$(dirname "${SCRIPT_DIR}")}"
REQUESTED_PORT="${2:-0}"
JOB_TIME="${3:-04:00:00}"
JOB_MEM="${4:-8G}"
SIF="${5:-${BULKANNEX_SIF:-}}"

RESULTS_DIR="$(realpath "${RESULTS_DIR}")"

# ---- Resolve Singularity image -----------------------------------------------
if [ -z "${SIF}" ]; then
    CANDIDATE="$(realpath "${SCRIPT_DIR}/../../containers/bulkannex_r/bulkannex_r_1.0.0.sif" 2>/dev/null || true)"
    if [ -f "${CANDIDATE}" ]; then
        SIF="${CANDIDATE}"
    else
        SIF="docker://damouzo/bulkannex_r:1.0.0"
        echo "  NOTE: local SIF not found — will pull from Docker registry on compute node."
        echo "  To use a pre-pulled SIF, set: export BULKANNEX_SIF=/path/to/bulkannex_r_1.0.0.sif"
        echo ""
    fi
fi

# ---- Auto-detect a free port -------------------------------------------------
if [ "${REQUESTED_PORT}" -eq 0 ]; then
    PORT=$(python3 -c "
import socket
s = socket.socket()
s.bind(('', 0))
print(s.getsockname()[1])
s.close()
" 2>/dev/null) || PORT=$(( RANDOM % 10000 + 20000 ))
else
    PORT="${REQUESTED_PORT}"
fi

FRONTEND="login.hpc.qmul.ac.uk"
USER_NAME="${USER}"

echo "========================================================"
echo "  bulkAnnex Dashboard — HPC Launcher (Apocrita/SLURM)"
echo "========================================================"
echo "  Results dir : ${RESULTS_DIR}"
echo "  App port    : ${PORT}"
echo "  SLURM time  : ${JOB_TIME}"
echo "  SLURM mem   : ${JOB_MEM}"
echo "  Container   : ${SIF}"
echo ""
echo "  Requesting SLURM allocation..."
echo "  SSH tunnel instructions will be printed once the job starts."
echo "========================================================"
echo ""

# ---- Launch via srun + singularity -------------------------------------------
srun \
    --job-name=bulkannex_dashboard \
    --cpus-per-task=2 \
    --mem="${JOB_MEM}" \
    --time="${JOB_TIME}" \
    --pty \
    bash -c "
        NODE=\$(hostname -s)
        FQDN=\$(hostname -f)

        echo ''
        echo '╔══════════════════════════════════════════════════════════════════════╗'
        echo '║  Dashboard running on: '\"\${NODE}\"'                                  ║'
        echo '╠══════════════════════════════════════════════════════════════════════╣'
        echo '║                                                                      ║'
        echo '║  Option A — VSCode Remote SSH (if you opened this via VSCode):       ║'
        echo '║  In a NEW VSCode terminal run:                                       ║'
        echo \"║    ssh -N -L ${PORT}:localhost:${PORT} \${FQDN}\"
        echo '║                                                                      ║'
        echo '║  Option B — Local terminal (MobaXterm / Mac / Linux):               ║'
        echo '║  On YOUR LAPTOP open a NEW terminal and run:                         ║'
        echo \"║    ssh -N -L ${PORT}:\${NODE}:${PORT} ${USER_NAME}@${FRONTEND}\"
        echo '║                                                                      ║'
        echo '║  Then open in browser: http://localhost:${PORT}                      ║'
        echo '║  Keep THIS terminal open. Ctrl-C to stop the dashboard.             ║'
        echo '╚══════════════════════════════════════════════════════════════════════╝'
        echo ''
        export BULKANNEX_RESULTS_DIR='${RESULTS_DIR}'
        singularity exec \
            --bind '${RESULTS_DIR}' \
            --bind '${SCRIPT_DIR}' \
            '${SIF}' \
            Rscript -e \"shiny::runApp('${SCRIPT_DIR}', port = ${PORT}, host = '0.0.0.0', launch.browser = FALSE)\"
    "

echo ""
echo "Dashboard session ended. To relaunch:"
echo "  BULKANNEX_SIF=${SIF} bash ${SCRIPT_DIR}/launch_dashboard_hpc.sh ${RESULTS_DIR}"
echo ""
