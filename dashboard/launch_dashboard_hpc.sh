#!/usr/bin/env bash
# launch_dashboard_hpc.sh — launch bulkAnnex dashboard on Apocrita (QMUL SLURM HPC)
#
# Requests an interactive SLURM allocation, starts the Shiny app, and prints
# SSH tunnel instructions so you can reach it from your local browser.
#
# Usage: bash launch_dashboard_hpc.sh [results_dir] [port] [time] [mem]
#
#   results_dir  Path to Nextflow outdir (default: parent of this script).
#   port         Port for the app; 0 = auto-detect a free port (default: 0).
#   time         SLURM time limit (default: 04:00:00).
#   mem          SLURM memory request (default: 8G).
#
# Example:
#   bash launch_dashboard_hpc.sh /data/projects/myproject/results 0 02:00:00 4G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-$(dirname "${SCRIPT_DIR}")}"
REQUESTED_PORT="${2:-0}"
JOB_TIME="${3:-04:00:00}"
JOB_MEM="${4:-8G}"

RESULTS_DIR="$(realpath "${RESULTS_DIR}")"

# ---- Auto-detect a free port if none requested ------------------------------
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

# ---- Print connection info --------------------------------------------------
echo "========================================================"
echo "  bulkAnnex Dashboard — HPC Launcher (Apocrita/SLURM)"
echo "========================================================"
echo ""
echo "  Results dir : ${RESULTS_DIR}"
echo "  App port    : ${PORT}"
echo "  SLURM time  : ${JOB_TIME}"
echo "  SLURM mem   : ${JOB_MEM}"
echo ""
echo "  Once the job starts and you see 'Listening on http://...',"
echo "  open a NEW terminal on your LOCAL machine and run:"
echo ""
echo "    ssh -L ${PORT}:\$(hostname -s):\${PORT} ${USER_NAME}@${FRONTEND}"
echo ""
echo "  (Replace \$(hostname -s) with the actual compute node name printed above)"
echo ""
echo "  Then open in your browser: http://localhost:${PORT}"
echo "========================================================"
echo ""

export BULKANNEX_RESULTS_DIR="${RESULTS_DIR}"

# ---- Launch via SLURM interactive session -----------------------------------
srun \
    --job-name=bulkannex_dashboard \
    --cpus-per-task=2 \
    --mem="${JOB_MEM}" \
    --time="${JOB_TIME}" \
    --pty \
    bash -c "
        echo \"Running on node: \$(hostname -s)\"
        echo \"Port           : ${PORT}\"
        echo \"Results dir    : ${RESULTS_DIR}\"
        echo \"\"
        export BULKANNEX_RESULTS_DIR='${RESULTS_DIR}'
        Rscript -e \"shiny::runApp('${SCRIPT_DIR}', port = ${PORT}, host = '0.0.0.0', launch.browser = FALSE)\"
    "
