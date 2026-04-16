#!/usr/bin/env bash
# launch_dashboard.sh — launch the bulkAnnex Shiny dashboard locally
#
# Usage: bash launch_dashboard.sh [results_dir] [port]
#
#   results_dir  Path to the Nextflow pipeline output directory.
#                Defaults to the parent of the directory containing this script
#                (i.e. <outdir>/ when the dashboard lives at <outdir>/dashboard/).
#   port         Port to serve the app on. Default: 3838.
#
# The script exports BULKANNEX_RESULTS_DIR so that global.R knows where to load data from.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-$(dirname "${SCRIPT_DIR}")}"
PORT="${2:-3838}"

RESULTS_DIR="$(realpath "${RESULTS_DIR}")"

if [ ! -f "${SCRIPT_DIR}/app.R" ]; then
    echo "ERROR: app.R not found in ${SCRIPT_DIR}"
    exit 1
fi

if [ ! -d "${RESULTS_DIR}" ]; then
    echo "ERROR: Results directory not found: ${RESULTS_DIR}"
    echo "Usage: bash launch_dashboard.sh [results_dir] [port]"
    exit 1
fi

export BULKANNEX_RESULTS_DIR="${RESULTS_DIR}"

echo "bulkAnnex dashboard"
echo "  App dir    : ${SCRIPT_DIR}"
echo "  Results dir: ${RESULTS_DIR}"
echo "  URL        : http://localhost:${PORT}"
echo ""

Rscript -e "shiny::runApp('${SCRIPT_DIR}', port = ${PORT}, launch.browser = TRUE)"
