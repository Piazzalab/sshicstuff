#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[12/12] Running the complete ssHiCstuff workflow on the test dataset..."

sshicstuff pipeline \
  -m "${COOL_INPUT}" \
  -c "${CAPTURE_OLIGOS}" \
  -C "${CHROM_COORDS}" \
  -a "${GROUPS_TABLE}" \
  -o "${OUTPUTS_DIR}/full_pipeline" \
  -b 1000 \
  -b 10000 \
  --window-cen 150000 \
  --window-telo 15000 \
  --bin-cen 10000 \
  --bin-telo 1000 \
  --balanced-stats \
  -r 50000 \
  -F \
  -N

echo
echo "Pipeline completed successfully."
echo "Main output directory:"
echo "  ${OUTPUTS_DIR}/full_pipeline"