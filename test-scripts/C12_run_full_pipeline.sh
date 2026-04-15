#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/00_set_variables.sh"

echo "[Full pipeline] Running the complete ssHiCstuff workflow on the test dataset..."

sshicstuff pipeline \
  -m "${GRAAL_MATRIX}" \
  -c "${CAPTURE_OLIGOS}" \
  -f "${FRAGMENTS_LIST}" \
  -C "${CHROM_COORDS}" \
  -a "${GROUPS_TABLE}" \
  -o "${OUTPUTS_DIR}/pipeline_run" \
  -b 1000 \
  -b 10000 \
  --window-size-cen 150000 \
  --window-size-telo 15000 \
  --binning-aggregate-cen 10000 \
  --binning-aggregate-telo 1000 \
  -r 50000 \
  -F \
  -N

echo
echo "Pipeline completed successfully."
echo "Main output directory:"
echo "  ${OUTPUTS_DIR}/pipeline_run"
