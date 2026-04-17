#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[10/12] Aggregating normalized 1 kb profiles around centromeres..."

sshicstuff aggregate \
  -p "${PROFILE_1KB_FREQ}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -C "${CHROM_COORDS}" \
  -o "${OUTPUTS_DIR}" \
  -N \
  --cen \
  -w 150000 \
  -F

echo "Aggregation written in:"
echo "  ${OUTPUTS_DIR}"