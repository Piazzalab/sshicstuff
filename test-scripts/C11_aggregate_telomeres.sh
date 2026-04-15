#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[11/11] Aggregating normalized 1 kb profiles around telomeres..."

sshicstuff aggregate \
  -p "${PROFILE_1KB_FREQ}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -C "${CHROM_COORDS}" \
  -o "${OUTPUTS_DIR}" \
  -N \
  --tel \
  -w 15000

echo "Output directory:"
echo "  ${OUTPUTS_DIR}/aggregated/telomeres/"
