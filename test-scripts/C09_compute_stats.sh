#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[9/12] Computing probe-level contact statistics..."

sshicstuff stats \
  -m "${COOL_INPUT}" \
  -p "${PROFILE_0KB_CONTACTS}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -C "${CHROM_COORDS}" \
  -o "${OUTPUTS_DIR}" \
  -r 50000 \
  -F

echo "Outputs:"
echo "  ${STATS_TABLE}"
echo "  ${CHR_FREQ_TABLE}"
echo "  ${INTER_CHR_FREQ_TABLE}"