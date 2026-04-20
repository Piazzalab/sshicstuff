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
  --balanced \
  -F

# Remove --balanced if ICE weights have not been computed yet (cooler balance).
# When --balanced is set, the denominator for capture-efficiency and
# ssDNA/dsDNA coverage ratios is derived from ICE-balanced pixel sums
# instead of raw counts, correcting for fragment-density and GC biases.

echo "Outputs:"
echo "  ${STATS_TABLE}"
echo "  ${CHR_FREQ_TABLE}"
echo "  ${INTER_CHR_FREQ_TABLE}"