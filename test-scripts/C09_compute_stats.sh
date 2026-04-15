#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[9/11] Computing probe-level contact statistics..."

sshicstuff stats \
  -m "${GRAAL_MATRIX}" \
  -p "${PROFILE_0KB_CONTACTS}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -C "${CHROM_COORDS}" \
  -o "${OUTPUTS_DIR}" \
  -r 50000 \
  -F

echo "Output:"
echo "  ${OUTPUTS_DIR}/Graal_sample_for_pipeline_statistics.tsv"
