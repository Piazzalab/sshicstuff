#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/00_set_variables.sh"

echo "[3/11] Extracting ssDNA-only sparse matrix..."

sshicstuff ssdnaonly \
  -m "${GRAAL_MATRIX}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}/Graal_sample_for_pipeline_ssdna_only.tsv" \
  -F

echo "Output:"
echo "  ${OUTPUTS_DIR}/Graal_sample_for_pipeline_ssdna_only.tsv"
