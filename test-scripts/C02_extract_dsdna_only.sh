#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[2/11] Extracting dsDNA-only sparse matrix..."

sshicstuff dsdnaonly \
  -m "${GRAAL_MATRIX}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -F

echo "Output:"
echo "  ${OUTPUTS_DIR}/Graal_sample_for_pipeline_dsdna_only.tsv"
