#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/00_set_variables.sh"

echo "[3/11] Extracting ssDNA-only sparse matrix..."

sshicstuff ssdnaonly \
  -m "${GRAAL_MATRIX}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -F

echo "Output:"
echo "  ${OUTPUTS_DIR}/AD433_sub4M_abs_graal_fragments_weighted_ssdna_to_ssdna_only.txt"
echo "  ${OUTPUTS_DIR}/AD433_sub4M_abs_graal_fragments_weighted_ssdna_only"
