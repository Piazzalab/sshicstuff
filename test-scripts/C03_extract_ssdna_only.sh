#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[3/12] Extracting ssDNA-only cools..."

sshicstuff ssdnaonly \
  -m "${COOL_INPUT}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -F

echo "Outputs:"
echo "  ${SSDNA_COOL}"
echo "  ${SSDNA_TO_SSDNA_COOL}"