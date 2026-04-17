#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[2/12] Extracting dsDNA-only cool..."

sshicstuff dsdnaonly \
  -m "${COOL_INPUT}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -F

echo "Output:"
echo "  ${DSDNA_COOL}"