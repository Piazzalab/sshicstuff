#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[5/12] Filtering probe-associated contacts..."

sshicstuff filter \
  -m "${COOL_INPUT}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${FILTERED_COOL}" \
  -t "${FILTERED_TSV}" \
  -F

echo "Outputs:"
echo "  ${FILTERED_COOL}"
echo "  ${FILTERED_TSV}"