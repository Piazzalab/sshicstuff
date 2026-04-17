#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[4/12] Computing fragment-level coverage..."

sshicstuff coverage \
  -m "${COOL_INPUT}" \
  -o "${OUTPUTS_DIR}" \
  -F

echo "[4/12] Computing binned coverage at 5 kb..."

sshicstuff coverage \
  -m "${COOL_INPUT}" \
  -c "${CHROM_COORDS}" \
  -b 5000 \
  -o "${OUTPUTS_DIR}" \
  -F

echo "Outputs:"
echo "  ${COVERAGE_FRAGMENT}"
echo "  ${COVERAGE_5KB}"