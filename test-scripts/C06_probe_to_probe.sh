#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[6/12] Building the probe-to-probe matrix..."

sshicstuff probe2probe \
  -c "${CAPTURE_ASSOCIATED}" \
  -f "${FILTERED_TSV}" \
  -o "${PROBE2PROBE_MATRIX}" \
  --normalize \
  -P \
  --plot-format png \
  --colormap YlOrBr \
  -F

echo "Outputs:"
echo "  ${PROBE2PROBE_MATRIX}"
echo "  ${PROBE2PROBE_PLOT}"
echo "  ${PROBE2PROBE_COOL}"