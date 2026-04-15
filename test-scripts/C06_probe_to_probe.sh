#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/00_set_variables.sh"

echo "[6/11] Building the probe-to-probe matrix, plotting a heatmap, and exporting a Cooler file..."

sshicstuff probe2probe \
  -c "${CAPTURE_ASSOCIATED}" \
  -f "${FILTERED_CONTACTS}" \
  -o "${PROBE2PROBE_MATRIX}" \
  --normalize \
  -P \
  --plot-format png \
  --colormap YlOrBr \
  --export-to-cooler \
  -F

echo "Outputs:"
echo "  ${PROBE2PROBE_MATRIX}"
echo "  ${PROBE2PROBE_MATRIX%.tsv}.png"
echo "  ${PROBE2PROBE_MATRIX%.tsv}.cool"
echo "  ${PROBE2PROBE_MATRIX%.tsv}_bins.tsv"
