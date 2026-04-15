#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[4/11] Computing fragment-level coverage..."

sshicstuff coverage \
  -m "${GRAAL_MATRIX}" \
  -f "${FRAGMENTS_LIST}" \
  --outdir "${OUTPUTS_DIR}" \
  -F

echo "[4/11] Computing binned coverage at 5 kb..."

sshicstuff coverage \
  -m "${GRAAL_MATRIX}" \
  -f "${FRAGMENTS_LIST}" \
  -c "${CHROM_COORDS}" \
  -b 5000 \
  --outdir "${OUTPUTS_DIR}" \
  -F

echo "Outputs:"
echo "  ${OUTPUTS_DIR}/Graal_sample_for_pipeline_contacts_coverage.bedgraph"
echo "  ${OUTPUTS_DIR}/Graal_sample_for_pipeline_contacts_coverage_5kb.bedgraph"
