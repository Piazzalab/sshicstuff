#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[7/12] Building the raw 4C-like profile in absolute contacts and normalized ..."

sshicstuff profile \
  -c "${CAPTURE_ASSOCIATED}" \
  -C "${CHROM_COORDS}" \
  -f "${FILTERED_TSV}" \
  -a "${GROUPS_TABLE}" \
  -o "${PROFILE_0KB_CONTACTS}" \
  -F \
  -N


echo "[7/12] Plotting the normalized 0 kb 4C-like profiles..."

sshicstuff plot4c \
  -p "${PROFILE_0KB_FREQ}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -C "${CHROM_COORDS}" \
  -e png \
  -H 600 \
  -W 1200 \
  -r 4 \
  -y 0 \
  -Y 0.01

echo "Outputs:"
echo "  ${PROFILE_0KB_CONTACTS}"
echo "  ${PROFILE_0KB_FREQ}"