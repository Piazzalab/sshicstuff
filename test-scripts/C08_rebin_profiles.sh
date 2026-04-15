#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[8/11] Rebinning the normalized profile at 1 kb..."

sshicstuff rebin \
  -p "${PROFILE_0KB_FREQ}" \
  -C "${CHROM_COORDS}" \
  -b 1000 \
  -F

echo "[8/11] Plotting the 1 kb profile..."

sshicstuff plot4c \
  -p "${PROFILE_1KB_FREQ}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -C "${CHROM_COORDS}" \
  -e png \
  -H 600 \
  -W 1200 \
  -r 4 \
  -y 0 \
  -Y 0.01

echo "[8/11] Rebinning the normalized profile at 10 kb..."

sshicstuff rebin \
  -p "${PROFILE_0KB_FREQ}" \
  -C "${CHROM_COORDS}" \
  -b 10000 \
  -F

echo "[8/11] Plotting the 10 kb profile..."

sshicstuff plot4c \
  -p "${PROFILE_10KB_FREQ}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -o "${OUTPUTS_DIR}" \
  -C "${CHROM_COORDS}" \
  -e png \
  -H 600 \
  -W 1200 \
  -r 4 \
  -y 0 \
  -Y 0.03

echo "Outputs:"
echo "  ${PROFILE_1KB_FREQ}"
echo "  ${PROFILE_10KB_FREQ}"
echo "  ${OUTPUTS_DIR}/plots/"
