#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[1/12] Associating probes with restriction fragments..."

sshicstuff associate \
  -m "${COOL_INPUT}" \
  -c "${CAPTURE_OLIGOS}" \
  -o "${CAPTURE_ASSOCIATED}"

echo "Output:"
echo "  ${CAPTURE_ASSOCIATED}"