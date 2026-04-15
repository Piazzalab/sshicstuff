#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/C00_set_variables.sh"

echo "[1/11] Associating probes with restriction fragments..."

sshicstuff associate \
  -c "${CAPTURE_OLIGOS}" \
  -f "${FRAGMENTS_LIST}" \
  -o "${CAPTURE_ASSOCIATED}"

echo "Output:"
echo "  ${CAPTURE_ASSOCIATED}"
