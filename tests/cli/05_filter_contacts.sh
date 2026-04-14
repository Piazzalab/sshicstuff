#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/00_set_variables.sh"

echo "[5/11] Filtering contacts to retain pairs involving at least one probe-associated fragment..."

sshicstuff filter \
  -m "${GRAAL_MATRIX}" \
  -c "${CAPTURE_ASSOCIATED}" \
  -f "${FRAGMENTS_LIST}" \
  -o "${FILTERED_CONTACTS}" \
  -F

echo "Output:"
echo "  ${FILTERED_CONTACTS}"
