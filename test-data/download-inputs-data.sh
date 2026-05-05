#!/bin/bash
set -euo pipefail
ZENODO_ID="19596177"
DEST_DIR="inputs"
API_URL="https://zenodo.org/api/records/${ZENODO_ID}/files"
mkdir -p "${DEST_DIR}"
curl -sSL "${API_URL}" \
  | python3 -c "import sys, json; [print(e['key'] + '\t' + e['links']['content']) for e in json.load(sys.stdin)['entries']]" \
  | while IFS=$'\t' read -r name url; do
        echo ">>> ${name}"
        wget -q --show-progress -c -O "${DEST_DIR}/${name}" "${url}"
    done
