#!/usr/bin/env bash
set -euo pipefail

# Directory containing this script:
# C-sshicstuff-pipeline/cli
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Example root:
# C-sshicstuff-pipeline
EXAMPLE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUTS_DIR="${EXAMPLE_DIR}/inputs"
OUTPUTS_DIR="${EXAMPLE_DIR}/outputs"

mkdir -p "${OUTPUTS_DIR}"

# Input files
GRAAL_MATRIX="${INPUTS_DIR}/AD433_abs_graal_fragments_weighted.txt"
CAPTURE_OLIGOS="${INPUTS_DIR}/capture_oligo_positions.csv"
FRAGMENTS_LIST="${INPUTS_DIR}/AD433_digested_fragments_list.txt"
CHROM_COORDS="${INPUTS_DIR}/chr_coordinates_for_pipeline.tsv"
GROUPS_TABLE="${INPUTS_DIR}/additional_probe_groups.tsv"

# Derived files
CAPTURE_ASSOCIATED="${OUTPUTS_DIR}/capture_oligo_positions_fragments_associated.csv"
FILTERED_CONTACTS="${OUTPUTS_DIR}/AD433_sub4M_filtered.tsv"

PROFILE_0KB_CONTACTS="${OUTPUTS_DIR}/AD433_sub4M_0kb_profile_contacts.tsv"
PROFILE_0KB_FREQ="${OUTPUTS_DIR}/AD433_sub4M_0kb_profile_frequencies.tsv"
PROFILE_1KB_FREQ="${OUTPUTS_DIR}/AD433_sub4M_1kb_profile_frequencies.tsv"
PROFILE_10KB_FREQ="${OUTPUTS_DIR}/AD433_sub4M_10kb_profile_frequencies.tsv"

PROBE2PROBE_MATRIX="${OUTPUTS_DIR}/AD433_sub4M_probe_to_probe_matrix.tsv"

echo "EXAMPLE_DIR=${EXAMPLE_DIR}"
echo "INPUTS_DIR=${INPUTS_DIR}"
echo "OUTPUTS_DIR=${OUTPUTS_DIR}"

for f in \
  "${GRAAL_MATRIX}" \
  "${CAPTURE_OLIGOS}" \
  "${FRAGMENTS_LIST}" \
  "${CHROM_COORDS}" \
  "${GROUPS_TABLE}" 
do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: missing input file: ${f}" >&2
    exit 1
  fi
done
