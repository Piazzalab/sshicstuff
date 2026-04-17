#!/usr/bin/env bash
set -euo pipefail

# Directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_ROOT="$(cd "${SCRIPT_DIR}/../test-data" && pwd)"

INPUTS_DIR="${TEST_ROOT}/inputs"
OUTPUTS_DIR="${TEST_ROOT}/C-output-sshicstuff"


mkdir -p "${OUTPUTS_DIR}"

# Shared sample name
SAMPLE="AD433_sub4M"

# Input files
GRAAL_MATRIX="${INPUTS_DIR}/${SAMPLE}_abs_graal_fragments_weighted_balanced.txt"
CAPTURE_OLIGOS="${INPUTS_DIR}/capture_oligo_positions.csv"
FRAGMENTS_LIST="${INPUTS_DIR}/${SAMPLE}_digested_fragments_list.txt"
INFO_CONTIGS="${INPUTS_DIR}/${SAMPLE}_info_contigs.txt"
CHROM_COORDS="${INPUTS_DIR}/chr_coordinates_for_pipeline.tsv"
GROUPS_TABLE="${INPUTS_DIR}/additional_probe_groups.tsv"

# Derived files
CAPTURE_ASSOCIATED="${OUTPUTS_DIR}/capture_oligo_positions_fragments_associated.csv"
FILTERED_CONTACTS="${OUTPUTS_DIR}/${SAMPLE}_filtered.tsv"

PROFILE_0KB_CONTACTS="${OUTPUTS_DIR}/${SAMPLE}_0kb_profile_contacts.tsv"
PROFILE_0KB_FREQ="${OUTPUTS_DIR}/${SAMPLE}_0kb_profile_frequencies.tsv"
PROFILE_1KB_FREQ="${OUTPUTS_DIR}/${SAMPLE}_1kb_profile_frequencies.tsv"
PROFILE_10KB_FREQ="${OUTPUTS_DIR}/${SAMPLE}_10kb_profile_frequencies.tsv"

PROBE2PROBE_MATRIX="${OUTPUTS_DIR}/${SAMPLE}_probe_to_probe_matrix.tsv"


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
