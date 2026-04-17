#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_ROOT="$(cd "${SCRIPT_DIR}/../test-data" && pwd)"

INPUTS_DIR="${TEST_ROOT}/inputs"


SAMPLE="AD433_sub4M"
OUTPUTS_DIR="${TEST_ROOT}/C-output-sshicstuff/${SAMPLE}"
mkdir -p "${OUTPUTS_DIR}"


# ---------------------------------------------------------------------
# Required inputs
# ---------------------------------------------------------------------
COOL_INPUT="${INPUTS_DIR}/${SAMPLE}.cool"
CAPTURE_OLIGOS="${INPUTS_DIR}/capture_oligo_positions.csv"
CHROM_COORDS="${INPUTS_DIR}/chr_coordinates_for_pipeline.tsv"
GROUPS_TABLE="${INPUTS_DIR}/additional_probe_groups.tsv"

COOL_BASENAME="$(basename "${COOL_INPUT}" .cool)"

# ---------------------------------------------------------------------
# Flat outputs
# ---------------------------------------------------------------------
CAPTURE_ASSOCIATED="${OUTPUTS_DIR}/capture_oligo_positions_fragments_associated.csv"

BALANCED_COOL="${COOL_INPUT}"  # same file if balance runs in-place

DSDNA_COOL="${OUTPUTS_DIR}/${COOL_BASENAME}_dsdna_only.cool"
SSDNA_COOL="${OUTPUTS_DIR}/${COOL_BASENAME}_ssdna_only.cool"
SSDNA_TO_SSDNA_COOL="${OUTPUTS_DIR}/${COOL_BASENAME}_ssdna_to_ssdna_only.cool"

FILTERED_COOL="${OUTPUTS_DIR}/${COOL_BASENAME}_filtered.cool"
FILTERED_TSV="${OUTPUTS_DIR}/${COOL_BASENAME}_filtered.tsv"

COVERAGE_FRAGMENT="${OUTPUTS_DIR}/${COOL_BASENAME}_coverage.counts.fragment_level.bedgraph"
COVERAGE_5KB="${OUTPUTS_DIR}/${COOL_BASENAME}_coverage.counts.5kb.bedgraph"

PROFILE_0KB_CONTACTS="${OUTPUTS_DIR}/${COOL_BASENAME}_0kb_profile_contacts.tsv"
PROFILE_0KB_FREQ="${OUTPUTS_DIR}/${COOL_BASENAME}_0kb_profile_frequencies.tsv"
PROFILE_1KB_FREQ="${OUTPUTS_DIR}/${COOL_BASENAME}_1kb_profile_frequencies.tsv"
PROFILE_10KB_FREQ="${OUTPUTS_DIR}/${COOL_BASENAME}_10kb_profile_frequencies.tsv"

PROBE2PROBE_MATRIX="${OUTPUTS_DIR}/${COOL_BASENAME}_probe_matrix.tsv"
PROBE2PROBE_PLOT="${OUTPUTS_DIR}/${COOL_BASENAME}_probe_matrix.png"
PROBE2PROBE_COOL="${OUTPUTS_DIR}/${COOL_BASENAME}_probe_matrix.cool"

STATS_TABLE="${OUTPUTS_DIR}/${COOL_BASENAME}_statistics.tsv"
CHR_FREQ_TABLE="${OUTPUTS_DIR}/${COOL_BASENAME}_norm_chr_freq.tsv"
INTER_CHR_FREQ_TABLE="${OUTPUTS_DIR}/${COOL_BASENAME}_norm_inter_chr_freq.tsv"

echo "INPUTS_DIR=${INPUTS_DIR}"
echo "OUTPUTS_DIR=${OUTPUTS_DIR}"
echo "COOL_INPUT=${COOL_INPUT}"

for f in \
  "${COOL_INPUT}" \
  "${CAPTURE_OLIGOS}" \
  "${CHROM_COORDS}" \
  "${GROUPS_TABLE}"
do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: missing input file: ${f}" >&2
    exit 1
  fi
done