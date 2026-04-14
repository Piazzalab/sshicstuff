#!/usr/bin/env bash
set -euo pipefail

# Example configuration for the ssHiCstuff CLI workflow.
# Edit these paths to match your local test data layout.


# Project root from this script location:
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../.."
DATA_DIR="${PROJECT_DIR}/test_data"
INPUTS_DIR="${DATA_DIR}/inputs"
OUTPUTS_DIR="${DATA_DIR}/outputs_cli"

mkdir -p "${OUTPUTS_DIR}"

# Input files
GRAAL_MATRIX="${INPUTS_DIR}/Graal_sample_for_pipeline.txt"
CAPTURE_OLIGOS="${INPUTS_DIR}/capture_oligo_for_pipeline.csv"
FRAGMENTS_LIST="${INPUTS_DIR}/fragments_list_for_pipeline.txt"
CHROM_COORDS="${INPUTS_DIR}/chr_coordinates_for_pipeline.tsv"
GROUPS_TABLE="${INPUTS_DIR}/probe_groups_for_pipeline.tsv"

# Derived files
CAPTURE_ASSOCIATED="${OUTPUTS_DIR}/capture_oligo_for_pipeline_fragments_associated.csv"
FILTERED_CONTACTS="${OUTPUTS_DIR}/Graal_sample_for_pipeline_filtered.tsv"
PROFILE_0KB_CONTACTS="${OUTPUTS_DIR}/Graal_sample_for_pipeline_0kb_profile_contacts.tsv"
PROFILE_0KB_FREQ="${OUTPUTS_DIR}/Graal_sample_for_pipeline_0kb_profile_frequencies.tsv"
PROFILE_1KB_FREQ="${OUTPUTS_DIR}/Graal_sample_for_pipeline_1kb_profile_frequencies.tsv"
PROFILE_10KB_FREQ="${OUTPUTS_DIR}/Graal_sample_for_pipeline_10kb_profile_frequencies.tsv"
PROBE2PROBE_MATRIX="${OUTPUTS_DIR}/Graal_sample_for_pipeline_probe_to_probe_matrix.tsv"

echo "PROJECT_DIR=${PROJECT_DIR}"
echo "INPUTS_DIR=${INPUTS_DIR}"
echo "OUTPUTS_DIR=${OUTPUTS_DIR}"
