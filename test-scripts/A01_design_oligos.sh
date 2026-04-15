#!/usr/bin/env bash
set -euo pipefail

# Directory containing this example
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Inputs
GENOME_NAME="S288c_DSB_chr3_rDNA"
FASTA="${SCRIPT_DIR}/${GENOME_NAME}.fa"

# Outputs
OUTDIR="${SCRIPT_DIR}/outputs"
mkdir -p "${OUTDIR}"

# Output filenames
O4S_OUTPUT_RAW="output_o4s_raw.fa"
O4S_OUTPUT_SNP="output_o4s_snp.fa"
ANNEALING_CSV="annealing_oligo_positions.csv"
CAPTURE_CSV="capture_oligo_positions.csv"

echo "SCRIPT_DIR=${SCRIPT_DIR}"
echo "FASTA=${FASTA}"
echo "OUTDIR=${OUTDIR}"

if [[ ! -f "${FASTA}" ]]; then
    echo "ERROR: FASTA file not found: ${FASTA}" >&2
    exit 1
fi

echo "[Design] Running oligo design on the test genome..."

sshicstuff design \
  -f "${FASTA}" \
  --forward-intervals "chr5:118710-133000" \
  --reverse-intervals "chr5:100000-118710" \
  --site "GATC" \
  --secondary-sites "CAATTG,AATATT" \
  --size "80" \
  --site-start "70" \
  --capture-size "60" \
  --outdir "${OUTDIR}" \
  --o4s-output-raw "${O4S_OUTPUT_RAW}" \
  --o4s-output-snp "${O4S_OUTPUT_SNP}" \
  --annealing-csv "${ANNEALING_CSV}" \
  --capture-csv "${CAPTURE_CSV}"

echo
echo "Design completed successfully."
echo "Outputs written to: ${OUTDIR}"
