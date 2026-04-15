#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FASTQ_DIR="${SCRIPT_DIR}/B-hicstuff-mapping/fastq"
PIPELINE_INPUTS_DIR="${SCRIPT_DIR}/C-sshicstuff-pipeline/inputs"

mkdir -p "${FASTQ_DIR}" "${PIPELINE_INPUTS_DIR}"

echo "Downloading Hi-C FASTQ files..."
wget -O "${FASTQ_DIR}/AD433_sub4M.end1.fastq.gz" "URL_FASTQ_R1"
wget -O "${FASTQ_DIR}/AD433_sub4M.end2.fastq.gz" "URL_FASTQ_R2"

echo "Downloading precomputed ssHiCstuff pipeline inputs..."
wget -O "${PIPELINE_INPUTS_DIR}/AD433_abs_graal_fragments_weighted.txt" "URL_GRAAL"
wget -O "${PIPELINE_INPUTS_DIR}/AD433_digested_fragments_list.txt" "URL_FRAGMENTS"
wget -O "${PIPELINE_INPUTS_DIR}/AD433_info_contigs.txt" "URL_CONTIGS"

echo "Done."
