#!/usr/bin/env sh

python3 ./pipeline.py  \
--samplesheet ../data/inputs/samplesheet.csv \
--fragments-list ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_v6_DpnIIHinfI.txt \
--outputs-dir ../data/outputs \
--chromosomes-arms-coordinates ../data/inputs/S288c_chr_centro_coordinates_S288c_DSB_LY_Capture_artificial_v6.tsv \
--oligos-capture ../data/inputs/capture_oligo_positions_v6.csv \
--additional-groups ../data/inputs/additional_probe_groups.tsv \
--binning-sizes 1000 2000 5000 10000 20000 40000 50000 80000 100000 \
--centromeres-aggregated-window-size 150000 \
--telomeres-aggregated-window-size 15000 \
--centromeres-aggregated-binning 10000 \
--telomeres-aggregated-binning 1000 \
--aggregate-by-arm-lengths \
--excluded-chr chr2 chr3 2_micron mitochondrion chr_artificial_donor chr_artificial_ssDNA \
--exclude-probe-chr