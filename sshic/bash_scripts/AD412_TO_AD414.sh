#!/usr/bin/bash

base_dir="/home/nicolas/Documents/Projects/ssHiC"
script="${base_dir}/hic_ssdna/sshic/pipeline.py"

fragments="${base_dir}/data/samples/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt"
oligos="${base_dir}/data/samples/inputs/capture_oligo_positions.csv"
additional="${base_dir}/data/samples/inputs/additional_probe_groups.tsv"
centromeres="${base_dir}/data/samples/inputs/S288c_chr_centro_coordinates.tsv"
binning="1000 2000 3000 5000 10000 20000 40000 50000 80000 100000"
ws_centros=150000
ws_telos=150000
excluded_chr="chr2 chr3 2_micron mitochondrion chr_artificial"

run_pipeline() {
    python3 "$script" -s "$sample" \
                      -f "$fragments" \
                      -o "$oligos" \
                      -c "$centromeres" \
                      -r "$reference" \
                      -b $binning \
                      -a "$additional" \
                      --window-size-centros $ws_centros \
                      --window-size-telos $ws_telos \
                      --excluded-chr $excluded_chr \
                      --exclude-probe-chr
}

samples_dir="${base_dir}/data/samples/AD412_AD414/pcrdupkept"
refs_dir="${base_dir}/data/samples/inputs/refs"

samples=(
    "${samples_dir}/AD412_S288c_DSB_LY_Capture_artificial_cutsite_q20_PCRdupkept.txt"
    "${samples_dir}/AD413_S288c_DSB_LY_Capture_artificial_cutsite_q20_PCRdupkept.txt"
    "${samples_dir}/AD414_S288c_DSB_LY_Capture_artificial_cutsite_q20_PCRdupkept.txt"

)

references=(
  "${refs_dir}/ref_WT4h_v3.tsv"
  "${refs_dir}/ref_WT4h_v3.tsv"
  "${refs_dir}/ref_WT4h_v3.tsv"
)



for ((i=0; i<${#samples[@]}; i++)); do
    sample="${samples[i]}"
    reference="${references[i]}"
    run_pipeline
done