#!/usr/bin/bash

base_dir="/home/nicolas/Documents/Projects/ssHiC"
script="${base_dir}/hic_ssdna/sshic/pipeline.py"

fragments="${base_dir}/data/samples/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt"
oligos="${base_dir}/data/samples/inputs/capture_oligo_positions.csv"
centromeres="${base_dir}/data/samples/inputs/S288c_chr_centro_coordinates.tsv"
binning="1000 2000 3000 5000 10000 20000 40000 50000 80000 100000"
ws_centros=150000
ws_telos=150000
excluded_chr="chr2 chr3 chr5 2_micron mitochondrion chr_artificial"

run_pipeline() {
    python3 "$script" -s "$sample" \
                      -f "$fragments" \
                      -o "$oligos" \
                      -c "$centromeres" \
                      -r "$reference" \
                      -b $binning \
                      --window-size-centros $ws_centros \
                      --window-size-telos $ws_telos \
                      --excluded-chr "$excluded_chr" \
                      --inter-norm
}

mode_pcr=$1

if [[ $mode_pcr == "pcrfree" ]]; then
    samples=(
        "${base_dir}/data/samples/pcrfree/AD403_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
        "${base_dir}/data/samples/pcrfree/AD404_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
        "${base_dir}/data/samples/pcrfree/AD401_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
        "${base_dir}/data/samples/pcrfree/AD402_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
        "${base_dir}/data/samples/pcrfree/AD405_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
        "${base_dir}/data/samples/pcrfree/AD406_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
        "${base_dir}/data/samples/pcrfree/AD407_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt"
    )
    references=(
        "${base_dir}/data/samples/inputs/wt4h_pcrfree.tsv"
        "${base_dir}/data/samples/inputs/wt2h_pcrfree.tsv"
        "${base_dir}/data/samples/pcrfree/AD403/AD403_global_statistics.tsv"
        "${base_dir}/data/samples/pcrfree/AD404/AD404_global_statistics.tsv"
        "${base_dir}/data/samples/pcrfree/AD404/AD404_global_statistics.tsv"
        "${base_dir}/data/samples/pcrfree/AD403/AD403_global_statistics.tsv"
        "${base_dir}/data/samples/pcrfree/AD402/AD402_global_statistics.tsv"
    )

elif [[ $mode_pcr == "pcrdupkept" ]]; then
    samples=(
        "${base_dir}/data/samples/pcrdupkept/AD403_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
        "${base_dir}/data/samples/pcrdupkept/AD404_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
        "${base_dir}/data/samples/pcrdupkept/AD401_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
        "${base_dir}/data/samples/pcrdupkept/AD402_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
        "${base_dir}/data/samples/pcrdupkept/AD405_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
        "${base_dir}/data/samples/pcrdupkept/AD406_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
        "${base_dir}/data/samples/pcrdupkept/AD407_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
    )
    references=(
        "${base_dir}/data/samples/inputs/wt4h_pcrdupkept.tsv"
        "${base_dir}/data/samples/inputs/wt2h_pcrdupkept.tsv"
        "${base_dir}/data/samples/pcrdupkept/AD403/AD403_global_statistics.tsv"
        "${base_dir}/data/samples/pcrdupkept/AD404/AD404_global_statistics.tsv"
        "${base_dir}/data/samples/pcrdupkept/AD404/AD404_global_statistics.tsv"
        "${base_dir}/data/samples/pcrdupkept/AD403/AD403_global_statistics.tsv"
        "${base_dir}/data/samples/pcrdupkept/AD404/AD404_global_statistics.tsv"
    )
else
    echo "Mode PCR invalide. Utilisez 'pcrfree' ou 'pcrdupkept'."
    exit 1
fi

for ((i=0; i<${#samples[@]}; i++)); do
    sample="${samples[i]}"
    reference="${references[i]}"
    run_pipeline
done