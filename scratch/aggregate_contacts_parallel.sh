#!/bin/bash


samples_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs"
outputs_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/ssHiC_filtered"
mkdir -p "${outputs_dir}/"

tasks=()
for samp in $samples_dir/10kb/*frequencies.tsv; do
    filename=$(basename "$samp")
    tasks+=("$filename")
done


parallel --jobs 1 '

	#### ----------- ####
	#### centromeres ####
	#### ----------- ####
	
    filename={}
    echo $filename

    contacts="'${samples_dir}'/10kb/${filename}"
    cen="'${samples_dir}'/S288c_chr_centro_coordinates.tsv"
    output="'${outputs_dir}'/"
    window=150000

    python3 /home/nicolas/Documents/Projects/ssHiC/hic_ssdna/src/hic_ssdna/aggregate/centromeres.py -b $contacts -c $cen -w $window -o $output 

    #### ----------- ####
	####  cohesins   ####
	#### ----------- ####

    scores_list=(50 100 200 500 1000 2000)
    for score in ${scores_list[@]}; do

		contacts="'${samples_dir}'/1kb/${filename}"
	    peaks="'${samples_dir}'/HB65_reference_peaks_score50min.bed"
	    output="'${outputs_dir}'/"
	    window=15000
	    
	    #python3 /home/nicolas/Documents/Projects/ssHiC/hic_ssdna/src/hic_ssdna/aggregate/cohesins.py -b $contacts -c $peaks  -w $window -s $score -o $output
    done

    #### ----------- ####
	#### telomeres	 ####
	#### ----------- ####

    contacts="'${samples_dir}'/10kb/${filename}"
	telo="'${samples_dir}'/S288c_chr_centro_coordinates.tsv"
	output="'${outputs_dir}'/"
	window=100000

	#python3 /home/nicolas/Documents/Projects/ssHiC/hic_ssdna/src/hic_ssdna/aggregate/telomeres.py -b $contacts -c $telo -w $window -o $output
	
' ::: "${tasks[@]}"




