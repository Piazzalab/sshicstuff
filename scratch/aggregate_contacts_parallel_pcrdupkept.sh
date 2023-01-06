#!/bin/bash


samples_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/10kb_(PCRduplicates)"


tasks=()
for samp in $samples_dir/*_frequencies_matrix.tsv; do
    filename=$(basename "$samp")
    tasks+=("$filename")
done


parallel --jobs 6 '
    filename={}
    echo $filename

    ############################################################################################################
    ############################################################################################################

    contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/10kb_(PCRduplicates)/${filename}"
    cen="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"
    output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/${filename}"
    window=150000

    python3 /home/nicolas/Documents/Projects/ssHiC/ssdna-hic/src/hic_ssdna/aggregate/centromeres.py -b $contacts -c $cen -w $window -o $output 

    ############################################################################################################
    ############################################################################################################

    scores_list=(50 100 200 500 1000 2000)
    for score in ${scores_list[@]}; do
	    contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/1kb_(PCRduplicates)/${filename}"
	    peaks="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/HB65_reference_peaks_score50min.bed"
	    output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/${filename}"
	    window=15000
	    
	    python3 /home/nicolas/Documents/Projects/ssHiC/ssdna-hic/src/hic_ssdna/aggregate/cohesins.py -b $contacts -c $peaks  -w $window -s $score -o $output
    done

    ############################################################################################################
    ############################################################################################################

    contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/10kb_(PCRduplicates)/${filename}"
	telo="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"
	output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/${filename}"
	window=100000

	python3 /home/nicolas/Documents/Projects/ssHiC/ssdna-hic/src/hic_ssdna/aggregate/telomeres.py -b $contacts -c $telo -w $window -o $output
	
' ::: "${tasks[@]}"




