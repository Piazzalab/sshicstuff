#!/bin/bash


samples_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/10kb_(PCRduplicates)"
scores_list=(50 100 200 500 1000 2000)

for samp in $samples_dir/*_frequencies_matrix.tsv; do
	filename=$(basename "$samp")
	echo $filename

	############################################################################################################
    ############################################################################################################
    echo "AGGREATE ON CENTROMERES"
	contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/10kb_(PCRduplicates)/${filename}"
	cen="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"
	output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/${filename}"
	window=150000
	

	python3 /home/nicolas/Documents/Projects/ssHiC/ssdna-hic/src/hic_ssdna/aggregate/centromeres.py -b $contacts -c $cen -w $window -o $output 
	
    ############################################################################################################
    ############################################################################################################
    echo "AGGREATE ON COHESINS PEAKS"
    for score in ${scores_list[@]}; do
        echo $score
	    contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/1kb_(PCRduplicates)/${filename}"
	    peaks="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/HB65_reference_peaks_score50min.bed"
	    output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/${filename}"
	    window=15000
	    
	    python3 /home/nicolas/Documents/Projects/ssHiC/ssdna-hic/src/hic_ssdna/aggregate/cohesins.py -b $contacts -c $peaks  -w $window -s $score -o $output
    done

    ############################################################################################################
    ############################################################################################################

    echo "AGGREATE ON TELOMERES"
    contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/10kb_(PCRduplicates)/${filename}"
	telo="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"
	output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregate_contacts/outputs/${filename}"
	window=100000

	python3 /home/nicolas/Documents/Projects/ssHiC/ssdna-hic/src/hic_ssdna/aggregate/telomeres.py -b $contacts -c $telo -w $window -o $output
	
done


