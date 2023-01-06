#!/bin/bash


samples_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/compute_ratio/inputs/0kb"


for samp in $samples_dir/*_contacts_matrix.tsv; do
	filename=$(basename "$samp")
	echo $filename
	
	contacts="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/compute_ratio/inputs/0kb/${filename}"
	output="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/compute_ratio/outputs/${filename}"
	cis_range=50000
	python3 /home/nicolas/Documents/Projects/ssHiC/hic_ssdna/src/hic_ssdna/compute_ratio.py -c $contacts -r $cis_range -O $output 
done


