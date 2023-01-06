#!/bin/bash


samples_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/contacts_binning/inputs/ssHiC_filtered"
outputs_dir="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/contacts_binning/outputs/ssHiC_filtered"
mkdir -p "${outputs_dir}/"


sizes=(0 1000 2000 5000 10000 20000 40000 80000 100000)

for bin_size in ${sizes[@]}; do
	echo "bin of size : ${bin_size} bp"
	for samp in $samples_dir/*.csv; do
		filename=$(basename "$samp")
		echo $filename
		
		genome="/home/nicolas/Documents/Projects/ssHiC/bash_scripts/contacts_binning/inputs/S288c_DSB_LY_capture_artificial.fa"
		contacts="${samples_dir}/${filename}"
		output="${outputs_dir}/"
		
		python3 /home/nicolas/Documents/Projects/ssHiC/hic_ssdna/src/hic_ssdna/contacts/binning.py -g $genome -c $contacts -b $bin_size -o $output 
	done
done

