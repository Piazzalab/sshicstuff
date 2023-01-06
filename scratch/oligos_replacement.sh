#!/bin/sh

genome="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/oligos_replacement/inputs/S288c_DSB_LY_Capture_original.fa"
oligos="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/oligos_replacement/inputs/oligo_positions.csv"
output="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/oligos_replacement/outputs/S288c_DSB_LY_Capture_original_artificial.fa"
bed="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/oligos_replacement/outputs/chr_artificial_coordinates.bed"

hic_ssdna.oligos_replacement -i $genome -o $output -c $oligos -b $bed -s 150