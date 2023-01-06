#!/bin/sh

# pip3 install oligos-replacement

fragments="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/contacts_filter/inputs/fragments_list.txt"
oligos="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/contacts_filter/inputs/capture_oligo_positions.csv"
contacts="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/contacts_filter/inputs/contacts.txt"
output="/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/contacts_filter/outputs/contacts_filtered.csv"

#oligos='/scratch/lanani/Stage_L3_cbp/Projet/Pycharm/outputs/new_capture_oligos.csv'
#fragments='/scratch/lanani/Stage_L3_cbp/Projet/Pycharm/inputs/fragments_list.txt'
#contacts='/scratch/lanani/Stage_L3_cbp/Projet/Pycharm/inputs/abs_fragments_contacts_weighted.txt'
#output='/scratch/lanani/Stage_L3_cbp/Projet/Pycharm/outputs/contacts_filtered.csv'

hic_ssdna.contacts_filter -o $oligos -f $fragments -c $contacts -O $output