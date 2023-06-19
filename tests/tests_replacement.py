from oligos_replacement import replacement

genome_path = '/scratch/lanani/Stage L3/' \
'Projet/Pycharm/inputs/S288c_DSB_LY_Capture_original.fa'

oligos_path = '/scratch/lanani/Stage L3/' \
'Projet/Pycharm/inputs/oligo_positions_modified2.csv'


output_path = '/scratch/lanani/Stage L3/' \
'Projet/Pycharm/outputs/new_genome.fa'

bed_paths = '/scratch/lanani/Stage L3/' \
'Projet/Pycharm/outputs/chr_articial_coordinates.bed'

reads_size = 150

replacement(input_genome=genome_path, input_oligos=oligos_path , output_genome=output_path, bed_path=bed_paths, flanking_size=reads_size)

