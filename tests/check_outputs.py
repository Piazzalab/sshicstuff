# test pour savoir si les oligos ont bien été remplacés correctement
import pandas as pd
from Bio import SeqIO

from Projet.hic_ssdna.src.hic_ssdna.oligos_replacement import oligo_correction

genome_path = "/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/" \
         "oligos_replacement/inputs/S288c_DSB_LY_Capture_original.fa"
oligos_path = "/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/" \
         "oligos_replacement/inputs/oligo_positions.csv"
output_path = "/Users/loqmenanani/OneDrive/ENS/L3_ENS/stage_l3/Projet/" \
         "oligos_replacement/outputs/S288c_DSB_LY_Capture_original_artificial.fa"


def startend(num_oligo, dataframe):
    start = dataframe['start'][num_oligo] - 1
    end = dataframe['end'][num_oligo]
    return start, end


ancien = SeqIO.parse(genome_path, 'fasta')
nouveau = SeqIO.parse(output_path, 'fasta')
oligos = oligo_correction(oligos_path)
reads_size = 150

old = []
for record in ancien:
    old.append(record)

new = []
for records in nouveau:
    new.append(records)

print(len(new))
print(len(old))
print(len(new[4].seq) == len(old[4].seq))
for i in range(21):
    start, end = startend(i, oligos)
    # print(i)
    print(new[4].seq[start:end] == oligos['sequence_modified'][i], 'nouveau = oligo ?')
    print(new[4].seq[start:end] != old[4].seq[start:end], 'nouveau != ancien ?')
    print(new[4].seq[start - 3:start] == old[4].seq[start - 3:start], ' pas de décalage ')
    # print(new[4].seq[start:end], 'nouveau')
    # print(oligos['sequence_modified'][i], 'oligo')
    # print(old[4].seq[start:end], 'ancienne sequence')
    # print('\n')
for k in range(len(old)):
    print(new[k].seq == old[k].seq, ' : meme sequence', new[k].id, old[k].id)
    print(new[k].id == old[k].id, ' : meme id')

artificial = str(new[-1].seq)

# print(artificial)
print(len(artificial))

for k in range(len(oligos)):
    start, end = startend(k, oligos)
    print(oligos['sequence_original'][k] in artificial)
    print(str(old[4].seq[start-reads_size:end+reads_size]) in artificial)
