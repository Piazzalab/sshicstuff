import numpy as np
import sys
from Bio.SeqIO.FastaIO import FastaIterator


def is_debug() -> bool:
    gettrace = getattr(sys, 'gettrace', None)

    if gettrace is None:
        return False
    else:
        v = gettrace()
        if v is None:
            return False
        else:
            return True


def build_bins_from_genome(path_to_genome: str,
                           bin_size: int):

    genome = open(path_to_genome, "r")
    nb_bins_per_chr = {}
    for record in FastaIterator(genome):
        chr_id = record.id
        nb_bins_per_chr[chr_id] = len(str(record.seq)) // bin_size + 1
    genome.close()
    total_nb_bins = np.sum(list(nb_bins_per_chr.values()))
    chr_bins = np.zeros(total_nb_bins, dtype='<U32')
    genome_bins = np.zeros(total_nb_bins, dtype='<U32')

    big_counter = 0
    counter = 0
    for ii_chr in nb_bins_per_chr:
        nb_bins = nb_bins_per_chr[ii_chr]
        for ii_bin in range(0, nb_bins, 1):
            start = ii_bin * bin_size
            stop = (ii_bin + 1) * bin_size
            chr_bins[counter] = ii_chr + '_' + str(start) + '_' + str(stop)
            genome_bins[counter] = str(start+big_counter*bin_size) + '_' + str(stop+big_counter*bin_size)
            counter += 1
        big_counter += nb_bins

    return genome_bins, chr_bins


if __name__ == "__main__":
    artificial_genome = "../../../oligos_replacement/outputs/S288c_DSB_LY_Capture_original_artificial_nicolas.fa"
    filtered_contacts = "../../../contacts_filter/outputs/contacts_filtered_nicolas.csv"
    genome_bins_names, chr_bins_names = build_bins_from_genome(artificial_genome, 10000)

    print('--- DONE ---')
