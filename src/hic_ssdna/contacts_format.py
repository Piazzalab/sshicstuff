#! /usr/bin/env python3

import numpy as np
import pandas as pd
import math
import sys
import getopt
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
    nb_bins_per_chr: dict = {}
    for record in FastaIterator(genome):
        chr_id = record.id
        nb_bins_per_chr[chr_id] = len(str(record.seq)) // bin_size + 1
    genome.close()
    total_nb_bins = np.sum(list(nb_bins_per_chr.values()))
    bed_array: dict = {'chr': np.zeros(total_nb_bins, dtype='<U16'),
                       'chr_and_bins': np.zeros(total_nb_bins, dtype='<U64'),
                       'bins_in_chr': np.zeros(total_nb_bins, dtype='<U64'),
                       'bins_in_genome': np.zeros(total_nb_bins, dtype='<U64'),
                       }

    big_counter = 0
    counter = 0
    for ii_chr, nb_bins in nb_bins_per_chr.items():
        for ii_bin in range(0, nb_bins, 1):
            start = ii_bin * bin_size
            # stop = (ii_bin + 1) * bin_size

            bed_array['chr'][counter] = ii_chr
            bed_array['chr_and_bins'][counter] = ii_chr + '_' + str(start)
            bed_array['bins_in_chr'][counter] = start
            bed_array['bins_in_genome'][counter] = start + big_counter*bin_size

            counter += 1
        big_counter += nb_bins

    return bed_array


def frag2(x):
    if x == 'a':
        y = 'b'
    else:
        y = 'a'
    return y


def count_occurrences(df: pd.DataFrame,
                      x: str,
                      res1: dict,
                      res2: dict,
                      bin_size):

    y = frag2(x)
    for ii_f, f in enumerate(df['frag_'+x].values):
        if not pd.isna(df['name_'+x][ii_f]):
            if f not in res1:
                res1[f] = {}

            if f not in res2:
                res2[f] = {'type': df['type_'+x][ii_f],
                           'name': df['name_'+x][ii_f],
                           'chr': df['chr_'+x][ii_f],
                           'start': df['start_'+x][ii_f],
                           'end': df['end_'+x][ii_f]}

            chr_id = df['chr_'+y][ii_f]
            start = math.floor(df['start_'+y][ii_f]/bin_size)*bin_size
            # end = start+bin_size
            # bin_id = chr_id + '_' + str(start) + '_' + str(end)
            bin_id = chr_id + '_' + str(start)

            if bin_id not in res1[f]:
                res1[f][bin_id] = df['contacts'][ii_f]
            else:
                res1[f][bin_id] += df['contacts'][ii_f]
    return res1, res2


def get_fragments_dict(contacts_path: str,
                       bin_size: int):

    df_contacts_filtered = pd.read_csv(contacts_path, sep=',')
    fragments_contacts: dict = {}
    fragments_infos: dict = {}
    fragments_contacts, fragments_infos = count_occurrences(df=df_contacts_filtered,
                                                            x='a',
                                                            res1=fragments_contacts,
                                                            res2=fragments_infos, bin_size=bin_size)

    fragments_contacts, fragments_infos = count_occurrences(df=df_contacts_filtered,
                                                            x='b',
                                                            res1=fragments_contacts,
                                                            res2=fragments_infos, bin_size=bin_size)

    return fragments_contacts, fragments_infos


def set_fragments_contacts_bins(bed_bins: dict,
                                bins_contacts_dict: dict,
                                fragment_infos_dict: dict,
                                output_path: str):

    chromosomes = bed_bins['chr']
    chr_and_bins = bed_bins['chr_and_bins']
    bins_in_chr = bed_bins['bins_in_chr']
    bins_in_genome = bed_bins['bins_in_genome']
    df = pd.DataFrame({'chr': chromosomes, 'chr_bins':  bins_in_chr, 'genome_bins': bins_in_genome})

    nb_bins = len(bins_in_genome)
    fragments = []
    types = []
    names = []
    for f in bins_contacts_dict:
        contacts = np.zeros(nb_bins, dtype=int)
        fragments.append(f)
        t = fragment_infos_dict[f]['type']
        n = fragment_infos_dict[f]['name']
        types.append(t)
        names.append(n)
        for _bin in bins_contacts_dict[f]:
            idx = np.where(chr_and_bins == _bin)[0]
            contacts[idx] = bins_contacts_dict[f][_bin]
        df[f] = contacts / np.sum(contacts)
        df = df.rename(columns={f: str(f) + '--' + str(n) + '--' + str(t)})

    df.to_csv(output_path+'_frequencies_matrix.csv')


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    artificial_genome_path, filtered_contacts_path, output_path, bin_size = ['' for _ in range(4)]

    try:
        opts, args = getopt.getopt(argv, "hg:c:b:o:", ["--help",
                                                       "--genome",
                                                       "--contacts",
                                                       "--bin_size"
                                                       "--output"])
    except getopt.GetoptError:
        print('contacts filter arguments :\n'
              '-g <fasta_genome_input> (artificially generated with oligos_replacement.py) \n'
              '-c <filtered_contacts_input.csv> (contacts filtered with contacts_filter.py) \n'
              '-b <bin_size> (size of a bin, in bp) \n'
              '-o <output_file_name.csv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('contacts filter arguments :\n'
                  '-g <fasta_genome_input> (artificially generated with oligos_replacement.py) \n'
                  '-c <filtered_contacts_input.csv> (contacts filtered with contacts_filter.py) \n'
                  '-b <bin_size> (size of a bin, in bp) \n'
                  '-o <output_file_name.csv>')
            sys.exit()
        elif opt in ("-g", "--genome"):
            artificial_genome_path = arg
        elif opt in ("-c", "--contacts"):
            filtered_contacts_path = arg
        elif opt in ("-b", "--bin_size"):
            bin_size = arg
        elif opt in ("-o", "--output"):
            output_path = arg.split('.csv')[0]

    bin_size = int(bin_size)
    bed_pos = build_bins_from_genome(artificial_genome_path, bin_size=bin_size)
    contacts_dict, infos_dict = get_fragments_dict(contacts_path=filtered_contacts_path, bin_size=bin_size)
    set_fragments_contacts_bins(bed_bins=bed_pos,
                                bins_contacts_dict=contacts_dict,
                                fragment_infos_dict=infos_dict,
                                output_path=output_path)


def debug(artificial_genome_path: str,
          filtered_contacts_path: str,
          bin_size: int,
          output_path: str):

    bed_pos = build_bins_from_genome(artificial_genome_path, bin_size=bin_size)
    contacts_dict, infos_dict = get_fragments_dict(contacts_path=filtered_contacts_path, bin_size=bin_size)
    set_fragments_contacts_bins(bed_bins=bed_pos,
                                bins_contacts_dict=contacts_dict,
                                fragment_infos_dict=infos_dict,
                                output_path=output_path)


if __name__ == "__main__":
    if is_debug():
        artificial_genome = "../../../contacts_format/inputs/S288c_DSB_LY_capture_artificial.fa"
        filtered_contacts = "../../../contacts_format/inputs/contacts_filtered_nicolas.csv"
        output = "../../../contacts_format/outputs/frequencies_per_bin_matrix.csv"
        debug(artificial_genome_path=artificial_genome,
              filtered_contacts_path=filtered_contacts,
              bin_size=10000,
              output_path=output)
    else:
        main(sys.argv[1:])

    print('--- DONE ---')
