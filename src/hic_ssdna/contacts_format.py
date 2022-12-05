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
    chr_bins = np.zeros(total_nb_bins, dtype='<U32')
    genome_bins = np.zeros(total_nb_bins, dtype='<U32')

    big_counter = 0
    counter = 0
    for ii_chr, nb_bins in nb_bins_per_chr.items():
        for ii_bin in range(0, nb_bins, 1):
            start = ii_bin * bin_size
            # stop = (ii_bin + 1) * bin_size
            # chr_bins[counter] = ii_chr + '_' + str(start) + '_' + str(stop)
            # genome_bins[counter] = str(start+big_counter*bin_size) + '_' + str(stop+big_counter*bin_size)
            chr_bins[counter] = ii_chr + '_' + str(start)
            genome_bins[counter] = str(start+big_counter*bin_size)
            counter += 1
        big_counter += nb_bins

    return genome_bins, chr_bins


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


def set_fragments_contacts_bins(bins_contacts_dict: dict,
                                fragment_infos_dict: dict,
                                chr_bins: np.ndarray,
                                genome_bins: np.ndarray,
                                output_path: str):

    merged_chr_and_genome_bins = np.zeros(len(chr_bins), dtype='<U64')
    for ii_b, (g_bin, ch_bin) in enumerate(zip(genome_bins, chr_bins)):
        merged_chr_and_genome_bins[ii_b] = ch_bin + ' / ' + g_bin

    # df = pd.DataFrame(columns=np.concatenate([['fragment', 'name', 'type'], merged_chr_and_genome_bins]))
    df_contacts = pd.DataFrame(columns=merged_chr_and_genome_bins)
    df_frequencies = df_contacts.copy(deep=True)
    nb_bins = len(chr_bins)
    fragments = []
    types = []
    names = []
    for f in bins_contacts_dict:
        contacts = np.zeros(nb_bins, dtype=int)
        fragments.append(f)
        types.append(fragment_infos_dict[f]['type'])
        names.append(fragment_infos_dict[f]['name'])
        for ctc in bins_contacts_dict[f]:
            idx = np.where(chr_bins == ctc)[0]
            contacts[idx] = bins_contacts_dict[f][ctc]
        df_contacts.loc[len(df_contacts)] = contacts
        df_frequencies.loc[len(df_frequencies)] = contacts / np.sum(contacts)

    df_contacts.insert(0, 'fragments', np.asarray(fragments))
    df_contacts.insert(1, 'names', np.asarray(names))
    df_contacts.insert(2, 'types', np.asarray(types))
    df_contacts.to_csv(output_path+'contacts_bins_matrix.csv')

    df_frequencies.insert(0, 'fragments', np.asarray(fragments))
    df_frequencies.insert(1, 'names', np.asarray(names))
    df_frequencies.insert(2, 'types', np.asarray(types))
    df_frequencies.to_csv(output_path+'frequencies_bins_matrix.csv')


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    artificial_genome_path, filtered_contacts_path, output_path = ['' for _ in range(3)]

    try:
        opts, args = getopt.getopt(argv, "hg:c:o:", ["--help",
                                                     "--genome",
                                                     "--contacts",
                                                     "--output"])
    except getopt.GetoptError:
        print('contacts filter arguments :\n'
              '-g <fasta_genome_input> (artificially generated with oligos_replacement.py) \n'
              '-c <filtered_contacts_input.csv> (contacts filtered with contacts_filter.py) \n'
              '-o <output>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('contacts filter arguments :\n'
                  '-g <fasta_genome_input> (artificially generated with oligos_replacement.py) \n'
                  '-c <filtered_contacts_input.csv> (contacts filtered with contacts_filter.py) \n'
                  '-o <output>')
            sys.exit()
        elif opt in ("-g", "--genome"):
            artificial_genome_path = arg
        elif opt in ("-c", "--contacts"):
            filtered_contacts_path = arg
        elif opt in ("-o", "--output"):
            output_path = arg

    genome_bins_names, chr_bins_names = build_bins_from_genome(artificial_genome_path, 10000)
    contacts_dict, infos_dict = get_fragments_dict(contacts_path=filtered_contacts_path, bin_size=10000)
    set_fragments_contacts_bins(bins_contacts_dict=contacts_dict,
                                fragment_infos_dict=infos_dict,
                                chr_bins=chr_bins_names,
                                genome_bins=genome_bins_names,
                                output_path=output_path)


def debug(artificial_genome_path: str,
          filtered_contacts_path: str,
          output_path: str):

    genome_bins_names, chr_bins_names = build_bins_from_genome(artificial_genome_path, 10000)
    contacts_dict, infos_dict = get_fragments_dict(contacts_path=filtered_contacts_path, bin_size=10000)
    set_fragments_contacts_bins(bins_contacts_dict=contacts_dict,
                                fragment_infos_dict=infos_dict,
                                chr_bins=chr_bins_names,
                                genome_bins=genome_bins_names,
                                output_path=output_path)


if __name__ == "__main__":
    if is_debug():
        artificial_genome = "../../../contacts_format/inputs/S288c_DSB_LY_capture_artificial_nicolas.fa"
        filtered_contacts = "../../../contacts_format/inputs/contacts_filtered_nicolas.csv"
        output = "../../../contacts_format/outputs/"
        debug(artificial_genome_path=artificial_genome,
              filtered_contacts_path=filtered_contacts,
              output_path=output)
    else:
        main(sys.argv[1:])

    print('--- DONE ---')
