#! /usr/bin/env python3

import re
import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import FastaIterator
from collections import Counter


def build_bins_from_genome(path_to_genome: str,
                           bin_size: int):
    """
    This function aims to parse the genome, and each chromosome into bins (regular range of bp)
    of size 'bin_size' a parameter given by the user as input.

    For the chr_bins, they will start at chrX_0, chrX_(0+bin_size) ... chrX_end
        idem for chrY, it starts again at chrX_0, and so on chrX_(0+bin_size) ... chrX_end
    For the genome they will start at 0, 0+bin_size, 0+2*bin_size ... genome_end.
        For the genome_bins, the positions are not reset when we go from the end of chrX to the start of chrY
    """
    genome = open(path_to_genome, "r")
    nb_bins_per_chr: dict = {}
    chr_length: dict = {}
    for record in FastaIterator(genome):
        chr_id = record.id
        nb_bins_per_chr[chr_id] = len(str(record.seq)) // bin_size + 1
        chr_length[chr_id] = len(record)
    genome.close()
    total_nb_bins = np.sum(list(nb_bins_per_chr.values()))
    #   Results are stored into a dictionary
    bed_array: dict = {'chr': np.zeros(total_nb_bins, dtype='<U16'),
                       'chr_and_bins': np.zeros(total_nb_bins, dtype='<U64'),
                       'bins_per_chr': {k: np.zeros(v, dtype=int) for k, v in nb_bins_per_chr.items()},
                       'bins_in_genome': np.zeros(total_nb_bins, dtype=int),
                       }

    #   Two counter are set as repair for creating bins
    #       big_counter is used as repair especially for genome_bins
    big_counter = 0
    counter = 0
    for ii_chr, nb_bins in nb_bins_per_chr.items():
        micro_counter = 0
        for ii_bin in range(0, nb_bins, 1):
            start = ii_bin * bin_size
            bed_array['chr'][counter] = ii_chr
            bed_array['chr_and_bins'][counter] = ii_chr + '_' + str(start)
            bed_array['bins_per_chr'][ii_chr][micro_counter] = start
            bed_array['bins_in_genome'][counter] = start + big_counter * bin_size
            micro_counter += 1
            counter += 1
        big_counter += nb_bins

    bed_array['all_bins_in_chr'] = np.concatenate(list(bed_array['bins_per_chr'].values()))
    return bed_array


def get_contacts_per_bin(
        genome_bins_dict: dict,
        formatted_contacts_path: str,
        bin_size: int,
        output_path: str
        ):

    unique_chr = pd.unique(genome_bins_dict['chr'])
    nb_bins = len(genome_bins_dict['bins_in_genome'])
    df_formatted_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_binned_contacts = pd.DataFrame({'chr': genome_bins_dict['chr'],
                                       'chr_bins': genome_bins_dict['all_bins_in_chr'],
                                       'genome_bins': genome_bins_dict['bins_in_genome']})

    df_binned_frequencies = df_binned_contacts.copy(deep=True)

    fragments = np.array([f for f in df_formatted_contacts.columns.values if re.match(r'\d+', f)])
    for frag in fragments:
        df_binned_contacts[frag] = np.zeros(nb_bins, dtype=int)
        for chrom in unique_chr:
            current_chr_bins_count = {k: 0 for k in genome_bins_dict['bins_per_chr'][chrom]}
            sub_formatted_df = df_formatted_contacts.loc[df_formatted_contacts['chr'] == chrom]
            contacted_positions = (sub_formatted_df['positions'][sub_formatted_df[frag] != 0] // bin_size) * bin_size
            if len(contacted_positions) > 0:
                contacts_count = dict(Counter(contacted_positions))
                for b, count in contacts_count.items():
                    current_chr_bins_count[b] = count

            df_binned_contacts.loc[df_binned_contacts['chr'] == chrom, frag] = \
                list(current_chr_bins_count.values())

        df_binned_frequencies[frag] = df_binned_contacts[frag] / sum(df_binned_contacts[frag])

        df_binned_contacts.to_csv(output_path + '_' + str(bin_size // 1000) + 'kb_contacts.tsv', sep='\t')
        df_binned_frequencies.to_csv(output_path + '_' + str(bin_size // 1000) + 'kb_frequencies.tsv', sep='\t')


def run(
        artificial_genome_path: str,
        formatted_contacts_path: str,
        bin_size: int,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    bed_bins_dict = build_bins_from_genome(
        artificial_genome_path,
        bin_size=bin_size)

    get_contacts_per_bin(formatted_contacts_path=formatted_contacts_path,
                         genome_bins_dict=bed_bins_dict,
                         bin_size=bin_size,
                         output_path=output_dir+sample_id)

    print('DONE: ', sample_id)
