#! /usr/bin/env python3

import numpy as np
import pandas as pd
import math
import sys
import getopt
from Bio.SeqIO.FastaIO import FastaIterator
import utils


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
    for record in FastaIterator(genome):
        chr_id = record.id
        nb_bins_per_chr[chr_id] = len(str(record.seq)) // bin_size + 1
    genome.close()
    total_nb_bins = np.sum(list(nb_bins_per_chr.values()))
    #   Results are stored into a dictionary
    bed_array: dict = {'chr': np.zeros(total_nb_bins, dtype='<U16'),
                       'chr_and_bins': np.zeros(total_nb_bins, dtype='<U64'),
                       'bins_in_chr': np.zeros(total_nb_bins, dtype='<U64'),
                       'bins_in_genome': np.zeros(total_nb_bins, dtype='<U64'),
                       }

    #   Two counter are set as repair for creating bins
    #       big_counter is used as repair especially for genome_bins
    big_counter = 0
    counter = 0
    for ii_chr, nb_bins in nb_bins_per_chr.items():
        for ii_bin in range(0, nb_bins, 1):
            start = ii_bin * bin_size
            bed_array['chr'][counter] = ii_chr
            bed_array['chr_and_bins'][counter] = ii_chr + '_' + str(start)
            bed_array['bins_in_chr'][counter] = start
            bed_array['bins_in_genome'][counter] = start + big_counter * bin_size
            counter += 1
        big_counter += nb_bins
    return bed_array


def count_occurrences(df: pd.DataFrame,
                      x: str,
                      contacts_res: dict,
                      infos_res: dict,
                      all_contacted_pos: dict,
                      bin_size):
    """
    This function will count the number of contacts for each read that comes from column 'frag_x' in each bin.
    The results are stored in three dictionaries, given as arguments.
        contacts_res of the form :  {oligoX : {chrX_binA : n ... } ...}
        infos_res of the form : {oligoX : {name: "probeX", type: "ds" , chr: 'chr1, ...} ...}
        all_contacted_pos of the form : {oligoX : {chrX_1456 : n, chrY_89445: m ... } ...}
    """
    #   if x = a get b, if x = b get a
    y = utils.frag2(x)
    for ii_f, f in enumerate(df['frag_' + x].values):
        if not pd.isna(df['name_' + x][ii_f]):
            if f not in infos_res:
                #   Nota Bene : A same read or fragment may contain two oligos because they are very close
                #       Thus for a same read we acn see two probe's names that are different
                #       For the moment the infos_res[f][names] is an array that may contain one (in most cases)
                #       or two probes (for two oligos in one read).
                infos_res[f] = {'type': df['type_' + x][ii_f],
                                'names': [df['name_' + x][ii_f]],
                                'chr': df['chr_' + x][ii_f],
                                'start': df['start_' + x][ii_f],
                                'end': df['end_' + x][ii_f]}

            elif f in infos_res:
                if df['name_' + x][ii_f] not in infos_res[f]['names']:
                    infos_res[f]['names'].append(df['name_' + x][ii_f])

            chr_id = df['chr_' + y][ii_f]
            start = df['start_' + y][ii_f]
            chr_and_pos = chr_id + '_' + str(start)

            if f not in contacts_res:
                contacts_res[f] = {}

            if chr_id not in all_contacted_pos:
                all_contacted_pos[chr_id] = set()

            if chr_and_pos not in all_contacted_pos[chr_id]:
                all_contacted_pos[chr_id].add(start)

            if bin_size > 0:
                start = math.floor(start / bin_size) * bin_size
                bin_id = chr_id + '_' + str(start)
            else:
                bin_id = chr_and_pos

            if bin_id not in contacts_res[f]:
                contacts_res[f][bin_id] = df['contacts'][ii_f]
            else:
                contacts_res[f][bin_id] += df['contacts'][ii_f]

    for f in infos_res:
        if len(infos_res[f]['names']) > 1:
            #   If the current fragment or read in the loop has multiple names i.e., contains two oligos
            #   we decided to merge the probe's names in a single one : 'name1-/-name2'
            infos_res[f]['uid'] = '-/-'.join(infos_res[f]['names'])
        else:
            infos_res[f]['uid'] = infos_res[f]['names'][0]

    return contacts_res, infos_res, all_contacted_pos


def get_fragments_dict(contacts_path: str,
                       bin_size: int):
    """
    From a filtered contacts file, count the number of contacts made for each read in every bin of size 'bin_size'

    """
    #   dataframe of the filtered contacts file created previously (contacts_filter)
    df_contacts_filtered = pd.read_csv(contacts_path, sep=',')
    #   dictionary that will store for each fragment (aka read) , for each bin, the number of contacts
    #       It will be in the form : {oligoX : {chrX_binA : n ... } ...}
    fragments_contacts: dict = {}
    #   dictionary that will store for each fragment (aka read) , its information
    #       It will be in the form : {oligoX : {name: "probeX", type: "ds" , chr: 'chr1, ...} ...}
    fragments_infos: dict = {}
    #   dictionary that will store for each fragment (aka read) , all the region it contacted, so NOT BINNED
    #       It will be in the form : {oligoX : {chrX_1456 : n, chrY_89445: m ... } ...}
    all_contacted_chr_pos: dict = {}

    #   At this moment these three dict are empty.
    #   We have to fill them with the function count_occurrences.

    #   There a first passage for reads that come from the 'frag_a' column in the filtered_contacts dataframe
    fragments_contacts, fragments_infos, all_contacted_chr_pos = \
        count_occurrences(df=df_contacts_filtered,
                          x='a',
                          contacts_res=fragments_contacts,
                          infos_res=fragments_infos,
                          all_contacted_pos=all_contacted_chr_pos,
                          bin_size=bin_size)

    #   Then a second passage for reads that come from the 'frag_b' column in the filtered_contacts dataframe
    fragments_contacts, fragments_infos, all_contacted_chr_pos = \
        count_occurrences(df=df_contacts_filtered,
                          x='b',
                          contacts_res=fragments_contacts,
                          infos_res=fragments_infos,
                          all_contacted_pos=all_contacted_chr_pos,
                          bin_size=bin_size)

    return fragments_contacts, fragments_infos, all_contacted_chr_pos


def concatenate_infos_and_contacts(df1: pd.DataFrame,
                                   df2: pd.DataFrame,
                                   headers: list):
    """
    concatenate dataframe DF2 containing information such as probe's name, type,
    location etc ... with the dataframe containing the amount of occurrence per bin DF1
    """
    #   transpose df2
    df2 = df2.T
    df2.columns = headers
    df3 = pd.concat([df2, df1])
    return df3


def set_fragments_contacts_bins(bed_bins: dict,
                                bins_contacts_dict: dict,
                                fragment_infos_dict: dict,
                                output_path: str):

    """
    Write a dataframe with contacts per bin, for a specified bin_size for each oligo.
    Formatted dataframe of the form :
    **********************************************************************************************************
                 chr    chr_bins    genome_bins      2630                5315                   8579
    names                                            Neg_chr2-199707     chr2-650593-MET8       chr4-64420-CDC13
    types                                            ds_neg              ds                     ds
    self_chr                                         chr2                chr2                   chr4
    self_start                                       199707	             650550	                64397
    self_end                                         199746	             650756	                64699
    0            1         0           0             0	                 0	                    0
    1            1         20000       20000         0                   3                      0
    2            1         40000       40000         1                   0                      2
    3            1         60000       60000         0                   0                      1
    |            |         |           |             |                   |                      |
    |            |         |           |             |                   |                      |
    36193    chr_art       0           1227000       2                   0                      0
    **********************************************************************************************************
    """
    chromosomes = bed_bins['chr']
    chr_and_bins = bed_bins['chr_and_bins']
    bins_in_chr = bed_bins['bins_in_chr']
    bins_in_genome = bed_bins['bins_in_genome']
    df_contc = pd.DataFrame({'chr': chromosomes, 'chr_bins': bins_in_chr, 'genome_bins': bins_in_genome})
    df_freq = pd.DataFrame({'chr': chromosomes, 'chr_bins': bins_in_chr, 'genome_bins': bins_in_genome})

    nb_bins = len(bins_in_genome)

    headers = list(df_contc.columns.values)
    fragments = []
    names = ['', '', '']
    types = ['', '', '']
    chrs = ['', '', '']
    starts = ['', '', '']
    ends = ['', '', '']
    for f in bins_contacts_dict:
        contacts = np.zeros(nb_bins, dtype=int)
        fragments.append(f)
        types.append(fragment_infos_dict[f]['type'])
        names.append(fragment_infos_dict[f]['uid'])
        chrs.append(fragment_infos_dict[f]['chr'])
        starts.append(fragment_infos_dict[f]['start'])
        ends.append(fragment_infos_dict[f]['end'])

        for _bin in bins_contacts_dict[f]:
            idx = np.where(chr_and_bins == _bin)[0]
            contacts[idx] = bins_contacts_dict[f][_bin]
        df_contc[f] = contacts
        df_freq[f] = contacts / np.sum(contacts)

    headers.extend(fragments)
    df_infos = pd.DataFrame(
        {'names': names, 'types': types, 'self_chr': chrs, 'self_start': starts, 'self_end': ends})
    df_contc = concatenate_infos_and_contacts(df1=df_contc, df2=df_infos, headers=headers)
    df_freq = concatenate_infos_and_contacts(df1=df_freq, df2=df_infos, headers=headers)
    df_contc.to_csv(output_path + '_contacts_matrix.tsv', sep='\t')
    df_freq.to_csv(output_path + '_frequencies_matrix.tsv', sep='\t')


def set_fragments_contacts_no_bin(contacts_pos_dict: dict,
                                  fragment_infos_dict: dict,
                                  all_chr_pos: dict,
                                  output_path: str):

    """
    Write a dataframe with all the contacts NOT BINNED here made for each read.
        Formatted dataframe of the form :
    **********************************************************************************************************
                 chr   positions     2630                    5315                    8579
    names                            Neg_chr2-199707        chr2-650593-MET8        chr4-64420-CDC13
    types                            ds_neg                 ds                      ds
    self_chr                         chr2                   chr2                    chr4
    self_start                       199707	                650550	                64397
    self_end                         199746	                650756	                64699
    0            1     0             0	                    0	                    0
    1            1     509           0                      3                       0
    2            1     1149          1                      0                       2
    3            1     1492          0                      0                       1
    |            |      |            |                      |                       |
    |            |      |            |                      |                       |
    36193    chr_art    7224          2                     0                       0
    **********************************************************************************************************
    """
    chr_and_pos = []
    chromosomes = []
    positions = []

    for chr_id, pos_list in all_chr_pos.items():
        all_chr_pos[chr_id] = sorted(pos_list)

    chr_unique_list = np.concatenate([['chr' + str(i) for i in range(1, 17)],
                                      ['2_micron', 'mitochondrion', 'chr_artificial']])

    for chr_id in chr_unique_list:
        if chr_id in all_chr_pos:
            new_pos = all_chr_pos[chr_id]
            positions.extend(new_pos)
            chromosomes.extend(np.repeat(chr_id, len(new_pos)))
            for npos in new_pos:
                chr_and_pos.append(chr_id + '_' + str(npos))

    chr_and_pos = np.asarray(chr_and_pos)
    df_contc = pd.DataFrame({'chr': chromosomes, 'positions': positions})
    df_freq = pd.DataFrame({'chr': chromosomes, 'positions': positions})

    headers = list(df_contc.columns.values)
    fragments = []
    names = ['', '']
    types = ['', '']
    chrs = ['', '']
    starts = ['', '']
    ends = ['', '']
    for f in contacts_pos_dict:
        contacts = np.zeros(len(chr_and_pos), dtype=int)
        fragments.append(f)
        types.append(fragment_infos_dict[f]['type'])
        names.append(fragment_infos_dict[f]['uid'])
        chrs.append(fragment_infos_dict[f]['chr'])
        starts.append(fragment_infos_dict[f]['start'])
        ends.append(fragment_infos_dict[f]['end'])

        for pos in contacts_pos_dict[f]:
            idx = np.argwhere(chr_and_pos == pos)[0]
            contacts[idx] = contacts_pos_dict[f][pos]

        df_contc[f] = contacts
        df_freq[f] = contacts / np.sum(contacts)

    headers.extend(fragments)
    df_infos = pd.DataFrame(
        {'names': names, 'types': types, 'self_chr': chrs, 'self_start': starts, 'self_end': ends})
    df_contc = concatenate_infos_and_contacts(df1=df_contc, df2=df_infos, headers=headers)
    df_freq = concatenate_infos_and_contacts(df1=df_freq, df2=df_infos, headers=headers)
    df_contc.to_csv(output_path + '_formatted_contacts_matrix.tsv', sep='\t')
    df_freq.to_csv(output_path + '_formatted_frequencies_matrix.tsv', sep='\t')


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
                                                       "--bin_size",
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
            output_path = arg.split('-filtered.csv')[0]

    bin_size = int(bin_size)
    contacts_dict, infos_dict, all_contacted_pos = get_fragments_dict(contacts_path=filtered_contacts_path,
                                                                      bin_size=bin_size)
    if bin_size > 0:
        bed_pos = build_bins_from_genome(artificial_genome_path, bin_size=bin_size)
        set_fragments_contacts_bins(bed_bins=bed_pos,
                                    bins_contacts_dict=contacts_dict,
                                    fragment_infos_dict=infos_dict,
                                    output_path=output_path)
    else:
        set_fragments_contacts_no_bin(contacts_pos_dict=contacts_dict,
                                      fragment_infos_dict=infos_dict,
                                      all_chr_pos=all_contacted_pos,
                                      output_path=output_path)


def debug(artificial_genome_path: str,
          filtered_contacts_path: str,
          bin_size: int,
          output_path: str):
    contacts_dict, infos_dict, all_contacted_pos = get_fragments_dict(contacts_path=filtered_contacts_path,
                                                                      bin_size=bin_size)
    if bin_size > 0:
        bed_pos = build_bins_from_genome(artificial_genome_path, bin_size=bin_size)
        set_fragments_contacts_bins(bed_bins=bed_pos,
                                    bins_contacts_dict=contacts_dict,
                                    fragment_infos_dict=infos_dict,
                                    output_path=output_path)
    else:
        set_fragments_contacts_no_bin(contacts_pos_dict=contacts_dict,
                                      fragment_infos_dict=infos_dict,
                                      all_chr_pos=all_contacted_pos,
                                      output_path=output_path)


if __name__ == "__main__":
    #   Go into debug function if debug mode is detected, else go for main script with sys arguments
    if utils.is_debug():
        #   Debug is mainly used for testing function of the script
        #   Parameters have to be declared here
        artificial_genome = "../../../bash_scripts/contacts_format/inputs/S288c_DSB_LY_capture_artificial.fa"
        filtered_contacts = "../../../bash_scripts/contacts_format/inputs/contacts_filtered_nicolas.csv"
        output = "../../../bash_scripts/contacts_format/outputs/frequencies_per_bin_matrix.csv"
        bin_size_value = 100000
        debug(artificial_genome_path=artificial_genome,
              filtered_contacts_path=filtered_contacts,
              bin_size=bin_size_value,
              output_path=output)
    else:
        main(sys.argv[1:])

    print('--- DONE ---')
