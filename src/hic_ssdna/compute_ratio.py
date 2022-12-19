#! /usr/bin/env python3

import numpy as np
import pandas as pd
import utils
import sys
import getopt


def fold_over(df_stats: pd.DataFrame):
    """
    Compute a fold over :
        number of contact for one oligo divided by the mean (or median) of all other 'ds' oligos in the genome
    """
    ds_average = np.mean(df_stats.loc[df_stats['types'] == 'ds', 'total_contacts'].values)
    df_stats['fold_over'] = df_stats['total_contacts'] / ds_average
    return df_stats


def compute_stats(formatted_contacts_path: str,
                  cis_range: int,
                  output_path: str):
    """
    After having formatted the contacts of each oligos in the genome (file with bin_size=0)
    We know cant to use this new formatted file to compute basics statistics like total number of contacts,
    frequencies of contacts intra .vs. inter chromosomes, cis .vs. trans oligos
     (with cis range as an input value given by the user, in bp)
    """

    chr_size = {'chr1': 230218, 'chr2': 813184, 'chr3': 316620, 'chr4': 1531933,
                'chr5': 576874, 'chr6': 270161, 'chr7': 1090940, 'chr8': 562643,
                'chr9': 439888, 'chr10': 745751, 'chr11': 666816, 'chr12': 1078177,
                'chr13': 924431, 'chr14': 784333, 'chr15': 1091291, 'chr16': 948066,
                'mitochondrion': 85779}

    genome_size = sum(chr_size.values())
    chr_size_normalized = {k: v/genome_size for k, v in chr_size.items()}

    #   dataframe of the formatted contacts csv file previously created,
    #   with DTYPE=object because multiple type are present in columns
    #   low_memory because df_all contains multiple dtypes within its columns,
    #   so pandas has to avoid to guess their types
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    #   It needs thus to split between numeric and not numeric data
    df_infos, df_contacts = utils.split_formatted_dataframe(df_all)

    #   df_stats : results dataframe containing all the measures and calculus we want.
    df_stats = pd.DataFrame(columns=['names', 'types', 'cis', 'trans', 'intra', 'inter', 'total_contacts'])
    #   df_chr_normalized_freq : frequencies of contact in each chr, for each oligo, normalized
    #   by the length of the chromosome, itself normalized by the length of the genome
    df_chr_normalized_freq = pd.DataFrame(columns=np.concatenate((['names', 'types'], list(chr_size.keys()))))

    #   df_infos is transposed here
    for index, row in df_infos.T.iterrows():
        name, ttype, current_chr, cis_left_boundary, cis_right_boundary = row
        #   set the cis position in downstream of the oligo/read : start position of the read - cis_range
        cis_left_boundary = int(cis_left_boundary) - cis_range
        #   set the cis position in upstream of the oligo/read : end position of the read + cis_range
        cis_right_boundary = int(cis_right_boundary) + cis_range

        #   Sum all contact made by the current oligo to get the total number of contacts
        all_contacts = np.sum(df_contacts.loc[:, index].values)

        #   subset dataframe that only contains the contacts made by the current oligo with the cis region
        sub_df_cis_contacts = df_contacts.loc[(df_contacts['positions'] > cis_left_boundary) &
                                              (cis_right_boundary > df_contacts['positions']) &
                                              (df_contacts['chr'] == current_chr)]

        #   Get the total amount of contacts made in cis region
        cis_contacts = np.sum(sub_df_cis_contacts.loc[:, index])

        #   To get the total amout of contacts made in trans region, just do overall total contacts - cis contacts
        trans_contacts = all_contacts - cis_contacts

        #   subset dataframe that only contains contacts made on the same chr (chr of the  current oligo)
        sub_df_intra_chromosome_contacts = df_contacts.loc[df_contacts['chr'] == current_chr]
        #   get total amount
        intra_chromosome_contacts = \
            np.sum(sub_df_intra_chromosome_contacts.loc[:, index].values)

        #   Compute the frequency of contact made with each chromosome
        #   and normalize this frequency to the size of the chromosome
        #   (itself normalized to the size of the genome)
        tmp_df_chr_normalized_freq = pd.DataFrame({'names': [name], 'types': [ttype]})
        for chrom, nrm_size in chr_size_normalized.items():
            if chrom == current_chr:
                continue
            sub_df_chrom_nrm = df_contacts.loc[df_contacts['chr'] == chrom]
            tmp_df_chr_normalized_freq[chrom] = \
                [(np.sum(sub_df_chrom_nrm.loc[:, index].values) / all_contacts) / nrm_size]

        #   Concatenated temporary dataframe with the result one
        df_chr_normalized_freq = pd.concat([df_chr_normalized_freq, tmp_df_chr_normalized_freq])

        #   subset dataframe that only contains contacts made on all other chr
        #   (other than the one of the current oligo)
        sub_df_inter_chromosome_contacts = df_contacts.loc[df_contacts['chr'] != current_chr]
        inter_chromosome_contacts = \
            np.sum(sub_df_inter_chromosome_contacts.loc[:, index].values)

        #   transform the contact number to frequencies by dividing
        #   them by the overall number of contacts of the current oligo
        cis_freq = cis_contacts / all_contacts
        trans_freq = trans_contacts / all_contacts
        inter_chromosome_freq = inter_chromosome_contacts / all_contacts
        intra_chromosome_freq = intra_chromosome_contacts / all_contacts

        #   temporary dataframe contains all the stats made for the current oligo/read in the loop
        #   and its information such as probe's name, type (ss, ds_neg etc ...) and inner location.
        tmp_df_stats = pd.DataFrame({'names': [name], 'types': [ttype], 'cis': [cis_freq],
                                     'trans': [trans_freq], 'intra': [intra_chromosome_freq],
                                     'inter': [inter_chromosome_freq], 'total_contacts': [all_contacts]})

        #   Concatenated temporary dataframe with the result one
        df_stats = pd.concat([df_stats, tmp_df_stats])
        #   Restore the index as a range of 0 to len(df)
        df_stats.index = range(len(df_stats))
        df_chr_normalized_freq.index = range(len(df_chr_normalized_freq))

    df_stats = fold_over(df_stats)
    df_stats.to_csv(output_path + 'global_statistics.tsv', sep='\t')
    df_chr_normalized_freq.to_csv(output_path + 'normalized_chr_frequencies.tsv', sep='\t')


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    cis_range, oligos_path, hic_contacts_list_path, formatted_contacts_path, output_path = ['' for _ in range(5)]
    try:
        opts, args = getopt.getopt(argv, "hc:r:O:", ["--help",
                                                     "--contacts",
                                                     "--cis_range",
                                                     "--output"])
    except getopt.GetoptError:
        print('compute ratios arguments :\n'
              '-c <formated_contacts.csv> (contacts filtered with contacts_format.py) \n'
              '-r <bin_size> (size of a bin, in bp) \n'
              '-O <output_file_name.csv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('compute ratios arguments :\n'
                  '-c <formated_contacts.csv> (contacts filtered with contacts_format.py) \n'
                  '-r <bin_size> (size of a bin, in bp) \n'
                  '-O <output_file_name.csv>')
            sys.exit()
        elif opt in ("-c", "--contacts"):
            formatted_contacts_path = arg
        elif opt in ("-r", "--cis_range"):
            cis_range = int(arg)
        elif opt in ("-O", "--output"):
            output_path = arg.split('contacts_matrix.tsv')[0]

    compute_stats(cis_range=cis_range,
                  formatted_contacts_path=formatted_contacts_path,
                  output_path=output_path)


def debug(cis_range: int,
          formatted_contacts_path: str,
          output_path: str):

    compute_stats(cis_range=cis_range,
                  formatted_contacts_path=formatted_contacts_path,
                  output_path=output_path)


if __name__ == "__main__":
    #   Go into debug function if debug mode is detected, else go for main script with sys arguments
    if utils.is_debug():
        #   Debug is mainly used for testing function of the script
        #   Parameters have to be declared here
        all_contacted_pos = "../../../bash_scripts/compute_ratio/inputs/" \
                            "AD162_test.csv"
        output = "../../../bash_scripts/compute_ratio/outputs/fragments_percentages.tsv"
        cis_range_value = 50000
        debug(cis_range=cis_range_value,
              formatted_contacts_path=all_contacted_pos,
              output_path=output)
    else:
        main()

    print('--- DONE ---')
