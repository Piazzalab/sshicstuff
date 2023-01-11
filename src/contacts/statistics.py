#! /usr/bin/env python3

import re
import numpy as np
import pandas as pd


def fold_over(df_stats: pd.DataFrame):
    """
    Compute a fold over :
        number of contact for one oligo divided by the mean (or median) of all other 'ds' oligos in the genome
    """
    ds_average = np.mean(df_stats.loc[df_stats['types'] == 'ds', 'total_contacts'].values)
    df_stats['fold_over'] = df_stats['total_contacts'] / ds_average
    return df_stats


def compute_stats(formatted_contacts_path: str,
                  fragments_to_oligos_path: str,
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
                'mitochondrion': 85779, '2_micron': 6318}

    genome_size = sum(chr_size.values())
    chr_size_normalized = {k: v/genome_size for k, v in chr_size.items()}


    #   dataframe of the formatted contacts tsv file previously created
    df_formatted_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_info = pd.read_csv(fragments_to_oligos_path, sep='\t', index_col=0)

    #   df_stats : results dataframe containing all the measures and calculus we want.
    df_stats = pd.DataFrame(columns=['names', 'types', 'cis', 'trans', 'intra', 'inter', 'total_contacts'])
    #   df_chr_normalized_freq : frequencies of contact in each chr, for each oligo, normalized
    #   by the length of the chromosome, itself normalized by the length of the genome
    df_chr_normalized_freq = pd.DataFrame(columns=np.concatenate((['names', 'types'], list(chr_size.keys()))))

    #   df_infos is transposed here
    for index, row in df_info.T.iterrows():
        name, ttype, current_chr, cis_left_boundary, cis_right_boundary = row
        #   set the cis position in downstream of the oligo/read : start position of the read - cis_range
        cis_left_boundary = int(cis_left_boundary) - cis_range
        #   set the cis position in upstream of the oligo/read : end position of the read + cis_range
        cis_right_boundary = int(cis_right_boundary) + cis_range

        #   Sum all contact made by the current oligo to get the total number of contacts
        all_contacts = np.sum(df_formatted_contacts.loc[:, index].values)

        #   subset dataframe that only contains the contacts made by the current oligo with the cis region
        sub_df_cis_contacts = df_formatted_contacts.loc[
            (df_formatted_contacts['positions'] > cis_left_boundary) &
            (cis_right_boundary > df_formatted_contacts['positions']) &
            (df_formatted_contacts['chr'] == current_chr)]

        #   Get the total amount of contacts made in cis region
        cis_contacts = np.sum(sub_df_cis_contacts.loc[:, index])

        #   To get the total amout of contacts made in trans region, just do overall total contacts - cis contacts
        trans_contacts = all_contacts - cis_contacts

        #   subset dataframe that only contains contacts made on the same chr (chr of the  current oligo)
        sub_df_intra_chromosome_contacts = df_formatted_contacts.loc[df_formatted_contacts['chr'] == current_chr]
        #   get total amount
        intra_chromosome_contacts = \
            np.sum(sub_df_intra_chromosome_contacts.loc[:, index].values)

        #   Compute the frequency of contact made with each chromosome
        #   and normalize this frequency to the size of the chromosome
        #   (itself normalized to the size of the genome)
        tmp_df_chr_normalized_freq = pd.DataFrame({'names': [name], 'types': [ttype]})
        for chrom, nrm_size in chr_size_normalized.items():
            sub_df_chrom_nrm = df_formatted_contacts.loc[df_formatted_contacts['chr'] == chrom]
            tmp_df_chr_normalized_freq[chrom] = \
                [(np.sum(sub_df_chrom_nrm.loc[:, index].values) / all_contacts) / nrm_size]

        #   Concatenated temporary dataframe with the result one
        df_chr_normalized_freq = pd.concat([df_chr_normalized_freq, tmp_df_chr_normalized_freq])

        #   subset dataframe that only contains contacts made on all other chr
        #   (other than the one of the current oligo)
        sub_df_inter_chromosome_contacts = df_formatted_contacts.loc[df_formatted_contacts['chr'] != current_chr]
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
    df_stats.to_csv(output_path + '_global_statistics.tsv', sep='\t')
    df_chr_normalized_freq.to_csv(output_path + '_normalized_chr_freq.tsv', sep='\t')


def run(
        cis_range: int,
        formatted_contacts_path: str,
        fragments_to_oligos_path: str,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_path = output_dir + sample_id
    compute_stats(cis_range=cis_range,
                  formatted_contacts_path=formatted_contacts_path,
                  fragments_to_oligos_path=fragments_to_oligos_path,
                  output_path=output_path)

    print('DONE: ', sample_id)

