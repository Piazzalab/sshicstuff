#! /usr/bin/env python3

import numpy as np
import pandas as pd
import utils
import sys
import getopt


def fold_over(df_stats: pd.DataFrame):
    ds_average = np.mean(df_stats.loc[df_stats['types'] == 'ds', 'total_contacts'].values)
    df_stats['fold_over'] = df_stats['total_contacts'] / ds_average
    return df_stats


def compute_stats(formated_contacts_path: str,
                  cis_range: int,
                  output_path: str):
    #   low_memory because df_all contains multiple dtypes within its columns,
    #   so pandas has to avoid to guess their types
    df_all = pd.read_csv(formated_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_infos, df_contacts = utils.split_formatted_dataframe(df_all)
    df_stats = pd.DataFrame(columns=['names', 'types', 'cis', 'trans', 'intra', 'inter', 'total_contacts'])

    for index, row in df_infos.T.iterrows():
        name, ttype, current_chr, cis_left_boundary, cis_right_boundary = row
        cis_left_boundary = int(cis_left_boundary) - cis_range
        cis_right_boundary = int(cis_right_boundary) + cis_range

        all_contacts = np.sum(df_contacts.loc[:, index].values)

        sub_df_cis_contacts = df_contacts.loc[(df_contacts['positions'] > cis_left_boundary) &
                                              (cis_right_boundary > df_contacts['positions']) &
                                              (df_contacts['chr'] == current_chr)]

        cis_contacts = np.sum(sub_df_cis_contacts.loc[:, index])
        trans_contacts = all_contacts - cis_contacts

        sub_df_intra_chromosome_contacts = df_contacts.loc[df_contacts['chr'] == current_chr]
        intra_chromosome_contacts = \
            np.sum(sub_df_intra_chromosome_contacts.loc[:, index].values)
        sub_df_inter_chromosome_contacts = df_contacts.loc[df_contacts['chr'] != current_chr]
        inter_chromosome_contacts = \
            np.sum(sub_df_inter_chromosome_contacts.loc[:, index].values)

        cis_freq = cis_contacts / all_contacts
        trans_freq = trans_contacts / all_contacts
        inter_chromosome_freq = inter_chromosome_contacts / all_contacts
        intra_chromosome_freq = intra_chromosome_contacts / all_contacts

        tmp_df_stats = pd.DataFrame({'names': [name], 'types': [ttype], 'cis': [cis_freq],
                                     'trans': [trans_freq], 'intra': [intra_chromosome_freq],
                                     'inter': [inter_chromosome_freq], 'total_contacts': [all_contacts]})

        df_stats = pd.concat([df_stats, tmp_df_stats])
        df_stats.index = range(len(df_stats))

    df_stats = fold_over(df_stats)
    df_stats.to_csv(output_path + 'statistics.csv')


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    cis_range, oligos_path, hic_contacts_list_path, formated_contacts_path, output_path = ['' for _ in range(5)]
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
            formated_contacts_path = arg
        elif opt in ("-r", "--cis_range"):
            cis_range = int(arg)
        elif opt in ("-O", "--output"):
            output_path = arg.split('contacts_matrix.csv')[0]

    compute_stats(cis_range=cis_range,
                  formated_contacts_path=formated_contacts_path,
                  output_path=output_path)


def debug(cis_range: int,
          formated_contacts_path: str,
          output_path: str):

    compute_stats(cis_range=cis_range,
                  formated_contacts_path=formated_contacts_path,
                  output_path=output_path)


if __name__ == "__main__":
    if utils.is_debug():
        all_contacted_pos = "../../../bash_scripts/compute_ratio/inputs/" \
                            "AD162_test.csv"
        output = "../../../bash_scripts/compute_ratio/outputs/fragments_percentages.csv"
        cis_range_value = 50000
        debug(cis_range=cis_range_value,
              formated_contacts_path=all_contacted_pos,
              output_path=output)
    else:
        main()

    print('--- DONE ---')
