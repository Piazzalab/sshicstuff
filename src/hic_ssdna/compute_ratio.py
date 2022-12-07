#! /usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import getopt


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


def oligos_correction(oligos_path):
    df_oligos = pd.read_csv(oligos_path, sep=",")
    df_oligos.columns = [df_oligos.columns[i].lower() for i in range(len(df_oligos.columns))]
    df_oligos.sort_values(by=['chr', 'start'], inplace=True)
    df_oligos.reset_index(drop=True, inplace=True)
    return df_oligos


def compute_stats(oligos_path: str,
                  formated_contacts_path: str,
                  cis_range: int,
                  output_path: str):
    df_contacts = pd.read_csv(formated_contacts_path, sep=',', index_col=0)
    df_oligos = oligos_correction(oligos_path)
    df_stats = pd.DataFrame(columns=['names', 'types', 'cis', 'trans', 'intra', 'inter'])
    for index, row in df_oligos.iterrows():
        name = row[4]
        ttype = row[3]
        cis_left_boundary = row[1] - cis_range
        cis_right_boundary = row[2] + cis_range
        current_chr = row[0]

        col_idx = 0
        for idx, colname in enumerate(df_contacts.columns.values):
            if name in colname.split('--'):
                col_idx = idx
                break

        if col_idx == 0:
            continue

        cis_contacts = df_contacts.loc[(df_contacts['positions'] > cis_left_boundary) &
                                       (cis_right_boundary > df_contacts['positions']) &
                                       (df_contacts['chr'] == current_chr)]

        cis_frequency = np.sum(cis_contacts.iloc[:, col_idx].values)
        trans_frequency = 1 - cis_frequency

        intra_chromosome_contacts = df_contacts.loc[df_contacts['chr'] == current_chr]
        intra_chromosome_frequency = np.sum(intra_chromosome_contacts.iloc[:, col_idx].values)
        inter_chromosome_contacts = df_contacts.loc[df_contacts['chr'] != current_chr]
        inter_chromosome_frequency = np.sum(inter_chromosome_contacts.iloc[:, col_idx].values)

        tmp_df_stats = pd.DataFrame({'names': [name], 'types': [ttype],
                                     'cis': [cis_frequency],
                                     'trans': [trans_frequency],
                                     'intra': [intra_chromosome_frequency],
                                     'inter': [inter_chromosome_frequency]})

        df_stats = pd.concat([df_stats, tmp_df_stats])
        df_stats.index = range(len(df_stats))
    df_stats.to_csv(output_path + 'percentage_cis_trans_intra_inter.csv')


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    cis_range, oligos_path, formated_contacts_path, output_path = ['' for _ in range(4)]
    try:
        opts, args = getopt.getopt(argv, "ho:c:r:O:", ["--help",
                                                       "--oligos",
                                                       "--contacts",
                                                       "--cis_range"
                                                       "--output"])
    except getopt.GetoptError:
        print('contacts filter arguments :\n'
              '-o <oligos_input.csv> \n'
              '-c <formated_contacts.csv> (contacts filtered with contacts_format.py) \n'
              '-r <bin_size> (size of a bin, in bp) \n'
              '-O <output_file_name.csv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('contacts filter arguments :\n'
                  '-o <oligos_input.csv> \n'
                  '-c <formated_contacts.csv> (contacts filtered with contacts_format.py) \n'
                  '-r <bin_size> (size of a bin, in bp) \n'
                  '-O <output_file_name.csv>')
            sys.exit()
        elif opt in ("-o", "--oligos"):
            oligos_path = arg
        elif opt in ("-c", "--contacts"):
            formated_contacts_path = arg
        elif opt in ("-r", "--cis_range"):
            cis_range = int(arg)
        elif opt in ("-O", "--output"):
            output_path = arg.split('frequencies_matrix.csv')[0]

    compute_stats(cis_range=cis_range,
                  oligos_path=oligos_path,
                  formated_contacts_path=formated_contacts_path,
                  output_path=output_path)


def debug(cis_range: int,
          oligos_path: str,
          formated_contacts_path: str,
          output_path: str):

    compute_stats(cis_range=cis_range,
                  oligos_path=oligos_path,
                  formated_contacts_path=formated_contacts_path,
                  output_path=output_path)


if __name__ == "__main__":
    if is_debug():
        oligos = "../../../compute_ratio/inputs/capture_oligo_positions.csv"
        all_contacted_pos = "../../../compute_ratio/inputs/fragments_frequencies_no_bin_matrix.csv"
        output = "../../../compute_ratio/outputs/fragments_percentages.csv"
        cis_range_value = 50000
        debug(cis_range=cis_range_value,
              oligos_path=oligos,
              formated_contacts_path=all_contacted_pos,
              output_path=output)
    else:
        main()

    print('--- DONE ---')
