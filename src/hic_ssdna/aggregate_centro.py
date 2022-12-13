#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
import sys
import os
import getopt
import utils

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def contacts_focus_around_centromeres(formatted_contacts_path: str,
                                      window_size: int,
                                      centros_infos_path: str):

    df_centros = pd.read_csv(centros_infos_path, sep='\t', index_col=None)
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_infos, df_contacts = utils.split_formatted_dataframe(df_all)

    df_res = pd.DataFrame()

    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]

    for index, row in df_centros.iterrows():
        current_chr = row[0]
        current_centros_pos = row[2]

        left_cutoff = current_centros_pos - window_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_centros_pos + window_size
        tmp_df = df_contacts.loc[(df_contacts['chr'] == current_chr) &
                                 (df_contacts['chr_bins'] > left_cutoff) &
                                 (df_contacts['chr_bins'] < right_cutoff)]

        tmp_df.index = range(len(tmp_df))
        current_centros_bin = utils.find_nearest(tmp_df['chr_bins'].values, current_centros_pos, mode='lower')

        for index2, row2 in tmp_df.iterrows():
            tmp_df.iloc[index2, 1] += bin_size/2 - current_centros_bin

        for c in tmp_df.columns[3:]:
            self_chr = df_infos.loc['self_chr', c]
            if self_chr == current_chr:
                tmp_df.loc[0:len(tmp_df), c] = np.nan

        df_res = pd.concat([df_res, tmp_df])
    df_res.index = range(len(df_res))
    return df_res, df_infos


def compute_mean_per_fragment(df_centros_bins: pd.DataFrame, output_file: str):
    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    bins_counter = dict(Counter(df_centros_bins['chr_bins'].values))
    for b in bins_counter:
        contacts_in_bin = df_centros_bins[df_centros_bins['chr_bins'] == b]
        tmp_df = contacts_in_bin.iloc[:, 3:]
        tmp_mean_df = pd.DataFrame(tmp_df.mean()).T
        tmp_std_df = pd.DataFrame(tmp_df.std()).T
        tmp_mean_df.index = [b]
        tmp_std_df.index = [b]
        df_mean = pd.concat([df_mean, tmp_mean_df])
        df_std = pd.concat([df_std, tmp_std_df])

    df_mean.to_csv(output_file + '_mean_on_cen.csv', sep='\t')
    df_std.to_csv(output_file + '_std_on_cen.csv', sep='\t')
    return df_mean, df_std


def plot_aggregated(mean_df: pd.DataFrame,
                    std_df: pd.DataFrame,
                    info_df: pd.DataFrame,
                    output_path: str):

    x = mean_df.index.tolist()
    for ii, oligo in enumerate(mean_df.columns):
        name = info_df.loc['names', oligo]
        y = mean_df[oligo]
        yerr = std_df[oligo]
        plt.figure(figsize=(14, 12))
        plt.bar(x, y)
        plt.errorbar(x, y, yerr=yerr, fmt="o", color='b', capsize=5)
        plt.title("Aggregated contacts for read {0} from probe {1} around chromosome's centromeres".format(oligo, name))
        plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
        plt.ylabel("Average contacts made and standard deviation")
        plt.savefig(output_path + "{0}-centromeres-aggregated_contacts_plot.{1}".format(oligo, 'jpg'), dpi=99)
        plt.close()


def debug(formatted_contacts_path: str,
          window_size: int,
          centros_coord_path: str,
          output_path: str):

    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    df_contacts_centros, df_infos = contacts_focus_around_centromeres(formatted_contacts_path=formatted_contacts_path,
                                                                      window_size=window_size,
                                                                      centros_infos_path=centros_coord_path)

    output_file = dir_res + output_path.split('/')[-2]
    df_mean, df_std = compute_mean_per_fragment(df_centros_bins=df_contacts_centros, output_file=output_file)
    plot_aggregated(df_mean, df_std, df_infos, dir_res)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    formatted_contacts_path, centros_coordinates_path, window_size, output_path, = ['' for _ in range(4)]

    try:
        opts, args = getopt.getopt(argv, "h:c:m:w:o:", ["--help",
                                                        "--contacts",
                                                        "--coordinates",
                                                        "--window",
                                                        "--output"])
    except getopt.GetoptError:
        print('aggregate centromeres arguments :\n'
              '-c <formatted_contacts_input.csv> (contacts filtered with contacts_filter.py) \n'
              '-m <chr_centros_coordinates.tsv>  \n'
              '-w <window> size at both side of the centromere to look around \n' 
              '-o <output_file_name.csv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromeres arguments :\n'
                  '-c <formatted_contacts_input.csv> (contacts filtered with contacts_filter.py) \n'
                  '-m <chr_centros_coordinates.tsv>  \n'
                  '-w <window> size at both side of the centromere to look around \n'
                  '-o <output_file_name.csv>')
            sys.exit()
        elif opt in ("-c", "--contacts"):
            formatted_contacts_path = arg
        elif opt in ("-m", "--coordinates"):
            centros_coordinates_path = arg
        elif opt in ("-w", "--window"):
            window_size = arg
        elif opt in ("-o", "--output"):
            output_path = arg.split('-filtered_contacts_matrix.csv')[0]

    window_size = int(window_size)
    dir_res = output_path + '/'
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    df_contacts_centros, df_infos = contacts_focus_around_centromeres(formatted_contacts_path=formatted_contacts_path,
                                                                      window_size=window_size,
                                                                      centros_infos_path=centros_coordinates_path)

    output_file = output_path + '/' + output_path.split('/')[-1]
    df_mean, df_std = compute_mean_per_fragment(df_centros_bins=df_contacts_centros, output_file=output_file)
    plot_aggregated(df_mean, df_std, df_infos, dir_res)


if __name__ == "__main__":
    if utils.is_debug():
        centros_coord = "../../../bash_scripts/aggregate_centro/inputs/S288c_chr_centro_coordinates.tsv"

        formatted_contacts = '../../../bash_scripts/aggregate_centro/inputs' \
                             '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
                             '-filtered_contacts_matrix.csv'

        output = "../../../bash_scripts/aggregate_centro/outputs/" \
                 "AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC-filtered_contacts_matrix.csv"

        oligos = "../../../bash_scripts/aggregate_centro/inputs/capture_oligo_positions.csv"
        window = 150000
        debug(formatted_contacts_path=formatted_contacts,
              window_size=window,
              centros_coord_path=centros_coord,
              output_path=output.split('_contacts_matrix.csv')[0]+'/')
    else:
        main()
    print('--- DONE ---')
