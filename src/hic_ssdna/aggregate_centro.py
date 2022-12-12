#! /usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
import pandas as pd
from typing import Optional
from collections import Counter

import sys
import os
import getopt

pd.options.mode.chained_assignment = None

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


def find_nearest(array: list | np.ndarray,
                 key: int | float,
                 mode: str,
                 maxi: Optional[int] = 10000,
                 mini: Optional[int] = 0,) -> int | float:
    array = np.asarray(array)
    if mode == 'upper':
        # the smallest element of array GREATER than key
        if key >= np.max(array):
            return maxi
        else:
            return array[array > key].min()
    elif mode == 'lower':
        # the largest element of array LESS than key
        if key <= np.min(array):
            return mini
        else:
            return array[array < key].max()
    else:
        return array[(np.abs(array - key)).argmin()]


def split_formated_dataframe(df0: pd.DataFrame):
    df0_a = df0.iloc[:5, :]
    df0_b = df0.iloc[5:, :]
    df1 = df0_a[[c for c in df0_a.columns if c not in ['chr', 'chr_bins', 'genome_bins', 'positions']]].astype(str)
    df2 = pd.DataFrame()
    df2['chr'] = df0_b.iloc[:, 0].astype(str)
    df2[df0_b.columns[1:]] = df0_b.iloc[:, 1:].astype(int)

    return df1, df2


def contacts_focus_around_centromeres(formated_contacts_path: str,
                                      window: int,
                                      centro_infos_path: str):

    df_centros = pd.read_csv(centro_infos_path, sep='\t', index_col=None)
    df_all = pd.read_csv(formated_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_infos, df_contacts = split_formated_dataframe(df_all)

    df_res = pd.DataFrame()

    for index, row in df_centros.iterrows():
        current_chr = row[0]
        current_centro_pos = row[2]

        left_cutoff = current_centro_pos - window
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_centro_pos + window
        tmp_df = df_contacts.loc[(df_contacts['chr'] == current_chr) &
                                 (df_contacts['chr_bins'] > left_cutoff) &
                                 (df_contacts['chr_bins'] < right_cutoff)]

        tmp_df.index = range(len(tmp_df))
        current_centro_bin = find_nearest(tmp_df['chr_bins'].values, current_centro_pos, mode='lower')

        for index2, row2 in tmp_df.iterrows():
            tmp_df.iloc[index2, 1] -= current_centro_bin

        for c in tmp_df.columns[3:]:
            self_chr = df_infos.loc['self_chr', c]
            if self_chr == current_chr:
                tmp_df.loc[0:len(tmp_df), c] = np.nan

        df_res = pd.concat([df_res, tmp_df])
    df_res.index = range(len(df_res))
    return df_res, df_infos


def compute_mean_per_fragment(df_centro_bins: pd.DataFrame, output_file: str):
    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    bins_counter = dict(Counter(df_centro_bins['chr_bins'].values))
    for b in bins_counter:
        contacts_in_bin = df_centro_bins[df_centro_bins['chr_bins'] == b]
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

    n = mean_df.shape[1]
    color_list = [colors.to_hex(colors.hsv_to_rgb([x * 1.0 / n, 1.0, 1.0])) for x in range(n)]

    x = mean_df.index.tolist()
    for ii, oligo in enumerate(mean_df.columns):
        name = info_df.loc['names', oligo]
        y = mean_df[oligo]
        yerr = std_df[oligo]
        plt.figure(figsize=(14, 12))
        plt.bar(x, y)
        plt.errorbar(x, y, yerr=yerr, fmt="o", color=color_list[ii], capsize=5)
        plt.title("Aggregated contacts for read {0} from probe {1} around chromosome's centromeres".format(oligo, name))
        plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
        plt.ylabel("Average contacts made and standard deviation")
        plt.savefig(output_path + "{0}-centromeres-aggregated_contacts_plot.{1}".format(oligo, 'jpg'), dpi=160)
        plt.close()


def debug(formated_contacts_path: str,
          window: int,
          centro_coord_path: str,
          output_path: str):

    dir_plots = output_path
    if not os.path.exists(dir_plots):
        os.makedirs(dir_plots)

    df_contacts_centro, df_infos = contacts_focus_around_centromeres(formated_contacts_path=formated_contacts_path,
                                                                     window=window,
                                                                     centro_infos_path=centro_coord_path)

    output_file = output_path + output_path.split('/')[-2]
    df_mean, df_std = compute_mean_per_fragment(df_centro_bins=df_contacts_centro, output_file=output_file)
    plot_aggregated(df_mean, df_std, df_infos, output_path)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    formated_contacts_path, centro_coordinates, window_size, output_path = ['' for _ in range(4)]
    try:
        opts, args = getopt.getopt(argv, "hc:m:w:O:", ["--help",
                                                       "--contacts",
                                                       "--coordinates"
                                                       "--window"
                                                       "--output"])
    except getopt.GetoptError:
        print('aggregate centromere arguments :\n'
              '-c <formated_contacts.csv> (contacts filtered with contacts_format.py) \n'
              '-m <chr_centromeres_coordinates.tsv> (positions of centromeres (CEN) for each chr (1 to 16) \n'
              '-w <window> (number of bp to look around the centromere, at both sides) \n'
              '-O <output_file_name.csv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromere arguments :\n'
                  '-c <formated_contacts.csv> (contacts filtered with contacts_format.py) \n'
                  '-m <chr_centromeres_coordinates.tsv> (positions of centromeres (CEN) for each chr (1 to 16) \n'
                  '-w <window> (number of bp to look around the centromere, at both sides) \n'
                  '-O <output_file_name.csv>')
            sys.exit()
        elif opt in ("-c", "--contacts"):
            formated_contacts_path = arg
        elif opt in ("-m", "--coordinates"):
            centro_coordinates = arg
        elif opt in ("-w", "--window"):
            window_size = int(arg)
        elif opt in ("-O", "--output"):
            output_path = arg.split('contacts_matrix.csv')[0]

    dir_plots = output_path + '/'
    if not os.path.exists(dir_plots):
        os.makedirs(dir_plots)

    df_contacts_centro, df_infos = contacts_focus_around_centromeres(formated_contacts_path=formated_contacts_path,
                                                                     window=window_size,
                                                                     centro_infos_path=centro_coordinates)

    output_file = output_path + output_path.split('/')[-2]
    df_mean, df_std = compute_mean_per_fragment(df_centro_bins=df_contacts_centro, output_file=output_file)
    plot_aggregated(df_mean, df_std, df_infos, output_path)


if __name__ == "__main__":
    if is_debug():
        centro_coord = "../../../bash_scripts/aggregate_centro/inputs/S288c_chr_centro_coordinates.tsv"

        formated_contacts = '../../../bash_scripts/aggregate_centro/inputs' \
                            '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
                            '-filtered_contacts_matrix.csv'

        output = "../../../bash_scripts/aggregate_centro/outputs/" \
                 "AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC-filtered_contacts_matrix.csv"

        oligos = "../../../bash_scripts/aggregate_centro/inputs/capture_oligo_positions.csv"
        window = 150000
        debug(formated_contacts_path=formated_contacts,
              window=window,
              centro_coord_path=centro_coord,
              output_path=output.split('_contacts_matrix.csv')[0]+'/')
    else:
        main()
    print('--- DONE ---')
