#! /usr/bin/env python3

import numpy as np
import pandas as pd
from collections import Counter
import sys
import getopt

from common import plot_aggregated, mkdir
from utils import tools

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


"""
**************************************************************************
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def compute_telomere_freq_per_oligo_per_chr(
        df_freq: pd.DataFrame,
        df_info: pd.DataFrame,
        dir_table: str):

    reads_array = df_info.columns.values
    chr_array = np.array(['chr'+str(i) for i in range(1, 17)])
    bin_size = df_freq.iloc[1, 2] - df_freq.iloc[0, 2]
    bins_array = pd.unique(df_freq['chr_bins'])

    for ol in reads_array:
        probe = df_info.loc['names', ol]
        if len(probe.split('-/-')) > 1:
            probe = '_&_'.join(probe.split('-/-'))

        df_freq_telo = pd.DataFrame(columns=chr_array, index=bins_array)
        grouped = df_freq.groupby(['chr', 'chr_bins'])
        for name, group in grouped:
            chr_name, bin_name = name
            df_freq_telo.loc[bin_name, chr_name] = group[ol].iloc[0]

        df_freq_telo = df_freq_telo.astype(float)

        #   Write to csv
        df_freq_telo.to_csv(dir_table + probe + '_chr1-16_freq_telo.tsv', sep='\t')


def freq_focus_around_centromeres(formatted_contacts_path: str,
                                  window_size: int,
                                  telomeres_coord_path: str):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    centromeres and for each of the 16 chromosomes of yeast genome
    """
    #   dataframe containing information about location of the centromere for each chr and the length of chr.
    df_centro = pd.read_csv(telomeres_coord_path, sep='\t', index_col=None)

    #   dataframe containing position of telemore for each chr
    #   i.e., pos 0 for left telomere and position length of chr for right telomere
    df_telo = pd.DataFrame({'chr': df_centro['Chr'], 'telo_l': 0, 'telo_r': df_centro['Length']})
    #   dataframe of the formatted contacts csv file previously created,
    #   with DTYPE=object because multiple type are present in columns
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    #   It needs thus to split between numeric and not numeric data
    df_info, df_contacts = tools.split_formatted_dataframe(df_all)

    #   result dataframe with bin around centromeres only
    df_res = pd.DataFrame()

    #   Size of a bin in our formatted file given as input
    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]

    for index, row in df_telo.iterrows():
        current_chr = row[0]
        current_telo_left = row[1]
        current_telo_right = row[2]

        tmp_df_left = df_contacts.loc[(df_contacts['chr'] == current_chr) &
                                      (df_contacts['chr_bins'] >= current_telo_left) &
                                      (df_contacts['chr_bins'] < window_size)]

        tmp_df_right = df_contacts.loc[(df_contacts['chr'] == current_chr) &
                                       (df_contacts['chr_bins'] > current_telo_right - window_size) &
                                       (df_contacts['chr_bins'] < current_telo_right)]

        right_telo_bin = tools.find_nearest(tmp_df_right['chr_bins'].values, current_telo_right, mode='lower')
        for index_r, row_r in tmp_df_right.iterrows():
            #   Indices shifting : bin of telomere becomes 0, bins in downstream becomes negative and bins
            #   in upstream becomes positive.
            tmp_df_right.loc[index_r, 'chr_bins'] = \
                "3'_" + str(abs(tmp_df_right.loc[index_r, 'chr_bins'] - right_telo_bin))

        for index_l, row_l in tmp_df_left.iterrows():
            #   Indices shifting : bin of telomere becomes 0, bins in downstream becomes negative and bins
            #   in upstream becomes positive.
            tmp_df_left.loc[index_l, 'chr_bins'] = \
                "5'_" + str(tmp_df_left.loc[index_l, 'chr_bins'])

        tmp_df = pd.concat((tmp_df_left, tmp_df_right))
        tmp_df.index = range(len(tmp_df))

        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for c in tmp_df.columns[3:]:
            self_chr = df_info.loc['self_chr', c]
            if self_chr == current_chr:
                tmp_df.loc[0:len(tmp_df), c] = np.nan

        #   Concatenate the temporary dataframe of the current chr with
        #   the results dataframe containing other chromosomes
        df_res = pd.concat([df_res, tmp_df])
    df_res.index = range(len(df_res))
    return df_res, df_info


def compute_aggregate_around_centromeres(
        df_centros_bins: pd.DataFrame,
        df_info: pd.DataFrame,
        output_file: str):
    """
    After fetching the contacts for each oligos around the centromere of the 16 chr,
    we need to make an average (and std) of the 16 chr.
    """

    #  df_mean :  dataframe with average contacts in the centromere areas (according to the window the user gives)
    #       for each oligo.
    #  df_std : same but with standard deviation/error instead of mean
    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    df_median = pd.DataFrame()
    bins_counter = dict(Counter(df_centros_bins['chr_bins'].values))
    for b in bins_counter:
        contacts_in_bin = df_centros_bins[df_centros_bins['chr_bins'] == b]
        tmp_df = contacts_in_bin.iloc[:, 3:]
        tmp_mean_df = pd.DataFrame(tmp_df.mean()).T
        tmp_std_df = pd.DataFrame(tmp_df.std()).T
        tmp_median_df = pd.DataFrame(tmp_df.median()).T
        tmp_mean_df.index = [b]
        tmp_std_df.index = [b]
        tmp_median_df.index = [b]
        df_mean = pd.concat([df_mean, tmp_mean_df])
        df_std = pd.concat([df_std, tmp_std_df])
        df_median = pd.concat([df_median, tmp_median_df])

    #   Concatenate with oligo names, types, locations ...
    df_mean_with_info = pd.concat([df_info, df_mean])
    df_std_with_info = pd.concat([df_info, df_std])
    df_median_with_info = pd.concat([df_info, df_median])

    #   Write to csv
    df_mean_with_info.to_csv(output_file + '_mean_on_cen.tsv', sep='\t')
    df_std_with_info.to_csv(output_file + '_std_on_cen.tsv', sep='\t')
    df_median_with_info.to_csv(output_file + '_median_on_cen.tsv', sep='\t')
    return df_mean, df_std, df_median


def debug(formatted_contacts_path: str,
          window_size: int,
          output_path: str,
          telomeres_coord_path: str):

    dir_table, dir_plot = mkdir(output_path=output_path, mode='telomeres')
    output_file = dir_table + output_path.split('/')[-2]

    df_contacts_centros, df_info = freq_focus_around_centromeres(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        telomeres_coord_path=telomeres_coord_path)

    compute_telomere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros, df_info=df_info, dir_table=dir_table)

    df_mean, df_std, df_median = compute_aggregate_around_centromeres(
        df_centros_bins=df_contacts_centros,
        df_info=df_info,
        output_file=output_file)

    plot_aggregated(
        mean_df=df_mean,
        std_df=df_std,
        info_df=df_info,
        mode='telomeres',
        output_path=dir_plot,
        pooled=False)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    binned_contacts_path, coordinates_path, window_size, output_path = [None for _ in range(4)]

    try:
        opts, args = getopt.getopt(argv, "h:b:c:w:o:", ["--help",
                                                        "--binning",
                                                        "--coordinates",
                                                        "--window",
                                                        "--output"])
    except getopt.GetoptError:
        print('aggregate centromeres arguments :\n'
              '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
              '-c <chr_centros_coordinates.tsv> \n'
              '-w <window> size at both side of the centromere to look around \n'
              '-o <output_file_name.tsv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromeres arguments :\n'
                  '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
                  '-c <chr_centros_coordinates.tsv> \n'
                  '-w <window> size at both side of the centromere to look around \n'
                  '-o <output_file_name.tsv>')
            sys.exit()
        elif opt in ("-b", "--binning"):
            binned_contacts_path = arg
        elif opt in ("-c", "--coordinates"):
            coordinates_path = arg
        elif opt in ("-w", "--window"):
            window_size = int(arg)
        elif opt in ("-o", "--output"):
            output_path = arg
            if 'formatted' in output_path:
                output_path = output_path.split('formatted_frequencies_matrix.tsv')[0]
            else:
                output_path = output_path.split('_frequencies_matrix.tsv')[0]

    dir_table, dir_plot = mkdir(output_path=output_path, mode='telomeres')
    output_file = dir_table + '/' + output_path.split('/')[-1]

    df_contacts_centros, df_info = freq_focus_around_centromeres(
        formatted_contacts_path=binned_contacts_path,
        window_size=window_size,
        telomeres_coord_path=coordinates_path)

    compute_telomere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros, df_info=df_info, dir_table=dir_table)

    df_mean, df_std, df_median = compute_aggregate_around_centromeres(
        df_centros_bins=df_contacts_centros,
        df_info=df_info,
        output_file=output_file)

    plot_aggregated(df_mean, df_std, df_info, 'telomeres', dir_plot)


if __name__ == "__main__":
    #   Go into debug function if debug mode is detected, else go for main script with sys arguments
    if tools.is_debug():
        #   Debug is mainly used for testing function of the script
        #   Parameters have to be declared here
        telo_coord = "../../../../bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"

        formatted_contacts_10kb = \
            '../../../../bash_scripts/aggregate_contacts/inputs' \
            '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
            '_10kb_frequencies_matrix.tsv'

        output = "../../../../bash_scripts/aggregate_contacts/outputs/" \
                 "AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC"

        oligos = "../../../../bash_scripts/aggregate_contacts/inputs/capture_oligo_positions.tsv"
        window = 100000

        debug(formatted_contacts_path=formatted_contacts_10kb,
              window_size=window,
              output_path=output.split('_frequencies_matrix.tsv')[0] + '/',
              telomeres_coord_path=telo_coord)

    else:
        main()

    print('--- DONE ---')
