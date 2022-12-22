#! /usr/bin/env python3

import numpy as np
import pandas as pd
from collections import Counter
import sys
import getopt

from utils import tools
from common import plot_aggregated, mkdir


#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


"""
**************************************************************************
*******   AGGREGATED CONTACTS AROUND COHESINS PEAKS FUNCTIONS   **********
**************************************************************************
"""


def freq_focus_around_cohesin_peaks(
        formatted_contacts_path: str,
        window_size: int,
        cohesins_peaks_path: str):

    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_info, df_contacts = tools.split_formatted_dataframe(df_all)
    df_res = pd.DataFrame()
    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]

    for index, row in df_peaks.iterrows():
        current_cohesin_chr = row[0]
        current_cohesin_peak = row[1]
        current_cohesin_score = row[4]

        left_cutoff = current_cohesin_peak - window_size - bin_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_cohesin_peak + window_size
        tmp_df = df_contacts.loc[(df_contacts['chr'] == current_cohesin_chr) &
                                 (df_contacts['chr_bins'] > left_cutoff) &
                                 (df_contacts['chr_bins'] < right_cutoff)]

        #   temporary dataframe containing the bins present in the windows for the current chr only
        tmp_df.index = range(len(tmp_df))
        current_cohesin_peak_bin = tools.find_nearest(tmp_df['chr_bins'].values, current_cohesin_peak, mode='lower')

        for index2, row2 in tmp_df.iterrows():
            #   Indices shifting : bin of centromere becomes 0, bins in downstream becomes negative and bins
            #   in upstream becomes positive.
            tmp_df.iloc[index2, 1] -= current_cohesin_peak_bin

        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for c in tmp_df.columns[3:]:
            self_chr = df_info.loc['self_chr', c]
            if self_chr == current_cohesin_chr:
                tmp_df.loc[0:len(tmp_df), c] = np.nan

        #   Concatenate the temporary dataframe of the current chr with
        #   the results dataframe containing other chromosomes
        df_res = pd.concat([df_res, tmp_df])
    df_res.index = range(len(df_res))
    return df_res, df_info


def compute_aggregate_around_cohesins_peaks(
        df_cohesins_peaks_bins: pd.DataFrame,
        df_info: pd.DataFrame,
        output_file: str):

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    bins_counter = dict(Counter(df_cohesins_peaks_bins['chr_bins'].values))
    for b in bins_counter:
        contacts_in_bin = df_cohesins_peaks_bins[df_cohesins_peaks_bins['chr_bins'] == b]
        tmp_df = contacts_in_bin.iloc[:, 3:]
        tmp_mean_df = pd.DataFrame(tmp_df.mean()).T
        tmp_std_df = pd.DataFrame(tmp_df.std()).T
        tmp_mean_df.index = [b]
        tmp_std_df.index = [b]
        df_mean = pd.concat([df_mean, tmp_mean_df])
        df_std = pd.concat([df_std, tmp_std_df])

    #   Sort the series according to index
    df_mean = df_mean.sort_index()
    df_std = df_std.sort_index()

    #   Concatenate with oligo names, types, locations ...
    df_mean_with_info = pd.concat([df_info, df_mean])
    df_std_with_info = pd.concat([df_info, df_std])

    #   Write to csv
    df_mean_with_info.to_csv(output_file + '_mean_on_cohesins_peaks.tsv', sep='\t')
    df_std_with_info.to_csv(output_file + '_std_on_cohesins_peaks.tsv', sep='\t')
    return df_mean, df_std


def debug(formatted_contacts_path: str,
          window_size: int,
          output_path: str,
          cohesins_peaks_path: str):

    dir_table, dir_plot = mkdir(output_path=output_path, mode='cohesins')
    output_file = dir_table + output_path.split('/')[-2]

    df_contacts_cohesins, df_info = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        cohesins_peaks_path=cohesins_peaks_path)

    df_mean, df_std = compute_aggregate_around_cohesins_peaks(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_info=df_info,
        output_file=output_file)

    plot_aggregated(df_mean, df_std, df_info, 'cohesins', dir_plot)


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
              '-c <chr_cohesins_peaks_coordinates.bed> \n'
              '-w <window> size at both side of the centromere to look around \n'
              '-o <output_file_name.tsv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromeres arguments :\n'
                  '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
                  '-c <chr_cohesins_peaks_coordinates.bed> \n'
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

    dir_table, dir_plot = mkdir(output_path=output_path, mode='cohesins')
    output_file = dir_table + '/' + output_path.split('/')[-1]

    df_contacts_cohesins, df_info = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=binned_contacts_path,
        window_size=window_size,
        cohesins_peaks_path=coordinates_path)

    df_mean, df_std = compute_aggregate_around_cohesins_peaks(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_info=df_info,
        output_file=output_file)

    plot_aggregated(df_mean, df_std, df_info, 'cohesins', dir_plot)


if __name__ == "__main__":
    #   Go into debug function if debug mode is detected, else go for main script with sys arguments
    if tools.is_debug():
        #   Debug is mainly used for testing function of the script
        #   Parameters have to be declared here
        cohesins_peaks = "../../../../bash_scripts/aggregate_contacts/inputs/HB65_reference_peaks_score50min.bed"

        formatted_contacts_1kb = \
            '../../../../bash_scripts/aggregate_contacts/inputs' \
            '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
            '_1kb_frequencies_matrix.tsv'

        output = "../../../../bash_scripts/aggregate_contacts/outputs/" \
                 "AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC"

        oligos = "../../../../bash_scripts/aggregate_contacts/inputs/capture_oligo_positions.tsv"
        window = 15000

        debug(formatted_contacts_path=formatted_contacts_1kb,
              window_size=window,
              output_path=output.split('_frequencies_matrix.tsv')[0] + '/',
              cohesins_peaks_path=cohesins_peaks)

    else:
        main()

    print('--- DONE ---')
