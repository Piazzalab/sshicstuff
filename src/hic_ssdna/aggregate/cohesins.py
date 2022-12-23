#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import getopt

from utils import tools


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
        cohesins_peaks_path: str,
        score_cutoff: int):

    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_peaks = df_peaks[df_peaks['score'] > score_cutoff]
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_info, df_contacts = tools.split_formatted_dataframe(df_all)
    df_res = pd.DataFrame()
    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]

    def process_row(row):
        current_chr = row[0]
        current_peak = row[1]

        left_cutoff = current_peak - window_size - bin_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_peak + window_size
        tmp_df = df_contacts.query("chr == @current_chr and chr_bins > @left_cutoff and chr_bins < @right_cutoff")
        #   temporary dataframe containing the bins present in the windows for the current chr only
        tmp_df.index = range(len(tmp_df))
        current_cohesin_peak_bin = tools.find_nearest(tmp_df['chr_bins'].values, current_peak, mode='lower')

        #   Indices shifting : bin of centromere becomes 0, bins in downstream becomes negative and bins
        #   in upstream becomes positive.
        tmp_df.iloc[:, 1] -= current_cohesin_peak_bin

        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        df_res.loc[:, df_res.columns[3:]] = \
            df_res.loc[:, df_res.columns[3:]].apply(
                lambda x: x.map(lambda y: np.nan if df_info.loc['self_chr', x.name].isin([current_chr]) else y)
            )

        return tmp_df
    df_res = pd.concat([process_row(row) for _, row in df_peaks.iterrows()])
    df_res.index = range(len(df_res))
    return df_res, df_info


def compute_average_aggregate(
        df_cohesins_peaks_bins: pd.DataFrame,
        df_info: pd.DataFrame,
        output_file: str):
    def process_bin(group):
        sub_group = group.iloc[:, 3:]
        tmp_mean_df = pd.DataFrame(sub_group.mean()).T
        tmp_std_df = pd.DataFrame(sub_group.std()).T
        tmp_mean_df.index = [group.name]
        tmp_std_df.index = [group.name]
        return tmp_mean_df, tmp_std_df

    df_cohesins_peaks_bins_grouped = df_cohesins_peaks_bins.groupby('chr_bins')
    df_mean, df_std = zip(*df_cohesins_peaks_bins_grouped.apply(process_bin))
    df_mean = pd.concat(df_mean)
    df_std = pd.concat(df_std)

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


def plot_aggregated(
        mean_df: pd.DataFrame,
        std_df: pd.DataFrame,
        info_df: pd.DataFrame,
        output_path: str):

    x = mean_df.index.tolist()
    for ii, oligo in enumerate(mean_df.columns):
        probe = info_df.loc['names', oligo]
        if len(probe.split('-/-')) > 1:
            probe = '_&_'.join(probe.split('-/-'))

        y = mean_df[oligo]
        yerr = std_df[oligo]
        ymin = -np.max((mean_df[oligo] + std_df[oligo])) * 0.01
        plt.figure(figsize=(18, 12))
        plt.bar(x, y)
        plt.errorbar(x, y, yerr=yerr, fmt="o", color='b', capsize=5)
        plt.ylim((ymin, None))
        plt.title("Aggregated frequencies for probe {0} cohesins peaks".format(probe))
        plt.xlabel("Bins around the cohesins peaks (in kb), 5' to 3'")
        plt.xticks(rotation=45)
        plt.ylabel("Average frequency made and standard deviation")
        plt.savefig(output_path + "{0}-cohesins-aggregated_frequencies_plot.{1}".format(probe, 'jpg'), dpi=99)
        plt.close()


def mkdir(output_path: str,
          score_h: int):
    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_type = dir_res + '/cohesins_peaks/'
    if not os.path.exists(dir_type):
        os.makedirs(dir_type)

    dir_score = dir_type + '/' + str(score_h) + '/'
    if not os.path.exists(dir_score):
        os.makedirs(dir_score)

    dir_plot = dir_score + 'plots/'
    if not os.path.exists(dir_plot):
        os.makedirs(dir_plot)

    dir_table = dir_score + 'tables/'
    if not os.path.exists(dir_table):
        os.makedirs(dir_table)

    return dir_table, dir_plot


def debug(formatted_contacts_path: str,
          window_size: int,
          output_path: str,
          cohesins_peaks_path: str,
          score_cutoff: int):

    dir_table, dir_plot = mkdir(output_path=output_path, score_h=score_cutoff)
    output_file = dir_table + output_path.split('/')[-2]

    df_contacts_cohesins, df_info = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        cohesins_peaks_path=cohesins_peaks_path,
        score_cutoff=score_cutoff)

    df_mean, df_std = compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_info=df_info,
        output_file=output_file)

    plot_aggregated(
        mean_df=df_mean,
        std_df=df_std,
        info_df=df_info,
        output_path=dir_plot)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    binned_contacts_path, coordinates_path, window_size, output_path, score_cutoff = [None for _ in range(5)]

    try:
        opts, args = getopt.getopt(argv, "h:b:c:w:s:o:", ["--help",
                                                          "--binning",
                                                          "--coordinates",
                                                          "--window",
                                                          "--score",
                                                          "--output"])
    except getopt.GetoptError:
        print('aggregate centromeres arguments :\n'
              '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
              '-c <chr_cohesins_peaks_coordinates.bed> \n'
              '-w <window> size at both side of the centromere to look around \n'
              '-s <score> select peak that have a score higher than s \n'
              '-o <output_file_name.tsv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromeres arguments :\n'
                  '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
                  '-c <chr_cohesins_peaks_coordinates.bed> \n'
                  '-w <window> size at both side of the centromere to look around \n'
                  '-s <score> select peak that have a score higher than s \n'
                  '-o <output_file_name.tsv>')
            sys.exit()
        elif opt in ("-b", "--binning"):
            binned_contacts_path = arg
        elif opt in ("-c", "--coordinates"):
            coordinates_path = arg
        elif opt in ("-w", "--window"):
            window_size = int(arg)
        elif opt in ("-s", "--score"):
            score_cutoff = int(arg)
        elif opt in ("-o", "--output"):
            output_path = arg
            if 'formatted' in output_path:
                output_path = output_path.split('formatted_frequencies_matrix.tsv')[0]
            else:
                output_path = output_path.split('_frequencies_matrix.tsv')[0]

    dir_table, dir_plot = mkdir(output_path=output_path, score_h=score_cutoff)
    output_file = dir_table + '/' + output_path.split('/')[-1]

    df_contacts_cohesins, df_info = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=binned_contacts_path,
        window_size=window_size,
        cohesins_peaks_path=coordinates_path,
        score_cutoff=score_cutoff)

    df_mean, df_std = compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_info=df_info,
        output_file=output_file)

    plot_aggregated(
        mean_df=df_mean,
        std_df=df_std,
        info_df=df_info,
        output_path=dir_plot)


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
        score = 1000

        debug(formatted_contacts_path=formatted_contacts_1kb,
              window_size=window,
              output_path=output.split('_frequencies_matrix.tsv')[0] + '/',
              cohesins_peaks_path=cohesins_peaks,
              score_cutoff=score)

    else:
        main()

    print('--- DONE ---')
