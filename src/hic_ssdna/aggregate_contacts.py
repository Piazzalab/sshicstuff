#! /usr/bin/env python3

from typing import Optional
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


"""
**************************************************************************
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def compute_centromere_freq_per_oligo_per_chr(
        df_freq: pd.DataFrame,
        df_info: pd.DataFrame,
        dir_table: str):

    reads_array = df_info.columns.values
    chr_array = np.array(['chr'+str(i) for i in range(1, 17)])
    bins_array = np.unique(df_freq['chr_bins'])

    for ol in reads_array:
        probe = df_info.loc['names', ol]
        if len(probe.split('-/-')) > 1:
            probe = '_&_'.join(probe.split('-/-'))

        df_freq_cen = pd.DataFrame(columns=chr_array, index=bins_array)
        grouped = df_freq.groupby(['chr', 'chr_bins'])
        for name, group in grouped:
            chr_name, bin_name = name
            df_freq_cen.loc[bin_name, chr_name] = group[ol].iloc[0]

        df_freq_cen = df_freq_cen.astype(float)

        #   Write to csv
        df_freq_cen.to_csv(dir_table + probe + '_chr1-16_freq_cen.tsv', sep='\t')


def freq_focus_around_centromeres(formatted_contacts_path: str,
                                  window_size: int,
                                  centros_infos_path: str):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    centromeres and for each of the 16 chromosomes of yeast genome
    """
    #   dataframe containing information about location of the centromere for each chr and the length of chr.
    df_centros = pd.read_csv(centros_infos_path, sep='\t', index_col=None)
    #   dataframe of the formatted contacts csv file previously created,
    #   with DTYPE=object because multiple type are present in columns
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    #   It needs thus to split between numeric and not numeric data
    df_info, df_contacts = utils.split_formatted_dataframe(df_all)

    #   result dataframe with bin around centromeres only
    df_res = pd.DataFrame()

    #   Size of a bin in our formatted file given as input
    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]

    for index, row in df_centros.iterrows():
        current_chr = row[0]
        current_centros_pos = row[2]

        left_cutoff = current_centros_pos - window_size - bin_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_centros_pos + window_size
        tmp_df = df_contacts.loc[(df_contacts['chr'] == current_chr) &
                                 (df_contacts['chr_bins'] > left_cutoff) &
                                 (df_contacts['chr_bins'] < right_cutoff)]

        #   temporary dataframe containing the bins present in the windows for the current chr only
        tmp_df.index = range(len(tmp_df))
        current_centros_bin = utils.find_nearest(tmp_df['chr_bins'].values, current_centros_pos, mode='lower')

        for index2, row2 in tmp_df.iterrows():
            #   Indices shifting : bin of centromere becomes 0, bins in downstream becomes negative and bins
            #   in upstream becomes positive.
            tmp_df.iloc[index2, 1] -= current_centros_bin

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
    #   df_std : same but with standard deviation/error instead of mean
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

    #   Sort the series according to index
    df_mean = df_mean.sort_index()
    df_std = df_std.sort_index()
    df_median = df_median.sort_index()

    #   Concatenate with oligo names, types, locations ...
    df_mean_with_info = pd.concat([df_info, df_mean])
    df_std_with_info = pd.concat([df_info, df_std])
    df_median_with_info = pd.concat([df_info, df_median])

    #   Write to csv
    df_mean_with_info.to_csv(output_file + '_mean_on_cen.tsv', sep='\t')
    df_std_with_info.to_csv(output_file + '_std_on_cen.tsv', sep='\t')
    df_median_with_info.to_csv(output_file + '_median_on_cen.tsv', sep='\t')
    return df_mean, df_std, df_median


"""
**************************************************************************
*******   AGGREGATED CONTACTS AROUND COHESINS PEAKS FUNCTIONS   **********
**************************************************************************
"""


def freq_focus_around_cohsin_peaks(
        formatted_contacts_path: str,
        window_size: int,
        cohesins_peaks_path: str):

    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_info, df_contacts = utils.split_formatted_dataframe(df_all)
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
        current_cohesin_peak_bin = utils.find_nearest(tmp_df['chr_bins'].values, current_cohesin_peak, mode='lower')

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


"""
**************************************************************************
************************   GENERAL FUNCTIONS   ***************************
**************************************************************************
"""


def pooled_stats(mean_df: pd.DataFrame,
                 std_df: pd.DataFrame):

    middle = int(np.where(mean_df.index.values == 0)[0])
    pooled_index = mean_df.index[middle:].values

    #   Pool the mean dataframe
    left_mean_df = mean_df.iloc[:middle+1]
    left_mean_df.index = pooled_index[::-1]
    left_mean_df = left_mean_df.sort_index()
    right_mean_df = mean_df.iloc[middle:]

    tmp_mean_df = pd.concat((left_mean_df, right_mean_df))
    pooled_mean_df = tmp_mean_df.groupby(tmp_mean_df.index).mean()

    #   Pool the std dataframe
    left_std_df = std_df.iloc[:middle + 1]
    left_std_df.index = pooled_index[::-1]
    left_std_df = left_std_df.sort_index()
    right_std_df = std_df.iloc[middle:]
    pooled_std_df = pd.DataFrame()

    for col in left_std_df.columns:
        n1 = left_std_df[col].shape[0]
        n2 = right_std_df[col].shape[0]
        std_pooled = np.sqrt(((n1 - 1) * left_std_df[col] ** 2 + (n2 - 1) * right_std_df[col] ** 2) / (n1 + n2 - 2))
        pooled_std_df[col] = std_pooled

    return pooled_mean_df, pooled_std_df


def plot_aggregated(mean_df: pd.DataFrame,
                    std_df: pd.DataFrame,
                    info_df: pd.DataFrame,
                    mode: str,
                    output_path: str):
    """
    Plot for each oligo/read a barplot of the average number of contacts
    made around the centromeres (average on the 16 chr of yeast).
    Gives also the standard deviation.
    """

    pooled_mean_df, pooled_std_df = pooled_stats(mean_df=mean_df, std_df=std_df)
    df_min = pooled_mean_df - pooled_std_df
    df_max = pooled_mean_df + pooled_std_df
    y_min = np.min(df_min.values) + 0.05 * np.min(df_min.values)
    y_max = np.max(df_max.values) + 0.05 * np.max(df_max.values)

    x = pooled_mean_df.index.tolist()
    for ii, oligo in enumerate(pooled_mean_df.columns):
        probe = info_df.loc['names', oligo]
        if len(probe.split('-/-')) > 1:
            probe = '_&_'.join(probe.split('-/-'))

        y = pooled_mean_df[oligo]
        yerr = pooled_std_df[oligo]
        plt.figure(figsize=(18, 12))
        plt.bar(x, y)
        plt.errorbar(x, y, yerr=yerr, fmt="o", color='b', capsize=5)
        plt.ylim((y_min, y_max))
        plt.title(
            "Aggregated frequencies for read {0} from probe {1} around {2}".format(oligo, probe, mode))
        plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
        plt.ylabel("Average frequency made and standard deviation")
        if mode == 'centromeres':
            plt.savefig(output_path + "{0}-centromeres-aggregated_frequencies_plot.{1}".format(probe, 'jpg'), dpi=99)
        elif mode == 'cohesins':
            plt.savefig(output_path + "{0}-cohesins_peaks-aggregated_frequencies_plot.{1}".format(probe, 'jpg'), dpi=99)
        plt.close()


def mkdir(output_path: str,
          mode: str):

    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_res2 = ''
    if mode == 'centromeres':
        dir_res2 = dir_res + '/centromeres/'
    elif mode == 'cohesins':
        dir_res2 = dir_res + '/cohesins_peaks/'
    if not os.path.exists(dir_res2):
        os.makedirs(dir_res2)

    dir_plot = dir_res2 + 'plots/'
    if not os.path.exists(dir_plot):
        os.makedirs(dir_plot)

    dir_table = dir_res2 + 'tables/'
    if not os.path.exists(dir_table):
        os.makedirs(dir_table)

    return dir_table, dir_plot


def debug(formatted_contacts_path: str,
          window_size: int,
          output_path: str,
          centros_coord_path: Optional[str] = None,
          cohesins_peaks_path: Optional[str] = None):

    if centros_coord_path is not None:

        dir_table, dir_plot = mkdir(output_path=output_path, mode='centromeres')
        output_file = dir_table + output_path.split('/')[-2]

        df_contacts_centros, df_info = freq_focus_around_centromeres(
            formatted_contacts_path=formatted_contacts_path,
            window_size=window_size,
            centros_infos_path=centros_coord_path)

        compute_centromere_freq_per_oligo_per_chr(
            df_freq=df_contacts_centros, df_info=df_info, dir_table=dir_table)

        df_mean, df_std, df_median = compute_aggregate_around_centromeres(
            df_centros_bins=df_contacts_centros,
            df_info=df_info,
            output_file=output_file)

        plot_aggregated(df_mean, df_std, df_info, 'centromeres', dir_plot)

    elif cohesins_peaks_path is not None:

        dir_table, dir_plot = mkdir(output_path=output_path, mode='cohesins')
        output_file = dir_table + output_path.split('/')[-2]

        df_contacts_cohesins, df_info = freq_focus_around_cohsin_peaks(
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

    binned_contacts_path, coordinates_path, window_size, output_path, mode = [None for _ in range(5)]

    try:
        opts, args = getopt.getopt(argv, "h:b:c:m:w:o:", ["--help",
                                                          "--binning",
                                                          "--coordinates",
                                                          "--mode",
                                                          "--window",
                                                          "--output"])
    except getopt.GetoptError:
        print('aggregate centromeres arguments :\n'
              '-b <binned_frequencies_matrix.csv> (contacts filtered with contacts_filter.py) \n'
              '-c <chr_centros_coordinates.tsv> or  <chr_cohesins_peaks_coordinates.bed> \n'
              '-m <mode> if we want to aggregate on chr centromeres or on cohesins peak positions'
              '-w <window> size at both side of the centromere to look around \n'
              '-o <output_file_name.tsv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromeres arguments :\n'
                  '-b <binned_frequencies_matrix.csv> (contacts filtered with contacts_filter.py) \n'
                  '-c <chr_centros_coordinates.tsv> or  <chr_cohesins_peaks_coordinates.bed> \n'
                  '-m <mode> if we want to aggregate on chr centromeres or on cohesins peak positions'
                  '-w <window> size at both side of the centromere to look around \n'
                  '-o <output_file_name.tsv>')
            sys.exit()
        elif opt in ("-b", "--binning"):
            binned_contacts_path = arg
        elif opt in ("-c", "--coordinates"):
            coordinates_path = arg
        elif opt in ("-m", "--mode"):
            mode = arg
        elif opt in ("-w", "--window"):
            window_size = int(arg)
        elif opt in ("-o", "--output"):
            output_path = arg
            if 'formatted' in output_path:
                output_path = output_path.split('formatted_frequencies_matrix.tsv')[0]
            else:
                output_path = output_path.split('_frequencies_matrix.tsv')[0]

    if mode == 'centromeres':

        dir_table, dir_plot = mkdir(output_path=output_path, mode='centromeres')
        output_file = dir_table + '/' + output_path.split('/')[-1]

        df_contacts_centros, df_info = freq_focus_around_centromeres(
            formatted_contacts_path=binned_contacts_path,
            window_size=window_size,
            centros_infos_path=coordinates_path)

        compute_centromere_freq_per_oligo_per_chr(
            df_freq=df_contacts_centros, df_info=df_info, dir_table=dir_table)

        df_mean, df_std, df_median = compute_aggregate_around_centromeres(
            df_centros_bins=df_contacts_centros,
            df_info=df_info,
            output_file=output_file)

        plot_aggregated(df_mean, df_std, df_info, 'centromeres', dir_plot)

    elif mode == 'cohesins':

        dir_table, dir_plot = mkdir(output_path=output_path, mode='cohesins')
        output_file = dir_table + '/' + output_path.split('/')[-1]

        df_contacts_cohesins, df_info = freq_focus_around_cohsin_peaks(
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
    if utils.is_debug():
        #   Debug is mainly used for testing function of the script
        #   Parameters have to be declared here
        cohesins_peaks = "../../../bash_scripts/aggregate_contacts/inputs/HB65_reference_peaks_score50min.bed"
        centros_coord = "../../../bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"

        formatted_contacts_10kb = \
            '../../../bash_scripts/aggregate_contacts/inputs' \
            '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
            '_10kb_frequencies_matrix.tsv'

        formatted_contacts_1kb = \
            '../../../bash_scripts/aggregate_contacts/inputs' \
            '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
            '_1kb_frequencies_matrix.tsv'

        output = "../../../bash_scripts/aggregate_contacts/outputs/" \
                 "AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC"

        oligos = "../../../bash_scripts/aggregate_contacts/inputs/capture_oligo_positions.tsv"
        window = 15000

        debug(formatted_contacts_path=formatted_contacts_1kb,
              window_size=window,
              output_path=output.split('_frequencies_matrix.tsv')[0] + '/',
              cohesins_peaks_path=cohesins_peaks)

    else:
        main()

    print('--- DONE ---')

