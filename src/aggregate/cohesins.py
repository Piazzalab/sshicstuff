#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
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
    excluded_chr = ['chr2', 'chr3']
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
        for c in tmp_df.columns[3:]:
            self_chr = df_info.loc['self_chr', c]
            if self_chr == current_chr:
                tmp_df.loc[:, c] = np.nan

        return tmp_df
    df_res = pd.concat([process_row(row) for _, row in df_peaks.iterrows()])
    df_res = df_res[~df_res['chr'].isin(excluded_chr)]

    df_res.index = range(len(df_res))
    return df_res, df_info


def compute_average_aggregate(
        df_cohesins_peaks_bins: pd.DataFrame,
        df_info: pd.DataFrame,
        table_path: str):
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

    probes = df_info.loc['names', :].values
    df_mean.columns = probes
    df_std.columns = probes

    df_mean.to_csv(table_path + '_mean_on_cohesins.tsv', sep='\t')
    df_std.to_csv(table_path + '_std_on_cohesins.tsv', sep='\t')

    return df_mean, df_std


def plot_aggregated(
        mean_df: pd.DataFrame,
        std_df: pd.DataFrame,
        plot_path: str):

    x = mean_df.index.tolist()
    for ii, probe in enumerate(mean_df.columns):
        if len(probe.split('-/-')) > 1:
            name = '_&_'.join(probe.split('-/-'))
        else:
            name = probe

        y = mean_df[probe]
        yerr = std_df[probe]
        ymin = -np.max((mean_df[probe] + std_df[probe])) * 0.01
        plt.figure(figsize=(18, 12))
        plt.bar(x, y)
        plt.errorbar(x, y, yerr=yerr, fmt="o", color='b', capsize=5)
        plt.ylim((ymin, None))
        plt.title("Aggregated frequencies for probe {0} cohesins peaks".format(name))
        plt.xlabel("Bins around the cohesins peaks (in kb), 5' to 3'")
        plt.xticks(rotation=45)
        plt.ylabel("Average frequency made and standard deviation")
        plt.savefig(plot_path + "_{0}_cohesins_aggregated_frequencies_plot.{1}".format(name, 'jpg'), dpi=99)
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


def run(
        formatted_contacts_path: str,
        window_size: int,
        output_dir: str,
        cohesins_peaks_path: str,
        score_cutoff: int):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(output_path=output_dir+sample_name, score_h=score_cutoff)

    df_contacts_cohesins, df_info = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        cohesins_peaks_path=cohesins_peaks_path,
        score_cutoff=score_cutoff)

    df_mean, df_std = compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_info=df_info,
        table_path=dir_table+sample_name)

    plot_aggregated(
        mean_df=df_mean,
        std_df=df_std,
        plot_path=dir_plot+sample_name)

    print('DONE: ', sample_name)
