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
        centros_info_path: str,
        window_size: int,
        cohesins_peaks_path: str,
        score_cutoff: int,
        filter_range: int,
        filter_mode: str | None):

    df_centros = pd.read_csv(centros_info_path, sep='\t', index_col=None)
    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_peaks = df_peaks[df_peaks['score'] > score_cutoff]
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_info, df_contacts = tools.split_formatted_dataframe(df_all)
    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]
    excluded_chr = ['chr2', 'chr3']

    def process_row(row):
        current_chr = row[0]
        current_peak = row[1]

        if current_chr in excluded_chr:
            return None

        left_cutoff = current_peak - window_size - bin_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_peak + window_size
        tmp_df = df_contacts.query("chr == @current_chr and chr_bins > @left_cutoff and chr_bins < @right_cutoff")
        #   temporary dataframe containing the bins present in the windows for the current chr only
        tmp_df.index = range(len(tmp_df))
        current_cohesin_peak_bin = tools.find_nearest(tmp_df['chr_bins'].values, current_peak, mode='lower')

        if filter_mode is not None:
            filtered_tmp_df = filter_peaks_around_centromeres(
                df_centros=df_centros,
                df_contacts_peaks=tmp_df,
                filter_range=filter_range,
                filter_mode=filter_mode,
                bin_size=bin_size)

            if filtered_tmp_df.shape[0] == 0:
                return None
        else:
            filtered_tmp_df = tmp_df

        #   Indices shifting : bin of centromere becomes 0, bins in downstream becomes negative and bins
        #   in upstream becomes positive.
        filtered_tmp_df.iloc[:, 1] -= current_cohesin_peak_bin
        filtered_tmp_df.index = range(len(filtered_tmp_df))
        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for c in filtered_tmp_df.columns[3:]:
            self_chr = df_info.loc['self_chr', c]
            if self_chr == current_chr:
                filtered_tmp_df.loc[:, c] = np.nan

        return filtered_tmp_df

    df_res = pd.concat([process_row(row) for _, row in df_peaks.iterrows()])
    df_res.index = range(len(df_res))
    return df_res, df_info


def filter_peaks_around_centromeres(
        df_centros: pd.DataFrame,
        df_contacts_peaks: pd.DataFrame,
        filter_range: int,
        filter_mode: str | None,
        bin_size: int):

    current_chr = pd.unique(df_contacts_peaks['chr'])[0]
    current_chr_cen = df_centros.loc[df_centros['Chr'] == current_chr, 'Left_arm_length'].values[0]
    left_cutoff = current_chr_cen - filter_range - bin_size
    right_cutoff = current_chr_cen + filter_range

    if filter_mode == 'outer':
        return df_contacts_peaks.query(
            "chr == @current_chr and (chr_bins < @left_cutoff or chr_bins > @right_cutoff)")
    elif filter_mode == 'inner':
        return df_contacts_peaks.query(
            "chr == @current_chr and chr_bins > @left_cutoff and chr_bins < @right_cutoff")
    elif filter_mode is None:
        return df_contacts_peaks


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


def pooled_stats(mean_df: pd.DataFrame,
                 std_df: pd.DataFrame,
                 table_path: str):

    middle = int(len(mean_df) / 2)
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

    pooled_mean_df.to_csv(table_path + '_pooled_mean_on_cohesins.tsv', sep='\t')
    pooled_mean_df.to_csv(table_path + '_pooled_std_on_cohesins.tsv', sep='\t')

    return pooled_mean_df, pooled_std_df


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
          score_h: int,
          filter_mode: None | str,
          filter_span: int):

    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_type = dir_res + '/cohesins_peaks/'
    if not os.path.exists(dir_type):
        os.makedirs(dir_type)

    if filter_mode is None:
        dir_mode = dir_type + 'all/'
    else:
        dir_mode = dir_type + filter_mode + '_' + str(filter_span // 1000) + 'kb' + '/'
    if not os.path.exists(dir_mode):
        os.makedirs(dir_mode)

    dir_score = dir_mode + str(score_h) + '/'
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
        centromere_info_path: str,
        score_cutoff: int,
        cen_filter_span: int,
        cen_filter_mode: str | None):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(
        output_path=output_dir+sample_name,
        score_h=score_cutoff,
        filter_mode=cen_filter_mode,
        filter_span=cen_filter_span)

    df_contacts_cohesins, df_info = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=formatted_contacts_path,
        centros_info_path=centromere_info_path,
        window_size=window_size,
        cohesins_peaks_path=cohesins_peaks_path,
        score_cutoff=score_cutoff,
        filter_range=cen_filter_span,
        filter_mode=cen_filter_mode)

    df_mean, df_std = compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_info=df_info,
        table_path=dir_table+sample_name)

    df_mean_pooled, df_std_pooled = pooled_stats(
        mean_df=df_mean,
        std_df=df_std,
        table_path=dir_table+sample_name)

    plot_aggregated(
        mean_df=df_mean_pooled,
        std_df=df_std_pooled,
        plot_path=dir_plot+sample_name)

    print('DONE: ', sample_name)
