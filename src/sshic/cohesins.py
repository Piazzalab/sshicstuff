#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from sshic import tools

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


"""
**************************************************************************
*******   AGGREGATED CONTACTS AROUND COHESINS PEAKS FUNCTIONS   **********
**************************************************************************
"""


def freq_focus_around_cohesin_peaks(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        centros_info_path: str,
        window_size: int,
        bin_size: int,
        cohesins_peaks_path: str,
        score_cutoff: int,
        filter_range: int,
        filter_mode: str | None,
        pooled: bool):

    excluded_chr = ['chr2', 'chr3']

    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_peaks = df_peaks[df_peaks['score'] > score_cutoff]
    df_peaks.sort_values(by='score', ascending=False, inplace=True)
    df_peaks.reset_index(drop=True, inplace=True)
    df_peaks = df_peaks.query("score > @score_cutoff and chr not in @excluded_chr")

    df_centros = pd.read_csv(centros_info_path, sep='\t', index_col=None)
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0)
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    df_probes_t = df_probes.transpose()
    unique_fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

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

        if pooled:
            filtered_tmp_df.loc[:, 'chr_bins'] = abs(filtered_tmp_df.loc[:, 'chr_bins'] - current_cohesin_peak_bin)
            filtered_tmp_df = filtered_tmp_df.groupby(['chr', 'chr_bins'], as_index=False).mean()
        else:
            filtered_tmp_df.loc[:, 'chr_bins'] -= current_cohesin_peak_bin

        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for f in unique_fragments:
            self_chr = df_probes_t.loc[df_probes_t['frag_id'] == f]['chr'][0]
            if self_chr == current_chr:
                tmp_df.loc[:, f] = np.nan

        return filtered_tmp_df

    df_res = pd.concat([process_row(row) for _, row in df_peaks.iterrows()])
    df_res.index = range(len(df_res))
    return df_res, df_probes


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
        df_probes: pd.DataFrame,
        table_path: str):

    all_probes = df_probes.columns.values
    unique_chr = pd.unique(df_cohesins_peaks_bins['chr'])
    bins_array = np.unique(df_cohesins_peaks_bins['chr_bins'])

    res: dict = {}
    for probe in all_probes:
        fragment = df_probes.loc['frag_id', probe]
        self_chr = df_probes.loc['chr', probe]
        if fragment not in df_cohesins_peaks_bins.columns:
            continue

        df_freq_cen = df_cohesins_peaks_bins.pivot_table(index='chr_bins', columns='chr', values=fragment, fill_value=np.nan)
        df_freq_cen[self_chr] = np.nan
        df_freq_cen = df_freq_cen[unique_chr].reindex(bins_array)

        res[probe] = df_freq_cen
        df_freq_cen.to_csv(table_path + probe + '_chr1-16_freq_cen.tsv', sep='\t')

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()

    for probe, df in res.items():
        df_mean[probe] = df.T.mean()
        df_std[probe] = df.T.std()

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
        plt.savefig(plot_path + "{0}_cohesins_aggregated_frequencies_plot.{1}".format(name, 'jpg'), dpi=99)
        plt.close()


def mkdir(output_path: str,
          score_h: int,
          filter_mode: None | str,
          filter_span: int):

    dir_res = output_path + '/'
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    if filter_mode is None:
        dir_mode = dir_res + 'all/'
    else:
        dir_mode = dir_res + filter_mode + '_' + str(filter_span // 1000) + 'kb' + '/'
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
        probes_to_fragments_path: str,
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

    df_contacts_cohesins, df_probes = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=formatted_contacts_path,
        probes_to_fragments_path=probes_to_fragments_path,
        centros_info_path=centromere_info_path,
        window_size=window_size,
        bin_size=1000,
        cohesins_peaks_path=cohesins_peaks_path,
        score_cutoff=score_cutoff,
        filter_range=cen_filter_span,
        filter_mode=cen_filter_mode,
        pooled=True)

    df_mean, df_std = compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_probes=df_probes,
        table_path=dir_table)

    plot_aggregated(
        mean_df=df_mean,
        std_df=df_std,
        plot_path=dir_plot)

    print('DONE: ', sample_name)
