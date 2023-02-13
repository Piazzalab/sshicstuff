#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from typing import Optional
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
        cohesins_peaks_path: str,
        score_cutoff: int,
        filter_range: int,
        filter_mode: str | None):

    excluded_chr = ['chr2', 'chr3', '2_micron', 'mitochondrion', 'chr_artificial']
    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_peaks = df_peaks[df_peaks['score'] > score_cutoff]
    df_peaks = df_peaks[~df_peaks['chr'].isin(excluded_chr)]
    df_peaks.sort_values(by='score', ascending=False, inplace=True)
    df_peaks.reset_index(drop=True, inplace=True)
    df_peaks = df_peaks.query("score > @score_cutoff and chr not in @excluded_chr")

    df_centros = pd.read_csv(centros_info_path, sep='\t', index_col=None)
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
    bin_size = df_contacts.loc[1, 'chr_bins'] - df_contacts.loc[0, 'chr_bins']
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    unique_fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    df_merged = pd.merge(df_contacts, df_peaks, on='chr')
    df_merged_cohesins_areas = df_merged[
        (df_merged.chr_bins > (df_merged.start-window_size)) &
        (df_merged.chr_bins < (df_merged.end+window_size))
    ]

    df_merged2 = pd.merge(df_merged_cohesins_areas, df_centros, on='chr')
    if filter_mode == 'inner':
        df_merged_cohesins_areas_filtered = df_merged2[
            (df_merged2.chr_bins > (df_merged2.left_arm_length-filter_range)) &
            (df_merged2.chr_bins < (df_merged2.left_arm_length+filter_range))
        ]
    elif filter_mode == 'outer':
        df_merged_cohesins_areas_filtered = df_merged2[
            (df_merged2.chr_bins < (df_merged2.left_arm_length - filter_range)) |
            (df_merged2.chr_bins > (df_merged2.left_arm_length + filter_range))
        ]
    else:
        df_merged_cohesins_areas_filtered = df_merged2.copy(deep=True)

    df_merged_cohesins_areas_filtered['chr_bins'] =\
        abs(df_merged_cohesins_areas_filtered['chr_bins'] -
            (df_merged_cohesins_areas_filtered['start'] // bin_size) * bin_size)

    df_merged_cohesins_areas_filtered.drop(columns=['start', 'end', 'score', 'length', 'uid',
                                                    'left_arm_length', 'right_arm_length'],
                                           axis=1, inplace=True)

    df_res = df_merged_cohesins_areas_filtered.groupby(['chr', 'chr_bins'], as_index=False).mean()
    df_res = tools.sort_by_chr(df_res, 'chr', 'chr_bins')

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in unique_fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_res.loc[df_res['chr'] == probe_chr, int(f)] = np.nan
        if df_res[f].sum() > 0:
            df_res[f] /= df_res[f].sum()

    return df_res, df_probes


def compute_average_aggregate(
        df_cohesins_peaks_bins: pd.DataFrame,
        df_probes: pd.DataFrame,
        table_path: str,
        plot: bool,
        plot_path: Optional[str]
):

    all_probes = df_probes.index.values
    unique_chr = pd.unique(df_cohesins_peaks_bins['chr'])
    bins_array = np.unique(df_cohesins_peaks_bins['chr_bins'])

    res: dict = {}
    for probe in all_probes:
        fragment = str(df_probes.loc[probe, 'frag_id'])
        self_chr = df_probes.loc[probe, 'chr']
        if fragment not in df_cohesins_peaks_bins.columns:
            continue

        df_freq_cen =\
            df_cohesins_peaks_bins.pivot_table(index='chr_bins', columns='chr', values=fragment, fill_value=np.nan)
        df_freq_cen[self_chr] = np.nan
        df_freq_cen = df_freq_cen[unique_chr].reindex(bins_array)

        res[probe] = df_freq_cen
        df_freq_cen.to_csv(table_path + probe + '_chr1-16_freq_cen.tsv', sep='\t')

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()

    for probe, df in res.items():
        mean = df.T.mean()
        std = df.T.std()

        if plot:
            pos = mean.index
            ymin = -np.max((mean + std)) * 0.01
            plt.figure(figsize=(16, 12))
            plt.bar(pos, mean)
            plt.errorbar(pos, mean, yerr=std, fmt="o", color='g', capsize=5, clip_on=True)
            plt.ylim((ymin, None))
            plt.title("Aggregated frequencies for probe {0} cohesins peaks".format(probe))
            plt.xlabel("Bins around the cohesins peaks (in kb), 5' to 3'")
            plt.xticks(rotation=45)
            plt.ylabel("Average frequency made and standard deviation")
            plt.savefig(plot_path + "{0}_cohesins_aggregated_frequencies_plot.{1}".format(probe, 'jpg'), dpi=96)
            plt.close()

        df_mean[probe] = mean
        df_std[probe] = std
    df_mean.to_csv(table_path + '_mean_on_cohesins.tsv', sep='\t')
    df_std.to_csv(table_path + '_std_on_cohesins.tsv', sep='\t')


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
        cohesins_peaks_path: str,
        centromere_info_path: str,
        score_cutoff: int,
        cen_filter_span: int,
        cen_filter_mode: str | None,
        output_dir: str,
        plot: bool = True
):

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
        cohesins_peaks_path=cohesins_peaks_path,
        score_cutoff=score_cutoff,
        filter_range=cen_filter_span,
        filter_mode=cen_filter_mode)

    compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_probes=df_probes,
        table_path=dir_table,
        plot=plot,
        plot_path=dir_plot)
