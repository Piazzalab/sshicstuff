#! /usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
from sshic import tools

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


"""
**************************************************************************
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def compute_telomere_freq_per_oligo_per_chr(
        df_freq: pd.DataFrame,
        df_probes: pd.DataFrame,
        table_path: str):

    all_probes = df_probes.columns.values
    res: dict = {}
    for probe in all_probes:
        fragment = df_probes.loc['frag_id', probe]
        if fragment not in df_freq.columns:
            continue
        df_freq_telo = df_freq.pivot_table(index='chr_bins', columns='chr', values=fragment, fill_value=np.nan)
        res[probe] = df_freq_telo
        df_freq_telo.to_csv(table_path + '_chr1-16_freq_on_telo.tsv', sep='\t')
    return res


def freq_focus_around_telomeres(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        window_size: int,
        telomeres_coord_path: str):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    telomeres and for each of the 16 chromosomes of yeast genome
    """

    df_centro = pd.read_csv(telomeres_coord_path, sep='\t', index_col=None)
    df_telos = pd.DataFrame({'chr': df_centro['chr'], 'telo_l': 0, 'telo_r': df_centro['length']})
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    df_probes_t = df_probes.transpose()
    bin_size = df_contacts.loc[1, 'chr_bins'] - df_contacts.loc[0, 'chr_bins']
    unique_fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])
    excluded_chr = ['chr2', 'chr3', '2_micron', 'mitochondrion', 'chr_artificial']

    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
    df_telos = df_telos[~df_telos['chr'].isin(excluded_chr)]

    df_merged = pd.merge(df_contacts, df_telos, on='chr')
    df_merged_telos_areas = df_merged[
        (df_merged.chr_bins <= (df_merged.telo_l+window_size)) |
        (df_merged.chr_bins >= (df_merged.telo_r-window_size))
    ]

    df_merged_telos_areas.loc[df_merged.chr_bins >= (df_merged.telo_r-window_size), 'chr_bins'] = \
        abs(df_merged_telos_areas['chr_bins'] - (df_merged_telos_areas['telo_r'] // bin_size)*bin_size)

    df_res = df_merged_telos_areas.groupby(['chr', 'chr_bins'], as_index=False).mean()
    df_res = tools.sort_by_chr(df_res, 'chr', 'chr_bins')
    df_res.drop(columns=['telo_l', 'telo_r'], axis=1, inplace=True)

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in unique_fragments:
        probe_chr = df_probes_t.loc[df_probes_t['frag_id'] == f, 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_res.loc[df_res['chr'] == probe_chr, f] = np.nan
        if df_res[f].sum() > 0:
            df_res[f] /= df_res[f].sum()

    return df_res, df_probes


def compute_average_aggregate(
        aggregated: dict[str: pd.DataFrame],
        table_path: str,
        plot_path: str):
    """
    After fetching the contacts for each oligos around the telomere of the 16 chr,
    we need to make an average (and std) of the 16 chr.
    """
    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()

    for probe, df in aggregated.items():
        mean = df.T.mean()
        std = df.T.std()

        ymin = -np.max((mean + std)) * 0.01
        pos = mean.index
        plt.figure(figsize=(18, 12))
        plt.bar(pos, mean)
        plt.errorbar(pos, mean, yerr=std, fmt="o", color='b', capsize=5, clip_on=True)
        plt.ylim((ymin, None))
        plt.title("Aggregated frequencies for probe {0} around telomeres".format(probe))
        plt.xlabel("Bins around the telomeres (in kb), 5' to 3'")
        plt.xticks(rotation=45)
        plt.ylabel("Average frequency made and standard deviation")
        plt.savefig(plot_path + "{0}_telomeres_aggregated_frequencies_plot.{1}".format(probe, 'jpg'), dpi=99)
        plt.close()

        df_mean[probe] = mean
        df_std[probe] = std

    #   Write to csv
    df_mean.to_csv(table_path + '_mean_on_telo.tsv', sep='\t')
    df_std.to_csv(table_path + '_std_on_telo.tsv', sep='\t')


def mkdir(output_path: str):
    dir_res = output_path + '/'
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_plot = dir_res + 'plots/'
    if not os.path.exists(dir_plot):
        os.makedirs(dir_plot)

    dir_table = dir_res + 'tables/'
    if not os.path.exists(dir_table):
        os.makedirs(dir_table)
    return dir_table, dir_plot


def run(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        telomeres_coord_path: str,
        window_size: int,
        output_path: str
):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(output_path=output_path+sample_name)

    df_contacts_centros, df_probes = freq_focus_around_telomeres(
        formatted_contacts_path=formatted_contacts_path,
        probes_to_fragments_path=probes_to_fragments_path,
        window_size=window_size,
        telomeres_coord_path=telomeres_coord_path)

    chr_aggregated_dict = compute_telomere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_probes=df_probes,
        table_path=dir_table+sample_name)

    compute_average_aggregate(
        aggregated=chr_aggregated_dict,
        table_path=dir_table,
        plot_path=dir_plot)
