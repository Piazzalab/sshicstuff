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
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def compute_centromere_freq_per_oligo_per_chr(
        df_freq: pd.DataFrame,
        df_probes: pd.DataFrame,
        table_path: str):

    all_probes = df_probes.columns.values
    unique_chr = pd.unique(df_freq['chr'])
    bins_array = np.unique(df_freq['chr_bins'])

    res: dict = {}
    for probe in all_probes:
        fragment = df_probes.loc['frag_id', probe]
        self_chr = df_probes.loc['chr', probe]
        if fragment not in df_freq.columns:
            continue

        df_freq_cen = df_freq.pivot_table(index='chr_bins', columns='chr', values=fragment, fill_value=np.nan)
        df_freq_cen[self_chr] = np.nan
        df_freq_cen = df_freq_cen[unique_chr].reindex(bins_array)

        res[probe] = df_freq_cen
        df_freq_cen.to_csv(table_path + probe + '_chr1-16_freq_cen.tsv', sep='\t')
    return res


def freq_focus_around_centromeres(
        formatted_contacts_path: str,
        fragments_to_oligos_path: str,
        window_size: int,
        bin_size: int,
        centros_infos_path: str,
        pooled: bool):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    centromeres and for each of the 16 chromosomes of yeast genome
    """

    df_centros = pd.read_csv(centros_infos_path, sep='\t', index_col=None)
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_probes = pd.read_csv(fragments_to_oligos_path, sep='\t', index_col=0)
    df_probes_t = df_probes.transpose()
    excluded_chr = ['chr2', 'chr3', '2_micron', 'mitochondrion']
    unique_fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    def process_row(row):
        current_chr = row[0]

        if current_chr in excluded_chr:
            return None

        current_centros_pos = row[2]
        left_cutoff = current_centros_pos - window_size - bin_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_centros_pos + window_size
        tmp_df = df_contacts.query("chr == @current_chr and chr_bins > @left_cutoff and chr_bins < @right_cutoff")

        #   temporary dataframe containing the bins present in the windows for the current chr only
        tmp_df.index = range(len(tmp_df))
        current_centros_bin = tools.find_nearest(tmp_df['chr_bins'].values, current_centros_pos, mode='lower')

        if pooled:
            tmp_df.loc[:, 'chr_bins'] = abs(tmp_df.loc[:, 'chr_bins'] - current_centros_bin)
            tmp_df = tmp_df.groupby(['chr', 'chr_bins'], as_index=False).mean()
        else:
            tmp_df.loc[:, 'chr_bins'] -= current_centros_bin

        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for f in unique_fragments:
            self_chr = df_probes_t.loc[df_probes_t['frag_id'] == f]['chr'][0]
            if self_chr == current_chr:
                tmp_df.loc[:, f] = np.nan

        return tmp_df

    df_res = pd.concat([process_row(row) for _, row in df_centros.iterrows()])
    df_res.index = range(len(df_res))
    return df_res, df_probes


def compute_average_aggregate(
        aggregated: dict[str: pd.DataFrame],
        table_path: str):
    """
    After fetching the contacts for each oligos around the centromere of the 16 chr,
    we need to make an average (and std) of the 16 chr.
    """

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    df_median = pd.DataFrame()

    for probe, df in aggregated.items():
        df_mean[probe] = df.T.mean()
        df_std[probe] = df.T.std()
        df_median[probe] = df.T.median()

    #   Write to csv
    df_mean.to_csv(table_path + 'mean_on_cen.tsv', sep='\t')
    df_std.to_csv(table_path + 'std_on_cen.tsv', sep='\t')
    df_median.to_csv(table_path + 'median_on_cen.tsv', sep='\t')

    return df_mean, df_std, df_median


def plot_aggregated(
        mean_df: pd.DataFrame,
        std_df: pd.DataFrame,
        plot_path: str):

    for probe in mean_df.columns.values:
        mean = mean_df[probe]
        std = std_df[probe]

        ymin = -np.max((mean + std)) * 0.01
        pos = mean.index
        plt.figure(figsize=(18, 12))
        plt.bar(pos, mean)
        plt.errorbar(pos, mean, yerr=std, fmt="o", color='b', capsize=5, clip_on=True)
        plt.ylim((ymin, None))
        plt.title("Aggregated frequencies for probe {0} around centromeres".format(probe))
        plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
        plt.xticks(rotation=45)
        plt.ylabel("Average frequency made and standard deviation")
        plt.savefig(plot_path + "{0}_centromeres_aggregated_freq_plot.{1}".format(probe, 'jpg'), dpi=99)
        plt.close()


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
        window_size: int,
        output_path: str,
        centros_coord_path: str):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(output_path=output_path+sample_name)

    df_contacts_centros, df_probes = freq_focus_around_centromeres(
        formatted_contacts_path=formatted_contacts_path,
        fragments_to_oligos_path=probes_to_fragments_path,
        window_size=window_size,
        bin_size=10000,
        centros_infos_path=centros_coord_path,
        pooled=True)

    chr_aggregated_dict = compute_centromere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_probes=df_probes,
        table_path=dir_table)

    df_mean, df_std, df_median = compute_average_aggregate(
        aggregated=chr_aggregated_dict,
        table_path=dir_table)

    plot_aggregated(
        mean_df=df_mean,
        std_df=df_std,
        plot_path=dir_plot)

    print('DONE: ', sample_name)
