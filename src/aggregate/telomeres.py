#! /usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
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
        table_path: str):

    reads_array = df_info.columns.values
    chr_array = np.array(['chr'+str(i) for i in range(1, 17)])
    bins_array = pd.unique(df_freq['chr_bins'])

    res: dict = {}
    for ol in reads_array:
        probe = df_info.loc['names', ol]
        self_chr = df_info.loc['self_chr', ol]
        if len(probe.split('-/-')) > 1:
            probe = '_&_'.join(probe.split('-/-'))

        df_freq_telo = pd.DataFrame(columns=chr_array, index=bins_array)
        grouped = df_freq.groupby(['chr', 'chr_bins'])
        for name, group in grouped:
            chr_name, bin_name = name
            df_freq_telo.loc[bin_name, chr_name] = group[ol].iloc[0]

        df_freq_telo = df_freq.pivot_table(index='chr_bins', columns='chr', values=ol, fill_value=np.nan)
        df_freq_telo[self_chr] = np.nan
        df_freq_telo = df_freq_telo[chr_array].reindex(bins_array)

        res[probe] = df_freq_telo
        df_freq_telo.to_csv(table_path + '_chr1-16_freq_on_telo.tsv', sep='\t')
    return res


def freq_focus_around_centromeres(formatted_contacts_path: str,
                                  window_size: int,
                                  telomeres_coord_path: str):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    centromeres and for each of the 16 chromosomes of yeast genome
    """

    df_centro = pd.read_csv(telomeres_coord_path, sep='\t', index_col=None)
    df_telo = pd.DataFrame({'chr': df_centro['Chr'], 'telo_l': 0, 'telo_r': df_centro['Length']})
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_info, df_contacts = tools.split_formatted_dataframe(df_all)
    df_res = pd.DataFrame()

    for index, row in df_telo.iterrows():
        current_chr = row[0]
        if current_chr == '2_micron' or current_chr == 'mitochondrion':
            continue

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


def compute_average_aggregate(
        aggregated: dict[str: pd.DataFrame],
        table_path: str):
    """
    After fetching the contacts for each oligos around the telomere of the 16 chr,
    we need to make an average (and std) of the 16 chr.
    """
    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()

    for probe, df in aggregated.items():
        df_mean[probe] = df.T.mean()
        df_std[probe] = df.T.std()

    #   Write to csv
    df_mean.to_csv(table_path + '_mean_on_telo.tsv', sep='\t')
    df_std.to_csv(table_path + '_std_on_telo.tsv', sep='\t')


def pooled_stats(mean_df: pd.DataFrame,
                 std_df: pd.DataFrame):

    mean_df = mean_df.sort_index()
    std_df = std_df.sort_index()

    index_float = [float(str(m).split('_')[-1]) for m in mean_df.index]
    mean_df.index = index_float
    std_df.index = index_float

    middle = np.where(mean_df.index == 0)[0][1]
    pooled_index = mean_df.index[middle:]

    #   Pool the mean dataframe
    left_mean_df = mean_df.iloc[:middle]
    right_mean_df = mean_df.iloc[middle:]

    tmp_mean_df = pd.concat((left_mean_df, right_mean_df))
    pooled_mean_df = pd.DataFrame(tmp_mean_df.groupby(tmp_mean_df.index).mean())

    #   Pool the std dataframe
    left_std_df = std_df.iloc[:middle]
    right_std_df = std_df.iloc[middle:]
    pooled_std_df = pd.DataFrame(index=pooled_index)
    n1 = left_std_df.values
    n2 = right_std_df.values
    s1 = len(n1)
    s2 = len(n2)
    std_pooled = np.sqrt(((s1 - 1) * n1 ** 2 + (s2 - 1) * n2 ** 2) / (s1 + s2 - 2))
    pooled_std_df[0] = std_pooled
    return pooled_mean_df, pooled_std_df


def plot_aggregated(
        aggregated: dict[str: pd.DataFrame],
        plot_path: str,
        pooled: bool = True):

    for probe, df in aggregated.items():
        mean = df.T.mean()
        std = df.T.std()

        if pooled:
            mean, std = pooled_stats(mean_df=mean, std_df=std)
            mean = mean.squeeze()
            std = std.squeeze()

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


def mkdir(output_path: str):
    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_type = dir_res + '/telomeres/'
    if not os.path.exists(dir_type):
        os.makedirs(dir_type)

    dir_plot = dir_type + 'plots/'
    if not os.path.exists(dir_plot):
        os.makedirs(dir_plot)

    dir_table = dir_type + 'tables/'
    if not os.path.exists(dir_table):
        os.makedirs(dir_table)
    return dir_table, dir_plot


def run(
        formatted_contacts_path: str,
        window_size: int,
        output_path: str,
        telomeres_coord_path: str):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(output_path=output_path+sample_name)

    df_contacts_centros, df_info = freq_focus_around_centromeres(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        telomeres_coord_path=telomeres_coord_path)

    chr_aggregated_dict = compute_telomere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_info=df_info,
        table_path=dir_table+sample_name)

    compute_average_aggregate(
        aggregated=chr_aggregated_dict,
        table_path=dir_table+sample_name)

    plot_aggregated(
        aggregated=chr_aggregated_dict,
        plot_path=dir_plot+sample_name,
        pooled=True)

    print('DONE: ', sample_name)
