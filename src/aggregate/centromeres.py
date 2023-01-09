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
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def compute_centromere_freq_per_oligo_per_chr(
        df_freq: pd.DataFrame,
        df_info: pd.DataFrame,
        table_path: str):

    reads_array = df_info.columns.values
    chr_array = np.array(['chr'+str(i) for i in range(1, 17)])
    bins_array = np.unique(df_freq['chr_bins'])

    res: dict = {}
    for ol in reads_array:
        probe = df_info.loc['names', ol]
        self_chr = df_info.loc['self_chr', ol]
        if len(probe.split('-/-')) > 1:
            probe = '_&_'.join(probe.split('-/-'))

        df_freq_cen = df_freq.pivot_table(index='chr_bins', columns='chr', values=ol, fill_value=np.nan)
        df_freq_cen[self_chr] = np.nan
        df_freq_cen = df_freq_cen[chr_array].reindex(bins_array)

        res[probe] = df_freq_cen
        df_freq_cen.to_csv(table_path + '_' + probe + '_chr1-16_freq_cen.tsv', sep='\t')
    return res


def freq_focus_around_centromeres(formatted_contacts_path: str,
                                  window_size: int,
                                  centros_infos_path: str):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    centromeres and for each of the 16 chromosomes of yeast genome
    """

    df_centros = pd.read_csv(centros_infos_path, sep='\t', index_col=None)
    df_all = pd.read_csv(formatted_contacts_path, sep='\t', index_col=0, low_memory=False)
    df_info, df_contacts = tools.split_formatted_dataframe(df_all)
    df_res = pd.DataFrame()
    bin_size = df_contacts.iloc[1, 1] - df_contacts.iloc[0, 1]

    def process_row(row):
        current_chr = row[0]
        if current_chr == '2_micron' or current_chr == 'mitochondrion':
            return pd.DataFrame()

        current_centros_pos = row[2]
        left_cutoff = current_centros_pos - window_size - bin_size
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_centros_pos + window_size
        tmp_df = df_contacts.query("chr == @current_chr and chr_bins > @left_cutoff and chr_bins < @right_cutoff")

        #   temporary dataframe containing the bins present in the windows for the current chr only
        tmp_df.index = range(len(tmp_df))
        current_centros_bin = tools.find_nearest(tmp_df['chr_bins'].values, current_centros_pos, mode='lower')

        tmp_df.iloc[:, 1] -= current_centros_bin

        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for c in tmp_df.columns[3:]:
            self_chr = df_info.loc['self_chr', c]
            if self_chr == current_chr:
                tmp_df.loc[:, c] = np.nan

        return tmp_df

    df_res = pd.concat([process_row(row) for _, row in df_centros.iterrows()])
    df_res.index = range(len(df_res))
    return df_res, df_info


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
    df_mean.to_csv(table_path + '_mean_on_cen.tsv', sep='\t')
    df_std.to_csv(table_path + '_std_on_cen.tsv', sep='\t')
    df_median.to_csv(table_path + '_median_on_cen.tsv', sep='\t')

    return df_mean, df_std, df_median


def pooled_stats(mean_df: pd.DataFrame,
                 std_df: pd.DataFrame,
                 table_path: str):

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

    pooled_mean_df.to_csv(table_path + '_pooled_mean_on_cen.tsv', sep='\t')
    pooled_mean_df.to_csv(table_path + '_pooled_std_on_cen.tsv', sep='\t')

    return pooled_mean_df, pooled_std_df


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
        plt.savefig(plot_path + '_' + "{0}_centromeres_aggregated_freq_plot.{1}".format(probe, 'jpg'), dpi=99)
        plt.close()


def mkdir(output_path: str):
    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_type = dir_res + '/centromeres/'

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
        centros_coord_path: str):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(output_path=output_path+sample_name)

    df_contacts_centros, df_info = freq_focus_around_centromeres(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        centros_infos_path=centros_coord_path)

    chr_aggregated_dict = compute_centromere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_info=df_info,
        table_path=dir_table+sample_name)

    df_mean, df_std, df_median = compute_average_aggregate(
        aggregated=chr_aggregated_dict,
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
