#! /usr/bin/env python3

import re
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
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def compute_centromere_freq_per_oligo_per_chr(
        df_freq: pd.DataFrame,
        df_info: pd.DataFrame,
        dir_table: str,
        sample_id: str):

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
        df_freq_cen.to_csv(dir_table + sample_id + '_' + probe + '_chr1-16_freq_cen.tsv', sep='\t')
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
        output_file: str):
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
    df_mean.to_csv(output_file + '_mean_on_cen.tsv', sep='\t')
    df_std.to_csv(output_file + '_std_on_cen.tsv', sep='\t')
    df_median.to_csv(output_file + '_median_on_cen.tsv', sep='\t')


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

def plot_aggregated(
        aggregated: dict[str: pd.DataFrame],
        dir_plot: str,
        sample_id: str,
        pooled: bool = True):

    for probe, df in aggregated.items():
        mean = df.T.mean()
        std = df.T.std()

        if pooled:
            mean, std = pooled_stats(mean_df=pd.DataFrame(mean), std_df=pd.DataFrame(std))
            mean = mean.squeeze()
            std = std.squeeze()

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
        plt.savefig(dir_plot + sample_id + '_' + "{0}_centromeres_aggregated_freq_plot.{1}".format(probe, 'jpg'), dpi=99)
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


def debug(formatted_contacts_path: str,
          window_size: int,
          output_path: str,
          centros_coord_path: str):

    dir_table, dir_plot = mkdir(output_path=output_path)
    sample_name = re.search(r"AD\d+", output_path).group()
    output_file = dir_table + output_path.split('/')[-2]

    df_contacts_centros, df_info = freq_focus_around_centromeres(
        formatted_contacts_path=formatted_contacts_path,
        window_size=window_size,
        centros_infos_path=centros_coord_path)

    chr_aggregated_dict = compute_centromere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_info=df_info,
        dir_table=dir_table,
        sample_id=sample_name)

    compute_average_aggregate(
        aggregated=chr_aggregated_dict,
        output_file=output_file)

    plot_aggregated(
        aggregated=chr_aggregated_dict,
        dir_plot=dir_plot,
        sample_id=sample_name,
        pooled=True)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    binned_contacts_path, coordinates_path, window_size, output_path = [None for _ in range(4)]

    try:
        opts, args = getopt.getopt(argv, "h:b:c:w:o:", ["--help",
                                                        "--binning",
                                                        "--coordinates",
                                                        "--window",
                                                        "--output"])
    except getopt.GetoptError:
        print('aggregate centromeres arguments :\n'
              '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
              '-c <chr_centros_coordinates.tsv> \n'
              '-w <window> size at both side of the centromere to look around \n'
              '-o <output_file_name.tsv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('aggregate centromeres arguments :\n'
                  '-b <binned_frequencies_matrix.csv> (contacts filtered with filter.py) \n'
                  '-c <chr_centros_coordinates.tsv> \n'
                  '-w <window> size at both side of the centromere to look around \n'
                  '-o <output_file_name.tsv>')
            sys.exit()
        elif opt in ("-b", "--binning"):
            binned_contacts_path = arg
        elif opt in ("-c", "--coordinates"):
            coordinates_path = arg
        elif opt in ("-w", "--window"):
            window_size = int(arg)
        elif opt in ("-o", "--output"):
            output_path = arg

    print(output_path)
    sample_id = re.search(r"AD\d+", binned_contacts_path).group()
    output_path += sample_id + '/'

    dir_table, dir_plot = mkdir(output_path=output_path)
    output_file = dir_table + '/' + sample_id

    df_contacts_centros, df_info = freq_focus_around_centromeres(
        formatted_contacts_path=binned_contacts_path,
        window_size=window_size,
        centros_infos_path=coordinates_path)

    chr_aggregated_dict = compute_centromere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_info=df_info,
        dir_table=dir_table,
        sample_id=sample_id)

    compute_average_aggregate(
        aggregated=chr_aggregated_dict,
        output_file=output_file)

    plot_aggregated(
        aggregated=chr_aggregated_dict,
        dir_plot=dir_plot,
        sample_id=sample_id,
        pooled=True)


if __name__ == "__main__":
    #   Go into debug function if debug mode is detected, else go for main script with sys arguments
    if tools.is_debug():
        #   Debug is mainly used for testing function of the script
        #   Parameters have to be declared here
        centros_coord = "../../../../bash_scripts/aggregate_contacts/inputs/S288c_chr_centro_coordinates.tsv"

        formatted_contacts = \
            '../../../../bash_scripts/aggregate_contacts/inputs' \
            '/AD162_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30_ssHiC' \
            '_10kb_frequencies_matrix.tsv'

        output = "/home/nicolas/Documents/Projects/ssHiC/bash_scripts/aggregated_contacts/outputs/"

        oligos = "../../../../bash_scripts/aggregate_contacts/inputs/capture_oligo_positions.tsv"
        window = 150000

        samp_name = re.search(r"AD\d+", formatted_contacts).group()
        full_output_path = output + samp_name + '_aggregated'

        debug(formatted_contacts_path=formatted_contacts,
              window_size=window,
              output_path=full_output_path,
              centros_coord_path=centros_coord)

    else:
        main()

    print('--- DONE ---')
