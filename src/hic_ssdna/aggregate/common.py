
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

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


def plot_aggregated(
        aggregated: dict[str: pd.DataFrame],
        mode: str,
        output_path: str,
        pooled: bool = True):

    for probe, df in aggregated.items():
        mean = df.T.mean()
        std = df.T.std()

        if pooled:
            mean, std = pooled_stats(mean_df=pd.DataFrame(mean), std_df=pd.DataFrame(std))
            mean = mean.squeeze()
            std = std.squeeze()

        pos = mean.index
        plt.figure(figsize=(18, 12))
        plt.bar(pos, mean)
        plt.errorbar(pos, mean, yerr=std, fmt="o", color='b', capsize=5, clip_on=True)
        plt.ylim((-0.0001, None))
        plt.title(
            "Aggregated frequencies for probe {0} around {1}".format(probe, mode))
        plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
        plt.xticks(rotation=45)
        plt.ylabel("Average frequency made and standard deviation")
        plt.savefig(output_path + "{0}-centromeres-aggregated_frequencies_plot.{1}".format(probe, 'jpg'), dpi=99)
        plt.close()

def mkdir(output_path: str,
          mode: str):
    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    dir_type = ''
    if mode == 'centromeres':
        dir_type = dir_res + '/centromeres/'
    if mode == 'cohesins':
        dir_type = dir_res + '/cohesins_peaks/'
    if mode == 'telomeres':
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