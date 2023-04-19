import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import re
from typing import Optional
from utils import is_debug, sort_by_chr

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
        inter_norm: bool,
        table_path: str):

    all_probes = df_probes.index.values
    res: dict = {}
    for probe in all_probes:
        fragment = str(df_probes.loc[probe, 'frag_id'])
        if fragment not in df_freq.columns:
            continue
        if df_freq[fragment].sum() == 0.:
            continue
        df_freq_cen = df_freq.pivot_table(index='chr_bins', columns='chr', values=fragment, fill_value=0.)
        if inter_norm:
            df_freq_cen.to_csv(table_path + probe + '_chr1-16_freq_cen_inter.tsv', sep='\t')
        else:
            df_freq_cen.to_csv(table_path + probe + '_chr1-16_freq_cen_absolute.tsv', sep='\t')
        res[probe] = df_freq_cen

    return res


def freq_focus_around_centromeres(
        formatted_contacts_path: str,
        fragments_to_oligos_path: str,
        window_size: int,
        bin_size: int,
        inter_norm: bool,
        centros_infos_path: str
):
    """
    Function to capture all the bins contained in a window in bp (specified by the user), at both side of the
    centromeres and for each of the 16 chromosomes of yeast genome
    """

    df_centros = pd.read_csv(centros_infos_path, sep='\t', index_col=None)
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_probes = pd.read_csv(fragments_to_oligos_path, sep='\t', index_col=0)
    excluded_chr = ['chr2', 'chr3', 'chr5', '2_micron', 'mitochondrion', 'chr_artificial']
    fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
    df_centros = df_centros[~df_centros['chr'].isin(excluded_chr)]

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_contacts.loc[df_contacts['chr'] == probe_chr, f] = np.nan

    if inter_norm:
        #   Inter normalization
        df_contacts[fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))

    df_merged = pd.merge(df_contacts, df_centros, on='chr')
    df_merged_cen_areas = df_merged[
        (df_merged.chr_bins > (df_merged.left_arm_length-window_size-bin_size)) &
        (df_merged.chr_bins < (df_merged.left_arm_length+window_size))
    ]

    df_merged_cen_areas['chr_bins'] = \
        abs(df_merged_cen_areas['chr_bins'] - (df_merged_cen_areas['left_arm_length'] // bin_size)*bin_size)

    df_res = df_merged_cen_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)
    df_res = sort_by_chr(df_res, 'chr', 'chr_bins')
    df_res.drop(columns=['length', 'left_arm_length', 'right_arm_length'], axis=1, inplace=True)

    return df_res, df_probes


def compute_average_aggregate(
        aggregated: dict[str: pd.DataFrame],
        table_path: str,
        inter_norm: bool,
        plot: bool,
        plot_path: Optional[str]):
    """
    After fetching the contacts for each oligos around the centromere of the 16 chr,
    we need to make an average (and std) of the 16 chr.
    """
    normalization = 'absolute'
    if inter_norm:
        normalization = 'inter'

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    df_median = pd.DataFrame()

    for probe, df in aggregated.items():
        mean = df.T.mean()
        std = df.T.std()
        median = df.T.median()

        if plot:
            ymin = -np.max((mean + std)) * 0.01
            pos = mean.index
            plt.figure(figsize=(16, 12))
            plt.bar(pos, mean)
            plt.errorbar(pos, mean, yerr=std, fmt="o", color='b', capsize=5, clip_on=True)
            plt.ylim((ymin, None))
            plt.title("Aggregated frequencies for probe {0} around centromeres {1}".format(probe, normalization))
            plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
            plt.xticks(rotation=45)
            plt.ylabel("Average frequency made and standard deviation")
            plt.savefig(plot_path + "{0}_centromeres_aggregated_freq_plot_{1}.{2}".format(
                probe, normalization, 'jpg'), dpi=96)
            plt.close()

        df_mean[probe] = mean
        df_std[probe] = std
        df_median[probe] = median

    #   Write to csv
    df_mean.to_csv(table_path + 'mean_on_cen_{0}.tsv'.format(normalization), sep='\t')
    df_std.to_csv(table_path + 'std_on_cen_{0}.tsv'.format(normalization), sep='\t')
    df_median.to_csv(table_path + 'median_on_cen_{0}.tsv'.format(normalization), sep='\t')


def aggregated(
        binned_contacts_path: str,
        centros_coord_path: str,
        window_size: int,
        inter_norm: bool,
        plot: bool = True
):

    sample_id = re.search(r"AD\d+", binned_contacts_path).group()
    sample_dir = os.path.dirname(binned_contacts_path)
    data_dir = os.path.dirname(sample_dir)

    aggregated_dir = os.path.join(sample_dir, 'centromeres')
    dir_tables, dir_plots = (os.path.join(aggregated_dir, 'tables'), os.path.join(aggregated_dir, 'plots'))

    os.makedirs(aggregated_dir, exist_ok=True)
    os.makedirs(dir_plots, exist_ok=True)
    os.makedirs(dir_tables, exist_ok=True)

    probes_to_fragments_path: str = os.path.join(data_dir, "probes_to_fragments.tsv")
    if not os.path.exists(probes_to_fragments_path):
        from probe2fragment import associate_probes_to_fragments
        associate_probes_to_fragments(
            fragments_list_path=os.path.join(data_dir, "fragments_list.txt"),
            oligos_capture_path=os.path.join(data_dir, "capture_oligo_positions.csv")
        )

    df_contacts_centros, df_probes = freq_focus_around_centromeres(
        formatted_contacts_path=binned_contacts_path,
        fragments_to_oligos_path=probes_to_fragments_path,
        window_size=window_size,
        bin_size=10000,
        inter_norm=inter_norm,
        centros_infos_path=centros_coord_path,)

    chr_aggregated_dict = compute_centromere_freq_per_oligo_per_chr(
        df_freq=df_contacts_centros,
        df_probes=df_probes,
        inter_norm=inter_norm,
        table_path=dir_table)

    compute_average_aggregate(
        aggregated=chr_aggregated_dict,
        table_path=dir_table,
        inter_norm=inter_norm,
        plot=plot,
        plot_path=dir_plot)




if __name__ == "__main__":
    pass