import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from typing import Optional, List
from utils import sort_by_chr

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


"""
**************************************************************************
**********   AGGREGATED CONTACTS AROUND CENTROMERE FUNCTIONS   **********
**************************************************************************
"""


def aggregated(
        binned_contacts_path: str,
        centros_coord_path: str,
        window_size: int,
        excluded_chr_list: List[str],
        exclude_probe_chr: bool = True,
        inter_normalization: bool = True,
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

    df_centros: pd.DataFrame = pd.read_csv(centros_coord_path, sep='\t', index_col=None)
    df_contacts: pd.DataFrame = pd.read_csv(binned_contacts_path, sep='\t')
    df_probes: pd.DataFrame = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    df_probes['frag_id'] = df_probes['frag_id'].astype(str)
    fragments: np.array = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    bin_size = df_contacts.loc[1, "chr_bins"] - df_contacts.loc[0, "chr_bins"]

    if len(excluded_chr_list) > 0:
        df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr_list)]
        df_centros = df_centros[~df_centros['chr'].isin(excluded_chr_list)]

    if exclude_probe_chr:
        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        for frag in fragments:
            probe_chr = df_probes.loc[df_probes['frag_id'] == frag, 'chr'].tolist()[0]
            if probe_chr not in excluded_chr_list:
                df_contacts.loc[df_contacts['chr'] == probe_chr, frag] = np.nan

    if inter_normalization:
        norm_suffix = "inter"
        #   Inter normalization
        df_contacts[fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
    else:
        norm_suffix = "absolute"

    df_merged: pd.DataFrame = pd.merge(df_contacts, df_centros, on='chr')
    df_merged_cen_areas: pd.DataFrame = df_merged[
        (df_merged.chr_bins > (df_merged.left_arm_length-window_size-bin_size)) &
        (df_merged.chr_bins < (df_merged.left_arm_length+window_size))
    ]

    df_merged_cen_areas['chr_bins'] = \
        abs(df_merged_cen_areas['chr_bins'] - (df_merged_cen_areas['left_arm_length'] // bin_size)*bin_size)

    df_grouped: pd.DataFrame = df_merged_cen_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)
    df_grouped = sort_by_chr(df_grouped, 'chr', 'chr_bins')
    df_grouped['chr_bins'] = df_grouped['chr_bins'].astype('int64')
    df_grouped.drop(columns=['length', 'left_arm_length', 'right_arm_length', 'genome_bins'], axis=1, inplace=True)

    df_aggregated_mean: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).mean(numeric_only=True)
    df_aggregated_mean.to_csv(
        os.path.join(dir_tables, f"aggregated_mean_contacts_around_centromeres_{norm_suffix}"), sep="\t")
    df_aggregated_std: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).std(numeric_only=True)
    df_aggregated_std.to_csv(
        os.path.join(dir_tables, f"aggregated_std_contacts_around_centromeres_{norm_suffix}"), sep="\t")
    df_aggregated_median: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).median(numeric_only=True)
    df_aggregated_median.to_csv(
        os.path.join(dir_tables, f"aggregated_median_contacts_around_centromeres_{norm_suffix}"), sep="\t")

    for probe, row in df_probes.iterrows():
        fragment = row['frag_id']
        df_chr_centros_pivot: pd.DataFrame = df_grouped.pivot_table(
            index='chr_bins', columns='chr', values=fragment, fill_value=0)
        df_chr_centros_pivot.to_csv(
            os.path.join(dir_tables, str(probe) + f"_contacts_around_centromeres_per_chr_{norm_suffix}.tsv"), sep='\t')

        if plot:
            mean = df_chr_centros_pivot.T.mean()
            std = df_chr_centros_pivot.T.std()
            median = df_chr_centros_pivot.T.median()

            ymin = -np.max((mean + std)) * 0.01
            pos = mean.index
            plt.figure(figsize=(16, 12))
            plt.bar(pos, mean)
            plt.errorbar(pos, mean, yerr=std, fmt="o", color='b', capsize=5, clip_on=True)
            plt.ylim((ymin, None))
            plt.title(f"Aggregated frequencies for probe {probe} around centromeres {norm_suffix} normalization")
            plt.xlabel("Bins around the centromeres (in kb), 5' to 3'")
            plt.xticks(rotation=45)
            plt.ylabel("Average frequency made and standard deviation")
            plt.savefig(
                os.path.join(dir_plots, f"{probe}_centromeres_aggregated_freq_plot_{norm_suffix}.jpg"), dpi=96)
            plt.close()


if __name__ == "__main__":

    aggregated(
        binned_contacts_path="../test_data/AD162/AD162_10kb_binned_frequencies.tsv",
        centros_coord_path="../test_data/S288c_chr_centro_coordinates.tsv",
        window_size=150000,
        excluded_chr_list=['chr2', 'chr3', 'chr5', '2_micron', 'mitochondrion', 'chr_artificial']
    )

    pass