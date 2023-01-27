#! /usr/bin/env python3
import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional
from scipy import stats
import re


def plot_size_distribution(
        df: pd.DataFrame,
        output_path: str,
        mode: str = 'all',
        bin_count: Optional[int] = None
):

    x = df['size'].values

    if bin_count is None:
        #   Freedman-Diaconis rule for optimal binning
        q1 = np.quantile(x, 0.25)
        q3 = np.quantile(x, 0.75)
        iqr = q3 - q1
        bin_width = (2 * iqr) / (len(x) ** (1 / 3))
        bin_count = int(np.ceil((x.max() - x.min()) / bin_width))

    xx = np.linspace(min(x), max(x), 10000)
    kde = stats.gaussian_kde(x)
    fig, ax = plt.subplots(figsize=(16, 14), dpi=300)
    ax.hist(x, density=True, bins=bin_count, alpha=0.3, linewidth=1.2, edgecolor='black')
    ax.plot(xx, kde(xx))
    plt.ylabel('Numbers')

    if mode == 'all':
        plt.xlabel('Sizes of NFR')
        plt.title("Distribution of NFR sizes")
        plt.savefig(output_path + 'nfr_sizes_distribution.jpg', dpi=300)
    else:
        plt.xlabel('Sizes of fragments')
        plt.title("Distribution of fragments sizes {0} NFR".format(mode))
        plt.savefig(output_path + 'fragments_sizes_{0}_nfr_distribution.jpg'.format(mode), dpi=300)
    # plt.show()
    plt.close()


def preprocess(
        fragments_list_path: str,
        fragments_nfr_filter: str,
        nucleosomes_path,
        output_dir: str
):
    """
    This function aims to part the fragments_list into two lists:
    . one with the fragments that are located at the level of a NFR
    . one with the  fragment that are not
    There is multiples ways to consider whereas a fragment overlaps correctly a NFR, that's why we try
    with multiple filters modes.

    Distribution of size sare then plot and stored :
    . histogram of fragments sizes within nfr
    . histogram of fragments sizes not in nfr
    . histogram on nfr sizes

    ARGUMENTS
    _____________

    fragments_list_path :  str
            path to the digested fragments list based on a restriction enzyme or a fixed chunk size.
    fragments_nfr_filter: str
        filter mode regarding how the fragment overlaps a nfr (the start of the fragment is
        contained into the nfr, or its end, or both, or its midpoint ..)
    nucleosomes_path : str
        path to file with all nucleosomes free regions positions and lengths
    output_dir  :  str
            the absolute path toward the output directory to save the results
    """

    df_fragments = pd.read_csv(fragments_list_path, sep='\t')
    df_fragments.insert(0, 'uid', df_fragments.index)
    df_nucleosomes = pd.read_csv(nucleosomes_path, sep='\t')
    df_nucleosomes.drop_duplicates(keep='first', inplace=True)
    df_nucleosomes.index = range(len(df_nucleosomes))
    excluded_chr = ['2_micron', 'mitochondrion', 'chr_artificial']
    df_fragments = df_fragments[~df_fragments['chrom'].isin(excluded_chr)]

    df_merged = pd.merge(df_fragments, df_nucleosomes, on='chrom')
    df_merged_in_nfr = pd.DataFrame()
    match fragments_nfr_filter:
        case 'start_only':
            df_merged_in_nfr = df_merged[
                (df_merged.start_pos > df_merged.start) & (df_merged.start_pos < df_merged.end)]
        case 'end_only':
            df_merged_in_nfr = df_merged[
                (df_merged.end_pos > df_merged.start) & (df_merged.end_pos < df_merged.end)]
        case 'middle':
            df_fragments["midpoint"] = (df_fragments["end_pos"] + df_fragments["start_pos"]) / 2
            df_merged = pd.merge(df_fragments, df_nucleosomes, on='chrom')
            df_merged_in_nfr = df_merged[
                (df_merged.midpoint >= df_merged.start) & (df_merged.midpoint <= df_merged.end)]
        case 'start_&_end':
            df_merged_in_nfr = df_merged[
                (df_merged.start_pos > df_merged.start) & (df_merged.end_pos < df_merged.end)]
    del df_merged
    index_in_nfr = np.unique(df_merged_in_nfr.uid)
    df_fragments_in_nfr = df_fragments.loc[index_in_nfr, :]
    df_fragments_out_nfr = df_fragments.drop(index_in_nfr)

    df_fragments.drop(columns='uid', inplace=True)
    df_fragments_in_nfr.to_csv(output_dir + 'fragments_list_in_nfr.tsv', sep='\t', index_label='fragments')
    df_fragments_out_nfr.to_csv(output_dir + 'fragments_list_out_nfr.tsv', sep='\t', index_label='fragments')
    df_nucleosomes = df_nucleosomes.rename(columns={'length': 'size'})

    plot_size_distribution(
        df=df_fragments_in_nfr,
        bin_count=120,
        mode='inside',
        output_path=output_dir
    )

    plot_size_distribution(
        df=df_fragments_out_nfr,
        bin_count=120,
        mode='outside',
        output_path=output_dir
    )

    plot_size_distribution(
        df=df_nucleosomes,
        bin_count=120,
        mode='all',
        output_path=output_dir
    )


def run(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        fragments_in_nfr_path: str,
        fragments_out_nfr_path: str,
        output_dir: str
):

    """
    This function's goal is to compute some statistics about the frequencies of contacts made by each
    probe, in NFR and non NFR regions.

    ARGUMENTS
    ______________
    formatted_contacts_path :  str
        path to the formatted contacts files of the current sample, basically the not_binned_contacts
        previously made with the function get_fragments_contacts in the binning script
    probes_to_fragments_path  :  str
        path to the table (.tsv file) that makes correspond the name of the probe with that of the fragment
        that contains it and complementary information, generated in the format script.
    fragments_in_nfr_path : str
        path to the list of fragments that are located on nfr
    fragments_in_nfr_path : str
        path to the list of fragments that aren't located on nfr
    output_dir  :  str
        the absolute path toward the output directory to save the results
    """

    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()

    df_fragments_in_nfr = pd.read_csv(fragments_in_nfr_path, sep='\t', index_col=0)
    df_fragments_out_nfr = pd.read_csv(fragments_out_nfr_path, sep='\t', index_col=0)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = output_dir+sample_id
    df_contacts_in = df_contacts[df_contacts['positions'].isin(df_fragments_in_nfr['start_pos'])]
    df_contacts_out = df_contacts[df_contacts['positions'].isin(df_fragments_out_nfr['start_pos'])]
    all_probes = df_probes.columns.tolist()

    probes = []
    fragments = []
    types = []
    nfr_in = []
    nfr_out = []
    total_sizes_in = sum(df_fragments_in_nfr['size'].values)
    total_sizes_out = sum(df_fragments_out_nfr['size'].values)
    total_sizes_all = total_sizes_in + total_sizes_out
    ii_probe = 0
    for probe in all_probes:
        probe_type, probe_start, probe_end, probe_chr, frag, frag_start, frag_end = df_probes[probe].tolist()
        if frag not in df_contacts.columns:
            continue
        probes.append(probe)
        fragments.append(frag)
        types.append(probe_type)

        sub_df_contacts_in = df_contacts_in.loc[df_contacts_in['chr'] != probe_chr]
        sub_df_contacts_out = df_contacts_out.loc[df_contacts_out['chr'] != probe_chr]

        #   cts_in:  sum of contacts made by the probe inside nfr
        #   cts_out: sum of contacts made by the probe outside nfr
        cts_in = np.sum(sub_df_contacts_in[frag].values)
        cts_out = np.sum(sub_df_contacts_out[frag].values)

        if cts_in + cts_out == 0:
            #   prevent from dividing by 0
            nfr_in.append(0)
            nfr_out.append(0)
        else:
            nfr_in.append(
                (cts_in / (cts_in + cts_out)) / (total_sizes_in / total_sizes_all)
            )

            nfr_out.append(
                (cts_out / (cts_in + cts_out)) / (total_sizes_out / total_sizes_all)
            )
        ii_probe += 1

    df_stats = pd.DataFrame({'probe': probes, 'fragment': fragments, 'type': types,
                             'contacts_in_nfr': nfr_in, 'contacts_out_nfr': nfr_out})

    df_stats.to_csv(output_path + '_statistics_nfr_per_probe.tsv', sep='\t')
