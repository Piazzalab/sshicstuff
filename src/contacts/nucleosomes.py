#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from utils.tools import find_nearest
import os


def get_nfr_contacts(
        fragments_path: str,
        nucleosomes_path: str,
        output_path: str):

    df_fragments = pd.read_csv(fragments_path, sep='\t')
    df_nucleosomes = pd.read_csv(nucleosomes_path, sep='\t')
    df_nucleosomes.drop_duplicates(keep='first', inplace=True)
    df_nucleosomes.index = range(len(df_nucleosomes))
    excluded_chr = ['2_micron', 'mitochondrion', 'chr_artificial']
    df_fragments = df_fragments[~df_fragments['chrom'].isin(excluded_chr)]

    fragments_in_nfr_index = []
    for index, row in df_fragments.iterrows():
        _, chrom, start, end, size, gc = row
        if chrom in excluded_chr:
            continue
        sub_df_nucleosomes = df_nucleosomes[df_nucleosomes['chrom'] == chrom]
        nearest_nfr_start = find_nearest(sub_df_nucleosomes['start'], start, mode='lower')
        nearest_nfr_end = sub_df_nucleosomes[sub_df_nucleosomes['start'] == nearest_nfr_start]['end'].values[0]
        if nearest_nfr_start <= start < nearest_nfr_end:
            fragments_in_nfr_index.append(index)

    df_fragments_in_nfr = df_fragments.iloc[fragments_in_nfr_index]
    df_fragments_out_nfr = df_fragments.drop(fragments_in_nfr_index)

    df_fragments_in_nfr.to_csv(output_path+'fragments_list_in_nfr.tsv', sep='\t', index_label='fragments')
    df_fragments_out_nfr.to_csv(output_path + 'fragments_list_out_nfr.tsv', sep='\t', index_label='fragments')

    return df_fragments_in_nfr, df_fragments_out_nfr


def nfr_statistics(
        df_contacts: pd.DataFrame,
        df_fragments_in: pd.DataFrame,
        df_fragments_out: pd.DataFrame,
        output_path: str):

    df_contacts_in = df_contacts[df_contacts['positions'].isin(df_fragments_in['start_pos'])]
    df_contacts_out = df_contacts[df_contacts['positions'].isin(df_fragments_out['start_pos'])]
    probes = [p for p in df_contacts.columns.values if p.isdigit()]
    nfr_in = []
    nfr_out = []
    total_sizes_in = sum(df_fragments_in['size'].values)
    total_sizes_out = sum(df_fragments_out['size'].values)
    total_sizes_all = total_sizes_in + total_sizes_out
    df_stats = pd.DataFrame()
    for p in probes:
        #   cts_in:  sum of contacts made by the probe inside nfr
        #   cts_out: sum of contacts made by the probe outside nfr
        cts_in = np.sum(df_contacts_in[p].values)
        cts_out = np.sum(df_contacts_out[p].values)

        nfr_in.append(
            (cts_in / (cts_in + cts_out)) / (total_sizes_in / total_sizes_all)
        )

        nfr_out.append(
            (cts_out / (cts_in + cts_out)) / (total_sizes_out / total_sizes_all)
        )

    pass


def plot_size_distribution(
        df_fragments: pd.DataFrame,
        mode: str,
        output_path: str):

    x = df_fragments['size'].values

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
    plt.xlabel('Sizes of fragments')
    plt.ylabel('Numbers of contacts')
    plt.title("Distribution of fragments sizes {0} NFR".format(mode))
    plt.savefig(output_path + '_fragments_sizes_{0}_nfr_distribution.jpg'.format(mode), dpi=300)
    # plt.show()
    plt.close()


def run(
        formatted_contacts_path: str,
        fragments_list_path: str,
        nucleosomes_path,
        output_dir: str):

    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
    files = [f for f in os.listdir(output_dir)]
    nfr_in = 'fragments_list_in_nfr.tsv'
    nfr_out = 'fragments_list_out_nfr.tsv'
    if np.sum(np.isin([nfr_in, nfr_out], files)) == 2:
        df_fragments_in_nfr = pd.read_csv(output_dir+nfr_in, sep='\t', index_col=0)
        df_fragments_out_nfr = pd.read_csv(output_dir+nfr_out, sep='\t', index_col=0)
    else:
        df_fragments_in_nfr, df_fragments_out_nfr = get_nfr_contacts(
            fragments_path=fragments_list_path,
            nucleosomes_path=nucleosomes_path,
            output_path=output_dir
        )

    nfr_statistics(
        df_contacts=df_contacts,
        df_fragments_in=df_fragments_in_nfr,
        df_fragments_out=df_fragments_out_nfr,
        output_path=output_dir
    )

    # plot_size_distribution(
    #     df_contacts=df_contacts_in_nfr,
    #     mode='inside',
    #     plot_path=dir_plot+sample_id
    # )
    #
    # plot_size_distribution(
    #     df_contacts=df_contacts_out_nfr,
    #     mode='outside',
    #     plot_path=dir_plot+sample_id
    # )

    print('DONE')
