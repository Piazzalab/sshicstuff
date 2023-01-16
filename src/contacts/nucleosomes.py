#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import re
import os


def get_fragments_sizes(fragments_to_probes_path: str):
    df_frag2oligo = pd.read_csv(fragments_to_probes_path, sep='\t', index_col=0)
    fragments_sizes = {}
    for frag in df_frag2oligo.columns.values:
        fragments_sizes[frag] = int(df_frag2oligo.loc['frag_end', frag]) - int(df_frag2oligo.loc['frag_start', frag])
    return fragments_sizes


def get_nfr_contacts(
        formatted_contacts_path: str,
        nucleosomes_path: str,
        table_path: str):

    df_nucleosomes_raw = pd.read_csv(nucleosomes_path, sep='\t')
    df_nucleosomes = df_nucleosomes_raw.drop_duplicates(keep='first')
    del df_nucleosomes_raw
    df_nucleosomes.index = range(len(df_nucleosomes))

    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
    unique_chr = pd.unique(df_contacts['chr'])
    excluded_chr = ['2_micron', 'mitochondrion']
    contacts_in_nfr: dict = {}
    df_contacts_in_nfr = pd.DataFrame()

    for chrom in unique_chr:
        if chrom in excluded_chr:
            continue
        contacts_in_nfr[chrom] = []
        positions = np.array(df_contacts.loc[df_contacts['chr'] == chrom, 'positions'])
        nucleosomes_starts = np.array(df_nucleosomes.loc[df_nucleosomes['chrom'] == chrom, 'start'])
        nucleosomes_ends = np.array(df_nucleosomes.loc[df_nucleosomes['chrom'] == chrom, 'end'])

        for pos in positions:
            for (start, end) in zip(nucleosomes_starts, nucleosomes_ends):
                if start < pos < end:
                    contacts_in_nfr[chrom].append(pos)
                    break

        df_contacts_in_nfr = pd.concat(
            [df_contacts_in_nfr, df_contacts[(
                    (df_contacts['positions'].isin(contacts_in_nfr[chrom])) &
                    (df_contacts['chr'] == chrom))]])

    df_contacts_out_nfr = df_contacts.merge(df_contacts_in_nfr, how='left', indicator=True)
    df_contacts_out_nfr = df_contacts_out_nfr[df_contacts_out_nfr['_merge'] == 'left_only']
    df_contacts_out_nfr = df_contacts_out_nfr.iloc[:, :-1]
    df_contacts_out_nfr.index = range(len(df_contacts_out_nfr))
    df_contacts_in_nfr.index = range(len(df_contacts_in_nfr))

    df_contacts_in_nfr.to_csv(table_path + '_contacts_in_nfr.tsv', sep='\t')
    df_contacts_out_nfr.to_csv(table_path + '_contacts_out_nfr.tsv', sep='\t')

    return df_contacts_in_nfr, df_contacts_out_nfr


def mkdir(output_path: str):
    dir_res = output_path
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)
    dir_plot = dir_res + '/plots/'
    if not os.path.exists(dir_plot):
        os.makedirs(dir_plot)
    dir_table = dir_res + '/tables/'
    if not os.path.exists(dir_table):
        os.makedirs(dir_table)
    return dir_table, dir_plot


def plot_size_distribution(
        df_contacts: pd.DataFrame,
        dict_sizes: dict,
        mode: str,
        plot_path: str):

    x = []
    for frag, size in dict_sizes.items():
        x.append([size]*sum(df_contacts[frag].values))
    x = np.concatenate(x)
    xx = np.linspace(min(x), max(x), 1000)
    kde = stats.gaussian_kde(x)

    fig, ax = plt.subplots(figsize=(16, 14), dpi=300)
    ax.hist(x, density=True, bins=100, alpha=0.3, linewidth=1.2, edgecolor='black')
    ax.plot(xx, kde(xx))
    plt.xlabel('Sizes of fragments')
    plt.ylabel('Numbers of contacts')
    plt.title("Distribution of fragment sizes {0} NFR".format(mode))
    plt.savefig(plot_path + '_fragment_size_{0}_nfr_distribution.jpg'.format(mode), dpi=300)
    # plt.show()
    plt.close()


def run(
        formatted_contacts_path: str,
        fragments_to_probes_path: str,
        nucleosomes_path,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_path = output_dir + sample_id

    frag_sizes = get_fragments_sizes(fragments_to_probes_path=fragments_to_probes_path)

    dir_table, dir_plot = mkdir(output_path=output_path)
    df_contacts_in_nfr, df_contacts_out_nfr = get_nfr_contacts(
        formatted_contacts_path=formatted_contacts_path,
        nucleosomes_path=nucleosomes_path,
        table_path=dir_table+sample_id
    )

    plot_size_distribution(
        df_contacts=df_contacts_in_nfr,
        dict_sizes=frag_sizes,
        mode='inside',
        plot_path=dir_plot+sample_id
    )

    plot_size_distribution(
        df_contacts=df_contacts_out_nfr,
        dict_sizes=frag_sizes,
        mode='outside',
        plot_path=dir_plot+sample_id
    )

    print('DONE: ', sample_id)
