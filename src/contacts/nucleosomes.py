#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import re
import os


def get_nfr_contacts(
        formatted_contacts_path: str,
        fragments_path: str,
        nucleosomes_path: str,
        table_path: str):

    df_fragments = pd.read_csv(fragments_path, sep='\t')
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


def nfr_statistics(
        df_contacts_in: pd.DataFrame,
        df_contacts_out: pd.DataFrame,
        output_file_name: str):

    if not os.path.exists(output_file_name):
        print("Please make sure that you already have run the statistics script")
        return None

    df_global_stats = pd.read_csv(output_file_name, sep='\t', index_col=0)
    fragments = df_global_stats['fragments'].astype(str).values

    nfr_in = []
    nfr_out = []
    nb_reads_in = df_contacts_in.shape[0]
    nb_reads_out = df_contacts_out.shape[0]
    for frag in fragments:
        nfr_in.append(
            np.sum(df_contacts_in[frag].values) / nb_reads_in
        )

        nfr_out.append(
            np.sum(df_contacts_out[frag].values) / nb_reads_out
        )

    df_global_stats['frac_nfr_in'] = nfr_in
    df_global_stats['frac_nfr_out'] = nfr_out



def fetch_fragments_sizes(
        df_fragments: pd.DataFrame,
        df_contacts: pd.DataFrame):

    chr_pos_fragments = (df_fragments['chrom'].astype(str) + '_' + df_fragments['start_pos'].astype(str)).to_numpy()
    chr_pos_contacts = (df_contacts['chr'].astype(str) + '_' + df_contacts['positions'].astype(str)).to_numpy()
    index_to_keep = np.where(np.isin(chr_pos_fragments, chr_pos_contacts))[0]
    df_fragments_filtered = df_fragments.iloc[index_to_keep]
    df_contacts.insert(2, 'size', df_fragments_filtered['size'].values)


def plot_size_distribution(
        df_contacts: pd.DataFrame,
        mode: str,
        plot_path: str):

    x = df_contacts['size'].values

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
    plt.savefig(plot_path + '_fragments_sizes_{0}_nfr_distribution.jpg'.format(mode), dpi=300)
    # plt.show()
    plt.close()


def run(
        formatted_contacts_path: str,
        statistics_path: str,
        fragments_list_path: str,
        nucleosomes_path,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_path = output_dir + sample_id

    df_fragments = pd.read_csv(fragments_list_path, sep='\t')
    dir_table, dir_plot = mkdir(output_path=output_path)
    df_contacts_in_nfr, df_contacts_out_nfr = get_nfr_contacts(
        formatted_contacts_path=formatted_contacts_path,
        fragments_path=fragments_list_path,
        nucleosomes_path=nucleosomes_path,
        table_path=dir_table+sample_id
    )

    nfr_statistics(
        df_contacts_in=df_contacts_in_nfr,
        df_contacts_out=df_contacts_out_nfr,
        output_file_name=statistics_path
    )



    # fetch_fragments_sizes(
    #     df_fragments=df_fragments,
    #     df_contacts=df_contacts_in_nfr
    # )
    #
    # fetch_fragments_sizes(
    #     df_fragments=df_fragments,
    #     df_contacts=df_contacts_out_nfr
    # )
    #
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

    print('DONE: ', sample_id)
