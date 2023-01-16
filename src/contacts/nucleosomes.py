#! /usr/bin/env python3

import numpy as np
import pandas as pd
import re
import os


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


def run(
        formatted_contacts_path: str,
        nucleosomes_path,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_path = output_dir + sample_id

    dir_table, dir_plot = mkdir(output_path=output_path)
    df_contacts_in_nfr, df_contacts_out_nfr = get_nfr_contacts(
        formatted_contacts_path=formatted_contacts_path,
        nucleosomes_path=nucleosomes_path,
        table_path=dir_table+sample_id
    )

    print('DONE: ', sample_id)
