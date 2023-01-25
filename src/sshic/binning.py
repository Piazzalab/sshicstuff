#! /usr/bin/env python3

import re
import pandas as pd


def run(
        formatted_contacts_path: str,
        bin_size: int,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()

    order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
             'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
             'chr15', 'chr16', '2_micron', 'mitochondrion', 'chr_artificial']

    df_formatted_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    fragments = [f for f in df_formatted_contacts if re.match(r'\d+', f)]
    df_formatted_contacts.insert(2, 'chr_bins', (df_formatted_contacts["positions"] // bin_size) * bin_size)
    df_binned_contacts = df_formatted_contacts.groupby(["chr", "chr_bins"], as_index=False).sum()
    df_binned_contacts.drop(['positions', 'sizes'], axis=1, inplace=True)

    df_binned_contacts['chr'] = df_binned_contacts['chr'].map(lambda x: order.index(x) if x in order else len(order))
    df_binned_contacts = df_binned_contacts.sort_values(by=['chr', 'chr_bins'])
    df_binned_contacts['chr'] = df_binned_contacts['chr'].map(lambda x: order[x])

    df_binned_frequencies = df_binned_contacts.copy(deep=True)
    for frag in fragments:
        df_binned_frequencies[frag] /= sum(df_binned_contacts[frag])

    df_binned_contacts.to_csv(output_dir+sample_id + '_contacts.tsv', sep='\t')
    df_binned_frequencies.to_csv(output_dir+sample_id + '_frequencies.tsv', sep='\t')
