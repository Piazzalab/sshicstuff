#! /usr/bin/env python3

import re
import pandas as pd
import sshic.tools as tl


def run(
        formatted_contacts_path: str,
        bin_size: int,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    df_formatted_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    fragments = [f for f in df_formatted_contacts if re.match(r'\d+', f)]
    df_formatted_contacts.insert(2, 'chr_bins', (df_formatted_contacts["positions"] // bin_size) * bin_size)
    df_binned_contacts = df_formatted_contacts.groupby(["chr", "chr_bins"], as_index=False).sum()
    df_binned_contacts.drop(['positions', 'sizes'], axis=1, inplace=True)

    df_binned_contacts = tl.sort_by_chr(df_binned_contacts, 'chr', 'chr_bins')
    df_binned_frequencies = df_binned_contacts.copy(deep=True)
    for frag in fragments:
        df_binned_frequencies[frag] /= sum(df_binned_contacts[frag])

    df_binned_contacts.to_csv(output_dir+sample_id + '_contacts.tsv', sep='\t')
    df_binned_frequencies.to_csv(output_dir+sample_id + '_frequencies.tsv', sep='\t')
