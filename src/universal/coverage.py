import os
import numpy as np
import pandas as pd

from universal.binning import rebin_contacts
from universal.utils import remove_columns

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def main(
        fragments_path: str,
        hic_contacts_path: str,
        output_dir: str
):

    df_fragments = pd.read_csv(fragments_path, sep='\t')
    df_fragments.rename(columns={'chrom': 'chr', 'start_pos': 'start', 'end_pos': 'end'}, inplace=True)
    df_fragments['id'] = df_fragments.index.values
    df_hic_contacts = pd.read_csv(hic_contacts_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])

    df_coverage = df_fragments[['chr', 'start', 'end']]
    df_coverage['contacts'] = np.nan

    df_merged_a = df_hic_contacts.merge(df_fragments[['id', 'chr', 'start', 'end']],
                                        left_on='frag_a',
                                        right_on='id',
                                        suffixes=('', '_a')).drop(columns=['frag_a', 'frag_b'])

    df_merged_b = df_hic_contacts.merge(df_fragments[['id', 'chr', 'start', 'end']],
                                        left_on='frag_b',
                                        right_on='id',
                                        suffixes=('', '_b')).drop(columns=['frag_a', 'frag_b'])

    df_grouped_a = df_merged_a.groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()
    df_grouped_b = df_merged_b.groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()

    df_grouped = pd.concat(
        (df_grouped_a, df_grouped_b)).groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()

    df_grouped.index = df_grouped.id
    df_grouped.drop(columns=['id'], inplace=True)
    df_grouped.to_csv(output_dir + 'coverage_per_fragment.tsv', sep='\t', index=False)
