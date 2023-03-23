#! /usr/bin/env python3
import numpy as np
import pandas as pd
from .utils import sort_by_chr, remove_columns

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def build_bins_from_genome(
        path_to_chr_coord: str,
        bin_size: int
):
    """
    This function aims to parse the genome, and each chromosome into bins (regular range of bp)
    of size 'bin_size' a parameter given by the user as input.
    For the chr_bins, they will start at chrX_0, chrX_(0+bin_size) ... chrX_end
        idem for chrY, it starts again at chrX_0, and so on chrX_(0+bin_size) ... chrX_end
    For the genome they will start at 0, 0+bin_size, 0+2*bin_size ... genome_end.
        For the genome_bins, the positions are not reset when we go from the end of chrX to the start of chrY
    """
    df = pd.read_csv(path_to_chr_coord, sep='\t', index_col=None)
    chr_sizes = dict(zip(df.chr, df.length))
    chr_list = []
    chr_bins = []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)

    df_res = pd.DataFrame({
        'chr': chr_list,
        'chr_bins': chr_bins,
        'genome_bins': np.arange(0, len(chr_bins)*bin_size, bin_size)
    })

    return df_res


def rebin_contacts(
        df_unbinned: pd.DataFrame,
        bin_size: int,
        chromosomes_coord_path: str,
):

    df_binned_template = build_bins_from_genome(
        path_to_chr_coord=chromosomes_coord_path,
        bin_size=bin_size
    )
    df = df_unbinned.copy(deep=True)
    df.insert(2, 'chr_bins', (df["start"] // bin_size) * bin_size)
    df_binned_contacts = df.groupby(["chr", "chr_bins"], as_index=False).sum()
    df_binned_contacts = sort_by_chr(df_binned_contacts, 'chr', 'chr_bins')
    df_binned_contacts = pd.merge(df_binned_template, df_binned_contacts,  on=['chr', 'chr_bins'], how='left')
    df_binned_contacts = remove_columns(df_binned_contacts, exclusion=['start', 'end', 'size'])
    df_binned_contacts.fillna(0, inplace=True)
    return df_binned_contacts
