#! /usr/bin/env python3
import os
import re
import numpy as np
import pandas as pd
from typing import Optional
from utils import sort_by_chr, make_groups_of_probes

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def build_bins_from_genome(path_to_chr_coord: str, bin_size: int):
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
        contacts_unbinned_path: str,
        chromosomes_coord_path: str,
        oligos_path: str,
        bin_size: int,
        output_dir: str,
        additional_path: Optional[str] = None,

):

    sample_filename = os.path.basename(contacts_unbinned_path)
    sample_id = re.search(r"AD\d+", sample_filename).group()
    bin_suffix = f'{bin_size // 1000}kb'
    output_path = os.path.join(output_dir, f'{sample_id}_{bin_suffix}_binned')

    df_binned_template = build_bins_from_genome(chromosomes_coord_path, bin_size)

    df_unbinned = pd.read_csv(contacts_unbinned_path, sep='\t')
    df_unbinned["end"] = df_unbinned["start"] + df_unbinned["sizes"]
    df_unbinned.drop(columns=["genome_start"], inplace=True)
    df_unbinned["start_bin"] = df_unbinned["start"] // bin_size * bin_size
    df_unbinned["end_bin"] = df_unbinned["end"] // bin_size * bin_size

    df_probes: pd.DataFrame = pd.read_csv(oligos_path, sep=',')
    probes = df_probes['name'].to_list()
    fragments = df_probes["fragment"].astype(str).tolist()
    if additional_path:
        df_additional: pd.DataFrame = pd.read_csv(additional_path, sep='\t')
        groups = df_additional['name'].to_list()
        df_unbinned.drop(columns=groups, inplace=True)
    else:
        df_additional: pd.DataFrame = pd.DataFrame()

    df_cross_bins = df_unbinned[df_unbinned["start_bin"] != df_unbinned["end_bin"]].copy()
    df_in_bin = df_unbinned.drop(df_cross_bins.index)
    df_in_bin["chr_bins"] = df_in_bin["start_bin"]

    df_cross_bins_a = df_cross_bins.copy()
    df_cross_bins_b = df_cross_bins.copy()
    df_cross_bins_a["chr_bins"] = df_cross_bins["start_bin"]
    df_cross_bins_b["chr_bins"] = df_cross_bins["end_bin"]

    fragments_columns = \
        [col for col in df_cross_bins_a.columns if col not in ["chr", "chr_bins", "start", "end", "sizes"]]

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned_contacts = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned_contacts.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned_contacts = df_binned_contacts.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned_contacts = sort_by_chr(df_binned_contacts, 'chr_bins')
    df_binned_contacts = pd.merge(df_binned_template, df_binned_contacts,  on=['chr', 'chr_bins'], how='left')
    df_binned_contacts.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned_contacts.fillna(0, inplace=True)

    df_binned_freq: pd.DataFrame = df_binned_contacts.copy(deep=True)
    for frag in fragments:
        df_binned_freq[frag] /= sum(df_binned_freq[frag])

    if additional_path:
        probes_to_fragments = dict(zip(probes, fragments))
        make_groups_of_probes(df_additional, df_binned_contacts, probes_to_fragments)
        make_groups_of_probes(df_additional, df_binned_freq, probes_to_fragments)

    df_binned_contacts.to_csv(f'{output_path}_contacts.tsv', sep='\t', index=False)
    df_binned_freq.to_csv(f'{output_path}_frequencies.tsv', sep='\t', index=False)

