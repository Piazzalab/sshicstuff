#! /usr/bin/env python3
import os
import logging
import numpy as np
import pandas as pd

import sshicstuff.utils as sshcu

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')


def rebin_dataframe(
        df: pd.DataFrame,
        bin_size: int,
        df_coords: pd.DataFrame
) -> pd.DataFrame:
    """
    Rebin function for the GUI to change resolution of contacts in live mode.
    """

    chr_sizes = dict(zip(df_coords.chr, df_coords.length))
    chr_list, chr_bins = [], []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)

    df_template = pd.DataFrame({
        'chr': chr_list,
        'chr_bins': chr_bins,
        'genome_bins': np.arange(0, len(chr_bins)*bin_size, bin_size)
    })

    df["end"] = df["start"] + df["sizes"]
    df["start_bin"] = df["start"] // bin_size * bin_size
    df["end_bin"] = df["end"] // bin_size * bin_size
    df.drop(columns=["genome_start"], inplace=True)

    df_cross_bins = df[df["start_bin"] != df["end_bin"]].copy()
    df_in_bin = df.drop(df_cross_bins.index)
    df_in_bin["chr_bins"] = df_in_bin["start_bin"]

    df_cross_bins_a = df_cross_bins.copy()
    df_cross_bins_b = df_cross_bins.copy()
    df_cross_bins_a["chr_bins"] = df_cross_bins["start_bin"]
    df_cross_bins_b["chr_bins"] = df_cross_bins["end_bin"]

    fragments_columns = df.filter(regex='^\d+$').columns.to_list()

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = sshcu.sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned,  on=['chr', 'chr_bins'], how='left')
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    return df_binned


def rebin_profile(
        contacts_unbinned_path: str,
        chromosomes_coord_path: str,
        bin_size: int,
        output_path: str = None,
        force: bool = False
) -> None:

    """
    Rebin the contacts from the unbinned contacts file to the binned contacts file.

    Parameters
    ----------
    contacts_unbinned_path : str
        Path to the unbinned contacts file.
    chromosomes_coord_path : str
        Path to the chromosomes coordinates file containing the length of each chromosome arms.
    bin_size : int
        Size of the bins (resolution).
    output_path : str
        Path to the output file.
    force : bool
        Force the overwriting of the output file if the file exists.
    """

    sshcu.check_if_exists(contacts_unbinned_path)
    sshcu.check_if_exists(chromosomes_coord_path)

    bin_suffix = f'{bin_size // 1000}kb'
    if not output_path:
        output_path = contacts_unbinned_path.replace("0kb_profile", f"{bin_suffix}_profile")

    if os.path.exists(output_path) and not force:
        logging.warning(f"Output file already exists: {output_path}")
        logging.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    df = pd.read_csv(contacts_unbinned_path, sep='\t')
    coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords: pd.DataFrame = pd.read_csv(chromosomes_coord_path, sep=coord_delim, index_col=None)

    chr_sizes = dict(zip(df_coords.chr, df_coords.length))
    chr_list, chr_bins = [], []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)

    df_template = pd.DataFrame({
        'chr': chr_list,
        'chr_bins': chr_bins,
        'genome_bins': np.arange(0, len(chr_bins)*bin_size, bin_size)
    })

    df["end"] = df["start"] + df["sizes"]
    df["start_bin"] = df["start"] // bin_size * bin_size
    df["end_bin"] = df["end"] // bin_size * bin_size
    df.drop(columns=["genome_start"], inplace=True)

    df_cross_bins = df[df["start_bin"] != df["end_bin"]].copy()
    df_in_bin = df.drop(df_cross_bins.index)
    df_in_bin["chr_bins"] = df_in_bin["start_bin"]

    df_cross_bins_a = df_cross_bins.copy()
    df_cross_bins_b = df_cross_bins.copy()
    df_cross_bins_a["chr_bins"] = df_cross_bins["start_bin"]
    df_cross_bins_b["chr_bins"] = df_cross_bins["end_bin"]

    fragments_columns = df.filter(regex='^\d+$').columns.to_list()

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = sshcu.sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned,  on=['chr', 'chr_bins'], how='left')
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    df_binned.to_csv(output_path, sep='\t', index=False)

    """
    Example of usage
    
    python3 ./main.py rebin
    ../data/sandbox//AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree_0kb_profile_frequencies.tsv
    ../data/sandbox/S288c_DSB_LY_Capture_artificial_coordinates.tsv \
    -b 10000 -F
    """


def profile_contacts(
        filtered_table_path: str,
        oligos_capture_path: str,
        chromosomes_coord_path: str,
        normalize: bool = False,
        output_path: str = None,
        additional_groups_path: str = None,
        force: bool = False
):
    """
    Organize the contacts made by each probe with the genome and save the results as two .tsv files:
    one for contacts and one for frequencies.

    Parameters
    ----------
    filtered_table_path : str
        Path to the filtered table (sshictuff filter script output).
    oligos_capture_path : str
        Path to the oligos capture file (table .csv or .tsv for oligos capture information).
    chromosomes_coord_path : str
        Path to the chromosomes coordinates file containing the length of each chromosome arms.
    normalize : bool
        Normalize the contacts by the total number of contacts.
    output_path : str
        Path to the output directory.
    additional_groups_path : str
        Path to the additional groups file (table .csv or .tsv for additional groups information).
    force : bool
        Force the overwriting of the output file if the file exists.
    """
    
    sshcu.check_if_exists(filtered_table_path)
    sshcu.check_if_exists(oligos_capture_path)
    sshcu.check_if_exists(chromosomes_coord_path)

    if not output_path:
        output_path = filtered_table_path.replace("filtered.tsv", "0kb_profile_contacts.tsv")

    basedir = os.path.dirname(output_path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    if os.path.exists(output_path) and not force:
        logging.warning(f"Output file already exists: {output_path}")
        logging.warning("Use the --force / -F flag to overwrite the existing file.")
        return
    
    chr_coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords: pd.DataFrame = pd.read_csv(chromosomes_coord_path, sep=chr_coord_delim, index_col=None)
    df_chr_len = df_coords[["chr", "length"]]
    chr_list = list(df_chr_len['chr'].unique())
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    oligos_delim = "," if oligos_capture_path.endswith(".csv") else "\t"
    df_oligos: pd.DataFrame = pd.read_csv(oligos_capture_path, sep=oligos_delim)
    probes = df_oligos['name'].to_list()
    fragments = df_oligos['fragment'].astype(str).to_list()

    df: pd.DataFrame = pd.read_csv(filtered_table_path, sep='\t')
    df_contacts: pd.DataFrame = pd.DataFrame(columns=['chr', 'start', 'sizes'])
    df_contacts: pd.DataFrame = df_contacts.astype(dtype={'chr': str, 'start': int, 'sizes': int})

    for x in ['a', 'b']:
        y = sshcu.frag2(x)
        df2 = df[~pd.isna(df['name_' + x])]

        for probe in probes:
            if probe not in pd.unique(df2['name_'+x]):
                tmp = pd.DataFrame({
                    'chr': [np.nan],
                    'start': [np.nan],
                    'sizes': [np.nan],
                    probe: [np.nan]})

            else:
                df3 = df2[df2['name_'+x] == probe]
                tmp = pd.DataFrame({
                    'chr': df3['chr_'+y],
                    'start': df3['start_'+y],
                    'sizes': df3['size_'+y],
                    probe: df3['contacts']})

            df_contacts = pd.concat([df_contacts, tmp])

    group = df_contacts.groupby(by=['chr', 'start', 'sizes'], as_index=False)
    df_contacts: pd.DataFrame = group.sum()
    df_contacts = sshcu.sort_by_chr(df_contacts, chr_list, 'chr', 'start')
    df_contacts.index = range(len(df_contacts))

    for probe, frag in zip(probes, fragments):
        df_contacts.rename(columns={probe: frag}, inplace=True)

    df_contacts: pd.DataFrame = df_contacts.loc[:, ~df_contacts.columns.duplicated()]

    df_merged: pd.DataFrame = df_contacts.merge(df_chr_len, on="chr")
    df_merged["genome_start"] = df_merged["cumu_start"] + df_merged["start"]
    df_contacts.insert(3, "genome_start", df_merged["genome_start"])

    if normalize:
        df_frequencies = df_contacts.copy(deep=True)
        for frag in fragments:
            frag_sum = df_frequencies[frag].sum()
            if frag_sum > 0:
                df_frequencies[frag] /= frag_sum

    if additional_groups_path:
        df_additional: pd.DataFrame = pd.read_csv(additional_groups_path, sep='\t')
        probes_to_fragments = dict(zip(probes, fragments))
        sshcu.make_groups_of_probes(df_additional, df_contacts, probes_to_fragments)
        if normalize:
            sshcu.make_groups_of_probes(df_additional, df_frequencies, probes_to_fragments)

    df_contacts.to_csv(output_path, sep='\t', index=False)
    if normalize:
        df_frequencies.to_csv(output_path.replace("contacts", "frequencies"), sep='\t', index=False)

    """
    Example of usage
    
    python3 ./main.py profile
    ../data/sandbox/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree_filtered.tsv \
    ../data/sandbox/capture_oligo_positions.csv \
    ../data/sandbox/S288c_DSB_LY_Capture_artificial_coordinates.tsv \
    -a ../data/sandbox/additional_probe_groups.tsv \
    -F -N
    
    """
