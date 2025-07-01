"""
This module contains functions to plot 4C-like profiles and to organize the contacts made by each probe with the genome.
"""

import os
import numpy as np
import pandas as pd

import sshicstuff.core.methods as methods
import sshicstuff.log as log

logger = log.logger


def profile_contacts(
    filtered_table_path: str,
    oligo_capture_with_frag_path: str,
    chromosomes_coord_path: str,
    normalize: bool = False,
    output_path: str = None,
    additional_groups_path: str = None,
    force: bool = False,
):
    """
    Organize the contacts made by each probe with the genome and save the results as two .tsv files:
    one for contacts and one for frequencies.

    Parameters
    ----------
    filtered_table_path : str
        Path to the filtered table (sshictuff filter script output).
    oligo_capture_with_frag_path : str
        Path to the oligo capture file (table .csv or .tsv for oligo capture information).
        Must be the file with the fragments associated made with the 'associate' command.
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

    methods.check_if_exists(filtered_table_path)
    methods.check_if_exists(oligo_capture_with_frag_path)
    methods.check_if_exists(chromosomes_coord_path)

    if not output_path:
        output_path = filtered_table_path.replace(
            "filtered.tsv", "0kb_profile_contacts.tsv"
        )

    basedir = os.path.dirname(output_path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    chr_coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords: pd.DataFrame = pd.read_csv(
        chromosomes_coord_path, sep=chr_coord_delim, index_col=None
    )
    df_coords.columns = [c.lower() for c in df_coords.columns]
    chrom_order: list[str] = df_coords["chr"].unique().tolist()

    df_chr_len = df_coords[["chr", "length"]]
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    probes = df_oligo["name"].to_list()
    fragments = df_oligo["fragment"].astype(str).to_list()

    df: pd.DataFrame = pd.read_csv(filtered_table_path, sep="\t")
    df_contacts: pd.DataFrame = pd.DataFrame(columns=["chr", "start", "sizes"])
    df_contacts: pd.DataFrame = df_contacts.astype(
        dtype={"chr": str, "start": int, "sizes": int}
    )

    for x in ["a", "b"]:
        y = methods.frag2(x)
        df2 = df[~pd.isna(df["name_" + x])]

        for probe in probes:
            if probe not in pd.unique(df2["name_" + x]):
                tmp = pd.DataFrame(
                    {
                        "chr": [np.nan],
                        "start": [np.nan],
                        "sizes": [np.nan],
                        probe: [np.nan],
                    }
                )

            else:
                df3 = df2[df2["name_" + x] == probe]
                tmp = pd.DataFrame(
                    {
                        "chr": df3["chr_" + y],
                        "start": df3["start_" + y],
                        "sizes": df3["size_" + y],
                        probe: df3["contacts"],
                    }
                )

            df_contacts = pd.concat([df_contacts, tmp])

    group = df_contacts.groupby(by=["chr", "start", "sizes"], as_index=False)
    df_contacts: pd.DataFrame = group.sum()

    # Sort the binned DataFrame
    df_contacts['chr'] = pd.Categorical(df_contacts['chr'], categories=chrom_order, ordered=True)
    df_contacts = df_contacts.sort_values(['chr', 'start']).reset_index(drop=True)

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
        df_additional: pd.DataFrame = pd.read_csv(additional_groups_path, sep="\t")
        probes_to_fragments = dict(zip(probes, fragments))
        methods.make_groups_of_probes(df_additional, df_contacts, probes_to_fragments)
        if normalize:
            methods.make_groups_of_probes(
                df_additional, df_frequencies, probes_to_fragments
            )

    df_contacts.to_csv(output_path, sep="\t", index=False)
    if normalize:
        df_frequencies.to_csv(
            output_path.replace("contacts", "frequencies"), sep="\t", index=False
        )


def profile_probes_only(
    filtered_table_path: str,
    oligo_capture_with_frag_path: str,
    output_path: str = None,
    force: bool = False,
):

    methods.check_if_exists(filtered_table_path)
    methods.check_if_exists(oligo_capture_with_frag_path)

    if not output_path:
        output_path = filtered_table_path.replace(
            "filtered.tsv", "probes_matrix.tsv"
        )

    basedir = os.path.dirname(output_path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    probes = df_oligo["name"].to_list()
    fragments = df_oligo["fragment"].astype(int).to_list()

    df: pd.DataFrame = pd.read_csv(filtered_table_path, sep="\t")
    df2 = df[df["frag_a"].isin(fragments) & df["frag_b"].isin(fragments)]

    square_matrix = np.zeros((len(probes), len(probes)))
    for i, probe_a in enumerate(probes):
        frag_a = df_oligo.loc[df_oligo["name"] == probe_a, "fragment"].values[0]
        for j, probe_b in enumerate(probes):
            frag_b = df_oligo.loc[df_oligo["name"] == probe_b, "fragment"].values[0]
            contacts = df2[(df2["frag_a"] == frag_a) & (df2["frag_b"] == frag_b)][
                "contacts"
            ].sum()
            square_matrix[i, j] = contacts

    # normalize the matrix
    total_contacts = square_matrix.sum()
    if total_contacts > 0:
        square_matrix /= total_contacts

    # make the matrix symmetric
    square_matrix = np.maximum(square_matrix, square_matrix.T)
    df_contacts = pd.DataFrame(square_matrix, columns=probes, index=probes)
    df_contacts.to_csv(output_path, sep="\t", index=True)


def rebin_profile(
    contacts_unbinned_path: str,
    chromosomes_coord_path: str,
    bin_size: int,
    output_path: str = None,
    force: bool = False,
) -> None:
    """
    Rebin the contacts profile from an unbinned contacts file into fixed-size bins.

    This function creates a template of bins based on chromosome lengths from the provided
    coordinates file and assigns contact values from the unbinned file into these bins.
    For contacts spanning multiple bins, the contact value is split proportionally based
    on the overlap. The resulting binned profile is then saved as a TSV file.

    Parameters
    ----------
    contacts_unbinned_path : str
        Path to the unbinned contacts file.
    chromosomes_coord_path : str
        Path to the chromosome coordinates file (CSV or TSV) with chromosome lengths.
    bin_size : int
        Size of the bins (resolution in bp).
    output_path : str, optional
        Path to the output file. If not provided, a default is derived from the input path.
    force : bool, optional
        If True, overwrite the output file if it exists.

    Returns
    -------
    None
    """

    # Check that input files exist (assumes methods.check_if_exists is defined)
    methods.check_if_exists(contacts_unbinned_path)
    methods.check_if_exists(chromosomes_coord_path)

    # Determine output file name if not provided
    bin_suffix = methods.get_bin_suffix(bin_size)
    if not output_path:
        output_path = contacts_unbinned_path.replace("0kb_profile", f"{bin_suffix}_profile")
    if os.path.exists(output_path) and not force:
        logger.warning("[Rebin] : Output file already exists: %s", output_path)
        logger.warning("[Rebin] : Use the --force / -F flag to overwrite the existing file.")
        return

    # Read the unbinned contacts and chromosome coordinates files
    df = pd.read_csv(contacts_unbinned_path, sep="\t")
    coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords = pd.read_csv(chromosomes_coord_path, sep=coord_delim)
    df_coords.columns = [c.lower() for c in df_coords.columns]
    chrom_order: list[str] = df_coords["chr"].unique().tolist()

    # Create a template DataFrame with all bins per chromosome
    chr_sizes = dict(zip(df_coords.chr, df_coords.length))
    chr_list_parts, chr_bins_parts = [], []
    for chrom, length in chr_sizes.items():
        n_bins = (length // bin_size) + 1
        chr_list_parts.append(np.full(n_bins, chrom))
        chr_bins_parts.append(np.arange(0, n_bins * bin_size, bin_size))
    chr_all = np.concatenate(chr_list_parts)
    bins_all = np.concatenate(chr_bins_parts)
    df_template = pd.DataFrame({
        "chr": chr_all,
        "chr_bins": bins_all,
        "genome_bins": np.arange(0, len(bins_all) * bin_size, bin_size)
    })

    # Compute the end coordinate for each contact and assign start and end bin positions
    df["end"] = df["start"] + df["sizes"]
    df["start_bin"] = (df["start"] // bin_size) * bin_size
    df["end_bin"] = (df["end"] // bin_size) * bin_size
    df.drop(columns=["genome_start"], inplace=True)

    # Separate contacts that fall entirely within one bin from those spanning bins
    df_in_bin = df.loc[df["start_bin"] == df["end_bin"]].copy()
    df_in_bin["chr_bins"] = df_in_bin["start_bin"]

    df_cross = df.loc[df["start_bin"] != df["end_bin"]].copy()
    # For cross-bin contacts, create two copies:
    #   - one for the first bin (using start_bin)
    #   - one for the second bin (using end_bin)
    df_cross_a = df_cross.copy()
    df_cross_b = df_cross.copy()
    df_cross_a["chr_bins"] = df_cross_a["start_bin"]
    df_cross_b["chr_bins"] = df_cross_b["end_bin"]

    # Identify columns corresponding to fragment data (numeric or starting with '$')
    fragments_columns = df.filter(regex=r"^\d+$|^\$").columns.tolist()
    # Compute the fraction that belongs to the second bin for cross-bin contacts
    correction = (df_cross_b["end"] - df_cross_b["chr_bins"]) / df_cross_b["sizes"]
    # Apply the proportional split in a vectorized way
    df_cross_a[fragments_columns] = df_cross_a[fragments_columns].multiply(1 - correction, axis=0)
    df_cross_b[fragments_columns] = df_cross_b[fragments_columns].multiply(correction, axis=0)

    # Combine in-bin and cross-bin fragments
    df_combined = pd.concat([df_in_bin, df_cross_a, df_cross_b], ignore_index=True)
    df_combined.drop(columns=["start_bin", "end_bin"], inplace=True)

    # Group by chromosome and bin to sum the contact values
    df_binned = df_combined.groupby(["chr", "chr_bins"], as_index=False).sum()

    # Sort the binned DataFrame
    df_binned['chr'] = pd.Categorical(df_binned['chr'], categories=chrom_order, ordered=True)
    df_binned = df_binned.sort_values(['chr', 'start']).reset_index(drop=True)

    # Merge with the template to ensure all bins are represented; missing bins are set to 0
    df_binned = pd.merge(df_template, df_binned, on=["chr", "chr_bins"], how="left")
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    # Save the final binned profile as a TSV file
    df_binned.to_csv(output_path, sep="\t", index=False)
    logger.info("[Rebin] : Binned profile saved to %s", output_path)