"""
Aggregate contact data around centromeres or telomeres.
"""

import os
import re
import numpy as np
import pandas as pd
from sshicstuff.log import logger
import sshicstuff.core.methods as methods

def aggregate(
    binned_contacts_path: str,
    chr_coord_path: str,
    oligo_capture_with_frag_path: str,
    window_size: int,
    telomeres: bool = False,
    centromeres: bool = False,
    output_dir: str = None,
    excluded_chr_list: list[str] = None,
    inter_only: bool = True,
    normalize: bool = True,
    arm_length_classification: bool = False,
) -> None:
    """
    Aggregate contact data around centromeres or telomeres.

    This function merges contact data with chromosome coordinates and, if needed,
    oligo capture fragment information. It aggregates contacts within a given window 
    around centromeres or telomeres, optionally excludes intra-chromosome contacts, 
    normalizes the contact frequencies, and (if specified) classifies contacts by 
    chromosome arm lengths.

    Parameters
    ----------
    binned_contacts_path : str
        Path to the binned_contacts file (tsv).
    chr_coord_path : str
        Path to the chromosome coordinates file (tsv or csv) containing centromere positions.
    oligo_capture_with_frag_path : str
        Path to the oligo_capture file (tsv or csv) containing the associated fragments.
    window_size : int
        Aggregation window size (in base pairs) around centromeres or telomeres.
    telomeres : bool, optional
        If True, aggregate contacts around telomeres. Default is False.
    centromeres : bool, optional
        If True, aggregate contacts around centromeres. Default is False.
    output_dir : str, optional
        Output directory path. Defaults to the directory of binned_contacts_path.
    excluded_chr_list : list[str], optional
        List of chromosomes to exclude from the analysis. Default is None.
    inter_only : bool, optional
        If True, exclude intra-chromosome contacts from the analysis. Default is True.
    normalize : bool, optional
        If True, normalize the contact values by the column sums. Default is True.
    arm_length_classification : bool, optional
        If True (and telomeres is True), classify contacts by chromosome arm lengths. Default is False.

    Returns
    -------
    None
    """

    # Set default output directory if not provided
    if output_dir is None:
        output_dir = os.path.dirname(binned_contacts_path)
    # Create a subdirectory for aggregated outputs based on region type (telomeres or centromeres)
    output_dir = os.path.join(output_dir, "aggregated", "telomeres" if telomeres else "centromeres")
    os.makedirs(output_dir, exist_ok=True)

    # Extract sample names from input file name using regex
    sample_name = re.search(r".+/([^/]+)_profile", binned_contacts_path).group(1)
    sample_short_name = sample_name.split("_")[0]
    output_prefix = os.path.join(output_dir, f"{sample_short_name}_agg_on_{'telo' if telomeres else 'cen'}")

    # Determine file delimiters based on file extensions
    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    coords_delim = "\t" if chr_coord_path.endswith(".tsv") else ","

    # Read input files
    df_coords = pd.read_csv(chr_coord_path, sep=coords_delim)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    df_contacts = pd.read_csv(binned_contacts_path, sep="\t")

    # Compute binsize based on consecutive rows in the contacts file
    binsize = int(df_contacts.loc[2, "chr_bins"] - df_contacts.loc[1, "chr_bins"])
    logger.info("[Aggregate] : Contacts binned profile resolution: %d bp", binsize)

    # Get unique chromosomes from coordinates
    chr_list = df_coords["chr"].unique().tolist()
    # Ensure excluded_chr_list is a list
    if excluded_chr_list is None:
        excluded_chr_list = []

    # Get fragment identifiers as strings
    fragments = df_oligo["fragment"].astype(str).tolist()
    # Identify additional group columns (prefixed with '$')
    groups = [col for col in df_contacts.columns if col.startswith("$")]

    # Exclude chromosomes if provided
    if excluded_chr_list:
        logger.info("[Aggregate] : Excluding chromosomes: %s", ", ".join(excluded_chr_list))
        df_contacts = df_contacts[~df_contacts["chr"].isin(excluded_chr_list)]
        df_coords = df_coords[~df_coords["chr"].isin(excluded_chr_list)]

    # Exclude intra-chromosome contacts if inter_only is True
    if inter_only:
        logger.info("[Aggregate] : Excluding intra-chromosome contacts")
        # Build a mapping from chromosome (chr_ori) to a list of fragments
        chr_to_frags = {}
        for row in df_oligo.itertuples(index=False):
            frag_id = str(row.fragment)
            chr_ori = row.chr_ori
            chr_to_frags.setdefault(chr_ori, []).append(frag_id)
        # For each chromosome, set contacts for its own fragments to NaN
        for chr_name, frags in chr_to_frags.items():
            if chr_name not in excluded_chr_list:
                df_contacts.loc[df_contacts["chr"] == chr_name, frags] = np.nan
        output_prefix += "_inter"

    # Normalize contact values by the total number of contacts per fragment if requested
    if normalize:
        logger.info("[Aggregate] : Normalizing contacts by column sums")
        df_contacts.loc[:, fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
        output_prefix += "_norm"

    # Aggregate based on centromere or telomere regions
    if centromeres:
        logger.info("[Aggregate] : Aggregating contacts around centromeres")
        logger.info("[Aggregate] : Using window size: %d bp on each side", window_size)
        # Merge contacts with chromosome coordinates
        df_merged = pd.merge(df_contacts, df_coords, on="chr")
        # Select bins within the specified window around the centromere
        cen_mask = (df_merged["chr_bins"] > (df_merged["left_arm_length"] - window_size - binsize)) & \
                   (df_merged["chr_bins"] < (df_merged["left_arm_length"] + window_size))
        df_cen = df_merged.loc[cen_mask].copy()
        # Recalculate bin positions relative to the centromere
        df_cen["chr_bins"] = (df_cen["chr_bins"] - (df_cen["left_arm_length"] // binsize) * binsize).abs()
        # Group by chromosome and bin position using vectorized mean calculation
        df_grouped = df_cen.groupby(["chr", "chr_bins"], as_index=False).mean(numeric_only=True)
        # Drop unnecessary columns
        df_grouped.drop(columns=["length", "left_arm_length", "right_arm_length", "genome_bins"], inplace=True)
    elif telomeres:
        logger.info("[Aggregate] : Aggregating contacts around telomeres")
        # Create telomere boundaries DataFrame
        df_telos = pd.DataFrame({
            "chr": df_coords["chr"],
            "telo_l": 0,
            "telo_r": df_coords["length"]
        })
        df_merged = pd.merge(df_contacts, df_telos, on="chr")
        # Define regions near the left and right telomeres
        left_mask = df_merged["chr_bins"] < (df_merged["telo_l"] + window_size + binsize)
        right_mask = df_merged["chr_bins"] > (df_merged["telo_r"] - window_size - binsize)
        df_telo_left = df_merged.loc[left_mask].copy()
        df_telo_right = df_merged.loc[right_mask].copy()
        # Adjust the right telomere bins to be relative to the chromosome end
        df_telo_right["chr_bins"] = (df_telo_right["chr_bins"] - (df_telo_right["telo_r"] // binsize) * binsize).abs()
        # Concatenate the two telomere regions
        df_telos_all = pd.concat([df_telo_left, df_telo_right])
        df_grouped = df_telos_all.groupby(["chr", "chr_bins"], as_index=False).mean(numeric_only=True)
        # Drop telomere-specific columns
        df_grouped.drop(columns=["telo_l", "telo_r", "genome_bins"], inplace=True)

        # Optionally classify contacts by chromosome arm lengths
        if arm_length_classification:
            if "category" not in df_coords.columns:
                logger.error(
                    "[Aggregate] : 'category' column missing in the coordinates file. "
                    "It must be in the form 'small_small' or 'long_middle' reflecting arm lengths."
                )
            else:
                logger.info("[Aggregate] : Classifying contacts by chromosome arm lengths")
                # Build a DataFrame of arm sizes via list comprehension
                arms_data = []
                for _, row in df_coords.iterrows():
                    if row["chr"] in excluded_chr_list:
                        continue
                    if pd.isna(row["left_arm_length"]) or pd.isna(row["right_arm_length"]) or pd.isna(row["category"]):
                        continue
                    left_cat, right_cat = row["category"].split("_")
                    arms_data.append((row["chr"], "left", row["left_arm_length"], left_cat))
                    arms_data.append((row["chr"], "right", row["right_arm_length"], right_cat))
                df_arms_size = pd.DataFrame(arms_data, columns=["chr", "arm", "size", "category"])
                # Merge and classify using fixed window parameters (e.g., 3000 and 1000)
                df_merged_tel = pd.merge(df_contacts, df_telos, on="chr")
                df_left = df_merged_tel.loc[df_merged_tel["chr_bins"] < (df_merged_tel["telo_l"] + 3000 + 1000)].copy()
                df_left.insert(2, "arm", "left")
                df_right = df_merged_tel.loc[df_merged_tel["chr_bins"] > (df_merged_tel["telo_r"] - 3000 - 1000)].copy()
                df_right.insert(2, "arm", "right")
                df_telo_freq = pd.concat([df_left, df_right])
                df_merged_arm = pd.merge(df_telo_freq, df_arms_size, on=["chr", "arm"])
                df_merged_arm.drop(columns=["telo_l", "telo_r", "size"], inplace=True)
                df_grouped_by_arm = df_merged_arm.groupby("category", as_index=False).mean(numeric_only=True)
                df_grouped_by_arm.drop(columns=["chr_bins", "genome_bins"], inplace=True)
                df_grouped_by_arm = df_grouped_by_arm.rename(columns={"category": "fragments"}).T
                df_grouped_by_arm.to_csv(f"{output_prefix}_by_arm_sizes.tsv", sep="\t", header=False)
    else:
        # If neither centromeres nor telomeres are specified, exit the function.
        return

    # Sort the grouped DataFrame by chromosome and bin position using the external helper
    df_grouped = methods.sort_by_chr(df_grouped, chr_list, "chr", "chr_bins")
    df_grouped["chr_bins"] = df_grouped["chr_bins"].astype("int64")

    logger.info("[Aggregate] : Computing aggregated mean, median, and standard deviation per bin.")
    # Calculate aggregated statistics across chromosomes
    df_aggregated_mean = df_grouped.groupby("chr_bins", as_index=False).mean(numeric_only=True)
    df_aggregated_mean.to_csv(f"{output_prefix}_mean.tsv", sep="\t")
    df_aggregated_std = df_grouped.groupby("chr_bins", as_index=False).std(numeric_only=True)
    df_aggregated_std.to_csv(f"{output_prefix}_std.tsv", sep="\t")
    df_aggregated_median = df_grouped.groupby("chr_bins", as_index=False).median(numeric_only=True)
    df_aggregated_median.to_csv(f"{output_prefix}_median.tsv", sep="\t")

    # For each fragment and group column, pivot the data by chromosome if there is non-zero signal
    for col in fragments + groups:
        col_name = col if col in fragments else col[1:]
        if df_grouped[col].sum() == 0:
            continue
        df_pivot = df_grouped.pivot_table(index="chr_bins", columns="chr", values=col, fill_value=0)
        df_pivot.to_csv(f"{output_prefix}_{col_name}_per_chr.tsv", sep="\t")
