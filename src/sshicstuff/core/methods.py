"""
Module containing functions to analyze the contacts and the capture efficiency of the oligos.
"""

import os
from os.path import join
import sys
import base64
import re
import shutil
import subprocess
import datetime

import numpy as np
import pandas as pd
import random as rd
import plotly.io as pio

import sshicstuff.log as log


logger = log.logger

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

pio.kaleido.scope.mathjax = None


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
    Aggregate the contacts around centromeres within defined regions.

    Parameters
    ----------

    binned_contacts_path : str
        Path to the binned_contacts.tsv file.
    chr_coord_path : str
        Path to the chr_centros_coordinates.tsv file containing the centromeres coordinates.
    oligo_capture_with_frag_path : str
        Path to the oligo_capture.tsv file. It must contain the fragments associated with the oligos.
        c.f associate_oligo_to_frag function.
    telomeres : bool, default=False
        Whether to aggregate the contacts around the telomeres.
    centromeres : bool, default=False
        Whether to aggregate the contacts around the centromeres.
    window_size : int
        Size of the region to aggregate around the centromeres.
    output_dir : str
        Path to the output directory.
    excluded_chr_list : list, default=None
        List of chromosome to exclude from the analysis.
    inter_only : bool, default=True
        Whether to exclude the chromosome of the probe from the analysis.
        (exclude intra-chr contacts)
    normalize : bool, default=True
        Whether to normalize the contacts by the total number of contacts remaining.
    arm_length_classification : bool, default=False
        Whether to classify the contacts by chromosome arm lengths.

    Returns
    -------
    None
    """

    if output_dir is None:
        output_dir = os.path.dirname(binned_contacts_path)

    output_dir = os.path.join(output_dir, "aggregated")
    output_dir = (
        os.path.join(output_dir, "telomeres")
        if telomeres
        else os.path.join(output_dir, "centromeres")
    )

    os.makedirs(output_dir, exist_ok=True)
    sample_name = re.match(r".+\/(.+)_profile", binned_contacts_path).group(1)
    sample_short_name = sample_name.split("_")[0]

    output_prefix = os.path.join(output_dir, sample_short_name) + "_agg_on_"
    output_prefix += "telo" if telomeres else "cen"

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    coords_delim = "\t" if chr_coord_path.endswith(".tsv") else ","
    df_coords: pd.DataFrame = pd.read_csv(
        chr_coord_path, sep=coords_delim, index_col=None
    )
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    df_contacts: pd.DataFrame = pd.read_csv(binned_contacts_path, sep="\t")

    binsize = int(df_contacts.loc[2, "chr_bins"] - df_contacts.loc[1, "chr_bins"])
    logger.info(
        "[Aggregate] : Contacts binned profile with resolution of : %d bp", binsize
    )

    chr_list = list(df_coords["chr"].unique())
    fragments = df_oligo["fragment"].astype(str).tolist()
    groups = [g for g in df_contacts.columns if g.startswith("$")]

    if len(excluded_chr_list) > 0:
        logger.info(
            "[Aggregate] : Excluding chromosomes:  %s", ", ".join(excluded_chr_list)
        )
        df_contacts = df_contacts[~df_contacts["chr"].isin(excluded_chr_list)]
        df_coords = df_coords[~df_coords["chr"].isin(excluded_chr_list)]

    if inter_only:
        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        logger.info("[Aggregate] : Excluding intra-chr contacts")
        for frag in fragments:
            ii_frag = df_oligo.loc[df_oligo["fragment"] == int(frag)].index[0]
            probe_chr_ori = df_oligo.loc[ii_frag, "chr_ori"]
            if probe_chr_ori not in excluded_chr_list:
                df_contacts.loc[df_contacts["chr"] == probe_chr_ori, frag] = np.nan

        output_prefix += "_inter"

    if normalize:
        logger.info("[Aggregate] : Normalizing the contacts")
        df_contacts.loc[:, fragments] = df_contacts[fragments].div(
            df_contacts[fragments].sum(axis=0)
        )
        output_prefix += "_norm"

    if centromeres:
        logger.info("[Aggregate] : Aggregating contacts around centromeres")
        logger.info(
            "[Aggregate] : Window size: %d bp on each side of the centromere",
            window_size,
        )

        df_merged: pd.DataFrame = pd.merge(df_contacts, df_coords, on="chr")
        df_merged_cen_areas: pd.DataFrame = df_merged[
            (df_merged.chr_bins > (df_merged.left_arm_length - window_size - binsize))
            & (df_merged.chr_bins < (df_merged.left_arm_length + window_size))
        ]
        df_merged_cen_areas["chr_bins"] = abs(
            df_merged_cen_areas["chr_bins"]
            - (df_merged_cen_areas["left_arm_length"] // binsize) * binsize
        )
        df_grouped: pd.DataFrame = df_merged_cen_areas.groupby(
            ["chr", "chr_bins"], as_index=False
        ).mean(numeric_only=True)
        df_grouped.drop(
            columns=["length", "left_arm_length", "right_arm_length", "genome_bins"],
            axis=1,
            inplace=True,
        )

    elif telomeres:
        df_telos: pd.DataFrame = pd.DataFrame(
            {"chr": df_coords["chr"], "telo_l": 0, "telo_r": df_coords["length"]}
        )
        df_merged: pd.DataFrame = pd.merge(df_contacts, df_telos, on="chr")
        df_merged_telos_areas_part_a: pd.DataFrame = df_merged[
            df_merged.chr_bins < (df_merged.telo_l + window_size + binsize)
        ]
        df_merged_telos_areas_part_b: pd.DataFrame = df_merged[
            df_merged.chr_bins > (df_merged.telo_r - window_size - binsize)
        ]
        df_merged_telos_areas_part_b["chr_bins"] = abs(
            df_merged_telos_areas_part_b["chr_bins"]
            - (df_merged_telos_areas_part_b["telo_r"] // binsize) * binsize
        )
        df_merged_telos_areas: pd.DataFrame = pd.concat(
            (df_merged_telos_areas_part_a, df_merged_telos_areas_part_b)
        )
        df_grouped: pd.DataFrame = df_merged_telos_areas.groupby(
            ["chr", "chr_bins"], as_index=False
        ).mean(numeric_only=True)
        df_grouped.drop(
            columns=["telo_l", "telo_r", "genome_bins"], axis=1, inplace=True
        )
        del df_merged, df_merged_telos_areas_part_a, df_merged_telos_areas_part_b

        if arm_length_classification:
            if "category" not in df_coords.columns:
                logger.error(
                    "[Aggregate] :"
                    "The 'category' column is missing in the centromeres file. "
                    "Must be in the form small_small or long_middle concerning lengths of left_right arms"
                )
            else:
                logger.info(
                    "[Aggregate] : Classifying the contacts by chromosome arm lengths"
                )

                df_arms_size: pd.DataFrame = pd.DataFrame(
                    columns=["chr", "arm", "size", "category"]
                )
                for _, row in df_coords.iterrows():
                    chr_ = row["chr"]
                    if chr_ not in excluded_chr_list:
                        left_, right_, category_ = (
                            row["left_arm_length"],
                            row["right_arm_length"],
                            row["category"],
                        )
                        if pd.isna(left_) or pd.isna(right_) or pd.isna(category_):
                            continue
                        df_arms_size.loc[len(df_arms_size)] = (
                            chr_,
                            "left",
                            left_,
                            category_.split("_")[0],
                        )
                        df_arms_size.loc[len(df_arms_size)] = (
                            chr_,
                            "right",
                            right_,
                            category_.split("_")[1],
                        )

                df_merged2 = pd.merge(df_contacts, df_telos, on="chr")
                df_merged_telos_areas_part_a = df_merged2[
                    df_merged2.chr_bins < (df_merged2.telo_l + 3000 + 1000)
                ]
                df_merged_telos_areas_part_a.insert(2, "arm", "left")
                df_merged_telos_areas_part_b = df_merged2[
                    df_merged2.chr_bins > (df_merged2.telo_r - 3000 - 1000)
                ]
                df_merged_telos_areas_part_b.insert(2, "arm", "right")

                df_telo_freq = pd.concat(
                    (df_merged_telos_areas_part_a, df_merged_telos_areas_part_b)
                )
                df_merged3 = pd.merge(df_telo_freq, df_arms_size, on=["chr", "arm"])
                df_merged3.drop(columns=["telo_l", "telo_r", "size"], inplace=True)

                df_grouped_by_arm = df_merged3.groupby(
                    by="category", as_index=False
                ).mean(numeric_only=True)
                df_grouped_by_arm.drop(
                    columns=["chr_bins", "genome_bins"], inplace=True
                )
                df_grouped_by_arm = df_grouped_by_arm.rename(
                    columns={"category": "fragments"}
                ).T
                df_grouped_by_arm.to_csv(
                    output_prefix + "_by_arm_sizes.tsv", sep="\t", header=False
                )

    else:
        return

    df_grouped = sort_by_chr(df_grouped, chr_list, "chr", "chr_bins")
    df_grouped["chr_bins"] = df_grouped["chr_bins"].astype("int64")

    logger.info(
        "[Aggregate] : Compute mean, median, std on the aggregated contacts per probe or group of probes, "
        "and per chromosome."
    )
    df_aggregated_mean: pd.DataFrame = df_grouped.groupby(
        by="chr_bins", as_index=False
    ).mean(numeric_only=True)
    df_aggregated_mean.to_csv(output_prefix + "_mean.tsv", sep="\t")
    df_aggregated_std: pd.DataFrame = df_grouped.groupby(
        by="chr_bins", as_index=False
    ).std(numeric_only=True)
    df_aggregated_std.to_csv(output_prefix + "_std.tsv", sep="\t")
    df_aggregated_median: pd.DataFrame = df_grouped.groupby(
        by="chr_bins", as_index=False
    ).median(numeric_only=True)
    df_aggregated_median.to_csv(output_prefix + "_median.tsv", sep="\t")

    for col in fragments + groups:
        if col in fragments:
            name = col
        else:
            name = col[1:]

        if df_grouped[col].sum() == 0:
            continue

        df_chr_centros_pivot = df_grouped.pivot_table(
            index="chr_bins", columns="chr", values=col, fill_value=0
        )
        df_chr_centros_pivot.to_csv(output_prefix + f"_{name}_per_chr.tsv", sep="\t")


def associate_oligo_to_frag(
    oligo_capture_path: str, fragments_path: str, force: bool = True
):
    """
    Associate oligo to fragments based on the fragment name.
    It adds 3 columns directly at the end of the oligo file :
    - fragment : id of the fragment, from fragment_list (hicstuff file output)
    - fragment_start : start position of the fragment
    - fragment_end : end position of the fragment

    Parameters
    ----------
    oligo_capture_path : str
        Path to the .csv file containing the oligo.
    fragments_path : str
        Path to the .csv file containing the fragments.
    force : bool
        If True, the function will overwrite the oligo file.
    Returns
    -------
    None
    """

    logger.info(
        "[Associate] : Associate oligo/probe name to fragment/read ID that contains it."
    )

    check_file_extension(fragments_path, ".txt")
    check_file_extension(oligo_capture_path, [".csv", ".tsv", ".txt"])

    output_path: str = oligo_capture_path.replace(".csv", "_fragments_associated.csv")
    logger.info(
        "[Associate] : Creating a new oligo_capture table : %s",
        output_path.split("/")[-1],
    )

    if os.path.exists(output_path) and not force:
        logger.info("[Associate] : Output file already exists: %s", output_path)
        logger.info(
            "[Associate] : Use the --force / -F flag to overwrite the existing file."
        )
        return

    # Read the oligo and fragments files
    oligo_delim = "," if oligo_capture_path.endswith(".csv") else "\t"
    df_oligo = pd.read_csv(oligo_capture_path, sep=oligo_delim)

    df_fragments = pd.read_csv(fragments_path, sep="\t")
    df_fragments["frag"] = [k for k in range(len(df_fragments))]

    fragments_id = []
    fragments_start = []
    fragments_end = []
    for _, row in df_oligo.iterrows():
        chr_ = row.iloc[0]
        probe_start = row.iloc[1]
        probe_end = row.iloc[2]
        df_sub_fragments = df_fragments[df_fragments["chrom"] == chr_]
        df_sub_fragment_sorted_start = np.sort(df_sub_fragments["start_pos"].to_numpy())

        probe_middle = int(probe_start + (probe_end - probe_start) / 2)

        idx = np.searchsorted(df_sub_fragment_sorted_start, probe_middle, side="left")
        nearest_frag_start = df_sub_fragment_sorted_start[idx - 1]

        frag_id = df_sub_fragments.index[
            df_sub_fragments["start_pos"] == nearest_frag_start
        ].tolist()[0]
        frag_start = df_sub_fragments.loc[frag_id, "start_pos"]
        frag_end = df_sub_fragments.loc[frag_id, "end_pos"]
        fragments_id.append(frag_id)
        fragments_start.append(frag_start)
        fragments_end.append(frag_end)

    df_oligo["fragment"] = fragments_id
    df_oligo["fragment_start"] = fragments_start
    df_oligo["fragment_end"] = fragments_end
    df_oligo.to_csv(output_path, sep=",", index=False)

    logger.info("[Associate] : oligos associated to fragments successfully.")


def check_file_extension(file_path: str, extension: str | list[str]):
    """
    Check if a file has the correct extension.

    Parameters
    ----------
    file_path : str
        Path to the file to check.
    extension : str
        Expected extension of the file.

    Returns
    -------
    bool
        True if the file has the correct extension, False otherwise.
    """
    if isinstance(extension, list):
        for ext in extension:
            if file_path.endswith(ext):
                return
        logger.error(f"File {file_path} does not have the correct extension {extension}.")

    else:
        if file_path.endswith(extension):
            return
        else:
            logger.error(f"File {file_path} does not have the correct extension {extension}.")


def check_gzip():
    """
    Check if gzip is installed and retrieve its version.
    """
    try:
        # Execute gzip with the --version argument to capture the version information
        result = subprocess.run(["gzip", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)
        # Parse the output to find the version number
        version_match = re.search(r"gzip (\d+\.\d+)", result.stdout)
        if version_match:
            version = version_match.group(1)
            logger.info(f"gzip version {version} is installed.")
            return version
        else:
            logger.error("Unable to determine gzip version from the output.")
            return None
    except subprocess.CalledProcessError:
        logger.error("gzip is not installed or not functioning correctly. "
                      "Please install or fix gzip before running this function.")
        return None
    except Exception as e:
        logger.error(f"Unexpected error when checking gzip version: {e}")
        return None


def check_if_exists(file_path: str):
    """
    Check if a file exists.

    Parameters
    ----------
    file_path : str
        Path to the file to check.

    Returns
    -------
    bool
        True if the file exists, False otherwise.
    """
    if os.path.exists(file_path):
        return
    else:
        logger.error(f"File {file_path} does not exist.")
        sys.exit(1)


def check_seqtk():
    """
    Check if seqtk is installed and retrieve its version.
    """
    try:
        # Execute seqtk without arguments to capture the usage information which contains the version
        result = subprocess.run("seqtk", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False, text=True)
        # Parse the output to find the version number
        version_match = re.search(r"Version: (\S+)", result.stdout)
        if version_match:
            version = version_match.group(1)
            logger.info(f"seqtk version {version} is installed.")
            return version
        else:
            logger.error("Unable to determine seqtk version from the output.")
            sys.exit(1)
    except subprocess.CalledProcessError:
        logger.error("seqtk is not installed or not functioning correctly. "
                      "Please install or fix seqtk before running this function.")
        sys.exit(1)


def compare_with_wt(
    stats1_path: str, stats2_path: str, ref_name: str, output_dir: str = None
) -> None:
    """
    Compare the capture efficiency of a sample with a wild-type reference.

    Gives a .csv file with the with the ratio of the capture efficiency
    of the sample over the wild-type reference.


    Parameters
    ----------
    stats1_path : str
        Path to the statistics file of the sample.
    stats2_path : str
        Path to the statistics file of the wild-type reference.
    ref_name : str
        Name of the wild-type reference.
    output_dir : str
        Path to the output directory.

    Returns
    -------
    None
    """

    logger.info(
        "Comparing the capture efficiency of a sample with a wild-type reference."
    )
    logger.info(
        "Be sure to have same number of reads for both samples. Otherwise use subsample function."
    )

    df_sample: pd.DataFrame = pd.read_csv(stats1_path, header=0, sep="\t")
    df_wt: pd.DataFrame = pd.read_csv(stats2_path, sep="\t")

    df_cap_eff = pd.DataFrame(
        columns=[
            "probe",
            "capture_efficiency",
            f"capture_efficiency_{ref_name}",
            "ratio_sample_vs_wt",
        ]
    )

    for index, row in df_sample.iterrows():
        probe = row["probe"]
        cap_eff = row["dsdna_norm_capture_efficiency"]

        if probe in df_wt["probe"].tolist():
            cap_eff_wt = df_wt.loc[
                df_wt["probe"] == probe, "dsdna_norm_capture_efficiency"
            ].tolist()[0]
            ratio = cap_eff / cap_eff_wt if cap_eff_wt != 0 else np.nan

        else:
            cap_eff_wt = np.nan
            ratio = np.nan

        df_cap_eff.loc[index, "probe"] = probe
        df_cap_eff.loc[index, "capture_efficiency"] = cap_eff
        df_cap_eff.loc[index, f"capture_efficiency_{ref_name}"] = cap_eff_wt
        df_cap_eff.loc[index, f"ratio_sample_vs_{ref_name}"] = ratio

    if output_dir is None:
        output_dir = os.path.dirname(stats1_path)

    output_path = os.path.join(
        output_dir, f"{os.path.basename(stats1_path).split('.')[0]}_vs_{ref_name}.csv"
    )
    df_cap_eff.to_csv(output_path, sep="\t", index=False)


def copy(source_path, destination_path):
    """
    Copy a file from source to destination.
    Useful to copy inputs files to the output directory in order to have a trace.

    Parameters
    ----------
    source_path : str
        Path to the file to copy.
    destination_path : str
        Path to the destination directory.
    """
    try:
        shutil.copy(source_path, destination_path)
        src_basename = source_path.split('/')[-1]
        logger.info(f"[Copy] : {src_basename} copied.")
    except IOError as e:
        logger.error(f"Unable to copy file. Error: {e}")
        sys.exit(1)


def coverage(
    sparse_mat_path: str,
    fragments_list_path: str,
    output_dir: str = None,
    normalize: bool = False,
    force: bool = False,
    bin_size: int = 0,
    chromosomes_coord_path: str = "",
) -> None:
    """
    Calculate the coverage per fragment and save the result to a bedgraph file in the output directory.

    Parameters
    ----------
    sparse_mat_path : str
        Path to the sparse_contacts_input.txt file (generated by hicstuff).
    fragments_list_path : str
        Path to the fragments_input.txt file (generated by hicstuff).
    output_dir : str
        Path to the output directory.
    normalize : bool
        Normalize the coverage by the total number of contacts.
    force : bool
        Force the overwriting of the output file if the file exists.
    bin_size : int
        Size of the bin to use for the coverage bedgraph.
    chromosomes_coord_path : str
        Path to the chromosomes coordinates file containing the length of each chromosome arms.
        Needed for the binning.

    Returns
    -------
    None
    """

    logger.info("[Coverage] : Calculating coverage per fragment into a bedgraph.")

    if output_dir is None:
        output_dir = os.path.dirname(sparse_mat_path)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(
        output_dir, os.path.basename(sparse_mat_path).split(".")[0]
    )

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    df_fragments: pd.DataFrame = pd.read_csv(fragments_list_path, sep="\t")
    df_fragments.rename(
        columns={"chrom": "chr", "start_pos": "start", "end_pos": "end"}, inplace=True
    )
    df_fragments["id"] = list(range(len(df_fragments)))

    df_hic_contacts: pd.DataFrame = pd.read_csv(
        sparse_mat_path, header=0, sep="\t", names=["frag_a", "frag_b", "contacts"]
    )

    df_merged_a: pd.DataFrame = df_hic_contacts.merge(
        df_fragments[["id", "chr", "start", "end"]],
        left_on="frag_a",
        right_on="id",
        suffixes=("", "_a"),
    ).drop(columns=["frag_a", "frag_b"])

    df_merged_b: pd.DataFrame = df_hic_contacts.merge(
        df_fragments[["id", "chr", "start", "end"]],
        left_on="frag_b",
        right_on="id",
        suffixes=("", "_b"),
    ).drop(columns=["frag_a", "frag_b"])

    df_grouped_a: pd.DataFrame = df_merged_a.groupby(
        by=["id", "chr", "start", "end"], as_index=False
    ).sum()
    df_grouped_b: pd.DataFrame = df_merged_b.groupby(
        by=["id", "chr", "start", "end"], as_index=False
    ).sum()

    df_contacts_cov: pd.DataFrame = (
        pd.concat((df_grouped_a, df_grouped_b))
        .groupby(by=["id", "chr", "start", "end"], as_index=False)
        .sum()
    )

    df_contacts_cov.index = df_contacts_cov.id
    df_contacts_cov.drop(columns=["id"], inplace=True)
    output_path = output_path + "_contacts_coverage.bedgraph"

    if bin_size > 0:
        # Define output file name with bin size suffix
        bin_suffix = get_bin_suffix(bin_size)
        output_path = output_path.replace(".bedgraph", f"_{bin_suffix}.bedgraph")
        logger.info("[Coverage] : Binning the bedgraph at %d resolution.", bin_size)

        # Load chromosome sizes
        check_if_exists(chromosomes_coord_path)
        delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
        df_chrom = pd.read_csv(chromosomes_coord_path, sep=delim)
        chrom_sizes = dict(zip(df_chrom.chr, df_chrom.length))

        # Create empty bins for all chromosomes
        chr_list, chr_bins = [], []
        for c, l in chrom_sizes.items():
            chr_list.append([c] * (l // bin_size + 1))
            chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))
        chr_list = np.concatenate(chr_list)
        chr_bins = np.concatenate(chr_bins)

        df_bins = pd.DataFrame(
            {
                "chr": chr_list,
                "start": chr_bins,
                "end": chr_bins + bin_size,
                "contacts": 0.0,
            }
        )

        df = df_contacts_cov.copy()
        df["size"] = df["end"] - df["start"]
        df["start_bin"] = df["start"] // bin_size * bin_size
        df["end_bin"] = df["end"] // bin_size * bin_size

        # Separate contacts spanning multiple bins
        df_cross = df[df["start_bin"] != df["end_bin"]].copy()
        df_in_bin = df.drop(df_cross.index)

        # Adjust contacts for bins that span two bins
        df_a, df_b = df_cross.copy(), df_cross.copy()
        df_a["start_bin"], df_b["start_bin"] = (
            df_cross["start_bin"],
            df_cross["end_bin"],
        )
        factor_b = (df_b["end"] - df_b["start_bin"]) / df_b["size"]
        df_a["contacts"] *= 1 - factor_b
        df_b["contacts"] *= factor_b

        # Merge corrected bins
        df_corrected = pd.concat([df_in_bin, df_a, df_b])
        df_corrected.drop(columns=["size", "start", "end", "end_bin"], inplace=True)
        df_corrected.rename(columns={"start_bin": "start"}, inplace=True)
        df_corrected = df_corrected.groupby(["chr", "start"]).sum().reset_index()
        df_corrected["end"] = df_corrected["start"] + bin_size
        df_corrected["contacts"] = np.round(df_corrected["contacts"], 4)

        # Merge with empty bins
        df_final = pd.concat([df_bins, df_corrected])
        df_final = df_final.groupby(["chr", "start", "end"]).sum().reset_index()
        df_final = sort_by_chr(df_final, chr_list, "chr", "start")
        df_final["contacts"] = df_final["contacts"].fillna(0)

        # Save output
        df_final.to_csv(output_path, sep="\t", index=False, header=False)
        logger.info(
            "[Coverage] : Contacts coverage binned file saved to %s", output_path
        )
        df_contacts_cov = df_final

    else:
        df_contacts_cov.to_csv(output_path, sep="\t", index=False, header=False)
        logger.info("[Coverage] : Contacts coverage file saved to %s", output_path)

    if normalize:
        output_path = output_path.replace("_contacts_", "_frequencies_")
        logger.info(
            "[Coverage] : Normalizing coverage by the total number of contacts."
        )
        df_frequencies_cov: pd.DataFrame = df_contacts_cov.copy(deep=True)
        df_frequencies_cov["contacts"] /= sum(df_frequencies_cov["contacts"])
        df_frequencies_cov.rename(columns={"contacts": "frequencies"})
        df_frequencies_cov.to_csv(output_path, sep="\t", index=False, header=False)

    logger.info("[Coverage] : Coverage calculation completed.")


def detect_delimiter(path: str):
    """
    Detect the delimiter of a file.
    The delimiter is detected by counting the number of tabs and commas in the file.

    Parameters
    ----------
    path : str
        Path to the file.

    Returns
    -------
    str
        Delimiter of the file.
    """

    with open(path, 'r', encoding='utf-8') as file:
        contents = file.read()
    tabs = contents.count('\t')
    commas = contents.count(',')
    if tabs > commas:
        return '\t'
    else:
        return ','


def edit_genome_ref(
    annealing_input: str,
    genome_input: str,
    enzyme: str,
    fragment_size: int = 150,
    fasta_spacer: str = "N",
    fasta_line_length: int = 80,
    additional_fasta_path: str = None,
):
    """
    Create an artificial chromosome that is the concatenation of the annealing oligo and the enzyme sequence.

    Insert it at the end of the original genome .FASTA file.

    Parameters
    ----------

    annealing_input : str
        Path to the annealing oligo input CSV file.
    genome_input : str
        Path to the original genome .FASTA file.
    enzyme : str
        Restriction Enzyme sequence (e.g., dpnII sequence : gatc).
    fragment_size : int, default=150
        Size of a digested fragment / read.
    fasta_spacer : str, default="N"
        Spacer character to insert between the enzyme and the annealing oligo.
    fasta_line_length : int, default=60
        Number of characters per line in the FASTA file.
    additional_fasta_path : str, default=None
        List of additional FASTA files to concatenate with the artificial chromosome ath
        the end of the genome reference .FASTA file.
    """

    basedir = os.path.dirname(genome_input)
    artificial_chr_path = os.path.join(basedir, "chr_artificial_ssDNA.fa")

    # Creating the artificial chromosome using annealing oligo sequences
    # and the enzyme sequence
    logger.info(
        "Creating the artificial chromosome with the annealing oligo and the enzyme %s",
        enzyme,
    )

    df = pd.read_csv(annealing_input, sep=",")
    ssdna_seq_series = df[df["type"] == "ss"]["sequence_modified"]
    ssdna_seq = [seq.lower() for seq in ssdna_seq_series.values]

    lg = fasta_line_length
    s = fragment_size - len(enzyme)
    p = fasta_spacer
    oneline = p * int(s / 2) + enzyme + p * s

    for seq in ssdna_seq:
        middle = len(seq) // 2
        enzyme_pos = seq.find(enzyme)
        if enzyme_pos < middle:
            seq2 = seq[enzyme_pos + len(enzyme) :].upper()
        else:
            seq2 = seq[:enzyme_pos].upper()

        oneline += seq2 + p * s + enzyme + p * s

    lines = "\n".join([oneline[i : i + lg] for i in range(0, len(oneline), lg)])
    fasta = f">chr_artificial_ssDNA\t ({len(oneline)} bp)\n{lines}"

    with open(artificial_chr_path, "w", encoding="utf-8") as f:
        f.write(fasta)

    # Inserting the artificial chromosome at the end of the genome .FASTA file
    genome_name = os.path.basename(genome_input)
    logger.info(
        "Inserting the artificial chromosome at the end of the original genome .FASTA file"
    )
    with open(genome_input, "r", encoding="utf-8") as f:
        genome = f.read()

    new_genome = genome + "\n" + fasta + "\n"

    # Concatenate with additional FASTA sequence(s), if any
    if additional_fasta_path:
        logger.info(
            "Looking for additional FASTA sequence(s) to concatenate with %s",
            genome_name,
        )
        add_fasta_name = os.path.basename(additional_fasta_path)
        logger.info("Concatenating %s with the genome .FASTA file", add_fasta_name)
        with open(additional_fasta_path, "r", encoding="utf-8") as f:
            add_fasta = f.read()

        # Check line length
        if len(add_fasta.split("\n")[1]) != lg:
            logger.warning("Line length of %s is not equal to %d", add_fasta_name, lg)

            # remove existing line breaks and add new ones
            add_fasta = add_fasta.replace("\n", "")
            add_fasta = "\n".join(
                [add_fasta[i : i + lg] for i in range(0, len(add_fasta), lg)]
            )

        # Concatenate the strings
        new_genome += "\n" + add_fasta

    new_genome_output = genome_input.replace(".fa", "_artificial.fa")
    with open(new_genome_output, "w", encoding="utf-8") as f:
        f.write(new_genome)

    logger.info(
        "Artificial chromosome created and inserted at the end of the genome .FASTA file"
    )


def frag2(x):
    """
    if x = a get b, if x = b get a

    Parameters
    ----------
    x : str
        String to invert.

    Returns
    -------
    str
        Inverted string.
    """
    if x == 'a':
        y = 'b'
    else:
        y = 'a'
    return y


def get_bin_suffix(bin_size: int):
    """
    Get the suffix of a bin size.
    The suffix is the bin size in kb or Mb if the bin size is greater than 1e3 or 1e6 respectively.
    """

    p = int(np.log10(bin_size))

    if p < 3:
        res = f"{bin_size}bp"
    elif p < 6:
        res = f"{bin_size // int(1e3)}kb"
    else:
        res = f"{bin_size // int(1e6)}Mb"

    return res


def generate_colors(color_type: str, n: int, a: float = 0.8, seed: int = 42) -> list[str]:
    """
    Generate random colors in hexadecimal or rgba format.
    """

    rd.seed(seed)

    if color_type == 'hex':
        white = '#FFFFFF'
        blue = '#0080FF'
        res = [f"#{rd.randint(0, 0xFFFFFF):06x}" for _ in range(n)]
        if white in res:
            res[res.index(white)] = blue

        return res

    elif color_type == 'rgba':
        white = 'rgba(255, 255, 255, 0.8)'
        blue = 'rgba(0, 128, 255, 0.8)'

        res = [
            (f"rgba({rd.randint(0, 255)}, "
             f"{rd.randint(0, 255)}, "
             f"{rd.randint(0, 255)}, {a})") for _ in range(n)
        ]

        if white in res:
            res[res.index(white)] = blue

        return res

    else:
        raise ValueError("type must be 'hex' or 'rgba'")


def is_debug() -> bool:
    """
    Check if the script is running in debug mode.

    Returns
    -------
    bool
        True if the script is running in debug mode, False otherwise.
    """
    gettrace = getattr(sys, 'gettrace', None)

    if gettrace is None:
        return False
    else:
        v = gettrace()
        if v is None:
            return False
        else:
            return True


def make_groups_of_probes(df_groups: pd.DataFrame, df: pd.DataFrame, prob2frag: dict):
    for _, row in df_groups.iterrows():
        group_probes = row["probes"].split(",")
        group_frags = np.unique([prob2frag[probe] for probe in group_probes])
        group_name = row["name"]
        group_name = "$" + group_name.lower()
        if row["action"] == "average":
            df[group_name] = df[group_frags].mean(axis=1)
        elif row["action"] == "sum":
            df[group_name] = df[group_frags].sum(axis=1)
        else:
            continue


def merge_sparse_mat(
    output_path: str = None, force: bool = False, matrices: list[str] = None
) -> None:
    if not matrices:
        logger.error("No sparse matrices provided")
        return

    n = len(matrices)
    logger.info("Merging %d sparse matrices into one", n)

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    for i in range(n):
        check_file_extension(matrices[i], ".txt")

    if not output_path:
        now_ = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        output_path = os.path.join(
            os.path.dirname(matrices[0]), f"{now_}_merged_sparse_contacts.tsv"
        )

    df_sparses: list[pd.DataFrame] = [
        pd.read_csv(matrix, sep="\t", header=0) for matrix in matrices
    ]

    n_frags = int(df_sparses[0].columns[0])
    if not all([int(df.columns[0]) == n_frags for df in df_sparses]):
        logger.error("All the sparse matrices must have the same number of fragments")
        return

    else:
        for df in df_sparses:
            df.columns = ["frag_a", "frag_b", "contacts"]

    df_concatenated: pd.DataFrame = pd.concat(df_sparses)
    df_merged = df_concatenated.groupby(["frag_a", "frag_b"], as_index=False)[
        "contacts"
    ].sum()

    df_merged.columns = [n_frags, n_frags, len(df_merged) + 1]
    df_merged.to_csv(output_path, sep="\t", index=False, header=True)
    logger.info("Merged sparse matrix saved to %s", output_path)


def save_file_cache(name, content, cache_dir):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(join(cache_dir, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def sort_by_chr(df: pd.DataFrame, chr_list: list[str], *args: str):
    """
    Sort a DataFrame by chromosome and then by other columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to sort.
    chr_list : List[str]
        List of chromosomes.
    args : str
        Columns to sort by after the chromosome.

    Returns
    -------
    pd.DataFrame
        Sorted DataFrame.
    """
    # use re to identify chromosomes of the form "chrX" with X being a number
    chr_with_number = [c for c in chr_list if re.match(r'chr\d+', c)]
    chr_with_number.sort(key=lambda x: int(x[3:]))
    chr_without_number = [c for c in chr_list if c not in chr_with_number]

    order = chr_with_number + chr_without_number
    df['chr'] = df['chr'].apply(lambda x: order.index(x) if x in order else len(order))

    if args:
        df = df.sort_values(by=['chr', *args])
    else:
        df = df.sort_values(by=['chr'])

    df['chr'] = df['chr'].map(lambda x: order[x])
    df.index = range(len(df))

    return df


def sparse_with_dsdna_only(
    sample_sparse_mat: str,
    oligo_capture_with_frag_path: str,
    n_flanking_dsdna: int = 2,
    output_path: str = None,
    force: bool = False,
) -> None:
    """
    Create a contact sparse matrix file with the same format as the input file, but with no ssDNA.
    i.e., it removes fragments containing a probe and the N fragment up/downstream of it.

    Parameters
    ----------
    sample_sparse_mat : str
        Path to the sparse matrix file (hicstuff given output).

    oligo_capture_with_frag_path : str
        Path to the oligo capture file (sshicstuff mandatory table).
        Must be the file with the fragments associated made with the 'associate' command

    n_flanking_dsdna : int
        Number of flanking fragments to remove around the probe fragment.
        Default is 2.

    output_path : str
        Path to the output file to be created.
        Default is None.

    force : bool
        Force the overwriting of the oligo file even if the columns are already present.
        Default is True.

    Returns
    -------
    None
    """

    if not output_path:
        output_path = sample_sparse_mat.replace(".txt", "_dsdna_only.txt")

    if not force and os.path.exists(output_path):
        logger.info(
            "[Sparse Matrix Graal (dsdna)] Output file already exists: %s", output_path
        )
        logger.warning(
            "[Sparse Matrix Graal (dsdna)] Use the --force / -F flag to overwrite the existing file."
        )
        return

    check_if_exists(sample_sparse_mat)
    check_if_exists(oligo_capture_with_frag_path)
    check_file_extension(sample_sparse_mat, ".txt")
    check_file_extension(oligo_capture_with_frag_path, [".csv", ".tsv"])

    oligo_capture_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep="\t", header=None)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_capture_delim)

    df_contacts_dsdna_only = df_sparse_mat.copy(deep=True)

    ssdna_frag = df_oligo.loc[df_oligo["type"] == "ss", "fragment"].tolist()
    df_ssdna = pd.DataFrame(ssdna_frag, columns=["fragments"])

    dsdna_frag = df_oligo.loc[df_oligo["type"] == "ds", "fragment"].tolist()
    dsdna_frag_flanking = []
    for f in dsdna_frag:
        for i in range(1, n_flanking_dsdna + 1):
            dsdna_frag_flanking.append(f + i)
            dsdna_frag_flanking.append(f - i)

    dsdna_frag_all = np.unique(dsdna_frag + dsdna_frag_flanking)
    df_dsdna = pd.DataFrame(dsdna_frag_all, columns=["fragments"])

    df_frag = pd.concat([df_ssdna, df_dsdna])
    del df_ssdna, df_dsdna

    df_sparse_mat["index"] = df_sparse_mat.index
    matches_a = pd.merge(
        df_sparse_mat,
        df_frag,
        left_on=0,
        right_on="fragments",
        how="inner",
        indicator=True,
    )
    matches_b = pd.merge(
        df_sparse_mat,
        df_frag,
        left_on=1,
        right_on="fragments",
        how="inner",
        indicator=True,
    )
    index_to_drop = np.unique(
        np.concatenate((matches_a["index"].to_numpy(), matches_b["index"].to_numpy()))
    )

    df_contacts_dsdna_only.drop(index_to_drop, inplace=True)

    df_contacts_dsdna_only.iloc[0, 0] -= len(df_frag)
    df_contacts_dsdna_only.iloc[0, 1] -= len(df_frag)
    df_contacts_dsdna_only.iloc[0, 2] -= len(index_to_drop)

    df_contacts_dsdna_only.to_csv(output_path, sep="\t", index=False, header=False)

    logger.info(
        "[Sparse Matrix Graal (dsdna)] : dsDNA only contacts saved to %s", output_path
    )


def sparse_with_ssdna_only(
    sample_sparse_mat: str,
    oligo_capture_with_frag_path: str,
    output_path: str = None,
    force: bool = False,
) -> None:
    """
    Create a contact sparse matrix file with the same format as the input file, but with only ssDNA.
    The idea is to make ssdna vs ssdna profile after that step

    Parameters
    ----------
    sample_sparse_mat : str
        Path to the sparse matrix file (hicstuff given output).

    oligo_capture_with_frag_path : str
        Path to the oligo capture file (sshicstuff mandatory table).
        Must be the file with the fragments associated made with the 'associate' command

    output_path : str
        Path to the output file to be created.
        Default is None.

    force : bool
        Force the overwriting of the oligo file even if the columns are already present.
        Default is True.

    Returns
    -------
    None
    """

    if not output_path:
        output_path = sample_sparse_mat.replace(".txt", "_ssdna_only.txt")

    if not force and os.path.exists(output_path):
        logger.info(
            "[Sparse Matrix Graal (ssdna)] : Output file already exists: %s",
            output_path,
        )
        logger.warning(
            "[Sparse Matrix Graal (ssdna)] : Use the --force / -F flag to overwrite the existing file."
        )
        return

    check_if_exists(sample_sparse_mat)
    check_if_exists(oligo_capture_with_frag_path)
    check_file_extension(sample_sparse_mat, ".txt")
    check_file_extension(oligo_capture_with_frag_path, [".csv", ".tsv"])

    oligo_capture_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep="\t", header=0)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_capture_delim)

    df_contacts_ssdna_only = df_sparse_mat.copy(deep=True)
    ssdna_frag = pd.unique(df_oligo.loc[df_oligo["type"] == "ss", "fragment"]).tolist()

    df_contacts_ssdna_only = df_contacts_ssdna_only[
        df_contacts_ssdna_only.iloc[:, 0].isin(ssdna_frag)
        & df_contacts_ssdna_only.iloc[:, 1].isin(ssdna_frag)
    ]

    df_contacts_ssdna_only.columns = [
        len(ssdna_frag),
        len(ssdna_frag),
        len(df_contacts_ssdna_only) + 1,
    ]
    df_contacts_ssdna_only.reset_index(drop=True, inplace=True)

    df_contacts_ssdna_only.to_csv(output_path, sep="\t", index=False, header=True)
    logger.info(
        "[Sparse Matrix Graal (ssdna)] : ssDNA only contacts saved to %s", output_path
    )


def subsample(
    input_path: str,
    seed: int = 100,
    size: int = 4000000,
    compress: bool = True,
    force: bool = False,
):
    """
    Subsample and compress FASTQ file using seqtk.

    Parameters
    ----------
    input_path : str
        Path to the input FASTQ file.
    seed : int
        Seed for random number generator.
    size : int
        Number of reads to subsample randomly.
    compress : bool
        Whether to compress the output file (with gzip).
    force : bool
        Whether to overwrite the output file if it already exists.
    """

    check_seqtk()

    # Determine the appropriate suffix for output file based on size
    if size >= 1000000000:
        suffix = f"sub{size // 1000000000}G"
    elif size >= 1000000:
        suffix = f"sub{size // 1000000}M"
    elif size >= 1000:
        suffix = f"sub{size // 1000}K"
    else:
        suffix = f"sub{size}"

    # Construct the output path
    pattern = r"(.+)(\.end[12]\.fastq)\.gz"
    match = re.match(pattern, input_path)
    if match:
        prefix = match.group(1)
        end_part = match.group(2)
        output_path = f"{prefix}_{suffix}{end_part}"
    else:
        logger.error("Input file does not match expected pattern")
        raise ValueError("Input file does not match expected pattern")

    if os.path.exists(output_path) and not force:
        logger.warning(
            "Output file %s already exists. skipping subsampling.", output_path
        )
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    # Log the start of subsampling
    logger.info(
        "Starting subsampling of %s to %s.gz with seed %d and size %d",
        input_path,
        output_path,
        seed,
        size,
    )

    # Run seqtk to subsample the FASTQ file
    seqtk_command = f"seqtk sample -s {seed} {input_path} {size} > {output_path}"
    try:
        subprocess.run(seqtk_command, shell=True, check=True)
        logger.info("Subsampling completed successfully for %s", output_path)
    except subprocess.CalledProcessError as e:
        logger.error("Error in subsampling %s: %s", input_path, e)
        raise

    # Optionally compress the output file
    if compress:
        check_gzip()
        if os.path.exists(output_path + ".gz"):
            logger.warning(
                "Output file %s.gz already exists. removing it.", output_path
            )
            os.remove(output_path + ".gz")

        gzip_command = f"gzip {output_path}"
        logger.info("Compressing %s", output_path)
        try:
            subprocess.run(gzip_command, shell=True, check=True)
            output_path += ".gz"
            logger.info("Compression completed successfully for %s", output_path)
        except subprocess.CalledProcessError as e:
            logger.error("Error in compressing %s: %s", output_path, e)
            raise


def transform_data(data: np.array, y_max: float, user_y_max: float, y_min: float, re_scale: bool):
    re_scale_output = ""
    if re_scale:
        if y_max <= 1.:
            # squared root transformation
            new_data = np.sqrt(data + 1e-8)
            y_max = np.sqrt(y_max) if not user_y_max else user_y_max
            y_min = np.sqrt(y_min) if y_min > 0 else 0
            re_scale_output = "sqrt"
        else:
            # log transformation
            new_data = np.log(data + 1)
            y_max = np.log(y_max) if not user_y_max else user_y_max
            y_min = 0
            re_scale_output = "log"
    else:
        new_data = data

    return new_data, y_max, y_min, re_scale_output


def uploaded_files_cache(cache_dir: str):
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(cache_dir):
        path = join(cache_dir, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files
