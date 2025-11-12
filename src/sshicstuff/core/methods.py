"""
Module containing functions to analyze the contacts and the capture efficiency of the oligos.
"""

import argparse
import base64
import datetime
import os
import random as rd
import re
import shutil
import subprocess
import sys
from os.path import join, dirname
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.io as pio
from Bio import SeqIO
from Bio.Seq import Seq

import sshicstuff.log as log

logger = log.logger

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

pio.kaleido.scope.mathjax = None

# cache dir
__CACHE_DIR__ = os.environ.get("SSHIC_CACHE_DIR", "/tmp/sshic_cache")
os.makedirs(__CACHE_DIR__, exist_ok=True)

def associate_oligo_to_frag(
        oligo_capture_path: str,
        fragments_path: str,
        output_path: str = None,
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
    output_path : str
        Path to the output file. If None, it will be created in the same directory as oligo_capture_path.
    Returns
    -------
    None
    """

    logger.info(
        "[Associate] : Associate oligo/probe name to fragment/read ID that contains it."
    )

    check_file_extension(fragments_path, ".txt")
    check_file_extension(oligo_capture_path, [".csv", ".tsv", ".txt"])

    if not output_path:
        output_path: str = oligo_capture_path.replace(".csv", "_fragments_associated.csv")

    os.makedirs(dirname(output_path), exist_ok=True)

    logger.info(
        "[Associate] : Creating a new oligo_capture table : %s",
        output_path.split("/")[-1],
    )

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

def annealing_to_capture(
    df_annealing: pd.DataFrame,
    enzyme: str,
    target_length: float,
):

    df_capture = df_annealing.copy()
    df_capture.drop(columns=["sequence_original", "sequence_modified", "length"], inplace=True)

    capture_oligos = []
    for _, row in df_annealing.iterrows():
        seq = row["sequence_modified"]
        pos = seq.lower().find(enzyme.lower())
        middle = len(seq) // 2
        if pos < middle:
            # cut from 5' end
            n_5_prime_deletion = pos + len(enzyme)
            n_3_prime_deletion = max(0, len(seq) - (n_5_prime_deletion + target_length))
        else:
            # cut from 3' end
            n_3_prime_deletion = len(seq) - pos
            n_5_prime_deletion = max(0, len(seq) - (n_3_prime_deletion + target_length))

        if seq is None or pd.isna(seq):
            new_seq = row["sequence_original"]
        else:
            new_seq = seq[n_5_prime_deletion: -n_3_prime_deletion]
        capture_oligos.append(new_seq)
    df_capture["sequence"] = capture_oligos

    logger.info(f"[Design/Capture] Creation of capture oligos from annealing oligos done and saved.")

    return df_capture



def check_file_extension(file_path: str | Path, extension: str | list[str]):
    """
    Check if a file has the correct extension.
    """

    file_path = str(file_path)  # force string to use endswith()

    if isinstance(extension, list):
        for ext in extension:
            if file_path.endswith(ext):
                return
        logger.error("File %s does not have the correct extension %s.", file_path, extension)
    else:
        if file_path.endswith(extension):
            return
        logger.error("File %s does not have the correct extension %s.", file_path, extension)


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
            logger.info("gzip version %s is installed.", version)
            return version
        else:
            logger.error("Unable to determine gzip version from the output.")
            return None
    except subprocess.CalledProcessError:
        logger.error("gzip is not installed or not functioning correctly. "
                      "Please install or fix gzip before running this function.")
        return None
    except FileNotFoundError as e:
        logger.error("gzip not found: %s", e)
        return None
    except subprocess.SubprocessError as e:
        logger.error("Subprocess error when checking gzip version: %s", e)
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
        logger.error("File %s does not exist.", file_path)
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
            logger.info("seqtk version %s is installed.", version)
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
        "[Compare] : Comparing the capture efficiency of a sample with a wild-type reference."
    )

    logger.info(
        "[Compare] : Be sure to have same number of reads for both samples. Otherwise use subsample sub-module."
    )

    df_sample: pd.DataFrame = pd.read_csv(stats1_path, header=0, sep="\t")
    df_wt: pd.DataFrame = pd.read_csv(stats2_path, sep="\t")

    df_cap_eff = pd.DataFrame(
        columns=[
            "probe",
            "capture_efficiency",
            f"capture_efficiency_{ref_name}",
            f"ratio_sample_vs_{ref_name}",
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
        output_dir, f"{os.path.basename(stats1_path).split('.')[0]}_vs_{ref_name}.tsv"
    )
    df_cap_eff.to_csv(output_path, sep="\t", index=False)

    logger.info(
        "[Compare] : Capture efficiency comparison file saved to %s", output_path
    )


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
        logger.info("[Copy] : %s copied.", src_basename)
    except IOError as e:
        logger.error("Unable to copy file. Error: %s", e)
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
    Compute fragment coverage from a sparse HiC contacts matrix and save a bedgraph file.

    This function computes per-fragment coverage by melting the sparse contacts matrix so that
    each contact is assigned to both fragments. It then aggregates these contacts based on
    fragment coordinates. Optionally, the coverage is binned into fixed-size bins (splitting
    fragments that span two bins proportionally) and normalized by the total number of contacts.

    Parameters
    ----------
    sparse_mat_path : str
        Path to the sparse contacts file.
    fragments_list_path : str
        Path to the fragments list file.
    output_dir : str, optional
        Directory to save the output file. Defaults to the directory of sparse_mat_path.
    normalize : bool, optional
        If True, normalize the coverage by the total number of contacts. Default is False.
    force : bool, optional
        If True, overwrite an existing output file. Default is False.
    bin_size : int, optional
        Size (in bp) of bins for coverage binning. If 0, no binning is performed.
    chromosomes_coord_path : str, optional
        Path to the chromosome sizes file (tsv or csv) required for binning.

    Returns
    -------
    None
    """

    # Set output directory and file name
    if output_dir is None:
        output_dir = os.path.dirname(sparse_mat_path)
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(sparse_mat_path))[0]
    output_path = os.path.join(output_dir, base_name + "_contacts_coverage.bedgraph")

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    # Read fragments list and prepare fragment information
    df_frag = pd.read_csv(fragments_list_path, sep="\t")
    df_frag.rename(
        columns={"chrom": "chr", "start_pos": "start", "end_pos": "end"}, inplace=True
    )
    df_frag["id"] = np.arange(len(df_frag))

    chrom_order = list(df_frag["chr"].unique())

    # Read sparse contacts matrix
    df_contacts = pd.read_csv(
        sparse_mat_path, header=0, sep="\t", names=["frag_a", "frag_b", "contacts"]
    )

    # Melt the contacts DataFrame so each fragment gets its share of the contact
    df_long = pd.concat([
        df_contacts[["frag_a", "contacts"]].rename(columns={"frag_a": "frag"}),
        df_contacts[["frag_b", "contacts"]].rename(columns={"frag_b": "frag"})
    ])

    # Merge with fragment info to retrieve chromosome coordinates
    df_cov = pd.merge(df_long, df_frag, left_on="frag", right_on="id")
    # Aggregate contact counts per fragment location
    df_cov = df_cov.groupby(["chr", "start", "end"], as_index=False)["contacts"].sum()

    if bin_size > 0:
        # --- BINNING SECTION ---
        # Load chromosome sizes and create a DataFrame with empty bins for each chromosome.
        delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
        df_chrom = pd.read_csv(chromosomes_coord_path, sep=delim)
        df_chrom.columns = [c.lower() for c in df_chrom.columns
                            ]
        chrom_sizes = df_chrom.set_index("chr")["length"].to_dict()

        bins_list = []
        for chr_, length in chrom_sizes.items():
            # Create bins from 0 to chromosome length in steps of bin_size
            starts = np.arange(0, length, bin_size)
            bins = pd.DataFrame({
                "chr": chr_,
                "start": starts,
                "end": starts + bin_size,
                "contacts": 0.0,
            })
            bins_list.append(bins)
        df_bins = pd.concat(bins_list, ignore_index=True)

        # Assign each fragment to bins by computing the bin start for its beginning and end.
        df_cov["start_bin"] = (df_cov["start"] // bin_size) * bin_size
        df_cov["end_bin"] = (df_cov["end"] // bin_size) * bin_size

        # Identify fragments that fall entirely within one bin.
        in_bin_mask = df_cov["start_bin"] == df_cov["end_bin"]
        df_in_bin = df_cov.loc[in_bin_mask].copy()
        df_in_bin["bin"] = df_in_bin["start_bin"]

        # For fragments that span two bins, split the contact counts proportionally.
        df_cross = df_cov.loc[~in_bin_mask].copy()
        if not df_cross.empty:
            first_bin_end = df_cross["start_bin"] + bin_size
            frag_length = df_cross["end"] - df_cross["start"]
            # Fraction of the fragment in the first bin
            frac_first = (first_bin_end - df_cross["start"]) / frag_length
            frac_second = 1 - frac_first

            df_cross_first = df_cross.copy()
            df_cross_first["bin"] = df_cross_first["start_bin"]
            df_cross_first["contacts"] *= frac_first

            df_cross_second = df_cross.copy()
            df_cross_second["bin"] = df_cross_second["end_bin"]
            df_cross_second["contacts"] *= frac_second

            df_cross_combined = pd.concat([df_cross_first, df_cross_second], ignore_index=True)
        else:
            df_cross_combined = pd.DataFrame(columns=df_cov.columns.tolist() + ["bin"])

        # Combine in-bin and split entries, then group by chromosome and bin start.
        df_binned = pd.concat([df_in_bin, df_cross_combined], ignore_index=True)
        df_binned = df_binned.groupby(["chr", "bin"], as_index=False)["contacts"].sum()
        df_binned["start"] = df_binned["bin"]
        df_binned["end"] = df_binned["bin"] + bin_size
        df_binned = df_binned[["chr", "start", "end", "contacts"]]

        # Merge with the complete bins to ensure no bin is missing
        df_final = pd.concat([df_bins, df_binned], ignore_index=True)
        df_final = df_final.groupby(["chr", "start", "end"], as_index=False)["contacts"].sum()

        df_final['chr'] = pd.Categorical(df_final['chr'], categories=chrom_order, ordered=True)
        df_final = df_final.sort_values(['chr', 'start']).reset_index(drop=True)
        df_final["contacts"] = df_final["contacts"].fillna(0).round(4)

        # Update output file name to include bin size suffix
        bin_suffix = get_bin_suffix(bin_size)
        output_path = output_path.replace(".bedgraph", f"_{bin_suffix}.bedgraph")
        df_final.to_csv(output_path, sep="\t", index=False, header=False)
        logger.info("[Coverage] Binned contacts coverage file saved to %s", output_path)
        result_df = df_final
    else:
        # No binning: save fragment-level coverage directly
        df_cov['chr'] = pd.Categorical(df_cov['chr'], categories=chrom_order, ordered=True)
        df_cov = df_cov.sort_values(['chr', 'start']).reset_index(drop=True)
        df_cov.to_csv(output_path, sep="\t", index=False, header=False)
        logger.info("[Coverage] Contacts coverage file saved to %s", output_path)
        result_df = df_cov

    if normalize:
        # Normalize contacts by the total sum
        total_contacts = result_df["contacts"].sum()
        if total_contacts:
            result_df["contacts"] /= total_contacts
        norm_output_path = output_path.replace("_contacts_", "_frequencies_")
        result_df.to_csv(norm_output_path, sep="\t", index=False, header=False)
        logger.info("[Coverage] Normalized coverage file saved to %s", norm_output_path)

    logger.info("[Coverage] Coverage calculation completed.")


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
    df_annealing: pd.DataFrame,
    genome_input: str,
    output_dir: str,
    enzyme: str,
    n_artificial_spacer: int = 150,
):
    # ── Stage 0: Inputs / setup
    logger.info("[Design/EditGenome] Start. Enzyme=%s | spacer N=%d | outdir=%s",
                enzyme, n_artificial_spacer, output_dir)

    # ── Stage 1: Collect sequences
    n_ss = df_annealing.query("type == 'ss'").shape[0]
    n_ds = df_annealing.query("type == 'ds'").shape[0]
    logger.info("[Design/EditGenome] Loaded oligos: ss=%d | ds=%d", n_ss, n_ds)

    # ── Stage 2: Build artificial chromosomes (strings)
    logger.info("[Design/EditGenome] Building artificial chromosomes (ssDNA & dsDNA) ...")

    chr_arti_dsdna_path = os.path.join(output_dir, "chr_artificial_dsDNA.fa")
    chr_arti_ssdna_path = os.path.join(output_dir, "chr_artificial_ssDNA.fa")

    enzyme = enzyme.upper()
    dsdna_seq_series = df_annealing[df_annealing["type"] == "ss"]["sequence_original"]
    dsdna_seq = [seq.upper() for seq in dsdna_seq_series.values]

    ssdna_seq_series = df_annealing[df_annealing["type"] == "ss"]["sequence_modified"]
    ssdna_seq = [seq.upper() for seq in ssdna_seq_series.values]

    s = n_artificial_spacer
    oneline_ssdna = "N" * s + enzyme + "N" * s
    oneline_dsdna = "N" * s + enzyme + "N" * s

    for ss, ds in zip(ssdna_seq, dsdna_seq):
        middle = len(ss) // 2
        enzyme_pos = ss.find(enzyme)
        if enzyme_pos < middle:
            ss2 = ss[enzyme_pos + len(enzyme) :].upper()
            ds2 = ds[enzyme_pos + len(enzyme) :].upper()
        else:
            ss2 = ss[:enzyme_pos].upper()
            ds2 = ds[:enzyme_pos].upper()

        oneline_ssdna += ss2 + "N" * s + enzyme + "N" * s
        oneline_dsdna += ds2 + "N" * s + enzyme + "N" * s

    logger.info("[Design/EditGenome] Built ssDNA (snp) and dsDNA (no snp) artificial, lengths = %d bp", len(oneline_ssdna))

    record_ssdna = SeqIO.SeqRecord(seq=Seq(oneline_ssdna), id="chr_artificial_ssDNA", description=f"({len(oneline_ssdna)} bp)")
    record_dsdna = SeqIO.SeqRecord(seq=Seq(oneline_dsdna), id="chr_artificial_dsDNA", description=f"({len(oneline_dsdna)} bp)")
    SeqIO.write(record_ssdna, chr_arti_ssdna_path, "fasta")
    SeqIO.write(record_dsdna, chr_arti_dsdna_path, "fasta")

    # ── Stage 3: Write the two artificial FASTA files
    logger.info("[Design/EditGenome] Writing artificial chromosomes to: %s, %s",
                chr_arti_ssdna_path, chr_arti_dsdna_path)

    # ── Stage 4: Mask native genome
    logger.info("[Design/EditGenome] Masking original oligo regions with N in %s ...", os.path.basename(genome_input))

    records = list(SeqIO.parse(genome_input, "fasta"))
    genome_name = os.path.basename(genome_input).split(".")[0]

    for _, row in df_annealing.iterrows():
        chr_    = row["chr"]
        seq_original = row["sequence_original"]
        start_ = row["start"]
        end_   = row["end"]

        for record in records:
            if record.id == chr_:
                seq = record.seq
                new_seq = seq[:start_] + "N" * len(seq_original) + seq[end_:]
                record.seq = new_seq

    # ── Stage 5: Append artificial chromosomes and write the edited genome
    genome_output_path = os.path.join(output_dir, f"{genome_name}_edited.fa")
    logger.info("[Design/EditGenome] Appending artificial chromosomes and writing %s", genome_output_path)

    records.append(record_dsdna)
    records.append(record_ssdna)
    SeqIO.write(records, genome_output_path, "fasta")


    pattern = re.compile(rf"N{{5,}}{enzyme}N{{5,}}", re.IGNORECASE)
    matches = list(pattern.finditer(oneline_ssdna))
    enzyme_positions = [m.start() + m.group().lower().find(enzyme.lower()) for m in matches]

    # ── Stage 6: Fragment map on chr_artificial_ssDNA
    logger.info("[Design/EditGenome] Computing fragment coordinates on chr_artificial_ssDNA (by enzyme sites)")
    logger.info("[Design/EditGenome] Found %d enzyme landmarks → %d fragments",
                len(enzyme_positions), max(0, len(enzyme_positions) - 1))

    # Now extract fragment coordinates based on the region between enzymes
    chr_arti_starts = []
    chr_arti_ends = []

    for i in range(len(enzyme_positions) - 1):
        start = enzyme_positions[i]
        end = enzyme_positions[i + 1]
        chr_arti_starts.append(start)
        chr_arti_ends.append(end)

    lengths = [end - start for start, end in zip(chr_arti_starts, chr_arti_ends)]

    df2 = df_annealing[df_annealing["type"] == "ss"].copy()
    df2.rename(columns={
        "chr": "chr_ori",
        "start": "start_ori",
        "end": "end_ori",
    }, inplace=True)

    df2["chr"]      = "chr_artificial_ssDNA"
    df2["start"]    = chr_arti_starts
    df2["end"]      = chr_arti_ends
    df2["length"]   = lengths

    order_col = [
        "chr",
        "start",
        "end",
        "length",
        "chr_ori",
        "start_ori",
        "end_ori",
        "orientation",
        "type",
        "name",
        "sequence_original",
        "sequence_modified",
    ]
    df2 = df2[order_col]

    # Add back ds control (not on artificial chromosome)
    df_dsdna = df_annealing[df_annealing["type"] == "ds"].copy()
    df_annealing_final = pd.concat([df2, df_dsdna], ignore_index=True)

    # ── Stage 7: Return dataframe
    logger.info("[Design/EditGenome] Final table: %d rows (ss artificial + ds controls)", len(df_annealing_final))
    logger.info("[Design/EditGenome] Done.")

    return df_annealing_final



def format_annealing_oligo_output(
        design_output_raw_path: str,
        design_output_snp_path: str,
):
    logger.info("[Design/Annealing] Formatting annealing oligo output from FASTA to CSV ...")
    df_raw = pd.read_csv(design_output_raw_path, sep="\t", header=None)
    df_snp = pd.read_csv(design_output_snp_path, sep="\t", header=None)

    chroms, starts, ends, lenghts, strands, types, names, raw_seqs, snp_seqs = [], [], [], [], [], [], [], [], []

    def parse_metadata(line: str):
        # Format: >chr4:848192:848271_mt0_+_raw
        parts = line[1:].split('_')
        coord_part = parts[0]
        strand = "w" if parts[2] == "+" else "c"  # Convert strand to W/C
        oligo_type = "ss"
        chrom, start, end = coord_part.split(':')

        name = "Probe_{}_{}_{}_{}".format(chrom, strand, start, end)
        return chrom, int(start), int(end), strand, oligo_type, name

    def search_snps(raw, mutated):
        res = list(mutated.upper())
        for p in range(len(raw)):
            if raw[p].upper() != mutated[p].upper():
                res[p] = mutated[p].lower()
        return "".join(res)


    for i in range(0, len(df_raw), 2):
        metadata_raw = df_raw.iloc[i, 0]
        seq_raw = df_raw.iloc[i + 1, 0]
        seq_snp = df_snp.iloc[i + 1, 0]

        seq_snp = search_snps(seq_raw, seq_snp)

        chrom, start, end, strand, oligo_type, name = parse_metadata(metadata_raw)

        chroms.append(chrom)
        starts.append(start)
        ends.append(end)
        lenghts.append(end - start + 1)  # Length of the oligo
        strands.append(strand)  # Convert strand to W/C
        types.append(oligo_type)
        names.append(name)
        raw_seqs.append(seq_raw)
        snp_seqs.append(seq_snp)

    df_final = pd.DataFrame({
        "chr": chroms,
        "start": starts,
        "end": ends,
        "length": lenghts,
        "orientation": strands,
        "type": types,
        "name": names,
        "sequence_original": raw_seqs,
        "sequence_modified": snp_seqs
    })

    return df_final

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
    """
    Aggregate probes into groups and add new columns to the DataFrame.

    For each group defined in df_groups, this function maps the comma-separated list of probes
    to their associated fragment identifiers (converted to strings) using prob2frag, and then computes
    the row-wise average or sum of the corresponding columns in df. If some fragments are missing,
    only the available ones are used. If none are found, the new group column is filled with NaN.

    Parameters
    ----------
    df_groups : pd.DataFrame
        DataFrame with group definitions. Must have columns "probes", "name", and "action".
    df : pd.DataFrame
        DataFrame containing contact data, where columns corresponding to probes (or fragments) are aggregated.
    prob2frag : dict
        Mapping from probe names to fragment identifiers (which may be int or str).

    Returns
    -------
    None
    """
    # Build a mapping from stringified column names to the actual column names in df
    col_map = {str(col): col for col in df.columns}

    for row in df_groups.itertuples(index=False):
        # Split the comma-separated probes and map each probe to its fragment identifier as string.
        group_probes = re.findall(r"[A-Za-z0-9_-]+", row.probes)
        group_frags = np.unique([str(prob2frag.get(probe, probe)) for probe in group_probes])
        group_name = "$" + row.name.lower()

        # Get only the columns that exist in df (using the stringified mapping)
        existing_frags = [col_map[frag] for frag in group_frags if frag in col_map]
        if not existing_frags:
            logger.warning("Group %s: none of the fragments %s are present in the DataFrame.",
                           group_name, group_frags)
            df[group_name] = np.nan
        else:
            if row.action.lower() == "average":
                df[group_name] = df[existing_frags].mean(axis=1)
            elif row.action.lower() == "sum":
                df[group_name] = df[existing_frags].sum(axis=1)


def merge_sparse_mat(
    output_path: str = None, force: bool = False, matrices: list[str] = None
) -> None:
    """
    Merge multiple sparse matrices into one.
    """

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

def namespace_to_args(ns: argparse.Namespace, flag_map: dict[str, str]) -> list[str]:
    """
    Convert an argparse namespace to a flat CLI list according to flag_map.
    Skips None; booleans become flags (present -> include flag).
    """
    args = []
    for field, flag in flag_map.items():
        val = getattr(ns, field, None)
        if isinstance(val, bool):
            if val:
                args.append(flag)
        elif val is None:
            continue
        else:
            args.extend([flag, str(val)])
    return args


def resolve_outpath(outdir: Path, name_or_path: str | None, default_name: str) -> Path:
    """
    If name_or_path is:
      - None           -> outdir/default_name
      - absolute path  -> as is
      - relative name  -> outdir/name_or_path
    """
    if name_or_path is None:
        return outdir / default_name
    p = Path(name_or_path)
    return p if p.is_absolute() else (outdir / p)


def save_file_cache(name, content, cache_dir):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(join(cache_dir, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def sparse_with_dsdna_only(
    sample_sparse_mat: str,
    oligo_capture_with_frag_path: str,
    n_flanking_dsdna: int = 2,
    output_path: str = None,
    force: bool = False,
) -> None:
    """
    Generate a sparse matrix file (with the same format as the input) that excludes single-stranded DNA (ssDNA)
    fragments. It removes fragments with a probe and the specified number of flanking fragments (both upstream
    and downstream) around each probe fragment.

    Parameters:
        sample_sparse_mat (str): Path to the input sparse matrix file (output from hicstuff).
        oligo_capture_with_frag_path (str): Path to the oligo capture file (required for sshicstuff),
            containing fragment associations from the 'associate' command.
        n_flanking_dsdna (int): Number of flanking fragments to remove around each probe fragment. Default is 2.
        output_path (str): File path to save the output. If None, the output filename is derived from sample_sparse_mat.
        force (bool): If False, the function will not overwrite an existing output file. Default is False.

    Returns:
        None
    """

    # Determine output path if not provided
    if not output_path:
        output_path = sample_sparse_mat.replace(".txt", "_dsdna_only.txt")

    # Do not overwrite existing file unless force is True
    if not force and os.path.exists(output_path):
        logger.info("[Sparse Matrix Graal (dsdna)] Output file already exists: %s", output_path)
        logger.warning("[Sparse Matrix Graal (dsdna)] Use the --force / -F flag to overwrite the existing file.")
        return

    # Check that input files exist and have correct file extensions
    check_if_exists(sample_sparse_mat)
    check_if_exists(oligo_capture_with_frag_path)
    check_file_extension(sample_sparse_mat, ".txt")
    check_file_extension(oligo_capture_with_frag_path, [".csv", ".tsv"])

    # Set delimiter based on file extension (.csv -> comma, .tsv -> tab)
    oligo_capture_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"

    # Read the input sparse matrix and the oligo capture table
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep="\t", header=None)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_capture_delim)

    # Create a working copy of the sparse matrix
    df_contacts_dsdna_only = df_sparse_mat.copy()

    # Extract ssDNA fragments
    ssdna_frag = df_oligo.loc[df_oligo["type"] == "ss", "fragment"].tolist()

    # Extract dsDNA fragments
    dsdna_frag = df_oligo.loc[df_oligo["type"] == "ds", "fragment"].tolist()

    # Compute flanking fragments (both upstream and downstream) using list comprehensions
    flanking_fragments = (
        [f + i for f in dsdna_frag for i in range(1, n_flanking_dsdna + 1)] +
        [f - i for f in dsdna_frag for i in range(1, n_flanking_dsdna + 1)]
    )
    # Get unique dsDNA fragments including their flanking fragments
    dsdna_frag_all = np.unique(dsdna_frag + flanking_fragments)

    # Build a DataFrame similar to the original implementation for header count adjustment
    df_ssdna = pd.DataFrame(ssdna_frag, columns=["fragments"])
    df_dsdna = pd.DataFrame(dsdna_frag_all, columns=["fragments"])
    df_frag = pd.concat([df_ssdna, df_dsdna])
    
    # For fast membership testing, create a set of fragments to remove
    fragments_to_remove = set(df_frag["fragments"])

    # Identify rows where either column 0 or 1 contains a fragment to remove
    mask = df_contacts_dsdna_only[0].isin(fragments_to_remove) | df_contacts_dsdna_only[1].isin(fragments_to_remove)
    index_to_drop = df_contacts_dsdna_only.index[mask]

    # Drop the identified rows from the matrix
    df_contacts_dsdna_only.drop(index_to_drop, inplace=True)

    # Adjust header counts (first row): subtract total number of fragments removed and dropped contacts
    df_contacts_dsdna_only.iloc[0, 0] -= len(df_frag)
    df_contacts_dsdna_only.iloc[0, 1] -= len(df_frag)
    df_contacts_dsdna_only.iloc[0, 2] -= len(index_to_drop)

    # Save the resulting sparse matrix to the output file
    df_contacts_dsdna_only.to_csv(output_path, sep="\t", index=False, header=False)
    logger.info("[Sparse Matrix Graal (dsdna)] dsDNA only contacts saved to %s", output_path)


def sparse_with_ssdna_only(
    sample_sparse_mat: str,
    oligo_capture_with_frag_path: str,
    output_path: str = None,
    force: bool = False,
) -> None:
    """
    Generate a sparse matrix file (with the same format as the input) containing only ssDNA contacts.
    This matrix is used for ssDNA vs ssDNA profiling.

    Parameters:
        sample_sparse_mat (str): Path to the input sparse matrix file (output from hicstuff).
        oligo_capture_with_frag_path (str): Path to the oligo capture file (required for sshicstuff),
            containing fragment associations from the 'associate' command.
        output_path (str): File path to save the output. If None, the filename is derived from sample_sparse_mat.
        force (bool): If False, the function will not overwrite an existing output file. Default is False.

    Returns:
        None
    """

    # Set the output file name if not provided
    if not output_path:
        output_path = sample_sparse_mat.replace(".txt", "_ssdna_only.txt")

    # Do not overwrite an existing file if force is not enabled
    if not force and os.path.exists(output_path):
        logger.info("[Sparse Matrix Graal (ssdna)] Output file already exists: %s", output_path)
        logger.warning("[Sparse Matrix Graal (ssdna)] Use the --force / -F flag to overwrite the existing file.")
        return

    # Check input file existence and file extension
    check_if_exists(sample_sparse_mat)
    check_if_exists(oligo_capture_with_frag_path)
    check_file_extension(sample_sparse_mat, ".txt")
    check_file_extension(oligo_capture_with_frag_path, [".csv", ".tsv"])

    # Set delimiter based on the oligo capture file extension
    oligo_capture_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"

    # Read the input sparse matrix and oligo capture files
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep="\t", header=0)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_capture_delim)

    # Extract unique ssDNA fragments from the oligo capture file
    ssdna_frag = pd.unique(df_oligo.loc[df_oligo["type"] == "ss", "fragment"]).tolist()

    # Filter the sparse matrix: keep rows where both fragment columns are ssDNA fragments
    mask = df_sparse_mat.iloc[:, 0].isin(ssdna_frag) & df_sparse_mat.iloc[:, 1].isin(ssdna_frag)
    df_contacts_ssdna_only = df_sparse_mat[mask].copy()

    # Update the header counts:
    # The first two columns become the count of ssDNA fragments,
    # and the third column is set to the number of rows in the filtered matrix plus one.
    df_contacts_ssdna_only.columns = [
        len(ssdna_frag),
        len(ssdna_frag),
        len(df_contacts_ssdna_only) + 1,
    ]
    df_contacts_ssdna_only.reset_index(drop=True, inplace=True)

    # Save the resulting matrix to the specified output file with headers
    df_contacts_ssdna_only.to_csv(output_path, sep="\t", index=False, header=True)
    logger.info("[Sparse Matrix Graal (ssdna)] ssDNA only contacts saved to %s", output_path)


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
    """
    Transform the data using a log or square root transformation if necessary.
    """
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
