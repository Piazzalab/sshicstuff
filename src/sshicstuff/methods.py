import os
import re
import subprocess
import datetime

import numpy as np
import pandas as pd

import plotly.graph_objs as go
import plotly.io as pio
from plotly.subplots import make_subplots

import sshicstuff.utils as utils
import sshicstuff.log as log
import sshicstuff.gui.graph as graph
import sshicstuff.colors as colors


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
        arm_length_classification: bool = False
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

    output_dir = os.path.join(output_dir, 'aggregated')
    output_dir = os.path.join(output_dir, 'telomeres') if telomeres else os.path.join(output_dir, 'centromeres')

    os.makedirs(output_dir, exist_ok=True)
    sample_name = re.match(r".+\/(.+)_profile", binned_contacts_path).group(1)
    sample_short_name = sample_name.split("_")[0]

    output_prefix = os.path.join(output_dir, sample_short_name) + f"_agg_on_"
    output_prefix += "telo" if telomeres else "cen"

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    coords_delim = "\t" if chr_coord_path.endswith(".tsv") else ","
    df_coords: pd.DataFrame = pd.read_csv(chr_coord_path, sep=coords_delim, index_col=None)
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    df_contacts: pd.DataFrame = pd.read_csv(binned_contacts_path, sep='\t')

    binsize = int(df_contacts.loc[2, 'chr_bins'] - df_contacts.loc[1, 'chr_bins'])
    logger.info(f"[Aggregate] : Contacts binned profile with resolution of : {binsize} bp")

    chr_list = list(df_coords['chr'].unique())
    fragments = df_oligo["fragment"].astype(str).tolist()
    groups = [g for g in df_contacts.columns if g.startswith("$")]

    if len(excluded_chr_list) > 0:
        logger.info(f"[Aggregate] : Excluding chromosomes:  {', '.join(excluded_chr_list)}")
        df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr_list)]
        df_coords = df_coords[~df_coords['chr'].isin(excluded_chr_list)]

    if inter_only:
        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        logger.info("[Aggregate] : Excluding intra-chr contacts")
        for frag in fragments:
            ii_frag = df_oligo.loc[df_oligo["fragment"] == int(frag)].index[0]
            probe_chr_ori = df_oligo.loc[ii_frag, 'chr_ori']
            if probe_chr_ori not in excluded_chr_list:
                df_contacts.loc[df_contacts['chr'] == probe_chr_ori, frag] = np.nan

        output_prefix += "_inter"

    if normalize:
        logger.info("[Aggregate] : Normalizing the contacts")
        df_contacts.loc[:, fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
        output_prefix += "_norm"

    if centromeres:
        logger.info(f"[Aggregate] : Aggregating contacts around centromeres")
        logger.info(f"[Aggregate] : Window size: {window_size} bp on each side of the centromere")

        df_merged: pd.DataFrame = pd.merge(df_contacts, df_coords, on='chr')
        df_merged_cen_areas: pd.DataFrame = df_merged[
            (df_merged.chr_bins > (df_merged.left_arm_length - window_size - binsize)) &
            (df_merged.chr_bins < (df_merged.left_arm_length + window_size))]
        df_merged_cen_areas['chr_bins'] = \
            abs(df_merged_cen_areas['chr_bins'] - (df_merged_cen_areas['left_arm_length'] // binsize) * binsize)
        df_grouped: pd.DataFrame = df_merged_cen_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(
            numeric_only=True)
        df_grouped.drop(columns=['length', 'left_arm_length', 'right_arm_length', 'genome_bins'], axis=1, inplace=True)

    elif telomeres:
        df_telos: pd.DataFrame = pd.DataFrame({'chr': df_coords['chr'], 'telo_l': 0, 'telo_r': df_coords['length']})
        df_merged: pd.DataFrame = pd.merge(df_contacts, df_telos, on='chr')
        df_merged_telos_areas_part_a: pd.DataFrame = \
            df_merged[df_merged.chr_bins < (df_merged.telo_l + window_size + binsize)]
        df_merged_telos_areas_part_b: pd.DataFrame = \
            df_merged[df_merged.chr_bins > (df_merged.telo_r - window_size - binsize)]
        df_merged_telos_areas_part_b['chr_bins'] = \
            abs(df_merged_telos_areas_part_b['chr_bins'] - (
                        df_merged_telos_areas_part_b['telo_r'] // binsize) * binsize)
        df_merged_telos_areas: pd.DataFrame = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
        df_grouped: pd.DataFrame = df_merged_telos_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(
            numeric_only=True)
        df_grouped.drop(columns=['telo_l', 'telo_r', 'genome_bins'], axis=1, inplace=True)
        del df_merged, df_merged_telos_areas_part_a, df_merged_telos_areas_part_b

        if arm_length_classification:
            if "category" not in df_coords.columns:
                logger.error(
                    "[Aggregate] :"
                    "The 'category' column is missing in the centromeres file. "
                    "Must be in the form small_small or long_middle concerning lengths of left_right arms")
            else:
                logger.info("[Aggregate] : Classifying the contacts by chromosome arm lengths")

                df_arms_size: pd.DataFrame = pd.DataFrame(columns=["chr", "arm", "size", "category"])
                for _, row in df_coords.iterrows():
                    chr_ = row["chr"]
                    if chr_ not in excluded_chr_list:
                        left_, right_, category_ = row["left_arm_length"], row["right_arm_length"], row["category"]
                        if pd.isna(left_) or pd.isna(right_) or pd.isna(category_):
                            continue
                        df_arms_size.loc[len(df_arms_size)] = chr_, "left", left_, category_.split("_")[0]
                        df_arms_size.loc[len(df_arms_size)] = chr_, "right", right_, category_.split("_")[1]

                df_merged2 = pd.merge(df_contacts, df_telos, on='chr')
                df_merged_telos_areas_part_a = df_merged2[df_merged2.chr_bins < (df_merged2.telo_l + 3000 + 1000)]
                df_merged_telos_areas_part_a.insert(2, 'arm', 'left')
                df_merged_telos_areas_part_b = df_merged2[df_merged2.chr_bins > (df_merged2.telo_r - 3000 - 1000)]
                df_merged_telos_areas_part_b.insert(2, 'arm', 'right')

                df_telo_freq = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
                df_merged3 = pd.merge(df_telo_freq, df_arms_size, on=['chr', 'arm'])
                df_merged3.drop(columns=['telo_l', 'telo_r', 'size'], inplace=True)

                df_grouped_by_arm = df_merged3.groupby(by='category', as_index=False).mean(numeric_only=True)
                df_grouped_by_arm.drop(columns=['chr_bins', 'genome_bins'], inplace=True)
                df_grouped_by_arm = df_grouped_by_arm.rename(columns={'category': 'fragments'}).T
                df_grouped_by_arm.to_csv(output_prefix + "_by_arm_sizes.tsv", sep='\t', header=False)

    else:
        return

    df_grouped = utils.sort_by_chr(df_grouped, chr_list, 'chr', 'chr_bins')
    df_grouped['chr_bins'] = df_grouped['chr_bins'].astype('int64')

    logger.info(f"[Aggregate] : Compute mean, median, std on the aggregated contacts per probe or group of probes, per chromosome")
    df_aggregated_mean: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).mean(numeric_only=True)
    df_aggregated_mean.to_csv(output_prefix + "_mean.tsv", sep="\t")
    df_aggregated_std: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).std(numeric_only=True)
    df_aggregated_std.to_csv(output_prefix + "_std.tsv", sep="\t")
    df_aggregated_median: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).median(numeric_only=True)
    df_aggregated_median.to_csv(output_prefix + "_median.tsv", sep="\t")

    for col in fragments + groups:
        if col in fragments:
            name = col
        else:
            name = col[1:]

        if df_grouped[col].sum() == 0:
            continue

        df_chr_centros_pivot = df_grouped.pivot_table(index='chr_bins', columns='chr', values=col, fill_value=0)
        df_chr_centros_pivot.to_csv(output_prefix + f"_{name}_per_chr.tsv", sep='\t')


def associate_oligo_to_frag(
        oligo_capture_path: str,
        fragments_path: str,
        force: bool = True
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

    logger.info("[Associate] : Associate oligo/probe name to fragment/read ID that contains it.")

    utils.check_file_extension(fragments_path, ".txt")
    utils.check_file_extension(oligo_capture_path, [".csv", ".tsv", ".txt"])

    output_path: str = oligo_capture_path.replace(".csv", "_fragments_associated.csv")
    logger.info(f"[Associate] : Creating a new oligo_capture table : {output_path.split('/')[-1]}")

    if os.path.exists(output_path) and not force:
        logger.info(f"[Associate] : Output file already exists: {output_path}")
        logger.info("[Associate] : Use the --force / -F flag to overwrite the existing file.")
        return

    # Read the oligo and fragments files
    oligo_delim = "," if oligo_capture_path.endswith(".csv") else "\t"
    df_oligo = pd.read_csv(oligo_capture_path, sep=oligo_delim)

    df_fragments = pd.read_csv(fragments_path, sep='\t')
    df_fragments['frag'] = [k for k in range(len(df_fragments))]

    fragments_id = []
    fragments_start = []
    fragments_end = []
    for index, row in df_oligo.iterrows():
        (chr_, probe_start, probe_end, probe_chr_ori, probe_start_ori,
         probe_end_ori, probe_type, probe, probe_seq) = row[:9]
        df_sub_fragments = df_fragments[df_fragments['chrom'] == chr_]
        df_sub_fragment_sorted_start = np.sort(df_sub_fragments['start_pos'].to_numpy())

        probe_middle = int(probe_start + (probe_end - probe_start) / 2)

        idx = np.searchsorted(df_sub_fragment_sorted_start, probe_middle, side="left")
        nearest_frag_start = df_sub_fragment_sorted_start[idx - 1]

        frag_id = df_sub_fragments.index[df_sub_fragments['start_pos'] == nearest_frag_start].tolist()[0]
        frag_start = df_sub_fragments.loc[frag_id, 'start_pos']
        frag_end = df_sub_fragments.loc[frag_id, 'end_pos']
        fragments_id.append(frag_id)
        fragments_start.append(frag_start)
        fragments_end.append(frag_end)

    df_oligo['fragment'] = fragments_id
    df_oligo['fragment_start'] = fragments_start
    df_oligo['fragment_end'] = fragments_end
    df_oligo.to_csv(output_path, sep=",", index=False)

    logger.info("[Associate] : oligos associated to fragments successfully.")


def compare_with_wt(
        stats1_path: str,
        stats2_path: str,
        ref_name: str,
        output_dir: str = None
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

    logger.info("Comparing the capture efficiency of a sample with a wild-type reference.")
    logger.info("Be sure to have same number of reads for both samples. Otherwise use subsample function.")

    df_sample: pd.DataFrame = pd.read_csv(stats1_path, header=0, sep="\t")
    df_wt: pd.DataFrame = pd.read_csv(stats2_path, sep='\t')

    df_cap_eff = pd.DataFrame(columns=[
        "probe",
        "capture_efficiency",
        f"capture_efficiency_{ref_name}",
        f"ratio_sample_vs_wt"
    ])

    for index, row in df_sample.iterrows():
        probe = row['probe']
        cap_eff = row['dsdna_norm_capture_efficiency']

        if probe in df_wt['probe'].tolist():
            cap_eff_wt = df_wt.loc[df_wt['probe'] == probe, "dsdna_norm_capture_efficiency"].tolist()[0]
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

    output_path = os.path.join(output_dir, f"{os.path.basename(stats1_path).split('.')[0]}_vs_{ref_name}.csv")
    df_cap_eff.to_csv(output_path, sep='\t', index=False)


def coverage(
        sparse_mat_path: str,
        fragments_list_path: str,
        output_dir: str = None,
        normalize: bool = False,
        force: bool = False,
        bin_size: int = 0
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

    Returns
    -------
    None
    """

    logger.info("[Coverage] : Calculating coverage per fragment into a bedgraph.")

    if output_dir is None:
        output_dir = os.path.dirname(sparse_mat_path)

    output_path = os.path.join(output_dir, os.path.basename(sparse_mat_path).split('.')[0])

    if os.path.exists(output_path) and not force:
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    df_fragments: pd.DataFrame = pd.read_csv(fragments_list_path, sep='\t')
    df_fragments.rename(columns={'chrom': 'chr', 'start_pos': 'start', 'end_pos': 'end'}, inplace=True)
    df_fragments['id'] = list(range(len(df_fragments)))

    df_hic_contacts: pd.DataFrame = pd.read_csv(
        sparse_mat_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])

    df_coverage: pd.DataFrame = df_fragments[['chr', 'start', 'end']]
    df_coverage['contacts'] = np.nan

    df_merged_a: pd.DataFrame = df_hic_contacts.merge(
        df_fragments[['id', 'chr', 'start', 'end']],
        left_on='frag_a',
        right_on='id',
        suffixes=('', '_a')).drop(columns=['frag_a', 'frag_b'])

    df_merged_b: pd.DataFrame = df_hic_contacts.merge(
        df_fragments[['id', 'chr', 'start', 'end']],
        left_on='frag_b',
        right_on='id',
        suffixes=('', '_b')).drop(columns=['frag_a', 'frag_b'])

    df_grouped_a: pd.DataFrame = df_merged_a.groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()
    df_grouped_b: pd.DataFrame = df_merged_b.groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()

    df_contacts_cov: pd.DataFrame = pd.concat(
        (df_grouped_a, df_grouped_b)).groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()

    df_contacts_cov.index = df_contacts_cov.id
    df_contacts_cov.drop(columns=['id'], inplace=True)
    output_path = output_path + "_contacts_coverage.bedgraph"

    if bin_size > 0:
        bin_suffix = str(bin_size // 1000) + "kb"
        output_path = output_path.replace(".bedgraph", f"_{bin_suffix}.bedgraph")
        logger.info(f"[Coverage] : Binning the bedgraph at {bin_suffix} resolution.")
        df_contacts_cov_bin: pd.DataFrame = df_contacts_cov.copy(deep=True)
        df_contacts_cov_bin['start'] = df_contacts_cov_bin['start'] // bin_size * bin_size
        df_contacts_cov_bin['end'] = df_contacts_cov_bin['end'] // bin_size * bin_size
        df_contacts_cov_bin = df_contacts_cov_bin.groupby(['chr', 'start', 'end'], as_index=False).sum()
        df_contacts_cov_bin.to_csv(output_path, sep='\t', index=False, header=False)
        logger.info(f"[Coverage] : Contacts coverage binned file saved to {output_path}")
        df_contacts_cov = df_contacts_cov_bin

    else:
        df_contacts_cov.to_csv(output_path, sep='\t', index=False, header=False)
        logger.info(f"[Coverage] : Contacts coverage file saved to {output_path}")

    if normalize:
        output_path = output_path.replace("_contacts_", "_frequencies_")
        logger.info("[Coverage] : Normalizing coverage by the total number of contacts.")
        df_frequencies_cov: pd.DataFrame = df_contacts_cov.copy(deep=True)
        df_frequencies_cov["contacts"] /= sum(df_frequencies_cov["contacts"])
        df_frequencies_cov.rename(columns={"contacts": "frequencies"})
        df_frequencies_cov.to_csv(output_path, sep='\t', index=False, header=False)

    logger.info("[Coverage] : Coverage calculation completed.")

def edit_genome_ref(
        annealing_input: str,
        genome_input: str,
        enzyme: str,
        fragment_size: int = 150,
        fasta_spacer: str = "N",
        fasta_line_length: int = 80,
        additional_fasta_path: str = None
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
    logger.info(f"Creating the artificial chromosome with the annealing oligo and the enzyme {enzyme}")

    df = pd.read_csv(annealing_input, sep=',')
    ssdna_seq_series = df[df["type"] == "ss"]['sequence_modified']
    ssdna_seq = [seq.lower() for seq in ssdna_seq_series.values]

    lg = fasta_line_length
    s = fragment_size - len(enzyme)
    p = fasta_spacer
    oneline = p * int(s/2) + enzyme + p * s

    for seq in ssdna_seq:
        middle = len(seq) // 2
        enzyme_pos = seq.find(enzyme)
        if enzyme_pos < middle:
            seq2 = seq[enzyme_pos+len(enzyme):].upper()
        else:
            seq2 = seq[:enzyme_pos].upper()

        oneline += seq2 + p * s + enzyme + p * s

    lines = "\n".join([oneline[i:i+lg] for i in range(0, len(oneline), lg)])
    fasta = f">chr_artificial_ssDNA\t ({len(oneline)} bp)\n{lines}"

    with open(artificial_chr_path, "w") as f:
        f.write(fasta)

    # Inserting the artificial chromosome at the end of the genome .FASTA file
    genome_name = os.path.basename(genome_input)
    logger.info(f"Inserting the artificial chromosome at the end of the original genome .FASTA file")
    with open(genome_input, "r") as f:
        genome = f.read()

    new_genome = genome + "\n" + fasta + "\n"

    # Concatenate with additional FASTA sequence(s), if any
    if additional_fasta_path:
        logger.info(f"Looking for additional FASTA sequence(s) to concatenate with {genome_name}")
        add_fasta_name = os.path.basename(additional_fasta_path)
        logger.info(f"Concatenating {add_fasta_name} with the genome .FASTA file")
        with open(additional_fasta_path, "r") as f:
            add_fasta = f.read()

        # Check line length
        if len(add_fasta.split("\n")[1]) != lg:
            logger.warning(f"Line length of {add_fasta_name} is not equal to {lg}")

            # remove existing line breaks and add new ones
            add_fasta = add_fasta.replace("\n", "")
            add_fasta = "\n".join([add_fasta[i:i+lg] for i in range(0, len(add_fasta), lg)])

        # Concatenate the strings
        new_genome += "\n" + add_fasta

    new_genome_output = genome_input.replace(".fa", "_artificial.fa")
    with open(new_genome_output, "w") as f:
        f.write(new_genome)

    logger.info(f"Artificial chromosome created and inserted at the end of the genome .FASTA file")

def get_stats(
        contacts_unbinned_path: str,
        sparse_mat_path: str,
        chr_coord_path: str,
        oligo_capture_with_frag_path: str,
        output_dir: str = None,
        cis_range: int = 50000,
        force: bool = False,
):
    """
    Generate statistics for contacts made by each probe.

    It generates 3 outputs file (.tsv):
    - contacts_statistics.tsv: contains different kinds of statistics for each probe.
    - norm_chr_freq.tsv: contains the normalized contacts for each probe on each chromosome.
    - norm_inter_chr_freq.tsv: contains the normalized contacts for each probe on each chromosome except its own.

    Parameters
    ----------

    contacts_unbinned_path : str
        Path to the unbinned_contacts.tsv file (generated by fragments).
    sparse_mat_path : str
        Path to the sparse_contacts_input.txt file (generated by hicstuff).
    chr_coord_path : str
        Path to the input chr_centros_coordinates.tsv file.
    oligo_capture_with_frag_path : str
        Path to the oligo input CSV file.
        Must contain fragments associated Made with the 'associate' command
    cis_range: int, default=50000
        Cis range to be considered around the probe.
    output_dir : str
        Path to the output directory.
    force : bool
        Force the overwriting of the output file if the file exists.

    Returns
    -------

    None
    """

    logger.info("[Stats] : Generating statistics for contacts made by each probe.")

    utils.check_if_exists(contacts_unbinned_path)
    utils.check_if_exists(sparse_mat_path)
    utils.check_if_exists(chr_coord_path)
    utils.check_if_exists(oligo_capture_with_frag_path)

    if output_dir is None:
        output_dir = os.path.dirname(contacts_unbinned_path)

    sample_name = os.path.basename(sparse_mat_path).split('.')[0]
    out_stats_path = os.path.join(output_dir, f"{sample_name}_statistics.tsv")
    out_chr_freq_path = os.path.join(output_dir, f"{sample_name}_norm_chr_freq.tsv")
    out_inter_chr_freq_path = os.path.join(output_dir, f"{sample_name}_norm_inter_chr_freq.tsv")

    if os.path.exists(out_stats_path) and not force:
        logger.warning(f"[Stats] : Output file already exists: {out_stats_path}")
        logger.warning("[Stats] : Use the --force / -F flag to overwrite the existing file.")
        return

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)

    coords_delim = "\t" if chr_coord_path.endswith(".tsv") else ","
    df_coords: pd.DataFrame = pd.read_csv(chr_coord_path, sep=coords_delim, index_col=None)

    chr_size_dict = {k: v for k, v in zip(df_coords['chr'], df_coords['length'])}
    chr_list = list(chr_size_dict.keys())

    df_unbinned_contacts: pd.DataFrame = pd.read_csv(contacts_unbinned_path, sep='\t')
    df_unbinned_contacts = df_unbinned_contacts.astype(dtype={'chr': str, 'start': int, 'sizes': int})

    df_sparse_contacts: pd.DataFrame = \
        pd.read_csv(sparse_mat_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])
    #   from sparse_matrix (hicstuff results): get total contacts from which probes enrichment is calculated
    total_sparse_contacts = sum(df_sparse_contacts["contacts"])

    chr_contacts_nrm = {k: [] for k in chr_size_dict}
    chr_inter_only_contacts_nrm = {k: [] for k in chr_size_dict}

    df_stats: pd.DataFrame = pd.DataFrame(columns=[
        "probe", "chr", "fragment", "type", "contacts",
        "coverage_over_hic_contacts", "cis", "trans",
        "intra_chr", "inter_chr"])

    probes = df_oligo['name'].to_list()
    fragments = df_oligo['fragment'].astype(str).to_list()
    for index, (probe, frag) in enumerate(zip(probes, fragments)):
        df_stats.loc[index, "probe"] = probe
        df_stats.loc[index, "fragment"] = frag
        df_stats.loc[index, "type"] = df_oligo.loc[index, "type"]

        #  get the probe's original coordinates
        self_chr_ori = df_oligo.loc[index, "chr_ori"]
        self_start_ori = df_oligo.loc[index, "start_ori"]
        self_stop_ori = df_oligo.loc[index, "stop_ori"]

        df_stats.loc[index, "chr"] = self_chr_ori

        sub_df = df_unbinned_contacts[['chr', 'start', 'sizes', frag]]
        sub_df.insert(3, 'end', sub_df['start'] + sub_df['sizes'])
        cis_start = self_start_ori - cis_range
        cis_stop = self_stop_ori + cis_range

        probe_contacts = sub_df[frag].sum()
        df_stats.loc[index, "contacts"] = probe_contacts
        df_stats.loc[index, 'coverage_over_hic_contacts'] = probe_contacts / total_sparse_contacts
        probes_contacts_inter = sub_df.query("chr != @self_chr_ori")[frag].sum()

        if probe_contacts > 0:
            cis_freq = sub_df.query("chr == @self_chr_ori & start >= @cis_start & end <= @cis_stop")[frag].sum()
            cis_freq /= probe_contacts

            trans_freq = 1 - cis_freq
            inter_chr_freq = probes_contacts_inter / probe_contacts
            intra_chr_freq = 1 - inter_chr_freq
        else:
            cis_freq = 0
            trans_freq = 0
            inter_chr_freq = 0
            intra_chr_freq = 0

        df_stats.loc[index, "cis"] = cis_freq
        df_stats.loc[index, "trans"] = trans_freq
        df_stats.loc[index, "intra_chr"] = intra_chr_freq
        df_stats.loc[index, "inter_chr"] = inter_chr_freq

        for chrom in chr_list:
            #   n1: sum contacts chr_i
            #   d1: sum contacts all chr
            #   chrom_size: chr_i's size
            #   genome_size: sum of sizes for all chr except frag_chr
            #   c1: normalized contacts on chr_i for frag_j
            chrom_size = chr_size_dict[chrom]
            genome_size = sum([s for c, s in chr_size_dict.items() if c != self_chr_ori])
            n1 = sub_df.loc[sub_df['chr'] == chrom, frag].sum()
            if n1 == 0:
                chr_contacts_nrm[chrom].append(0)
            else:
                d1 = probe_contacts
                c1 = (n1 / d1) / (chrom_size / genome_size)
                chr_contacts_nrm[chrom].append(c1)

            #   n2: sum contacts chr_i if chr_i != probe_chr
            #   d2: sum contacts all inter chr (exclude the probe_chr)
            #   c2: normalized inter chr contacts on chr_i for frag_j
            n2 = sub_df.loc[
                (sub_df['chr'] == chrom) &
                (sub_df['chr'] != self_chr_ori), frag].sum()

            if n2 == 0:
                chr_inter_only_contacts_nrm[chrom].append(0)
            else:
                d2 = probes_contacts_inter
                c2 = (n2 / d2) / (chrom_size / genome_size)
                chr_inter_only_contacts_nrm[chrom].append(c2)

    #  capture_efficiency_vs_dsdna: amount of contact for one oligo divided
    #  by the mean of all other 'ds' oligo in the genome
    n3 = df_stats.loc[:, 'contacts']
    d3 = np.mean(df_stats.loc[df_stats['type'] == 'ds', 'contacts'])
    df_stats['dsdna_norm_capture_efficiency'] = n3 / d3

    df_chr_nrm = pd.DataFrame({
        "probe": probes, "fragment": fragments, "type": df_oligo["type"].values
    })

    df_chr_inter_only_nrm = df_chr_nrm.copy(deep=True)

    for chr_id in chr_list:
        df_chr_nrm[chr_id] = chr_contacts_nrm[chr_id]
        df_chr_inter_only_nrm[chr_id] = chr_inter_only_contacts_nrm[chr_id]

    df_stats.to_csv(out_stats_path, sep='\t', index=False)
    df_chr_nrm.to_csv(out_chr_freq_path, sep='\t', index=False)
    df_chr_inter_only_nrm.to_csv(out_inter_chr_freq_path, sep='\t', index=False)

    logger.info(f"[Stats] : Statistics saved to {out_stats_path}")
    logger.info(f"[Stats] : Normalized chr contacts saved to {out_chr_freq_path}")
    logger.info(f"[Stats] : Normalized inter-only chr contacts saved to {out_inter_chr_freq_path}")

def filter_contacts(
        sparse_mat_path: str,
        oligo_capture_path: str,
        fragments_list_path: str,
        output_path: str = None,
        force: bool = False
) -> None:

    """
    Filter the sparse matrix by creating a nex table that contains only pairs of fragments that have an oligo
    either on frag_a or frag_b.

    The output table will contain the following columns:
    - frag_a: fragment id of the first fragment
    - frag_b: fragment id of the second fragment
    - contacts: number of contacts between the two fragments
    - chr_a: chromosome of the first fragment
    - start_a: start position of the first fragment
    - end_a: end position of the first fragment
    - size_a: size of the first fragment
    - gc_content_a: gc content of the first fragment
    - type_a: type of the oligo on the first fragment
    - name_a: name of the oligo on the first fragment
    - sequence_a: sequence of the oligo on the first fragment
    - chr_b: chromosome of the second fragment
    - start_b: start position of the second fragment
    - end_b: end position of the second fragment
    - size_b: size of the second fragment
    - gc_content_b: gc content of the second fragment
    - type_b: type of the oligo on the second fragment
    - name_b: name of the oligo on the second fragment
    - sequence_b: sequence of the oligo on the second fragment

    Parameters
    ----------
    sparse_mat_path : str
        Path to the sparse matrix file (hicstuff given output).

    oligo_capture_path : str
        Path to the oligo capture file (sshicstuff mandatory table).

    fragments_list_path : str
        Path to the fragments list file (hicstuff given output).

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
        output_path = sparse_mat_path.replace(".txt", "_filtered.tsv")

    out_basedir = os.path.dirname(output_path)
    if not os.path.exists(out_basedir):
        os.makedirs(out_basedir)

    if not force and os.path.exists(output_path):
        logger.warning(f"[Filter] : Output file already exists: {output_path}")
        logger.warning("[Filter] : Use the --force / -F flag to overwrite the existing file.")
        return

    utils.check_if_exists(sparse_mat_path)
    utils.check_if_exists(oligo_capture_path)
    utils.check_if_exists(fragments_list_path)

    utils.check_file_extension(sparse_mat_path, ".txt")
    utils.check_file_extension(oligo_capture_path, [".csv", ".tsv"])
    utils.check_file_extension(fragments_list_path, ".txt")

    df_fragments: pd.DataFrame = fragments_correction(fragments_list_path)
    df_oligo: pd.DataFrame = oligo_correction(oligo_capture_path)
    df_contacts: pd.DataFrame = sparse_mat_correction(sparse_mat_path)

    """
    Joining of the 3 dataframes
    """

    df_oligo_fragments = oligo_fragments_joining(df_fragments, df_oligo)
    df1 = second_join('a', df_fragments, df_oligo_fragments, df_contacts)
    df2 = second_join('b', df_fragments, df_oligo_fragments, df_contacts)
    df_contacts_joined = pd.concat([df1, df2])
    df_contacts_joined.drop("frag", axis=1, inplace=True)
    df_contacts_joined.sort_values(by=['frag_a', 'frag_b', 'start_a', 'start_b'], inplace=True)
    df_contacts_filtered = df_contacts_joined.convert_dtypes().reset_index(drop=True)

    df_contacts_filtered.to_csv(output_path, sep='\t', index=False)

    logger.info(f"[Filter] : Filtered contacts saved to {output_path}")


def first_join(x: str, oligo_fragments: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Join the contacts and oligo_fragments DataFrames, keeping only the rows that have their 'x' fragment
    (either 'frag_a' or 'frag_b', see contacts_correction function).

    Parameters
    ----------
    x : str
        Either 'a' or 'b', indicating whether to join on 'frag_a' or 'frag_b'.
    oligo_fragments : pd.DataFrame
        The joined oligo and fragments DataFrame.
    contacts : pd.DataFrame
        The corrected contacts DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined contacts and oligo_fragments DataFrame.
    """

    joined = contacts.merge(oligo_fragments, left_on='frag_'+x, right_on='frag', how='inner')
    return joined


def fragments_correction(fragments_path):
    fragments = pd.read_csv(fragments_path, sep='\t')
    fragments = pd.DataFrame(
        {'frag': [k for k in range(len(fragments))],
         'chr': fragments['chrom'],
         'start': fragments['start_pos'],
         'end': fragments['end_pos'],
         'size': fragments['size'],
         'gc_content': fragments['gc_content']
         })

    return fragments


def merge_sparse_mat(output_path: str = None, force: bool = False, matrices: list[str] = None) -> None:
    if not matrices:
        logger.error("No sparse matrices provided")
        return

    N = len(matrices)
    logger.info("Merging {0} sparse matrices into one".format(N))

    if os.path.exists(output_path) and not force:
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    for i in range(N):
        utils.check_file_extension(matrices[i], ".txt")

    if not output_path:
        now_ = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        output_path = os.path.join(os.path.dirname(matrices[0]), f"{now_}_merged_sparse_contacts.tsv")

    df_sparses: list[pd.DataFrame] = [
        pd.read_csv(matrix, sep='\t', header=0) for matrix in matrices
    ]

    n_frags = int(df_sparses[0].columns[0])
    if not all([int(df.columns[0]) == n_frags for df in df_sparses]):
        logger.error("All the sparse matrices must have the same number of fragments")
        return

    else:
        for df in df_sparses:
            df.columns = ["frag_a", "frag_b", "contacts"]

    df_concatenated: pd.DataFrame = pd.concat(df_sparses)
    df_merged = df_concatenated.groupby(['frag_a', 'frag_b'], as_index=False)['contacts'].sum()

    df_merged.columns = [n_frags, n_frags, len(df_merged)+1]
    df_merged.to_csv(output_path, sep='\t', index=False, header=True)
    logger.info(f"Merged sparse matrix saved to {output_path}")


def oligo_correction(oligo_path):
    delim = "," if oligo_path.endswith(".csv") else "\t"
    oligo = pd.read_csv(oligo_path, sep=delim)
    column_to_keep = ['chr', 'start', 'end', 'name', 'type', 'sequence']
    oligo = oligo[column_to_keep]
    oligo.columns = [oligo.columns[i].lower() for i in range(len(oligo.columns))]
    oligo.sort_values(by=['chr', 'start'], inplace=True)
    oligo.reset_index(drop=True, inplace=True)

    return oligo


def oligo_fragments_joining(fragments: pd.DataFrame, oligo: pd.DataFrame) -> pd.DataFrame:
    """
    Join the oligo and fragments DataFrames, removing fragments that do not contain an oligo region.

    Updates the start and end columns with the corresponding fragment positions.

    Parameters
    ----------
    fragments : pd.DataFrame
        The fragments DataFrame.
    oligo : pd.DataFrame
        The oligo DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined oligo and fragments DataFrame.
    """
    oligo = starts_match(fragments, oligo)
    oligo.set_index(['chr', 'start'])
    oligo.pop("end")
    fragments.set_index(['chr', 'start'])
    oligo_fragments = fragments.merge(oligo, on=['chr', 'start'])
    oligo_fragments.sort_values(by=['chr', 'start'])
    return oligo_fragments


def plot_profiles(
        profile_contacts_path: str,
        oligo_capture_path: str,
        chr_coord_path: str,
        output_dir: str = None,
        extension: str = "pdf",
        rolling_window: int = 1,
        region: str = None,
        log_scale: bool = False,
        user_y_min: float = None,
        user_y_max: float = None,
        width: int = 1200,
        height: int = 600
):

    profile_type = 'contacts'
    if 'frequencies' in profile_contacts_path:
        profile_type = 'frequencies'

    df: pd.DataFrame = pd.read_csv(profile_contacts_path, sep='\t')
    frags_col = df.filter(regex=r'^\d+$|^\$').columns.to_list()
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_path, sep=',')
    probes_to_frag = dict(zip(df_oligo['fragment'].astype(str), df_oligo['name'].astype(str)))
    df_coords = pd.read_csv(chr_coord_path, sep='\t')


    full_genome_size = df_coords.loc[:, 'length'].sum()
    x_min = 0
    x_max = full_genome_size

    b = re.search(r'_(\d+)kb_profile_', profile_contacts_path).group(1)
    if b == 0:
        binsize = 0
        bin_suffix = "0kb"
    else:
        binsize = int(b) * 1000  # kb to bp
        bin_suffix = f"{b}kb"

    sample_name = os.path.basename(profile_contacts_path).split(f'_{bin_suffix}_')[0]

    # Create the output directory
    if not output_dir:
        output_dir = os.path.dirname(profile_contacts_path)

    output_dir = os.path.join(output_dir, "plots")
    output_dir = os.path.join(output_dir, bin_suffix)
    if log_scale:
        output_dir = os.path.join(output_dir, "log")
    else:
        output_dir = os.path.join(output_dir, "raw")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    chr_list_unique = pd.unique(df.chr)
    n_chr = len(chr_list_unique)

    # Set coordinates preference
    coord_mode = ["genomic", "unbinned"]
    x_col = "genome_start" if binsize == 0 else "genome_bins"
    if region:
        coord_mode[0] = "chromosomal"
        region_list = region.split("-")
        if len(region_list) == 1:
            chr_, start_, end_ = region_list[0], "", ""
        else:
            chr_, start_, end_ = region_list
        max_chr_region_len = df_coords.loc[df_coords.chr == chr_, 'length'].values[0]
        x_col = "start" if binsize == 0 else "chr_bins"

        if not start_ and not end_:
            df = df[df['chr'] == chr_]
            df_coords = df_coords[df_coords['chr'] == chr_]
            x_min = 0
            x_max = max_chr_region_len
        else:
            start_, end_ = int(start_), int(end_)
            x_min = start_
            x_max = end_

        df = df[(df['chr'] == chr_) & (df[x_col] >= x_min) & (df[x_col] <= x_max)]
        x_col = "start" if binsize == 0 else "chr_bins"

    if binsize > 0:
        df_bins = graph.build_bins_template(df_coords, binsize)
        chr_bins = df_bins.chr_bins.values
        genome_bins = df_bins.genome_bins.values
        n_bins = len(chr_bins)
        if region:
            df_bins = df_bins[df_bins['chr'] == chr_]

        coord_mode[1] = "binned"
        x_min = x_min // binsize * binsize
        x_max = x_max // binsize * binsize + binsize
        df_bins = df_bins[(df_bins['chr_bins'] >= x_min) & (df_bins['chr_bins'] <= x_max)]

        if rolling_window > 1:
            for chr_ in chr_list_unique:
                df.loc[df['chr'] == chr_, frags_col] = (
                    df.loc[df['chr'] == chr_, frags_col].rolling(window=rolling_window, min_periods=1).mean())

    y_min = float(user_y_min) if user_y_min else 0.
    y_max = float(user_y_max) if user_y_max else df[frags_col].max().max()

    log_status = ""
    if log_scale:
        log_status = "log_"
        data = df[frags_col].values
        data[data == 0] = np.nan
        new_data = np.log10(data)
        y_min = np.nanmin(new_data) if not user_y_max else float(user_y_max)
        y_max = np.nanmax(new_data) if not user_y_min else float(user_y_min)
        df[frags_col] = new_data

    y_ticks = np.linspace(y_min, y_max, 5)
    y_tick_text = [f"{tick:.3f}" for tick in y_ticks]

    colors_rgba = colors.generate('rgba', len(frags_col))

    if region:
        for ii_f, frag in enumerate(frags_col):
            fig = go.Figure()
            probe = probes_to_frag.get(frag, "")
            output_path = os.path.join(
                output_dir,
                f"{sample_name}_{frag}_{probe}_{profile_type}_{bin_suffix}_{log_status}{region}.{extension}"
            )

            fig.add_trace(
                go.Scattergl(
                    x=df[x_col],
                    y=df[frag],
                    name=frag,
                    mode='lines',
                    line=dict(width=1, color=colors_rgba[ii_f]),
                    marker=dict(size=4),
                    showlegend=False
                )
            )

            fig.update_layout(
                title=f"{sample_name}",
                xaxis=dict(
                    title=dict(text=f"{region} coordinates", standoff=80),
                    tickformat='d',
                    range=[x_min, x_max],
                    showgrid=False,
                ),
                yaxis=dict(
                    title=f"{profile_type.capitalize()}{' - log' if log_scale else ''}",
                    tickvals=y_ticks,
                    ticktext=y_tick_text,
                    range=[y_min, y_max],
                    showgrid=False,
                ),
                xaxis_type='linear',
                xaxis_tickformat="d",
                yaxis_tickformat='%.4e' if log_scale else '%.4d',
                plot_bgcolor='white',
                paper_bgcolor='white',
                width=width,
                height=height,
            )

            pio.write_image(fig, output_path, engine="kaleido")

    else:
        df_10kb_tmp = graph.build_bins_template(df_coords=df_coords, bin_size=10000)
        colorbar, chr_ticks_pos = graph.colorbar_maker(df_10kb_tmp)

        # plot for each prob or group of probes
        for ii_f, frag in enumerate(frags_col):
            probe = probes_to_frag.get(frag, "")
            output_path = os.path.join(
                output_dir, f"{sample_name}_{frag}_{probe}_{profile_type}_{bin_suffix}_{log_status}.{extension}"
            )

            fig = make_subplots(
                rows=2, cols=1, row_heights=[0.94, 0.06], vertical_spacing=0.06,
                specs=[[{'type': 'scatter'}], [{'type': 'bar'}]]
            )

            fig.add_trace(
                go.Scatter(
                    x=df[x_col],
                    y=df[frag],
                    mode='lines',
                    line=dict(width=1, color=colors_rgba[ii_f]),
                    marker=dict(size=4),
                    showlegend=False
                ),
                row=1, col=1
            )

            fig.add_trace(colorbar, row=2, col=1)
            title = f" {sample_name} <br> {frag} {'- ' + probe}"
            fig.update_layout(
                title=title,
                xaxis=dict(
                    title=dict(text="Genomic coordinates", standoff=80),
                    tickformat='d',
                    range=[x_min, x_max],
                    showgrid=False,
                ),
                xaxis2=dict(
                    tickmode='array',
                    tickvals=chr_ticks_pos,
                    ticktext=df['chr'].unique(),
                    tickfont=dict(size=12),
                ),
                yaxis=dict(
                    title=f"{profile_type.capitalize()}{' - log' if log_scale else ''}",
                    tickvals=y_ticks,
                    ticktext=y_tick_text,
                    range=[y_min, y_max],
                    showgrid=False,
                ),
                yaxis2=dict(
                    showticklabels=False,
                ),

                xaxis_showgrid=False,
                yaxis_showgrid=False,
                xaxis_type='linear',
                xaxis_tickformat="d",
                yaxis_tickformat='%.4e' if log_scale else '%.4d',
                xaxis_range=[x_min, x_max],
                yaxis_range=[y_min, y_max],
                hovermode='closest',
                plot_bgcolor='white',
                paper_bgcolor='white',
                width=width,
                height=height,
            )

            pio.write_image(fig, output_path, engine="kaleido")

def profile_contacts(
        filtered_table_path: str,
        oligo_capture_with_frag_path: str,
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

    utils.check_if_exists(filtered_table_path)
    utils.check_if_exists(oligo_capture_with_frag_path)
    utils.check_if_exists(chromosomes_coord_path)

    if not output_path:
        output_path = filtered_table_path.replace("filtered.tsv", "0kb_profile_contacts.tsv")

    basedir = os.path.dirname(output_path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    if os.path.exists(output_path) and not force:
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    chr_coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords: pd.DataFrame = pd.read_csv(chromosomes_coord_path, sep=chr_coord_delim, index_col=None)
    df_chr_len = df_coords[["chr", "length"]]
    chr_list = list(df_chr_len['chr'].unique())
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    probes = df_oligo['name'].to_list()
    fragments = df_oligo['fragment'].astype(str).to_list()

    df: pd.DataFrame = pd.read_csv(filtered_table_path, sep='\t')
    df_contacts: pd.DataFrame = pd.DataFrame(columns=['chr', 'start', 'sizes'])
    df_contacts: pd.DataFrame = df_contacts.astype(dtype={'chr': str, 'start': int, 'sizes': int})

    for x in ['a', 'b']:
        y = utils.frag2(x)
        df2 = df[~pd.isna(df['name_' + x])]

        for probe in probes:
            if probe not in pd.unique(df2['name_' + x]):
                tmp = pd.DataFrame({
                    'chr': [np.nan],
                    'start': [np.nan],
                    'sizes': [np.nan],
                    probe: [np.nan]})

            else:
                df3 = df2[df2['name_' + x] == probe]
                tmp = pd.DataFrame({
                    'chr': df3['chr_' + y],
                    'start': df3['start_' + y],
                    'sizes': df3['size_' + y],
                    probe: df3['contacts']})

            df_contacts = pd.concat([df_contacts, tmp])

    group = df_contacts.groupby(by=['chr', 'start', 'sizes'], as_index=False)
    df_contacts: pd.DataFrame = group.sum()
    df_contacts = utils.sort_by_chr(df_contacts, chr_list, 'chr', 'start')
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
        utils.make_groups_of_probes(df_additional, df_contacts, probes_to_fragments)
        if normalize:
            utils.make_groups_of_probes(df_additional, df_frequencies, probes_to_fragments)

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

    utils.check_if_exists(contacts_unbinned_path)
    utils.check_if_exists(chromosomes_coord_path)

    bin_suffix = f'{bin_size // 1000}kb'
    if not output_path:
        output_path = contacts_unbinned_path.replace("0kb_profile", f"{bin_suffix}_profile")

    if os.path.exists(output_path) and not force:
        logger.warning(f"[Rebin] : Output file already exists: {output_path}")
        logger.warning("[Rebin] : Use the --force / -F flag to overwrite the existing file.")
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
        'genome_bins': np.arange(0, len(chr_bins) * bin_size, bin_size)
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

    fragments_columns = df.filter(regex=r'^\d+$|^\$').columns.to_list()

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = utils.sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned, on=['chr', 'chr_bins'], how='left')
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    df_binned.to_csv(output_path, sep='\t', index=False)


def second_join(
        x: str, fragments: pd.DataFrame, oligo_fragments: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Add the fragments DataFrame information (=columns) for the y fragment after the first join
    (see first_join function). This is only for the y fragment, because the x fragments already have their
    information in the oligo_fragments DataFrame.

    Parameters
    ----------
    x : str
        Either 'a' or 'b', indicating which fragment corresponds to an oligo.
    fragments : pd.DataFrame
        The corrected fragments DataFrame.
    oligo_fragments : pd.DataFrame
        The joined oligo and fragments DataFrame.
    contacts : pd.DataFrame
        The corrected contacts DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined DataFrame with added fragment information for the y fragment.
    """
    new_contacts = first_join(x, oligo_fragments, contacts)
    y = utils.frag2(x)
    joined = new_contacts.join(fragments.drop("frag", axis=1),
                               on='frag_'+y,
                               lsuffix='_' + x[-1],
                               rsuffix='_' + y[-1], how='left')

    # puts a suffix to know what fragment corresponds to an oligo
    joined.rename(columns={"type": "type_" + x[-1],
                           "name": "name_" + x[-1],
                           "sequence": "sequence_" + x[-1]
                           },
                  inplace=True)
    return joined


def sparse_mat_correction(sparse_mat_path):
    """
    Re-organizes the sparse contacts matrix file
    """
    contacts = pd.read_csv(sparse_mat_path, sep='\t', header=None)
    contacts.drop([0], inplace=True)
    contacts.reset_index(drop=True, inplace=True)
    contacts.columns = ['frag_a', 'frag_b', 'contacts']

    return contacts


def sparse_with_dsdna_only(
        sample_sparse_mat: str,
        oligo_capture_with_frag_path: str,
        n_flanking_dsdna: int = 2,
        output_path: str = None,
        force: bool = False
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
        logger.info(f"[Sparse Matrix Graal (dsdna)] Output file already exists: {output_path}")
        logger.warning("[Sparse Matrix Graal (dsdna)] Use the --force / -F flag to overwrite the existing file.")
        return

    utils.check_if_exists(sample_sparse_mat)
    utils.check_if_exists(oligo_capture_with_frag_path)
    utils.check_file_extension(sample_sparse_mat, ".txt")
    utils.check_file_extension(oligo_capture_with_frag_path, [".csv", ".tsv"])

    oligo_capture_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep='\t', header=None)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_capture_delim)

    df_contacts_dsdna_only = df_sparse_mat.copy(deep=True)

    ssdna_frag = df_oligo.loc[df_oligo["type"] == "ss", "fragment"].tolist()
    df_ssdna = pd.DataFrame(ssdna_frag, columns=['fragments'])

    dsdna_frag = df_oligo.loc[df_oligo["type"] == "ds", "fragment"].tolist()
    dsdna_frag_flanking = []
    for f in dsdna_frag:
        for i in range(1, n_flanking_dsdna + 1):
            dsdna_frag_flanking.append(f + i)
            dsdna_frag_flanking.append(f - i)

    dsdna_frag_all = np.unique(dsdna_frag + dsdna_frag_flanking)
    df_dsdna = pd.DataFrame(dsdna_frag_all, columns=['fragments'])

    df_frag = pd.concat([df_ssdna, df_dsdna])
    del df_ssdna, df_dsdna

    df_sparse_mat["index"] = df_sparse_mat.index
    matches_a = pd.merge(df_sparse_mat, df_frag, left_on=0, right_on='fragments', how='inner', indicator=True)
    matches_b = pd.merge(df_sparse_mat, df_frag, left_on=1, right_on='fragments', how='inner', indicator=True)
    index_to_drop = np.unique(np.concatenate((matches_a['index'].to_numpy(), matches_b['index'].to_numpy())))

    df_contacts_dsdna_only.drop(index_to_drop, inplace=True)

    df_contacts_dsdna_only.iloc[0, 0] -= len(df_frag)
    df_contacts_dsdna_only.iloc[0, 1] -= len(df_frag)
    df_contacts_dsdna_only.iloc[0, 2] -= len(index_to_drop)

    df_contacts_dsdna_only.to_csv(output_path, sep='\t', index=False, header=False)

    logger.info(f"[Sparse Matrix Graal (dsdna)] : dsDNA only contacts saved to {output_path}")


def sparse_with_ssdna_only(
        sample_sparse_mat: str,
        oligo_capture_with_frag_path: str,
        output_path: str = None,
        force: bool = False
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
        logger.info(f"[Sparse Matrix Graal (ssdna)] : Output file already exists: {output_path}")
        logger.warning("[Sparse Matrix Graal (ssdna)] : Use the --force / -F flag to overwrite the existing file.")
        return

    utils.check_if_exists(sample_sparse_mat)
    utils.check_if_exists(oligo_capture_with_frag_path)
    utils.check_file_extension(sample_sparse_mat, ".txt")
    utils.check_file_extension(oligo_capture_with_frag_path, [".csv", ".tsv"])

    oligo_capture_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep='\t', header=0)
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_capture_delim)

    df_contacts_ssdna_only = df_sparse_mat.copy(deep=True)
    ssdna_frag = pd.unique(df_oligo.loc[df_oligo["type"] == "ss", "fragment"]).tolist()

    df_contacts_ssdna_only = df_contacts_ssdna_only[
        df_contacts_ssdna_only.iloc[:, 0].isin(ssdna_frag) &
        df_contacts_ssdna_only.iloc[:, 1].isin(ssdna_frag)
    ]

    df_contacts_ssdna_only.columns = [len(ssdna_frag), len(ssdna_frag), len(df_contacts_ssdna_only)+1]
    df_contacts_ssdna_only.reset_index(drop=True, inplace=True)

    df_contacts_ssdna_only.to_csv(output_path, sep='\t', index=False, header=True)
    logger.info(f"[Sparse Matrix Graal (ssdna)] : ssDNA only contacts saved to {output_path}")



def starts_match(fragments: pd.DataFrame, oligo: pd.DataFrame) -> pd.DataFrame:
    """
    Update the start positions of the oligo DataFrame based on the corresponding fragment positions.

    If the capture oligo is inside a fragment, update the start position of the oligo DataFrame with the start
    position of the fragment.

    Parameters
    ----------
    fragments : pd.DataFrame
        The fragments DataFrame.
    oligo : pd.DataFrame
        The oligo DataFrame.

    Returns
    -------
    pd.DataFrame
        The updated oligo DataFrame.
    """
    l_starts = []
    for i in range(len(oligo)):
        oligo_chr = oligo['chr'][i]
        middle = int((oligo['end'][i] - oligo['start'][i] - 1) / 2 + oligo['start'][i] - 1)
        if oligo_chr == 'chr_artificial':
            for k in reversed(range(len(fragments))):
                interval = range(fragments['start'][k], fragments['end'][k])
                fragments_chr = fragments['chr'][k]
                if middle in interval and fragments_chr == oligo_chr:
                    l_starts.append(fragments['start'][k])
                    break
        else:
            for k in range(len(fragments)):
                interval = range(fragments['start'][k], fragments['end'][k] + 1)
                fragments_chr = fragments['chr'][k]

                if middle in interval and fragments_chr == oligo_chr:
                    l_starts.append(fragments['start'][k])
                    break
    oligo['start'] = list(l_starts)
    return oligo


def subsample(
        input_path: str,
        seed: int = 100,
        size: int = 4000000,
        compress: bool = True,
        force: bool = False
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

    utils.check_seqtk()

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
        logger.warning(f"Output file {output_path} already exists. skipping subsampling.")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    # Log the start of subsampling
    logger.info(f"Starting subsampling of {input_path} to {output_path}.gz with seed {seed} and size {size}")

    # Run seqtk to subsample the FASTQ file
    seqtk_command = f"seqtk sample -s {seed} {input_path} {size} > {output_path}"
    try:
        subprocess.run(seqtk_command, shell=True, check=True)
        logger.info(f"Subsampling completed successfully for {output_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in subsampling {input_path}: {e}")
        raise

    # Optionally compress the output file
    if compress:
        utils.check_gzip()
        if os.path.exists(output_path + ".gz"):
            logger.warning(f"Output file {output_path}.gz already exists. removing it.")
            os.remove(output_path + ".gz")

        gzip_command = f"gzip {output_path}"
        logger.info(f"Compressing {output_path}")
        try:
            subprocess.run(gzip_command, shell=True, check=True)
            output_path += ".gz"
            logger.info(f"Compression completed successfully for {output_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error in compressing {output_path}: {e}")
            raise


