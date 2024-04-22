import os
import re
import subprocess
import numpy as np
import pandas as pd

import sshicstuff.utils as sshcu
from hicstuff.log import logger

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def aggregate(
        binned_contacts_path: str,
        chr_coord_path: str,
        oligos_capture_path: str,
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
    oligos_capture_path : str
        Path to the oligos_capture.tsv file.
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

    oligos_delim = "," if oligos_capture_path.endswith(".csv") else "\t"
    coords_delim = "\t" if chr_coord_path.endswith(".tsv") else ","
    df_coords: pd.DataFrame = pd.read_csv(chr_coord_path, sep=coords_delim, index_col=None)
    df_oligos: pd.DataFrame = pd.read_csv(oligos_capture_path, sep=oligos_delim)
    df_contacts: pd.DataFrame = pd.read_csv(binned_contacts_path, sep='\t')

    binsize = int(df_contacts.loc[2, 'chr_bins'] - df_contacts.loc[1, 'chr_bins']) * 1000
    logger.info(f"Contacts binned profile with resolution of : {binsize} bp")

    chr_list = list(df_coords['chr'].unique())
    fragments = df_oligos["fragment"].astype(str).tolist()
    groups = [g for g in df_contacts.columns if g.startswith("$")]

    if len(excluded_chr_list) > 0:
        logger.info(f"Excluding chromosomes:  {', '.join(excluded_chr_list)}")
        df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr_list)]
        df_coords = df_coords[~df_coords['chr'].isin(excluded_chr_list)]

    if inter_only:
        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        logger.info("Excluding intra-chr contacts")
        for frag in fragments:
            ii_frag = df_oligos.loc[df_oligos["fragment"] == int(frag)].index[0]
            probe_chr_ori = df_oligos.loc[ii_frag, 'chr_ori']
            if probe_chr_ori not in excluded_chr_list:
                df_contacts.loc[df_contacts['chr'] == probe_chr_ori, frag] = np.nan

        output_prefix += "_inter"

    if normalize:
        logger.info("Normalizing the contacts")
        df_contacts.loc[:, fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
        output_prefix += "_norm"

    if centromeres:
        logger.info("Aggregating contacts around centromeres")
        logger.info(f"Window size: {window_size} bp on each side of the centromere")

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

        if arm_length_classification:
            if "category" not in df_coords.columns:
                logger.error(
                    "The 'category' column is missing in the centromeres file. "
                    "Must be in the form small_small or long_middle concerning lengths of left_right arms")
            else:
                logger.info("Classifying the contacts by chromosome arm lengths")
                aggregated_dir = os.path.join(output_dir, "aggregated_by_arm_sizes")

                df_arms_size: pd.DataFrame = pd.DataFrame(columns=["chr", "arm", "size", "category"])
                for _, row in df_coords.iterrows():
                    chr_ = row["chr"]
                    if chr_ not in excluded_chr_list:
                        left_, right_, category_ = row["left_arm_length"], row["right_arm_length"], row["category"]
                        if pd.isna(left_) or pd.isna(right_) or pd.isna(category_):
                            continue
                        df_arms_size.loc[len(df_arms_size)] = chr_, "left", left_, category_.split("_")[0]
                        df_arms_size.loc[len(df_arms_size)] = chr_, "right", right_, category_.split("_")[1]

                df_merged = pd.merge(df_contacts, df_telos, on='chr')
                df_merged_telos_areas_part_a = df_merged[
                    df_merged.chr_bins < (df_merged.telo_l + 3000 + 1000)]
                df_merged_telos_areas_part_a.insert(2, 'arm', 'left')
                df_merged_telos_areas_part_b = df_merged[
                    df_merged.chr_bins > (df_merged.telo_r - 3000 - 1000)]
                df_merged_telos_areas_part_b.insert(2, 'arm', 'right')

                df_telo_freq = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
                df_merged2 = pd.merge(df_telo_freq, df_arms_size, on=['chr', 'arm'])
                df_merged2.drop(columns=['telo_l', 'telo_r', 'size'], inplace=True)

                df_grouped = df_merged2.groupby(by='category', as_index=False).mean(numeric_only=True)
                df_grouped.drop(columns=['chr_bins', 'genome_bins'], inplace=True)
                df_grouped = df_grouped.rename(columns={'category': 'fragments'}).T
                df_grouped.to_csv(output_prefix + "_by_arm_sizes.tsv", sep='\t', header=False)

    else:
        return

    df_grouped = sshcu.sort_by_chr(df_grouped, chr_list, 'chr', 'chr_bins')
    df_grouped['chr_bins'] = df_grouped['chr_bins'].astype('int64')

    logger.info(f"Compute mean, median, std on the aggregated contacts per probe or group of probes, per chromosome")
    df_aggregated_mean: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).mean(numeric_only=True)
    df_aggregated_mean.to_csv(output_prefix + "_mean.tsv", sep="\t")
    df_aggregated_std: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).std(numeric_only=True)
    df_aggregated_std.to_csv(output_prefix + "_std.tsv", sep="\t")
    df_aggregated_median: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).median(numeric_only=True)
    df_aggregated_median.to_csv(output_prefix + "_median.tsv", sep="\t")

    for col in fragments + groups:
        if col in fragments:
            name = df_oligos.loc[df_oligos["fragment"] == int(col), 'name'].values[0]
        else:
            name = col[1:]

        if df_grouped[col].sum() == 0:
            continue

        df_chr_centros_pivot: pd.DataFrame = df_grouped.pivot_table(index='chr_bins',
                                                                    columns='chr', values=col, fill_value=0)

        df_chr_centros_pivot.to_csv(output_prefix + f"_{name}_per_chr.tsv", sep='\t', index=False)

    """
    Example of usage

    python3 ./main.py aggregate 
    ../data/sandbox/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree_10kb_profile_frequencies.tsv
    ../data/sandbox/S288c_DSB_LY_Capture_artificial_coordinates.tsv
    ../data/sandbox/capture_oligo_positions.csv
    -w 150000 -E chr3 -E chr2 -E 2_micron -E mitochondrion -E chr_artificial_donor -E chr_artificial_ssDNA
    -I -N -C
    """


def associate_oligo_to_frag(
        oligos_capture_path: str,
        fragments_path: str,
        frag_id_shift: int = 0,
        force: bool = True
):
    """
    Associate oligos to fragments based on the fragment name.
    It adds 3 columns directly at the end of the oligos file :
    - fragment : id of the fragment, from fragment_list (hicstuff file output)
    - fragment_start : start position of the fragment
    - fragment_end : end position of the fragment

    Parameters
    ----------
    oligos_capture_path : str
        Path to the .csv file containing the oligos.
    fragments_path : str
        Path to the .csv file containing the fragments.
    frag_id_shift : int
        Shift to apply to the fragment ids.
    force : bool
        If True, the function will overwrite the oligos file.
    Returns
    -------
    None
    """

    logger.info("Associating oligos to fragments based on the fragment id, start and end positions.")

    sshcu.check_file_extension(fragments_path, ".txt")
    sshcu.check_file_extension(oligos_capture_path, [".csv", ".tsv", ".txt"])

    # Read the oligos and fragments files
    oligos_delim = "," if oligos_capture_path.endswith(".csv") else "\t"
    df_oligos = pd.read_csv(oligos_capture_path, sep=oligos_delim)

    if "fragment" in df_oligos.columns and not force:
        logger.info("Oligos already associated to fragments. Use --force=True to overwrite.")
        return

    df_fragments = pd.read_csv(fragments_path, sep='\t')
    df_fragments['frag'] = [k for k in range(len(df_fragments))]
    df_fragments["frag"] = df_fragments["frag"] + frag_id_shift

    fragments_id = []
    fragments_start = []
    fragments_end = []
    for index, row in df_oligos.iterrows():
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

    df_oligos['fragment'] = fragments_id
    df_oligos['fragment_start'] = fragments_start
    df_oligos['fragment_end'] = fragments_end
    df_oligos.to_csv(oligos_capture_path, sep=",", index=False)

    logger.info("Oligos associated to fragments successfully.")

    """
    Example of usage:

    python3 ./main.py associate \
    ../data/inputs/capture_oligo_positions.csv \
    ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt \
    -F
    """


# def compare_with_wt(
#         stats1_path: str,
#         stats2_path: str,
#         output_dir: str = None,
# ):
#     """
#     Compare the capture efficiency of a sample with a wild-type reference.
#
#     Gives a .csv file with the with the ratio of the capture efficiency
#     of the sample over the wild-type reference.
#
#
#     Parameters
#     ----------
#     stats1_path : str
#         Path to the statistics file of the sample.
#     stats2_path : str
#         Path to the statistics file of the wild-type reference.
#     output_dir : str
#         Path to the output directory.
#
#     Returns
#     -------
#     None
#     """
#     df_sample: pd.DataFrame = pd.read_csv(stats1_path, header=0, sep="\t")
#     df_wt: pd.DataFrame = pd.read_csv(stats2_path, sep='\t')
#
#     df_cap_eff = pd.DataFrame(columns=[
#         "probe",
#         "capture_efficiency",
#         f"dsdna_norm_capture_efficiency_{wt_name}",
#         f"ratio"
#     ])
#
#     df_stats[f"capture_efficiency_vs_{wt_ref_name}"] = np.nan
#     for index, row in df_stats.iterrows():
#         probe = row['probe']
#         wt_capture_eff = df_wt.loc[df_wt['probe'] == probe, "dsdna_norm_capture_efficiency"].tolist()[0]
#
#         if wt_capture_eff > 0:
#             df_stats.loc[index, f"capture_efficiency_vs_{wt_ref_name}"] = \
#                 df_stats.loc[index, 'dsdna_norm_capture_efficiency'] / wt_capture_eff
#
#     df_stats.to_csv(statistics_path, sep='\t')

    pass


def coverage(
        sparse_mat_path: str,
        fragments_list_path: str,
        output_path: str = None,
        frag_id_shift: int = 0,
        normalize: bool = False,
        force: bool = False
) -> None:
    """
    Calculate the coverage per fragment and save the result to a bedgraph file in the output directory.

    Parameters
    ----------
    sparse_mat_path : str
        Path to the sparse_contacts_input.txt file (generated by hicstuff).
    fragments_list_path : str
        Path to the fragments_input.txt file (generated by hicstuff).
    output_path : str
        Path to the output directory.
    frag_id_shift : int
        Shift the fragment id by this value.
    normalize : bool
        Normalize the coverage by the total number of contacts.
    force : bool
        Force the overwriting of the output file if the file exists.

    Returns
    -------
    None
    """

    logger.info("Calculating coverage per fragment into a bedgraph.")

    if output_path is None:
        output_path = sparse_mat_path.replace(".txt", "_coverage.bedgraph")

    out_basedir = os.path.dirname(output_path)
    if not os.path.exists(out_basedir):
        os.makedirs(out_basedir)

    if os.path.exists(output_path) and not force:
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    df_fragments: pd.DataFrame = pd.read_csv(fragments_list_path, sep='\t')
    df_fragments.rename(columns={'chrom': 'chr', 'start_pos': 'start', 'end_pos': 'end'}, inplace=True)
    df_fragments['id'] = list(range(len(df_fragments)))
    df_fragments["id"] = df_fragments["id"] + frag_id_shift

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

    if normalize:
        logger.info("Normalizing coverage by the total number of contacts.")
        df_frequencies_cov: pd.DataFrame = df_contacts_cov.copy(deep=True)
        df_frequencies_cov["contacts"] /= sum(df_frequencies_cov["contacts"])
        df_frequencies_cov.rename(columns={"contacts": "frequencies"})
        df_frequencies_cov.to_csv(output_path, sep='\t', index=False, header=False)

    else:
        df_contacts_cov.to_csv(output_path, sep='\t', index=False, header=False)
        logger.info(f"Coverage file saved to {output_path}")

    logger.info("Coverage calculation completed.")

    """
    Example of usage:

    python3 ./main.py coverage \
    ../data/samples/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree.txt \
    ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt \
     -s 0 -F -N
    """


def get_stats(
        contacts_unbinned_path: str,
        sparse_mat_path: str,
        chr_coord_path: str,
        oligos_path: str,
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
    oligos_path : str
        Path to the oligos input CSV file.
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

    logger.info("Generating statistics for contacts made by each probe.")

    sshcu.check_if_exists(contacts_unbinned_path)
    sshcu.check_if_exists(sparse_mat_path)
    sshcu.check_if_exists(chr_coord_path)
    sshcu.check_if_exists(oligos_path)

    if output_dir is None:
        output_dir = os.path.dirname(contacts_unbinned_path)

    sample_name = os.path.basename(sparse_mat_path).split('.')[0]
    out_stats_path = os.path.join(output_dir, f"{sample_name}_statistics.tsv")
    out_chr_freq_path = os.path.join(output_dir, f"{sample_name}_norm_chr_freq.tsv")
    out_inter_chr_freq_path = os.path.join(output_dir, f"{sample_name}_norm_inter_chr_freq.tsv")

    if os.path.exists(out_stats_path) and not force:
        logger.warning(f"Output file already exists: {out_stats_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    oligos_delim = "," if oligos_path.endswith(".csv") else "\t"
    df_oligos: pd.DataFrame = pd.read_csv(oligos_path, sep=oligos_delim)

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

    probes = df_oligos['name'].to_list()
    fragments = df_oligos['fragment'].astype(str).to_list()
    for index, (probe, frag) in enumerate(zip(probes, fragments)):
        df_stats.loc[index, "probe"] = probe
        df_stats.loc[index, "fragment"] = frag
        df_stats.loc[index, "type"] = df_oligos.loc[index, "type"]

        #  get the probe's original coordinates
        self_chr_ori = df_oligos.loc[index, "chr_ori"]
        self_start_ori = df_oligos.loc[index, "start_ori"]
        self_stop_ori = df_oligos.loc[index, "stop_ori"]

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
    #  by the mean of all other 'ds' oligos in the genome
    n3 = df_stats.loc[:, 'contacts']
    d3 = np.mean(df_stats.loc[df_stats['type'] == 'ds', 'contacts'])
    df_stats['dsdna_norm_capture_efficiency'] = n3 / d3

    df_chr_nrm = pd.DataFrame({
        "probe": probes, "fragment": fragments, "type": df_oligos["type"].values
    })

    df_chr_inter_only_nrm = df_chr_nrm.copy(deep=True)

    for chr_id in chr_list:
        df_chr_nrm[chr_id] = chr_contacts_nrm[chr_id]
        df_chr_inter_only_nrm[chr_id] = chr_inter_only_contacts_nrm[chr_id]

    df_stats.to_csv(out_stats_path, sep='\t', index=False)
    df_chr_nrm.to_csv(out_chr_freq_path, sep='\t', index=False)
    df_chr_inter_only_nrm.to_csv(out_inter_chr_freq_path, sep='\t', index=False)

    logger.info(f"Statistics saved to {out_stats_path}")
    logger.info(f"Normalized chr contacts saved to {out_chr_freq_path}")
    logger.info(f"Normalized inter-only chr contacts saved to {out_inter_chr_freq_path}")

    """
    Example of usage:

    python3 ./main.py stats
    ../data/sandbox/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree_0kb_profile_contacts.tsv
    ../data/sandbox/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree.txt
    ../data/sandbox/S288c_DSB_LY_Capture_artificial_coordinates.tsv
    ../data/sandbox/capture_oligo_positions.csv
    -F
    """


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
    Create an artificial chromosome that is the concatenation of the annealing oligos and the enzyme sequence.

    Insert it at the end of the original genome .FASTA file.

    Parameters
    ----------

    annealing_input : str
        Path to the annealing oligos input CSV file.
    genome_input : str
        Path to the original genome .FASTA file.
    enzyme : str
        Restriction Enzyme sequence (e.g., dpnII sequence : gatc).
    fragment_size : int, default=150
        Size of a digested fragment / read.
    fasta_spacer : str, default="N"
        Spacer character to insert between the enzyme and the annealing oligos.
    fasta_line_length : int, default=60
        Number of characters per line in the FASTA file.
    additional_fasta_path : str, default=None
        List of additional FASTA files to concatenate with the artificial chromosome ath
        the end of the genome reference .FASTA file.
    """

    basedir = os.path.dirname(genome_input)
    artificial_chr_path = os.path.join(basedir, "chr_artificial_ssDNA.fa")

    # Creating the artificial chromosome using annealing oligos sequences
    # and the enzyme sequence
    logger.info(f"Creating the artificial chromosome with the annealing oligos and the enzyme {enzyme}")

    df = pd.read_csv(annealing_input, sep=',')
    ssDNA_seq_series = df[df["type"] == "ss"]['sequence_modified']
    ssDNA_seq = [seq.lower() for seq in ssDNA_seq_series.values]

    lg = fasta_line_length
    s = fragment_size - len(enzyme)
    p = fasta_spacer
    oneline = p * int(s/2) + enzyme + p * s

    for seq in ssDNA_seq:
        middle = len(seq) // 2
        dpnII_pos = seq.find(enzyme)
        if dpnII_pos < middle:
            seq2 = seq[dpnII_pos+len(enzyme):].upper()
        else:
            seq2 = seq[:dpnII_pos].upper()

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
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
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
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
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

    fragments_columns = df.filter(regex='^\d+$').columns.to_list()

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = sshcu.sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned, on=['chr', 'chr_bins'], how='left')
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

    sshcu.check_seqtk()

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
        sshcu.check_gzip()
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


