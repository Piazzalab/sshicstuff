import os
import re
import logging
import numpy as np
import pandas as pd

import sshicstuff.utils as sshcu

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

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

    Returns
    -------
    None
    """

    if telomeres == centromeres:
        logging.error("You must specify either telomeres or centromeres. Not both")
        logging.error("Exiting...")
        return

    sshcu.check_if_exists(binned_contacts_path)
    sshcu.check_if_exists(chr_coord_path)
    sshcu.check_if_exists(oligos_capture_path)
    
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
    logging.info(f"Contacts binned profile with resolution of : {binsize} bp")

    chr_list = list(df_coords['chr'].unique())
    fragments = df_oligos["fragment"].astype(str).tolist()
    groups = [g for g in df_contacts.columns if g.startswith("$")]

    if len(excluded_chr_list) > 0:
        logging.info(f"Excluding chromosomes:  {', '.join(excluded_chr_list)}")
        df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr_list)]
        df_coords = df_coords[~df_coords['chr'].isin(excluded_chr_list)]

    if inter_only:
        #   We need to remove for each oligo the number of contact it makes with its own chr.
        #   Because we know that the frequency of intra-chr contact is higher than inter-chr
        #   We have to set them as NaN to not bias the average
        logging.info("Excluding intra-chr contacts")
        for frag in fragments:
            ii_frag = df_oligos.loc[df_oligos["fragment"] == int(frag)].index[0]
            probe_chr_ori = df_oligos.loc[ii_frag, 'chr_ori']
            if probe_chr_ori not in excluded_chr_list:
                df_contacts.loc[df_contacts['chr'] == probe_chr_ori, frag] = np.nan

        output_prefix += "_inter"

    if normalize:
        logging.info("Normalizing the contacts")
        df_contacts.loc[:, fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
        output_prefix += "_norm"

    if centromeres:
        logging.info("Aggregating contacts around centromeres")
        logging.info(f"Window size: {window_size} bp on each side of the centromere")

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
            abs(df_merged_telos_areas_part_b['chr_bins'] - (df_merged_telos_areas_part_b['telo_r'] // binsize) * binsize)
        df_merged_telos_areas: pd.DataFrame = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
        df_grouped: pd.DataFrame = df_merged_telos_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(
            numeric_only=True)
        df_grouped.drop(columns=['telo_l', 'telo_r', 'genome_bins'], axis=1, inplace=True)

        if arm_length_classification:
            if "category" not in df_coords.columns:
                logging.error(
                    "The 'category' column is missing in the centromeres file. "
                    "Must be in the form small_small or long_middle concerning lengths of left_right arms")
            else:
                logging.info("Classifying the contacts by chromosome arm lengths")
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
                chr_arm(
                    df_chr_arm=df_arms_size,
                    df_telos=df_telos,
                    df_contacts=df_contacts,
                    telomeres_size=30000,
                    output_path=output_prefix+"_by_arm_sizes.tsv"
                )

    else:
        return

    df_grouped = sshcu.sort_by_chr(df_grouped, chr_list, 'chr', 'chr_bins')
    df_grouped['chr_bins'] = df_grouped['chr_bins'].astype('int64')

    logging.info(f"Compute mean, median, std on the aggregated contacts per probe or group of probes, per chromosome")
    df_aggregated_mean: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).mean(numeric_only=True)
    df_aggregated_mean.to_csv(output_prefix+"_mean.tsv", sep="\t")
    df_aggregated_std: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).std(numeric_only=True)
    df_aggregated_std.to_csv(output_prefix+"_std.tsv", sep="\t")
    df_aggregated_median: pd.DataFrame = df_grouped.groupby(by="chr_bins", as_index=False).median(numeric_only=True)
    df_aggregated_median.to_csv(output_prefix+"_median.tsv", sep="\t")

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


def chr_arm(
        df_chr_arm: pd.DataFrame,
        df_telos: pd.DataFrame,
        df_contacts: pd.DataFrame,
        telomeres_size: int,
        output_path: str
):
    """
    Classify the contacts aggregate around the telomeres by chromosome arm lengths.
    """

    df_merged = pd.merge(df_contacts, df_telos, on='chr')
    df_merged_telos_areas_part_a = df_merged[df_merged.chr_bins < (df_merged.telo_l + telomeres_size + 1000)]
    df_merged_telos_areas_part_a.insert(2, 'arm', 'left')
    df_merged_telos_areas_part_b = df_merged[df_merged.chr_bins > (df_merged.telo_r - telomeres_size - 1000)]
    df_merged_telos_areas_part_b.insert(2, 'arm', 'right')

    df_telo_freq = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
    df_merged2 = pd.merge(df_telo_freq, df_chr_arm, on=['chr', 'arm'])
    df_merged2.drop(columns=['telo_l', 'telo_r', 'size'], inplace=True)

    df_grouped = df_merged2.groupby(by='category', as_index=False).mean(numeric_only=True)
    df_grouped.drop(columns=['chr_bins', 'genome_bins'], inplace=True)
    df_grouped = df_grouped.rename(columns={'category': 'fragments'}).T
    df_grouped.to_csv(output_path, sep='\t', header=False)
