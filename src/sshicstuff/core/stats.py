"""
This module contains the functions to generate statistics for contacts made by each probe.
"""
import pandas as pd
from pathlib import Path
import logging

import sshicstuff.core.methods as methods

logger = logging.getLogger(__name__)



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

    It generates three output files (.tsv):
    - *_statistics.tsv: contact statistics per probe.
    - *_norm_chr_freq.tsv: normalized contact frequencies per chromosome.
    - *_norm_inter_chr_freq.tsv: normalized inter-chromosomal contact frequencies.

    Parameters
    ----------
    contacts_unbinned_path : str
        Path to unbinned_contacts.tsv (from fragments).
    sparse_mat_path : str
        Path to sparse_contacts_input.txt (from hicstuff).
    chr_coord_path : str
        Path to chr_centros_coordinates.tsv.
    oligo_capture_with_frag_path : str
        Path to oligo input file (CSV/TSV) with fragment associations.
    cis_range : int, optional
        Range (bp) to consider for cis contacts, by default 50000.
    output_dir : str, optional
        Output directory. Defaults to the input contacts path.
    force : bool, optional
        Overwrite existing output files, by default False.

    Returns
    -------
    None
    """

    logger.info("[Stats] : Generating statistics for contacts made by each probe.")

    methods.check_if_exists(contacts_unbinned_path)
    methods.check_if_exists(sparse_mat_path)
    methods.check_if_exists(chr_coord_path)
    methods.check_if_exists(oligo_capture_with_frag_path)

    contacts_path = Path(contacts_unbinned_path)
    output_dir = Path(output_dir) if output_dir else contacts_path.parent
    sample_name = Path(sparse_mat_path).stem

    out_stats_path = output_dir / f"{sample_name}_statistics.tsv"
    out_chr_freq_path = output_dir / f"{sample_name}_norm_chr_freq.tsv"
    out_inter_chr_freq_path = output_dir / f"{sample_name}_norm_inter_chr_freq.tsv"

    if out_stats_path.exists() and not force:
        logger.warning("[Stats] : Output file already exists: %s", out_stats_path)
        logger.warning("[Stats] : Use the --force / -F flag to overwrite the existing file.")
        return

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)

    coords_delim = "\t" if chr_coord_path.endswith(".tsv") else ","
    df_chr_coords = pd.read_csv(chr_coord_path, sep=coords_delim)
    df_chr_coords.columns = [c.lower() for c in df_chr_coords.columns]

    chr_sizes = dict(zip(df_chr_coords["chr"], df_chr_coords["length"]))
    chromosomes = list(chr_sizes.keys())

    df_contacts = pd.read_csv(contacts_unbinned_path, sep="\t").astype({"chr": str, "start": int, "sizes": int})
    df_contacts_no_artificial = df_contacts[~df_contacts["chr"].isin(["chr_artificial_ssDNA", "chr_artificial_dsDNA"])]
    df_contacts_artificial = df_contacts[df_contacts["chr"].isin(["chr_artificial_ssDNA", "chr_artificial_dsDNA"])]

    df_total_contacts = pd.read_csv(sparse_mat_path, sep="\t", header=0, names=["frag_a", "frag_b", "contacts"])
    total_contacts_count = df_total_contacts["contacts"].sum()

    normalized_contacts_per_chr = {chrom: [] for chrom in chromosomes}
    normalized_inter_contacts_per_chr = {chrom: [] for chrom in chromosomes}

    stats_records = []
    probes = df_oligo["name"].tolist()
    fragment_ids = df_oligo["fragment"].astype(str).tolist()

    for idx, (probe, fragment_id) in enumerate(zip(probes, fragment_ids)):
        probe_type = df_oligo.loc[idx, "type"]
        probe_chr = df_oligo.loc[idx, "chr"]
        probe_chr_ori = df_oligo.loc[idx, "chr_ori"]
        probe_start = df_oligo.loc[idx, "start_ori"]
        probe_end = df_oligo.loc[idx, "stop_ori"]

        record = {
            "probe": probe,
            "fragment": fragment_id,
            "type": probe_type,
            "chr": probe_chr_ori
        }

        df_probe = df_contacts[["chr", "start", "sizes", fragment_id]]
        probe_contact_sum = df_probe[fragment_id].sum()

        if probe_contact_sum > 0:
            record["contacts"] = probe_contact_sum
            record["coverage_over_hic_contacts"] = probe_contact_sum / total_contacts_count

            inter_contact_sum = df_probe.query("chr != @probe_chr_ori and chr != @probe_chr")[fragment_id].sum()
            record["inter_chr"] = inter_contact_sum / probe_contact_sum
            record["intra_chr"] = 1 - record["inter_chr"]

            # cis without artificial chromosomes
            df_probe_no_art = df_contacts_no_artificial[["chr", "start", "sizes", fragment_id]].copy()
            df_probe_no_art["end"] = df_probe_no_art["start"] + df_probe_no_art["sizes"]

            cis_start = probe_start - cis_range
            cis_end = probe_end + cis_range

            cis_contacts = df_probe_no_art.query(
                "chr == @probe_chr_ori and start >= @cis_start and end <= @cis_end"
            )[fragment_id].sum()

            cis_artificial_contacts = df_contacts_artificial[fragment_id].sum()
            total_cis = cis_contacts + cis_artificial_contacts

            record["cis"] = cis_contacts / probe_contact_sum
            record["trans"] = 1 - record["cis"]
            record["cis_with_artificial"] = total_cis / probe_contact_sum
            record["trans_with_artificial"] = 1 - record["cis_with_artificial"]

        else:
            for key in ["contacts", "coverage_over_hic_contacts", "cis", "trans",
                        "cis_with_artificial", "trans_with_artificial", "intra_chr", "inter_chr"]:
                record[key] = 0

        for chrom in chromosomes:
            chrom_size = chr_sizes[chrom]
            genome_size = sum(size for chr_, size in chr_sizes.items() if chr_ != probe_chr_ori)

            contacts_on_chrom = df_probe.loc[df_probe["chr"] == chrom, fragment_id].sum()
            inter_contacts = df_probe.loc[(df_probe["chr"] == chrom) & (df_probe["chr"] != probe_chr_ori), fragment_id].sum()

            norm_contact = (contacts_on_chrom / probe_contact_sum) / (chrom_size / genome_size) if contacts_on_chrom else 0
            norm_inter_contact = (inter_contacts / inter_contact_sum) / (chrom_size / genome_size) if inter_contacts else 0

            normalized_contacts_per_chr[chrom].append(norm_contact)
            normalized_inter_contacts_per_chr[chrom].append(norm_inter_contact)

        stats_records.append(record)

    df_stats = pd.DataFrame(stats_records)
    mean_ds_contacts = df_stats.loc[df_stats["type"] == "ds", "contacts"].mean()
    df_stats["dsdna_norm_capture_efficiency"] = df_stats["contacts"] / mean_ds_contacts

    df_chr_norm = pd.DataFrame({"probe": probes, "fragment": fragment_ids, "type": df_oligo["type"]})
    df_chr_norm_inter = df_chr_norm.copy()

    for chrom in chromosomes:
        df_chr_norm[chrom] = normalized_contacts_per_chr[chrom]
        df_chr_norm_inter[chrom] = normalized_inter_contacts_per_chr[chrom]

    df_stats.to_csv(out_stats_path, sep="\t", index=False)
    df_chr_norm.to_csv(out_chr_freq_path, sep="\t", index=False)
    df_chr_norm_inter.to_csv(out_inter_chr_freq_path, sep="\t", index=False)

    logger.info("[Stats] : Statistics saved to %s", out_stats_path)
    logger.info("[Stats] : Normalized chr contacts saved to %s", out_chr_freq_path)
    logger.info("[Stats] : Normalized inter-only chr contacts saved to %s", out_inter_chr_freq_path)
