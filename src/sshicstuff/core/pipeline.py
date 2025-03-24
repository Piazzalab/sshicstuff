"""
This module contains the full pipeline for the analysis of a single sample.
The pipeline is composed of the following steps:
- Create a new sparse matrix with only dsDNA reads
- Create a new sparse matrix with only ssDNA reads
- Filter the contacts to keep only the ones with at least one oligo/probe
- Generate a 4C-like profile for each ssDNA oligo
- Generate a profile containing only contacts frequencies between oligos
- Make basic statistics on the contacts (inter/intra chr, cis)
- Change bin resolution of the 4-C like profile (unbinned -> binned)
- Aggregate all 4C-like profiles on centromeric regions
- Aggregate all 4C-like profiles on telomeric regions
"""

import os
from os.path import join
from datetime import datetime

import sshicstuff.core.methods as methods
import sshicstuff.core.filter as filt
import sshicstuff.core.profile as prof
import sshicstuff.core.stats as stats
import sshicstuff.core.aggregate as agg
import sshicstuff.log as log

logger = log.logger

SEED = 1999
CWD = os.getcwd()


def full_pipeline(
    sample_sparse_mat: str,
    oligo_capture: str,
    fragments_list: str,
    chr_coordinates: str,
    output_dir: str = None,
    additional_groups: str = None,
    bin_sizes: list[int] = None,
    cen_agg_window_size: int = 15000,
    cen_aggregated_binning: int = 10000,
    telo_agg_window_size: int = 15000,
    telo_agg_binning: int = 10000,
    arm_length_classification: bool = False,
    excluded_chr: list[str] = None,
    cis_region_size: int = 50000,
    n_flanking_dsdna: int = 2,
    inter_chr_only: bool = False,
    copy_inputs: bool = True,
    force: bool = False,
    normalize: bool = False,
):
    """ "
    Run the full pipeline for a given sample
    """

    # Files and path alias names
    sample_name = os.path.basename(sample_sparse_mat).split(".")[0]
    input_basedir = os.path.dirname(sample_sparse_mat)
    if not output_dir:
        output_dir = join(input_basedir, sample_name)

    copy_dir = join(output_dir, "inputs")
    dsdnaonly_name = sample_name + "_dsdna_only.txt"
    ssdnaonly_name = sample_name + "_ssdna_only.txt"
    filtered_name = sample_name + "_filtered.tsv"
    profile_0kb_contacts_name = sample_name + "_0kb_profile_contacts.tsv"
    profile_0kb_frequencies_name = sample_name + "_0kb_profile_frequencies.tsv"

    oligo_capture_with_frag = oligo_capture.replace(".csv", "_fragments_associated.csv")

    now = datetime.now()
    now_string = now.strftime("%Y-%m-%d %H:%M:%S")

    logger.info("[START] : %s", now_string)
    logger.info("[Pipeline] : %s ", sample_name)
    os.makedirs(output_dir, exist_ok=True)

    if copy_inputs:
        os.makedirs(copy_dir, exist_ok=True)
        methods.copy(sample_sparse_mat, copy_dir)
        methods.copy(oligo_capture, copy_dir)
        methods.copy(fragments_list, copy_dir)
        methods.copy(chr_coordinates, copy_dir)
        if additional_groups:
            methods.copy(additional_groups, copy_dir)

    methods.associate_oligo_to_frag(
        oligo_capture_path=oligo_capture,
        fragments_path=fragments_list,
        force=force,
    )
    if copy_inputs:
        methods.copy(oligo_capture_with_frag, copy_dir)

    # dsDNA reads only
    logger.info(
        "[Sparse Matrix Graal (dsdna)] : creating a new sparse matrix with only dsDNA reads"
    )
    dsdna_sparse_mat = join(output_dir, dsdnaonly_name)
    methods.sparse_with_dsdna_only(
        sample_sparse_mat=sample_sparse_mat,
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        n_flanking_dsdna=n_flanking_dsdna,
        output_path=dsdna_sparse_mat,
        force=force,
    )

    logger.info("[Coverage] : Calculate the coverage for dsDNA reads only")
    methods.coverage(
        sparse_mat_path=dsdna_sparse_mat,
        fragments_list_path=fragments_list,
        normalize=normalize,
        output_dir=output_dir,
        force=force,
    )

    # ssDNA reads only
    logger.info(
        "[Sparse Matrix Graal (ssdna)] : creating a new sparse matrix with only ssDNA reads"
    )
    ssdna_sparse_mat = join(output_dir, ssdnaonly_name)
    methods.sparse_with_ssdna_only(
        sample_sparse_mat=sample_sparse_mat,
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        output_path=ssdna_sparse_mat,
        force=force,
    )

    logger.info(
        "[Coverage] : Calculate the coverage per ssDNA fragment and save the result to a bedgraph"
    )
    methods.coverage(
        sparse_mat_path=ssdna_sparse_mat,
        fragments_list_path=fragments_list,
        normalize=normalize,
        output_dir=output_dir,
        force=force,
    )

    # All reads
    logger.info(
        "[Filter] : Only keep pairs of reads that contain at least one oligo/probe"
    )
    filt.filter_contacts(
        sparse_mat_path=sample_sparse_mat,
        oligo_capture_path=oligo_capture,
        fragments_list_path=fragments_list,
        output_path=join(output_dir, filtered_name),
        force=force,
    )

    logger.info(
        "[Coverage] : Calculate the coverage per fragment and save the result to a bedgraph"
    )
    for bn in [0] + bin_sizes:
        methods.coverage(
            sparse_mat_path=sample_sparse_mat,
            fragments_list_path=fragments_list,
            normalize=normalize,
            output_dir=output_dir,
            force=force,
            bin_size=bn,
            chromosomes_coord_path=chr_coordinates,
        )

    logger.info("[Profile] : Generate a 4C-like profile for each ssDNA oligo")
    logger.info("[Profile] : Basal résolution : 0 kb (max resolution)")
    prof.profile_contacts(
        filtered_table_path=join(output_dir, filtered_name),
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        chromosomes_coord_path=chr_coordinates,
        normalize=normalize,
        force=force,
        additional_groups_path=additional_groups,
    )

    logger.info(
        "[Profile] : Generate a profile conaining only contacts frequencie between oligos"
    )
    prof.profile_probes_only(
        filtered_table_path=join(output_dir, filtered_name),
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        force=force,
    )

    logger.info(
        "[Stats] : Make basic statistics on the contacts (inter/intra chr, cis/trans, ssdna/dsdna etc ...)"
    )
    stats.get_stats(
        contacts_unbinned_path=join(output_dir, profile_0kb_contacts_name),
        sparse_mat_path=sample_sparse_mat,
        chr_coord_path=chr_coordinates,
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        output_dir=output_dir,
        cis_range=cis_region_size,
        force=force,
    )

    logger.info(
        "[Rebin] : Change bin resolution of the 4-C like profile (unbinned -> binned)"
    )
    for bn in bin_sizes:
        bin_suffix = methods.get_bin_suffix(bn)
        logger.info("[Rebin] : %s", bin_suffix)

        prof.rebin_profile(
            contacts_unbinned_path=join(output_dir, profile_0kb_contacts_name),
            chromosomes_coord_path=chr_coordinates,
            bin_size=bn,
            force=force,
        )

        prof.rebin_profile(
            contacts_unbinned_path=join(output_dir, profile_0kb_frequencies_name),
            chromosomes_coord_path=chr_coordinates,
            bin_size=bn,
            force=force,
        )

    logger.info("[Aggregate] : Aggregate all 4C-like profiles on centromeric regions")

    binsize_for_cen = cen_aggregated_binning
    correct_profile = profile_0kb_frequencies_name.replace(
        "_0kb_profile_", f"_{binsize_for_cen // 1000}kb_profile_"
    )

    agg.aggregate(
        binned_contacts_path=join(output_dir, correct_profile),
        chr_coord_path=chr_coordinates,
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        window_size=cen_agg_window_size,
        centromeres=True,
        output_dir=output_dir,
        excluded_chr_list=excluded_chr,
        inter_only=inter_chr_only,
        normalize=normalize,
    )

    logger.info("[Aggregate] : Aggregate all 4C-like profiles on telomeric regions")
    binsize_for_telo = telo_agg_binning
    correct_profile = profile_0kb_frequencies_name.replace(
        "_0kb_profile_", f"_{binsize_for_telo // 1000}kb_profile_"
    )

    agg.aggregate(
        binned_contacts_path=join(output_dir, correct_profile),
        chr_coord_path=chr_coordinates,
        oligo_capture_with_frag_path=oligo_capture_with_frag,
        window_size=telo_agg_window_size,
        telomeres=True,
        output_dir=output_dir,
        excluded_chr_list=excluded_chr,
        inter_only=inter_chr_only,
        normalize=normalize,
        arm_length_classification=arm_length_classification,
    )

    now = datetime.now()
    now_string = now.strftime("%Y-%m-%d %H:%M:%S")
    logger.info("[Pipeline] : %s  Done", {sample_name})
    logger.info("[END] : %s", now_string)
