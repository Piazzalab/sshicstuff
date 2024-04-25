import os
from os.path import join

import sshicstuff.sshicstuff as sshic
import sshicstuff.utils as shcu
import sshicstuff.log as log

logger = log.logger

SEED = 1999
CWD = os.getcwd()
DEFAULT_DATA_DIR = join(CWD, "data")

"""
Example of usage :

python3 -m sshicstuff.main pipeline 
-c /home/nicolas/Documents/projects/sshicstuff/test_data/capture_oligo_positions.csv
-C /home/nicolas/Documents/projects/sshicstuff/test_data/chr_coords.tsv
-f /home/nicolas/Documents/projects/sshicstuff/test_data/AD162/fragments_list.txt
-m /home/nicolas/Documents/projects/sshicstuff/test_data/AD162/AD162_pcrdupkept.txt
-a /home/nicolas/Documents/projects/sshicstuff/test_data/additional_probe_groups.tsv
-b 1000 -b 2000 -b 5000 -b 10000
-E chr2 -E chr3 -E 2_micron -E mitochondrion -E chr_artificial_donor -E chr_artificial_ssDNA
-F -I -L -N 
-n 2
-o /home/nicolas/Documents/projects/sshicstuff/test_data/AD162-TEST-PIPELINE
--binning-aggregate-cen 10000
--binning-aggregate-telo 1000
--window-size-cen 150000
--window-size-telo 15000
--copy-inputs
"""


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
        frag_id_shift: int = 0,
        copy_inputs: bool = True,
        force: bool = False,
        normalize: bool = False,
):

    # Files and path alias names
    sample_name = os.path.basename(sample_sparse_mat).split('.')[0]
    if not output_dir:
        output_dir = join(DEFAULT_DATA_DIR, sample_name)

    copy_dir = join(output_dir, "inputs")
    hiconly_name = sample_name + "_hic_only.txt"
    filtered_name = sample_name + "_filtered.tsv"
    profile_0kb_contacts_name = sample_name + "_0kb_profile_contacts.tsv"
    profile_0kb_frequencies_name = sample_name + "_0kb_profile_frequencies.tsv"

    oligo_capture_with_frag = oligo_capture.replace(".csv", "_fragments_associated.csv")

    logger.info(f" -- Sample {sample_name} -- ")

    os.makedirs(output_dir, exist_ok=True)

    if copy_inputs:
        logger.info("Copying inputs file for reproducibility")
        os.makedirs(copy_dir, exist_ok=True)
        shcu.copy(sample_sparse_mat, copy_dir)
        shcu.copy(oligo_capture, copy_dir)
        shcu.copy(fragments_list, copy_dir)
        shcu.copy(chr_coordinates, copy_dir)
        if additional_groups:
            shcu.copy(additional_groups, copy_dir)

    logger.info("Associate : Associate oligo/probe name to fragment/read ID that contains it")
    sshic.associate_oligo_to_frag(
        oligo_capture_path=oligo_capture,
        fragments_path=fragments_list,
        force=force,
        frag_id_shift=frag_id_shift
    )

    logger.info("HiC only : keep only Hi-C reads, create a new sparse matrix file 'hic_only'")
    sshic.hic_only(
        sample_sparse_mat=sample_sparse_mat,
        oligo_capture_path=oligo_capture_with_frag,
        n_flanking_dsdna=n_flanking_dsdna,
        output_path=join(output_dir, hiconly_name),
        force=force
    )

    logger.info("Filter : Only keep pairs of reads that contain at least one oligo/probe")
    sshic.filter_contacts(
        sparse_mat_path=sample_sparse_mat,
        oligo_capture_path=oligo_capture,
        fragments_list_path=fragments_list,
        output_path=join(output_dir, filtered_name),
        frag_id_shift=frag_id_shift,
        force=force
    )

    logger.info("Coverage : Calculate the coverage per fragment and save the result to a bedgraph")

    sshic.coverage(
        sparse_mat_path=sample_sparse_mat,
        fragments_list_path=fragments_list,
        frag_id_shift=frag_id_shift,
        normalize=normalize,
        force=force
    )

    logger.info("Coverage : Calculate the coverage per fragment but from the")
    logger.info("'hic_only' sparse matrix and save the result to a bedgraph")
    sshic.coverage(
        sparse_mat_path=join(output_dir, hiconly_name),
        fragments_list_path=fragments_list,
        frag_id_shift=frag_id_shift,
        normalize=normalize,
        force=force
    )

    logger.info("Profile : Generate a 4C-like profile for each ssDNA oligo")
    logger.info("Basal résolution : 0 kb (unbinned)")
    sshic.profile_contacts(
        filtered_table_path=join(output_dir, filtered_name),
        oligo_capture_path=oligo_capture_with_frag,
        chromosomes_coord_path=chr_coordinates,
        normalize=normalize,
        force=force,
        additional_groups_path=additional_groups
    )

    logger.info("Stats : Make basic statistics on the contacts (inter/intra chr, cis/trans, ssdna/dsdna etc ...)")
    sshic.get_stats(
        contacts_unbinned_path=join(output_dir, profile_0kb_contacts_name),
        sparse_mat_path=sample_sparse_mat,
        chr_coord_path=chr_coordinates,
        oligo_path=oligo_capture_with_frag,
        output_dir=output_dir,
        cis_range=cis_region_size
    )

    logger.info(f"Rebin : Change bin resolution of the 4-C like profile (unbinned -> binned)")
    for bn in bin_sizes:
        bin_suffix = str(bn // 1000) + "kb"
        logger.info(f"Rebinning at {bin_suffix}")

        sshic.rebin_profile(
            contacts_unbinned_path=join(output_dir, profile_0kb_contacts_name),
            chromosomes_coord_path=chr_coordinates,
            bin_size=bn,
            force=force
        )

        sshic.rebin_profile(
            contacts_unbinned_path=join(output_dir, profile_0kb_frequencies_name),
            chromosomes_coord_path=chr_coordinates,
            bin_size=bn,
            force=force
        )

    logger.info("Aggregate : Aggregate all 4C-like profiles on centromeric regions")

    binsize_for_cen = cen_aggregated_binning
    correct_profile = profile_0kb_frequencies_name.replace(
        "_0kb_profile_", f"_{binsize_for_cen // 1000}kb_profile_"
    )

    sshic.aggregate(
        binned_contacts_path=join(output_dir, correct_profile),
        chr_coord_path=chr_coordinates,
        oligo_capture_path=oligo_capture_with_frag,
        window_size=cen_agg_window_size,
        centromeres=True,
        output_dir=output_dir,
        excluded_chr_list=excluded_chr,
        inter_only=inter_chr_only,
        normalize=normalize
    )

    logger.info("Aggregate : Aggregate all 4C-like profiles on telomeric regions")
    binsize_for_telo = telo_agg_binning
    correct_profile = profile_0kb_frequencies_name.replace(
        "_0kb_profile_", f"_{binsize_for_telo // 1000}kb_profile_"
    )

    sshic.aggregate(
        binned_contacts_path=join(output_dir, correct_profile),
        chr_coord_path=chr_coordinates,
        oligo_capture_path=oligo_capture_with_frag,
        window_size=telo_agg_window_size,
        telomeres=True,
        output_dir=output_dir,
        excluded_chr_list=excluded_chr,
        inter_only=inter_chr_only,
        normalize=normalize,
        arm_length_classification=arm_length_classification
    )

    logger.info(f"--- {sample_name} DONE --- \n\n")
