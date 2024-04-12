import os
from os.path import join
import itertools
import logging
from typing import List

import sshicstuff.filter as shcf
import sshicstuff.coverage as shcc
import sshicstuff.weight as shcw
import sshicstuff.aggregated as shca
import sshicstuff.statistics as shcs
import sshicstuff.binning as shcb
import sshicstuff.utils as shcu


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("sshicstuff_pipeline.log"),
        logging.StreamHandler()
    ]
)


def full_pipeline(
        sample_sparse_mat: str,
        oligos_capture: str,
        fragments_list: str,
        chromosomes_arms_coordinates: str,
        out_dir: str = None,
        additional_groups: str = None,
        reference: List[str] = None,
        binning_sizes=None,
        centromeres_aggregated_window_size: int = 15000,
        centromeres_aggregated_binning: int = 10000,
        telomeres_aggregated_window_size: int = 15000,
        telomeres_aggregated_binning: int = 10000,
        aggregate_by_arm_lengths: bool = False,
        excluded_chr: List[str] = None,
        cis_region_size: int = 50000,
        hic_only: bool = False,
        exclude_probe_chr: bool = False,
        psmn_shift: bool = False,
        plot: bool = False,
        copy_inputs: bool = True,
):
    sample_name = os.path.basename(sample_sparse_mat).split('.')[0]
    logging.info(f" -- Sample {sample_name} -- ")

    if not out_dir:
        out_dir = os.getcwd()

    sample_output_dir = join(out_dir, sample_name)
    del out_dir
    sample_copy_input_dir = join(sample_output_dir, 'inputs')
    no_weight_dir = join(sample_output_dir, 'classic')
    sample_sparse_no_probe_file_path = join(sample_output_dir, sample_name + "_hic.txt")
    filter_contacts = join(sample_output_dir, f"{sample_name}_contacts_filtered.tsv")
    unbinned_contacts = join(no_weight_dir, f"{sample_name}_profile_contacts.tsv")
    unbinned_frequencies = join(no_weight_dir, f"{sample_name}_profile_frequencies.tsv")
    global_statistics = join(sample_output_dir, f"{sample_name}_contacts_statistics.tsv")

    os.makedirs(sample_output_dir, exist_ok=True)
    os.makedirs(no_weight_dir, exist_ok=True)
    os.makedirs(sample_copy_input_dir, exist_ok=True)

    logging.info("Copying inputs file for reproducibility")
    if copy_inputs:
        os.makedirs(sample_copy_input_dir, exist_ok=True)
        shcu.copy(sample_sparse_mat, sample_copy_input_dir)
        shcu.copy(oligos_capture, sample_copy_input_dir)
        shcu.copy(fragments_list, sample_copy_input_dir)
        shcu.copy(chromosomes_arms_coordinates, sample_copy_input_dir)
        if additional_groups:
            shcu.copy(additional_groups, sample_copy_input_dir)
        if reference:
            for ref in reference:
                shcu.copy(ref, sample_copy_input_dir)

    logging.info("Filter contacts")
    shcf.filter_contacts(
        sample_name=sample_name,
        oligos_path=oligos_capture,
        fragments_path=fragments_list,
        contacts_path=sample_sparse_mat,
        output_dir=sample_output_dir,
        hic_only=hic_only,
        psmn_shift=psmn_shift
    )

    logging.info("Make the coverages")
    shcc.coverage(
        sshic_contacts_path=sample_sparse_mat,
        fragments_path=fragments_list,
        output_dir=sample_output_dir,
        psmn_shift=psmn_shift
    )

    if hic_only:
        shcc.coverage(
            sshic_contacts_path=sample_sparse_no_probe_file_path,
            fragments_path=fragments_list,
            output_dir=sample_output_dir,
            psmn_shift=psmn_shift
        )

    logging.info("Organize the contacts between probe fragments and the rest of the genome")
    shcb.profile_contacts(
        sample_name=sample_name,
        filtered_contacts_path=filter_contacts,
        output_dir=no_weight_dir,
        oligos_path=oligos_capture,
        additional_path=additional_groups,
        chromosomes_coord_path=chromosomes_arms_coordinates,
    )

    logging.info("Make basic statistics on the contacts (inter/intra chr, cis/trans, ssdna/dsdna etc ...)")
    shcs.get_stats(
        sample_name=sample_name,
        contacts_unbinned_path=unbinned_contacts,
        sparse_contacts_path=sample_sparse_mat,
        centros_coord_path=chromosomes_arms_coordinates,
        oligos_path=oligos_capture,
        output_dir=sample_output_dir,
        cis_range=cis_region_size
    )

    if reference:
        for ref in reference:
            ref_name = os.path.basename(ref).split('.')[0]
            ref_output_dir = join(sample_output_dir, f"vs_{ref_name}")
            logging.info(f"Compare the capture efficiency with that of a wild type (may be another sample) \n")
            shcs.compare_to_wt(statistics_path=global_statistics, reference_path=ref, wt_ref_name=ref_name)

            logging.info(f"Weight the unbinned contacts and frequencies tables by the efficiency score got on step ahead \n")
            shcw.weight_mutant(
                statistics_path=global_statistics,
                wt_ref_name=ref_name,
                contacts_path=unbinned_contacts,
                frequencies_path=unbinned_frequencies,
                output_dir=ref_output_dir,
                additional_path=additional_groups
            )

    logging.info(f"Rebin and weight the unbinned tables (contacts and frequencies) at : \n")
    if binning_sizes:
        for bn in binning_sizes:
            bin_suffix = str(bn // 1000) + "kb"
            logging.info(f"Rebinning at {bin_suffix}")

            shcb.rebin_contacts(
                sample_name=sample_name,
                contacts_unbinned_path=unbinned_contacts,
                chromosomes_coord_path=chromosomes_arms_coordinates,
                bin_size=bn,
                oligos_path=oligos_capture,
                output_dir=no_weight_dir,
                additional_path=additional_groups,
            )

            current_binned_contacts = join(no_weight_dir, f"{sample_name}_{bin_suffix}_profile_contacts.tsv")
            current_binned_frequencies = join(no_weight_dir, f"{sample_name}_{bin_suffix}_profile_frequencies.tsv")

            if reference:
                for ref in reference:
                    ref_name = os.path.basename(ref).split('.')[0]
                    ref_output_dir = join(sample_output_dir, f"vs_{ref_name}")
                    shcs.compare_to_wt(statistics_path=global_statistics, reference_path=ref, wt_ref_name=ref_name)

                    shcw.weight_mutant(
                        statistics_path=global_statistics,
                        wt_ref_name=ref_name,
                        contacts_path=current_binned_contacts,
                        frequencies_path=current_binned_frequencies,
                        output_dir=ref_output_dir,
                        additional_path=additional_groups
                    )

    regions = ["centromeres", "telomeres"]
    normalization = [True, False]
    weights_dir = [no_weight_dir]
    if reference:
        for r in reference:
            weights_dir.append(join(sample_output_dir, f"vs_{os.path.basename(r).split('.')[0]}"))

    param_combinations = list(itertools.product(regions, weights_dir, normalization))
    for region, weight_dir, is_normalized in param_combinations:
        weight_suffix = "_" + weight_dir.split("/")[-1] if "vs" in weight_dir else ""
        if region == "centromeres":
            binning_suffix = str(centromeres_aggregated_binning // 1000) + "kb"
        elif region == "telomeres":
            binning_suffix = str(telomeres_aggregated_binning // 1000) + "kb"
        else:
            continue

        full_name = f"{sample_name}_{binning_suffix}_profile{weight_suffix}_contacts.tsv"
        binned_contacts_path = join(weight_dir, full_name)
        output_dir = weight_dir
        ws = centromeres_aggregated_window_size if region == "centromeres" else telomeres_aggregated_window_size

        logging.info(
            f"Make an aggregated of contacts around {region} ({weight_dir.split('/')[-1]}, "
            f"{'with' if is_normalized else 'no'} normalization)")

        shca.aggregate(
            binned_contacts_path=binned_contacts_path,
            centros_coord_path=chromosomes_arms_coordinates,
            oligos_path=oligos_capture,
            window_size=ws,
            on=region,
            output_dir=output_dir,
            aggregate_by_arm_sizes=aggregate_by_arm_lengths,
            exclude_probe_chr=exclude_probe_chr,
            excluded_chr_list=excluded_chr,
            additional_path=additional_groups,
            inter_normalization=is_normalized,
            plot=plot
        )

    logging.info(f"--- {sample_name} DONE --- \n\n")

