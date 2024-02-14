import os
from os.path import join
import itertools
import argparse
import shutil
import pandas as pd
from typing import List, Optional

from core.filter import filter_contacts
from core.coverage import coverage
from core.weight import weight_mutant
from core.aggregated import aggregate
from core.statistics import get_stats, compare_to_wt
from core.binning import rebin_contacts, unbinned_contacts


class PathBundle:
    def __init__(self, sample_sparse_file_path: str, outputs_dir: str, reference_path_list: List[str] = None, ):

        self.sample_sparse_file_path = sample_sparse_file_path
        self.samp_name = sample_sparse_file_path.split("/")[-1].split(".")[0]

        self.sample_outputs_dir = join(outputs_dir, self.samp_name)
        self.sample_copy_inputs_dir = join(self.sample_outputs_dir, "inputs")
        self.not_weighted_dir = join(self.sample_outputs_dir, "not_weighted")

        self.sample_sparse_no_probe_file_path = join(self.sample_outputs_dir, self.samp_name + "_hic.txt")

        os.makedirs(self.sample_outputs_dir, exist_ok=True)
        os.makedirs(self.sample_copy_inputs_dir, exist_ok=True)
        os.makedirs(self.not_weighted_dir, exist_ok=True)

        self.filtered_contacts_input = join(self.sample_outputs_dir, "contacts_filtered.tsv")
        self.unbinned_contacts_input = join(self.not_weighted_dir, "unbinned_contacts.tsv")
        self.unbinned_frequencies_input = join(self.not_weighted_dir, "unbinned_frequencies.tsv")
        self.global_statistics_input = join(self.sample_outputs_dir, "contacts_statistics.tsv")

        self.wt_references_path = []
        self.wt_references_name = []
        self.weighted_dirs = []
        if reference_path_list:
            for ref_path in reference_path_list:
                ref_name = ref_path.split("/")[-1].split(".")[0]
                self.wt_references_path.append(ref_path)
                self.wt_references_name.append(ref_name)
                weighted_dir = join(self.sample_outputs_dir, f"weighted_{ref_name}")
                self.weighted_dirs.append(weighted_dir)
                os.makedirs(weighted_dir, exist_ok=True)


class AggregateParams:
    def __init__(self, window_size_centro, window_size_telos,
                 binning_size_centro, binning_size_telo, aggregate_by_arm_lengths,
                 excluded_chr_list, excluded_probe_chr, ):

        self.window_size_centromeres = window_size_centro
        self.window_size_telomeres = window_size_telos
        self.binning_centromeres = binning_size_centro
        self.binning_telomeres = binning_size_telo
        self.aggregate_by_arm_lengths = aggregate_by_arm_lengths
        self.excluded_probe_chr = excluded_probe_chr
        self.excluded_chr_list = excluded_chr_list


def check_and_run(output_path, func, *args):
    if not os.path.exists(output_path):
        func(*args)


def copy_file(source_path, destination_path):
    try:
        shutil.copy(source_path, destination_path)
        print(f"File {source_path.split('/')[-1]} copied successfully.")
    except IOError as e:
        print(f"Unable to copy file. Error: {e}")


def pipeline(
    path_bundle: PathBundle,
    oligos_path: str,
    fragments_list_path: str,
    centromeres_coordinates_path: str,
    binning_size_list: List[int],
    aggregate_params: AggregateParams,
    additional_groups: Optional[str] = None,
    hic_only: Optional[bool] = False
):
    print(f" -- Sample {path_bundle.samp_name} -- \n")

    copy_file(fragments_list_path, path_bundle.sample_copy_inputs_dir)
    copy_file(centromeres_coordinates_path, path_bundle.sample_copy_inputs_dir)
    copy_file(additional_groups, path_bundle.sample_copy_inputs_dir)
    copy_file(oligos_path, path_bundle.sample_copy_inputs_dir)
    copy_file(path_bundle.sample_sparse_file_path, path_bundle.sample_copy_inputs_dir)
    for rp in path_bundle.wt_references_path:
        copy_file(rp, path_bundle.sample_copy_inputs_dir)

    print("\n")

    print(f"Filter contacts \n")
    check_and_run(
        path_bundle.filtered_contacts_input, filter_contacts, oligos_path,
        fragments_list_path, path_bundle.sample_sparse_file_path, path_bundle.sample_outputs_dir, hic_only)

    print(f"Make the coverages\n")
    coverage(path_bundle.sample_sparse_file_path, fragments_list_path, path_bundle.sample_outputs_dir)
    if hic_only:
        coverage(path_bundle.sample_sparse_no_probe_file_path, fragments_list_path, path_bundle.sample_outputs_dir)

    print(f"Organize the contacts between probe fragments and the rest of the genome 'unbinned tables' \n")
    check_and_run(
        path_bundle.unbinned_contacts_input, unbinned_contacts, path_bundle.filtered_contacts_input,
        oligos_path, centromeres_coordinates_path, path_bundle.not_weighted_dir, additional_groups)

    print(f"Make basic statistics on the contacts (inter/intra chr, cis/trans, ssdna/dsdna etc ...) \n")
    check_and_run(
        path_bundle.global_statistics_input, get_stats, path_bundle.unbinned_contacts_input,
        path_bundle.sample_sparse_file_path, oligos_path, path_bundle.sample_outputs_dir)

    for rp, rn, rd in zip(path_bundle.wt_references_path, path_bundle.wt_references_name, path_bundle.weighted_dirs):
        print(f"Compare the capture efficiency with that of a wild type (may be another sample) \n")
        compare_to_wt(
            statistics_path=path_bundle.global_statistics_input,
            reference_path=rp,
            wt_ref_name=rn)

        print(f"Weight the unbinned contacts and frequencies tables by the efficiency score got on step ahead \n")
        weight_mutant(
            statistics_path=path_bundle.global_statistics_input, wt_ref_name=rn,
            contacts_path=path_bundle.unbinned_contacts_input, frequencies_path=path_bundle.unbinned_frequencies_input,
            binned_type="unbinned", output_dir=rd, additional_path=additional_groups)

    print(f"Rebin and weight the unbinned tables (contacts and frequencies) at : \n")
    for bn in binning_size_list:
        bin_suffix = str(bn // 1000) + "kb"
        print(bin_suffix)
        rebin_contacts(
            contacts_unbinned_path=path_bundle.unbinned_contacts_input,
            chromosomes_coord_path=centromeres_coordinates_path, oligos_path=oligos_path, bin_size=bn,
            output_dir=path_bundle.not_weighted_dir, additional_path=additional_groups)

        binned_contacts_input = join(path_bundle.not_weighted_dir, f"{bin_suffix}_binned_contacts.tsv")
        binned_frequencies_input = join(path_bundle.not_weighted_dir, f"{bin_suffix}_binned_frequencies.tsv")

        for rn, rd in zip(path_bundle.wt_references_name, path_bundle.weighted_dirs):
            weight_mutant(
                statistics_path=path_bundle.global_statistics_input, wt_ref_name=rn,
                contacts_path=binned_contacts_input, frequencies_path=binned_frequencies_input,
                binned_type=f"{bin_suffix}_binned", output_dir=rd,
                additional_path=additional_groups)

    print("\n")

    regions = ["centromeres", "telomeres"]
    weights_dir = [rd for rd in path_bundle.weighted_dirs] + [path_bundle.not_weighted_dir]
    normalization = [True, False]

    param_combinations = list(itertools.product(regions, weights_dir, normalization))
    for region, weight_dir, is_normalized in param_combinations:
        if region == "centromeres":
            binning_suffix = str(aggregate_params.binning_centromeres // 1000) + "kb"
            binned_contacts_path = join(weight_dir, f"{binning_suffix}_binned_frequencies.tsv")
        elif region == "telomeres":
            binning_suffix = str(aggregate_params.binning_telomeres // 1000) + "kb"
            binned_contacts_path = join(weight_dir, f"{binning_suffix}_binned_frequencies.tsv")
        else:
            continue

        output_dir = weight_dir
        ws = aggregate_params.window_size_centromeres \
            if region == "centromeres" else aggregate_params.window_size_telomeres

        print(
            f"Make an aggregated of contacts around {region} ({weight_dir.split('/')[-1]}, "
            f"{'with' if is_normalized else 'no'} normalization)")

        aggregate(
            binned_contacts_path=binned_contacts_path,
            centros_coord_path=centromeres_coordinates_path,
            oligos_path=oligos_path,
            window_size=ws,
            on=region,
            output_dir=output_dir,
            aggregate_by_arm_sizes=aggregate_params.aggregate_by_arm_lengths,
            exclude_probe_chr=aggregate_params.excluded_probe_chr,
            excluded_chr_list=aggregate_params.excluded_chr_list,
            additional_path=additional_groups,
            inter_normalization=is_normalized,
            plot=False
        )

    print(f"--- {path_bundle.samp_name} DONE --- \n\n")


def check_nan(str_):
    return str_ != str_


if __name__ == "__main__":
    #   Example command to enter for parameters (parse)
    """
    --samplesheet ../data/inputs/samplesheet.csv
    --fragments-list ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt
    --outputs-dir ../data/outputs
    --chromosomes-arms-coordinates ../data/inputs/S288c_chr_centro_coordinates.tsv 
    --oligos-capture ../data/inputs/capture_oligo_positions_v2.csv
    --additional-groups ../data/inputs/additional_probe_groups.tsv
    --binning-sizes 1000 2000 5000 10000 20000 40000 50000 80000 100000
    --centromeres-aggregated-window-size 150000  
    --telomeres-aggregated-window-size 15000
    --centromeres-aggregated-binning 10000
    --telomeres-aggregated-binning 1000
    --aggregate-by-arm-lengths 
    --excluded-chr chr2 chr3 2_micron mitochondrion chr_artificial
    --exclude-probe-chr 
    """

    """
    --samplesheet ../data/inputs/samplesheet.csv
    --fragments-list ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_v6_DpnIIHinfI.txt
    --outputs-dir ../data/outputs
    --chromosomes-arms-coordinates ../data/inputs/S288c_chr_centro_coordinates_S288c_DSB_LY_Capture_artificial_v6.tsv
    --oligos-capture ../data/inputs/capture_oligo_positions_v6.csv
    --additional-groups ../data/inputs/additional_probe_groups.tsv
    --binning-sizes 1000 2000 5000 10000 20000 40000 50000 80000 100000
    --centromeres-aggregated-window-size 150000  
    --telomeres-aggregated-window-size 15000
    --centromeres-aggregated-binning 10000
    --telomeres-aggregated-binning 1000
    --aggregate-by-arm-lengths 
    --excluded-chr chr2 chr3 2_micron mitochondrion chr_artificial_donor chr_artificial_ssDNA
    --exclude-probe-chr 
    """

    # default folders paths
    samples_dir = "../data/samples"
    references_dir = "../data/references"
    inputs_dir = "../data/inputs"
    outputs_dir = "../data/outputs"

    parser = argparse.ArgumentParser(
        description="Script that processes sshic samples data.")

    parser.add_argument('--samplesheet', type=str, required=True,
                        help='Path to the samplesheet (.csv) that contains samples and their respective references ')

    parser.add_argument('--oligos-capture', type=str, required=True,
                        help='Path to the file that contains positions of oligos')

    parser.add_argument('--fragments-list', type=str, required=True,
                        help='Path to the file fragments_list (hic_stuff output)')

    parser.add_argument('--outputs-dir', type=str, required=True,
                        help='Path to the output directory that will contain the results for each samples')

    parser.add_argument('--chromosomes-arms-coordinates', type=str, required=True,
                        help='Path to the file containing centromeres coordinates and chromosomes arms lengths')

    parser.add_argument('--binning-sizes', nargs='+', type=int, required=True,
                        help='desired bin size for the rebin step')

    parser.add_argument('--additional-groups', type=str, required=False,
                        help='Path to additional groups of probes table')

    parser.add_argument('--centromeres-aggregated-window-size', type=int, required=True,
                        help="window (in bp) that defines a focus region to aggregated centromeres")

    parser.add_argument('--centromeres-aggregated-binning', type=int, required=True,
                        help="bin size (in bp) to use for the aggregated centromeres contacts")

    parser.add_argument('--telomeres-aggregated-window-size', type=int, required=True,
                        help="window (in bp) that defines a focus region to aggregated telomeres")

    parser.add_argument('--telomeres-aggregated-binning', type=int, required=True,
                        help="bin size (in bp) to use for the aggregated telomeres contacts")

    parser.add_argument('--aggregate-by-arm-lengths', action='store_true', required=False,
                        help="aggregate contacts by arm lengths")

    parser.add_argument('--excluded-chr', nargs='+', type=str, required=False,
                        help='list of chromosomes to excludes to prevent bias of contacts')

    parser.add_argument('--hic-only', action='store_true', required=False,
                        help="remove from sparse the fragment that contains oligo sshic position")

    parser.add_argument('--exclude-probe-chr', action='store_true', required=False,
                        help="exclude the chromosome where the probe comes from (oligo's chromosome)")

    args = parser.parse_args()

    if len(args.samplesheet.split("/")) == 1:
        args.samplesheet = join(inputs_dir, args.samplesheet)
    if len(args.oligos_capture.split("/")) == 1:
        args.oligos_capture = join(inputs_dir, args.oligos_capture)
    if len(args.fragments_list.split("/")) == 1:
        args.fragments_list = join(inputs_dir, args.fragments_list)
    if len(args.chromosomes_arms_coordinates.split("/")) == 1:
        args.chromosomes_arms_coordinates = join(inputs_dir, args.chromosomes_arms_coordinates)
    if args.additional_groups and len(args.additional_groups.split("/")) == 1:
        args.additional_groups = join(inputs_dir, args.additional_groups)

    df_samplesheet: pd.DataFrame = pd.read_csv(args.samplesheet, sep=",")
    samples_and_refs_paths = {}
    for _, row in df_samplesheet.iterrows():
        samp = row.loc["sample"]
        if len(samp.split("/")) == 1:
            samp_path = join(samples_dir, samp)
        else:
            samp_path = samp
        samples_and_refs_paths[samp_path] = []

        if len(row) > 1:
            for i in range(1, len(row)):
                ref = row.iloc[i]
                if not check_nan(ref):
                    if len(ref.split("/")) == 1:
                        ref_path = join(references_dir, ref)
                    else:
                        ref_path = ref
                    samples_and_refs_paths[samp_path].append(ref_path)

    for samp in samples_and_refs_paths:
        refs = samples_and_refs_paths[samp]
        sample_path_bundle = PathBundle(samp, args.outputs_dir, refs)

        sample_aggregate_params_centros = AggregateParams(
            args.centromeres_aggregated_window_size, args.telomeres_aggregated_window_size,
            args.centromeres_aggregated_binning, args.telomeres_aggregated_binning, args.aggregate_by_arm_lengths,
            args.excluded_chr, args.exclude_probe_chr)

        sample_data = [
            sample_path_bundle, args.oligos_capture, args.fragments_list, args.chromosomes_arms_coordinates,
            args.binning_sizes, sample_aggregate_params_centros, args.additional_groups, args.hic_only]
        pipeline(*sample_data)
