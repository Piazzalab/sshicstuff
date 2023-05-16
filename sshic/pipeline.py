import re
import os
import argparse

import pandas as pd
from typing import List

from filter import filter_contacts
from probe2fragment import associate_probes_to_fragments
from coverage import coverage
from fragments import organize_contacts
from statistics import get_stats, compare_to_wt
from binning import rebin_contacts
from ponder import ponder_mutant
from aggregated import aggregate
from utils import is_debug


class PathBundle:
    def __init__(self, sample_sparse_file_path, samples_wt_table_path: str):
        self.sample_sparse_file_path = sample_sparse_file_path
        self.sample_input_dir = os.path.dirname(sample_sparse_file_path)
        self.samp_id = re.match(r'^AD\d+', sample_sparse_file_path.split("/")[-1]).group()
        self.sample_output_dir = os.path.join(self.sample_input_dir, self.samp_id)
        os.makedirs(self.sample_output_dir, exist_ok=True)

        self.filtered_contacts_input = os.path.join(self.sample_output_dir, self.samp_id + "_filtered.tsv")
        self.cover = os.path.join(self.sample_output_dir, self.samp_id + "_coverage_per_fragment.bedgraph")
        self.unbinned_contacts_input = os.path.join(self.sample_output_dir, self.samp_id+"_unbinned_contacts.tsv")
        self.unbinned_frequencies_input = os.path.join(self.sample_output_dir, self.samp_id+"_unbinned_frequencies.tsv")
        self.global_statistics_input = os.path.join(self.sample_output_dir, self.samp_id+"_global_statistics.tsv")

        mutants_vs_wt: pd.DataFrame = pd.read_csv(samples_wt_table_path, sep="\t")
        ref = mutants_vs_wt[mutants_vs_wt["sample"] == self.samp_id]["ref"].tolist()[0]
        if re.match(r'^AD\d+', ref):
            self.wt_to_compare_path = \
                os.path.join(self.sample_input_dir, os.path.join(ref, os.path.join(f"{ref}_global_statistics.tsv")))
        else:
            self.wt_to_compare_path = os.path.join(self.sample_input_dir, f"wt/{ref}.tsv")


class AggregateParams:
    def __init__(self, window_size, excluded_probe_chr, excluded_chr_list, inter_normalization):
        self.window_size = window_size
        self.excluded_probe_chr = excluded_probe_chr
        self.excluded_chr_list = excluded_chr_list
        self.inter_normalization = inter_normalization


def check_and_run(output_path, func, *args):
    if not os.path.exists(output_path):
        func(*args)


def step_one(
        path_bundle: PathBundle,
        oligos_path: str,
        fragments_list_path: str
):

    check_and_run(path_bundle.filtered_contacts_input,
                  filter_contacts,
                  oligos_path,
                  fragments_list_path,
                  path_bundle.sample_sparse_file_path,
                  path_bundle.sample_output_dir)

    associate_probes_to_fragments(fragments_list_path, oligos_path)

    check_and_run(path_bundle.cover, coverage,
                  path_bundle.sample_sparse_file_path,
                  fragments_list_path, path_bundle.sample_output_dir)

    organize_contacts(filtered_contacts_path=path_bundle.filtered_contacts_input, oligos_path=oligos_path)

    get_stats(contacts_unbinned_path=path_bundle.unbinned_contacts_input,
              sparse_contacts_path=path_bundle.sample_sparse_file_path,
              oligos_path=oligos_path)


def step_two(
        path_bundle: PathBundle,
        oligos_path: str,
        centromeres_coordinates_path: str,
        binning_size_list: List[int],
        aggregate_params_centros: AggregateParams,
        aggregate_params_telos: AggregateParams,
):

    compare_to_wt(path_bundle.global_statistics_input, path_bundle.wt_to_compare_path)

    ponder_mutant(statistics_path=path_bundle.global_statistics_input,
                  contacts_path=path_bundle.unbinned_contacts_input,
                  frequencies_path=path_bundle.unbinned_frequencies_input,
                  binned_type="unbinned")

    for bn in binning_size_list:
        rebin_contacts(contacts_unbinned_path=path_bundle.unbinned_contacts_input,
                       chromosomes_coord_path=centromeres_coordinates_path,
                       bin_size=bn)
        bin_suffix = str(bn // 1000) + "kb"
        binned_contacts_input = \
            os.path.join(path_bundle.sample_output_dir, path_bundle.samp_id + f"_{bin_suffix}_binned_contacts.tsv")
        binned_frequencies_input = \
            os.path.join(path_bundle.sample_output_dir, path_bundle.samp_id + f"_{bin_suffix}_binned_frequencies.tsv")

        ponder_mutant(
            statistics_path=path_bundle.global_statistics_input,
            contacts_path=binned_contacts_input,
            frequencies_path=binned_frequencies_input,
            binned_type=f"{bin_suffix}_binned"
        )

    aggregate(
        binned_contacts_path=os.path.join(
            path_bundle.sample_output_dir,
            path_bundle.samp_id+"_10kb_binned_frequencies.tsv"),
        centros_coord_path=centromeres_coordinates_path,
        oligos_path=oligos_path,
        window_size=aggregate_params_centros.window_size,
        on="centromeres",
        exclude_probe_chr=aggregate_params_centros.excluded_probe_chr,
        excluded_chr_list=aggregate_params_centros.excluded_chr_list,
        inter_normalization=aggregate_params_centros.inter_normalization,
        plot=True)

    aggregate(
        binned_contacts_path=os.path.join(path_bundle.sample_output_dir,
                                          path_bundle.samp_id+"_10kb_binned_frequencies.tsv"),
        centros_coord_path=centromeres_coordinates_path,
        oligos_path=oligos_path,
        window_size=aggregate_params_telos.window_size,
        on="telomeres",
        exclude_probe_chr=aggregate_params_telos.excluded_probe_chr,
        excluded_chr_list=aggregate_params_telos.excluded_chr_list,
        inter_normalization=aggregate_params_telos.inter_normalization,
        plot=True)


if __name__ == "__main__":

    #   Example :

    """
    -s ../data/AD401_AD407_classic
    -f ../data/AD401_AD407_classic/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt
    -c ../data/AD401_AD407_classic/S288c_chr_centro_coordinates.tsv 
    -o ../data/AD401_AD407_classic/capture_oligo_positions.csv
    --samples-vs-wt ../data/AD401_AD407_classic/mutants_vs_ref.csv
    -t 15
    -b 1000 2000 3000 5000 10000 20000 40000 50000 80000 10000
    --window-size-centros 150000  
    --window-size-telos 150000 
    --excluded-chr chr2 chr3 chr5 2_micron mitochondrion, chr_artificial
    --exclude-probe-chr 
    --inter-norm
    """

    parser = argparse.ArgumentParser(
        description="Script that processes sshic samples data.")

    parser.add_argument('-s', '--samples-dir', type=str, required=True,
                        help='Path to inputs directory that contains samples files')

    parser.add_argument('-o', '--oligos-input', type=str, required=True,
                        help='Path to the file that contains positions of oligos')

    parser.add_argument('-f', '--fragments-list-input', type=str, required=True,
                        help='Path to the file fragments_list (hic_stuff output)')

    parser.add_argument('-c', '--centromeres-coordinates-input', type=str, required=True,
                        help='Path to the file centromeres_coordinates')

    parser.add_argument('-b', '--binning-sizes-list', nargs='+', type=int, required=True,
                        help='desired bin size for the rebin step')

    parser.add_argument('-t', '--threads', type=int, required=True,
                        help='desired number of thread to parallelize')

    parser.add_argument('--samples-vs-wt', type=str, required=True,
                        help="Path to the table that associates for each sample a WT dir.")

    parser.add_argument('--window-size-centros', type=int, required=True,
                        help="window (in bp) that defines a focus region to aggregated centromeres")

    parser.add_argument('--window-size-telos', type=int, required=True,
                        help="window (in bp) that defines a focus region to aggregated telomeres")

    parser.add_argument('--excluded-chr', nargs='+', type=str, required=False,
                        help='list of chromosomes to excludes to prevent bias of contacts')

    parser.add_argument('--exclude-probe-chr', action='store_true', required=False,
                        help="exclude the chromosome where the probe comes from (oligo's chromosome)")

    parser.add_argument('--inter-norm', action='store_true', required=False,
                        help="normalize the contacts only on contacts made "
                             "on chromosomes that have not been excluded (inter)")

    args = parser.parse_args()

    if is_debug():
        args.threads = 1

    list_files = [f for f in os.listdir(args.samples_dir) if os.path.isfile(os.path.join(args.samples_dir, f))]
    my_samples = [os.path.join(args.samples_dir, s) for s in list_files if re.match(r'^AD\d+', s)]

    step_one_list = []
    step_two_list = []
    for sample_file in my_samples:
        sample_path_bundle = PathBundle(sample_file, args.samples_vs_wt)

        sample_aggregate_params_centros = AggregateParams(args.window_size_centros, args.exclude_probe_chr,
                                                          args.excluded_chr, args.inter_norm)
        sample_aggregate_params_telos = AggregateParams(args.window_size_telos, args.exclude_probe_chr,
                                                        args.excluded_chr, args.inter_norm)

        step_one_list.append((sample_path_bundle, args.oligos_input, args.fragments_list_input))

        step_two_list.append((sample_path_bundle, args.oligos_input, args.centromeres_coordinates_input,
                              args.binning_sizes_list, sample_aggregate_params_centros, sample_aggregate_params_telos))

    print("step 1 : ")
    for single_sample_data in step_one_list:
        print(single_sample_data[0].samp_id)
        step_one(*single_sample_data)

    print("step 2 : ")
    for single_sample_data in step_two_list:
        print(single_sample_data[0].samp_id)
        step_two(*single_sample_data)
