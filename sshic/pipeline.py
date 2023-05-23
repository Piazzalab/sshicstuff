import re
import os
import argparse
from typing import List, Optional

from filter import filter_contacts
from probe2fragment import associate_probes_to_fragments
from coverage import coverage
from fragments import organize_contacts
from statistics import get_stats, compare_to_wt
from binning import rebin_contacts
from ponder import ponder_mutant
from aggregated import aggregate


class PathBundle:
    def __init__(self, sample_sparse_file_path: str, reference_path):

        ref_name = reference_path.split("/")[-1].split(".")[0]

        self.sample_sparse_file_path = sample_sparse_file_path
        self.samp_id = re.match(r'^AD\d+', sample_sparse_file_path.split("/")[-1]).group()
        parent_dir = os.path.dirname(sample_sparse_file_path)
        self.sample_dir = os.path.join(parent_dir, self.samp_id)
        self.not_pondered_dir = os.path.join(self.sample_dir, "not_pondered")
        self.pondered_dir = os.path.join(self.sample_dir, f"pondered_{ref_name}")
        os.makedirs(self.sample_dir, exist_ok=True)
        os.makedirs(self.pondered_dir, exist_ok=True)
        os.makedirs(self.not_pondered_dir, exist_ok=True)

        self.filtered_contacts_input = os.path.join(self.sample_dir, self.samp_id + "_filtered.tsv")
        self.cover = os.path.join(self.sample_dir, self.samp_id + "_coverage_per_fragment.bedgraph")
        self.unbinned_contacts_input = os.path.join(self.not_pondered_dir, self.samp_id+"_unbinned_contacts.tsv")
        self.unbinned_frequencies_input = os.path.join(self.not_pondered_dir, self.samp_id+"_unbinned_frequencies.tsv")
        self.global_statistics_input = os.path.join(self.sample_dir, self.samp_id+"_global_statistics.tsv")

        if not os.path.exists(reference_path):
            raise ValueError(f"file {reference_path} doesnt exist. "
                             f"Please be sure to first run the pipeline on the reference sample before")
        else:
            self.wt_to_compare_path = reference_path


class AggregateParams:
    def __init__(self, window_size, excluded_probe_chr, excluded_chr_list, inter_normalization):
        self.window_size = window_size
        self.excluded_probe_chr = excluded_probe_chr
        self.excluded_chr_list = excluded_chr_list
        self.inter_normalization = inter_normalization


def check_and_run(output_path, func, *args):
    if not os.path.exists(output_path):
        func(*args)


def do_it(
    path_bundle: PathBundle,
    oligos_path: str,
    fragments_list_path: str,
    centromeres_coordinates_path: str,
    binning_size_list: List[int],
    aggregate_params_centros: AggregateParams,
    aggregate_params_telos: AggregateParams,
    additional_groups: Optional[str] = None
):
    print(f" -- Sample {path_bundle.samp_id} -- \n")
    print(f"Filter contacts \n")
    check_and_run(
        path_bundle.filtered_contacts_input, filter_contacts, oligos_path,
        fragments_list_path, path_bundle.sample_sparse_file_path, path_bundle.sample_dir)

    print(f"Associate the fragment name to probe where it is located \n")
    associate_probes_to_fragments(fragments_list_path, oligos_path)

    print(f"Make the coverage \n")
    check_and_run(
        path_bundle.cover, coverage, path_bundle.sample_sparse_file_path, fragments_list_path, path_bundle.sample_dir)

    print(f"Organize the contacts between probe fragments and the rest of the genome 'unbinned tables' \n")
    organize_contacts(filtered_contacts_path=path_bundle.filtered_contacts_input,
                      oligos_path=oligos_path, chromosomes_coord_path=centromeres_coordinates_path,
                      output_dir=path_bundle.not_pondered_dir, additional_path=additional_groups)

    print(f"Make basic statistics on the contacts (inter/intra chr, cis/trans, ssdna/dsdna etc ...) \n")
    get_stats(
        contacts_unbinned_path=path_bundle.unbinned_contacts_input,
        sparse_contacts_path=path_bundle.sample_sparse_file_path,
        oligos_path=oligos_path, output_dir=path_bundle.sample_dir)

    print(f"Compare the capture efficiency with that of a wild type (may be another sample) \n")
    compare_to_wt(statistics_path=path_bundle.global_statistics_input, reference_path=path_bundle.wt_to_compare_path)

    print(f"Ponder the unbinned contacts and frequencies tables by the efficiency score got on step ahead \n")
    ponder_mutant(
        statistics_path=path_bundle.global_statistics_input, contacts_path=path_bundle.unbinned_contacts_input,
        frequencies_path=path_bundle.unbinned_frequencies_input, binned_type="unbinned",
        output_dir=path_bundle.pondered_dir, additional_path=additional_groups)

    print(f"Rebin and ponder the unbinned tables (contacts and frequencies) at : \n")
    for bn in binning_size_list:
        bin_suffix = str(bn // 1000) + "kb"
        print(bin_suffix)
        rebin_contacts(contacts_unbinned_path=path_bundle.unbinned_contacts_input,
                       chromosomes_coord_path=centromeres_coordinates_path, bin_size=bn,
                       output_dir=path_bundle.not_pondered_dir)

        binned_contacts_input = \
            os.path.join(path_bundle.not_pondered_dir, path_bundle.samp_id + f"_{bin_suffix}_binned_contacts.tsv")
        binned_frequencies_input = \
            os.path.join(path_bundle.not_pondered_dir, path_bundle.samp_id + f"_{bin_suffix}_binned_frequencies.tsv")

        ponder_mutant(
            statistics_path=path_bundle.global_statistics_input, contacts_path=binned_contacts_input,
            frequencies_path=binned_frequencies_input, binned_type=f"{bin_suffix}_binned",
            output_dir=path_bundle.pondered_dir, additional_path=additional_groups)
    print("\n")

    print("Make an aggregated of contacts around centromeres \n")
    aggregate(
        binned_contacts_path=os.path.join(path_bundle.not_pondered_dir,
                                          path_bundle.samp_id+"_10kb_binned_frequencies.tsv"),
        centros_coord_path=centromeres_coordinates_path, oligos_path=oligos_path,
        window_size=aggregate_params_centros.window_size, on="centromeres",
        output_dir=path_bundle.sample_dir,
        exclude_probe_chr=aggregate_params_centros.excluded_probe_chr,
        excluded_chr_list=aggregate_params_centros.excluded_chr_list,
        inter_normalization=aggregate_params_centros.inter_normalization, plot=True)

    print("Make an aggregated of contacts around telomeres \n")
    aggregate(
        binned_contacts_path=os.path.join(path_bundle.not_pondered_dir,
                                          path_bundle.samp_id+"_10kb_binned_frequencies.tsv"),
        centros_coord_path=centromeres_coordinates_path, oligos_path=oligos_path,
        window_size=aggregate_params_telos.window_size, on="telomeres",
        output_dir=path_bundle.sample_dir,
        exclude_probe_chr=aggregate_params_telos.excluded_probe_chr,
        excluded_chr_list=aggregate_params_telos.excluded_chr_list,
        inter_normalization=aggregate_params_telos.inter_normalization, plot=True)

    print(f"--- {path_bundle.samp_id} DONE --- \n\n")


if __name__ == "__main__":

    #   Command to enter for parameters (parse)
    """
    -s ../../data/samples/AD162_AD407/AD404_S288c_DSB_LY_Capture_artificial_cutsite_PCRfree_q20.txt
    -f ../../data/samples/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt
    -c ../../data/samples/inputs/S288c_chr_centro_coordinates.tsv 
    -o ../../data/samples/inputs/capture_oligo_positions.csv
    -r ../../data/samples/inputs/refs/ref_WT4h_v2.tsv
    -a ../../data/samples/inputs/additional_probe_groups.tsv
    -b 1000 2000 3000 5000 10000 20000 40000 50000 80000 10000
    --window-size-centros 150000  
    --window-size-telos 150000 
    --excluded-chr chr2 chr3 chr5 2_micron mitochondrion, chr_artificial
    --exclude-probe-chr 
    --inter-norm
    """

    parser = argparse.ArgumentParser(
        description="Script that processes sshic samples data.")

    parser.add_argument('-s', '--sparse', type=str, required=True,
                        help='Path hicstuff sparse matrix output ')

    parser.add_argument('-o', '--oligos-input', type=str, required=True,
                        help='Path to the file that contains positions of oligos')

    parser.add_argument('-f', '--fragments-list', type=str, required=True,
                        help='Path to the file fragments_list (hic_stuff output)')

    parser.add_argument('-c', '--centromeres-coordinates-input', type=str, required=True,
                        help='Path to the file centromeres_coordinates')

    parser.add_argument('-b', '--binning-sizes-list', nargs='+', type=int, required=True,
                        help='desired bin size for the rebin step')

    parser.add_argument('-r', '--reference', type=str, required=False,
                        help="Path to the reference WT to ponder the sample.")

    parser.add_argument('-a', '--additional', type=str, required=False,
                        help='Path to additional groups of probes table')

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

    sample_path_bundle = PathBundle(args.sparse, args.reference)
    sample_aggregate_params_centros = AggregateParams(
        args.window_size_centros, args.exclude_probe_chr, args.excluded_chr, args.inter_norm)
    sample_aggregate_params_telos = AggregateParams(
        args.window_size_telos, args.exclude_probe_chr, args.excluded_chr, args.inter_norm)

    sample_data = [
        sample_path_bundle, args.oligos_input, args.fragments_list, args.centromeres_coordinates_input,
        args.binning_sizes_list, sample_aggregate_params_centros, sample_aggregate_params_telos, args.additional]
    do_it(*sample_data)
