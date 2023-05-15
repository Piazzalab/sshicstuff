import re
import os
import argparse

import pandas as pd

from filter import filter_contacts
from coverage import coverage
from fragments import organize_contacts
from statistics import get_stats
from binning import rebin_contacts
from ponder import ponder_mutant
from aggregated import aggregate
from utils import is_debug


def main(
        my_sample_sparse_file_path: str,
        oligos_path: str,
        fragments_list_path: str,
        centromeres_coordinates_path: str,
        binning_size_list: list,
        df_mutant_vs_ref: pd.DataFrame,
        window_size_centromeres: int,
        window_size_telomeres: int,
        excluded_chr_list: list,
        excluded_probe_chr: bool = True,
        inter_normalization: bool = True,
):
    my_sample_input_dir = os.path.dirname(my_sample_sparse_file_path)
    samp_id = re.match(r'^AD\d+', my_sample_sparse_file_path.split("/")[-1]).group()
    my_sample_output_dir = os.path.join(my_sample_input_dir, samp_id)
    os.makedirs(my_sample_output_dir, exist_ok=True)

    print(samp_id)

    mutants = df_mutant_vs_ref['sample'].values.tolist()
    wt_to_compare_path = None
    if samp_id in mutants:
        ref = df_mutant_vs_ref.loc[df_mutant_vs_ref["sample"] == samp_id]["ref"].tolist()[0]
        wt_to_compare_path = os.path.join(args.wt_dir, ref+".tsv")

    """
    Filtering Sparse matrix
    """
    filtered_contacts_input = os.path.join(my_sample_output_dir, samp_id + "_filtered.tsv")
    if not os.path.exists(filtered_contacts_input):
        filter_contacts(
            oligos_path=oligos_path,
            fragments_path=fragments_list_path,
            contacts_path=my_sample_sparse_file_path,
            output_dir=my_sample_output_dir)

    """
    Computing coverage per digested fragments
    """
    cover = os.path.join(my_sample_output_dir, samp_id + "_coverage_per_fragment.bedgraph")
    if not os.path.exists(cover):
        coverage(
            hic_contacts_path=my_sample_sparse_file_path,
            fragments_path=fragments_list_path,
            output_dir=my_sample_output_dir)

    """
    Tidying contacts between probe and the rest of the genome (unbinned table)
    """
    organize_contacts(
        filtered_contacts_path=filtered_contacts_input,
        oligos_path=oligos_path)
    unbinned_contacts_input = os.path.join(my_sample_output_dir, samp_id+"_unbinned_contacts.tsv")
    unbinned_frequencies_input = os.path.join(my_sample_output_dir, samp_id+"_unbinned_frequencies.tsv")

    """
    Computing some basic statistic about the contacts made
    """
    get_stats(
        contacts_unbinned_path=unbinned_contacts_input,
        sparse_contacts_path=my_sample_sparse_file_path,
        oligos_path=oligos_path,
        wt_reference=wt_to_compare_path,
    )
    global_statistics_input = os.path.join(my_sample_output_dir, samp_id+"_global_statistics.tsv")

    ponder_mutant(
        statistics_path=global_statistics_input,
        contacts_path=unbinned_contacts_input,
        frequencies_path=unbinned_frequencies_input,
        binned_type="unbinned"
    )

    """
    Rebinning the unbinned table at n kb (aggregates contacts on regular range of bp)
    """
    for bn in binning_size_list:
        rebin_contacts(
            contacts_unbinned_path=unbinned_contacts_input,
            chromosomes_coord_path=centromeres_coordinates_path,
            bin_size=bn)

        bin_suffix = str(bn // 1000) + "kb"
        binned_contacts_input = os.path.join(my_sample_output_dir, samp_id + f"_{bin_suffix}_binned_contacts.tsv")
        binned_frequencies_input = os.path.join(my_sample_output_dir, samp_id + f"_{bin_suffix}_binned_frequencies.tsv")
        ponder_mutant(
            statistics_path=global_statistics_input,
            contacts_path=binned_contacts_input,
            frequencies_path=binned_frequencies_input,
            binned_type=f"{bin_suffix}_binned"
        )

    """
    Aggregating contacts around centromeres
    """
    aggregate(
        binned_contacts_path=os.path.join(my_sample_output_dir, samp_id+"_10kb_binned_frequencies.tsv"),
        centros_coord_path=centromeres_coordinates_path,
        oligos_path=oligos_path,
        window_size=window_size_centromeres,
        on="centromeres",
        exclude_probe_chr=excluded_probe_chr,
        excluded_chr_list=excluded_chr_list,
        inter_normalization=inter_normalization,
        plot=True)

    """
    Aggregating contacts around telomeres
    """
    aggregate(
        binned_contacts_path=os.path.join(my_sample_output_dir, samp_id+"_10kb_binned_frequencies.tsv"),
        centros_coord_path=centromeres_coordinates_path,
        oligos_path=oligos_path,
        window_size=window_size_telomeres,
        on="telomeres",
        exclude_probe_chr=excluded_probe_chr,
        excluded_chr_list=excluded_chr_list,
        inter_normalization=inter_normalization,
        plot=True)


if __name__ == "__main__":

    #   Example :

    """
    -s ../test_data/AD401_AD407_classic
    -f ../test_data/AD401_AD407_classic/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt
    -c ../test_data/AD401_AD407_classic/S288c_chr_centro_coordinates.tsv 
    -o ../test_data/AD401_AD407_classic/capture_oligo_positions.csv
    --wt-dir ../test_data/AD401_AD407_classic/capture_efficiencies/
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

    parser.add_argument('--wt-dir', type=str, required=True,
                        help="Path to the wt_capture_efficiency dir.")

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

    df_mutant_vs_wt: pd.DataFrame = pd.read_csv(os.path.join(args.wt_dir, "mutants_vs_ref.csv"), sep="\t")

    if args.threads > 1:
        import multiprocessing as mp
        with mp.Pool(args.threads) as p:
            p.starmap(main, [(
                samp,
                args.oligos_input,
                args.fragments_list_input,
                args.centromeres_coordinates_input,
                args.binning_sizes_list,
                df_mutant_vs_wt,
                args.window_size_centros,
                args.window_size_telos,
                args.excluded_chr,
                args.exclude_probe_chr,
                args.inter_norm) for samp in my_samples])

    else:
        for samp in my_samples:
            main(
                my_sample_sparse_file_path=samp,
                oligos_path=args.oligos_input,
                fragments_list_path=args.fragments_list_input,
                centromeres_coordinates_path=args.centromeres_coordinates_input,
                binning_size_list=args.binning_sizes_list,
                df_mutant_vs_ref=df_mutant_vs_wt,
                excluded_chr_list=args.excluded_chr,
                excluded_probe_chr=args.exclude_probe_chr,
                window_size_centromeres=args.window_size_centros,
                window_size_telomeres=args.window_size_telos
            )

