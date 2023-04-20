import re
import os
import argparse

from utils import list_folders
from filter import filter_contacts
from coverage import coverage
from probe2fragment import associate_probes_to_fragments as p2f
from fragments import organize_contacts
from statistics import get_stats
from binning import rebin_contacts
from aggregated import aggregate


def main(
        inputs_dir_path: str,
        samples_dir_path: str,
        oligos_path: str,
        fragments_list_path: str,
        centromeres_coordinates_path: str,
        binning_size_list: list,
        window_size_centromeres: int,
        window_size_telomeres: int,
        excluded_chr_list: list,
        excluded_probe_chr: bool = True,
        inter_normalization: bool = True
):

    listed_dirs = list_folders(inputs_dir_path)
    my_samples = [d for d in listed_dirs if re.match(r'^AD\d+$', d)]
    my_samples_absolute = [os.path.join(inputs_dir_path, samp) for samp in my_samples]
    for sample_dir in my_samples_absolute:
        sample_id = re.search(r"AD\d+", sample_dir).group()
        sparse_contacts_input = os.path.join(sample_dir, sample_id+"_S288c_DSB_LY_Capture_artificial_cutsite_q30.txt")

        filter_contacts(
            oligos_path=oligos_path,
            fragments_path=fragments_list_path,
            contacts_path=sparse_contacts_input)
        filtered_contacts_input = os.path.join(sample_dir, sample_id+"_filtered.tsv")

        coverage(hic_contacts_path=sparse_contacts_input, fragments_path=fragments_list_path)
        p2f(fragments_list_path=fragments_list_path, oligos_capture_path=oligos_path)
        organize_contacts(filtered_contacts_path=filtered_contacts_input)
        unbinned_contacts_input = os.path.join(sample_dir, sample_id+"_unbinned_contacts.tsv")
        unbinned_frequencies_input = os.path.join(sample_dir, sample_id+"_unbinned_frequencies.tsv")

        get_stats(contacts_unbinned_path=unbinned_contacts_input, sparse_contacts_path=sparse_contacts_input)

        for bn in binning_size_list:
            rebin_contacts(
                contacts_unbinned_path=unbinned_contacts_input,
                chromosomes_coord_path=centromeres_coordinates_path,
                bin_size=bn)

        aggregate(
            binned_contacts_path=os.path.join(sample_dir, sample_id+"_1kb_binned_frequencies.tsv"),
            centros_coord_path=centromeres_coordinates_path,
            window_size=window_size_centromeres,
            on="centromeres",
            exclude_probe_chr=excluded_probe_chr,
            excluded_chr_list=excluded_chr_list,
            inter_normalization=inter_normalization,
            plot=True)

        aggregate(
            binned_contacts_path=os.path.join(sample_dir, sample_id+"_1kb_binned_frequencies.tsv"),
            centros_coord_path=centromeres_coordinates_path,
            window_size=window_size_telomeres,
            on="telomeres",
            exclude_probe_chr=excluded_probe_chr,
            excluded_chr_list=excluded_chr_list,
            inter_normalization=inter_normalization,
            plot=True)
            

if __name__ == "__main__":
    """
    -i ../test_data/ 
    -o ../test_data/outputs
    -s ../test_data/S288c_DSB_chrIV845464_Capture_APO1345
    -oligos-input ../test_data/capture_oligo_positions.csv
    -f ../test_data/fragments_list.txt
    -c ../test_data/S288c_chr_centro_coordinates.tsv 
    -b 1000 10000 10000
    --window-size-centros 150000  
    --window-size-telos 150000 
    --excluded-chr chr2 chr3 chr5 2_micron mitochondrion, chr_artificial
    --exclude-probe-chr 
    --inter-norm
    """

    parser = argparse.ArgumentParser(
        description="Script that processes sshic samples data.")

    parser.add_argument('-i', '--inputs-dir', type=str, required=True,
                        help='Path to inputs directory that contains inputs files')
    parser.add_argument('-o', '--outputs-dir', type=str, required=True,
                        help='Path to outputs directory that contains outputs files for each sample')
    parser.add_argument('-s', '--samples-dir', type=str, required=True,
                        help='Path to inputs directory that contains samples files')
    parser.add_argument('--oligos-input', type=str, required=True,
                        help='name of the file that contains positions of oligos')
    parser.add_argument('-f', '--fragments-list-input', type=str, required=True,
                        help='name of the fragments_list file (hic_stuff output)')
    parser.add_argument('-c', '--centromeres-coordinates-input', type=str, required=True,
                        help='name of the centromeres_coordinates file')
    parser.add_argument('-b', '--binning-sizes-list', nargs='+', type=int, required=True,
                        help='desired bin size for the rebin step')
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
    main(
        inputs_dir_path=args.inputs_dir,
        samples_dir_path=args.outputs_dir,
        oligos_path=args.oligos_input,
        fragments_list_path=args.fragments_list_input,
        centromeres_coordinates_path=args.centromeres_coordinates_input,
        binning_size_list=args.binning_sizes_list,
        excluded_chr_list=args.excluded_chr,
        excluded_probe_chr=args.exclude_probe_chr,
        window_size_centromeres=args.window_size_centros,
        window_size_telomeres=args.window_size_telos
    )

