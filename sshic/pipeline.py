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


def main(
        inputs_dir_path: str,
        oligos_path: str,
        fragments_list_path: str,
        centromeres_coordinates_path: str,
        binning_size_list: list
):

    my_samples = list_folders(inputs_dir_path)
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

if __name__ == "__main__":

    # data_dir = "../test_data"
    # oligos_input = os.path.join(data_dir, "capture_oligo_positions.csv")
    # fragments_list_input = os.path.join(data_dir, "fragments_list.txt")
    # centromeres_coordinates_input = os.path.join(data_dir, "S288c_chr_centro_coordinates.tsv")
    # binning_sizes_list = [1000, 5000, 10000, 50000, 100000]

    parser = argparse.ArgumentParser(
        description="Script that processes sshic samples data.")

    parser.add_argument('--inputs-dir', type=str, required=True,
                        help='Path to inputs directory that contains samples and inputs files')
    parser.add_argument('--oligos-input', type=str, required=True,
                        help='name of the file that contains positions of oligos')
    parser.add_argument('--fragments-list-input', type=str, required=True,
                        help='name of the fragments_list file (hic_stuff output)')
    parser.add_argument('--centromeres-coordinates-input', type=str, required=True,
                        help='name of the centromeres_coordinates file')
    parser.add_argument('--binning-sizes-list', nargs='+', type=int, required=True,
                        help='desired bin size for the rebin step')

    args = parser.parse_args()
    main(
        args.inputs_dir,
        args.oligos_input,
        args.fragments_list_input,
        args.centromeres_coordinates_input,
        args.binning_sizes_list)

