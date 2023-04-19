import os
import re
import sys
import getopt

from universal.utils import list_folders
from universal.filter import filter_contacts
from universal.coverage import coverage
from universal.probe2fragment import associate_probes_to_fragments as p2f
from universal.fragments import organize_contacts
from universal.statistics import get_stats
from universal.binning import rebin_contacts


if __name__ == "__main__":

    data_dir = "../../test_data/"
    oligos_input = data_dir + "capture_oligo_positions.csv"
    fragments_list_input = data_dir + "fragments_list.txt"
    centromeres_coordinates_input = data_dir + "S288c_chr_centro_coordinates.tsv"

    binning_sizes_list = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
    my_samples = list_folders(data_dir)
    for sample_dir in my_samples:
        sample_id = re.search(r"AD\d+", sample_dir).group()
        sparse_contacts_input = sample_dir + "AD162_S288c_DSB_LY_Capture_artificial_cutsite_q30.txt"

        filter_contacts(
            oligos_path=oligos_input,
            fragments_path=fragments_list_input,
            contacts_path=sparse_contacts_input
        )
        filtered_contacts_input = sample_dir+sample_id+"_filtered.tsv"

        coverage(
            hic_contacts_path=sparse_contacts_input,
            fragments_path=fragments_list_input
        )

        p2f(
            fragments_list_path=fragments_list_input,
            oligos_capture_path=oligos_input
        )

        organize_contacts(
            filtered_contacts_path=filtered_contacts_input
        )
        unbinned_contacts_input = sample_dir+sample_id+"_unbinned_contacts.tsv"
        unbinned_frequencies_input = sample_dir+sample_id+"_unbinned_frequencies.tsv"

        get_stats(
            contacts_unbinned_path=unbinned_contacts_input,
            sparse_contacts_path=sparse_contacts_input
        )

        for bn in binning_sizes_list:
            rebin_contacts(
                contacts_unbinned_path=sample_dir+sample_id+"_unbinned_contacts.tsv",
                chromosomes_coord_path=centromeres_coordinates_input,
                bin_size=bn
            )
