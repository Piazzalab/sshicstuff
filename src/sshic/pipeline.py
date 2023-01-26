import os
import re
from typing import Optional
from itertools import chain
import multiprocessing as mp
import numpy as np
from sshic import binning, nucleosomes, statistics, filter, format, \
    ponder_mutants, telomeres, centromeres, cohesins


def run_single(
        fragment_list_path: str,
        oligos_positions_path: str,
        probes_to_fragments_path: str,
        centromeres_positions_path: str,
        cohesins_peaks_path: str,
        wt_references_dir: str,
        samples_to_compare_wt: dict,
        nfr_list_path: str,
        outputs_dir: str,
        operations: dict,
        sshic_pcrdupt_dir: str,
):

    #   DEFAULT ARGUMENTS
    bins_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]
    fragments_nfr_filter_list = ['start_only', 'end_only', 'middle', 'start_&_end']
    statistics_cis_region_span = 50000
    centromere_filter_window = 40000
    centromeres_filter_list = ['inner', 'outer', None]
    cohesins_filter_scores_list = [100, 200, 500, 1000, 2000]

    #   OUTPUTS
    hicstuff_dir = outputs_dir + "hicstuff/"
    filter_dir = outputs_dir + "filtered/"
    pondered_dir = outputs_dir + "pondered/"
    binning_dir = outputs_dir + "binned/"
    statistics_dir = outputs_dir + "statistics/"
    nucleosomes_dir = outputs_dir + "nucleosomes/"
    centromeres_dir = outputs_dir + "centromeres/"
    telomeres_dir = outputs_dir + "telomeres/"
    cohesins_dir = outputs_dir + "cohesins/"

    #   FILTER
    if operations['filter'] == 1:
        samples_dir = hicstuff_dir+sshic_pcrdupt_dir
        samples = np.unique(os.listdir(samples_dir))
        if not os.path.exists(filter_dir+sshic_pcrdupt_dir):
            os.makedirs(filter_dir+sshic_pcrdupt_dir)
        for samp in samples:
            samp_id = re.search(r"AD\d+", samp).group()
            filter.run(
                oligos_input_path=oligos_positions_path,
                fragments_input_path=fragment_list_path,
                contacts_input=samples_dir+samp,
                output_path=filter_dir+sshic_pcrdupt_dir+samp_id)

    #   FORMAT
    if operations['format'] == 1:
        format.run(
            fragments_list_path=fragment_list_path,
            oligos_capture_path=oligos_positions_path,
            output_path=probes_to_fragments_path
        )

    #   BINNING
    if operations['binning'] == 1:
        samples_dir = filter_dir+sshic_pcrdupt_dir
        samples = os.listdir(samples_dir)
        for bs in bins_sizes_list:
            print('bin of size: ', bs)

            for samp in samples:
                binning.run(
                    filtered_contacts_path=samples_dir+samp,
                    bin_size=bs,
                    output_dir=binning_dir+sshic_pcrdupt_dir
                )

    if operations['statistics'] == 1:
        not_binned_dir = binning_dir+sshic_pcrdupt_dir+'0kb/'
        sparse_matrix_list = sorted(os.listdir(hicstuff_dir+sshic_pcrdupt_dir))

        samples = [f for f in sorted(os.listdir(not_binned_dir)) if '_contacts' in f]
        samples_id = sorted([re.search(r"AD\d+", f).group() for f in samples])

        for ii_samp, samp in enumerate(samples_id):
            statistics.run(
                cis_range=statistics_cis_region_span,
                sparse_mat_path=hicstuff_dir+sshic_pcrdupt_dir+sparse_matrix_list[ii_samp],
                wt_references_dir=wt_references_dir+sshic_pcrdupt_dir,
                samples_vs_wt=samples_to_compare_wt,
                formatted_contacts_path=not_binned_dir+samples[ii_samp],
                probes_to_fragments_path=probes_to_fragments_path,
                output_dir=statistics_dir+sshic_pcrdupt_dir
            )

    if operations['ponder'] == 1:
        binned_dir_list = os.listdir(binning_dir+sshic_pcrdupt_dir)
        statistics_tables_list = [s for s in sorted(os.listdir(statistics_dir+sshic_pcrdupt_dir)) if 'global' in s]
        samples_id = sorted([re.search(r"AD\d+", f).group() for f in statistics_tables_list])
        for bin_dir in binned_dir_list:
            bin_dir += '/'
            binned_contacts_list = \
                [f for f in sorted(os.listdir(binning_dir+sshic_pcrdupt_dir+bin_dir)) if 'frequencies' in f]

            if not os.path.exists(pondered_dir+sshic_pcrdupt_dir+bin_dir):
                os.makedirs(pondered_dir+sshic_pcrdupt_dir+bin_dir)

            for ii_samp, samp in enumerate(samples_id):
                if samp not in list(chain(*samples_to_compare_wt.values())):
                    continue
                else:
                    binned_contacts_sample = binning_dir+sshic_pcrdupt_dir+bin_dir+binned_contacts_list[ii_samp]
                    stats_table_sample = statistics_dir+sshic_pcrdupt_dir+statistics_tables_list[ii_samp]

                    ponder_mutants.run(
                        samples_vs_wt=samples_to_compare_wt,
                        binned_contacts_path=binned_contacts_sample,
                        statistics_path=stats_table_sample,
                        output_dir=pondered_dir+sshic_pcrdupt_dir+bin_dir,
                    )
        pass

    if operations['nucleosomes'] == 1:
        pass

    if operations['centromeres'] == 1:
        pass

    if operations['telomeres'] == 1:
        pass

    if operations['cohesins'] == 1:
        pass

    pass


