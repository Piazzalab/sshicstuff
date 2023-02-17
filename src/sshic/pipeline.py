import os
import re
from itertools import chain
import multiprocessing as mp
import numpy as np
from sshic import binning, nucleosomes, statistics, filter, format, \
    ponder_mutants, telomeres, centromeres, cohesins


def run(
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
        parallel: bool = True,
):
    threads = mp.cpu_count()

    #   DEFAULT ARGUMENTS
    bins_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]
    fragments_nfr_filter_list = ['start_only', 'end_only', 'middle', 'start_&_end']

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
    not_binned_dir = binning_dir + sshic_pcrdupt_dir + '0kb/'

    #################################
    #   FILTERING
    #################################
    if operations['filter'] == 1:
        print("Filtered hicstuff sparse matrix by keeping only contacts made by probe")
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
                output_path=filter_dir+sshic_pcrdupt_dir+samp_id
            )

    #################################
    #   FORMATTING
    #################################
    if operations['format'] == 1:
        print("Creating a table to make correspond probe with fragments that contains it")
        format.run(
            fragments_list_path=fragment_list_path,
            oligos_capture_path=oligos_positions_path,
            output_path=probes_to_fragments_path
        )

    #################################
    #   BINNING
    #################################
    if operations['binning'] == 1:
        if not os.path.exists(not_binned_dir):
            print("Organizing the contacts for each probe in the genome ('0kb' binning)")
            samples_dir = filter_dir + sshic_pcrdupt_dir
            samples = os.listdir(samples_dir)
            os.makedirs(not_binned_dir)
            if parallel:
                with mp.Pool(threads) as p:
                    p.starmap(binning.get_fragments_contacts, [(
                        samples_dir + samp,
                        not_binned_dir) for samp in samples]
                        )
            else:
                for samp in samples:
                    binning.get_fragments_contacts(
                        filtered_contacts_path=samples_dir+samp,
                        output_dir=not_binned_dir
                    )

        samples_dir = not_binned_dir
        samples = [f for f in os.listdir(samples_dir) if 'contacts.tsv' in f]
        for bs in bins_sizes_list:
            this_bin_dir = binning_dir + sshic_pcrdupt_dir + str(bs // 1000) + 'kb/'
            print('Rebinning the contacts on bins of size {0} bp'.format(bs))
            if not os.path.exists(this_bin_dir):
                os.makedirs(this_bin_dir)
                if bs == 1000:
                    os.makedirs(this_bin_dir+'probes_centered/')

            if parallel:
                with mp.Pool(threads) as p:
                    p.starmap(binning.rebin_contacts, [(
                        samples_dir+samp,
                        bs,
                        this_bin_dir,
                        probes_to_fragments_path,
                        centromeres_positions_path) for samp in samples]
                    )

            else:
                for samp in samples:
                    binning.rebin_contacts(
                        not_binned_samp_path=samples_dir+samp,
                        bin_size=bs,
                        output_dir=this_bin_dir,
                        probes_to_fragments_path=probes_to_fragments_path,
                        chromosomes_coord_path=centromeres_positions_path
                    )

    #################################
    #   STATISTICS
    #################################
    if operations['statistics'] == 1:
        print('Compute some general statistics on sshic contacts')
        sparse_matrix_list = sorted(os.listdir(hicstuff_dir+sshic_pcrdupt_dir))
        samples = [f for f in sorted(os.listdir(not_binned_dir)) if '_contacts' in f]
        samples_id = sorted([re.search(r"AD\d+", f).group() for f in samples])

        if not os.path.exists(statistics_dir+sshic_pcrdupt_dir):
            os.makedirs(statistics_dir+sshic_pcrdupt_dir)

        if parallel:
            with mp.Pool(threads) as p:
                p.starmap(statistics.run, [(
                    50000,
                    hicstuff_dir+sshic_pcrdupt_dir+sparse_matrix_list[ii_samp],
                    wt_references_dir + sshic_pcrdupt_dir,
                    samples_to_compare_wt,
                    not_binned_dir+samples[ii_samp],
                    probes_to_fragments_path,
                    statistics_dir+sshic_pcrdupt_dir) for ii_samp, samp in enumerate(samples_id)]
                          )
        else:
            for ii_samp, samp in enumerate(samples_id):
                statistics.run(
                    cis_range=50000,
                    sparse_mat_path=hicstuff_dir+sshic_pcrdupt_dir+sparse_matrix_list[ii_samp],
                    wt_references_dir=wt_references_dir+sshic_pcrdupt_dir,
                    samples_vs_wt=samples_to_compare_wt,
                    formatted_contacts_path=not_binned_dir+samples[ii_samp],
                    probes_to_fragments_path=probes_to_fragments_path,
                    output_dir=statistics_dir+sshic_pcrdupt_dir
                )

    #################################
    #   PONDERING
    #################################
    if operations['ponder'] == 1:
        binned_dir_list = os.listdir(binning_dir+sshic_pcrdupt_dir)
        statistics_tables_list = [s for s in sorted(os.listdir(statistics_dir+sshic_pcrdupt_dir)) if 'global' in s]
        samples_id = sorted([re.search(r"AD\d+", f).group() for f in statistics_tables_list])
        for bin_dir in binned_dir_list:
            print('Ponder mutant contacts (rebinned at {0} over WT references)'.format(bin_dir))
            bin_dir += '/'
            binned_contacts_list = \
                [f for f in sorted(os.listdir(binning_dir+sshic_pcrdupt_dir+bin_dir)) if 'frequencies' in f]

            if not os.path.exists(pondered_dir+sshic_pcrdupt_dir+bin_dir):
                os.makedirs(pondered_dir+sshic_pcrdupt_dir+bin_dir)

            if parallel:
                with mp.Pool(threads) as p:
                    p.starmap(ponder_mutants.run, [(
                        samples_to_compare_wt,
                        binning_dir+sshic_pcrdupt_dir+bin_dir+binned_contacts_list[ii_samp],
                        statistics_dir+sshic_pcrdupt_dir+statistics_tables_list[ii_samp],
                        pondered_dir+sshic_pcrdupt_dir+bin_dir) for ii_samp, samp in enumerate(samples_id)]
                              )

            else:
                for ii_samp, samp in enumerate(samples_id):
                    if samp in list(chain(*samples_to_compare_wt.values())):
                        binned_contacts_sample = binning_dir+sshic_pcrdupt_dir+bin_dir+binned_contacts_list[ii_samp]
                        stats_table_sample = statistics_dir+sshic_pcrdupt_dir+statistics_tables_list[ii_samp]
                        ponder_mutants.run(
                            samples_vs_wt=samples_to_compare_wt,
                            binned_contacts_path=binned_contacts_sample,
                            statistics_path=stats_table_sample,
                            output_dir=pondered_dir+sshic_pcrdupt_dir+bin_dir,
                        )

    #################################
    #   NUCLEOSOMES
    #################################
    if operations['nucleosomes'] == 1:
        print('look for fragments inside and outside NFR')

        nfr_in_file = 'fragments_list_in_nfr.tsv'
        nfr_out_file = 'fragments_list_out_nfr.tsv'
        samples = sorted([f for f in os.listdir(not_binned_dir) if 'contacts.tsv' in f])
        for f_filter in fragments_nfr_filter_list:
            print(f_filter)
            filter_dir = f_filter + '/'
            if not os.path.exists(nucleosomes_dir+filter_dir):
                os.makedirs(nucleosomes_dir+filter_dir)
                nucleosomes.preprocess(
                    fragments_list_path=fragment_list_path,
                    fragments_nfr_filter=f_filter,
                    nucleosomes_path=nfr_list_path,
                    output_dir=nucleosomes_dir+filter_dir
                )
            if parallel:
                with mp.Pool(threads) as p:
                    p.starmap(nucleosomes.run, [(
                        not_binned_dir+samp,
                        probes_to_fragments_path,
                        nucleosomes_dir+filter_dir+nfr_in_file,
                        nucleosomes_dir+filter_dir+nfr_out_file,
                        nucleosomes_dir+filter_dir+sshic_pcrdupt_dir) for samp in samples]
                              )
            else:
                for samp in samples:
                    nucleosomes.run(
                        formatted_contacts_path=not_binned_dir+samp,
                        probes_to_fragments_path=probes_to_fragments_path,
                        fragments_in_nfr_path=nucleosomes_dir+filter_dir+nfr_in_file,
                        fragments_out_nfr_path=nucleosomes_dir+filter_dir+nfr_out_file,
                        output_dir=nucleosomes_dir+filter_dir+sshic_pcrdupt_dir
                    )

    #################################
    #   CENTROMERES
    #################################
    if operations['centromeres'] == 1:
        print('aggregated on centromeres positions')
        print('\n')
        print('raw binned tables')
        samples_not_pondered = \
            sorted([f for f in os.listdir(binning_dir+sshic_pcrdupt_dir+'10kb/') if 'contacts.tsv' in f])
        if parallel:
            with mp.Pool(threads) as p:
                p.starmap(centromeres.run, [(
                    binning_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path,
                    centromeres_positions_path,
                    150000,
                    centromeres_dir+'not_pondered/'+sshic_pcrdupt_dir) for samp in samples_not_pondered]
                )
        else:
            for samp in samples_not_pondered:
                centromeres.run(
                    formatted_contacts_path=binning_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path=probes_to_fragments_path,
                    centros_coord_path=centromeres_positions_path,
                    window_size=150000,
                    output_path=centromeres_dir+'not_pondered/'+sshic_pcrdupt_dir
                )
        print('\n')
        print('pondered binned tables')
        samples_pondered = sorted(os.listdir(pondered_dir+sshic_pcrdupt_dir+'10kb/'))
        if parallel:
            with mp.Pool(threads) as p:
                p.starmap(centromeres.run, [(
                    pondered_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path,
                    centromeres_positions_path,
                    150000,
                    centromeres_dir+'pondered/'+sshic_pcrdupt_dir) for samp in samples_pondered]
                )
        else:
            for samp in samples_pondered:
                centromeres.run(
                    formatted_contacts_path=pondered_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path=probes_to_fragments_path,
                    centros_coord_path=centromeres_positions_path,
                    window_size=150000,
                    output_path=centromeres_dir+'pondered/'+sshic_pcrdupt_dir
                )

    #################################
    #   TELOMERES
    #################################
    if operations['telomeres'] == 1:
        print('aggregated on telomeres positions')
        print('\n')
        print('raw binned tables')
        samples_not_pondered = \
            sorted([f for f in os.listdir(binning_dir + sshic_pcrdupt_dir + '10kb/') if 'frequencies.tsv' in f])
        if parallel:
            with mp.Pool(threads) as p:
                p.starmap(telomeres.run, [(
                    binning_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path,
                    centromeres_positions_path,
                    150000,
                    telomeres_dir+'not_pondered/'+sshic_pcrdupt_dir) for samp in samples_not_pondered]
                )
        else:
            for samp in samples_not_pondered:
                telomeres.run(
                    formatted_contacts_path=binning_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path=probes_to_fragments_path,
                    window_size=150000,
                    telomeres_coord_path=centromeres_positions_path,
                    output_path=telomeres_dir+'not_pondered/'+sshic_pcrdupt_dir,
                )
        print('\n')
        print('pondered binned tables')
        samples_pondered = sorted(os.listdir(pondered_dir+sshic_pcrdupt_dir+'10kb/'))
        if parallel:
            with mp.Pool(threads) as p:
                p.starmap(telomeres.run, [(
                    pondered_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path,
                    centromeres_positions_path,
                    150000,
                    telomeres_dir+'pondered/'+sshic_pcrdupt_dir) for samp in samples_pondered]
                )
        else:
            for samp in samples_pondered:
                telomeres.run(
                    formatted_contacts_path=pondered_dir+sshic_pcrdupt_dir+'10kb/'+samp,
                    probes_to_fragments_path=probes_to_fragments_path,
                    window_size=150000,
                    telomeres_coord_path=centromeres_positions_path,
                    output_path=telomeres_dir+'pondered/'+sshic_pcrdupt_dir,
                )

    #################################
    #   COHESINS PEAKS
    #################################
    if operations['cohesins'] == 1:
        cohesins_filter_list = ['inner', 'outer', None]
        cohesins_filter_span = 40000
        cohesins_filter_scores_list = [100, 200, 600, 1000, 2000, 3000]
        print('\n')
        print('raw binned tables')
        samples_not_pondered = \
            sorted([f for f in os.listdir(binning_dir+sshic_pcrdupt_dir+'1kb/') if 'frequencies.tsv' in f])
        for m in cohesins_filter_list:
            if m is not None:
                print('aggregated on cohesins peaks, {1} {0} '
                      'filtered around the chr centromeres'.format(cohesins_filter_span, m))
            else:
                print('aggregated on cohesins peaks')
            for sc in cohesins_filter_scores_list:
                print('peak scores higher than {0}'.format(sc))
                if parallel:
                    with mp.Pool(threads) as p:
                        p.starmap(cohesins.run, [(
                            binning_dir+sshic_pcrdupt_dir+'1kb/'+samp,
                            probes_to_fragments_path,
                            15000,
                            cohesins_peaks_path,
                            centromeres_positions_path,
                            sc,
                            cohesins_filter_span,
                            m,
                            cohesins_dir+'not_pondered/'+sshic_pcrdupt_dir,
                            False) for samp in samples_not_pondered]
                        )
                else:
                    for samp in samples_not_pondered:
                        cohesins.run(
                            formatted_contacts_path=binning_dir+sshic_pcrdupt_dir+'1kb/'+samp,
                            probes_to_fragments_path=probes_to_fragments_path,
                            window_size=15000,
                            cohesins_peaks_path=cohesins_peaks_path,
                            centromere_info_path=centromeres_positions_path,
                            score_cutoff=sc,
                            cen_filter_span=40000,
                            cen_filter_mode=m,
                            output_dir=cohesins_dir+'not_pondered/'+sshic_pcrdupt_dir,
                            plot=False
                        )
        print('\n')
        print('pondered binned tables')
        samples_pondered = sorted(os.listdir(pondered_dir+sshic_pcrdupt_dir+'1kb/'))
        for m in cohesins_filter_list:
            if m is not None:
                print('aggregated on cohesins peaks, {1} {0} '
                      'filtered around the chr centromeres'.format(cohesins_filter_span, m))
            else:
                print('aggregated on cohesins peaks')
            for sc in cohesins_filter_scores_list:
                print('peak scores higher than {0}'.format(sc))
                if parallel:
                    with mp.Pool(threads) as p:
                        p.starmap(cohesins.run, [(
                            pondered_dir+sshic_pcrdupt_dir+'1kb/'+samp,
                            probes_to_fragments_path,
                            15000,
                            cohesins_peaks_path,
                            centromeres_positions_path,
                            sc,
                            cohesins_filter_span,
                            m,
                            cohesins_dir+'pondered/'+sshic_pcrdupt_dir,
                            False) for samp in samples_pondered]
                        )
                else:
                    for samp in samples_pondered:
                        cohesins.run(
                            formatted_contacts_path=pondered_dir+sshic_pcrdupt_dir+'1kb/'+samp,
                            probes_to_fragments_path=probes_to_fragments_path,
                            window_size=15000,
                            cohesins_peaks_path=cohesins_peaks_path,
                            centromere_info_path=centromeres_positions_path,
                            score_cutoff=sc,
                            cen_filter_span=40000,
                            cen_filter_mode=m,
                            output_dir=cohesins_dir+'pondered/'+sshic_pcrdupt_dir,
                            plot=False
                        )
