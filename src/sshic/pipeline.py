import os
import re
from typing import Optional
from itertools import chain
import multiprocessing as mp
import numpy as np
from sshic import binning, nucleosomes, statistics, filter, format, \
    ponder_mutants, telomeres, centromeres, cohesins


def do_filter(
        fragments: str,
        oligos: str,
        samples_dir: str,
        output_dir: str):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique(os.listdir(samples_dir))
    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        filter.run(
            oligos_input_path=oligos,
            fragments_input_path=fragments,
            contacts_input=samples_dir+samp,
            output_path=output_dir+samp_id)


def do_format(
        fragments: str,
        oligos: str,
        probes2frag: str,
        samples_dir: str,
        output_dir: str,
        parallel: bool = True):

    if not os.path.exists(probes2frag):
        format.fragments_to_oligos(
            fragments_list_path=fragments,
            oligos_capture_path=oligos,
            output_path=probes2frag
        )

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique(os.listdir(samples_dir))

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(format.run, [(samples_dir + samp,
                                    output_dir) for samp in samples])

    else:
        for samp in samples:
            format.run(
                filtered_contacts_path=samples_dir+samp,
                output_dir=output_dir
            )


def do_binning(
        artificial_genome: str,
        samples_dir: str,
        output_dir: str,
        parallel: bool = True):
    bin_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique(
        [f for f in os.listdir(samples_dir) if '_frequencies.tsv' in f])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            for bs in bin_sizes_list:
                print('bin of size: ', bs)

                new_output_dir = output_dir + str(bs // 1000) + 'kb/'
                if not os.path.exists(new_output_dir):
                    os.makedirs(new_output_dir)

                p.starmap(binning.run, [(artificial_genome,
                                         samples_dir + samp,
                                         bs,
                                         new_output_dir) for samp in samples])
    else:
        for bs in bin_sizes_list:
            print('bin of size: ', bs)
            new_output_dir = output_dir + str(bs // 1000) + 'kb/'
            if not os.path.exists(new_output_dir):
                os.makedirs(new_output_dir)
            for samp in samples:
                binning.run(
                    artificial_genome_path=artificial_genome,
                    formatted_contacts_path=samples_dir+samp,
                    bin_size=bs,
                    output_dir=new_output_dir
                )


def do_stats(
        hicstuff_dir: str,
        samples_dir: str,
        probes2frag: str,
        wt_references: Optional[str],
        samples_vs_wt: Optional[dict],
        output_dir: str,
        cis_span: int,
        parallel: bool = True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sparse_mat_list = sorted(os.listdir(hicstuff_dir))
    formatted_contacts_list = [f for f in sorted(os.listdir(samples_dir)) if '_contacts' in f]
    samples_id = sorted([re.search(r"AD\d+", f).group() for f in sparse_mat_list])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(statistics.run, [(
                cis_span,
                hicstuff_dir + sparse_mat_list[ii_samp],
                wt_references,
                samples_vs_wt,
                samples_dir+formatted_contacts_list[ii_samp],
                probes2frag,
                output_dir) for ii_samp, samp in enumerate(samples_id)])

    else:
        for ii_samp, samp in enumerate(samples_id):
            statistics.run(
                cis_range=cis_span,
                sparse_mat_path=hicstuff_dir+sparse_mat_list[ii_samp],
                wt_references_dir=wt_references,
                samples_vs_wt=samples_vs_wt,
                formatted_contacts_path=samples_dir+formatted_contacts_list[ii_samp],
                probes_to_fragments_path=probes2frag,
                output_dir=output_dir
            )


def do_ponder(
        samples_vs_wt: dict,
        binned_contacts_dir: str,
        statistics_contacts_dir: str,
        output_dir: str,
        parallel: bool = True):

    bins_dir_list = os.listdir(binned_contacts_dir)
    statistics_contacts_list = [s for s in sorted(os.listdir(statistics_contacts_dir)) if 'global' in s]
    samples_id = sorted([re.search(r"AD\d+", f).group() for f in statistics_contacts_list])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            for bin_dir in bins_dir_list:
                print('bin of size: ', bin_dir)
                bin_dir += '/'
                binned_contacts_list = \
                    [f for f in sorted(os.listdir(binned_contacts_dir + bin_dir)) if 'frequencies' in f]

                new_output_dir = output_dir + bin_dir
                if not os.path.exists(new_output_dir):
                    os.makedirs(new_output_dir)

                p.starmap(ponder_mutants.run, [
                    (samples_vs_wt,
                     binned_contacts_dir+bin_dir+binned_contacts_list[ii_samp],
                     statistics_contacts_dir+statistics_contacts_list[ii_samp],
                     new_output_dir) for ii_samp, samp in enumerate(samples_id)])

    else:
        for bin_dir in bins_dir_list:
            bin_dir += '/'
            binned_contacts_list = \
                [f for f in sorted(os.listdir(binned_contacts_dir + bin_dir)) if 'frequencies' in f]

            new_output_dir = output_dir + bin_dir
            if not os.path.exists(new_output_dir):
                os.makedirs(new_output_dir)

            for ii_samp, samp in enumerate(samples_id):
                if samp not in list(chain(*samples_vs_wt.values())):
                    continue
                else:
                    binned_contacts = binned_contacts_list[ii_samp]
                    stats_contacts = statistics_contacts_list[ii_samp]

                    ponder_mutants.run(
                        samples_vs_wt=samples_vs_wt,
                        binned_contacts_path=binned_contacts_dir+bin_dir+binned_contacts,
                        statistics_path=statistics_contacts_dir+stats_contacts,
                        output_dir=new_output_dir,
                    )


def do_nucleo(
        samples_dir: str,
        fragments: str,
        probe2frag: str,
        nucleosomes_path: str,
        output_dir: str,
        parallel: bool = True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_parent_dir = os.path.dirname(os.path.dirname(output_dir))+'/'
    files = os.listdir(output_parent_dir)
    nfr_in = 'fragments_list_in_nfr.tsv'
    nfr_out = 'fragments_list_out_nfr.tsv'

    if np.sum(np.isin([nfr_in, nfr_out], files)) != 2:
        nucleosomes.preprocess(
            fragments_list_path=fragments,
            nucleosomes_path=nucleosomes_path,
            output_dir=output_parent_dir
        )

    samples = np.unique(
        [re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir)])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(nucleosomes.run, [(samples_dir + samp + '_contacts.tsv',
                                         probe2frag,
                                         output_parent_dir + nfr_in,
                                         output_parent_dir + nfr_out,
                                         output_dir) for samp in samples])

    else:
        for samp in samples:
            nucleosomes.run(
                formatted_contacts_path=samples_dir+samp+'_contacts.tsv',
                probes_to_fragments_path=probe2frag,
                fragments_in_nfr_path=output_parent_dir+nfr_in,
                fragments_out_nfr_path=output_parent_dir+nfr_out,
                output_dir=output_dir
            )


def do_centro(
        centromeres_coordinates: str,
        probes2frag: str,
        samples_dir: str,
        output_dir: str,
        span: int,
        parallel: bool = True):

    samples = np.unique(
        [re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir+'10kb/')])
    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(centromeres.run, [(samples_dir + '10kb/' + samp + '_frequencies.tsv',
                                         probes2frag,
                                         span,
                                         output_dir,
                                         centromeres_coordinates) for samp in samples])
    else:
        for samp in samples:
            centromeres.run(
                formatted_contacts_path=samples_dir + '10kb/' + samp + '_frequencies.tsv',
                probes_to_fragments_path=probes2frag,
                window_size=span,
                output_path=output_dir,
                centros_coord_path=centromeres_coordinates
            )


def do_telo(
        centromeres_coordinates: str,
        probes2frag: str,
        samples_dir: str,
        output_dir: str,
        span: int,
        parallel: bool = True):

    samples = np.unique(
        [re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir+'10kb/')])
    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(telomeres.run, [(samples_dir + '10kb/' + samp + '_frequencies.tsv',
                                       probes2frag,
                                       span,
                                       output_dir,
                                       centromeres_coordinates) for samp in samples])
    else:
        for samp in samples:
            telomeres.run(
                formatted_contacts_path=samples_dir + '10kb/' + samp + '_frequencies.tsv',
                probes_to_fragments_path=probes2frag,
                window_size=span,
                output_path=output_dir,
                telomeres_coord_path=centromeres_coordinates
            )


def do_cohesins(
        samples_dir: str,
        centromeres_coordinates: str,
        probes2frag: str,
        cohesins_peaks: str,
        output_dir: str,
        span: int,
        cen_filter_operations: list[str | None],
        cen_filter_span: int,
        scores: list[int],
        parallel: bool = True):

    samples = np.unique(
        [re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir+'1kb/')])
    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            for m in cen_filter_operations:
                for sc in scores:
                    print('score higher than: ', sc)
                    p.starmap(cohesins.run, [(samples_dir + '1kb/' + samp + '_frequencies.tsv',
                                              probes2frag,
                                              span,
                                              output_dir,
                                              cohesins_peaks,
                                              centromeres_coordinates,
                                              sc,
                                              cen_filter_span,
                                              m) for samp in samples])
    else:
        for m in cen_filter_operations:
            for sc in scores:
                for samp in samples:
                    cohesins.run(
                        formatted_contacts_path=samples_dir + '1kb/' + samp + '_frequencies.tsv',
                        probes_to_fragments_path=probes2frag,
                        window_size=span,
                        output_dir=output_dir,
                        cohesins_peaks_path=cohesins_peaks,
                        centromere_info_path=centromeres_coordinates,
                        score_cutoff=sc,
                        cen_filter_span=cen_filter_span,
                        cen_filter_mode=m)
