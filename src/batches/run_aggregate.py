import os
import re
import numpy as np
import multiprocessing as mp
import utils.tools as tools
from aggregate import centromeres, telomeres, cohesins


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
            for sc in scores_list:
                for samp in samples:
                    cohesins.run(
                        formatted_contacts_path=samples_dir + '1kb/' + samp + '_frequencies.tsv',
                        probes_to_fragments_path=probes2frag,
                        window_size=span,
                        output_dir=output_dir,
                        cohesins_peaks_path=cohesins_peaks,
                        centromere_info_path=centromeres_coordinates,
                        score_cutoff=sc,
                        cen_filter_span=cen_filter_window,
                        cen_filter_mode=m)


if __name__ == "__main__":

    centromeres_positions = "../../data/inputs/S288c_chr_centro_coordinates.tsv"
    probes_and_fragments = "../../data/inputs/probes_to_fragments.tsv"
    cohesins_peaks_bed = "../../data/inputs/HB65_reference_peaks_score50min.bed"
    scores_list = [100, 200, 500, 1000, 2000]
    cen_filter_window = 40000
    cen_filter_modes = ['inner', 'outer', None]

    parallel_state: bool = True
    if tools.is_debug():
        parallel_state = False

    sshic_dir = ['sshic/', 'sshic_pcrdupkept/']
    modes = ['cohesins']

    for hicd in sshic_dir:
        print(hicd)

        binned_contacts = "../../data/outputs/binning/" + hicd
        aggregated_output_dir = "../../data/outputs/aggregated/" + hicd

        if not os.path.exists(aggregated_output_dir):
            os.makedirs(aggregated_output_dir)

        if 'centromeres' in modes:
            print('Centromeres')
            do_centro(
                centromeres_coordinates=centromeres_positions,
                probes2frag=probes_and_fragments,
                samples_dir=binned_contacts,
                span=150000,
                output_dir=aggregated_output_dir,
                parallel=parallel_state
            )

        if 'telomeres' in modes:
            print('Telomeres')
            do_telo(
                centromeres_coordinates=centromeres_positions,
                probes2frag=probes_and_fragments,
                samples_dir=binned_contacts,
                span=100000,
                output_dir=aggregated_output_dir,
                parallel=parallel_state
            )

        if 'cohesins' in modes:
            print('Cohesins Peaks')
            do_cohesins(
                samples_dir=binned_contacts,
                centromeres_coordinates=centromeres_positions,
                probes2frag=probes_and_fragments,
                cohesins_peaks=cohesins_peaks_bed,
                output_dir=aggregated_output_dir,
                span=15000,
                scores=scores_list,
                cen_filter_operations=cen_filter_modes,
                cen_filter_span=cen_filter_window,
                parallel=parallel_state
            )

    print('--- DONE ---')
