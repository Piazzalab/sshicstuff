import os
import re
import numpy as np
import multiprocessing as mp

from aggregate import centromeres, telomeres, cohesins


def do_centro(parallel: bool = True):
    span = 150000
    samples = np.unique([re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir+'10kb/')])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(centromeres.run, [(samples_dir + '10kb/' + samp + '_frequencies.tsv',
                                         frag2prob_dir + samp + '_frag_to_prob.tsv',
                                         span,
                                         output_dir,
                                         centromeres_coordinates) for samp in samples])
    else:
        for samp in samples:
            centromeres.run(
                formatted_contacts_path=samples_dir + '10kb/' + samp + '_frequencies.tsv',
                fragments_to_oligos_path=frag2prob_dir + samp + '_frag_to_prob.tsv',
                window_size=span,
                output_path=output_dir,
                centros_coord_path=centromeres_coordinates
            )


def do_telo(parallel: bool = True):
    span = 100000
    samples = np.unique([re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir+'10kb/')])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(telomeres.run, [(samples_dir + '10kb/' + samp + '_frequencies.tsv',
                                       frag2prob_dir + samp + '_frag_to_prob.tsv',
                                       span,
                                       output_dir,
                                       centromeres_coordinates) for samp in samples])
    else:
        for samp in samples:
            telomeres.run(
                formatted_contacts_path=samples_dir + '10kb/' + samp + '_frequencies.tsv',
                fragments_to_oligos_path=frag2prob_dir + samp + '_frag_to_prob.tsv',
                window_size=span,
                output_path=output_dir,
                telomeres_coord_path=centromeres_coordinates
            )


def do_cohesins(parallel: bool = True):
    scores_list = [100, 200, 500, 1000, 2000]
    span = 15000
    cen_filter_window = 40000
    cen_filter_mode = ['inner', 'outer', None]
    cohesins_peaks = "../../data/inputs/HB65_reference_peaks_score50min.bed"
    samples = np.unique([re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir+'10kb/')])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            for m in cen_filter_mode:
                for sc in scores_list:
                    print('score higher than: ', sc)
                    p.starmap(cohesins.run, [(samples_dir + '1kb/' + samp + '_frequencies.tsv',
                                              frag2prob_dir + samp + '_frag_to_prob.tsv',
                                              span,
                                              output_dir,
                                              cohesins_peaks,
                                              centromeres_coordinates,
                                              sc,
                                              cen_filter_window,
                                              m) for samp in samples])
    else:
        for m in cen_filter_mode:
            for sc in scores_list:
                for samp in samples:
                    cohesins.run(
                        formatted_contacts_path=samples_dir + '1kb/' + samp + '_frequencies.tsv',
                        fragments_to_oligos_path=frag2prob_dir + samp + '_frag_to_prob.tsv',
                        window_size=span,
                        output_dir=output_dir,
                        cohesins_peaks_path=cohesins_peaks,
                        centromere_info_path=centromeres_coordinates,
                        score_cutoff=sc,
                        cen_filter_span=cen_filter_window,
                        cen_filter_mode=m)


if __name__ == "__main__":
    samples_dir = "../../data/outputs/binning/sshic_pcrdupkept/"
    frag2prob_dir = "../../data/outputs/formatted/sshic_pcrdupkept/"
    output_dir = "../../data/outputs/aggregated/sshic_pcrdupkept/"
    centromeres_coordinates = "../../data/inputs/S288c_chr_centro_coordinates.tsv"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    modes = ['cohesins']

    if 'centromeres' in modes:
        do_centro()
    if 'telomeres' in modes:
        do_telo()
    if 'cohesins' in modes:
        do_cohesins(parallel=False)

    print('--- DONE ---')
