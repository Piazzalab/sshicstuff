import os
import re
import multiprocessing as mp

from aggregate import centromeres, telomeres, cohesins


def do_centro(parallel: bool = True):
    span = 150000
    centromeres_coordinates = "../../data/inputs/S288c_chr_centro_coordinates.tsv"
    samples = os.listdir(samples_dir + '10kb/')

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(centromeres.run, [(samples_dir+'10kb/'+samp,
                                         span,
                                         output_dir,
                                         centromeres_coordinates) for samp in samples])
    else:
        for samp in samples:
            centromeres.run(
                formatted_contacts_path=samples_dir+'10kb/'+samp,
                window_size=span,
                output_path=output_dir,
                centros_coord_path=centromeres_coordinates
            )


def do_telo(parallel: bool = True):
    span = 100000
    telomeres_coordinates = "../../data/inputs/S288c_chr_centro_coordinates.tsv"
    samples = os.listdir(samples_dir + '10kb/')

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(telomeres.run, [(samples_dir+'10kb/'+samp,
                                       span,
                                       output_dir,
                                       telomeres_coordinates) for samp in samples])
    else:
        for samp in samples:
            telomeres.run(
                formatted_contacts_path=samples_dir+'10kb/'+samp,
                window_size=span,
                output_path=output_dir,
                telomeres_coord_path=telomeres_coordinates
            )


def do_cohesins():
    scores_list = [50, 100, 200, 600, 1000, 2000]
    for sc in scores_list:
        span = 15000
        cohesins_peaks = "../../data/inputs/HB65_reference_peaks_score50min.bed"
        samples = os.listdir(samples_dir + '1kb/')
        for samp in samples:
            samp_id = re.search(r"AD\d+", samp).group()
            cohesins.run(
                formatted_contacts_path=samples_dir+'1kb/'+samp,
                window_size=span,
                output_path=output_dir,
                sample_name=samp_id,
                cohesins_peaks_path=cohesins_peaks,
                score_cutoff=sc)


if __name__ == "__main__":
    samples_dir = "../../data/outputs/binning/sshic/"
    output_dir = "../../data/outputs/aggregated/sshic/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    modes = ['telomeres']

    if 'centromeres' in modes:
        do_centro()
    if 'telomeres' in modes:
        do_telo()
    if 'cohesins' in modes:
        do_cohesins()

    print('--- DONE ---')
