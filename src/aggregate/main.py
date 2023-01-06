import os
import re

from aggregate import centromeres, telomeres, cohesins


def do_centro():
    span = 150000
    centromeres_coordinates = "../../data/inputs/S288c_chr_centro_coordinates.tsv"
    samples = os.listdir(samples_dir + '10kb/')
    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        centromeres.run(
            formatted_contacts_path=samples_dir+'10kb/'+samp,
            window_size=span,
            output_path=output_dir,
            sample_name=samp_id,
            centros_coord_path=centromeres_coordinates
        )


def do_telo():
    span = 100000
    telomeres_coordinates = "../../data/inputs/S288c_chr_centro_coordinates.tsv"
    samples = os.listdir(samples_dir + '10kb/')
    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        telomeres.run(
            formatted_contacts_path=samples_dir+'10kb/'+samp,
            window_size=span,
            output_path=output_dir,
            sample_name=samp_id,
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
    samples_dir = "../../data/outputs/binning/ssHiC_filtered/"
    output_dir = "../../data/outputs/aggregated/ssHiC_filtered/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    modes = ['cohesins']

    if 'centromeres' in modes:
        do_centro()
    if 'telomeres' in modes:
        do_telo()
    if 'cohesins' in modes:
        do_cohesins()

    print('--- DONE ---')
