import os
import re

from contacts import binning, statistics


def do_binning():
    bin_sizes_list = [0, 1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]
    artificial_genome = "../../data/inputs/S288c_DSB_LY_capture_artificial.fa"
    samples_dir = "../../data/outputs/filter/ssHiC_filtered/"

    for bs in bin_sizes_list:

        output_dir = "../../data/outputs/binning/ssHiC_filtered/" + str(bs // 1000) + 'kb/'

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        samples = os.listdir(samples_dir)
        for samp in samples:
            samp_id = re.search(r"AD\d+", samp).group()
            binning.run(
                artificial_genome_path=artificial_genome,
                filtered_contacts_path=samples_dir + samp,
                bin_size=bs,
                output_path=output_dir + samp_id
            )


def do_stats():
    cis_range = 50000
    samples_dir = "../../data/outputs/binning/ssHiC_filtered_PCRduplicateskept/0kb/"
    output_dir = "../../data/outputs/statistics/ssHiC_filtered_PCRduplicateskept/"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    samples = os.listdir(samples_dir)
    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        statistics.run(
            cis_range=cis_range,
            binned_contacts_path=samples_dir + samp,
            output_path=output_dir + samp_id
        )


if __name__ == "__main__":
    modes = ['binning', 'statistics']
    if 'binning' in modes:
        do_binning()
    if 'statistics' in modes:
        do_stats()

    print('--- DONE ---')
