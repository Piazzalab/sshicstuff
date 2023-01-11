import os
import re
import multiprocessing as mp
import numpy as np
from contacts import filter, format, binning, statistics


def do_filter():
    fragments = "../../data/inputs/fragments_list.txt"
    oligos = "../../data/inputs/capture_oligo_positions.csv"
    samples_dir = "../../data/outputs/hicstuff/sshic/"
    output_dir = "../../data/outputs/filtered/sshic/"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = os.listdir(samples_dir)
    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        filter.run(
            oligos_input_path=oligos,
            fragments_input_path=fragments,
            contacts_input=samples_dir+samp,
            output_path=output_dir+samp_id)


def do_format(parallel: bool = True):
    samples_dir = "../../data/outputs/filtered/sshic/"
    output_dir = "../../data/outputs/formatted/sshic/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = os.listdir(samples_dir)

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(format.run, [(samples_dir+samp,
                                    output_dir) for samp in samples])

    else:
        for samp in samples:
            format.run(
                filtered_contacts_path=samples_dir+samp,
                output_dir=output_dir
            )


def do_binning(parallel: bool = True):
    bin_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000, 10**9]

    artificial_genome = "../../data/inputs/S288c_DSB_LY_capture_artificial.fa"
    samples_dir = "../../data/outputs/formatted/sshic/"
    output_dir = "../../data/outputs/binning/sshic/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique([f for f in os.listdir(samples_dir) if '_frequencies.tsv' in f])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            for bs in bin_sizes_list:
                print('bin of size: ', bs)

                output_dir += str(bs // 1000) + 'kb/'
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                p.starmap(binning.run, [(artificial_genome,
                                         samples_dir+samp,
                                         bs,
                                         output_dir) for samp in samples])
    else:
        for bs in bin_sizes_list:
            print('bin of size: ', bs)
            output_dir += str(bs // 1000) + 'kb/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            for samp in samples:
                binning.run(
                    artificial_genome_path=artificial_genome,
                    formatted_contacts_path=samples_dir+samp,
                    bin_size=bs,
                    output_dir=output_dir
                )


def do_stats(parallel: bool = True):
    cis_range = 50000
    samples_dir = "../../data/outputs/binning/sshic/0kb/"
    output_dir = "../../data/outputs/statistics/sshic/"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    samples = os.listdir(samples_dir)

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(statistics.run, [(cis_range,
                                        samples_dir+samp,
                                        output_dir) for samp in samples])

    else:
        for samp in samples:
            statistics.run(
                cis_range=cis_range,
                binned_contacts_path=samples_dir+samp,
                output_dir=output_dir
            )


if __name__ == "__main__":
    modes = ['binning']
    if 'filter' in modes:
        do_filter()
    if 'format' in modes:
        do_format()
    if 'binning' in modes:
        do_binning()
    if 'statistics' in modes:
        do_stats()

    print('--- DONE ---')
