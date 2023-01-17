import os
import re
import multiprocessing as mp
import numpy as np
from contacts import filter, format, binning, statistics, nucleosomes
import utils.tools as tools


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
        samples_dir: str,
        output_dir: str,
        parallel: bool = True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique(os.listdir(samples_dir))

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


def do_binning(
        artificial_genome: str,
        samples_dir: str,
        output_dir: str,
        parallel: bool = True):
    bin_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000, 10**9]

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
                                         samples_dir+samp,
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
        samples_dir: str,
        output_dir: str,
        cis_span: int,
        parallel: bool = True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique(
        [re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir)])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(statistics.run, [(cis_span,
                                        samples_dir+samp+'_contacts.tsv',
                                        samples_dir + samp + '_frag_to_prob.tsv',
                                        output_dir) for samp in samples])

    else:
        for samp in samples:
            statistics.run(
                cis_range=cis_span,
                formatted_contacts_path=samples_dir+samp+'_contacts.tsv',
                fragments_to_oligos_path=samples_dir+samp+'_frag_to_prob.tsv',
                output_dir=output_dir
            )


def do_nucleo(
        samples_dir: str,
        fragments: str,
        statistics_dir: str,
        nucleosomes_path: str,
        output_dir: str,
        parallel: bool = True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = np.unique(
        [re.search(r"AD\d+", f).group() for f in os.listdir(samples_dir)])

    if parallel:
        with mp.Pool(mp.cpu_count()) as p:
            p.starmap(nucleosomes.run, [(samples_dir+samp+'_contacts.tsv',
                                         statistics_dir+samp+'_global_statistics.tsv',
                                         fragments,
                                         nucleosomes_path,
                                         output_dir) for samp in samples])

    else:
        for samp in samples:
            nucleosomes.run(
                formatted_contacts_path=samples_dir+samp+'_contacts.tsv',
                statistics_path=statistics_dir+samp+'_global_statistics.tsv',
                fragments_list_path=fragments,
                nucleosomes_path=nucleosomes_path,
                output_dir=output_dir
            )


if __name__ == "__main__":
    sshic_dir = ['sshic/', 'sshic_pcrdupkept/']
    modes = ['nucleosomes']

    for hicd in sshic_dir:
        print(hicd)

        #  ARGUMENTS
        #######################################################
        fragments_list = "../../data/inputs/fragments_list.txt"
        artificial_genome_fa = "../../data/inputs/S288c_DSB_LY_capture_artificial.fa"
        oligos_positions = "../../data/inputs/capture_oligo_positions.csv"
        nucleosomes_free_regions = "../../data/inputs/Chereji_Henikoff_genome_research_NFR.bed"
        hicstuff_dir = "../../data/outputs/hicstuff/" + hicd
        filter_output_dir = "../../data/outputs/filtered/" + hicd
        format_output_dir = "../../data/outputs/formatted/" + hicd
        binning_output_dir = "../../data/outputs/binning/" + hicd
        statistics_output_dir = "../../data/outputs/statistics/" + hicd
        nfr_output_dir = "../../data/outputs/nucleosomes/" + hicd

        parallel_state: bool = True
        if tools.is_debug():
            parallel_state = False
        #######################################################

        if 'filter' in modes:
            print('Filtering')
            do_filter(
                fragments=fragments_list,
                oligos=oligos_positions,
                samples_dir=hicstuff_dir,
                output_dir=filter_output_dir
            )

        if 'format' in modes:
            print('Formatting')
            do_format(
                samples_dir=filter_output_dir,
                output_dir=format_output_dir,
                parallel=parallel_state
            )

        if 'binning' in modes:
            print('Binning')
            do_binning(
                artificial_genome=artificial_genome_fa,
                samples_dir=format_output_dir,
                output_dir=binning_output_dir,
                parallel=parallel_state
            )

        if 'statistics' in modes:
            print('Statistics')
            do_stats(
                samples_dir=format_output_dir,
                output_dir=statistics_output_dir,
                cis_span=50000,
                parallel=parallel_state
            )

        if 'nucleosomes' in modes:
            print('nucleosomes')
            do_nucleo(
                samples_dir=format_output_dir,
                fragments=fragments_list,
                statistics_dir=statistics_output_dir,
                nucleosomes_path=nucleosomes_free_regions,
                output_dir=nfr_output_dir,
                parallel=parallel_state
            )

    print('--- DONE ---')
