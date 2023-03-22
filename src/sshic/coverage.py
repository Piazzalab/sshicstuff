import os
import re
import numpy as np
import pandas as pd
from tools import is_debug
import multiprocessing as mp

from universal.binning import rebin_contacts
from universal.utils import remove_columns

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def main(
        fragments_path: str,
        hic_contacts_path: str,
        output_dir: str
):

    sample_id = re.search(r"AD\d+", hic_contacts_path).group()
    df_fragments = pd.read_csv(fragments_path, sep='\t')
    df_fragments.rename(columns={'chrom': 'chr', 'start_pos': 'start', 'end_pos': 'end'}, inplace=True)
    df_fragments['id'] = df_fragments.index.values
    df_hic_contacts = pd.read_csv(hic_contacts_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])

    df_coverage = df_fragments[['chr', 'start', 'end']]
    df_coverage['contacts'] = np.nan

    df_merged_a = df_hic_contacts.merge(df_fragments[['id', 'chr', 'start', 'end']],
                                        left_on='frag_a',
                                        right_on='id',
                                        suffixes=('', '_a')).drop(columns=['frag_a', 'frag_b'])

    df_merged_b = df_hic_contacts.merge(df_fragments[['id', 'chr', 'start', 'end']],
                                        left_on='frag_b',
                                        right_on='id',
                                        suffixes=('', '_b')).drop(columns=['frag_a', 'frag_b'])

    df_grouped_a = df_merged_a.groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()
    df_grouped_b = df_merged_b.groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()

    df_grouped = pd.concat(
        (df_grouped_a, df_grouped_b)).groupby(by=['id', 'chr', 'start', 'end'], as_index=False).sum()

    df_grouped.index = df_grouped.id
    df_grouped.drop(columns=['id'], inplace=True)

    for bs in bins_sizes_list:
        bin_dir = str(bs // 1000) + 'kb/'
        if not os.path.exists(output_dir+bin_dir):
            os.makedirs(output_dir+bin_dir)
        if not os.path.exists(output_dir+'bedgraph/'+bin_dir):
            os.makedirs(output_dir+'bedgraph/'+bin_dir)
        if bs == 0:
            df_grouped.to_csv(
                output_dir+'bedgraph/'+bin_dir+sample_id+'_coverage_per_fragment.bedgraph',
                sep='\t', index=False, header=False)
            df_grouped.to_csv(output_dir+bin_dir+sample_id+'_coverage_per_fragment.tsv', sep='\t', index=False)
            continue

        df_rebinned = rebin_contacts(
            df_unbinned=df_grouped,
            bin_size=bs,
            chromosomes_coord_path=centromeres_positions,
        )

        df_rebinned_bed = remove_columns(df_rebinned, exclusion=['start', 'genome_bins', 'size'])
        df_rebinned_bed.rename(columns={'chr_bins': 'start'}, inplace=True)
        df_rebinned_bed.insert(2, 'end', df_rebinned_bed['start']+bs)
        df_rebinned_bed.to_csv(
            output_dir+'bedgraph/'+bin_dir+sample_id+'_coverage_per_fragment_rebinned.bedgraph',
            sep='\t', index=False, header=False)

        df_rebinned.to_csv(output_dir+bin_dir+sample_id+'_coverage_per_fragment_rebinned.tsv', sep='\t', index=False)

    print(sample_id)


if __name__ == "__main__":

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    inputs_dir = data_dir + 'inputs/'
    outputs_dir = data_dir + 'outputs/'
    fragments_list = inputs_dir + "fragments_list.txt"
    hicstuff_dir = outputs_dir + "hicstuff/"
    coverage_dir = outputs_dir + "coverage/"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"

    bins_sizes_list = [0, 1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]

    parallel = True
    if is_debug():
        parallel = False

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        print("Making coverage for each digested fragment in the genome")
        samples_dir = hicstuff_dir + sshic_dir
        samples = np.unique(os.listdir(samples_dir))
        if not os.path.exists(coverage_dir+sshic_dir+'bedgraph/'):
            os.makedirs(coverage_dir+sshic_dir+'bedgraph/')

        if parallel:
            with mp.Pool(int(mp.cpu_count()*0.75)) as p:
                p.starmap(main, [(
                    fragments_list,
                    samples_dir+samp,
                    coverage_dir+sshic_dir) for samp in samples]
                )
        else:
            for samp in samples:
                main(
                    fragments_path=fragments_list,
                    hic_contacts_path=samples_dir+samp,
                    output_dir=coverage_dir+sshic_dir
                )

    print('-- DONE --')
