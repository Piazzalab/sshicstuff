import os
import re
import numpy as np
import pandas as pd
from tools import is_debug
import multiprocessing as mp

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

    # Write .bedgraph files :
    if not os.path.exists(output_dir+'bedgraph/'):
        os.makedirs(output_dir+'bedgraph/')
    df_grouped.to_csv(output_dir+'bedgraph/'+sample_id+'_coverage_per_fragment.bedgraph',
                      sep='\t',
                      index=False,
                      header=False)

    print(sample_id)


if __name__ == "__main__":

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    inputs_dir = data_dir + 'inputs/'
    outputs_dir = data_dir + 'outputs/'
    fragments_list = inputs_dir + "fragments_list.txt"
    hicstuff_dir = outputs_dir + "hicstuff/"
    coverage_dir = outputs_dir + "coverage/"

    parallel = True
    if is_debug():
        parallel = False

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        print("Making coverage for each digested fragment in the genome")
        samples_dir = hicstuff_dir + sshic_dir
        samples = np.unique(os.listdir(samples_dir))
        if not os.path.exists(coverage_dir + sshic_dir):
            os.makedirs(coverage_dir + sshic_dir)

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
