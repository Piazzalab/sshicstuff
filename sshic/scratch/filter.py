import pandas as pd
import re
import os
import numpy as np
import multiprocessing as mp
from utils import is_debug


def oligos_correction(oligos_path):
    oligos = pd.read_csv(oligos_path, sep=",")
    oligos.columns = [oligos.columns[i].lower() for i in range(len(oligos.columns))]
    oligos.sort_values(by=['chr', 'start'], inplace=True)
    oligos.reset_index(drop=True, inplace=True)

    return oligos


def fragments_correction(fragments_path):
    fragments = pd.read_csv(fragments_path, sep='\t')
    fragments = pd.DataFrame({'frag': [k for k in range(len(fragments))],
                              'chr': fragments['chrom'],
                              'start': fragments['start_pos'],
                              'end': fragments['end_pos'],
                              'size': fragments['size'],
                              'gc_content': fragments['gc_content']
                              })
    return fragments


def starts_match(fragments, oligos):
    """
    If the capture oligo is inside a fragment, it changes the start of
    the oligos dataframe with the fragments starts.
    """
    l_starts = []
    for i in range(len(oligos)):
        oligos_chr = oligos['chr'][i]
        middle = int((oligos['end'][i] - oligos['start'][i] - 1) / 2 + oligos['start'][i] - 1)
        if oligos_chr == 'chr_artificial':
            for k in reversed(range(len(fragments))):
                interval = range(fragments['start'][k], fragments['end'][k])
                fragments_chr = fragments['chr'][k]
                if middle in interval and fragments_chr == oligos_chr:
                    l_starts.append(fragments['start'][k])
                    break
        else:
            for k in range(len(fragments)):
                interval = range(fragments['start'][k], fragments['end'][k] + 1)
                fragments_chr = fragments['chr'][k]

                if middle in interval and fragments_chr == oligos_chr:
                    l_starts.append(fragments['start'][k])
                    break
    oligos['start'] = list(l_starts)
    return oligos


def oligos_fragments_joining(fragments, oligos):
    """
    Removes the fragments that does not contain an oligo region, puts all the columns of fragments_list
    that corresponds. It also changes the starts and ends columns by the fragments ones.
    """
    oligos = starts_match(fragments, oligos)
    oligos.set_index(['chr', 'start'])
    oligos.pop("end")
    fragments.set_index(['chr', 'start'])
    oligos_fragments = fragments.merge(oligos, on=['chr', 'start'])
    oligos_fragments.sort_values(by=['chr', 'start'])
    return oligos_fragments


def contacts_correction(contacts_path):
    """
    Re-organizes the contacts file
    """
    contacts = pd.read_csv(contacts_path, sep='\t', header=None)
    contacts.drop([0], inplace=True)
    contacts.reset_index(drop=True, inplace=True)
    contacts.columns = ['frag_a', 'frag_b', 'contacts']

    return contacts


def first_join(x, oligos_fragments, contacts):
    """
    Join the contacts and the oligos_fragments dataframes keeping only
    the rows that have their 'x' frag (frag_a or frag_b, see contacts_correction function)
    """
    joined = contacts.merge(oligos_fragments, left_on=x, right_on='frag', how='inner')
    return joined


def frag2(x):
    if x == 'frag_a':
        y = 'frag_b'
    else:
        y = 'frag_a'
    return y


def second_join(x, fragments, oligos_fragments, contacts):
    """
    Adds the fragments file information (=columns) for the y fragment after the first join (see first_join function)
    Only for y because the x fragments have already their information because it is the oligos_fragments.
    """
    new_contacts = first_join(x, oligos_fragments, contacts)
    y = frag2(x)
    joined = new_contacts.join(fragments.drop("frag", axis=1),
                               on=y,
                               lsuffix='_' + x[-1],
                               rsuffix='_' + y[-1], how='left')

    # puts a suffix to know what fragment corresponds to an oligo
    joined.rename(columns={"type": "type_" + x[-1],
                           "name": "name_" + x[-1],
                           "sequence": "sequence_" + x[-1]
                           },
                  inplace=True)
    return joined


def main(
        oligos_path: str,
        fragments_path: str,
        contacts_path: str,
        output_path: str
):
    """
    Does the two joining (for 'frag_a' and 'frag_b') and then concatenate the two results
    """

    sample_id = re.search(r"AD\d+", contacts_path).group()
    fragments = fragments_correction(fragments_path)
    oligos = oligos_correction(oligos_path)
    contacts = contacts_correction(contacts_path)
    oligos_fragments = oligos_fragments_joining(fragments, oligos)
    df1 = second_join('frag_a', fragments, oligos_fragments, contacts)
    df2 = second_join('frag_b', fragments, oligos_fragments, contacts)

    contacts_joined = pd.concat([df1, df2])
    contacts_joined.drop("frag", axis=1, inplace=True)
    contacts_joined.sort_values(by=['frag_a', 'frag_b', 'start_a', 'start_b'], inplace=True)
    contacts_filtered = contacts_joined.convert_dtypes().reset_index(drop=True)
    contacts_filtered.to_csv(output_path+sample_id+'_filtered.tsv', sep='\t', index=False)


if __name__ == "__main__":

    data_dir = os.path.dirname(os.getcwd()) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    hicstuff_dir = outputs_dir + "sparse/"
    filter_dir = outputs_dir + "filtered/"
    fragments_list = inputs_dir + "fragments_list.txt"
    oligos_positions = inputs_dir + "capture_oligo_positions.csv"
    parallel = True
    if is_debug():
        parallel = False

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        print("Filtered hicstuff sparse matrix by keeping only contacts made by probe")
        samples_dir = hicstuff_dir + sshic_dir
        samples = np.unique(os.listdir(samples_dir))
        if not os.path.exists(filter_dir + sshic_dir):
            os.makedirs(filter_dir + sshic_dir)
        if parallel:
            with mp.Pool(mp.cpu_count()) as p:
                p.starmap(main, [(
                    oligos_positions,
                    fragments_list,
                    samples_dir+samp,
                    filter_dir+sshic_dir) for samp in samples]
                )
        else:
            for samp in samples:
                main(
                    oligos_path=oligos_positions,
                    fragments_path=fragments_list,
                    contacts_path=samples_dir+samp,
                    output_path=filter_dir+sshic_dir
                )

        print('-- DONE --')
