#! /usr/bin/env python3
import os.path
import re
import numpy as np
import pandas as pd
import sshic.tools as tl


def get_fragments_contacts(
        filtered_contacts_path: str,
        output_dir: str):
    """
    This function will count the number of contacts for each read that comes from column 'frag_x' in each bin.
    The results are stored in three dictionaries, given as arguments.
        contacts_res of the form :  {oligoX : {chrX_binA : n ... } ...}
        all_contacted_pos of the form : {oligoX : {chrX_1456 : n, chrY_89445: m ... } ...}
    """

    samp_id = re.search(r"AD\d+", filtered_contacts_path).group()
    output_path = output_dir + samp_id

    df = pd.read_csv(filtered_contacts_path, sep='\t')
    contacts = pd.DataFrame(columns=['chr', 'positions', 'sizes'])
    contacts = contacts.astype(dtype={'chr': str, 'positions': int, 'sizes': int})
    frequencies = contacts.copy(deep=True)

    for x in ['a', 'b']:
        #   if x = a get b, if x = b get a
        y = tl.frag2(x)
        df2 = df[~pd.isna(df['name_' + x])]
        unique_frag = pd.unique(df2['frag_'+x])
        for frag in unique_frag:
            df3 = df2[df2['frag_'+x] == frag]

            tmp_c = pd.DataFrame({'chr': df3['chr_'+y], 'positions': df3['start_'+y],
                                  'sizes': df3['size_'+y], frag: df3['contacts']})

            tmp_f = tmp_c.copy(deep=True)
            tmp_f[frag] /= np.sum(tmp_f[frag])

            contacts = pd.concat([contacts, tmp_c])
            frequencies = pd.concat([frequencies, tmp_f])

    group_c = contacts.groupby(by=['chr', 'positions', 'sizes'], as_index=False)
    group_f = frequencies.groupby(by=['chr', 'positions', 'sizes'], as_index=False)

    res_c = group_c.sum()
    res_f = group_f.sum()

    res_c = tl.sort_by_chr(res_c, 'chr', 'positions')
    res_f = tl.sort_by_chr(res_f, 'chr', 'positions')

    res_c.index = range(len(res_c))
    res_c.index = range(len(res_f))

    res_c.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    res_f.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)

    print('DONE : ', samp_id)


def rebin_contacts(
        not_binned_samp_path: str,
        bin_size: int,
        output_dir: str
):

    samp_id = re.search(r"AD\d+", not_binned_samp_path).group()
    df = pd.read_csv(not_binned_samp_path, sep='\t')
    bin_dir = output_dir + str(bin_size // 1000) + 'kb/'
    if not os.path.exists(bin_dir):
        os.makedirs(bin_dir)
    output_path = bin_dir+samp_id
    fragments = [f for f in df.columns if re.match(r'\d+', str(f))]
    df.insert(2, 'chr_bins', (df["positions"] // bin_size) * bin_size)
    df_binned_contacts = df.groupby(["chr", "chr_bins"], as_index=False).sum()
    df_binned_contacts.drop(['positions', 'sizes'], axis=1, inplace=True)

    df_binned_contacts = tl.sort_by_chr(df_binned_contacts, 'chr', 'chr_bins')
    df_binned_frequencies = df_binned_contacts.copy(deep=True)
    for frag in fragments:
        df_binned_frequencies[frag] /= sum(df_binned_contacts[frag])

    df_binned_contacts.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    df_binned_frequencies.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)
