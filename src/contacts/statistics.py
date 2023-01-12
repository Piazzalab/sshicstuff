#! /usr/bin/env python3

import re
import numpy as np
import pandas as pd


def compute_stats(formatted_contacts_path: str,
                  fragments_to_oligos_path: str,
                  cis_range: int,
                  output_path: str):

    """
    After having formatted the contacts of each oligos in the genome (file with bin_size=0)
    We know cant to use this new formatted file to compute basics statistics like total number of contacts,
    frequencies of contacts intra .vs. inter chromosomes, cis .vs. trans oligos
     (with cis range as an input value given by the user, in bp)
    """

    chr_size = {
        'chr1': 230218, 'chr2': 813184, 'chr3': 316620, 'chr4': 1531933, 'chr5': 576874, 'chr6': 270161,
        'chr7': 1090940, 'chr8': 562643, 'chr9': 439888, 'chr10': 745751, 'chr11': 666816, 'chr12': 1078177,
        'chr13': 924431, 'chr14': 784333, 'chr15': 1091291, 'chr16': 948066, 'mitochondrion': 85779, '2_micron': 6318}

    genome_size = sum(chr_size.values())
    chr_size_normalized = {k: v/genome_size for k, v in chr_size.items()}
    df_formatted_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_info = pd.read_csv(fragments_to_oligos_path, sep='\t', index_col=0)

    fragments = [f for f in df_formatted_contacts.columns.values if re.search(r"\d+", f) is not None]
    probes = np.asarray(df_info.loc['oligo', :].values, dtype='<U64')
    types = np.asarray(df_info.loc['type', :].values, dtype='<U8')

    cis_contacts = []
    trans_contacts = []
    intra_chr_contacts = []
    inter_chr_contacts = []
    total_contacts = []
    chr_contacts_nrm = {k: [] for k in chr_size}

    for ii_f, frag in enumerate(fragments):
        sub_df = df_formatted_contacts[['chr', 'positions', frag]]
        cis_limits = [int(df_info.loc['frag_start', frag]) - cis_range, int(df_info.loc['frag_end', frag]) + cis_range]
        frag_chr = df_info.loc['frag_chr', frag]
        total_contacts.append(np.sum(sub_df[frag].values))
        cis_contacts.append(
            np.sum(
                sub_df.query(
                    "chr == @frag_chr and positions > @cis_limits[0] and positions <@cis_limits[1]")[frag].values) /
            total_contacts[ii_f]
        )
        trans_contacts.append(1 - cis_contacts[ii_f])
        intra_chr_contacts.append(
            np.sum(sub_df.query("chr == @frag_chr")[frag].values) / total_contacts[ii_f]
        )
        inter_chr_contacts.append(
            np.sum(sub_df.query("chr != @frag_chr")[frag].values) / total_contacts[ii_f]
        )

        for chrom, size in chr_size_normalized.items():
            chr_contacts_nrm[chrom].append(
                (np.sum(sub_df.query("chr == @chrom")[frag].values) / total_contacts[ii_f]) / size
            )

    df_global = pd.DataFrame({'fragments': fragments, 'probes': probes, 'types': types, 'total': total_contacts,
                              'cis': cis_contacts, 'trans': trans_contacts, 'intra_chr': intra_chr_contacts,
                              'inter_chr': inter_chr_contacts})

    #  fold over : number of contact for one oligo divided by the mean (or median)
    #  of all other 'ds' oligos in the genome
    df_global['fold_over'] = \
        df_global.loc[:, 'total'] / np.mean(df_global.loc[df_global['types'] == 'ds', 'total'].values)

    df_chr_nrm = pd.DataFrame({'fragments': fragments, 'probes': probes, 'types': types})
    for chr_id, freqs in chr_contacts_nrm.items():
        df_chr_nrm[chr_id] = freqs

    df_global.to_csv(output_path + '_global_statistics.tsv', sep='\t')
    df_chr_nrm.to_csv(output_path + '_normalized_chr_freq.tsv', sep='\t')


def run(
        cis_range: int,
        formatted_contacts_path: str,
        fragments_to_oligos_path: str,
        output_dir: str):

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_path = output_dir + sample_id
    compute_stats(
        cis_range=cis_range,
        formatted_contacts_path=formatted_contacts_path,
        fragments_to_oligos_path=fragments_to_oligos_path,
        output_path=output_path)

    print('DONE: ', sample_id)

