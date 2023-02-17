#! /usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
from typing import Optional


def run(
        cis_range: int,
        sparse_mat_path: Optional[str],
        wt_references_dir: Optional[str],
        samples_vs_wt: Optional[dict],
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        output_dir: str
):

    """
    This function's purpose is to compute some basic statistics about the contacts made by the oligo probes
    to get a better insight of the single stranded filament behaviour.
    Among these statistics we mainly do :
        . the total number of contacts made by the probe, or at least the read that contains it
        . that amount of contacts normalized by the total number of contacts measured
        in the sparse matrix (whole hic experiment)
        . the number of 'cis' contacts i.e., contacts made by probe with very close (defined areas as an argument)
        space
        . the number of 'trans' contacts, basically what is not cis contacts, so 1 - cis
        . the intra chromosomal contacts i.e., contacts made by the probe on the same chromosome as its
        . the inter chromosomal contacts i.e., contacts by the probe made on other chromosome, so 1 - intra
        . the dsdna_norm_capture_efficiency, number of contacts for one oligo divided by the mean (or median)
        of all other 'ds' oligos in the genome
        . (for mutants only) the capture efficiency by probe normalized over that of the wt

    These statistics are stored as 'global' statistics.
    Moreover, we also generate tables for :
        . the contacts made by a probe with each chromosome, normalized by the all the contacts made by
        the probe, itself normalized by the length of the chromosome, itself normalized by the size of the genome
        . the contacts made by a probe with only inter chromosome (chromosome of the probe is excluded)
        normalized by the sum of contacts made on all inter chromosomes, itself normalized by chromosomes length,
        itself normalized by genome size.

    ARGUMENTS
    ______________

    cis_range : int
        number of bases (bp) around each probe where we consider that if the probe did a contact in there,
        it is a 'cis' contact. Outside this window is will be a 'trans' contacts
    sparse_mat_path : str
        raw matrix (sparse matrix) given by hicstuff with all contacts made in hic. The file has
        three columns : frag_a | frag_b | contacts
    wt_references_dir : str
        the path to the directory containing the capture efficiency averaged on the wild types
    samples_vs_wt : dict
        dictionary of samples that need to be weighted over the WT references.
        keys are the wt time point like 2h, 4h, 6h etc ...
        values are lists of samples names to be pondered using the key reference wt
    formatted_contacts_path :  str
        path to the formatted contacts files of the current sample, basically the not_binned_contacts
        previously made with the function get_fragments_contacts in the binning script
    probes_to_fragments_path  :  str
        path to the table (.tsv file) that makes correspond the name of the probe with that of the fragment
        that contains it and complementary information, generated in the format script.
    output_dir  :  str
                the absolute path toward the output directory to save the results
    """
    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_path = output_dir + sample_id

    """
    After having formatted the contacts of each oligos in the genome (file with bin_size=0)
    We know cant to use this new formatted file to compute basics statistics like total number of contacts,
    frequencies of contacts intra .vs. inter chromosomes, cis .vs. trans oligos
     (with cis range as an input value given by the user, in bp)
    """

    chr_size_dict = {
        'chr1': 230218, 'chr2': 813184, 'chr3': 316620, 'chr4': 1531933, 'chr5': 576874, 'chr6': 270161,
        'chr7': 1090940, 'chr8': 562643, 'chr9': 439888, 'chr10': 745751, 'chr11': 666816, 'chr12': 1078177,
        'chr13': 924431, 'chr14': 784333, 'chr15': 1091291, 'chr16': 948066, 'mitochondrion': 85779, '2_micron': 6318}

    unique_chr = list(chr_size_dict.keys())
    df_formatted_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    all_probes = np.asarray(df_probes.index, dtype='<U64')

    df_sparse_mat = pd.read_csv(sparse_mat_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])
    #   from sparse_matrix (hicstuff results): get total contacts from which probes enrichment is calculated
    total_sparse_contacts = sum(df_sparse_mat["contacts"])

    wt_capture_ref = False
    if os.path.exists(wt_references_dir) and os.listdir(wt_references_dir) is not None:
        wt_capture_ref = True
        ref_wt: dict = \
            {k.lower().split('.')[0]: pd.read_csv(wt_references_dir + k, sep='\t')
             for k in os.listdir(wt_references_dir)}
    else:
        ref_wt: dict = {}

    probes = []
    fragments = []
    types = []
    cis_contacts = []
    trans_contacts = []
    intra_chr_contacts = []
    inter_chr_contacts = []
    total_contacts = []
    total_contacts_inter = []

    chr_contacts_nrm = {k: [] for k in chr_size_dict}
    chr_inter_only_contacts_nrm = {k: [] for k in chr_size_dict}

    ii_probe = 0
    dsdna_counter = 0
    for probe in all_probes:
        probe_type, probe_start, probe_end, probe_chr, frag_id, frag_start, frag_end = df_probes.loc[probe].tolist()
        frag_id = str(frag_id)

        if frag_id not in df_formatted_contacts.columns.tolist():
            continue

        probes.append(probe)
        fragments.append(frag_id)
        types.append(probe_type)
        sub_df = df_formatted_contacts[['chr', 'positions', frag_id]]
        cis_limits = \
            [int(df_probes.loc[probe, 'probe_start']) - cis_range, int(df_probes.loc[probe, 'probe_end']) + cis_range]
        probe_contacts = np.sum(sub_df[frag_id].values)
        total_contacts.append(probe_contacts)
        total_contacts_inter.append(np.sum(sub_df.query("chr != @probe_chr")[frag_id].values))

        cis_contacts.append(
            np.sum(
                sub_df.query(
                    "chr == @probe_chr and positions > @cis_limits[0] and positions <@cis_limits[1]")[frag_id].values) /
            total_contacts[ii_probe]
        )
        trans_contacts.append(1 - cis_contacts[ii_probe])
        intra_chr_contacts.append(
            np.sum(sub_df.query("chr == @probe_chr")[frag_id].values) / total_contacts[ii_probe]
        )
        inter_chr_contacts.append(
            np.sum(sub_df.query("chr != @probe_chr")[frag_id].values) / total_contacts[ii_probe]
        )

        for chrom in chr_size_dict:

            #   n1: sum contacts chr_i
            #   d1: sum contacts all chr
            #   chrom_size: chr_i's size
            #   genome_size: sum of sizes for all chr except frag_chr
            #   c1 : normalized contacts on chr_i for frag_j
            chrom_size = chr_size_dict[chrom]
            genome_size = sum([s for c, s in chr_size_dict.items() if c != probe_chr])
            n1 = np.sum(sub_df.query("chr == @chrom")[frag_id].values)
            if n1 == 0:
                chr_contacts_nrm[chrom].append(0)
            else:
                d1 = total_contacts[ii_probe]
                c1 = (n1/d1) / (chrom_size/genome_size)
                chr_contacts_nrm[chrom].append(c1)

            #   n2: sum contacts chr_i if chr_i != frag_chr
            #   d2: sum contacts all inter chr (exclude the frag_chr)
            #   c2 : normalized inter chr contacts on chr_i for frag_j
            n2 = np.sum(sub_df.query("chr == @chrom and chr != @probe_chr")[frag_id].values)
            if n2 == 0:
                chr_inter_only_contacts_nrm[chrom].append(0)
            else:
                d2 = total_contacts_inter[ii_probe]
                c2 = (n2 / d2) / (chrom_size / genome_size)
                chr_inter_only_contacts_nrm[chrom].append(c2)
        ii_probe += 1

    df_global = pd.DataFrame({'probes': probes, 'fragments': fragments, 'types': types,
                              'total_contacts': total_contacts, 'cis': cis_contacts, 'trans': trans_contacts,
                              'intra_chr': intra_chr_contacts, 'inter_chr': inter_chr_contacts})

    df_global['total_contacts_pondered'] = df_global['total_contacts'] / total_sparse_contacts

    if wt_capture_ref:
        #  dsdna_norm_capture_efficiency : number of contact for one oligo divided by the mean (or median)
        #  of all other 'ds' oligos in the genome
        df_global['dsdna_norm_capture_efficiency'] = \
            df_global.loc[:, 'total_contacts'] / np.mean(df_global.loc[df_global['types'] == 'ds',
                                                                       'total_contacts'].values)

        #   capture_efficiency_norm_'+wt : divide the dsDNA-normalized contacts to
        #   the WT_capture efficiency for each probe to get the correction factor for each probe
        for wt in ref_wt:
            if sample_id in samples_vs_wt[wt]:
                wt_capture_eff_values = \
                    df_global.merge(ref_wt[wt], on='probes')['Capture_efficiency_WT'].values
                df_global['capture_efficiency_norm_'+wt] = \
                    df_global['dsdna_norm_capture_efficiency'] / wt_capture_eff_values
            else:
                continue

    df_chr_nrm = pd.DataFrame({'probes': probes, 'fragments': fragments, 'types': types})
    df_chr_inter_only_nrm = df_chr_nrm.copy(deep=True)

    for chr_id in unique_chr:
        df_chr_nrm[chr_id] = chr_contacts_nrm[chr_id]
        df_chr_inter_only_nrm[chr_id] = chr_inter_only_contacts_nrm[chr_id]

    df_global.to_csv(output_path + '_global_statistics.tsv', sep='\t')
    df_chr_nrm.to_csv(output_path + '_normalized_chr_freq.tsv', sep='\t')
    df_chr_inter_only_nrm.to_csv(output_path + '_normalized_inter_chr_only_freq.tsv', sep='\t')
