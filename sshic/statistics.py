import os
import sys
import getopt
import re
import numpy as np
import pandas as pd
from typing import Optional


def get_stats(
        contacts_unbinned_path: str,
        sparse_contacts_path: str,
        wt_reference: Optional[str] = None,
        cis_range: int = 5000
):

    sample_id = re.search(r"AD\d+", contacts_unbinned_path).group()
    sample_dir = os.path.dirname(contacts_unbinned_path)
    data_dir = os.path.dirname(sample_dir)
    output_path = os.path.join(sample_dir, sample_id)

    probes_to_fragments_path: str = os.path.join(data_dir, "probes_to_fragments.tsv")
    if not os.path.exists(probes_to_fragments_path):
        from probe2fragment import associate_probes_to_fragments
        associate_probes_to_fragments(
            fragments_list_path=os.path.join(data_dir, "fragments_list.txt"),
            oligos_capture_path=os.path.join(data_dir, "capture_oligo_positions.csv")
        )
        
    df_probes: pd.DataFrame = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)

    chr_size_dict: dict = {
        'chr1': 230218, 'chr2': 813184, 'chr3': 316620, 'chr4': 1531933, 'chr5': 576874, 'chr6': 270161,
        'chr7': 1090940, 'chr8': 562643, 'chr9': 439888, 'chr10': 745751, 'chr11': 666816, 'chr12': 1078177,
        'chr13': 924431, 'chr14': 784333, 'chr15': 1091291, 'chr16': 948066, 'mitochondrion': 85779, '2_micron': 6318}

    chr_list = list(chr_size_dict.keys())

    df_unbinned_contacts: pd.DataFrame = pd.read_csv(contacts_unbinned_path, sep='\t')
    df_unbinned_contacts = df_unbinned_contacts.astype(dtype={'chr': str, 'start': int, 'sizes': int})

    df_sparse_contacts: pd.DataFrame = \
        pd.read_csv(sparse_contacts_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])
    #   from sparse_matrix (hicstuff results): get total contacts from which probes enrichment is calculated
    total_sparse_contacts = sum(df_sparse_contacts["contacts"])

    cis_contacts = []
    trans_contacts = []
    intra_chr_contacts = []
    inter_chr_contacts = []
    total_contacts = []
    total_contacts_inter = []

    chr_contacts_nrm = {k: [] for k in chr_size_dict}
    chr_inter_only_contacts_nrm = {k: [] for k in chr_size_dict}

    index = 0
    for probe, row in df_probes.iterrows():
        frag_id = str(row['frag_id'])
        sub_df = df_unbinned_contacts[['chr', 'start', frag_id]]
        cis_limits = [
            int(df_probes.loc[probe, 'probe_start']) - cis_range,
            int(df_probes.loc[probe, 'probe_end']) + cis_range
        ]
        probe_contacts = np.sum(sub_df[frag_id].values)
        total_contacts.append(probe_contacts)
        total_contacts_inter.append(np.sum(sub_df.query("chr != @row['chr']")[frag_id].values))

        if total_contacts[index] > 0:
            cis_contacts.append(
                sub_df.loc[(sub_df["chr"] == row['chr']) &
                           (sub_df['start'] > cis_limits[0]) &
                           (sub_df['start'] < cis_limits[1]), frag_id].sum() / total_contacts[index])

            trans_contacts.append(1 - cis_contacts[index])
            intra_chr_contacts.append(sub_df.loc[sub_df['chr'] == row['chr'], frag_id].sum() / total_contacts[index])
            inter_chr_contacts.append(sub_df.loc[sub_df['chr'] != row['chr'], frag_id].sum() / total_contacts[index])
        else:
            cis_contacts.append(0.)
            trans_contacts.append(0.)
            intra_chr_contacts.append(0.)
            inter_chr_contacts.append(0.)

        for chrom in chr_list:
            #   n1: sum contacts chr_i
            #   d1: sum contacts all chr
            #   chrom_size: chr_i's size
            #   genome_size: sum of sizes for all chr except frag_chr
            #   c1 : normalized contacts on chr_i for frag_j
            chrom_size = chr_size_dict[chrom]
            genome_size = sum([s for c, s in chr_size_dict.items() if c != row['chr']])
            n1 = sub_df.loc[sub_df['chr'] == chrom, frag_id].sum()
            if n1 == 0:
                chr_contacts_nrm[chrom].append(0)
            else:
                d1 = total_contacts[index]
                c1 = (n1/d1) / (chrom_size/genome_size)
                chr_contacts_nrm[chrom].append(c1)

            #   n2: sum contacts chr_i if chr_i != probe_chr
            #   d2: sum contacts all inter chr (exclude the probe_chr)
            #   c2 : normalized inter chr contacts on chr_i for frag_j
            n2 = sub_df.loc[
                (sub_df['chr'] == chrom) &
                (sub_df['chr'] != row['chr']), frag_id].sum()

            if n2 == 0:
                chr_inter_only_contacts_nrm[chrom].append(0)
            else:
                d2 = total_contacts_inter[index]
                c2 = (n2 / d2) / (chrom_size / genome_size)
                chr_inter_only_contacts_nrm[chrom].append(c2)
        index += 1

    df_stats: pd.DataFrame = pd.DataFrame({
        'probes': df_probes.index.values, 'fragments': df_probes["frag_id"].values,
        'types': df_probes["type"].values, 'contacts': total_contacts, 'cis': cis_contacts,
        'trans': trans_contacts, 'intra_chr': intra_chr_contacts, 'inter_chr': inter_chr_contacts
    })

    df_stats['contacts_over_hic_contacts'] = df_stats['contacts'] / total_sparse_contacts

    #  capture_efficiency_vs_dsdna : number of contact for one oligo divided
    #  by the mean of all other 'ds' oligos in the genome
    n3 = df_stats.loc[:, 'contacts']
    d3 = np.mean(df_stats.loc[df_stats['types'] == 'ds', 'contacts'])
    df_stats['capture_efficiency_vs_dsdna'] = n3 / d3

    if wt_reference is not None:
        df_wt: pd.DataFrame = pd.read_csv(wt_reference, sep='\t')
        df_stats[f"capture_efficiency_vs_wt"] = np.nan
        for index, row in df_stats.iterrows():
            probe = row['probe']
            wt_capture_eff = df_wt.loc[df_wt['probes'] == probe, "Capture_efficiency_WT"]
            if wt_capture_eff > 0:
                df_stats.loc[index, f"capture_efficiency_vs_wt"] = wt_capture_eff

    df_chr_nrm = pd.DataFrame({
        'probes': df_probes.index.values, 'fragments': df_probes["frag_id"].values, 'types': df_probes["type"].values
    })

    df_chr_inter_only_nrm = df_chr_nrm.copy(deep=True)

    for chr_id in chr_list:
        df_chr_nrm[chr_id] = chr_contacts_nrm[chr_id]
        df_chr_inter_only_nrm[chr_id] = chr_inter_only_contacts_nrm[chr_id]

    df_stats.to_csv(output_path + '_global_statistics.tsv', sep='\t')
    df_chr_nrm.to_csv(output_path + '_normalized_chr_freq.tsv', sep='\t')
    df_chr_inter_only_nrm.to_csv(output_path + '_normalized_inter_chr_freq.tsv', sep='\t')


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    try:
        opts, args = getopt.getopt(
            argv,
            "hu:s:w:", [
                "--help",
                "--contacts",
                "--sparse",
                "--wildtype"]
        )

    except getopt.GetoptError:
        print('Making some statistics and normalization around the contacts make for each probe:\n'
              '-u <unbinned_contacts.tsv> (generated by fragments) \n'
              '-s <sparse_contacts_input.txt> (generated by hicstuff) \n'
              '-w <wt_capture_efficiency> (Optional, if you want to pondered sample) \n'
              )
        sys.exit(2)

    contacts_unbinned_input, sparse_contacts_input = ['' for _ in range(2)]
    wt_reference_input = None
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('Making some statistics and normalization around the contacts make for each probe:\n'
                  '-u <unbinned_contacts.tsv> (generated by fragments) \n'
                  '-s <sparse_contacts_input.txt> (generated by hicstuff) \n'
                  '-w <wt_capture_efficiency> (Optional, if you want to pondered sample) \n'
                  )
            sys.exit()
        elif opt in ("-u", "--contacts"):
            contacts_unbinned_input = arg
        elif opt in ("-s", "--sparse"):
            sparse_contacts_input = arg
        elif opt in ("-w", "--wildtype"):
            wt_reference_input = arg

    get_stats(
        contacts_unbinned_path=contacts_unbinned_input,
        sparse_contacts_path=sparse_contacts_input,
        wt_reference=wt_reference_input,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
