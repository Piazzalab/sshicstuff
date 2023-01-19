import re
import numpy as np
import pandas as pd
from utils import tools


def fragments_to_oligos(
        fragments_list_path: str,
        oligos_capture_path: str,
        output_path: str):

    df_fragments = pd.read_csv(fragments_list_path, sep='\t')
    df_oligos = pd.read_csv(oligos_capture_path, sep=",")
    df_probes_in_frag = pd.DataFrame()
    df_probes_in_frag.index = ['type', 'probe_start', 'probe_end', 'chr', 'frag_id', 'frag_start', 'frag_end']

    for index, row in df_oligos.iterrows():
        chrom, probe_start, probe_end, probe_type, probe, probe_seq = row
        sub_df_fragments = df_fragments[df_fragments['chrom'] == chrom]
        oligo_middle = int(probe_start + (probe_end-probe_start)/2)
        nearest_frag_start = tools.find_nearest(
            array=sub_df_fragments['start_pos'], key=oligo_middle, mode='lower'
        )
        frag_id = sub_df_fragments.index[sub_df_fragments['start_pos'] == nearest_frag_start].tolist()[0]
        frag_start = sub_df_fragments.loc[frag_id, 'start_pos']
        frag_end = sub_df_fragments.loc[frag_id, 'end_pos']
        df_probes_in_frag[probe] = [probe_type, probe_start, probe_end, chrom, frag_id, frag_start, frag_end]

    df_probes_in_frag.to_csv(output_path, sep='\t', index_label='probe')


def format_fragments_contacts(
        df: pd.DataFrame,
        output_path: str):
    """
    This function will count the number of contacts for each read that comes from column 'frag_x' in each bin.
    The results are stored in three dictionaries, given as arguments.
        contacts_res of the form :  {oligoX : {chrX_binA : n ... } ...}
        all_contacted_pos of the form : {oligoX : {chrX_1456 : n, chrY_89445: m ... } ...}
    """
    contacts = pd.DataFrame(columns=['chr', 'positions', 'sizes'])
    contacts = contacts.astype(dtype={'chr': str, 'positions': int, 'sizes': int})
    frequencies = contacts.copy(deep=True)

    for x in ['a', 'b']:
        #   if x = a get b, if x = b get a
        y = tools.frag2(x)
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

    res_c.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    res_f.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)


def run(
        filtered_contacts_path: str,
        output_dir: str):

    sample_id = re.search(r"AD\d+", filtered_contacts_path).group()
    df_contacts_filtered = pd.read_csv(filtered_contacts_path, sep=',')

    format_fragments_contacts(
        df=df_contacts_filtered,
        output_path=output_dir+sample_id)

    print('DONE: ', sample_id)
