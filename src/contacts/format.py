import re
import numpy as np
import pandas as pd
from utils import tools


def fragments_to_oligos(
        df: pd.DataFrame,
        output_path: str):

    df_frag2oligo = pd.DataFrame()
    df_frag2oligo.index = ['oligo', 'type', 'frag_chr', 'frag_start', 'frag_end']
    #   res of the form : {oligoX : {name: "probeX", type: "ds" , chr: 'chr1, ...} ...}
    res: dict = {}
    for x in ['a', 'b']:
        sub_df = df[~pd.isna(df['name_' + x])]
        unique_frag = pd.unique(sub_df['frag_'+x])
        for frag in unique_frag:
            #   Nota Bene : A same read or fragment may contain two oligos because they are very close
            #       Thus for a same read we acn see two probe's names that are different
            #       For the moment the infos_res[f][names] is an array that may contain one (in most cases)
            #       or two probes (for two oligos in one read).
            if frag not in res:
                res[frag] = {
                    'probes': pd.unique(sub_df.loc[sub_df['frag_'+x] == frag, 'name_'+x]),
                    'type': sub_df.loc[sub_df['frag_'+x] == frag, 'type_'+x].values[0],
                    'frag_chr': sub_df.loc[sub_df['frag_'+x] == frag, 'chr_'+x].values[0],
                    'frag_start': sub_df.loc[sub_df['frag_' + x] == frag, 'start_' + x].values[0],
                    'frag_end': sub_df.loc[sub_df['frag_' + x] == frag, 'end_' + x].values[0]
                }

        for frag, val in res.items():
            if len(val['probes']) > 1:
                #   If the current fragment or read in the loop has multiple names i.e., contains two oligos
                #   we decided to merge the probe's names in a single one : 'probe1_&_probe2'
                uid = '_&_'.join(res[frag]['probes'])
            else:
                uid = res[frag]['probes'][0]

            df_frag2oligo[frag] = np.array([uid,
                                           val['type'],
                                           val['frag_chr'],
                                           val['frag_start'],
                                           val['frag_end']])

        df_frag2oligo.to_csv(output_path + '_frag_to_prob.tsv', sep='\t')
        return res


def format_fragments_contacts(
        df: pd.DataFrame,
        output_path: str):

    """
    This function will count the number of contacts for each read that comes from column 'frag_x' in each bin.
    The results are stored in three dictionaries, given as arguments.
        contacts_res of the form :  {oligoX : {chrX_binA : n ... } ...}
        all_contacted_pos of the form : {oligoX : {chrX_1456 : n, chrY_89445: m ... } ...}
    """
    contacted_pos_per_chr: dict = {}
    contacted_bins_count_per_frag: dict = {}
    
    for x in ['a', 'b']:
        #   if x = a get b, if x = b get a
        y = tools.frag2(x)
        for ii_f, f in enumerate(df['frag_' + x].values):
            if not pd.isna(df['name_' + x][ii_f]):
                chr_id = df['chr_' + y][ii_f]
                start = df['start_' + y][ii_f]
                chr_and_pos = chr_id + '_' + str(start)

                if f not in contacted_bins_count_per_frag:
                    contacted_bins_count_per_frag[f] = {}

                if chr_id not in contacted_pos_per_chr:
                    contacted_pos_per_chr[chr_id] = set()

                if chr_and_pos not in contacted_pos_per_chr[chr_id]:
                    contacted_pos_per_chr[chr_id].add(start)

                bin_id = chr_and_pos

                if bin_id not in contacted_bins_count_per_frag[f]:
                    contacted_bins_count_per_frag[f][bin_id] = df['contacts'][ii_f]
                else:
                    contacted_bins_count_per_frag[f][bin_id] += df['contacts'][ii_f]

    for chr_id, pos_list in contacted_pos_per_chr.items():
        contacted_pos_per_chr[chr_id] = sorted(pos_list)

    chr_unique_list = np.concatenate([['chr' + str(i) for i in range(1, 17)],
                                      ['2_micron', 'mitochondrion', 'chr_artificial']])

    chr_and_pos = []
    chromosomes = []
    positions = []
    for chr_id in chr_unique_list:
        if chr_id in contacted_pos_per_chr:
            new_pos = contacted_pos_per_chr[chr_id]
            positions.extend(new_pos)
            chromosomes.extend(np.repeat(chr_id, len(new_pos)))
            for npos in new_pos:
                chr_and_pos.append(chr_id + '_' + str(npos))

    chr_and_pos = np.asarray(chr_and_pos)
    df_formatted_contacts = pd.DataFrame({'chr': chromosomes, 'positions': positions})
    df_formatted_frequencies = pd.DataFrame({'chr': chromosomes, 'positions': positions})

    for f in contacted_bins_count_per_frag:
        contacts = np.zeros(len(chr_and_pos), dtype=int)
        for pos in contacted_bins_count_per_frag[f]:
            idx = np.argwhere(chr_and_pos == pos)[0]
            contacts[idx] = contacted_bins_count_per_frag[f][pos]

        df_formatted_contacts[f] = contacts
        df_formatted_frequencies[f] = contacts / np.sum(contacts)

    df_formatted_contacts.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    df_formatted_frequencies.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)


def run(
        filtered_contacts_path: str,
        output_dir: str):

    sample_id = re.search(r"AD\d+", filtered_contacts_path).group()
    df_contacts_filtered = pd.read_csv(filtered_contacts_path, sep=',')

    fragments_to_oligos(
        df=df_contacts_filtered,
        output_path=output_dir+sample_id)

    format_fragments_contacts(
        df=df_contacts_filtered,
        output_path=output_dir+sample_id)

    print('DONE: ', sample_id)
