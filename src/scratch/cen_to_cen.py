import os
import re
import numpy as np
import pandas as pd


if __name__ == "__main__":
    bin_size = 1000
    window_size = 60000
    excluded_chr = ['chr3', 'chr5']
    data_dir = '../../data/tmp/'
    samples_dir = data_dir + 'HiC_WT_2h_4h/'
    samples = sorted(os.listdir(samples_dir))

    output_dir = data_dir + 'outputs/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_centros = pd.read_csv(data_dir+'S288c_chr_centro_coordinates.tsv', sep='\t', index_col=None)
    df_fragments = pd.read_csv(data_dir+'AD154to160_S288c_DSB_cutsite_q20_chrs_1kb.frag.tsv', sep='\t', index_col=None)

    df_merged_frag_centros = pd.merge(df_fragments, df_centros, on='chr')
    df_fragments_filtered = df_merged_frag_centros[
        (df_merged_frag_centros.start_pos > (df_merged_frag_centros.left_arm_length-window_size-bin_size)) &
        (df_merged_frag_centros.start_pos < (df_merged_frag_centros.left_arm_length+window_size)) &
        (~df_merged_frag_centros.chr.isin(excluded_chr))
    ]
    df_fragments_filtered_bis = df_fragments_filtered[
        (df_fragments_filtered.start_pos == (df_fragments_filtered.left_arm_length // bin_size)*bin_size)
    ]

    fragid_to_chr = dict(pd.Series(df_fragments_filtered_bis.iloc[:, 1]))

    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        df1 = pd.read_csv(samples_dir+samp, sep=' ', header=None)
        #   remove fragments on row that are not in the windows [-nkb -- centromere -- + nkb]
        df2 = df1.filter(items=df_fragments_filtered.index.tolist(), axis=0)
        #   only keep on columns fragments that are on the centromere's bin
        df3 = df2.filter(items=df_fragments_filtered_bis.index.tolist(), axis=1)
        #   add columns with chr ID for each fragment on row
        df3.insert(0, 'chr', df_fragments_filtered.chr)
        #   add columns with bin for each fragment on row
        df3.insert(1, 'chr_bins', df_fragments_filtered.start_pos)
        #   indices shifting for bins in 'chr_bins' column
        #   use absolute value to allow a groupby method in a further step
        df3['chr_bins'] = abs(df3['chr_bins']-(df_fragments_filtered['left_arm_length'] // bin_size)*bin_size)
        #   replace columns name by chr names
        df3.rename(columns=fragid_to_chr, inplace=True)

        #   always favour deepcopy for dataframe
        df3_inter = df3.copy(deep=True)
        #   set n.a.n for all inter contacts i.e., contacts made by chr X with chr X
        for c in fragid_to_chr.values():
            df3_inter.loc[df3_inter['chr'] == c, c] = np.nan

        df3['mean'] = df3.loc[:, fragid_to_chr.values()].mean(axis=1)
        df3_inter['mean'] = df3_inter.loc[:, fragid_to_chr.values()].mean(axis=1)

        df4 = df3.loc[:, ['chr', 'chr_bins', 'mean']]
        df4_inter = df3_inter.loc[:, ['chr', 'chr_bins', 'mean']]

        #   pivot the table
        #   left and right bins are merged such as a groupby method, followed by a mean method
        df5 = df4.pivot_table(index='chr_bins', columns=['chr'],
                              values='mean', aggfunc=np.mean, fill_value=0)
        df5_inter = df4_inter.pivot_table(index='chr_bins', columns=['chr'],
                                          values='mean', aggfunc=np.mean, fill_value=0)

        res = df5.copy(deep=True)
        res_inter = df5_inter.copy(deep=True)

        for b in res.index:
            #   1 : freq divided by the row's sum
            #   2 : freq of inter contacts divided by row's sum of inter contacts
            res.loc[b, :] /= res.loc[b, :].sum()
            res_inter.loc[b, :] /= res_inter.loc[b, :].sum()

        output_path = output_dir + samp_id
        res.to_csv(output_path + '_freq.tsv', sep='\t')
        res_inter.to_csv(output_path + '_freq_inter.tsv', sep='\t')

