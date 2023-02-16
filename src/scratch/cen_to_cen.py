import os
import re
import numpy as np
import pandas as pd


def compute_sums(df: pd.DataFrame,
                 chr2frag: dict):
    sum1 = df.sum(axis=0)
    matrix = df.to_numpy(dtype=float)

    for chrom, frags in chr2frag.items():
        start = frags[0]
        end = frags[-1]+1
        matrix[start:end, start:end] = np.nan

    sum2 = pd.DataFrame(matrix).sum(axis=0)
    return sum1, sum2


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

    chr_to_fragid = {}
    for c in np.unique(df_fragments.chr):
        chr_to_fragid[c] = df_fragments[df_fragments['chr'] == c].index.tolist()


    df_merged_frag_centros = pd.merge(df_fragments, df_centros, on='chr')
    df_fragments_filtered = df_merged_frag_centros[
        (df_merged_frag_centros.start_pos > (df_merged_frag_centros.left_arm_length-window_size-bin_size)) &
        (df_merged_frag_centros.start_pos < (df_merged_frag_centros.left_arm_length+window_size)) &
        (~df_merged_frag_centros.chr.isin(excluded_chr))
    ]
    df_fragments_filtered_bis = df_fragments_filtered[
        (df_fragments_filtered.start_pos == (df_fragments_filtered.left_arm_length // bin_size)*bin_size)
    ]

    cen_bin_to_chr = dict(pd.Series(df_fragments_filtered_bis.iloc[:, 1]))

    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        df1 = pd.read_csv(samples_dir+samp, sep=' ', header=None)
        sum_, sum_inter = compute_sums(df1, chr_to_fragid)

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
        df3.rename(columns=cen_bin_to_chr, inplace=True)

        df3['mean'] = df3.loc[:, cen_bin_to_chr.values()].mean(axis=1) / sum_
        df3['mean_inter'] = df3.loc[:, cen_bin_to_chr.values()].mean(axis=1) / sum_inter

        df4 = df3.loc[:, ['chr', 'chr_bins', 'mean', 'mean_inter']].fillna(0.)

        #   pivot the table
        #   left and right bins are merged such as a groupby method, followed by a mean method
        #   res1 : normalized over all contacts in df1
        #   res2 : normalized over all inter contacts in df1
        res1 = df4.pivot_table(index='chr_bins', columns=['chr'], values='mean', aggfunc=np.mean, fill_value=0)
        res2 = df4.pivot_table(index='chr_bins', columns=['chr'], values='mean_inter', aggfunc=np.mean, fill_value=0)

        output_path = output_dir + samp_id
        res1.to_csv(output_path + '_freq.tsv', sep='\t')
        res2.to_csv(output_path + '_freq_inter.tsv', sep='\t')

        print(samp_id)

