import os
import re
import numpy as np
import pandas as pd


if __name__ == "__main__":
    bin_size = 1000
    window_size = 100000
    excluded_chr = ['chr3', 'chr2', 'chr5']
    data_dir = '../../data/'
    samples_dir = data_dir + 'inputs/HiC_WT_2h_4h/samples/'
    fragments_dir = data_dir + 'inputs/HiC_WT_2h_4h/'
    samples = sorted(os.listdir(samples_dir))

    output_dir = '../../data/outputs/hic/' + 'cen2rdna/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_centros = pd.read_csv(data_dir+'inputs/S288c_chr_centro_coordinates.tsv', sep='\t', index_col=None)
    df_fragments = pd.read_csv(fragments_dir+'AD154to160_S288c_DSB_cutsite_q20_chrs_1kb.frag.tsv',
                               sep='\t', index_col=None)

    rdna_regions = ['chr12:'+str(b) for b in range(451000, 468000, 1000)]
    chr_to_fragid = {}
    for c in np.unique(df_fragments.chr):
        chr_to_fragid[c] = df_fragments[df_fragments['chr'] == c].index.tolist()

    df_merged_frag_centros = pd.merge(df_fragments, df_centros, on='chr')

    #   filter fragments that belongs to excluded chromosomes lits
    df_fragments_filtered1 = df_merged_frag_centros[
        (~df_merged_frag_centros.chr.isin(excluded_chr))
    ]

    #   filter fragments that are outside the windows of 60kb left-right the centromeres
    df_fragments_filtered2 = df_fragments_filtered1[
        (df_fragments_filtered1.start_pos > (df_fragments_filtered1.left_arm_length-window_size-bin_size)) &
        (df_fragments_filtered1.start_pos < (df_fragments_filtered1.left_arm_length+window_size))
    ]

    #   filter fragments that are not in the chr12 rDNA region 451000:468000
    df_fragments_filtered3 = df_fragments_filtered1[
        (df_fragments_filtered1.start_pos >= 451000) &
        (df_fragments_filtered1.end_pos <= 468000) &
        (df_fragments_filtered1.chr == 'chr12')
    ]

    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        df0 = pd.read_csv(samples_dir+samp, sep=' ', header=None)
        mat = df0.to_numpy(dtype=float)

        for chrom, frags in chr_to_fragid.items():
            start = frags[0]
            end = frags[-1] + 1
            mat[start:end, start:end] = np.nan

        df1 = pd.DataFrame(mat)
        #   remove excluded chromosomes
        df1b = df1.filter(items=df_fragments_filtered1.index.tolist(), axis=0)
        df1c = df1b.filter(items=df_fragments_filtered1.index.tolist(), axis=1)
        #   inter normalization
        #   add 1e-9 to prevent from dividing by zero
        df1d = df1c.div(df1c.sum(axis=0)+1e-9)

        #   remove fragments on row that are not in the windows [-nkb -- centromere -- + nkb]
        df2 = df1d.filter(items=df_fragments_filtered2.index.tolist(), axis=0)


        #   only keep on columns fragments that are on rDNA region
        df3 = df2.filter(items=df_fragments_filtered3.index.tolist(), axis=1)
        #   replace columns name by chr names
        df3.columns = rdna_regions
        #   add columns with chr ID for each fragment on row
        df3.insert(0, 'chr', df_fragments_filtered2.chr)
        #   add columns with bin for each fragment on row
        df3.insert(1, 'chr_bins', df_fragments_filtered2.start_pos)
        #   indices shifting for bins in 'chr_bins' column
        #   use absolute value to allow a groupby method in a further step
        df3['chr_bins'] = abs(df3['chr_bins'] - (df_fragments_filtered2['left_arm_length'] // bin_size) * bin_size)
        df3.insert(2, 'mean', df3[rdna_regions].mean(axis=1))

        df4 = df3.loc[:, ['chr', 'chr_bins', 'mean']]

        #   pivot the table
        #   left and right bins are merged such as a groupby method, followed by a mean method
        #   res : normalized over all inter contacts in df1
        res = df4.pivot_table(index='chr_bins', columns=['chr'],
                              values='mean', aggfunc=np.mean, fill_value=np.nan)

        res['mean'] = res.mean(axis=1)

        output_path = output_dir + samp_id
        res.to_csv(output_path + '_freq_pericentro_to_rdna.tsv', sep='\t')

        print(samp_id)

