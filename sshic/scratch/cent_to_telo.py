import os
import re
import numpy as np
import pandas as pd

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

if __name__ == "__main__":
    bin_size = 1000
    centro_region_size = 150000
    telo_region_size = 30000
    excluded_chr = ['chr3', 'chr2', 'chr5']
    data_dir = '../../data/'
    samples_dir = data_dir + 'inputs/HiC_WT_2h_4h/samples/'
    fragments_dir = data_dir + 'inputs/HiC_WT_2h_4h/'
    chr_arm = data_dir + "inputs/S288c_chr_arm_sorted_by_length.csv"
    samples = sorted(os.listdir(samples_dir))

    output_dir = '../../data/outputs/hic/' + 'cen2telo/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_centros = pd.read_csv(data_dir+'inputs/S288c_chr_centro_coordinates.tsv', sep='\t', index_col=None)
    df_telos = pd.DataFrame({'chr': df_centros['chr'], 'telo_l': 0, 'telo_r': df_centros['length']})
    df_chr_arm = pd.read_csv(chr_arm, sep='\t', header=None)
    df_chr_arm.columns = ['chr', 'arm', 'size', 'category']

    df_fragments = pd.read_csv(fragments_dir+'AD154to160_S288c_DSB_cutsite_q20_chrs_1kb.frag.tsv',
                               sep='\t', index_col=None)

    chr_to_fragid = {}
    for c in np.unique(df_fragments.chr):
        chr_to_fragid[c] = df_fragments[df_fragments['chr'] == c].index.tolist()

    df_merged1 = pd.merge(df_fragments, df_centros, on='chr')

    #   filter fragments that belongs to excluded chromosomes lits
    df_fragments_filtered1 = df_merged1[
        (~df_merged1.chr.isin(excluded_chr))
    ]

    #   filter fragments that are outside the windows of n kb left-right the centromeres
    df_fragments_filtered2 = df_fragments_filtered1[
        (df_fragments_filtered1.start_pos > (df_fragments_filtered1.left_arm_length - centro_region_size - bin_size)) &
        (df_fragments_filtered1.start_pos < (df_fragments_filtered1.left_arm_length + centro_region_size))
    ]

    df_merged2 = pd.merge(df_fragments_filtered1, df_telos, on='chr')

    #   filter fragments that are not in the telomeres regions
    df_fragments_filtered3l = df_merged2[
        (df_merged2.start_pos < df_merged2.telo_l + telo_region_size + bin_size)
    ]
    df_fragments_filtered3l['arm'] = 'left'

    df_fragments_filtered3r = df_merged2[
        (df_merged2.start_pos > df_merged2.telo_r - telo_region_size - bin_size)
    ]
    df_fragments_filtered3r['arm'] = 'right'
    df_fragments_filtered4 = pd.concat((df_fragments_filtered3l, df_fragments_filtered3r))

    df_merged3 = pd.merge(df_fragments_filtered4, df_chr_arm, on=['chr', 'arm'])
    df_fragments_filtered5 = df_merged3.drop(columns=[
        'id', 'telo_l', 'telo_r', 'left_arm_length', 'right_arm_length', 'length', 'size_x', 'size_y', 'gc_content'])
    df_fragments_filtered5.index = df_fragments_filtered4.index

    small_mask = set(df_fragments_filtered1.index)
    small_all = set(df_fragments_filtered5.loc[df_fragments_filtered5['category'] == 'small'].index)
    small_intersect = small_all.intersection(small_mask)

    middle_mask = set(df_fragments_filtered1.index)
    middle_all = set(df_fragments_filtered5.loc[df_fragments_filtered5['category'] == 'middle'].index)
    middle_intersect = middle_all.intersection(middle_mask)

    long_mask = set(df_fragments_filtered1.index)
    long_all = set(df_fragments_filtered5.loc[df_fragments_filtered5['category'] == 'long'].index)
    long_intersect = long_all.intersection(long_mask)

    categories_of_arm_idx = {
        'small': list(small_intersect),
        'middle': list(middle_intersect),
        'long': list(long_intersect)
    }

    del df_merged1, df_merged2, df_merged3, df_fragments_filtered3r, df_fragments_filtered3l
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
        #   remove fragments on row that are not in the windows [-nkb -- telomeres -- + nkb]
        df3 = df2.filter(items=df_fragments_filtered5.index.tolist(), axis=1)

        df3['small'] = df3.loc[:, categories_of_arm_idx['small']].mean(axis=1)
        df3['middle'] = df3.loc[:, categories_of_arm_idx['middle']].mean(axis=1)
        df3['long'] = df3.loc[:, categories_of_arm_idx['long']].mean(axis=1)

        df4 = df3[list(categories_of_arm_idx.keys())]
        #   add columns with chr ID for each fragment on row
        df4.insert(0, 'chr', df_fragments_filtered2.chr)
        #   add columns with bin for each fragment on row
        df4.insert(1, 'chr_bins', df_fragments_filtered2.start_pos)
        #   indices shifting for bins in 'chr_bins' column
        #   use absolute value to allow a groupby method in a further step
        df4['chr_bins'] = abs(df4['chr_bins'] - (df_fragments_filtered2['left_arm_length'] // bin_size) * bin_size)

        df5 = df4.groupby(by='chr_bins', as_index=False).mean(numeric_only=True)
        output_path = output_dir + samp_id
        df5.to_csv(output_path + 'cen_to_telo_over_arm_length_average.tsv', sep='\t', index=False)

        print(samp_id)

