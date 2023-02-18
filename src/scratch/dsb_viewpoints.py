import os
import re

import numpy as np
import pandas as pd

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

if __name__ == "__main__":
    excluded_chr = ['chr3', 'chr5']
    data_dir = '../../data/'
    samples_dir = data_dir + 'inputs/HiC_WT_2h_4h/samples/'
    fragments_dir = data_dir + 'inputs/HiC_WT_2h_4h/'
    samples = sorted(os.listdir(samples_dir))

    #   position on both end of the breaks, 6 bins of 1000 bp
    chr5_dsb_left_bins = [105000, 106000, 107000, 108000, 109000, 110000]
    chr5_dsb_right_bins = [125000, 126000, 127000, 128000, 129000, 130000]
    #   same but for the control
    chr8_ctrl_bins = [140000, 141000, 142000, 143000, 144000, 145000]

    #   read the centromeres coordinates table
    #   from the centromere df create a telomeres one.
    #       telo_left -> 0
    #       telo_right -> chr length
    #   read the cohesins peaks table and keeps only peaks with score >= 100
    #   read the fragment library file
    df_centros = pd.read_csv(data_dir+'inputs/S288c_chr_centro_coordinates.tsv', sep='\t', index_col=None)
    df_telos = pd.DataFrame({'chr': df_centros['chr'], 'telo_l': 0, 'telo_r': df_centros['length']})
    df_cohesins = pd.read_csv(data_dir+'inputs/HB65_reference_peaks_score50min.bed',
                              sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'uid', 'score'])
    df_cohesins = df_cohesins.loc[df_cohesins['score'] >= 100, :]
    df_fragments = pd.read_csv(fragments_dir+'AD154to160_S288c_DSB_cutsite_q20_chrs_1kb.frag.tsv',
                               sep='\t', index_col=None)

    chr_to_fragid = {}
    for c in np.unique(df_fragments.chr):
        chr_to_fragid[c] = df_fragments[df_fragments['chr'] == c].index.tolist()

    #   Only keep fragments that belong to the bins of interest
    df_fragments_filtered = df_fragments[
        ((df_fragments.chr == 'chr5') &
         ((df_fragments.start_pos.isin(chr5_dsb_left_bins)) |
         (df_fragments.start_pos.isin(chr5_dsb_right_bins))) |
         (df_fragments.chr == 'chr8') &
         (df_fragments.start_pos.isin(chr8_ctrl_bins)))
    ]

    for samp in samples:
        #   just the name of the sample (AD157, AD254, ...)
        samp_id = re.search(r"AD\d+", samp).group()
        print(samp_id)
        #   df0: raw dataframe, dense matrix
        #   may require a lot of time to import the table as it is heavy table
        #   NB : df1 is a square matrix
        df0 = pd.read_csv(samples_dir+samp, sep=' ', header=None)
        mat = df0.to_numpy(dtype=float)
        for chrom, frags in chr_to_fragid.items():
            start = frags[0]
            end = frags[-1] + 1
            mat[start:end, start:end] = np.nan

        df1 = pd.DataFrame(mat)
        #   df2: only keep columns corresponding to fragment index present in df_fragment_filtered,
        #       i.e.,  bins of interest
        #   Inter normalization
        df1 = df1.div(df1.sum(axis=0))
        df2 = df1.filter(items=df_fragments_filtered.index.tolist(), axis=1)
        #   group the 6 bins at the left of the breaks site by mean and add columns for result
        df2['chr5_dsb_left'] = df2.iloc[:, 0:6].mean(axis=1)
        #   group the 6 bins on the rigth by mean and add columns for result
        df2['chr5_dsb_right'] = df2.iloc[:, 6:12].mean(axis=1)
        #   same for the 6 bins of the control on chr8
        df2['chr8_ctrl'] = df2.iloc[:, 12:18].mean(axis=1)
        #   df3: df2 but with only the 3 averaged columns
        df3 = df2[['chr5_dsb_left', 'chr5_dsb_right', 'chr8_ctrl']]
        #   add chromosome number column at hte beginning of th df3
        df3.insert(0, 'chr', df_fragments.chr)
        #   df4: create a deep copy of df3.
        #       df3 is binned at 1kb
        #       df4 will be binned at 10 kb
        df4 = df3.copy(deep=True)
        #   add chr_bins columns for both df3 and df4
        df3.insert(1, 'chr_bins', df_fragments.start_pos)
        #   modify the binning of df4, from 1kb to 10kb
        df4.insert(1, 'chr_bins', (df_fragments.start_pos // 10000)*10000)
        #   df5 : df4 but we grouped the row sharing the same bin
        df5 = df4.groupby(['chr', 'chr_bins'], as_index=False).sum()

        """
        ##########################
        #   CENTROMERES AGGREGATED
        ##########################
        """
        #   output dir for centromeres aggregated
        cen_dir = '../../data/outputs/hic/' + 'centromeres/' + samp_id + '/'
        if not os.path.exists(cen_dir):
            os.makedirs(cen_dir)

        #   merge centromere df and df5
        #   only keep contacts with the 150000 around each centromere
        df_merged1 = pd.merge(df5, df_centros, on='chr')
        df_merged_cen_areas = df_merged1[
            (df_merged1.chr_bins > (df_merged1.left_arm_length - 150000 - 10000)) &
            (df_merged1.chr_bins < (df_merged1.left_arm_length + 150000))
        ]

        #   indices shifting i.e., bin of the centromere becomes bin 0 etc ...
        #   NB : we use absolute value to make easy the groupby method in a further step
        df_merged_cen_areas['chr_bins'] = \
            abs(df_merged_cen_areas['chr_bins'] - (df_merged_cen_areas['left_arm_length'] // 10000) * 10000)

        #   df6: groupby chromosomes AND chr_bins, with an average, the rows of df_merged_cen_areas
        df6 = df_merged_cen_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)
        #   remove useless columns
        df6.drop(columns=['length', 'left_arm_length', 'right_arm_length'], axis=1, inplace=True)
        #   df7: df6 but without chr3 and chr5
        df7 = df6[~df6.chr.isin(excluded_chr)]

        #   Initialize df for mean, std and median
        df_cen_mean = pd.DataFrame()
        df_cen_std = pd.DataFrame()
        df_cen_median = pd.DataFrame()
        for col in ['chr5_dsb_left', 'chr5_dsb_right', 'chr8_ctrl']:
            #   make a pivot of the table :
            #       chromosomes becoming the columns
            #       chr_bins the row index
            #       applied for each column ('chr5_dsb_left', 'chr5_dsb_right', 'chr8_ctrl') one by one
            #       we go from one dataframe to 3 new ones
            df_cen_chr = df7.pivot_table(index='chr_bins', columns='chr', values=col, fill_value=0.)
            df_cen_chr.to_csv(cen_dir + col + '_chr1-16_freq_cen.tsv', sep='\t')

            #   compute the mean, std and median for each column
            #   add the result the corresponding tables
            df_cen_mean[col] = df_cen_chr.mean(axis=1)
            df_cen_median[col] = df_cen_chr.median(axis=1)
            df_cen_std[col] = df_cen_chr.std(axis=1)

        #   save the result in tsv tables
        df_cen_mean.to_csv(cen_dir + 'mean_on_cen.tsv', sep='\t')
        df_cen_std.to_csv(cen_dir + 'std_on_cen.tsv', sep='\t')
        df_cen_median.to_csv(cen_dir + 'median_on_cen.tsv', sep='\t')

        print("\t DONE : centromeres")

        """
        ##########################
        #   TELOMERES AGGREGATED
        ##########################
        """

        #   No comment here, same process as centromeres.
        telo_dir = '../../data/outputs/hic/' + 'telomeres/' + samp_id + '/'
        if not os.path.exists(telo_dir):
            os.makedirs(telo_dir)

        df_merged2 = pd.merge(df5, df_telos, on='chr')
        df_merged_telo_areas = df_merged2[
            (df_merged2.chr_bins < (df_merged2.telo_l + 150000 + 10000)) |
            (df_merged2.chr_bins > (df_merged2.telo_r - 150000))
        ]

        df_merged_telo_areas.loc[df_merged2.chr_bins >= (df_merged2.telo_r - 150000), 'chr_bins'] = \
            abs(df_merged_telo_areas['chr_bins'] - (df_merged_telo_areas['telo_r'] // 10000) * 10000)

        df8 = df_merged_telo_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)
        df8.drop(columns=['telo_l', 'telo_r'], axis=1, inplace=True)
        df9 = df8[~df8.chr.isin(excluded_chr)]

        for c in df9.columns[2:]:
            df9[c] /= df9[c].sum()

        df_telo_mean = pd.DataFrame()
        df_telo_std = pd.DataFrame()
        df_telo_median = pd.DataFrame()
        for col in ['chr5_dsb_left', 'chr5_dsb_right', 'chr8_ctrl']:
            df_telo_chr = df9.pivot_table(index='chr_bins', columns='chr', values=col, fill_value=0.)
            df_telo_chr.to_csv(telo_dir + col + '_chr1-16_freq_telo.tsv', sep='\t')

            df_telo_mean[col] = df_telo_chr.mean(axis=1)
            df_telo_median[col] = df_telo_chr.median(axis=1)
            df_telo_std[col] = df_telo_chr.std(axis=1)

        df_telo_mean.to_csv(telo_dir + 'mean_on_telo.tsv', sep='\t')
        df_telo_std.to_csv(telo_dir + 'std_on_telo.tsv', sep='\t')
        df_telo_median.to_csv(telo_dir + 'median_on_telo.tsv', sep='\t')

        print("\t DONE : telomeres")

        """
        ##########################
        #   COHESINS AGGREGATED
        ##########################
        """

        cohesins_dir = '../../data/outputs/hic/' + 'cohesins/' + samp_id + '/'
        if not os.path.exists(cohesins_dir):
            os.makedirs(cohesins_dir)

        df_merged3 = pd.merge(df3, df_cohesins, on='chr')
        df_merged_cohesins_areas = df_merged3[
            (df_merged3.chr_bins > (df_merged3.start - 15000 - 1000)) &
            (df_merged3.chr_bins < (df_merged3.end + 15000))
        ]

        df_merged_cohesins_areas.sort_values(by=['start', 'chr_bins'], inplace=True)
        df_merged_cohesins_areas['chr_bins'] = \
            abs(df_merged_cohesins_areas['chr_bins'] - (df_merged_cohesins_areas['start'] // 1000) * 1000)

        df10 = df_merged_cohesins_areas.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)
        df10.drop(columns=['start', 'end', 'score'], axis=1, inplace=True)
        df11 = df10[~df10.chr.isin(excluded_chr)]

        for c in df11.columns[2:]:
            df11[c] /= df11[c].sum()

        df_cohesins_mean = pd.DataFrame()
        df_cohesins_std = pd.DataFrame()
        df_cohesins_median = pd.DataFrame()
        for col in ['chr5_dsb_left', 'chr5_dsb_right', 'chr8_ctrl']:
            df_cohesins_chr = df11.pivot_table(index='chr_bins', columns='chr', values=col, fill_value=0.)
            df_cohesins_chr.to_csv(cohesins_dir + col + '_chr1-16_freq_telo.tsv', sep='\t')

            df_cohesins_mean[col] = df_cohesins_chr.mean(axis=1)
            df_cohesins_median[col] = df_cohesins_chr.median(axis=1)
            df_cohesins_std[col] = df_cohesins_chr.std(axis=1)

        df_cohesins_mean.to_csv(cohesins_dir + 'mean_on_cohesins.tsv', sep='\t')
        df_cohesins_std.to_csv(cohesins_dir + 'std_on_cohesins.tsv', sep='\t')
        df_cohesins_median.to_csv(cohesins_dir + 'median_on_cohesins.tsv', sep='\t')

        print("\t DONE : cohesins")
        print("\n")

