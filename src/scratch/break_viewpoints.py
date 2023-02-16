import os
import re
import numpy as np
import pandas as pd

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

if __name__ == "__main__":
    excluded_chr = ['chr3', 'chr5']
    data_dir = '../../data/tmp/'
    samples_dir = data_dir + 'HiC_WT_2h_4h/'
    samples = sorted(os.listdir(samples_dir))

    chr5_dsb_left_bins = [105000, 106000, 107000, 108000, 109000, 110000]
    chr5_dsb_right_bins = [125000, 126000, 127000, 128000, 129000, 130000]
    chr8_ctrl_bins = [140000, 141000, 142000, 143000, 144000, 145000]

    df_centros = pd.read_csv(data_dir+'S288c_chr_centro_coordinates.tsv', sep='\t', index_col=None)
    df_fragments = pd.read_csv(data_dir+'AD154to160_S288c_DSB_cutsite_q20_chrs_1kb.frag.tsv', sep='\t', index_col=None)

    df_fragments_filtered = df_fragments[
        ((df_fragments.chr == 'chr5') &
         ((df_fragments.start_pos.isin(chr5_dsb_left_bins)) |
         (df_fragments.start_pos.isin(chr5_dsb_right_bins))) |
         (df_fragments.chr == 'chr8') &
         (df_fragments.start_pos.isin(chr8_ctrl_bins)))
    ]

    for samp in samples:
        samp_id = re.search(r"AD\d+", samp).group()
        df1 = pd.read_csv(samples_dir+samp, sep=' ', header=None)
        df2 = df1.filter(items=df_fragments_filtered.index.tolist(), axis=1)
        df2['chr5_dsb_left'] = df2.iloc[:, 0:6].mean(axis=1)
        df2['chr5_dsb_right'] = df2.iloc[:, 6:12].mean(axis=1)
        df2['chr8_ctrl'] = df2.iloc[:, 12:18].mean(axis=1)
        df3 = df2[['chr5_dsb_left', 'chr5_dsb_right', 'chr8_ctrl']]

        df3.insert(0, 'chr', df_fragments.chr)
        df4 = df3.copy(deep=True)
        df3.insert(1, 'chr_bins', df_fragments.start_pos)
        df4.insert(1, 'chr_bins', (df_fragments.start_pos // 10000)*10000)
        df5 = df4.groupby(['chr', 'chr_bins'], as_index=False).sum()

        df_merged = pd.merge(df5, df_centros, on='chr')

        """
        ##########################
        #   CENTROMERES AGGREGATED
        ##########################
        """
        df_merged_cen_areas = df_merged[
            (df_merged.chr_bins > (df_merged.left_arm_length - 150000 - 10000)) &
            (df_merged.chr_bins < (df_merged.left_arm_length + 150000))
        ]

        df_merged_cen_areas['chr_bins'] = \
            abs(df_merged_cen_areas['chr_bins'] - (df_merged_cen_areas['left_arm_length'] // 10000) * 10000)

        df_cen = df_merged_cen_areas.groupby(['chr', 'chr_bins'], as_index=False).mean()
        df_cen.drop(columns=['length', 'left_arm_length', 'right_arm_length'], axis=1, inplace=True)

        """
        ##########################
        #   TELOMERES AGGREGATED
        ##########################
        """

        """
        ##########################
        #   COHESINS AGGREGATED
        ##########################
        """

        pass

