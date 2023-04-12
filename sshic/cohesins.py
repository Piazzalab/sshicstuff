#! /usr/bin/env python3
import numpy as np
import pandas as pd
import os
import re

import universal.utils as tools

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def main(
        df_peaks: pd.DataFrame,
        df_probes: pd.DataFrame,
        df_centro: pd.DataFrame,
        fragments_path: str,
        output_dir: str
):
    samples_id = re.search(r"AD\d+", fragments_path).group()
    fragments = pd.unique(df_probes['frag_id'].astype(str))
    probes_high_quality_sum = ["18535", "18589", "18605", "18611", "18614", "18666", "18694"]
    probes_to_average = {
        'Average_left': ['18535', '18589', '18605', '18611', '18614', '18616'],
        'Average_right': ['18621', '18632', '18634', '18666', '18694'],
        'Average_3Left_(2599-3081-3728)': ['18605', '18611', '18614'],
        'Average_4left_(less_4kb)': ['18605', '18611', '18614', '18616'],
        'Average_left_(2599-3081-3728-6065)': ['18605', '18611', '18614', '18589'],
        'Average_3right_(1439-2715-2954)': ['18621', '18632', '18634'],
        'Average_4right_(1429-2715-2954-8072)': ['18621', '18632', '18634', '18666'],
        'Average_2right_(2954-8072)': ['18634', '18666'],
        'Average_right_(1439-2715)': ['18621', '18632']
    }

    control_probes = ['8579', '32542', '68339', '5315', '65930', '38864', '30750']
    df_contacts: pd.DataFrame = pd.read_csv(fragments_path, sep='\t')
    df_contacts["Sum_7high"] = df_contacts[probes_high_quality_sum].sum(axis=1)

    for colname, colfrag in probes_to_average.items():
        df_contacts[colname] = df_contacts[colfrag].mean(axis=1)

    df_contacts.drop(columns=[
        c for c in fragments if c not in control_probes
    ], inplace=True)

    col_of_interest = [c for c in df_contacts.columns if c not in ['chr', 'positions', 'sizes']]

    df_merged: pd.DataFrame = pd.merge(df_contacts, df_peaks, on='chr')
    df_filtered: pd.DataFrame = df_merged.loc[
        (df_merged['positions'] >= df_merged['interval_start']) &
        (df_merged['positions'] + df_merged['sizes'] < df_merged['interval_end'])
    ]

    del df_merged

    df_grouped1: pd.DataFrame = df_filtered.groupby(
        by=['chr', 'start', 'end', 'interval_start', 'interval_end', 'interval_size', 'interval_shifted'],
        as_index=False).sum(numeric_only=True).drop(columns=['positions', 'sizes'])

    df_grouped2: pd.DataFrame = df_grouped1.groupby(
        by=['chr', 'start', 'end', 'interval_size'], as_index=False).sum(numeric_only=True)
    df_grouped2.drop(columns=['interval_start', 'interval_end', 'interval_size', 'interval_shifted'], inplace=True)

    for c in col_of_interest:
        df_grouped1[c] = df_grouped1[c].div(df_grouped1['interval_size'])
        df_grouped2[c] = df_grouped2[c].div(df_grouped2['end'] - df_grouped2['start'] + 1)

    df_merged2: pd.DataFrame = df_grouped1.merge(df_grouped2, on=['chr', 'start', 'end'])
    df_res = df_merged2.iloc[:, :7]

    for c in col_of_interest:
        df_res[c] = df_merged2[c+'_x'] / df_merged2[c+'_y']

    df_aggregated: pd.DataFrame = df_res.groupby(by=['interval_shifted'], as_index=False).mean(numeric_only=True)
    df_aggregated.drop(
        columns=['start', 'end', 'interval_start', 'interval_end', 'interval_size', 'interval_shifted'], inplace=True)

    df_aggregated.to_csv(output_dir + "enrichment_cohesins_peaks_intervals_"+samples_id+".tsv", sep='\t')
    return df_aggregated


def merge(
        wt_res: dict,
        output_dir: str
):
    df_merged = pd.DataFrame()
    for k, v in wt_res.items():
        df_merged = pd.concat((df_merged, v))
    df_merged = df_merged.groupby(level=0).mean(numeric_only=True)
    df_merged.to_csv(
        output_dir + "average_enrichment_cohesins_peaks_intervals_" + '-'.join(list(wt_res.keys())) + '.tsv', sep='\t')


if __name__ == "__main__":

    inputs_dir = os.path.dirname(os.getcwd()) + '/inputs/'
    outputs_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/outputs/'
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    genes_list_path = inputs_dir + "genes_sacCer3.tsv"
    binning_dir = outputs_dir + "binned/"
    pondered_dir = outputs_dir + "pondered/"
    cohesins_dir = outputs_dir + "cohesins/"
    cohesins_peaks_path = inputs_dir + "HB65_reference_peaks_score50min.bed"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    df_cohesins_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                                    names=['chr', 'start', 'end', 'uid', 'score'])
    df_cohesins_peaks.drop(columns='end', inplace=True)
    df_cohesins_peaks = tools.sort_by_chr(df_cohesins_peaks, col1='chr', col2='start')

    df_cohesins_peaks.drop(['uid', 'score'], axis=1, inplace=True)
    df_inter_peaks_10kb = pd.DataFrame()

    for chrom in pd.unique(df_cohesins_peaks['chr']):
        df = df_cohesins_peaks.loc[df_cohesins_peaks['chr'] == chrom]
        df['end'] = df['start'].shift(-1) - 1
        df = df.iloc[:-1, :]
        df['size'] = df['end'] - df['start'] + 1
        df_filtered = df.loc[df['size'] >= 10000]
        df_inter_peaks_10kb = pd.concat((df_inter_peaks_10kb, df_filtered))

    df_inter_peaks_10kb.index = range(len(df_inter_peaks_10kb))
    df_inter_peaks_10kb[['start', 'end', 'size']] = df_inter_peaks_10kb[['start', 'end', 'size']].astype('int64')

    rebin = df_inter_peaks_10kb['size'].min() // 1000 + 1
    interval_linspace = np.hstack(
        [np.linspace(row['start'], row['end']+1, rebin, endpoint=False).astype('int64')
         for _, row in df_inter_peaks_10kb.iterrows()]
    )

    interval_sizes = np.hstack(
        [np.diff(np.linspace(row['start'], row['end']+1, rebin, endpoint=False).astype('int64')).mean().astype(int)
         for _, row in df_inter_peaks_10kb.iterrows()]
    )
    interval_sizes = interval_sizes.repeat(rebin)

    new_df = df_inter_peaks_10kb.loc[df_inter_peaks_10kb.index.values.repeat(rebin)].reset_index(drop=True)
    new_df['interval_start'] = interval_linspace
    new_df['interval_end'] = interval_linspace + interval_sizes
    new_df['interval_size'] = interval_sizes
    new_df['interval_shifted'] = (new_df['interval_start'] - new_df['start']) // interval_sizes

    new_df.drop(columns=['size'], inplace=True)
    df_probes2frag = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)
    df_centro_pos = pd.read_csv(centromeres_positions, sep='\t', index_col=None)

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        not_binned_dir = binning_dir + sshic_dir + '0kb/'
        samples = [f for f in sorted(os.listdir(not_binned_dir)) if 'frequencies' in f]

        if not os.path.exists(cohesins_dir+sshic_dir):
            os.makedirs(cohesins_dir+sshic_dir)

        wt_df = {}
        for samp in samples:
            samp_id = re.search(r"AD\d+", samp).group()
            if samp_id in ["AD162", "AD242", "AD296", "AD300"]:
                print(samp_id)
                df_aggregated = main(
                    df_peaks=new_df,
                    df_probes=df_probes2frag,
                    df_centro=df_centro_pos,
                    fragments_path=not_binned_dir+samp,
                    output_dir=cohesins_dir+sshic_dir
                )
                wt_df[samp_id] = df_aggregated

        merge(wt_res=wt_df, output_dir=cohesins_dir+sshic_dir)
