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
    fragments_of_interest = ["18535", "18589", "18605", "18611", "18614", "18666", "18694"]

    df_contacts = pd.read_csv(fragments_path, sep='\t')
    bin_size = int(df_contacts.chr_bins[1] - df_contacts.chr_bins[0])
    df_contacts = df_contacts.loc[:, (df_contacts.columns.isin(fragments_of_interest)) |
                                     (df_contacts.columns.isin(['chr', 'chr_bins', 'sizes']))]
    df_contacts.insert(2, 'chr_bins_end', df_contacts.chr_bins + bin_size)
    df_contacts["sum"] = df_contacts[fragments_of_interest].sum(axis=1)
    df_contacts.drop(columns=fragments_of_interest, inplace=True)

    df_merged = pd.merge(df_peaks, df_contacts, on='chr')
    df_filtered = df_merged.loc[
        (df_merged['inter'] >= df_merged['chr_bins']) &
        (df_merged['inter'] < df_merged['chr_bins_end'])
    ]

    pass


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
        while True:
            if not any(df['start'].diff() < 10000):
                break
            df['end'] = df['start'].shift(-1) - 1
            df = df.iloc[:-1, :]
            df['size'] = df['end'] - df['start'] + 1
            df_filtered = df.loc[df['size'] >= 10000]
            df = df_filtered[['chr', 'start', 'end', 'size']]
            df = pd.concat(
                (df, pd.DataFrame({'chr': [chrom], 'start': [df_filtered.iloc[-1, 2]+1]})), axis=0, ignore_index=True)
        df_inter_peaks_10kb = pd.concat((df_inter_peaks_10kb, df))
    del df_filtered, df

    df_inter_peaks_10kb.dropna(axis=0, inplace=True)
    df_inter_peaks_10kb.index = range(len(df_inter_peaks_10kb))
    df_inter_peaks_10kb[['start', 'end', 'size']] = df_inter_peaks_10kb[['start', 'end', 'size']].astype('int64')
    df_inter_peaks_10kb['bin_start'] = (df_inter_peaks_10kb['start'] // 1000) * 1000
    df_inter_peaks_10kb['bin_end'] = (df_inter_peaks_10kb['end'] // 1000) * 1000
    df_inter_peaks_10kb['bin_diff'] = df_inter_peaks_10kb['bin_end'] - df_inter_peaks_10kb['bin_start']

    rebin = df_inter_peaks_10kb['bin_diff'].min() // 1000 + 1
    new_df = pd.DataFrame()
    for index, row in df_inter_peaks_10kb.iterrows():
        interval_linspace = np.linspace(row['bin_start'], row['bin_end'], rebin).astype('int64')
        interval_size = int(np.diff(interval_linspace).mean())
        duplicated_rows = pd.concat([row] * len(interval_linspace), axis=1).T.reset_index(drop=True)
        duplicated_rows['inter'] = interval_linspace
        duplicated_rows['inter_size'] = interval_size
        duplicated_rows['inter_shifted'] = (duplicated_rows['inter'] - row['bin_start']) // interval_size

        new_df = pd.concat([new_df, duplicated_rows], ignore_index=False)

    new_df.drop(columns=['bin_start', 'bin_end', 'bin_diff', 'size'], inplace=True)
    df_probes2frag = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)
    df_centro_pos = pd.read_csv(centromeres_positions, sep='\t', index_col=None)

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        not_binned_dir = binning_dir + sshic_dir + '1kb/'
        samples = [f for f in sorted(os.listdir(not_binned_dir)) if 'frequencies' in f]

        if not os.path.exists(cohesins_dir+sshic_dir):
            os.makedirs(cohesins_dir+sshic_dir)

        wt_df = {}
        for samp in samples:
            samp_id = re.search(r"AD\d+", samp).group()
            if samp_id in ["AD162", "AD242", "AD296", "AD300"]:
                print(samp_id)
                main(
                    df_peaks=new_df,
                    df_probes=df_probes2frag,
                    df_centro=df_centro_pos,
                    fragments_path=not_binned_dir+samp,
                    output_dir=cohesins_dir+sshic_dir
                )

