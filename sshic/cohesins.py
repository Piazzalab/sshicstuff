#! /usr/bin/env python3
import numpy as np
import pandas as pd
import os
import re

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
    bin_size = int(df_contacts.chr_bins[1]-df_contacts.chr_bins[0])
    df_contacts = df_contacts.loc[:, (df_contacts.columns.isin(fragments_of_interest)) |
                                     (df_contacts.columns.isin(['chr', 'chr_bins', 'sizes']))]
    df_contacts.insert(2, 'end', df_contacts.chr_bins+bin_size)
    df_contacts["sum"] = df_contacts[fragments_of_interest].sum(axis=1)
    df_contacts.drop(columns=fragments_of_interest, inplace=True)

    df_merged1 = df_contacts.merge(df_centro, on='chr')
    df_merged_cen = df_merged1[
        (df_merged1.chr_bins > (df_merged1.left_arm_length-80000-bin_size)) &
        (df_merged1.chr_bins < (df_merged1.left_arm_length+80000))
    ]

    df_contacts2 = df_merged_cen.drop(columns=["length", "left_arm_length", "right_arm_length"])
    df_contacts2.index = range(len(df_contacts2))

    df_merged2 = df_contacts2.merge(df_peaks, on='chr')
    df_merged_cohesins_areas = df_merged2[
        (df_merged2.start >= df_merged2.chr_bins) &
        (df_merged2.start + 1 <= df_merged2.end)
    ]
    df_merged_cohesins_areas.drop(columns=['end', 'uid'], axis=1, inplace=True)
    df_merged_cohesins_areas.rename(columns={'start': 'cohesin'}, inplace=True)

    df_enrichment_cohesins = df_merged_cohesins_areas[['chr', 'chr_bins', 'cohesin', 'score']]
    df_enrichment_cohesins['enrichment'] = np.nan
    for index, row in df_merged_cohesins_areas.iterrows():
        cohesins_bin_contact_per_bp = row['sum'] / 1000
        df_extended_region = df_contacts2.loc[
            (df_contacts2['chr_bins'] >= row['chr_bins'] - 3000) &
            (df_contacts2['end'] <= row['chr_bins'] + 4000) &
            (df_contacts2['chr'] == row['chr'])
        ]

        extended_region_contacts_per_bp = df_extended_region['sum'].sum() / 7000
        df_enrichment_cohesins.loc[index, 'enrichment'] = cohesins_bin_contact_per_bp / extended_region_contacts_per_bp

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
                    df_peaks=df_cohesins_peaks,
                    df_probes=df_probes2frag,
                    df_centro=df_centro_pos,
                    fragments_path=not_binned_dir+samp,
                    output_dir=cohesins_dir+sshic_dir
                )

