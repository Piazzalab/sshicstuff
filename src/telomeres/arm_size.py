#! /usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import multiprocessing as mp
from typing import Optional
from universal.utils import is_debug, sort_by_chr

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def main(
        binned_contacts_path: str,
        telomeres_size: int,
        output_dir: str
):

    samp_id = re.search(r"AD\d+", binned_contacts_path).group()
    df_contacts = pd.read_csv(binned_contacts_path, sep='\t')
    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
    fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_contacts.loc[df_contacts['chr'] == probe_chr, f] = np.nan

    #   Inter normalization
    df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
    df_merged = pd.merge(df_contacts, df_telo, on='chr')

    df_merged_telos_areas_part_a = df_merged[
        df_merged.chr_bins < (df_merged.telo_l + telomeres_size + 1000)
    ]

    df_merged_telos_areas_part_b = df_merged[
        df_merged.chr_bins > (df_merged.telo_r - telomeres_size - 1000)
    ]

    df_merged_telos_areas = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
    pass


if __name__ == "__main__":

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    binning_dir = outputs_dir + "binned/"
    pondered_dir = outputs_dir + "pondered/"
    arm_length_dir = outputs_dir + "telomeres/"
    chr_arm = inputs_dir + "S288c_chr_arm_sorted_by_length.csv"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    telomeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"

    parallel = True
    if is_debug():
        parallel = False

    excluded_chr = ['chr2', 'chr3', '2_micron', 'mitochondrion', 'chr_artificial']
    df_probes = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)
    df_chr_arm = pd.read_csv(chr_arm, sep='\t', header=None)
    df_chr_arm.columns = ['chr', 'arm', 'size', 'category']
    df_centro = pd.read_csv(telomeres_positions, sep='\t', index_col=None)
    df_telo = pd.DataFrame({'chr': df_centro['chr'], 'telo_l': 0, 'telo_r': df_centro['length']})
    df_telo = df_telo[~df_telo['chr'].isin(excluded_chr)]


    print('aggregated on telomeres positions')
    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        print('\n')
        print('raw binned tables')
        samples_not_pondered = \
            sorted([f for f in os.listdir(binning_dir+sshic_dir+'1kb/') if 'frequencies.tsv' in f])
        for samp in samples_not_pondered:
            main(
                binned_contacts_path=binning_dir+sshic_dir+'1kb/'+samp,
                telomeres_size=30000,
                output_dir=arm_length_dir
            )

    print('-- DONE --')
