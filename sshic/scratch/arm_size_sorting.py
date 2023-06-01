#! /usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def chr_arm(
        df_chr_arm_sorted: pd.DataFrame,
        binned_contacts_path: str,
        telomeres_size: int,
        output_dir: str
):

    samp_id = re.search(r"AD\d+", binned_contacts_path).group()
    df_contacts = pd.read_csv(binned_contacts_path, sep='\t')
    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)].drop(columns=['genome_bins'])
    fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_contacts.loc[df_contacts['chr'] == probe_chr, f] = np.nan

    #   Inter normalization
    df_contacts[fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
    #   add additional subset of probes where contacts are averaged
    for colname, colfrag in probes_averages.items():
        df_contacts[colname] = df_contacts[colfrag].mean(axis=1)

    df_merged = pd.merge(df_contacts, df_telo, on='chr')
    df_merged_telos_areas_part_a = df_merged[df_merged.chr_bins < (df_merged.telo_l + telomeres_size + 1000)]
    df_merged_telos_areas_part_a.insert(2, 'arm', 'left')
    df_merged_telos_areas_part_b = df_merged[df_merged.chr_bins > (df_merged.telo_r - telomeres_size - 1000)]
    df_merged_telos_areas_part_b.insert(2, 'arm', 'right')

    df_telo_freq = pd.concat((df_merged_telos_areas_part_a, df_merged_telos_areas_part_b))
    df_merged2 = pd.merge(df_telo_freq, df_chr_arm_sorted, on=['chr', 'arm'])
    df_merged2.drop(columns=['telo_l', 'telo_r', 'size'], inplace=True)

    df_grouped = df_merged2.groupby(by='category', as_index=False).mean(numeric_only=True).drop(columns=['chr_bins'])
    df_grouped = df_grouped.rename(columns={'category': 'fragments'}).T
    df_grouped.to_csv(output_dir+samp_id+'_arm_length_telomeres_frequencies.tsv', sep='\t', header=False)
    print(samp_id)


if __name__ == "__main__":

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    pondering_mode = ['not_pondered', 'pondered']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    binning_dir = outputs_dir + "binned/"
    pondered_dir = outputs_dir + "pondered/"
    arm_length_dir = outputs_dir + "arm_sizes_telo/"
    chr_arm = inputs_dir + "S288c_chr_arm_sorted_by_length.csv"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    telomeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"

    excluded_chr = ['chr2', 'chr3', '2_micron', 'mitochondrion', 'chr_artificial']
    chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                'chr10', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16']

    probes_averages = {
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

    df_probes = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)

    df_chr_arm = pd.read_csv(chr_arm, sep='\t', header=None)
    df_chr_arm.columns = ['chr', 'arm', 'size', 'category']

    df_centro = pd.read_csv(telomeres_positions, sep='\t', index_col=None)
    df_telo = pd.DataFrame({'chr': df_centro['chr'], 'telo_l': 0, 'telo_r': df_centro['length']})
    df_telo = df_telo[~df_telo['chr'].isin(excluded_chr)]

    print("look if the contacts frequencies on telomeres region is related to the chr arm's length")
    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        for ponder in pondering_mode:
            output_dir = ''
            samples_dir = ''
            samples = []
            if ponder == 'not_pondered':
                print('not pondered contacts')
                output_dir = arm_length_dir+'not_pondered/'+sshic_dir
                samples_dir = binning_dir + sshic_dir + '1kb/'
                samples = sorted([f for f in os.listdir(samples_dir) if 'frequencies' in f])

            elif ponder == 'pondered':
                print('pondered contacts')
                output_dir = arm_length_dir+'pondered/'+sshic_dir
                samples_dir = pondered_dir + sshic_dir + '1kb/'
                samples = sorted([f for f in os.listdir(samples_dir) if 'frequencies' in f])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            for samp in samples:
                main(
                    df_chr_arm_sorted=df_chr_arm,
                    binned_contacts_path=samples_dir+samp,
                    telomeres_size=30000,
                    output_dir=output_dir
                )

    print('-- DONE --')
