#! /usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

if __name__ == "__main__":
    data_dir = os.path.dirname(os.getcwd()) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    binning_dir = outputs_dir + "binned/"
    rdna_dir = outputs_dir + "rdna/"

    rdna_regions = pd.DataFrame({'chr': ['chr12'] * 17, 'chr_bins': np.arange(451000, 468000, 1000)})
    rdna_flanking_left = pd.DataFrame({'chr': ['chr12'] * 6, 'chr_bins': np.arange(445000, 451000, 1000)})
    rdna_flanking_right = pd.DataFrame({'chr': ['chr12'] * 6, 'chr_bins': np.arange(470000, 476000, 1000)})

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

    print("Look at contacts between probes and rDNA regions")
    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        samples = sorted([f for f in os.listdir(binning_dir+sshic_dir+'1kb/') if 'frequencies.tsv' in f])
        if not os.path.exists(rdna_dir+sshic_dir):
            os.makedirs(rdna_dir+sshic_dir)
        for samp in samples:
            samp_id = re.search(r"AD\d+", binning_dir+sshic_dir+'1kb/'+samp).group()
            print(samp_id)
            df_contacts = pd.read_csv(binning_dir+sshic_dir+'1kb/'+samp, sep='\t')
            for colname, colfrag in probes_averages.items():
                df_contacts[colname] = df_contacts[colfrag].mean(axis=1)

            df_rdna_focus = pd.DataFrame()
            df_rdna_focus['avg_flanking_left \n (chr12-445000:450000)'] = \
                pd.merge(rdna_flanking_left, df_contacts, on=['chr', 'chr_bins']).mean(axis=0, numeric_only=True)

            df_rdna_focus['avg_rdna_region \n (chr12-451000:468000)'] = \
                pd.merge(rdna_regions, df_contacts, on=['chr', 'chr_bins']).mean(axis=0, numeric_only=True)

            df_rdna_focus['avg_flanking_right \n (chr12-470000:475000)'] = \
                pd.merge(rdna_flanking_right, df_contacts, on=['chr', 'chr_bins']).mean(axis=0, numeric_only=True)

            df_rdna_focus.drop(['chr_bins', 'genome_bins'], inplace=True)

            df_rdna_focus.to_csv(rdna_dir+sshic_dir+samp_id+'_rDNA_regions_frequencies.tsv', sep='\t')

    print('--DONE--')

