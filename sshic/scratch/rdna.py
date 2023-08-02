#! /usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def main(
        binned_contacts_path: str,
        output_dir: str
):

    samp_id = re.search(r"AD\d+", binned_contacts_path).group()
    print(samp_id)
    df_contacts = pd.read_csv(binned_contacts_path, sep='\t')
    fragments = pd.unique(df_probes['frag_id'].astype(str))

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        df_contacts.loc[df_contacts['chr'] == probe_chr, f] = np.nan
    #   Inter normalization
    df_contacts[fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))

    for colname, colfrag in probes_averages.items():
        df_contacts[colname] = df_contacts[colfrag].mean(axis=1)

    df_rdna_focus = pd.DataFrame()
    df_rdna_focus['avg_flanking_left_(chr12-440000:450000)'] = \
        pd.merge(rdna_flanking_left, df_contacts, on=['chr', 'chr_bins']).mean(axis=0, numeric_only=True)

    df_rdna_focus['avg_rdna_region_(chr12-451000:467000)'] = \
        pd.merge(rdna_regions, df_contacts, on=['chr', 'chr_bins']).mean(axis=0, numeric_only=True)

    df_rdna_focus['avg_flanking_right_(chr12-490000:500000)'] = \
        pd.merge(rdna_flanking_right, df_contacts, on=['chr', 'chr_bins']).mean(axis=0, numeric_only=True)

    df_rdna_focus.drop(['chr_bins', 'genome_bins'], inplace=True)
    df_rdna_focus.to_csv(output_dir + samp_id + '_rDNA_regions_frequencies.tsv', sep='\t')

    del df_contacts, df_rdna_focus


if __name__ == "__main__":
    data_dir = os.path.dirname(os.getcwd()) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    weighting_mode = ['not_weighted', 'weighted']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    binning_dir = outputs_dir + "binned/"
    weighted_dir = outputs_dir + "weighted/"
    rdna_dir = outputs_dir + "rdna/"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"

    df_probes = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)
    rdna_regions = pd.DataFrame({'chr': ['chr12'] * 17, 'chr_bins': np.arange(451000, 468000, 1000)})
    rdna_flanking_left = pd.DataFrame({'chr': ['chr12'] * 11, 'chr_bins': np.arange(440000, 451000, 1000)})
    rdna_flanking_right = pd.DataFrame({'chr': ['chr12'] * 10, 'chr_bins': np.arange(490000, 500000, 1000)})

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
        for weight in weighting_mode:
            output_dir = ''
            samples_dir = ''
            samples = []
            if weight == 'not_weighted':
                print('not weighted contacts')
                output_dir = rdna_dir+'not_weighted/'+sshic_dir
                samples_dir = binning_dir + sshic_dir + '1kb/'
                samples = sorted([f for f in os.listdir(samples_dir) if 'frequencies' in f])

            elif weight == 'weighted':
                print('weighted contacts')
                output_dir = rdna_dir+'weighted/'+sshic_dir
                samples_dir = weighted_dir + sshic_dir + '1kb/'
                samples = sorted([f for f in os.listdir(samples_dir) if 'frequencies' in f])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            for samp in samples:
                main(
                    binned_contacts_path=samples_dir+samp,
                    output_dir=output_dir
                )

    print('--DONE--')
