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
        df_chr_arms: pd.DataFrame,
        binned_contacts_path: str,
        telomeres_size: int,
        output_dir: str
):
    samp_id = re.search(r"AD\d+", binned_contacts_path).group()
    df_contacts = pd.read_csv(binned_contacts_path, sep='\t')
    fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])


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

    parallel = True
    if is_debug():
        parallel = False

    df_chr_arm_sorted = pd.read_csv(chr_arm, sep='\t', header=None)
    df_chr_arm_sorted.columns = ['chr', 'arm', 'size', 'category']

    print('aggregated on telomeres positions')
    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        print('\n')
        print('raw binned tables')
        samples_not_pondered = \
            sorted([f for f in os.listdir(binning_dir+sshic_dir+'1kb/') if 'frequencies.tsv' in f])
        for samp in samples_not_pondered:
            main(
                df_chr_arms=df_chr_arm_sorted,
                binned_contacts_path=binning_dir+sshic_dir+'1kb/'+samp,
                telomeres_size=30000,
                output_dir=arm_length_dir
            )

    print('-- DONE --')
