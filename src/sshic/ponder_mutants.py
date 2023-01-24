import numpy as np
import pandas as pd
import os


def run(
        probes_to_fragments_path: str,
        binned_contacts_path: str,
        statistics_path: str,
        output_dir: str):

    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0).T
    df_binned_contacts = pd.read_csv(binned_contacts_path, header=0, sep="\t", index_col=0)
    df_stats = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)


    dsDNA_correct_factor: dict = {}

    pass