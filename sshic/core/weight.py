import re
import os
import pandas as pd
import numpy as np
from typing import Optional
from utils import make_groups_of_probes


def weight_mutant(
        sample_name: str,
        statistics_path: str,
        wt_ref_name: str,
        contacts_path: str,
        frequencies_path: str,
        binned_type: str,
        output_dir: str,
        additional_path: Optional[str] = None
):

    """
    This function allows to weight / normalize every sample that are a mutant (present in
    the dict samples_vs_wt) contacts by the normalized capture efficiency for each probe
    by using the newly made statistics table.

    Do it for each bin.

    ARGUMENTS
    ________________

    """

    df_stats: pd.DataFrame = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)
    df_stats["fragment"] = df_stats["fragment"].astype(str)
    df_contacts: pd.DataFrame = pd.read_csv(contacts_path, header=0, sep="\t")
    df_frequencies: pd.DataFrame = pd.read_csv(frequencies_path, header=0, sep="\t")

    probes = df_stats['probe'].tolist()
    fragments = df_stats['fragment'].astype(str).tolist()
    if additional_path:
        df_additional: pd.DataFrame = pd.read_csv(additional_path, sep='\t')
        groups = df_additional['name'].to_list()
        df_contacts.drop(columns=groups, inplace=True)
        df_frequencies.drop(columns=groups, inplace=True)
    else:
        df_additional: pd.DataFrame = pd.DataFrame()

    wt_colname = f"capture_efficiency_vs_{wt_ref_name}"
    for frag in fragments:
        weight_coefficient = df_stats.loc[df_stats["fragment"] == frag, wt_colname].tolist()[0]
        if weight_coefficient == np.nan:
            weight_coefficient = 0
        df_contacts.loc[:, frag] = df_contacts.loc[:, frag] * weight_coefficient
        df_frequencies.loc[:, frag] = df_frequencies.loc[:, frag] * weight_coefficient

    if additional_path:
        probes_to_fragments = dict(zip(probes, fragments))
        make_groups_of_probes(df_additional, df_contacts, probes_to_fragments)
        make_groups_of_probes(df_additional, df_frequencies, probes_to_fragments)

    output_path = os.path.join(output_dir, f"{sample_name}_{binned_type}_vs_{wt_ref_name}")
    df_contacts.to_csv(f"{output_path}_contacts.tsv", sep='\t', index=False)
    df_frequencies.to_csv(f"{output_path}_frequencies.tsv", sep='\t', index=False)
