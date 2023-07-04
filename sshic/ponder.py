import re
import os
import pandas as pd
import numpy as np
from typing import Optional
from utils import make_groups_of_probes


def ponder_mutant(
        statistics_path: str,
        wt_ref_name: str,
        contacts_path: str,
        frequencies_path: str,
        binned_type: str,
        output_dir: str,
        additional_path: Optional[str] = None
):

    """
    This function allows to ponder / normalize every sample that are a mutant (present in
    the dict samples_vs_wt) contacts by the normalized capture efficiency for each probe
    by using the newly made statistics table.

    Do it for each bin.

    ARGUMENTS
    ________________

    """

    sample_id = re.search(r"AD\d+", contacts_path.split("/")[-1]).group()
    output_path = os.path.join(output_dir, sample_id)

    df_stats: pd.DataFrame = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)
    df_stats["fragment"] = df_stats["fragment"].astype(str)
    df_contacts: pd.DataFrame = pd.read_csv(contacts_path, header=0, sep="\t")
    df_frequencies: pd.DataFrame = pd.read_csv(frequencies_path, header=0, sep="\t")

    probes = df_stats['probe'].tolist()
    fragments = df_stats['fragment'].astype(str).tolist()

    wt_colname = f"capture_efficiency_vs_{wt_ref_name}"
    for frag in fragments:
        ponder_coefficient = df_stats.loc[df_stats["fragment"] == frag, wt_colname].tolist()[0]
        if ponder_coefficient == np.nan:
            ponder_coefficient = 0
        df_contacts.loc[:, frag] = df_contacts.loc[:, frag] * ponder_coefficient
        df_frequencies.loc[:, frag] = df_frequencies.loc[:, frag] * ponder_coefficient

    if additional_path:
        df_additional: pd.DataFrame = pd.read_csv(additional_path, sep='\t')
        probes_to_fragments = dict(zip(probes, fragments))
        make_groups_of_probes(df_additional, df_contacts, probes_to_fragments)
        make_groups_of_probes(df_additional, df_frequencies, probes_to_fragments)

    df_contacts.to_csv(output_path+f"_{binned_type}_pondered_contacts.tsv", sep='\t', index=False)
    df_frequencies.to_csv(output_path + f"_{binned_type}_pondered_frequencies.tsv", sep='\t', index=False)

