import re
import os
import sys
import argparse
import pandas as pd
import numpy as np
from typing import Optional


def ponder_mutant(
        statistics_path: str,
        contacts_path: str,
        frequencies_path: str,
        binned_type: str,
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
    output_dir = os.path.dirname(contacts_path)
    output_path = os.path.join(output_dir, sample_id)

    df_stats: pd.DataFrame = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)
    df_stats["fragment"] = df_stats["fragment"].astype(str)
    df_contacts: pd.DataFrame = pd.read_csv(contacts_path, header=0, sep="\t")
    df_frequencies: pd.DataFrame = pd.read_csv(frequencies_path, header=0, sep="\t")

    probes = df_stats['probe'].tolist()
    fragments = df_stats['fragment'].astype(str).tolist()

    for frag in fragments:
        ponder_coefficient = df_stats.loc[df_stats["fragment"] == frag, "capture_efficiency_vs_wt"].tolist()[0]
        if ponder_coefficient == np.nan:
            ponder_coefficient = 0
        df_contacts.loc[:, frag] = df_contacts.loc[:, frag] * ponder_coefficient
        df_frequencies.loc[:, frag] = df_frequencies.loc[:, frag] * ponder_coefficient

    if additional_path:
        probe_to_frag = dict(zip(probes, fragments))
        df_additional: pd.DataFrame = pd.read_csv(additional_path, sep='\t')
        for index, row in df_additional.iterrows():
            group_probes = row["probes"].split(",")
            group_frags = pd.unique([probe_to_frag[probe] for probe in group_probes])
            group_name = row["name"]
            if row["action"] == "average":
                df_contacts[group_name] = df_contacts[group_frags].mean(axis=1)
            elif row["action"] == "sum":
                df_contacts[group_name] = df_contacts[group_frags].sum(axis=1)
            else:
                continue
            df_frequencies[group_name] = df_contacts[group_name].div(sum(df_contacts[group_name]))

    df_contacts.to_csv(output_path+f"_{binned_type}_pondered_contacts.tsv", sep='\t', index=False)
    df_frequencies.to_csv(output_path + f"_{binned_type}_pondered_frequencies.tsv", sep='\t', index=False)


def main(argv):
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    parser = argparse.ArgumentParser(
        description='Ponder a mutant sample by an efficiency coefficient over the wt')
    parser.add_argument('-s', '--sample-dir', type=str, required=True,
                        help='Path to the sample_directory (containing binned and unbinned contacts/frequencies tsv)')

    parser.add_argument('-b', '--binning', type=int, required=True, help='Binning size (in bp)')

    parser.add_argument('-a', '--additional', type=str, required=False,
                        help='Path to additional groups of probes table')

    args = parser.parse_args(argv)


if __name__ == "__main__":
    main(sys.argv[1:])

