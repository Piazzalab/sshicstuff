import re
import os
import sys
import argparse
import pandas as pd
import numpy as np


def ponder_mutant(
        sample_dir: str,
        bins_size: int = 0
):

    """
    This function allows to ponder / normalize every sample that are a mutant (present in
    the dict samples_vs_wt) contacts by the normalized capture efficiency for each probe
    by using the newly made statistics table.

    Do it for each bin.

    ARGUMENTS
    ________________

    """

    bin_prefix = str(bins_size // 1000) + "kb"
    sample_id = re.search(r"AD\d+", sample_dir).group()
    output_path = os.path.join(sample_dir, sample_id)

    df_stats: pd.DataFrame = pd.read_csv(
        os.path.join(sample_dir, sample_id+"_global_statistics.tsv"), header=0, sep="\t", index_col=0)

    if bin_prefix == "0kb":
        bin_prefix = "unbinned"
    else:
        bin_prefix += "_binned"
    df_contacts: pd.DataFrame = pd.read_csv(
        os.path.join(sample_dir, sample_id+f"_{bin_prefix}_contacts.tsv"), header=0, sep="\t")
    df_frequencies: pd.DataFrame = pd.read_csv(
        os.path.join(sample_dir, sample_id+f"_{bin_prefix}_frequencies.tsv"), header=0, sep="\t")

    fragments_list = df_stats['fragments'].unique()
    for frag in fragments_list:
        ponder_coefficient = df_stats.loc[frag, f"capture_efficiency_vs_wt"]
        if ponder_coefficient == np.nan:
            ponder_coefficient = 0
        df_contacts.loc[:, frag] = df_contacts.loc[:, frag] * ponder_coefficient
        df_frequencies.loc[:, frag] = df_frequencies.loc[:, frag] * ponder_coefficient

    df_contacts.to_csv(output_path+f"_{bin_prefix}_contacts.tsv")
    df_frequencies.to_csv(output_path + f"_{bin_prefix}_frequencies.tsv")


def main(argv):
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    parser = argparse.ArgumentParser(
        description='Ponder a mutant sample by an efficiency coefficient over the wt')
    parser.add_argument('-s', '--sample-dir', type=str, required=True,
                        help='Path to the sample_directory (containing binned and unbinned contacts/frequencies tsv)')
    parser.add_argument('-b', '--binning', type=int, required=True, help='Binning size (in bp)')

    args = parser.parse_args(argv)
    ponder_mutant(
        sample_dir=args.sample_dir,
        bins_size=args.binning
    )


if __name__ == "__main__":
    main(sys.argv[1:])

