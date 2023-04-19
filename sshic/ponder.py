import re
import os
import sys
import getopt
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


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    try:
        opts, args = getopt.getopt(
            argv,
            "hs:b:", [
                "--help",
                "--sample_dir",
                "--binning"]
        )

    except getopt.GetoptError:
        print('Ponder a mutant sample by an efficiency coefficient over the wt:\n'
              '-s <sample_directory> (contained binned and unbinned contacts/frequencies tsv) \n'
              '-b <binning size> (in bp) \n'
              )
        sys.exit(2)

    sample_dir_input = ""
    binning = 0
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('Ponder a mutant sample by an efficiency coefficient over the wt:\n'
                  '-s <sample_directory> (contained binned and unbinned contacts/frequencies tsv) \n'
                  '-b <binning size> (in bp) \n'
                  )
            sys.exit()
        elif opt in ("-s", "--sample_dir"):
            sample_dir_input = arg
        elif opt in ("-b", "--binning"):
            binning = int(arg)

    ponder_mutant(
        sample_dir=sample_dir_input,
        bins_size=binning
    )


if __name__ == "__main__":
    main(sys.argv[1:])

