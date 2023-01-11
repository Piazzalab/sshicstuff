#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from utils import tools


if __name__ == "__main__":

    df_nucleosomes_raw = pd.read_csv(
        '../../data/inputs/Chereji_Henikoff_genome_research_NFR.bed', sep='\t')
    df_nucleosomes = df_nucleosomes_raw.drop_duplicates(keep='first')
    del df_nucleosomes_raw
    df_nucleosomes.index = range(len(df_nucleosomes))

    df_fragments_all = pd.read_csv(
        "../../data/outputs/binning/sshic/0kb/AD162_0kb_frequencies.tsv", sep='\t', index_col=0, low_memory=False)

    df_info, df_fragments = tools.split_formatted_dataframe(df_fragments_all)
    pass