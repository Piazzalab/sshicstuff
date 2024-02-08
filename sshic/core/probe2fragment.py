import numpy as np
import pandas as pd
from utils import find_nearest


def associate_probes_to_fragments(fragments_list_path: str, oligos_capture_path: str):
    """
    This function aims at formatting and creating a correspondence between each probe (oligo capture) and
    the fragments aka the read that contains it.
    Complementary information are given such as chromosomes of the probe (basically the same as the fragment),
    the position on the chromosome of the fragment and the probe, and the type (ss, ds, ds_neg etc ..) of the probe

    The resulting dataframe is written in a tsv file, in the same location of that of the fragments_list.txt

    ARGUMENTS
    _______________________
    fragments_list_path : str
        path to the digested fragments list based on a restriction enzyme or a fixed chunk size.
    oligos_capture_path : str
        path to the file containing the oligo-nucleotides capture information
    """

    df_fragments = pd.read_csv(fragments_list_path, sep='\t')
    df_oligos = pd.read_csv(oligos_capture_path, sep=",")
    if "fragment" in df_oligos.columns:
        return
    else:
        fragment_id = []
        fragment_start = []
        fragment_end = []

        for index, row in df_oligos.iterrows():
            chrom, probe_start, probe_end, probe_type, probe, probe_seq = row
            sub_df_fragments = df_fragments[df_fragments['chrom'] == chrom]
            oligo_middle = int(probe_start + (probe_end-probe_start)/2)
            nearest_frag_start = find_nearest(
                array=sub_df_fragments['start_pos'], key=oligo_middle, mode='lower'
            )

            id_ = sub_df_fragments.index[sub_df_fragments['start_pos'] == nearest_frag_start].tolist()[0]
            start = sub_df_fragments.loc[id_, 'start_pos']
            end = sub_df_fragments.loc[id_, 'end_pos']
            fragment_id.append(id_)
            fragment_start.append(start)
            fragment_end.append(end)

        df_oligos.insert(5, "fragment", np.array(fragment_id))
        df_oligos.insert(6, "fragment_start", np.array(fragment_start))
        df_oligos.insert(7, "fragment_end", np.array(fragment_end))

    df_oligos.to_csv(oligos_capture_path, sep=",", index=False)





