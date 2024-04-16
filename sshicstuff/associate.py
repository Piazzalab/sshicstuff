import logging
import pandas as pd
import numpy as np

import sshicstuff.utils as sshcu

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')



def associate_oligo_to_frag(
        oligos_capture_path: str,
        fragments_path: str,
        frag_id_shift: int = 0,
        force: bool = True
):
    """
    Associate oligos to fragments based on the fragment name.
    It adds 3 columns directly at the end of the oligos file :
    - fragment : id of the fragment, from fragment_list (hicstuff file output)
    - fragment_start : start position of the fragment
    - fragment_end : end position of the fragment

    Parameters
    ----------
    oligos_capture_path : str
        Path to the .csv file containing the oligos.
    fragments_path : str
        Path to the .csv file containing the fragments.
    frag_id_shift : int
        Shift to apply to the fragment ids.
    force : bool
        If True, the function will overwrite the oligos file.
    Returns
    -------
    None
    """
    
    logging.info("Associating oligos to fragments based on the fragment id, start and end positions.")
    
    sshcu.check_if_exists(oligos_capture_path)
    sshcu.check_if_exists(fragments_path)
    sshcu.check_file_extension(fragments_path, ".txt")
    sshcu.check_file_extension(oligos_capture_path, [".csv", ".tsv", ".txt"])

    # Read the oligos and fragments files
    oligos_delim = "," if oligos_capture_path.endswith(".csv") else "\t"
    df_oligos = pd.read_csv(oligos_capture_path, sep=oligos_delim)

    if "fragment" in df_oligos.columns and not force:
        logging.info("Oligos already associated to fragments. Use --force=True to overwrite.")
        return

    df_fragments = pd.read_csv(fragments_path, sep='\t')
    df_fragments['frag'] = [k for k in range(len(df_fragments))]
    df_fragments["frag"] = df_fragments["frag"] + frag_id_shift

    fragments_id = []
    fragments_start = []
    fragments_end = []
    for index, row in df_oligos.iterrows():
        (chr_, probe_start, probe_end, probe_chr_ori, probe_start_ori,
         probe_end_ori, probe_type, probe, probe_seq) = row[:9]
        df_sub_fragments = df_fragments[df_fragments['chrom'] == chr_]
        df_sub_fragment_sorted_start = np.sort(df_sub_fragments['start_pos'].to_numpy())

        probe_middle = int(probe_start + (probe_end - probe_start) / 2)

        idx = np.searchsorted(df_sub_fragment_sorted_start, probe_middle, side="left")
        nearest_frag_start = df_sub_fragment_sorted_start[idx - 1]

        frag_id = df_sub_fragments.index[df_sub_fragments['start_pos'] == nearest_frag_start].tolist()[0]
        frag_start = df_sub_fragments.loc[frag_id, 'start_pos']
        frag_end = df_sub_fragments.loc[frag_id, 'end_pos']
        fragments_id.append(frag_id)
        fragments_start.append(frag_start)
        fragments_end.append(frag_end)

    df_oligos['fragment'] = fragments_id
    df_oligos['fragment_start'] = fragments_start
    df_oligos['fragment_end'] = fragments_end
    df_oligos.to_csv(oligos_capture_path, sep=",", index=False)

    logging.info("Oligos associated to fragments successfully.")

    """
    Example of usage:

    python3 ./main.py associate \
    ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt \
    ../data/inputs/capture_oligo_positions.csv \
    -F
    """