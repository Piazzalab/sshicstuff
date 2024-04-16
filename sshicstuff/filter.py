import os
import logging
import pandas as pd
import numpy as np

import sshicstuff.utils as sshcu

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')


def starts_match(fragments: pd.DataFrame, oligos: pd.DataFrame) -> pd.DataFrame:
    """
    Update the start positions of the oligos DataFrame based on the corresponding fragment positions.

    If the capture oligo is inside a fragment, update the start position of the oligos DataFrame with the start
    position of the fragment.

    Parameters
    ----------
    fragments : pd.DataFrame
        The fragments DataFrame.
    oligos : pd.DataFrame
        The oligos DataFrame.

    Returns
    -------
    pd.DataFrame
        The updated oligos DataFrame.
    """
    l_starts = []
    for i in range(len(oligos)):
        oligos_chr = oligos['chr'][i]
        middle = int((oligos['end'][i] - oligos['start'][i] - 1) / 2 + oligos['start'][i] - 1)
        if oligos_chr == 'chr_artificial':
            for k in reversed(range(len(fragments))):
                interval = range(fragments['start'][k], fragments['end'][k])
                fragments_chr = fragments['chr'][k]
                if middle in interval and fragments_chr == oligos_chr:
                    l_starts.append(fragments['start'][k])
                    break
        else:
            for k in range(len(fragments)):
                interval = range(fragments['start'][k], fragments['end'][k] + 1)
                fragments_chr = fragments['chr'][k]

                if middle in interval and fragments_chr == oligos_chr:
                    l_starts.append(fragments['start'][k])
                    break
    oligos['start'] = list(l_starts)
    return oligos


def oligos_fragments_joining(fragments: pd.DataFrame, oligos: pd.DataFrame) -> pd.DataFrame:
    """
    Join the oligos and fragments DataFrames, removing fragments that do not contain an oligo region.

    Updates the start and end columns with the corresponding fragment positions.

    Parameters
    ----------
    fragments : pd.DataFrame
        The fragments DataFrame.
    oligos : pd.DataFrame
        The oligos DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined oligos and fragments DataFrame.
    """
    oligos = starts_match(fragments, oligos)
    oligos.set_index(['chr', 'start'])
    oligos.pop("end")
    fragments.set_index(['chr', 'start'])
    oligos_fragments = fragments.merge(oligos, on=['chr', 'start'])
    oligos_fragments.sort_values(by=['chr', 'start'])
    return oligos_fragments


def first_join(x: str, oligos_fragments: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Join the contacts and oligos_fragments DataFrames, keeping only the rows that have their 'x' fragment
    (either 'frag_a' or 'frag_b', see contacts_correction function).

    Parameters
    ----------
    x : str
        Either 'a' or 'b', indicating whether to join on 'frag_a' or 'frag_b'.
    oligos_fragments : pd.DataFrame
        The joined oligos and fragments DataFrame.
    contacts : pd.DataFrame
        The corrected contacts DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined contacts and oligos_fragments DataFrame.
    """

    joined = contacts.merge(oligos_fragments, left_on='frag_'+x, right_on='frag', how='inner')
    return joined


def second_join(
        x: str, fragments: pd.DataFrame, oligos_fragments: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Add the fragments DataFrame information (=columns) for the y fragment after the first join
    (see first_join function). This is only for the y fragment, because the x fragments already have their
    information in the oligos_fragments DataFrame.

    Parameters
    ----------
    x : str
        Either 'a' or 'b', indicating which fragment corresponds to an oligo.
    fragments : pd.DataFrame
        The corrected fragments DataFrame.
    oligos_fragments : pd.DataFrame
        The joined oligos and fragments DataFrame.
    contacts : pd.DataFrame
        The corrected contacts DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined DataFrame with added fragment information for the y fragment.
    """
    new_contacts = first_join(x, oligos_fragments, contacts)
    y = sshcu.frag2(x)
    joined = new_contacts.join(fragments.drop("frag", axis=1),
                               on='frag_'+y,
                               lsuffix='_' + x[-1],
                               rsuffix='_' + y[-1], how='left')

    # puts a suffix to know what fragment corresponds to an oligo
    joined.rename(columns={"type": "type_" + x[-1],
                           "name": "name_" + x[-1],
                           "sequence": "sequence_" + x[-1]
                           },
                  inplace=True)
    return joined


def filter_contacts(
        sparse_mat_path: str,
        oligos_capture_path: str,
        fragments_list_path: str,
        output_path: str = None,
        frag_id_shift: int = 0,
        force: bool = False
) -> None:

    """
    Filter the sparse matrix by creating a nex table that contains only pairs of fragments that have an oligo
    either on frag_a or frag_b.

    The output table will contain the following columns:
    - frag_a: fragment id of the first fragment
    - frag_b: fragment id of the second fragment
    - contacts: number of contacts between the two fragments
    - chr_a: chromosome of the first fragment
    - start_a: start position of the first fragment
    - end_a: end position of the first fragment
    - size_a: size of the first fragment
    - gc_content_a: gc content of the first fragment
    - type_a: type of the oligo on the first fragment
    - name_a: name of the oligo on the first fragment
    - sequence_a: sequence of the oligo on the first fragment
    - chr_b: chromosome of the second fragment
    - start_b: start position of the second fragment
    - end_b: end position of the second fragment
    - size_b: size of the second fragment
    - gc_content_b: gc content of the second fragment
    - type_b: type of the oligo on the second fragment
    - name_b: name of the oligo on the second fragment
    - sequence_b: sequence of the oligo on the second fragment

    Parameters
    ----------
    sparse_mat_path : str
        Path to the sparse matrix file (hicstuff given output).

    oligos_capture_path : str
        Path to the oligo capture file (sshicstuff mandatory table).

    fragments_list_path : str
        Path to the fragments list file (hicstuff given output).

    output_path : str
        Path to the output file to be created.
        Default is None.

    frag_id_shift : int
        Shift the fragment id by this value.
        Default is 0.

    force : bool
        Force the overwriting of the oligos file even if the columns are already present.
        Default is True.

    Returns
    -------
    None
    """

    if not output_path:
        output_path = sparse_mat_path.replace(".txt", "_filtered.txt")

    if not force and os.path.exists(output_path):
        logging.info(f"Output file already exists: {output_path}")
        logging.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    sshcu.check_if_exists(sparse_mat_path)
    sshcu.check_if_exists(oligos_capture_path)
    sshcu.check_if_exists(fragments_list_path)

    sshcu.check_file_extension(sparse_mat_path, ".txt")
    sshcu.check_file_extension(oligos_capture_path, [".csv", ".tsv"])
    sshcu.check_file_extension(fragments_list_path, ".txt")

    sparse_delim = "\t"
    oligos_capture_delim = "," if oligos_capture_path.endswith(".csv") else "\t"
    fragment_delim = "\t"

    df_fragments_raw = pd.read_csv(fragments_list_path, sep='\t')
    df_fragments = pd.DataFrame(
        {'frag': [k for k in range(len(df_fragments_raw))],
         'chr': df_fragments_raw['chrom'],
         'start': df_fragments_raw['start_pos'],
         'end': df_fragments_raw['end_pos'],
         'size': df_fragments_raw['size'],
         'gc_content': df_fragments_raw['gc_content']
         }
    )

    df_fragments["frag"] = df_fragments["frag"] + frag_id_shift
    df_oligos = pd.read_csv(oligos_capture_path, sep=",")

    """
    Contacts import and correction of col names
    """

    df_contacts = pd.read_csv(sparse_mat_path, sep='\t', header=None)
    df_contacts = df_contacts.drop([0])
    df_contacts.reset_index(drop=True, inplace=True)
    df_contacts.columns = ['frag_a', 'frag_b', 'contacts']

    """
    Joining of the 3 dataframes
    """

    df_oligos_fragments = oligos_fragments_joining(df_fragments, df_oligos)
    df1 = second_join('a', df_fragments, df_oligos_fragments, df_contacts)
    df2 = second_join('b', df_fragments, df_oligos_fragments, df_contacts)
    df_contacts_joined = pd.concat([df1, df2])
    df_contacts_joined.drop("frag", axis=1, inplace=True)
    df_contacts_joined.sort_values(by=['frag_a', 'frag_b', 'start_a', 'start_b'], inplace=True)
    df_contacts_filtered = df_contacts_joined.convert_dtypes().reset_index(drop=True)

    df_contacts_filtered.to_csv(output_path, sep='\t', index=False)


def onlyhic(
        sample_sparse_mat: str,
        oligos_capture_path: str,
        n_flanking_fragment: int = 2,
        output_path: str = None,
        force: bool = False
):

    """
    Create a contact sparse matrix file with the same format as the input file, but with no ssDNA.
    i.e., it removes fragments containing a probe and the N fragment up/downstream of it.

    Parameters
    ----------
    sample_sparse_mat : str
        Path to the sparse matrix file (hicstuff given output).

    oligos_capture_path : str
        Path to the oligo capture file (sshicstuff mandatory table).

    n_flanking_fragment : int
        Number of flanking fragments to remove around the probe fragment.
        Default is 2.

    output_path : str
        Path to the output file to be created.
        Default is None.
    """

    if not output_path:
        output_path = sample_sparse_mat.replace(".txt", "_HiC_only.txt")

    if not force and os.path.exists(output_path):
        logging.info(f"Output file already exists: {output_path}")
        logging.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    sshcu.check_if_exists(sample_sparse_mat)
    sshcu.check_if_exists(oligos_capture_path)
    sshcu.check_file_extension(sample_sparse_mat, ".txt")
    sshcu.check_file_extension(oligos_capture_path, [".csv", ".tsv"])

    sparse_delim = "\t"
    oligos_capture_delim = "," if oligos_capture_path.endswith(".csv") else "\t"
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep=sparse_delim, header=None)
    df_oligos = pd.read_csv(oligos_capture_path, sep=oligos_capture_delim)

    df_contacts_hic_only = df_sparse_mat.copy(deep=True)

    ssdna_frag = df_oligos['fragment'].tolist()
    ssdna_frag_flanking = []
    for f in ssdna_frag:
        for i in range(1, n_flanking_fragment + 1):
            ssdna_frag_flanking.append(f + i)
            ssdna_frag_flanking.append(f - i)

    ssdna_frag_all = np.unique(ssdna_frag + ssdna_frag_flanking)
    df_ssdna = pd.DataFrame(ssdna_frag_all, columns=['fragments'])

    df_sparse_mat["index"] = df_sparse_mat.index
    matches_a = pd.merge(df_sparse_mat, df_ssdna, left_on=0, right_on='fragments', how='inner', indicator=True)
    matches_b = pd.merge(df_sparse_mat, df_ssdna, left_on=1, right_on='fragments', how='inner', indicator=True)
    index_to_drop = np.unique(np.concatenate((matches_a['index'].to_numpy(), matches_b['index'].to_numpy())))

    df_contacts_hic_only.drop(index_to_drop, inplace=True)

    df_contacts_hic_only.iloc[0, 0] -= len(df_ssdna)
    df_contacts_hic_only.iloc[0, 1] -= len(df_ssdna)
    df_contacts_hic_only.iloc[0, 2] -= len(index_to_drop)

    df_contacts_hic_only.to_csv(output_path, sep='\t', index=False, header=False)

    """
    Example of usage:
    
    python3 ./main.py hiconly \
      ../data/samples/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree.txt \
      ../data/inputs/capture_oligo_positions.csv \
      -o ../data/outputs/AD241_S288c_Hic_only.txt \
      -n 2 \
      -F
    """