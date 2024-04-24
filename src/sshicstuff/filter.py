import os
import pandas as pd
import numpy as np

import sshicstuff.utils as utils
import sshicstuff.log as log

logger = log.logger


def oligo_correction(oligo_path):
    delim = "," if oligo_path.endswith(".csv") else "\t"
    oligo = pd.read_csv(oligo_path, sep=delim)
    column_to_keep = ['chr', 'start', 'end', 'name', 'type', 'sequence']
    oligo = oligo[column_to_keep]
    oligo.columns = [oligo.columns[i].lower() for i in range(len(oligo.columns))]
    oligo.sort_values(by=['chr', 'start'], inplace=True)
    oligo.reset_index(drop=True, inplace=True)

    return oligo


def fragments_correction(fragments_path, shift=0):
    fragments = pd.read_csv(fragments_path, sep='\t')
    fragments = pd.DataFrame({'frag': [k for k in range(len(fragments))],
                              'chr': fragments['chrom'],
                              'start': fragments['start_pos'],
                              'end': fragments['end_pos'],
                              'size': fragments['size'],
                              'gc_content': fragments['gc_content']
                              })

    fragments["frag"] = fragments["frag"] + shift
    return fragments


def sparse_mat_correction(sparse_mat_path):
    """
    Re-organizes the sparse contacts matrix file
    """
    contacts = pd.read_csv(sparse_mat_path, sep='\t', header=None)
    contacts.drop([0], inplace=True)
    contacts.reset_index(drop=True, inplace=True)
    contacts.columns = ['frag_a', 'frag_b', 'contacts']

    return contacts


def starts_match(fragments: pd.DataFrame, oligo: pd.DataFrame) -> pd.DataFrame:
    """
    Update the start positions of the oligo DataFrame based on the corresponding fragment positions.

    If the capture oligo is inside a fragment, update the start position of the oligo DataFrame with the start
    position of the fragment.

    Parameters
    ----------
    fragments : pd.DataFrame
        The fragments DataFrame.
    oligo : pd.DataFrame
        The oligo DataFrame.

    Returns
    -------
    pd.DataFrame
        The updated oligo DataFrame.
    """
    l_starts = []
    for i in range(len(oligo)):
        oligo_chr = oligo['chr'][i]
        middle = int((oligo['end'][i] - oligo['start'][i] - 1) / 2 + oligo['start'][i] - 1)
        if oligo_chr == 'chr_artificial':
            for k in reversed(range(len(fragments))):
                interval = range(fragments['start'][k], fragments['end'][k])
                fragments_chr = fragments['chr'][k]
                if middle in interval and fragments_chr == oligo_chr:
                    l_starts.append(fragments['start'][k])
                    break
        else:
            for k in range(len(fragments)):
                interval = range(fragments['start'][k], fragments['end'][k] + 1)
                fragments_chr = fragments['chr'][k]

                if middle in interval and fragments_chr == oligo_chr:
                    l_starts.append(fragments['start'][k])
                    break
    oligo['start'] = list(l_starts)
    return oligo


def oligo_fragments_joining(fragments: pd.DataFrame, oligo: pd.DataFrame) -> pd.DataFrame:
    """
    Join the oligo and fragments DataFrames, removing fragments that do not contain an oligo region.

    Updates the start and end columns with the corresponding fragment positions.

    Parameters
    ----------
    fragments : pd.DataFrame
        The fragments DataFrame.
    oligo : pd.DataFrame
        The oligo DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined oligo and fragments DataFrame.
    """
    oligo = starts_match(fragments, oligo)
    oligo.set_index(['chr', 'start'])
    oligo.pop("end")
    fragments.set_index(['chr', 'start'])
    oligo_fragments = fragments.merge(oligo, on=['chr', 'start'])
    oligo_fragments.sort_values(by=['chr', 'start'])
    return oligo_fragments


def first_join(x: str, oligo_fragments: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Join the contacts and oligo_fragments DataFrames, keeping only the rows that have their 'x' fragment
    (either 'frag_a' or 'frag_b', see contacts_correction function).

    Parameters
    ----------
    x : str
        Either 'a' or 'b', indicating whether to join on 'frag_a' or 'frag_b'.
    oligo_fragments : pd.DataFrame
        The joined oligo and fragments DataFrame.
    contacts : pd.DataFrame
        The corrected contacts DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined contacts and oligo_fragments DataFrame.
    """

    joined = contacts.merge(oligo_fragments, left_on='frag_'+x, right_on='frag', how='inner')
    return joined


def second_join(
        x: str, fragments: pd.DataFrame, oligo_fragments: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Add the fragments DataFrame information (=columns) for the y fragment after the first join
    (see first_join function). This is only for the y fragment, because the x fragments already have their
    information in the oligo_fragments DataFrame.

    Parameters
    ----------
    x : str
        Either 'a' or 'b', indicating which fragment corresponds to an oligo.
    fragments : pd.DataFrame
        The corrected fragments DataFrame.
    oligo_fragments : pd.DataFrame
        The joined oligo and fragments DataFrame.
    contacts : pd.DataFrame
        The corrected contacts DataFrame.

    Returns
    -------
    pd.DataFrame
        The joined DataFrame with added fragment information for the y fragment.
    """
    new_contacts = first_join(x, oligo_fragments, contacts)
    y = utils.frag2(x)
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
        oligo_capture_path: str,
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

    oligo_capture_path : str
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
        Force the overwriting of the oligo file even if the columns are already present.
        Default is True.

    Returns
    -------
    None
    """

    if not output_path:
        output_path = sparse_mat_path.replace(".txt", "_filtered.tsv")

    out_basedir = os.path.dirname(output_path)
    if not os.path.exists(out_basedir):
        os.makedirs(out_basedir)

    if not force and os.path.exists(output_path):
        logger.warning(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    utils.check_if_exists(sparse_mat_path)
    utils.check_if_exists(oligo_capture_path)
    utils.check_if_exists(fragments_list_path)

    utils.check_file_extension(sparse_mat_path, ".txt")
    utils.check_file_extension(oligo_capture_path, [".csv", ".tsv"])
    utils.check_file_extension(fragments_list_path, ".txt")

    df_fragments: pd.DataFrame = fragments_correction(fragments_list_path, frag_id_shift)
    df_oligo: pd.DataFrame = oligo_correction(oligo_capture_path)
    df_contacts: pd.DataFrame = sparse_mat_correction(sparse_mat_path)

    """
    Joining of the 3 dataframes
    """

    df_oligo_fragments = oligo_fragments_joining(df_fragments, df_oligo)
    df1 = second_join('a', df_fragments, df_oligo_fragments, df_contacts)
    df2 = second_join('b', df_fragments, df_oligo_fragments, df_contacts)
    df_contacts_joined = pd.concat([df1, df2])
    df_contacts_joined.drop("frag", axis=1, inplace=True)
    df_contacts_joined.sort_values(by=['frag_a', 'frag_b', 'start_a', 'start_b'], inplace=True)
    df_contacts_filtered = df_contacts_joined.convert_dtypes().reset_index(drop=True)

    df_contacts_filtered.to_csv(output_path, sep='\t', index=False)

    logger.info(f"Filtered contacts saved to {output_path}")


def onlyhic(
        sample_sparse_mat: str,
        oligo_capture_path: str,
        n_flanking_fragment: int = 2,
        output_path: str = None,
        force: bool = False
) -> None:

    """
    Create a contact sparse matrix file with the same format as the input file, but with no ssDNA.
    i.e., it removes fragments containing a probe and the N fragment up/downstream of it.

    Parameters
    ----------
    sample_sparse_mat : str
        Path to the sparse matrix file (hicstuff given output).

    oligo_capture_path : str
        Path to the oligo capture file (sshicstuff mandatory table).

    n_flanking_fragment : int
        Number of flanking fragments to remove around the probe fragment.
        Default is 2.

    output_path : str
        Path to the output file to be created.
        Default is None.

    force : bool
        Force the overwriting of the oligo file even if the columns are already present.
        Default is True.

    Returns
    -------
    None
    """

    if not output_path:
        output_path = sample_sparse_mat.replace(".txt", "_hic_only.txt")

    if not force and os.path.exists(output_path):
        logger.info(f"Output file already exists: {output_path}")
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    utils.check_if_exists(sample_sparse_mat)
    utils.check_if_exists(oligo_capture_path)
    utils.check_file_extension(sample_sparse_mat, ".txt")
    utils.check_file_extension(oligo_capture_path, [".csv", ".tsv"])

    oligo_capture_delim = "," if oligo_capture_path.endswith(".csv") else "\t"
    df_sparse_mat = pd.read_csv(sample_sparse_mat, sep='\t', header=None)
    df_oligo = pd.read_csv(oligo_capture_path, sep=oligo_capture_delim)

    df_contacts_hic_only = df_sparse_mat.copy(deep=True)

    ssdna_frag = df_oligo['fragment'].tolist()
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

    logger.info(f"Hi-C only contacts saved to {output_path}")
