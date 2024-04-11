import sys
import re
import shutil
import numpy as np
import pandas as pd

from typing import List


def is_debug() -> bool:
    """
    Check if the script is running in debug mode.

    Returns
    -------
    bool
        True if the script is running in debug mode, False otherwise.
    """
    gettrace = getattr(sys, 'gettrace', None)

    if gettrace is None:
        return False
    else:
        v = gettrace()
        if v is None:
            return False
        else:
            return True


def detect_delimiter(path: str):
    """
    Detect the delimiter of a file.
    The delimiter is detected by counting the number of tabs and commas in the file.

    Parameters
    ----------
    path : str
        Path to the file.

    Returns
    -------
    str
        Delimiter of the file.
    """
    with open(path, 'r') as file:
        contents = file.read()
    tabs = contents.count('\t')
    commas = contents.count(',')
    if tabs > commas:
        return '\t'
    else:
        return ','


def frag2(x):
    """
    if x = a get b, if x = b get a

    Parameters
    ----------
    x : str
        String to invert.

    Returns
    -------
    str
        Inverted string.
    """
    if x == 'a':
        y = 'b'
    else:
        y = 'a'
    return y


def sort_by_chr(df: pd.DataFrame, chr_list: List[str], *args: str):
    """
    Sort a DataFrame by chromosome and then by other columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to sort.
    chr_list : List[str]
        List of chromosomes.
    args : str
        Columns to sort by after the chromosome.

    Returns
    -------
    pd.DataFrame
        Sorted DataFrame.
    """
    # use re to identify chromosomes of the form "chrX" with X being a number
    chr_with_number = [c for c in chr_list if re.match(r'chr\d+', c)]
    chr_with_number.sort(key=lambda x: int(x[3:]))
    chr_without_number = [c for c in chr_list if c not in chr_with_number]

    order = chr_with_number + chr_without_number
    df['chr'] = df['chr'].apply(lambda x: order.index(x) if x in order else len(order))

    if args:
        df = df.sort_values(by=['chr', *args])
    else:
        df = df.sort_values(by=['chr'])

    df['chr'] = df['chr'].map(lambda x: order[x])
    df.index = range(len(df))

    return df


def make_groups_of_probes(df_groups: pd.DataFrame, df: pd.DataFrame, prob2frag: dict):
    for index, row in df_groups.iterrows():
        group_probes = row["probes"].split(",")
        group_frags = np.unique([prob2frag[probe] for probe in group_probes])
        group_name = row["name"]
        if row["action"] == "average":
            df[group_name] = df[group_frags].mean(axis=1)
        elif row["action"] == "sum":
            df[group_name] = df[group_frags].sum(axis=1)
        else:
            continue


def copy(source_path, destination_path):
    """
    Copy a file from source to destination.
    Useful to copy inputs files to the output directory in order to have a trace.

    Parameters
    ----------
    source_path : str
        Path to the file to copy.
    destination_path : str
        Path to the destination directory.
    """
    try:
        shutil.copy(source_path, destination_path)
        print(f"File {source_path.split('/')[-1]} copied successfully.")
    except IOError as e:
        print(f"Unable to copy file. Error: {e}")


def check_if_exists(file_path: str):
    """
    Check if a file exists.

    Parameters
    ----------
    file_path : str
        Path to the file to check.

    Returns
    -------
    bool
        True if the file exists, False otherwise.
    """
    try:
        with open(file_path, 'r'):
            return True
    except FileNotFoundError:
        return False

