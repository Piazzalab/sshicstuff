import sys
import numpy as np
from typing import Optional

import pandas as pd


def is_debug() -> bool:
    """
    Function to see if the script is running in debug mode.
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


def find_nearest(array: list | np.ndarray,
                 key: int | float,
                 mode: Optional[str] = None) -> int | float:
    """
    Find the nearest element from key into a list of numerics (integer or float).
    Possibility to specify if we want the nearest element in upstream or downstream from the key.
    Optionally we can specify manually boundaries that are not in the list like '0' for the case where the key hasn't
    an element in downstream.
    """

    array = np.asarray(array)
    if mode == 'upper':
        # the smallest element of array GREATER than key
        if key >= np.max(array):
            return np.max(array)
        else:
            return array[array > key].min()
    elif mode == 'lower':
        # the largest element of array LESS than key
        if key <= np.min(array):
            return np.min(array)
        else:
            return array[array < key].max()
    else:
        return array[(np.abs(array - key)).argmin()]


def frag2(x):
    """
    if x = a get b, if x = b get a
    """
    if x == 'a':
        y = 'b'
    else:
        y = 'a'
    return y


def sort_by_chr(
        df: pd.DataFrame,
        col1: str = 'chr',
        col2: str = 'positions'):

    order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
             'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
             'chr15', 'chr16', '2_micron', 'mitochondrion', 'chr_artificial']

    df[col1] = df[col1].map(lambda x: order.index(x) if x in order else len(order))
    df = df.sort_values(by=[col1, col2])
    df[col1] = df[col1].map(lambda x: order[x])
    df.index = range(len(df))
    return df
