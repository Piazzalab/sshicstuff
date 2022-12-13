import sys
import pandas as pd
import numpy as np
from typing import Optional


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


def split_formatted_dataframe(df0: pd.DataFrame):
    """
    Take a dataframe created in contacts_formats that contains within it infomations
    about each oligo (probe's name, type, location etc ...) and the number of contacts per bins.
    This function will separate the 'info' from the contacts in the df0 into a df1 and df2.
    This allows to use dtype=str for df1 (information about oligos) and dtype=int for df2
    (containing contacts per bins).
    """
    df0_a = df0.iloc[:5, :]
    df0_b = df0.iloc[5:, :]
    df1 = df0_a[[c for c in df0_a.columns if c not in ['chr', 'chr_bins', 'genome_bins', 'positions']]].astype(str)
    df2 = pd.DataFrame()
    df2['chr'] = df0_b.iloc[:, 0].astype(str)
    df2[df0_b.columns[1:]] = df0_b.iloc[:, 1:].astype(int)

    return df1, df2


def find_nearest(array: list | np.ndarray,
                 key: int | float,
                 mode: str,
                 maxi: Optional[int] = 10000,
                 mini: Optional[int] = 0,) -> int | float:
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
            return maxi
        else:
            return array[array > key].min()
    elif mode == 'lower':
        # the largest element of array LESS than key
        if key <= np.min(array):
            return mini
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

if __name__ == "__main__":
    pass
