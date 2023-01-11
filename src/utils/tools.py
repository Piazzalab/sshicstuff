import sys
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
