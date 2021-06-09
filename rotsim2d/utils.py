import collections
import numpy as np
import scipy.constants as C


def wn2nu(wn):
    return wn*C.c*1e2


def nu2wn(nu):
    return nu/C.c*1e-2


def find_index(arr, val, axis=None):
    """Returns index of an `arr` value that is closest to `val`."""
    return np.argmin(np.abs(arr-val), axis=axis)


def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d
