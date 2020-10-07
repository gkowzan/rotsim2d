"""Calculate interaction matrix elements between rotational states."""
import numpy as np
import pywigxjpf.pywigxjpf as wig


def jm_dipole(j1: int, m1: int, j2:int, m2: int) -> float:
    """Geometric part of the element <j2,m2|mu|j1,m1>."""
    val = 0.0
    for q in (-1, 0, 1):
        val += wig.cg(j1, m2, 1, q, j2, m2)*\
            wig.cg(j1, 0, 1, 0, j2, 0)

    return val*np.sqrt((2*j1+1)/(2*j2+1))
