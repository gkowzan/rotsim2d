"""Straightforward implementation of Felker's rotational coherence spectroscopy equations.

Units:
- Rotational constants (A, B, C) in cm-1,
- Time variables in s,
- Temperature in Kelvin.
"""
from typing import Union
import numpy as np
import scipy.constants as C
from math import isclose

k_cm = (C.k/C.h)/C.c*1e-2


def normalized(arr):
    return arr/np.max(np.abs(arr))


def est_jmax(B, T):
    return np.ceil(np.sqrt(25*k_cm*T/B)).astype(np.int)


def alpha(Ji: int, Ki: int, parallel: bool=True):
    F = 2 if parallel else -1
    return (2*Ji+1)/9 + F/45*( ( (Ji*(Ji-1) - 3*Ki**2)*(Ji**2-Ki**2) )/( Ji**2*(2*Ji+1) ) +
                               ( (3*Ki**2-Ji*(Ji+1) )*(2*Ji+1)*Ki**2 )/( Ji**2*(2*Ji+1)**2 ) +
                               ( ((Ji+1)*(Ji+2)-3*Ki**2)*((Ji+1)**2-Ki**2) )/( (Ji+1)**2*(2*Ji+1) ) )


def beta(Ji: int, Ki: int, parallel: bool=True):
    F = 2 if parallel else -1
    return ( 2*F*Ki**2*(Ji**2 - Ki**2) )/( 15*Ji**2*(Ji+1) )


def gamma(Ji: int, Ki: int, parallel: bool=True):
    F = 2 if parallel else -1
    return ( 2*F*Ki**2*((Ji+1)**2-Ki**2) )/( 15*Ji*(Ji+1)**2 )


def delta(Ji: int, Ki: int, parallel: bool=True):
    F = 2 if parallel else -1
    return ( 2*F*((Ji+1)**2 - Ki**2)*(Ji**2-Ki**2) )/(15*Ji*(Ji+1)*(2*Ji+1))


def decay_single(Ji: int, Ki: int, B: float, t: Union[float, np.ndarray], parallel: bool=True):
    """Fluorescence decay function without coherence/population decay.

    `B` is in cm-1 but `t` is in s, convert `B` to Hz."""
    B = B*C.c*1e2
    nu1 = 2*B*Ji
    nu2 = 2*B*(Ji+1)
    return alpha(Ji,Ki,parallel) + beta(Ji,Ki,parallel)*np.cos(2*np.pi*nu1*t) +\
        gamma(Ji,Ki,parallel)*np.cos(2*np.pi*nu2*t) +\
        delta(Ji,Ki,parallel)*np.cos(2*np.pi*(nu1+nu2)*t)


def decay_total(Jmax: int, A: float, B: float, T: float, t: Union[float, np.ndarray], parallel: bool=True):
    """Thermally-averaged fluorescence decay."""
    t = np.atleast_1d(t)
    ret = np.zeros(t.shape)
    if not isclose(A, 0.0):     # symmetric top case
        for j in range(1, Jmax+1):
            ks = np.arange(-j, j+1)
            # print(j, (j*(j+1)*B + (A-B)*ks[:, None].T**2))
            ret[...] += np.sum(decay_single(j, ks[:, None], B, t[None, :], parallel)*
                               np.exp(-(j*(j+1)*B + (A-B)*ks[:, None]**2)/k_cm/T), axis=0)
    else:                       # linear case
        js = np.arange(1, Jmax+1)
        ret[...] += np.sum(decay_single(js[:, None], 0, B, t[None, :])*\
                           np.exp(-js[:, None]*(js[:, None]+1)*B/k_cm/T), axis=0)

    return ret


def decay_total_anisotropy(Jmax: int, A: float, B: float, T: float, t: Union[float, np.ndarray]):
    """Anistropy of thermally-averaged fluorescence decay."""
    par = decay_total(Jmax, A, B, T, t, parallel=True)
    perp = decay_total(Jmax, A, B, T, t, parallel=False)

    return (par-perp)/(par+2*perp)
