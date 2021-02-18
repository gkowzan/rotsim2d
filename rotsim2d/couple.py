"""Calculate interaction matrix elements between rotational states."""
from typing import Sequence, List, Union, Tuple
from numbers import Number
import numpy as np
import pywigxjpf.pywigxjpf as wig


def jm_dipole(j1: int, m1: int, j2:int, m2: int) -> float:
    """Geometric part of the element <j2,m2|mu|j1,m1>."""
    val = 0.0
    for q in (-1, 0, 1):
        val += wig.cg(j1, m2, 1, q, j2, m2)*\
            wig.cg(j1, 0, 1, 0, j2, 0)

    return val*np.sqrt((2*j1+1)/(2*j2+1))


def wigner6j0(a: int, b: int, c: int) -> float:
    """Value of Wigner6j(a,b,c;0,c,b)"""
    s = a+b+c
    return (-1)**s/np.sqrt((2*b+1)*(2*c+1))


def G(ji: int, jj: int, jk: int, jl: int, k: int):
    """Wigner6j part of four-fold reduced matrix element."""
    return (2*k+1)*wigner6j0(jk, ji, k)*wig.wig6j(1, 1, k, jk, ji, jj)*\
        wig.wig6j(1, 1, k, ji, jk, jl)


def T00(phi: float, phj: float, phk: float, phl: float, k: int):
    """Recoupling of four collinear beams with total Q=K=0.

    Only linear polarization.
    """
    if k==0:
        return np.cos(phi-phj)*np.cos(phk-phl)/3.0
    elif k==1:
        return np.sin(phi-phk)*np.sin(phk-phl)*np.sqrt(3)/6
    elif k==2:
        return np.sqrt(5)/60*(np.cos(phi-phj-phk+phl)+np.cos(phi-phj+phk-phl)+6*np.cos(phi+phj-phk-phl))
    else:
        return 0.0

def T00_circ(phi: float, phj: float, phk: float, phl: float,
             delti: float, deltj: float, deltk: float, deltl: float, k: int):
    """Recoupling of four collinear beams with total Q=K=0.

    Possibly circular polarization.
    """
    if k==0:
        return (1/3*(np.cos(phi)*np.cos(phj)*np.cos(phk)*np.cos(phl)+\
                     np.exp(1.0j*(deltl-deltk))*\
                     np.sin(phk)*np.sin(phl)*np.cos(phi)*np.cos(phj)+\
                     np.exp(1.0j*(deltj-delti))*\
                     np.sin(phi)*np.sin(phj)*np.cos(phk)*np.cos(phl)+\
                     np.exp(1.0j*(deltj+deltl-delti-deltk))*\
                     np.sin(phi)*np.sin(phj)*np.sin(phk)*np.sin(phl)))
    elif k==1:
        return (1/6*np.sqrt(3)*\
                (np.exp(1.0j*(deltj+deltl))*np.sin(phj)*np.sin(phl)*np.cos(phi)*np.cos(phk)-\
                 np.exp(1.0j*(deltj-deltk))*np.sin(phj)*np.sin(phk)*np.cos(phi)*np.cos(phl)-\
                 np.exp(1.0j*(deltl-delti))*np.sin(phi)*np.sin(phl)*np.cos(phj)*np.cos(phk)+\
                 np.exp(-1.0j*(delti+deltk))*np.sin(phi)*np.sin(phk)*np.cos(phj)*np.cos(phl)))
    elif k==2:
        return (np.sqrt(5)*\
                (0.1*np.exp(1.0j*(deltj+deltl))*np.sin(phj)*np.sin(phl)*np.cos(phi)*np.cos(phk)+\
                 0.1*np.exp(1.0j*(deltj-deltk))*np.sin(phj)*np.sin(phk)*np.cos(phi)*np.cos(phl)+\
                 1/7.5*np.cos(phi)*np.cos(phj)*np.cos(phk)*np.cos(phl)-\
                 1/15*np.exp(1.0j*(deltl-deltk))*np.sin(phk)*np.sin(phl)*np.cos(phi)*np.cos(phj)-\
                 1/15*np.exp(1.0j*(deltj-delti))*np.sin(phi)*np.sin(phj)*np.cos(phk)*np.cos(phl)+\
                 1/7.5*np.exp(1.0j*(deltj+deltl-delti-deltk))*\
                 np.sin(phi)*np.sin(phj)*np.sin(phk)*np.sin(phl)+\
                 0.1*np.exp(1.0j*(deltl-delti))*np.sin(phi)*np.sin(phl)*np.cos(phj)*np.cos(phk)+\
                 0.1*np.exp(-1.0j*(delti+deltk))*np.sin(phi)*np.sin(phk)*np.cos(phj)*np.cos(phl)))
    else:
        return 0.0


def four_couple(js: List[int], angles: Union[List[float], List[Tuple[float]]]):
    """j values in `js` correspond to bras associated with dipole operators"""
    if isinstance(angles[0], Number):
        return four_couple_linear(js, angles)
    else:
        return four_couple_circ(js, angles)


def four_couple_linear(js: List[int], angles: List[float]):
    return G(js[0], js[1], js[2], js[3], 0)*T00(angles[0], angles[1], angles[2], angles[3], 0)+\
        G(js[0], js[1], js[2], js[3], 1)*T00(angles[0], angles[1], angles[2], angles[3], 1)+\
        G(js[0], js[1], js[2], js[3], 2)*T00(angles[0], angles[1], angles[2], angles[3], 2)


def four_couple_circ(js: List[int], angles: List[Tuple[float]]):
    """`angles` is a list of four tuples with angles."""
    return G(js[0], js[1], js[2], js[3], 0)*T00_circ(angles[0][0], angles[1][0], angles[2][0], angles[3][0],
                                                     angles[0][1], angles[1][1], angles[2][1], angles[3][1], 0)+\
        G(js[0], js[1], js[2], js[3], 1)*T00_circ(angles[0][0], angles[1][0], angles[2][0], angles[3][0],
                                                  angles[0][1], angles[1][1], angles[2][1], angles[3][1], 1)+\
        G(js[0], js[1], js[2], js[3], 2)*T00_circ(angles[0][0], angles[1][0], angles[2][0], angles[3][0],
                                                  angles[0][1], angles[1][1], angles[2][1], angles[3][1], 2)
