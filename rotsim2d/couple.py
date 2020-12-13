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


def wigner6j0(a: int, b: int, c: int) -> float:
    """Value of Wigner6j(a,b,c;0,c,b)"""
    s = a+b+c
    return (-1)**s/np.sqrt((2*b+1)*(2*c+1))


def G(ji, j1, j2, j3, k):
    """Wigner6j part of four-fold reduced matrix element."""
    return (2*k+1)*wigner6j0(j2, ji, k)*wig.wig6j(1, 1, k, j2, ji, j1)*\
        wig.wig6j(1, 1, k, ji, j2, j3)


def T00(phi, phj, phk, phl, k):
    """Recoupling of four collinear beams with total Q=K=0."""
    if k==0:
        return np.cos(phi-phj)*np.cos(phk-phl)/3.0
    elif k==1:
        return np.sin(phi-phk)*np.sin(phk-phl)*np.sqrt(3)/6
    elif k==2:
        return np.sqrt(5)/60*(np.cos(phi-phj-phk+phl)+np.cos(phi-phj+phk-phl)+6*np.cos(phi+phj-phk-phl))
    else:
        return 0.0


def four_couple(js, angles):
    return G(js[0], js[1], js[2], js[3], 0)*T00(angles[0], angles[1], angles[2], angles[3], 0)+\
        G(js[0], js[1], js[2], js[3], 1)*T00(angles[0], angles[1], angles[2], angles[3], 1)+\
        G(js[0], js[1], js[2], js[3], 2)*T00(angles[0], angles[1], angles[2], angles[3], 2)
