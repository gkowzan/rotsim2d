"""Calculate interaction matrix elements between rotational states."""
import itertools as it
from numbers import Number
from typing import Callable, List, Optional, Sequence, Tuple, Union

import numpy as np
from sympy.physics.wigner import wigner_3j


def brute_geometric_factor(js, p):
    """Brute force sum over 4 dipole matrix elements without RME part for only
    one polarization component.

    This is mostly for testing.
    """
    res = 0.0
    for m0, m1, m2, m3 in it.product(
            range(-js[0], js[0]+1),
            range(-js[1], js[1]+1),
            range(-js[2], js[2]+1),
            range(-js[3], js[3]+1)):
        res += (-1)**(js[0]-m0)*wigner_3j(js[0], 1, js[1], -m0, p, m1)*\
            (-1)**(js[1]-m1)*wigner_3j(js[1], 1, js[2], -m1, p, m2)*\
            (-1)**(js[2]-m2)*wigner_3j(js[2], 1, js[3], -m2, p, m3)*\
            (-1)**(js[3]-m3)*wigner_3j(js[3], 1, js[0], -m3, p, m0)
    res /= (2*js[0]+1)**0.5

    return float(res)


def brute_rme(js, k):
    """This should just be four-fold, signed Honl-London factor.

    This is mostly for testing."""
    res = (-1)**(js[0]-k)*wigner_3j(js[0], 1, js[1], -k, 0, k)*\
        ((2*js[0]+1)*(2*js[1]+1))**(1/2)*\
        (-1)**(js[1]-k)*wigner_3j(js[1], 1, js[2], -k, 0, k)*\
        ((2*js[1]+1)*(2*js[2]+1))**(1/2)*\
        (-1)**(js[2]-k)*wigner_3j(js[2], 1, js[3], -k, 0, k)*\
        ((2*js[2]+1)*(2*js[3]+1))**(1/2)*\
        (-1)**(js[3]-k)*wigner_3j(js[3], 1, js[0], -k, 0, k)*\
        ((2*js[3]+1)*(2*js[0]+1))**(1/2)
    return res


def wigner6j0(a: int, b: int, c: int, sqrt: Callable=np.sqrt) -> float:
    """Value of Wigner6j(a,b,c;0,c,b)"""
    s = a+b+c
    return (-1)**s/sqrt((2*b+1)*(2*c+1))


def wigner6j1_nn(j1: int, j2: int, j3: int, sqrt: Callable=np.sqrt) -> float:
    """{j1, j2, j3}
       { 1, j3, j2}"""
    s = j1+j2+j3
    nom = 2*(j1*(j1+1)-j2*(j2+1)-j3*(j3+1))
    denom = sqrt( 2*j2*(2*j2+1)*(2*j2+2)*2*j3*(2*j3+1)*(2*j3+2) )
    if denom == 0.0:
        return 0.0
    return (-1)**s*nom/denom


def wigner6j1_nm(j1: int, j2: int, j3: int, sqrt: Callable=np.sqrt) -> float:
    """{j1,   j2, j3}
       { 1, j3-1, j2}"""
    s = j1+j2+j3
    nom = 2*(s+1)*(s-2*j1)*(s-2*j2)*(s-2*j3+1)
    denom = 2*j2*(2*j2+1)*(2*j2+2)*(2*j3-1)*(2*j3)*(2*j3+1)
    if denom == 0.0:
        return 0.0
    return (-1)**s*sqrt(nom/denom)


def wigner6j1_mm(j1: int, j2: int, j3: int, sqrt: Callable=np.sqrt) -> float:
    """{j1,   j2,   j3}
       { 1, j3-1, j2-1}"""
    s = j1+j2+j3
    nom = s*(s+1)*(s-2*j1-1)*(s-2*j1)
    denom = (2*j2-1)*2*j2*(2*j2+1)*(2*j3-1)*2*j3*(2*j3+1)
    if denom == 0.0:
        return 0.0
    return (-1)**s*sqrt(nom/denom)


def wigner6j1_pm(j1: int, j2: int, j3: int, sqrt: Callable=np.sqrt) -> float:
    """{j1,   j2,   j3}
       { 1, j3-1, j2+1}"""
    s = j1+j2+j3
    nom = (s-2*j2-1)*(s-2*j2)*(s-2*j3+1)*(s-2*j3+2)
    denom = (2*j2+1)*(2*j2+2)*(2*j2+3)*(2*j3-1)*2*j3*(2*j3+1)
    if denom == 0.0:
        return 0.0
    return (-1)**s*sqrt(nom/denom)


wigner6j1_map = {
    (0,  0,  0): wigner6j0,
    (1,  0,  0): wigner6j1_nn,
    (1, -1,  0): wigner6j1_nm,
    (1, -1, -1): wigner6j1_mm,
    (1, -1,  1): wigner6j1_pm,
}


def w6j_args_match(args: Sequence[int]) -> Optional[Callable]:
    for pat, func in wigner6j1_map.items():
        if args[2]+pat[1] == args[4] and args[1]+pat[2] == args[5]\
           and args[3] == pat[0]:
            return func


def w6j_equiv_args(args: Sequence[int]) -> List[Tuple[int, ...]]:
    """Return all equivalent lists of Wigner-6j arguments.

    `args` has the form: (j1, j2, j3, j4, j5, j6), which correspond to the
    following Wigner-6j coefficient::

      {j1, j2, j3}
      {j4, j5, j6}
    """
    cols = [(args[0], args[3]), (args[1], args[4]), (args[2], args[5])]
    cols = it.permutations(cols)
    new_cols = []
    for arg_list in cols:
        new_cols.append(arg_list)
        new_cols.append((arg_list[0][::-1], arg_list[1][::-1], arg_list[2]))
        new_cols.append((arg_list[0][::-1], arg_list[1], arg_list[2][::-1]))
        new_cols.append((arg_list[0], arg_list[1][::-1], arg_list[2][::-1]))

    return [(alist[0][0], alist[1][0], alist[2][0], alist[0][1], alist[1][1], alist[2][1])
            for alist in new_cols]


def w6j_special(*args: int, sqrt: Callable=np.sqrt) -> float:
    """Wigner-6j coefficient from analytical expressions."""
    args_list = w6j_equiv_args(args)
    for args_cand in args_list:
        func = w6j_args_match(args_cand)
        if func is not None:
            args_cand = list(args_cand[:3]) + [sqrt]
            return func(*args_cand)


def G(ji: int, jj: int, jk: int, jl: int, k: int, sqrt: Callable=np.sqrt) -> float:
    """Wigner-6j part of four-fold reduced matrix element."""
    # return (2*k+1)*(-1)**(jj+jk+jl-ji)*\
    return (2*k+1)*\
        w6j_special(k, k, 0, ji, ji, jk, sqrt=sqrt)*\
        w6j_special(1, 1, k, jk, ji, jj, sqrt=sqrt)*\
        w6j_special(1, 1, k, ji, jk, jl, sqrt=sqrt)


def T00(phi: float, phj: float, phk: float, phl: float, k: int):
    """Recoupling of four collinear beams with total Q=K=0.

    Only linear polarization.
    """
    if k==0:
        return np.cos(phi-phj)*np.cos(phk-phl)/3.0
    elif k==1:
        return np.sin(phi-phj)*np.sin(phk-phl)*np.sqrt(3)/6
    elif k==2:
        return np.sqrt(5)/60*\
            (np.cos(phi-phj-phk+phl)+\
             np.cos(phi-phj+phk-phl)+\
             6*np.cos(phi+phj-phk-phl))
    else:
        raise ValueError(f'No spherical component for k={k!r}')


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
        raise ValueError(f'No spherical component for k={k!r}')


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
    return G(js[0], js[1], js[2], js[3], 0)*\
        T00_circ(angles[0][0], angles[1][0], angles[2][0], angles[3][0],
                 angles[0][1], angles[1][1], angles[2][1], angles[3][1], 0)+\
        G(js[0], js[1], js[2], js[3], 1)*\
        T00_circ(angles[0][0], angles[1][0], angles[2][0], angles[3][0],
                 angles[0][1], angles[1][1], angles[2][1], angles[3][1], 1)+\
        G(js[0], js[1], js[2], js[3], 2)*\
        T00_circ(angles[0][0], angles[1][0], angles[2][0], angles[3][0],
                 angles[0][1], angles[1][1], angles[2][1], angles[3][1], 2)
