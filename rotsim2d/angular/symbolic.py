# pylint: disable=line-too-long
"""This module contains functions to derive and manipulate SymPy expressions
related to polarization dependence and angular momentum dependence of four-fold
dipole interaction operator. The derived expressions are in
:mod:`rotsim2d.angular.symbolic_results`."""
# * Imports
from typing import Sequence, Dict, Tuple, Optional 
from sympy import *
import sympy.physics.quantum.cg as cg
from spectroscopy.molecule import SymTopState, DiatomState

J_i = symbols("J_i", integer=True, nonnegative=True)
phi, phj, phk, phl = symbols(r"\phi_i \phi_j \phi_k \phi_l", real=True)
thetas = symbols(r"\theta_i \theta_j \theta_k \theta_l", real=True)
theta_i, theta_j, theta_k, theta_l = thetas
x0, x1, x2, x3 = symbols("x0 x1 x2 x3", real=True)

# ** Utilities
def inf_syms():
    i = 0
    while True:
        yield Symbol("x"+str(i))
        i += 1

# * Analytical Wigner-6j functions
# Analytical formulas for Wigner-6j coefficients for arguments differing by
# small factors. Based on formulas from "Microwave molecular spectra" by Gordy
# and Cook.
def wigner6j1_nn(j1, j2, j3):
    """{j1, j2, j3}
       { 1, j3, j2}"""
    s = j1+j2+j3
    nom = j1*(j1+1)-j2*(j2+1)-j3*(j3+1)
    denom = sqrt( 2*j2*(2*j2+1)*(2*j2+2)*2*j3*(2*j3+1)*(2*j3+2) )
    return (-1)**s*nom/denom


def wigner6j1_nm(j1, j2, j3):
    """{j1,   j2, j3}
       { 1, j3-1, j2}"""
    s = j1+j2+j3
    nom = 2*(s+1)*(s-2*j1)*(s-2*j2)*(s-2*j3+1)
    denom = 2*j2*(2*j2+1)*(2*j2+2)*(2*j3-1)*(2*j3)*(2*j3+1)
    return (-1)**s*sqrt(nom/denom)


def wigner6j1_mn(j1, j2, j3):
    """{j1,   j2, j3}
       { 1, j3, j2-1}"""
    return wigner6j1_nm(j1, j3, j2)


def wigner6j1_pn(j1, j2, j3):
    """{j1,   j2, j3}
       { 1, j3, j2+1}"""
    return wigner6j1_mn(j1, j2+1, j3)


def wigner6j1_np(j1, j2, j3):
    """{j1,   j2, j3}
       { 1, j3+, j2}"""
    return wigner6j1_pn(j1, j3, j2)


def wigner6j1_mm(j1, j2, j3):
    """{j1,   j2,   j3}
       { 1, j3-1, j2-1}"""
    s = j1+j2+j3
    nom = s*(s+1)*(s-2*j1-1)*(s-2*j1)
    denom = (2*j2-1)*2*j2*(2*j2+1)*(2*j3-1)*2*j3*(2*j3+1)
    return (-1)**s*sqrt(nom/denom)


def wigner6j1_pp(j1, j2, j3):
    """{j1,   j2,   j3}
       { 1, j3-1, j2-1}"""
    return wigner6j1_mm(j1, j2+1, j3+1)


def wigner6j1_pm(j1, j2, j3):
    """{j1,   j2,   j3}
       { 1, j3-1, j2+1}"""
    s = j1+j2+j3
    nom = (s-2*j2-1)*(s-2*j2)*(s-2*j3+1)*(s-2*j3+2)
    denom = (2*j2+1)*(2*j2+2)*(2*j2+3)*(2*j3-1)*2*j3*(2*j3+1)
    return (-1)**s*sqrt(nom/denom)


def wigner6j1_mp(j1, j2, j3):
    """{j1,   j2,   j3}
       { 1, j3+1, j2-1}"""
    return wigner6j1_pm(j1, j3, j2)

# * Polarization
# ** Linear polarization
def T1q(q, ph):
    """Linear polarization spherical tensor component `q`."""
    if q == 1:
        return -1/sqrt(2)*exp(1j*ph)
    elif q == -1:
        return 1/sqrt(2)*exp(-1j*ph)
    else:
        return 0


def T00_phis(k, phis):
    """k-th component of isotropic four-fold polarization tensor."""
    pre = 2*k+1
    ret = 0
    for q in range(-k, k+1):
        for qp in (-1, 0, 1):
            for qpp in (-1, 0 ,1):
                ret += cg.wigner_3j(k, k, 0, q, -q, 0)*\
                    cg.wigner_3j(1,1,k,qp,q-qp,-q)*\
                    cg.wigner_3j(1,1,k,qpp,-q-qpp,q)*\
                    T1q(qp,phis[0])*T1q(q-qp,phis[1])*T1q(qpp,phis[2])*T1q(-q-qpp,phis[3])
    return pre*ret

#: Cosines terms present in :var:`T00_exprs`
T00_trigs = [cos(phi + phj - phk - phl), cos(phi - phj + phk - phl), cos(phi - phj - phk + phl)]

# * R-factors
# Combined polarization-angular momentum factors for four-fold dipole
# interaction operator.
# ** Simplify R-factor
def rfactorize(gterm: Sequence, pterms: Sequence, cfac: bool=True, relative: bool=False, coeffs: bool=False):
    """Produce a nice factorized expression for the R coefficient.

    Parameters
    ----------
    gterm: sequence
        Spherical components of G-factor.
    pterms: sequence
        Spherical components of four-fold polarization tensor.
    cfac: bool
        Separate the 1/(60*(2*J_i+1)**1.5) factor.
    relative: bool
        Return polarization-dependent R-factor relative to XXXX polarization.
    coeffs: bool
         Return dict of coefficients instead of a SymPy expression.
    """
    # factor out a term common to all J labels and write polarization in terms
    # of four-angles cosines
    if cfac:
        cfac = S(1)/(60*(2*J_i+1)**(S(3)/S(2)))
    else:
        cfac = S(1)
    rfactor = sum(gterm[i]*FU['TR8'](pterms[i]) for i in range(3))/cfac
    if relative:
        rfactor = rfactor/simplify(rfactor.subs({phi: S(0), phj: S(0), phk: S(0), phl: S(0)}))
    
    # find common multiplicative factor for this specific J label
    # substitution necessary to identify (un)common factor by checking free
    # symbols
    subs1 = dict(zip(T00_trigs, [x1, x2, x3]))
    subs2 = dict(zip([x1, x2, x3], T00_trigs))
    rfactor = rfactor.subs(subs1)
    rfactor_dict = {k: powdenest(factor(powdenest(v, force=True), deep=True), force=True)
                    for k, v in collect(expand(rfactor), [x1, x2, x3], evaluate=False).items()}
    rfactor = collect(factor(sum(k*v for k, v in rfactor_dict.items())), [x1, x2, x3])

    if not coeffs:
        return rfactor.subs(subs2)*cfac

    common, uncommon = S(1), []
    for term in rfactor.args:
        if {x1, x2, x3} & term.free_symbols:
            uncommon.append(term)
        else:
            common = common*term
    # print(common, uncommon)
    # special case when common == 1
    # could also be dealt with by checking if top-level operator is Mul or Add
    if len(uncommon) > 1 and common == S(1):
        uncommon = [Add(*uncommon)]
    common = cfac*common

    # `uncommon` is now a sum of >=3 cosines with individual non-factorable coeffs
    # collect terms corresponding to each of three cosines
    uncommon_dict = collect(uncommon[0], [x1, x2, x3], evaluate=False)
    return {'c0': common, 'c12': uncommon_dict[x1], 'c13': uncommon_dict[x2],
            'c14': uncommon_dict[x3]}

# * Classify pathways with regards to R-factor
# ** Functions
def dl_to_rfactor(dl, rfactors):
    """Return R-factor corresponding to DressedLeaf.
    
    DressedLeaf polarization angles need to be indices (integers) from 0 to 3.
    `rfactors` is one of dicts from this modules containing SymPy expressions.
    """
    return rfactors[dl.geo_label].subs(
        {phi: thetas[dl.angles[0]],
         phj: thetas[dl.angles[1]],
         phk: thetas[dl.angles[2]],
         phl: thetas[dl.angles[3]]})


def classify_dls(dressed_pws: Sequence, rfactors, states=False) -> dict:
    """Return a mapping between polarization expressions and pathways in dressed_pws.

    If `states` is True, convert Pathways to lists of states."""
    classified = {}
    for dressed_leaf, cur_rfac in zip(dressed_pws, rfactors):
        found = False
        rfacs = list(classified.keys())
        for rfac in rfacs:
            if factor(rfac-cur_rfac) == S(0):
                found = rfac
                break
        if found:
            classified[found].append(dressed_leaf)
        else:
            classified[cur_rfac] = [dressed_leaf]

    if states:
        return {k: [pw.leaf.to_statelist(diatom=True, normalize=True) for pw in v]
                for k, v in classified.items()}
            
    return classified


def suppression_angles(exprs: Sequence, angles: Sequence) -> dict:
    """Return roots of expressions in `exprs` with respect to detection angle.

    `angles` contains linear polarization angles of three pulses in sequence.
    """
    return {k: solve(expand_trig(k.subs({theta_i: angles[0],
                                         theta_j: angles[1],
                                         theta_k: angles[2]})), theta_l)[0]
            for k in exprs}


def classify_suppression(classified: dict, angles: dict):
    """Return a map between suppression angles and pathways."""
    angles_to_pws = {}
    for k, v in angles.items():
        angles_to_pws.setdefault(v, []).extend(classified[k])

    return angles_to_pws
