# pylint: disable=line-too-long
"""This module contains functions to derive and manipulate SymPy expressions
related to polarization dependence and angular momentum dependence of four-fold
dipole interaction operator. The derived expressions are in
:mod:`rotsim2d.symbolic.results`."""
# * Imports
from typing import Sequence, Dict, Tuple, Optional
from collections.abc import Mapping, Sequence
import itertools as it
from sympy import *
import sympy.physics.quantum.cg as cg
from molspecutils.molecule import SymTopState, DiatomState
import rotsim2d.dressedleaf as dl
import rotsim2d.pathways as pw
from rotsim2d.symbolic.common import *
from rotsim2d.symbolic.results import rfactors_highj

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
def wigner6j0(j1, j2, j3):
    s = j1+j2+j3
    return (-1)**s/sqrt((2*j2+1)*(2*j3+1))


def wigner6j1_nn(j1, j2, j3):
    """{j1, j2, j3}
       { 1, j3, j2}"""
    s = j1+j2+j3
    nom = 2*(j1*(j1+1)-j2*(j2+1)-j3*(j3+1))
    denom = sqrt( 2*j2*(2*j2+1)*(2*j2+2)*2*j3*(2*j3+1)*(2*j3+2) )
    return (-1)**s*nom/denom


def wigner6j1_nm(j1, j2, j3):
    """{j1,   j2, j3}
       { 1, j3-1, j2}"""
    s = j1+j2+j3
    nom = 2*(s+1)*(s-2*j1)*(s-2*j2)*(s-2*j3+1)
    denom = 2*j2*(2*j2+1)*(2*j2+2)*(2*j3-1)*(2*j3)*(2*j3+1)
    return (-1)**s*sqrt(nom/denom)


# def wigner6j1_mn(j1, j2, j3):
#     """{j1,   j2, j3}
#        { 1, j3, j2-1}"""
#     return wigner6j1_nm(j1, j3, j2)


# def wigner6j1_pn(j1, j2, j3):
#     """{j1,   j2, j3}
#        { 1, j3, j2+1}"""
#     return wigner6j1_mn(j1, j2+1, j3)


# def wigner6j1_np(j1, j2, j3):
#     """{j1,   j2, j3}
#        { 1, j3+1, j2}"""
#     return wigner6j1_pn(j1, j3, j2)


def wigner6j1_mm(j1, j2, j3):
    """{j1,   j2,   j3}
       { 1, j3-1, j2-1}"""
    s = j1+j2+j3
    nom = s*(s+1)*(s-2*j1-1)*(s-2*j1)
    denom = (2*j2-1)*2*j2*(2*j2+1)*(2*j3-1)*2*j3*(2*j3+1)
    return (-1)**s*sqrt(nom/denom)


# def wigner6j1_pp(j1, j2, j3):
#     """{j1,   j2,   j3}
#        { 1, j3+1, j2+1}"""
#     return wigner6j1_mm(j1, j2+1, j3+1)


def wigner6j1_pm(j1, j2, j3):
    """{j1,   j2,   j3}
       { 1, j3-1, j2+1}"""
    s = j1+j2+j3
    nom = (s-2*j2-1)*(s-2*j2)*(s-2*j3+1)*(s-2*j3+2)
    denom = (2*j2+1)*(2*j2+2)*(2*j2+3)*(2*j3-1)*2*j3*(2*j3+1)
    return (-1)**s*sqrt(nom/denom)


# def wigner6j1_mp(j1, j2, j3):
#     """{j1,   j2,   j3}
#        { 1, j3+1, j2-1}"""
#     return wigner6j1_pm(j1, j3, j2)


wigner6j1_map = {
    (0,  0,  0): wigner6j0,
    (1,  0,  0): wigner6j1_nn,
    (1, -1,  0): wigner6j1_nm,
    # (1,  0, -1): wigner6j1_mn,
    # (1,  0,  1): wigner6j1_pn,
    # (1,  1,  0): wigner6j1_np,
    (1, -1, -1): wigner6j1_mm,
    # (1,  1,  1): wigner6j1_pp,
    (1, -1,  1): wigner6j1_pm,
    # (1,  1, -1): wigner6j1_mp
}


def w6j_args_match(args):
    for pat, func in wigner6j1_map.items():
        if args[2]+pat[1] == args[4] and args[1]+pat[2] == args[5]\
           and args[3] == pat[0]:
            return func


def w6j_equiv_args(args):
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
    

def w6j_expr(*args):
    """Return analytical expression for Wigner-6j coefficient."""
    args_list = w6j_equiv_args(args)
    for args_cand in args_list:
        func = w6j_args_match(args_cand)
        if func is not None:
            return func(*args_cand[:3])


def gfactor_expr(ji, jj, jk, jl, k):
    return (2*k+1)*(-1)**(jj+jk+jl-ji)*w6j_expr(k, k, 0, ji, ji, jk)*\
        w6j_expr(1, 1, k, jk, ji, jj)*w6j_expr(1, 1, k, ji, jk, jl)


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

#: Cosines terms present in `T00_exprs` (dummy angles)
T00_trigs = [cos(phi + phj - phk - phl), cos(phi - phj + phk - phl), cos(phi - phj - phk + phl)]
#: Cosines terms present in `T00_exprs` (experimental angles)
T00_theta_trigs = [cos(theta_i + theta_j - theta_k - theta_l),
                   cos(theta_i - theta_j + theta_k - theta_l),
                   cos(theta_i - theta_j - theta_k + theta_l)]

# * R-factors
# Combined polarization-angular momentum factors for four-fold dipole
# interaction operator.
# ** Simplify R-factor
def rfactorize(gterm: Sequence, pterms: Sequence, relative: bool=False, coeffs: bool=False):
    """Produce a nice factorized expression for the R coefficient.

    Parameters
    ----------
    gterm: sequence
        Spherical components of G-factor.
    pterms: sequence
        Spherical components of four-fold polarization tensor.
    relative: bool
        Return polarization-dependent R-factor relative to XXXX polarization.
    coeffs: bool
         Return dict of coefficients instead of a SymPy expression.
    """
    # factor out a term common to all J labels and write polarization in terms
    # of four-angles cosines
    rfactor = sum(gterm[i]*FU['TR8'](pterms[i]) for i in range(3))
    if relative:
        rfactor = rfactor/simplify(rfactor.subs({phi: S(0), phj: S(0), phk: S(0), phl: S(0)}))

    return rfactor_simplify(rfactor, T00_trigs, coeffs=coeffs)


def rfactor_simplify(expr, subterms: Sequence, coeffs: bool=False):
    # Find common multiplicative factor for this specific J label.
    # Substitution necessary to identify (un)common factor by checking free
    # symbols.
    subs1 = dict(zip(subterms, [x1, x2, x3]))
    subs2 = dict(zip([x1, x2, x3], subterms))
    expr = expr.subs(subs1)
    expr_dict = {k: powdenest(factor(powdenest(v, force=True), deep=True), force=True)
                 for k, v in collect(expand(expr), [x1, x2, x3], evaluate=False).items()}
    expr = collect(factor(sum(k*v for k, v in expr_dict.items())), [x1, x2, x3])

    if not coeffs:              # back to cosines
        return expr.subs(subs2)

    common, uncommon = S(1), []
    for term in expr.args:
        if {x1, x2, x3} & term.free_symbols:
            uncommon.append(term)
        else:
            common = common*term
    # special case when common == 1
    # could also be dealt with by checking if top-level operator is Mul or Add
    if len(uncommon) > 1 and common == S(1):
        uncommon = [Add(*uncommon)]

    # `uncommon` is now a sum of >=3 cosines with individual non-factorable coeffs
    # collect terms corresponding to each of three cosines
    uncommon_dict = collect(uncommon[0], [x1, x2, x3], evaluate=False)
    return {'c0': common, 'c12': uncommon_dict[x1], 'c13': uncommon_dict[x2],
            'c14': uncommon_dict[x3]}


def coeff_dict_to_tuple(cdict):
    return (cdict['c0'], cdict['c12'], cdict['c13'], cdict['c14'])


def coeffs_to_expr(coeffs, T00_exprs):
    """Convert dict or tuple of coeffs to SymPy expression."""
    if isinstance(coeffs, Mapping):
        return coeffs['c0']*(coeffs['c12']*T00_exprs[0] +
                             coeffs['c13']*T00_exprs[1] +
                             coeffs['c14']*T00_exprs[2])
    elif len(coeffs) == 4:
        return coeffs[0]*(coeffs[0]*T00_exprs[0] +
                          coeffs[1]*T00_exprs[1] +
                          coeffs[2]*T00_exprs[2])
    elif len(coeffs) == 3:
        return (coeffs[0]*T00_exprs[0] +
                coeffs[1]*T00_exprs[1] +
                coeffs[2]*T00_exprs[2])
    else:
        return ValueError("`coeffs` is not a valid data structure of R-factor coefficients.")


    
# * Classify pathways with regards to R-factor
# TODO Define RFactor class with different representation and high-J variant
# ** Functions
def dl_to_rfactor(dl, rfactors: dict, normalize=False, coeffs=False):
    """Return R-factor corresponding to :class:`rotsim2d.dressedleaf.Pathway`.
    
    Pathway polarization angles need to be indices (integers) from 0 to 3.
    `rfactors` is one of dicts from this module containing SymPy expressions.

    If `normalize` is True, then make the sign of the common factor positive and
    flip the signs of the coefficients if the `c12` coefficient is negative.
    """
    phis = [phi, phj, phk, phl]
    phi_angles = dl._phi_angles(thetas)
    rfac = rfactors[dl.geo_label].subs(dict(zip(phis, phi_angles)))
    # rfac = rfactors[dl.geo_label].subs(
    #     {phi: thetas[dl.angles[0]],
    #      phj: thetas[dl.angles[1]],
    #      phk: thetas[dl.angles[2]],
    #      phl: thetas[dl.angles[3]]})

    if normalize or coeffs:
        rfac_dict = rfactor_simplify(rfac, T00_theta_trigs, coeffs=True)
        if normalize:
            rfac_dict['c0'] = abs(rfac_dict['c0'])
            if rfac_dict['c12'] < 0:
                rfac_dict['c12'] = -rfac_dict['c12']
                rfac_dict['c13'] = -rfac_dict['c13']
                rfac_dict['c14'] = -rfac_dict['c14']

        if not coeffs:
            rfac = coeffs_to_expr(rfac_dict, T00_theta_trigs)
        else:
            rfac = rfac_dict

    return rfac


def dl_gfactors(dl):
    """Return G-factors for a :class:`rotsim2d.dressedleaf.Pathway`."""
    js = list(dl.js)
    return (gfactor_expr(*(js + [0])).evalf(),
            gfactor_expr(*(js + [1])).evalf(),
            gfactor_expr(*(js + [2])).evalf())


def dl_T00s(dl, angles):
    """Return polarization components a :class:`rotsim2d.dressedleaf.Pathway`."""
    phis = [phi, phj, phk, phl]
    phi_angles = dl._phi_angles(angles)
    T00_exprs = [
        cos(phi - phj)*cos(phk - phl)/3,
        sqrt(3)*sin(phi - phj)*sin(phk - phl)/6,
        sqrt(5)*(cos(phi - phj - phk + phl) +
                 cos(phi - phj + phk - phl) +
                 6*cos(phi + phj - phk - phl))/60]

    return (T00_exprs[0].subs(dict(zip(phis, phi_angles))).evalf(),
            T00_exprs[1].subs(dict(zip(phis, phi_angles))).evalf(),
            T00_exprs[2].subs(dict(zip(phis, phi_angles))).evalf())
    

def classify_dls(dressed_pws: Sequence, rfactors: Sequence, states=False) -> dict:
    """Return a mapping between polarization expressions and pathways in `dressed_pws`.

    If `states` is True, convert Pathways to lists of states. `rfactors` are
    R-factors as functions of experimental angles (thetas) not dummy angles
    (phis), i.e. obtained by applying :func:`dl_to_rfactor`.
    """
    if isinstance(rfactors[0], Basic): # SymPy base class
        cmp = lambda x, y: factor(x-y) == S(0)
    else:
        cmp = lambda x, y: x==y

    classified = {}
    for dressed_leaf, cur_rfac in zip(dressed_pws, rfactors):
        found = False
        rfacs = list(classified.keys())
        for rfac in rfacs:
            if cmp(rfac, cur_rfac):
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


def solve_det_angle(rexpr, angles: Optional[Sequence]=None):
    """Return expr. for `tan(theta_l)` which zeroes `rexpr`.

    `angles` contains linear polarization angles of up to three pulses in
    sequence. If `angles` is None, then the first angle is set to 0.
    """
    if angles is None:
        angles = [0]
    expr = expand_trig(rexpr).subs(dict(zip(thetas, angles))).subs(
        {sin(theta_l): x1, cos(theta_l): x2}).subs({x1: x1x2*x2, x2: x1/x1x2})

    return atan(factor(expand_trig(solve(expr, x1x2)[0]), deep=True))


def suppression_angles(exprs: Sequence, angles: Optional[Sequence]=None) -> dict:
    """Return roots of expressions in `exprs` with respect to detection angle.

    `angles` contains linear polarization angles of up to three pulses in
    sequence. If `angles` is None, then the first angle is set to 0.
    """
    return {k: solve_det_angle(k, angles=angles) for k in exprs}


def dummify_angle_expr(angle_expr):
    """Substitute dummy variables for tangents of angles.

    SymPy struggles with solving equations involving several trigonometric
    functions.  Simplify the task by solving for tangents of angles.
    """
    angle_expr = factor(angle_expr.subs({sin(theta_j): x1, cos(theta_j): x2}).subs({x1: x1x2*x2, x2: x1/x1x2}))
    angle_expr = factor(angle_expr.subs({sin(theta_k): x3, cos(theta_k): x4}).subs({x3: x3x4*x4, x4: x3/x3x4}))

    return angle_expr


def common_angles(exprs: Sequence) -> dict:
    """Find angles simultaneously zeroing all in `exprs`.

    Uses :func:`sympy.solve`.
    """
    back_subs = {x0: tan(theta_l), x1x2: tan(theta_j), x3x4: tan(theta_k)}
    system = [dummify_angle_expr(v)-x0 for v in exprs]
    sols = solve(system, [x1x2, x3x4, x0], dict=True)
    sols = [{k.subs(back_subs): v.subs(back_subs) for k, v in s.items()} for s in sols]

    return sols


def classify_suppression(classified: dict, angles: dict) -> dict:
    """Return a map between suppression angles and pathways.

    `classified` is the map between expressions and pathways returned by
    :func:`classify_dls`.
    """
    angles_to_pws = {}
    for k, v in angles.items():
        angles_to_pws.setdefault(v, []).extend(classified[k])

    return angles_to_pws


def pathway_angles(pws: Sequence, angles: Sequence) -> dict:
    """Return a map between detection angles and elements of `pws` they suppress.

    Parameters
    ----------
    pws
        Sequence of class:`dressedleaf.Pathway`.
    angles
        Linear polarization angles of three interacting pulses in sequence.
    """
    pws_rfactors = [dl_to_rfactor(pw, rfactors_highj) for pw in pws]
    classified = classify_dls(pws, pws_rfactors)
    zeroing_angles = suppression_angles(classified.keys(), angles)
    ret = classify_suppression(classified, zeroing_angles)

    return ret


def detection_angles(angles: Sequence, meths=None):
    """Zeroing det. angles for a minimal complete set of pathways.

    Generate a minimal complete set of non-R-factor-equivalent pathways
    including Q-branch transitions. The resulting dict can be used with
    :func:`rotsim2d.dressedleaf.equiv_peaks` or
    :func:`rotsim2d.dressedleaf.split_by_equiv_peaks` to identify other pathways that are
    also zeroed by any of the detection angles.

    Parameters
    ----------
    angles
        Light beam polarizations.
    meths
        Filters for KetBra excitation tree, see
        :func:`rotsim2d.pathways.gen_pathways`.

    Returns
    -------
    det_angles
        Dict from detection angles to symmetric top pathways.
    """
    kbs = pw.gen_pathways([5], rotor='symmetric', kiter_func=lambda x: [1],
                          meths=meths)
    pws = dl.Pathway.from_kb_list(kbs)
    det_angles = pathway_angles(pws, angles)

    return det_angles
