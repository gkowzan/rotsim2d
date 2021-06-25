# pylint: disable=line-too-long
"""This module contains functions to derive and manipulate SymPy expressions
related to polarization dependence and angular momentum dependence of four-fold
dipole interaction operator. The derived expressions are in
:mod:`rotsim2d.symbolic.results`."""
# * Imports
from collections.abc import Mapping
import itertools as it
from sympy import *
import sympy.physics.quantum.cg as cg
from molspecutils.molecule import SymTopState, DiatomState
import rotsim2d.dressedleaf as dl
import rotsim2d.pathways as pw
from rotsim2d.symbolic.common import *
from rotsim2d.symbolic.results import rfactors_highj, gfactors, gfactors_highj, T00_exprs
from typing import Sequence, Dict, List, Tuple, Optional, NewType, Any, Union

#: Dummy type for any SymPy expressions, since SymPy is not annotated
# Expr = NewType('Expr', Any)

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
class RFactor:
    def __init__(self, coeffs: Union[Mapping, Sequence], angles: str='dummy'):
        """Expression describing polarization and angular momentum dependence of a pathway.

        Parameters
        ----------
        coeffs
            Mapping or a sequence of `c0`, `c12`, `c13` and `c14` coefficients.
        angles, optional
            R-factor as a function of 'dummy' (phis) or 'experimental` (thetas)
            angles.
        """
        if isinstance(coeffs, Sequence):
            if len(coeffs) == 3:
                coeffs = (1,) + tuple(coeffs)
            if len(coeffs) == 4:
                self.dict = dict(zip(('c0', 'c12', 'c13', 'c14'),
                                     coeffs))
            else:
                raise ValueError("coeffs has to be a Mapping or a Sequence with len=3 or len=4")
        else:
            self.dict = coeffs

        if angles == 'dummy':
            self.angles = (phi, phj, phk, phl)
            self.trigs = T00_trigs
        elif angles == 'experimental':
            self.trigs = T00_theta_trigs
            self.angles = thetas
        else:
            raise ValueError("angles has to be either 'dummy' or 'experimental'")

        self.dict_to_expr()

    def __repr__(self):
        return f"RFactor(coeffs={self.tuple!r}, angles={self.angles!r})"

    @classmethod
    def from_gterms(cls, gterms: Sequence[Basic], pterms: Sequence[Basic]=T00_exprs,
                    angles: str='dummy'):
        """Make :class:`RFactor` from G-factor and polarization terms.

        Parameters
        ----------
        gterms
            Spherical tensor components of G-factor.
        pterms
            Spherical tensor components of four-fold polarization tensor.
        angles, optional
            R-factor as a function of 'dummy' (phis) or 'experimental` (thetas)
            angles, has to match `pterms`.
        """
        rexpr = sum(gterms[i]*FU['TR8'](pterms[i]) for i in range(3))
        return cls.from_expr(rexpr, angles)

    @classmethod
    def from_expr(cls, rexpr: Basic, angles: str='dummy'):
        """Make :class:`RFactor` from R-factor expression.

        Parameters
        ----------
        rexpr
            SymPy expression for R-factor.
        angles, optional
            R-factor as a function of 'dummy' (phis) or 'experimental` (thetas)
            angles, has to match `rexpr`.
        """
        if angles == 'dummy':
            trigs = T00_trigs
        elif angles == 'experimental':
            trigs = T00_theta_trigs
        _, d = cls._simplify(rexpr, trigs)

        return cls(d, angles=angles)


    @classmethod
    def from_pathway(cls, pw: dl.Pathway, highj: bool=False, normalize: bool=False):
        """Return R-factor corresponding to :class:`rotsim2d.dressedleaf.Pathway`.

        Parameters
        ----------
        pw
            Pathway or DressedPathway.
        highj
            Use high-J limit versions of G-factors.
        normalize
            Make the sign of the common factor positive and flip the signs of
            the coefficients if the `c12` coefficient is negative.
        """
        if highj:
            gterms = gfactors_highj[pw.geo_label]
        else:
            gterms = gfactors[pw.geo_label]

        subs_dict = dict(zip([phi, phj, phk, phl], pw._phi_angles(thetas)))
        pterms = [e.subs(subs_dict) for e in T00_exprs[:]]
        rfac = cls.from_gterms(gterms, pterms, 'experimental')
        if normalize:
            rfac.dict['c0'] = abs(rfac.dict['c0'])
            if rfac.dict['c12'] < 0:
                rfac.dict['c12'] = -rfac.dict['c12']
                rfac.dict['c13'] = -rfac.dict['c13']
                rfac.dict['c14'] = -rfac.dict['c14']
            rfac.dict_to_expr()

        return rfac

    def dict_to_expr(self):
        """Regenerate expression from dict."""
        self.expr = self.dict['c0']*(self.dict['c12']*self.trigs[0] +
                                     self.dict['c13']*self.trigs[1] +
                                     self.dict['c14']*self.trigs[2])

    @property
    def tuple(self):
        """Return coefficients as a tuple."""
        return (self.dict['c0'], self.dict['c12'], self.dict['c13'], self.dict['c14'])

    @staticmethod
    def _simplify(expr: Basic, subterms: Sequence[Basic]) -> Tuple[Basic, Dict]:
        """Simplify RFactor SymPy expression.

        Parameters
        ----------
        expr
            Expression to simplify,
        subterms
            Collect coefficients for expressions in `subterms`

        Returns
        -------
        Tuple factorized expression for R-factor and dict of coefficients.
        """
        subs1 = dict(zip(subterms, [x1, x2, x3]))
        subs2 = dict(zip([x1, x2, x3], subterms))
        expr = expr.subs(subs1)
        expr_dict = {k: powdenest(factor(powdenest(v, force=True), deep=True), force=True)
                     for k, v in collect(expand(expr), [x1, x2, x3], evaluate=False).items()}
        expr = collect(factor(sum(k*v for k, v in expr_dict.items())), [x1, x2, x3])
        ret_expr  = expr.subs(subs2)

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
        ret_dict = {'c0': common, 'c12': uncommon_dict[x1], 'c13': uncommon_dict[x2],
                    'c14': uncommon_dict[x3]}

        return ret_expr, ret_dict

    def expr_xxxx(self):
        """Return R-factor expression for XXXX polarization."""
        expr = self.expr.subs(dict(zip(self.angles, [0]*4)))

        return simplify(factor(powdenest(expr, force=True), deep=True))

    def expr_relative(self):
        """Return R-factor expressions relative to XXXX polarization."""
        return self._simplify(self.expr/self.expr_xxxx(), self.trigs)[0]

    def __eq__(self, o):
        if isinstance(o, RFactor):
            return self == o.tuple
        elif isinstance(o, Mapping):
            return all(factor(self.dict[k]-o[k], deep=True) == S(0) for k in self.dict)
        else:
            return all(factor(st-ot, deep=True) == S(0) for st, ot in zip(self.tuple, o))

    def __hash__(self):
        # could also convert expressions to str but I'm not sure if the strings
        # will be unambiguous
        return hash(tuple([float(x.subs(J_i, 100) if isinstance(x, Basic) else x)
                           for x in self.tuple]))


# * Classify pathways with regards to R-factor
# ** Functions
def dl_gfactors(dl: dl.Pathway):
    """Return G-factors for a :class:`rotsim2d.dressedleaf.Pathway`."""
    js = list(dl.js)
    return (gfactor_expr(*(js + [0])).evalf(),
            gfactor_expr(*(js + [1])).evalf(),
            gfactor_expr(*(js + [2])).evalf())


def dl_T00s(dl: dl.Pathway, angles):
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
    

def classify_dls(dressed_pws: Sequence[dl.Pathway], rfactors: Sequence[RFactor],
                 states: bool=False) -> Dict[RFactor, dl.Pathway]:
    """Return a mapping between polarization expressions and pathways in `dressed_pws`.

    If `states` is True, convert Pathways to lists of states. `rfactors` are
    R-factors as functions of experimental angles (thetas) not dummy angles
    (phis).
    """
    classified = {}
    for dressed_leaf, cur_rfac in zip(dressed_pws, rfactors):
        found = False
        rfacs = list(classified.keys())
        for rfac in rfacs:
            if rfac == cur_rfac:
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


def solve_det_angle(rexpr: RFactor, angles: Optional[Sequence]=None) -> Basic:
    """Return expr. for `tan(theta_l)` which zeroes `rexpr`.

    `angles` contains linear polarization angles of up to three pulses in
    sequence. If `angles` is None, then the first angle is set to 0.
    """
    if angles is None:
        angles = [0]
    expr = expand_trig(rexpr.expr).subs(dict(zip(thetas, angles))).subs(
        {sin(theta_l): x1, cos(theta_l): x2}).subs({x1: x1x2*x2, x2: x1/x1x2})

    return atan(factor(expand_trig(solve(expr, x1x2)[0]), deep=True))


def suppression_angles(exprs: Sequence[RFactor],
                       angles: Optional[Sequence[Basic]]=None) -> Dict[Basic, Basic]:
    """Return roots of expressions in `exprs` with respect to detection angle.

    `angles` contains linear polarization angles of up to three pulses in
    sequence. If `angles` is None, then the first angle is set to 0.
    """
    return {k: solve_det_angle(k, angles=angles) for k in exprs}


def dummify_angle_expr(angle_expr: Basic) -> Basic:
    """Substitute dummy variables for tangents of angles.

    SymPy struggles with solving equations involving several trigonometric
    functions.  Simplify the task by solving for tangents of angles.
    """
    angle_expr = factor(angle_expr.subs({sin(theta_j): x1, cos(theta_j): x2}).subs({x1: x1x2*x2, x2: x1/x1x2}))
    angle_expr = factor(angle_expr.subs({sin(theta_k): x3, cos(theta_k): x4}).subs({x3: x3x4*x4, x4: x3/x3x4}))

    return angle_expr


def common_angles(exprs: Sequence[Basic]) -> Dict[Basic, Basic]:
    """Find angles simultaneously zeroing all in `exprs`.

    Uses :func:`sympy.solve`.
    """
    back_subs = {x0: tan(theta_l), x1x2: tan(theta_j), x3x4: tan(theta_k)}
    system = [dummify_angle_expr(v)-x0 for v in exprs]
    sols = solve(system, [x1x2, x3x4, x0], dict=True)
    sols = [{k.subs(back_subs): v.subs(back_subs) for k, v in s.items()} for s in sols]

    return sols


def classify_suppression(classified: Dict[RFactor, dl.Pathway],
                         angles: Dict[RFactor, Basic]) -> Dict[Basic, dl.Pathway]:
    """Return a map between suppression angles and pathways.

    `classified` is the map between expressions and pathways returned by
    :func:`classify_dls`.
    """
    angles_to_pws = {}
    for k, v in angles.items():
        angles_to_pws.setdefault(v, []).extend(classified[k])

    return angles_to_pws


def pathway_angles(pws: Sequence[dl.Pathway], angles: Sequence) -> Dict[Basic, List[dl.Pathway]]:
    """Return a map between detection angles and elements of `pws` they suppress.

    Parameters
    ----------
    pws
        Sequence of class:`dressedleaf.Pathway`.
    angles
        Linear polarization angles of three interacting pulses in sequence.
    """
    pws_rfactors = [RFactor.from_pathway(pw, True, True) for pw in pws]
    # pws_rfactors = [dl_to_rfactor(pw, rfactors_highj) for pw in pws]
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
