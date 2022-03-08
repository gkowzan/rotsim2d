# pylint: disable=line-too-long
"""This module contains functions to derive and manipulate SymPy expressions
related to polarization dependence and angular momentum dependence of four-fold
dipole interaction operator. The derived expressions are in
:mod:`rotsim2d.symbolic.results`."""
# * Imports
import inspect
import itertools as it
import re
from collections.abc import Mapping
from functools import reduce
from operator import mul

import numpy as np
import rotsim2d.dressedleaf as dl
import rotsim2d.pathways as pw
import sympy.physics.quantum.cg as cg
from molspecutils.molecule import RotState
from rotsim2d.couple import T00, G
from rotsim2d.symbolic.common import *
from rotsim2d.symbolic.results import (T00_exprs, gfactors, gfactors_highj,
                                       gfactors_highj_numeric)
from scipy.optimize import minimize
from sympy import *
from typing import (Any, Callable, Dict, Iterable, List, NewType, Optional,
                    Sequence, Tuple, Union)

#: Dummy type for any SymPy expressions, since SymPy is not annotated
# Expr = NewType('Expr', Any)

# * Utilities
def inf_syms():
    i = 0
    while True:
        yield Symbol("x"+str(i))
        i += 1

# * Analytical Wigner-6j functions
def gfactor_expr(ji, jj, jk, jl, k):
    """Call :func:`rotsim2d.couple.G` with symbolic `sqrt`."""
    return G(ji, jj, jk, jl, k, sqrt=sqrt)


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
        """Expression describing polarization and angular momentum dependence of a
        pathway.

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
            self.trigs = T00_trigs
            self.angles = (phi, phj, phk, phl)
            self.angles_type = angles
        elif angles == 'experimental':
            self.trigs = T00_theta_trigs
            self.angles = thetas
            self.angles_type = angles
        else:
            raise ValueError("angles has to be either 'dummy' or 'experimental'")

        self._numeric = None
        self._numeric_rel = None
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
        else:
            raise ValueError("`trigs` has to be either 'experimental' or 'dummy'")
        _, d = cls._simplify(rexpr, trigs) # 2nd run gets rid of last two dups
        _, d = cls._simplify(_, trigs)

        return cls(d, angles=angles)


    @classmethod
    def from_pathway(cls, pw: dl.Pathway, highj: bool=False, normalize: bool=False,
                     only_angles: bool=False, rotated: Optional[str]=None):
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
        only_angles
            Set 'c00' to 1.
        """
        if highj:
            gterms = gfactors_highj[pw.geo_label]
        else:
            gterms = gfactors[pw.geo_label]
        if rotated is not None:
            if rotated == 'forward':
                deltaj = pw.js[0]-pw.js[-1] # Jj-Ji, Jj = Ji + (Jj-Ji)
                gterms = tuple([(S(gterm).subs(J_i, J_j)*sqrt((2*J_j+1)/(2*J_i+1)))\
                                .subs(J_j, J_i+deltaj) for gterm in gterms])
            elif rotated == 'backward':
                deltaj = pw.js[0]-pw.js[1] # Jl-Ji, Jl = Ji + (Jl-Ji)
                gterms = tuple([(S(gterm).subs(J_i, J_l)*sqrt((2*J_l+1)/(2*J_i+1)))\
                                .subs(J_l, J_i+deltaj) for gterm in gterms])
            else:
                raise ValueError("`rotated` needs to be either 'forward' or 'backward'")

        subs_dict = dict(zip([phi, phj, phk, phl], pw._phi_angles(thetas)))
        pterms = [e.subs(subs_dict) for e in T00_exprs[:]]
        rfac = cls.from_gterms(gterms, pterms, 'experimental')
        if only_angles:
            rfac.dict['c0'] = 1
            rfac.dict_to_expr()
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
        expr_dict = {k: powdenest(factor(powdenest(v, force=True), deep=True),
                                  force=True)
                     for k, v in collect(expand(expr), [x1, x2, x3],
                                         evaluate=False).items()}
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

        # `uncommon` is now a sum of >=3 cosines with individual non-factorable
        # coeffs
        # collect terms corresponding to each of three cosines
        uncommon_dict = collect(uncommon[0], [x1, x2, x3], evaluate=False)
        ret_dict = {'c0': common,
                    'c12': uncommon_dict[x1],
                    'c13': uncommon_dict[x2],
                    'c14': uncommon_dict[x3]}

        return ret_expr, ret_dict

    def expr_xxxx(self):
        """Return R-factor expression for XXXX polarization."""
        expr = self.expr.subs(dict(zip(self.angles, [0]*4)))

        return simplify(factor(powdenest(expr, force=True), deep=True))

    def expr_relative(self):
        """Return R-factor expressions relative to XXXX polarization."""
        return self._simplify(self.expr/self.expr_xxxx(), self.trigs)[0]

    def numeric(self, thetai, thetaj, thetak, thetal, Ji=None):
        """Try to efficiently evaluate R-factor with NumPy.

        It is always safe to pass `Ji` argument. It will be ignored if it's not
        needed. This is because even fully general R-factors sometimes do not
        depend on `J_i`.
        """
        if self._numeric is None:
            args = list(self.angles)
            if J_i in self.expr.free_symbols:
                args = args + [J_i]
            self._numeric = lambdify(args, self.expr)
            self._numeric_nargs = len(inspect.signature(
                self._numeric).parameters)

        args = (thetai, thetaj, thetak, thetal, Ji)

        return self._numeric(*args[:self._numeric_nargs])

    def numeric_rel(self, thetai, thetaj, thetak, thetal, Ji=None):
        """Try to efficiently evaluate relative R-factor with NumPy.

        It is always safe to pass `Ji` argument. It will be ignored if it's not
        needed. This is because even fully general R-factors sometimes do not
        depend on `J_i`.
        """
        if self._numeric_rel is None:
            args = list(self.angles)
            expr = self.expr_relative()
            if J_i in expr.free_symbols:
                args = args + [J_i]
            self._numeric_rel = lambdify(args, expr)
            self._numeric_rel_nargs = len(inspect.signature(
                self._numeric_rel).parameters)

        args = (thetai, thetaj, thetak, thetal, Ji)

        return self._numeric_rel(*args[:self._numeric_rel_nargs])

    def det_angle(self, angles: Optional[Sequence]=None) -> Basic:
        """Return expr. for `tan(theta_l)` which zeroes `rexpr`.

        `angles` contains linear polarization angles of up to three pulses in
        sequence. If `angles` is None, then the first angle is set to 0.
        """
        return solve_det_angle(self, angles)
    
    def __eq__(self, o):
        if isinstance(o, RFactor):
            return self == o.tuple
        elif isinstance(o, Mapping):
            return all(factor(self.dict[k]-o[k], deep=True) == S(0)
                       for k in self.dict)
        else:
            return all(factor(st-ot, deep=True) == S(0)
                       for st, ot in zip(self.tuple, o))

    def __hash__(self):
        # could also convert expressions to str but I'm not sure if the strings
        # will be unambiguous
        return hash(tuple([float(x.subs(J_i, 100) if isinstance(x, Basic) else x)
                           for x in self.tuple]))


# * Classify pathways with regards to R-factor
# ** Functions
def dl_gfactors(dl: dl.Pathway, evalf: bool=True, abstract: bool=False):
    """Return G-factors for a :class:`rotsim2d.dressedleaf.Pathway`."""
    js = list(dl.js)
    if abstract:
        evalf = False
        js = [J_i+(j-js[0]) for j in js]
    gfactors = (gfactor_expr(*(js + [0])),
                gfactor_expr(*(js + [1])),
                gfactor_expr(*(js + [2])))
    if evalf:
        gfactors = tuple(gfactor.evalf() for gfactor in gfactors)

    return gfactors


def dl_T00s(dl: dl.Pathway, angles, evalf: bool=True):
    """Return polarization components a :class:`rotsim2d.dressedleaf.Pathway`."""
    phis = [phi, phj, phk, phl]
    phi_angles = dl._phi_angles(angles)
    T00_exprs = [
        cos(phi - phj)*cos(phk - phl)/3,
        sqrt(3)*sin(phi - phj)*sin(phk - phl)/6,
        sqrt(5)*(cos(phi - phj - phk + phl) +
                 cos(phi - phj + phk - phl) +
                 6*cos(phi + phj - phk - phl))/60]
    T00_exprs = tuple(T00_expr.subs(dict(zip(phis, phi_angles)))
                      for T00_expr in T00_exprs)

    if evalf:
        T00_exprs = tuple(T00_expr.evalf() for T00_expr in T00_exprs)

    return T00_exprs
    

def classify_dls(dressed_pws: Sequence[dl.Pathway], rfactors: Sequence[RFactor],
                 states: bool=False) -> Dict[RFactor, List[dl.Pathway]]:
    """Return a mapping between R-factors and pathways in `dressed_pws`.

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
    expr = expand_trig(rexpr.expr).subs(dict(zip(rexpr.angles, angles))).subs(
        {sin(rexpr.angles[3]): x1, cos(rexpr.angles[3]): x2}).subs({x1: x1x2*x2, x2: x1/x1x2})

    sol = solve(expr, x1x2)
    if sol:
        return atan(factor(expand_trig(sol[0]), deep=True))

    sol = solve(expr, x1)
    if sol:
        return asin(factor(expand_trig(sol[0]), deep=True))

    sol = solve(expr, x2)
    if sol:
        return acos(factor(expand_trig(sol[0]), deep=True))


def suppression_angles(exprs: Sequence[RFactor],
                       angles: Optional[Sequence[Basic]]=None) -> Dict[Basic, Basic]:
    """Return roots of expressions in `exprs` with respect to detection angle.

    `angles` contains linear polarization angles of up to three pulses in
    sequence. If `angles` is None, then the first angle is set to 0.
    """
    return {k: solve_det_angle(k, angles=angles)
            for k in exprs}


def dummify_angle_expr(angle_expr: Basic,
                       angles_vars: Optional[Sequence[Basic]]=thetas) -> Basic:
    """Substitute dummy variables for tangents of angles.

    SymPy struggles with solving equations involving several trigonometric
    functions.  Simplify the task by solving for tangents of angles.
    """
    angle_expr = factor(
        angle_expr.\
        subs({sin(angles_vars[1]): x1, cos(angles_vars[1]): x2}).\
        subs({x1: x1x2*x2, x2: x1/x1x2}))
    angle_expr = factor(
        angle_expr.\
        subs({sin(angles_vars[2]): x3, cos(angles_vars[2]): x4}).\
        subs({x3: x3x4*x4, x4: x3/x3x4}))

    return angle_expr


def common_angles(exprs: Sequence[Basic],
                  angles_vars: Optional[Sequence[Basic]]=thetas) -> Dict[Basic, Basic]:
    r"""Find angles simultaneously zeroing all in `exprs`.

    The expressions in `exprs` should be formulas for the tan of detection
    angle, :math:`\theta_l`, zeroing some pathways, as obtained by calling
    :func:`solve_det_angle`. Uses :func:`sympy.solve`.
    """
    back_subs = {x0: tan(angles_vars[3]), x1x2: tan(angles_vars[1]), x3x4: tan(angles_vars[2])}
    system = [dummify_angle_expr(v, angles_vars)-x0 for v in exprs]
    sols = solve(system, [x1x2, x3x4, x0], dict=True)
    sols = [{k.subs(back_subs): v.subs(back_subs) for k, v in s.items()} for s in sols]

    return sols


def classify_suppression(classified: Dict[RFactor, Sequence[dl.Pathway]],
                         angles: Dict[RFactor, Basic]) -> Dict[Basic, List[dl.Pathway]]:
    """Return a map between suppression angles and pathways.

    `classified` is the map between expressions and pathways returned by
    :func:`classify_dls`.
    """
    angles_to_pws = {}
    for k, v in angles.items():
        angles_to_pws.setdefault(v, []).extend(classified[k])

    return angles_to_pws


def pathway_angles(pws: Sequence[dl.Pathway], angles: Sequence,
                   angles_vars: Optional[Sequence[Basic]]=thetas) -> Dict[Basic, List[dl.Pathway]]:
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
    zeroing_angles = suppression_angles(classified.keys(), angles,
                                        angles_vars)
    ret = classify_suppression(classified, zeroing_angles)

    return ret


def detection_angles(angles: Sequence, meths: Optional[Sequence[Callable]]=None,
                     angles_vars: Optional[Sequence[Basic]]=thetas):
    """Zeroing det. angles for a minimal complete set of pathways.

    Generate a minimal complete set of non-R-factor-equivalent pathways
    including Q-branch transitions. The resulting dict can be used with
    :func:`rotsim2d.dressedleaf.equiv_peaks` or
    :func:`rotsim2d.dressedleaf.split_by_equiv_peaks` to identify other pathways
    that are also zeroed by any of the detection angles.

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
    det_angles = pathway_angles(pws, angles, angles_vars)

    return det_angles


# ** RFactorPathways
class RFactorPathways:
    """Container for :class:`RFactor` and associated
    :class:`rotsim2d.dressedleaf.Pathway`s.

    Attributes
    ----------
    rfactor: RFactor
        R-factor associated with these pathways.
    pws: list of :class:`rotsim2d.dressedleaf.Pathway`
        List of pathways.
    trans_labels: list of str
        Combined transition labels for all included pathways.
    trans_labels_deg: list of str
        Combined degenerate transition labels for all included pathways.
    peak_labels: list of str
        Combined 2D peak labels for all included pathways.
    """
    def __init__(self, rfactor: RFactor, pws: List[dl.Pathway]):
        self.rfactor = rfactor
        self.pws = pws
        self.props = {}
        self.det_angle = None

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "RFactorPathways(rfactor={rfac!r}, peak_labels={pl!r}, trans_labels_deg={tld!r})".format(
            rfac=self.rfactor.tuple, tld=self.trans_labels_deg, pl=self.peak_labels
        )
    
    @classmethod
    def from_pwlist(cls, pwlist: Sequence[dl.Pathway], highj: bool=False,
                    normalize: bool=False) -> List['RFactorPathways']:
        rfactors = [RFactor.from_pathway(pw, highj, normalize) for pw in pwlist]
        classified_pws = classify_dls(pwlist, rfactors)

        return [cls(rfactor, pws) for rfactor, pws in classified_pws.items()]

    def calc_det_angle(self, angles: Optional[Sequence]=None):
        """Calculate zeroing detection angle for this R-factor."""
        self.det_angle = solve_det_angle(self.rfactor, angles)

    @property
    def trans_labels(self) -> List[str]:
        ":meta private:"
        labels = [pw.trans_label for pw in self.pws]
        labels.sort(key=lambda x: re.sub(r'[)(0-9)]', '', x))

        return labels

    @property
    def trans_labels_deg(self) -> List[str]:
        ":meta private:"
        labels = list(set([pw.trans_label_deg for pw in self.pws]))
        labels.sort()

        return labels

    @property
    def peak_labels(self) -> List[str]:
        ":meta private:"
        labels = list(set([pw.peak_label for pw in self.pws]))
        labels.sort()

        return labels

    def add_property(self, name: str, d: Mapping):
        assert len(d) == len(self.pws)
        self.props[name] = d


def calc_rfactors(rf_pw: RFactorPathways, *args, **kwargs):
    """Calculate J-dependent R-factors and uniquify them."""
    rf_pw.add_property(
        'rfactors', [RFactor.from_pathway(pw, *args, **kwargs)
                     for pw in rf_pw.pws])


def calc_angle_funcs(rf_pw: RFactorPathways, angles: Sequence):
    """Generate J-dependent zeroing angle functions."""
    angle_funcs, angle_exprs = [], []
    for rf in rf_pw.props['rfactors']:
        rf_det_angle = tan(solve_det_angle(rf, angles=angles[:3]))
        angle_exprs.append(rf_det_angle)
        if J_i in rf_det_angle.free_symbols:
            angle_funcs.append(
                lambdify(J_i, simplify(rf_det_angle)))
        else:
                angle_funcs.append(None)

    rf_pw.add_property('angle_funcs', angle_funcs)
    rf_pw.add_property('angle_exprs', angle_exprs)


def calc_angle_amps(rf_pw: RFactorPathways, angles_num: Sequence[float],
                    js: Sequence[int]):
    """Calculate geometric amplitudes for specific angles and J_i values."""
    amps = []
    for rf in rf_pw.props['rfactors']:
        d = np.atleast_1d(rf.numeric_rel(*(angles_num+[js])))
        if len(d) == 1:
            d = np.full(len(js), d[0])
        amps.append(d)

    rf_pw.add_property('angle_amps', amps)


def calc_relative_exprs(rf_pw: RFactorPathways, angles: Sequence):
    """Calculate relative amplitudes for specified angles."""
    exprs = []
    for rf in rf_pw.props['rfactors']:
        expr = rf.expr_relative()
        expr = expr.subs(dict(zip(rf.angles, angles)))
        expr = fu(expr)     # sym.FU['TR11'] might be enough

        exprs.append(expr)

    rf_pw.add_property('exprs_relative', exprs)


def calc_amplitude_exprs(rf_pw: RFactorPathways, angles: Sequence,
                         norm: bool=True):
    """Calculate pathway amplitude expressions."""
    exprs = []
    for pw1 in rf_pw.pws:
        rfactor = RFactor.from_pathway(pw1)
        exprs.append(
            (amplitude_expr_pw(pw1, norm=norm)*\
             rfactor.expr.subs(
                 dict(zip(rfactor.angles, angles)))).\
            replace(Abs, Id))

    rf_pw.add_property('amplitude_exprs', exprs)

    
def optimize_contrast(rfpws_min: Sequence[RFactorPathways],
                      rfpws_max: Sequence[RFactorPathways],
                      initial_guess: Sequence[float],
                      **min_kwargs):
    """Optimize contrast between two sets of R-factors.

    This is useful in case there are no analytical angles that simultaneously
    zero all R-factors in `rfpws_min`. Keyword arguments are passed to
    scipy.optimize.minimize. Returns OptimizeResult object.
    """
    min_pws = sum((rfpw.pws for rfpw in rfpws_min), [])
    min_pws = list(dl.split_by_peaks(min_pws).values())
    max_pws = sum((rfpw.pws for rfpw in rfpws_max), [])
    max_pws = list(dl.split_by_peaks(max_pws).values())

    def sum_rfactors(peaks: Iterable[dl.Pathway], angles: Sequence[float]):
        angles = list(angles)
        sum_geo = 0.0
        for peak_pws in peaks:
            peak_sum = 0.0
            for pw in peak_pws:
                thetas = pw._phi_angles([0.0] + angles)
                gfacs = gfactors_highj_numeric[pw.geo_label]
                for k in (0, 1, 2):
                    peak_sum += gfacs[k]*T00(*(thetas + [k]))
            sum_geo += abs(peak_sum)

        return sum_geo

    def opt_func(thetas):
        min_geo = sum_rfactors(min_pws, thetas)
        max_geo = sum_rfactors(max_pws, thetas)
        print("max_geo={:f}, min_geo={:f}".format(max_geo, min_geo))

        return -(max_geo-min_geo)

    return minimize(opt_func, initial_guess,
                    bounds=[(-np.pi/2, np.pi/2)]*3,
                    **min_kwargs)

    
    
# * Rotational coherence
# Retrieve rotational coherence expressions from Pathways.
F = Function('F')

def rot_expression(state: RotState, jref: int, higher_order: str=None) -> Basic:
    """Rotational energy expression relative to ``jref``."""
    dj = state.j-jref

    globs = globals()
    nu = globs['nu'+str(state.nu)]
    is_symtop = getattr(state, "k", None)
    if higher_order == 'exact':
        if is_symtop:
            expr = F(nu, J_i+dj, K)
        else:
            expr = F(nu, J_i+dj)
    else:
        expr = nu + globs['B'+str(state.nu)]*(J_i+dj)*(J_i+dj+1)
        if is_symtop:
            expr += globs['BK'+str(state.nu)]*K**2
        if higher_order == 'DJ':
            expr += globs['DJ'+str(state.nu)]*((J_i+dj)*(J_i+dj+1))**2

    return expr


def rcs_expression(pair: Tuple[RotState, RotState],
                   jref: int, higher_order: str=None) -> Basic:
    """Return rotational beat expression."""
    return factor(rot_expression(pair[1], jref, higher_order)-
                  rot_expression(pair[0], jref, higher_order), deep=True)

honl_london_factors = {0: {-1: (Jpp+Kpp)*(Jpp-Kpp)/Jpp,
                            0: (2*Jpp+1)*Kpp**2/Jpp/(Jpp+1),
                            1: (Jpp+1+Kpp)*(Jpp+1-Kpp)/(Jpp+1)}}


def honl_london(pair: Tuple[RotState, RotState], jref: int) -> Basic:
    """Return Honl-London factor in terms of J_i = jref."""
    if 'k' in pair[0]._fields:
        K = Kpp
    else:
        K = 0

    dj = pair[0].j-jref

    return honl_london_factors[0][pair[1].j-pair[0].j].\
        subs(Kpp, K).subs(Jpp, J_i+dj)


def vib_expr(pair: Tuple[RotState, RotState]) -> Basic:
    """Return mu_ij symbol for `pair`."""
    nus = sorted((pair[0].nu, pair[1].nu))
    globs = globals()
    return globs['mu'+str(nus[0])+str(nus[1])]


def honl_london_pw(transitions: Sequence[Tuple[RotState, RotState]],
                   jref: int) -> Basic:
    """Four-fold Honl-London factor for a pathway."""
    return factor(reduce(mul, [honl_london(t, jref) for t in transitions], 1), deep=True)


def amplitude_expr_pw(pw1, norm: bool=False) -> Basic:
    """Pathway amplitude expression.

    Without equilibrium population (Boltzmann) factor and R-factor. Includes
    isotropy factor. For a quantity roughly equal for each J_i state, multiply
    by (2*J_i+1)**(1/2).

    Parameters
    ----------
    pw1 : Pathway
    norm : bool
        Flip the sign for negative probe pathway.

    Returns
    -------
    Basic
        Sympy pathway amplitude expression.
    """
    # Need to include the sign of rmu.
    const = -(I/hbar)**(len(pw1.transitions)-1)*pw1.leaf.total_side()/sqrt(2*J_i+1)
    vibs = factor(reduce(mul, [vib_expr(t) for t in pw1.transitions], 1), deep=True)
    rots = sqrt(honl_london_pw(pw1.transitions, pw1.leaf.root.ket.j))

    rot_sign = 1
    sides = [li.side for li in pw1.leaf.interactions()]
    for pair, side in zip(pw1.transitions, sides):
        if side == pw.Side.BRA:
            pair = pair[::-1]
        rot_sign *= -1 if pair[0].j < pair[1].j else 1

    if norm:
        probe = pw1.coherences[2]
        if (probe[0].nu < probe[1].nu) or\
           (probe[0].nu == probe[1].nu and probe[0].j < probe[1].j):
            rot_sign *= -1

    return const*vibs*rots*rot_sign
