r"""Associate polarizations, angular momentum factors and pathways amplitudes
with pathways. The time- and frequency-dependent signals are calculated in
:mod:`rotsim2d.propagate`.

:class:`Pathway` represents a double-sided Feynmann pathway for third-order
rovibrational excitation, without associating it to any specific
molecule. :class:`DressedPathway` specializes the object to specific vibrational
mode at a given temperature.  The main quantity, calculable with
:meth:`DressedPathway.amplitude`, is the pathway amplitude given by:

.. math::

    A(\widetilde{\theta}, \widetilde{J}) = (-1)^\kappa \frac{i}{\hbar^3} \langle O_{ijkl}\rangle
    = (-1)^\kappa \frac{i}{\hbar^3} \frac{N_{\eta_i J_i}}{N} \langle T^{(0)}_0(\eta_iJ_i)^\dagger\rangle R^{(0)}_0(\widetilde{\theta}, \widetilde{J}) \langle \nu_i J_i\|T^{(0)}(\tilde{\mu})\|\nu_i J_i\rangle,

where :math:`\widetilde{\theta}` and :math:`\widetilde{J}` are the sequences of
polarization angles and J values for the pathway. :math:`N_{\eta_i J_i}`
includes all relevant degeneracy factors (nuclear spin, rovibrational symmetry)
and :math:`\langle T^{(0)}_0(\eta_iJ_i)^\dagger\rangle=(2J_i+1)^{-1/2}`.

The macroscopic polarization is related to the pathway amplitude by:

.. math::

    \vec{\epsilon}_4\cdot P^{(n)} =  \frac{1}{8} N \mathcal{E}_1 \mathcal{E}_2 \mathcal{E}_3 A(\widetilde{\theta}, \widetilde{J}) \mathcal{I}(\widetilde{\omega}\text{ or } \widetilde{t}),

where :math:`N` is the number density and
:math:`\mathcal{I}(\widetilde{\omega}\text{ or } \widetilde{t})` determines the
n-dimensional frequency- or time-domain response. The absorption coefficient for the probe (in the Lambert-Beer law sense) is then given by:

.. math::

    \alpha_{\mathrm{probe}} = \frac{1}{4}N \frac{k^{\mathrm{probe}}_0}{\epsilon_0} \mathcal{E}_1 \mathcal{E}_2 A(\widetilde{\theta}, \widetilde{J}) \mathcal{I}(\widetilde{\omega}\text{ or } \widetilde{t}).

"""
import collections.abc as abc
import json
from abc import ABCMeta, abstractmethod
from collections import namedtuple
from math import isclose
from pathlib import Path
from typing import (Iterable, List, Mapping, Optional, Sequence, Tuple, Union, Dict, Any)

import h5py
import molspecutils.molecule as mol
import numpy as np
import scipy.constants as C
from molspecutils.molecule import CH3ClAlchemyMode, COAlchemyMode

import rotsim2d.couple as cp
import rotsim2d.pathways as pw
import rotsim2d.utils as u
from rotsim2d.couple import four_couple
from rotsim2d.pathways import KetBra, KSign, Side, gen_pathways

#: Spectroscopic notation for transitions
dj_to_letter = {-2: "O", -1: "P", 0: "Q", 1: "R", 2: "S"}
def abstract_line_label(pair: Tuple[mol.RotState], vib=False) -> str:
    ":meta private:"
    pair = sorted(pair, key=lambda x: x.nu)
    dnu = abs(pair[1].nu-pair[0].nu)
    label = str(dnu) if dnu>1 else ""
    label += dj_to_letter[pair[1].j-pair[0].j]
    if vib and (pair[1].nu > 1 or pair[0].nu > 1):
        label = str(max(pair[1].nu, pair[0].nu))+label

    return label


def abstract_format(dnu: int, dj: int):
    ":meta private:"
    if dnu==0 and dj==0:
        return "0"
    return str(dnu)+dj_to_letter[dj]


def abstract_state_label(state: mol.RotState, ref_state: mol.RotState) -> str:
    ":meta private:"
    return abstract_format(state.nu-ref_state.nu, state.j-ref_state.j)


def abstract_pair_label(pair: Tuple[mol.RotState], ref_state: mol.RotState) -> str:
    ":meta private:"
    return "|{:s}><{:s}|".format(
        abstract_state_label(pair[0], ref_state),
        abstract_state_label(pair[1], ref_state))

class Pathway:
    r"""Collect information on a pathway based on KetBra tree leaf without
    specializing it to any specific vibrational mode of a molecule.

    Attributes
    ----------
    leaf: rotsim2d.pathways.KetBra
        The leaf used to construct this pathway.
    coherences: list of tuple
        Coherences created by light-matter interactions as a list of pairs of
        :class:`molspecutils.molecule.RotState`.
    transitions: list of tuple
        List of pairs of :class:`molspecutils.molecule.RotState` representing
        transitions between states.
    js: tuple of int
        Arguments of G factor.
    light_inds: tuple of int
        Ordering of polarization vectors in four-fold dipole operator.
    const: complex
        :math:`(-1)^\kappa\left(\frac{i}{2\hbar}\right)^n`, where `n` is the
        order of interaction (usually 3) and :math:`\kappa` is the sign factor
        due to multiple interactions on either ket- or bra-side of the density
        matrix.
    geo_label: str
        Short hand notation for `js`, see module-level description.
    trans_label: str
        Short hand notation for series of transitions in the pathway.
    trans_label_deg: str
        Ambiguous version of `trans_label`.
    tw_coherence: bool
        Whether the molecule is in coherent state after second interaction.
    peak: tuple of str
        Pair of strings representing coherences created by first and third
        interaction (2D peak label).
    peak_label: str
        Two-letter identifier of the 2D peak.
    abstract_peak: tuple of str
        Same as `peak` but using P-, Q-, R-branch notation instead of absolute
        `J` numbers.
    experimental_label: str
        The kind of interaction: ground-state hole-burning, excited-states
        absorption, etc.
    colors: int
        Number of different optical frequencies requried to produce this
        pathway.
    """
    fields = ['leaf', 'coherences', 'transitions', 'js', 'light_inds', 'const',
              'geo_label', 'trans_label', 'trans_label_deg', 'tw_coherence',
              'peak', 'peak_label', 'abstract_peak', 'experimental_label',
              'colors']

    def __init__(self, leaf: KetBra):
        self.leaf = leaf
        self.isotropy = 1/np.sqrt(2*leaf.root.ket.j+1)
        self.coherences, self.transitions, wbras = leaf._pathway_info()
        self.js = tuple(x[0] for x in wbras)
        self.light_inds = tuple(x[1] for x in wbras)
        self.const = -(1.0j/C.hbar)**(len(self.transitions)-1)*\
            self.leaf.total_side()
        self.tw_coherence = not KetBra(*self.coherences[1]).is_diagonal()
        self.peak = ("|{:s}><{:s}|".format(self.coherences[0][0].name,
                                           self.coherences[0][1].name),
                     "|{:s}><{:s}|".format(self.coherences[2][0].name,
                                           self.coherences[2][1].name))
        self.abstract_peak = (abstract_pair_label(self.coherences[0],
                                                  self.leaf.root.ket),
                              abstract_pair_label(self.coherences[2],
                                                  self.leaf.root.ket))

    def __eq__(self, o):
        if not isinstance(o, Pathway):
            return NotImplemented
        return tuple(self.leaf.to_statelist()) == tuple(o.leaf.to_statelist())

    def __hash__(self):
        return hash(tuple(self.leaf.to_statelist()))

    @property
    def geo_label(self):
        ":meta private:"
        return geometric_labels[tuple(j-self.js[0] for j in self.js)]

    @property
    def trans_label(self):
        ":meta private:"
        return ''.join((self._trans_label(i) for i in (0, 1, 2)))

    @property
    def trans_label_deg(self):
        ":meta private:"
        return ''.join((self._trans_label(i, False) for i in (0, 1, 2)))

    @property
    def peak_label(self):
        ":meta private:"
        return abstract_line_label(self.coherences[0], True)+\
            '-' + abstract_line_label(self.coherences[2], True)

    @property
    def experimental_label(self):
        ":meta private:"
        if self.leaf.is_esa():
            return "Excited-state pump-probe"
        if self.leaf.is_gshb():
            return "Ground-state hole-burning"
        if self.leaf.is_sep():
            return "Stimulated-emission pumping"
        if self.leaf.is_doublequantum():
            return "Double quantum"

    @property
    def colors(self):
        ":meta private:"
        return self.leaf.color_tier()

    def _trans_label(self, i: int, unique: bool=True):
        labels = {-1: 'P', 0: 'Q', 1: 'R'}
        trans = sorted(self.transitions[i], key=lambda x: x.nu)
        side = self.leaf.sides()[i]

        label = labels[trans[1].j-trans[0].j]
        if unique:
            if trans[1].nu > 1:
                label = str(trans[1].nu)+label
            if side == Side.BRA:
                label = '('+ label +')'

        return label

    def _phi_angles(self, theta_angles: Union[Mapping, Sequence]) -> List[float]:
        """Order pulse/detection angles to evaluate R-factor."""
        if isinstance(theta_angles, abc.Sequence):
            theta_angles = dict(zip(('omg1', 'omg2', 'omg3', 'mu'),
                                    theta_angles))
        ints = self.leaf.interactions()

        return [theta_angles[ints[i].name] for i in self.light_inds]

    def geometric_factor(self, relative: bool=False,
                         angles: Optional[Union[Sequence, Dict]]=None) -> float:
        r"""Geometric R-factor for pathway intensity for isotropic initial density
        matrix:

        .. math::

            R^{(0)}_0(\epsilon_i^\ast, \epsilon_j, \epsilon_k^\ast, \epsilon_l; J_i, J_j, J_k, J_l) = \sum_{k=0}^2  T^{(0)}_0(\epsilon_i^\ast, \epsilon_j, \epsilon_k^\ast, \epsilon_l; k)G(J_i, J_j, J_k, J_l; k)

        Parameters
        ----------
        relative : bool
            Gives value relative to XXXX polarization.

        Returns
        -------
        float
            Purely J- and polarization-dependent part of the response.
        """
        if angles is None:
            angles = [0.0]*4
        ret = four_couple(self.js, self._phi_angles(angles))
        if relative:
            ret /= four_couple(self.js, [0.0]*4)

        return ret

    def gfactors(self):
        r"""Geometric factors for `k=0,1,2`:

        .. math::

            G(J_i, J_j, J_k, J_l; k) = (2k+1)\begin{Bmatrix} k & k & 0\\ J_i & J_i & J_k \end{Bmatrix} \begin{Bmatrix} 1 & 1 & k\\ J_k & J_i & J_j \end{Bmatrix} \begin{Bmatrix} 1 & 1 & k\\ J_k & J_i & J_l \end{Bmatrix}.
        """
        js = list(self.js)
        return (cp.G(*(js + [0])), cp.G(*(js + [1])), cp.G(*(js + [2])))

    def T00s(self, angles):
        r"""Polarization tensor components :math:`T^{(0)}_0(\epsilon_i^{\ast}, \epsilon_j, \epsilon_k^{\ast}, \epsilon_l; k)` for `k=0,1,2`
        """
        angles = self._phi_angles(angles)
        return (cp.T00(*(angles + [0])), cp.T00(*(angles + [1])), cp.T00(*(angles + [2])))

    @classmethod
    def from_kb_tree(cls, kb_tree: KetBra) -> List["Pathway"]:
        """Make a list of Pathway's from KetBra tree."""
        return [cls(leaf) for leaf in kb_tree.root.leaves]

    @classmethod
    def from_kb_list(cls, kb_list: Sequence[KetBra]) -> List["Pathway"]:
        """Make a list of Pathway's from KetBra list."""
        return sum((cls.from_kb_tree(kb_tree) for kb_tree in kb_list), [])

    def custom_str(self, fields=None):
        fields = fields if fields is not None else self.fields
        s = ', '.join(["{:s}={!s}".format(f, getattr(self, f)) for f in fields])

        return "Pathway({:s})".format(s)

    @staticmethod
    def _ketbra_symbols(kb: KetBra) -> Tuple[str]:
        l, r = '  ', '  '
        if kb.parent:
            arr = '->' if kb.parent.sign == KSign.POS else '<-'
            if kb.parent.side == Side.KET:
                l = arr
            else:
                r = arr
        return l, r

    def print_diagram(self, abstract=False, print=print):
        """Pretty print double-sided Feynmann diagram.

        Parameters
        ----------
        abstract: bool
            Use spectroscopic transition notation (P, Q, R) instead of absolute
            J values.
        """
        for kb in self.leaf.ketbras()[::-1]:
            l, r = Pathway._ketbra_symbols(kb)
            if abstract:
                print("{l:s}|{ket:2s}><{bra:>2s}|{r:s}".format(
                    l=l, r=r, ket=abstract_state_label(kb.ket, kb.root.ket),
                    bra=abstract_state_label(kb.bra, kb.root.ket)))
            else:
                print(f"{l}|{kb.ket.name}><{kb.bra.name}|{r}")

    def _tw_pprint(self, end='\n', print=print):
        if self.tw_coherence:
            import rotsim2d.symbolic.functions as sym
            tw_coherence = sym.rcs_expression(
                self.coherences[1], self.leaf.root.ket.j)
        else:
            tw_coherence = False
        print(f"Coherence during waiting time: {tw_coherence!r}", end=end)

    def pprint(self, abstract=False, angles=None, print=print):
        """Pretty print this pathway.

        Parameters
        ----------
        abstract: bool
            Use spectroscopic transition notation (P, Q, R) instead of absolute
            J values.
        """
        print("diagram:")
        self.print_diagram(abstract=abstract, print=print)
        print(f"G-factor label: {self.geo_label}")
        print("G-factors: {!s}".format(self.gfactors()))
        print("Total side: {!s}".format(self.leaf.total_side()))
        print("4-fold amplitude: {!s}".format(self.const))
        print(f"Transition chain label: {self.trans_label}")
        if angles is not None:
            print("Polarizations components: {!s}".format(self.T00s(angles=angles)))
            print("R-factor value: {:f}".format(self.geometric_factor(angles=angles)))
        print(f"Experimental label: {self.experimental_label}")
        print("Colors: {:d}".format(self.leaf.color_tier()))
        self._tw_pprint(print=print)

    def __repr__(self):
        return f"Pathway(leaf={self.leaf!r})"


class NDResonance(metaclass=ABCMeta):
    @abstractmethod
    def nu(self, i: int) -> float:
        """Frequency of `i` resonance."""
        pass

    @abstractmethod
    def gamma(self, i: int) -> float:
        """Decay rate of `i` resonance."""
        pass

    @abstractmethod
    def intensity(self, tw: Optional[float]=None, angles=None) -> float:
        """Amplitude of the whole resonance."""
        pass


class DressedPathway(Pathway, NDResonance):
    r"""Excitation pathway specialized to a vibrational mode.

    Parameters
    ----------
    leaf
        Leaf of :class:`rotsim2d.pathways.KetBra` excitation tree.
    vib_mode
        Object implementing :class:`molspecutils.molecule.VibrationalMode`
        interface.
    T
        Temperature in Kelvin.

    Attributes
    ----------
    const: complex
        The `const` factor of :class:`Pathway` multiplied by the fractional
        population of the initial state of the pathway in thermal equilibrium,
        :meth:`molspecutils.molecule.VibrationalMode.equilibrium_pop`, and by the
        four-fold reduced matrix element, :math:`\langle \nu_i J_i\|T^{(0)}(\tilde{\mu})\|\nu_i J_i\rangle`:

        .. math::

            \langle \nu_i J_i\|T^{(0)}(\tilde{\mu})\|\nu_i J_i\rangle = \langle \nu_i J_i\|\mu_1\|\nu_1 J_1\rangle \langle \nu_1 J_1\|\mu_2\|\nu_2 J_2\rangle \langle \nu_2 J_2\|\mu_3\|\nu_3 J_3\rangle \langle \nu_3 J_3\|\mu_4\|\nu_i J_i\rangle
    """
    def __init__(self, leaf: KetBra, vib_mode: mol.VibrationalMode, T: float):
        Pathway.__init__(self, leaf)
        self.vib_mode = vib_mode
        self.T = T
        sides = [li.side for li in leaf.interactions()]
        for pair, side in zip(self.transitions, sides):
            if side == Side.BRA:
                pair = pair[::-1]
            self.const *= vib_mode.mu(pair)
        self.const *= vib_mode.equilibrium_pop(self.leaf.root.ket, T)

    def __eq__(self, o):
        if not isinstance(o, DressedPathway):
            return NotImplemented
        return Pathway.__eq__(self, o) and isclose(self.T, o.T) and self.vib_mode == o.vib_mode

    def __hash__(self):
        return hash((
            tuple(self.leaf.to_statelist()),
            self.T,
            id(self.vib_mode)))

    def nu(self, i: int) -> float:
        """Frequency of `i`-th coherence."""
        return self.vib_mode.nu(self.coherences[i])

    def pump_fraction(self, E12: float) -> float:
        """Fraction of initial population excited by pump pulses.

        Parameters
        ----------
        E12 : float
           Electric field integral for the pump pulses (same for both pulses).
        """
        rmu2 = self.vib_mode.mu(self.transitions[0])*\
            self.vib_mode.mu(self.transitions[1])
        deg_fac = np.sqrt(2*self.leaf.root.ket.j+1)

        return 3/4/C.hbar**2*deg_fac*rmu2*E12**2*\
            self.geometric_factor()

    def gamma(self, i: int) -> float:
        """Pressure-broadening coefficient of `i`-th coherence."""
        return self.vib_mode.gamma(self.coherences[i])

    def intensity(self, tw: Optional[float]=None, angles=None) -> float:
        r"""Amplitude of the pathway, given by the product of isotropic coefficient of
        the initial density matrix, :math:`\langle T^{(0)}_0(\eta_i
        J_i)^\dagger\rangle=(2J_i+1)^{-1/2}`, :attr:`const` and :meth:`Pathway.geometric_factor`.
        """
        ret = self.isotropy*self.const*self.geometric_factor(angles=angles)
        if tw is not None:
            ret *= np.exp(-2.0*np.pi*tw*(1.0j*self.nu(1)))

        return ret

    def _tw_pprint(self, print=print):
        Pathway._tw_pprint(self, end='', print=print)
        if self.tw_coherence:
            print(", {:.2f} cm-1".format(u.nu2wn(self.nu(1))))
        else:
            print()

    def pprint(self, abstract=False, angles=None, print=print):
        """Pretty print this pathway.

        Parameters
        ----------
        abstract: bool
            Use spectroscopic transition notation (P, Q, R) instead of absolute
            J values.
        """
        Pathway.pprint(self, abstract=abstract, angles=angles, print=print)
        print("Gammas: {:.3e}, {:.3e}, {:.3e} cm-1".format(
            self.gamma(0)/C.c/100.0,
            self.gamma(1)/C.c/100.0,
            self.gamma(2)/C.c/100.0))

    @classmethod
    def from_kb_tree(cls, kb_tree: KetBra, vib_mode: mol.VibrationalMode,
                     T: float) -> List["DressedPathway"]:
        """Make a list of DressedPathway's from KetBra tree."""
        dp_list = []
        for leaf in kb_tree.root.leaves:
            try:
                dp_list.append(cls(leaf, vib_mode, T))
            except mol.MissingStateError:
                continue

        return dp_list

    @classmethod
    def from_kb_list(cls, kb_list: Sequence[KetBra], vib_mode: mol.VibrationalMode,
                     T: float) -> List["DressedPathway"]:
        """Make a list of DressedPathway's from KetBra list."""
        return sum((cls.from_kb_tree(kb_tree, vib_mode, T) for kb_tree in kb_list), [])

    base_params_dict: Dict[str, Any] = {
        'isotopologue': 1,
    }

    @classmethod
    def from_params_dict(cls, params: Mapping) -> List["DressedPathway"]:
        """Make a list of DressedPathway's from dict of parameters."""
        fparams = cls.base_params_dict.copy()
        fparams.update(params)
        if fparams['molecule'] == 'CH3Cl':
            mode = CH3ClAlchemyMode(iso=fparams['isotopologue'])
            rotor = 'symmetric'
        elif fparams['molecule'] == 'CO':
            mode = COAlchemyMode(iso=fparams['isotopologue'])
            rotor = 'linear'
        else:
            raise ValueError("Invalid molecule")

        meths = []
        if "direction" in fparams:
            meths.append(getattr(pw, "only_"+fparams["direction"]))
        meths.extend([getattr(pw, meth) for meth in fparams['filters']])
        kbs = gen_pathways(
            range(fparams['jmax']), meths=meths, rotor=rotor,
            kiter_func=fparams['kiter'])

        return cls.from_kb_list(kbs, mode, fparams['T'])

    def __repr__(self):
        return f"DressedPathway(leaf={self.leaf!r}, vib_mode={self.vib_mode!r}, T={self.T!r})"


class DressedLeaf:
    """Assign polarizations in proper order, line intensities, coherence frequencies
    and allow for calculating the intensity of the pathway. Do everything except
    calculating the associated line shape/response time dependence."""
    fields = ['peak', 'const', 'cohs', 'js', 'pols', 'tw_coherence', 'tw_peak', 'diagonal', 'geo_label', 'leaf']

    def __init__(self, leaf: KetBra, vib_mode: mol.VibrationalMode, T: float):
        self.leaf = leaf
        kb_series = leaf.ketbras()
        wkets, wbras = [], []
        self.cohs = []
        self.const = np.complex(1.0)
        for i in range(1, len(kb_series)):
            kb, kbp = kb_series[i], kb_series[i-1]

            # dipole interaction reduced matrix element
            if kb.parent.side is Side.KET:
                pair = (kbp.ket, kb.ket)
                wkets.insert(0, (kb.ket.j,  kb.parent.angle))
            else:
                pair = (kbp.bra, kb.bra)
                wbras.append((kb.parent.parent.bra.j, kb.parent.angle))
            mu = vib_mode.mu(pair) # do permutation in vib_mode.mu

            if kb.parent.readout:
                self.const *= mu*vib_mode.equilibrium_pop(leaf.root.ket, T)
                wbras.extend(wkets)
                self.js = tuple(x[0] for x in wbras)
                self.pols = tuple(x[1] for x in wbras)
                break
            self.const *= 1.0j/C.hbar*kb.parent.side*mu # this probably should be divided by half

            # collect pairs
            self.cohs.append((kb.bra, kb.ket))

        self.vib_mode = vib_mode
        self.tw_coherence = not self.nu(1) == 0.0
        self.tw_peak = kb_series[2].name
        self.peak = (kb_series[1].name, kb_series[3].name)
        self.diagonal = self.peak[0] == self.peak[1]
        self.geo_label = geometric_label(self)
        self.isotropy = 1/np.sqrt(2*leaf.root.ket.j+1)

    def nu(self, i):
        return self.vib_mode.nu(self.cohs[i])

    def gamma(self, i):
        return self.vib_mode.gamma(self.cohs[i])

    def intensity(self, tw: Optional[float]=None) -> float:
        """Intensity of the pathway."""
        ret = self.isotropy*self.const*four_couple(self.js, self.pols)
        if tw is not None and self.tw_coherence:
            ret *= np.exp(-2.0*np.pi*tw*(1.0j*self.nu(1)))

        return ret

    def geometric_factor(self, relative: bool=False) -> float:
        """Geometric factor for pathway intensity.

        Parameters
        ----------
        relative : bool
            Gives value relative to XXXX polarization.

        Returns
        -------
        float
            Purely J- and polarization-dependent part of the response.
        """
        ret = four_couple(self.js, self.pols)
        if relative:
            ret /= four_couple(self.js, [0.0]*4)

        return ret

    def custom_str(self, fields=None):
        fields = fields if fields is not None else self.fields
        s = ', '.join(["{:s}={!s}".format(f, getattr(self, f)) for f in fields])

        return "DressedLeaf({:s})".format(s)

    def print_diagram(self):
        for ket, bra in self.cohs[::-1]:
            print(f"|{ket.name}><{bra.name}|")

    def pprint(self):
        print('diagram:')
        self.print_diagram()
        print(f'geometric factor: {self.geo_label}')

    def __str__(self):
        return self.custom_str()

    def __repr__(self):
        return self.custom_str()


def split_by_js(kbl: Iterable[DressedLeaf]):
    ret = {}
    for dl in kbl:
        ret.setdefault(dl.js, []).append(dl)
    return ret


def perm_pols(pols: Tuple) -> Tuple:
    return (pols[2], pols[3], pols[0], pols[1])


def perm_js(js: Tuple) -> Tuple:
    return (js[0], js[3], js[2], js[1])


def undegenerate_js(ljs: Sequence[Tuple]) -> List[Tuple]:
    nondeg = []
    for jset in ljs:
        if jset not in nondeg and perm_js(jset) not in nondeg:
            nondeg.append(jset)

    return nondeg


def split_by_pols(kbl: Iterable[DressedLeaf]):
    ret = {}
    for dl in kbl:
        pols_set = frozenset((dl.pols, perm_pols(dl.pols)))
        ret.setdefault(pols_set, []).append(dl)
    return ret


def split_by_pols_js(kbl: Iterable[DressedLeaf]):
    ret = {}
    for dl in kbl:
        pols_set = frozenset((dl.pols, perm_pols(dl.pols)))
        ret.setdefault((dl.js, pols_set), []).append(dl)
    return ret


geometric_labels = {
    # linear labels
    # (0, -1, -2, -1): 'O',
    # (0,  1,  2,  1): 'S',
    # (0, -1,  0, -1): 'PP',
    # (0,  1,  0,  1): 'RR',
    # (0,  1,  0, -1): 'RP',
    # (0, -1,  0,  1): 'PR',      # RP
    (0, -1, -2, -1): 'PPR',
    (0,  1,  2,  1): 'RRP',
    (0, -1,  0, -1): 'PRP',
    (0,  1,  0,  1): 'RPR',
    (0,  1,  0, -1): 'RPP',
    (0, -1,  0,  1): 'PRR',      # RP
    # symmetric top labels
    (0, -1, -1, -1): 'PQQ',
    (0,  0, -1, -1): 'QPQ',
    (0,  0, -1,  0): 'QPR',
    (0,  0,  0, -1): 'QQP',
    (0,  0,  0,  0): 'QQQ',
    (0,  0,  1,  0): 'QRP',
    (0,  1,  0,  0): 'RPQ',
    (0,  1,  1,  0): 'RQP',
    (0,  1,  1,  1): 'RQQ',
    (0, -1, -1,  0): 'PQR',     # QPQ
    (0, -1,  0,  0): 'PRQ',     # QQP
    (0,  0,  0,  1): 'QQR',     # RPQ
    (0,  0,  1,  1): 'QRQ',     # RQP
}

def geometric_label(dl: DressedLeaf):
    js = dl.js; ji = js[0]
    js = tuple(j-ji for j in js)
    return geometric_labels[js]


def split_by_pols_highjs(kbl: Iterable[DressedLeaf]):
    sopr = {'RRP', 'PPR', 'PRR', 'RPP'}
    rrpp = {'RPR', 'PRP'}
    ret = {}
    for dl in kbl:
        pols_set = frozenset((dl.pols, perm_pols(dl.pols)))
        if geometric_label(dl) in sopr:
            cat = 'sopr'
        elif geometric_label(dl) in rrpp:
            cat = 'rrpp'
        ret.setdefault((cat, pols_set), []).append(dl)
    return ret


def split_by_peaks(kbl: Iterable[Pathway], abstract: bool=False)\
    -> Dict[str, List[Pathway]]:
    ret = {}
    for dl in kbl:
        if abstract:
            ret.setdefault(dl.abstract_peak, []).append(dl)
        else:
            ret.setdefault(dl.peak, []).append(dl)

    return ret


def pprint_dllist(dllist, abstract=False, angles=None, print=print):
    for i, dl in enumerate(dllist):
        if i == 0:
            print('-'*10)
            if isinstance(dl, DressedPathway):
                print('pump = {:.2f} cm-1, probe = {:.2f} cm-1'.format(
                    u.nu2wn(dl.nu(0)), u.nu2wn(dl.nu(2))))
            dl.pprint(abstract=abstract, angles=angles, print=print)
        else:
            print()
            dl.pprint(abstract=abstract, angles=angles, print=print)
        if i == len(dllist)-1:
            print('-'*10)
            print()


def print_dl_dict(dldict, fields=None):
    for k in dldict:
        print(k)
        for dl in dldict[k]:
            print("   ", dl.custom_str(fields=fields))


def print_dl_tuple_dict(dldict, fields=None):
    for k in dldict:
        print(k)
        for dl, a in dldict[k]:
            print("   {:s}, {:f}".format(dl.custom_str(fields=fields), a))


def dress_pws(pws, vib_mode, T):
    return sum(([DressedLeaf(l, vib_mode, T) for l in root.leaves] for root in pws), [])


# * Peaks without line shapes
Peak2D = namedtuple("Peak2D", "pump_wl probe_wl sig peak")

class Peak2DList(list):
    """List of 2D peaks with easy access to pump, probe frequencies and peak
    intensities.

    Calling :meth:`copy` or slicing returns a regular list.
    """
    def __init__(self, iterable=(), quantity='amplitudes'):
        super().__init__(iterable)
        self.quantity = quantity

    @property
    def pumps(self):
        """List of pump frequencies."""
        return [peak.pump_wl for peak in self]

    @property
    def probes(self):
        """List of probe frequencies."""
        return [peak.probe_wl for peak in self]

    @property
    def sigs(self):
        """Peak amplitude--sum of :meth:`DressedPathway.intensity` over all pathways
        contributing to a 2D peak.
        """
        return [peak.sig for peak in self]

    @property
    def peaks(self):
        """Peak strings."""
        return [peak.peak for peak in self]

    @staticmethod
    def _sort_func(peak):
        return abs(peak.sig)

    def sort_by_sigs(self):
        """Sort peaks by amplitude.

        Ensures that strong peaks are not covered by weak ones in scatter plot.
        """
        self.sort(key=self._sort_func)

    def normalize_sig_signs(self) -> 'Peak2DList':
        """Flip signs of amplitudes with negative probe frequencies.

        The Maxwell wave equation flips the sign of the negative frequency
        material polarization component, such that the positive and negative
        components add up to a real cosine wave. At the level of TD perturbation
        theory the signs are oppposite and summing corresponding rephasing and
        non-rephasing pathway amplitudes (or a pathway with its complex
        conjugate) gives zero. This means that if you want to estimate some
        signal amplitude/intensity based on :class:`Peak2DList` amplitudes, you
        might get a zero unless you normalize the signs.
        """
        pnorm = Peak2DList(quantity=self.quantity)
        for p in self:
            if p.probe_wl < 0.0:
                pnorm.append(p._replace(sig=-p.sig))
            else:
                pnorm.append(p)

        return pnorm

    def get_by_peak(self, peak: Tuple[str, str]) -> Peak2D:
        """Return peak with `peak` attribute equal to `peak` argument."""
        for p in self:
            if p.peak == peak:
                return p
        else:
            raise IndexError("'{!s}' not found in the list".format(peak))

    def to_abs_coeff(self, E12: float, conc: float) -> "Peak2DList":
        """Convert peaks of 'peak_intensity' type to abs. coefficients.

        Assumes self-heterodyne probe detection.

        Parameters
        ----------
        E12 : float
            Product of field integrals of pump pulses.
        conc : float
            Concentration in 1/m**3.
        """
        if self.quantity != 'peak_intensity':
            raise ValueError("This method only works for 'peak_intensity' type peaks.")

        pabs = Peak2DList(quantity=self.quantity)
        for p in self:
            pabs.append(p._replace(sig=p.sig/np.pi*conc*E12))

        return pabs

    def to_file(self, path: Union[str, Path], metadata: Optional[Dict]=None):
        """Save peak list to HDF5 file."""
        with h5py.File(path, mode='w') as f:
            f.create_dataset("pumps", data=self.pumps)
            f.create_dataset("probes", data=self.probes)
            f.create_dataset("sigs", data=self.sigs)
            f.create_dataset(
                "peaks", data=[json.dumps(peak) for peak in self.peaks])
            if metadata:
                f.attrs['metadata'] = json.dumps(metadata)

    @classmethod
    def from_file(cls, path: Union[str, Path]):
        """Read peak list from HDF5 file."""
        with h5py.File(path, mode='r') as f:
            pl = cls()
            for pu, pr, sig, peak in zip(
                    f['pumps'], f['probes'], f['sigs'], f['peaks']):
                pl.append(Peak2D(pu, pr, sig, tuple(json.loads(peak))))

        return pl


def run_peak_list(params: Mapping) -> Peak2DList:
    """Calculate list of peaks based on toml input data."""
    if params['spectrum']['type'] != 'peaks':
        raise ValueError("Wrong spectrum type requested.")

    dpws = DressedPathway.from_params_dict(params['pathways'])
    pl = peak_list(dpws, tw=params['spectrum']['tw']*1e-12,
                   angles=params['spectrum']['angles'])

    return pl


def peak_list(ll: List[DressedPathway], tw: Optional[float]=0.0,
              angles=None, return_dls: bool=False,
              quantity: str='amplitude', p: float=1.0) -> \
              Union[Peak2DList, Tuple[Peak2DList, List[List[DressedPathway]]]]:
    """Create a list of 2D peaks from a list of :class:`DressedPathway`.

    Optionally return sorted list of :class:`DressedPathway` corresponding to peaks.

    `quantity` can be:

    - 'amplitude', pathway amplitude as returned by
      :meth:`DressedPathway.intensity`,
    - 'line_intensity', HITRAN-analogous line intensity,
    - 'peak_intensity', 'line_intensity' divided by pressure width. Peak
       absorption value.
    """
    ll: dict = split_by_peaks(ll)
    pl = Peak2DList(quantity=quantity)
    dls = []
    for peak, dll in ll.items():
        pu, pr = u.nu2wn(dll[0].nu(0)), u.nu2wn(dll[0].nu(2))
        sig = 0.0
        for dl in dll:
            # with wigxjpf(300, 6):
            pre_sig = np.imag(dl.intensity(tw=tw, angles=angles))
            if quantity == 'line_intensity':
                pre_sig *= np.pi*2*np.pi*dll[0].nu(2)/4/C.epsilon_0/C.c
            elif quantity == 'peak_intensity':
                pre_sig *= np.pi*2*np.pi*dll[0].nu(2)/4/C.epsilon_0/C.c/dll[0].gamma(2)/p
            sig += pre_sig
        pl.append(Peak2D(pu, pr, sig, peak))
        if return_dls:
            dls.append(dll)
    if return_dls:
        pairs = sorted(zip(pl, dls), key=lambda x: abs(x[0].sig))
        return Peak2DList([x[0] for x in pairs]), [x[1] for x in pairs]

    pl.sort_by_sigs()
    return pl


def equiv_peaks(pw: Pathway, pl: Peak2DList, dll: Sequence[Sequence[Pathway]]) -> Peak2DList:
    """Return peaks from `pl` which are polarization-equivalent to `pw`."""
    new_pl = Peak2DList(quantity=pl.quantity)
    for peak, dp_list in zip(pl, dll):
        if any(dp.leaf.is_equiv_pathway(pw) for dp in dp_list):
            new_pl.append(peak)
    new_pl.sort_by_sigs()

    return new_pl


def split_by_equiv_peaks(det_angles: dict, pl: Peak2DList, dll: Sequence[Pathway]) -> dict:
    """Map zeroing angles from `det_angles` to peaks from `pl`.

    If `det_angles` is complete, e.g. a result of
    :func:`rotsim2d.symbolic.functions.detection_angles`, then this function
    will split *all* peaks from `pl` into disjoint sets distinguishable by
    polarization (under the polarization scheme assumed when constructing
    `det_angles`).

    Parameters
    ----------
    det_angles
        Map between detection angles and a list of pathways.
    pl
        Any peak list.
    dll
        A list of lists of pathways associated with each peak in `pl`.

    Returns
    -------
    equiv_peaks_dict
        Map between `det_angles` dict keys and peak lists constructed from `pl`.
    """
    equiv_peaks_dict = {}
    for angle in det_angles:
        equiv_peaks_dict[angle] = Peak2DList()
        for pw1 in det_angles[angle]:
            equiv_peaks_dict[angle].extend(equiv_peaks(pw1.leaf, pl, dll))
        equiv_peaks_dict[angle].sort_by_sigs()

    return equiv_peaks_dict
