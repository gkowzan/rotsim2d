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
import re
import json
from abc import ABCMeta, abstractmethod
from math import isclose
from pathlib import Path
from typing import (Iterable, List, Mapping, Optional, Sequence, Tuple, Union,
                    Dict, Any, Callable)
from attr import define, field

import h5py
import molspecutils.molecule as mol
import numpy as np
import scipy.constants as C
from molspecutils.molecule import CH3ClAlchemyMode, COAlchemyMode, RotState

import rotsim2d.couple as cp
import rotsim2d.pathways as pw
import rotsim2d.utils as u
from rotsim2d.couple import four_couple
from rotsim2d.pathways import KetBra, KSign, Side, gen_pathways

#: Spectroscopic notation for transitions
dj_to_letter = {-2: "O", -1: "P", 0: "Q", 1: "R", 2: "S"}


def abstract_line_label(pair: Tuple[mol.RotState, mol.RotState],
                        vib=False) -> str:
    "Assign branch label, e.g: P, 2P, Q, to ``pair`` representing a transition."
    pair = sorted(pair, key=lambda x: x.nu)
    dnu = abs(pair[1].nu-pair[0].nu)
    label = str(dnu) if dnu>1 else ""
    label += dj_to_letter[pair[1].j-pair[0].j]
    if vib and (pair[1].nu > 1 or pair[0].nu > 1):
        label = str(max(pair[1].nu, pair[0].nu))+label

    return label


def abstract_format(dnu: int, dj: int) -> str:
    "Abstract ket (or bra) label, e.g: 1P, 0Q."
    if dnu==0 and dj==0:
        return "0"
    return str(dnu)+dj_to_letter[dj]


def abstract_state_label(state: mol.RotState, ref_state: mol.RotState) -> str:
    "Abstract label for ``state`` relative to ``ref_state``, e.g: 1P, 0Q."
    return abstract_format(state.nu-ref_state.nu, state.j-ref_state.j)


def abstract_pair_label(pair: Tuple[mol.RotState, mol.RotState],
                        ref_state: mol.RotState) -> str:
    "Abstract label for ``pair`` relative to ``ref_state``, e.g: |1P><0|. "
    return "|{:s}><{:s}|".format(
        abstract_state_label(pair[0], ref_state),
        abstract_state_label(pair[1], ref_state))


_coh_conj_regex = re.compile(r"\|(.+)><(.+)\|")
def coherence_conj(coh_string: str) -> str:
    """Conjugate coherence string, `|ket><bra|`->`|bra><ket|`."""
    match = re.fullmatch(_coh_conj_regex, coh_string)
    if match:
        return "|{:s}><{:s}|".format(match.group(2), match.group(1))
    else:
        raise ValueError("Malformed peak string '{:s}'".format(coh_string))


def peak_conj(peak: Tuple[str, str]) -> Tuple[str, str]:
    """Conjugate last coherence string in peak identifier.

    See also
    --------
    Pathway.peak, Pathway.abstract_peak
    """
    return (peak[0], coherence_conj(peak[1]))


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

AnglesT = Union[Mapping[str, Any], Sequence]


class Pathway:
    r"""Collect information on a pathway based on KetBra tree leaf without
    specializing it to any specific vibrational mode of a molecule.

    Parameters
    ----------
    leaf : KetBra
        Leaf of a KetBra excitation tree.
    """
    fields = ['leaf', 'coherences', 'transitions', 'js', 'light_inds', 'const',
              'geo_label', 'trans_label', 'trans_label_deg', 'tw_coherence',
              'peak', 'peak_label', 'abstract_peak', 'experimental_label',
              'colors']

    def __init__(self, leaf: KetBra):
        self.leaf: KetBra = leaf
        "The leaf used to construct this pathway."
        self.isotropy: float = 1/np.sqrt(2*leaf.root.ket.j+1)
        c, t, wbras = leaf._pathway_info()
        self.coherences: List[Tuple[RotState, RotState]] = c
        """Coherences created by light-matter interactions."""
        self.transitions: List[Tuple[RotState, RotState]] = t
        """Time-ordered transitions between states for the pathway."""
        self.js: Tuple[int, ...] = tuple(x[0] for x in wbras)
        """Arguments for the G-factor."""
        self.light_inds: Tuple[int, ...] = tuple(x[1] for x in wbras)
        """Ordering of polarization vectors in four-fold dipole-operator."""
        self.const: complex = -(1.0j/C.hbar)**(len(self.transitions)-1)*\
            self.leaf.total_side()
        r""":math:`(-1)^\kappa\left(\frac{i}{2\hbar}\right)^n`, where `n` is the
        order of interaction (usually 3) and :math:`\kappa` is the sign factor
        due to multiple interactions on either ket- or bra-side of the density
        matrix."""
        self.tw_coherence: bool = not KetBra(*self.coherences[1]).is_diagonal()
        "Whether the molecule is in coherent state after second interaction."
        self.peak: Tuple[str, str] = (
            "|{:s}><{:s}|".format(self.coherences[0][0].name,
                                  self.coherences[0][1].name),
            "|{:s}><{:s}|".format(self.coherences[2][0].name,
                                  self.coherences[2][1].name))
        """Pair of strings representing coherences created by first and third
        interaction (2D peak label)."""
        self.abstract_peak: Tuple[str, str] = (
            abstract_pair_label(self.coherences[0],
                                self.leaf.root.ket),
            abstract_pair_label(self.coherences[2],
                                self.leaf.root.ket))
        """Same as `peak` but using P-, Q-, R-branch notation instead of absolute
        `J` numbers."""

    def __eq__(self, o):
        if not isinstance(o, Pathway):
            return NotImplemented
        return tuple(self.leaf.to_statelist()) == tuple(o.leaf.to_statelist())

    def __hash__(self):
        return hash(tuple(self.leaf.to_statelist()))

    @property
    def geo_label(self) -> str:
        """Shorthand notation for :attr:`js`."""
        return geometric_labels[tuple(j-self.js[0] for j in self.js)]

    @property
    def trans_label(self) -> str:
        "Short hand notation for series of transitions in the pathway."
        return ''.join((self._trans_label(i) for i in (0, 1, 2)))

    @property
    def trans_label_deg(self) -> str:
        "Ambiguous version of :attr:`trans_label`."
        return ''.join((self._trans_label(i, False) for i in (0, 1, 2)))

    @property
    def peak_label(self) -> str:
        "Two-letter identifier for the 2D peak."
        return abstract_line_label(self.coherences[0], True)+\
            '-' + abstract_line_label(self.coherences[2], True)

    @property
    def experimental_label(self) -> str:
        """The kind of interaction: ground-state hole-burning, excited-states
        absorption, etc."""
        if self.leaf.is_esa():
            return "Excited-state pump-probe"
        if self.leaf.is_gshb():
            return "Ground-state hole-burning"
        if self.leaf.is_sep():
            return "Stimulated-emission pumping"
        if self.leaf.is_doublequantum():
            return "Double quantum"

    @property
    def colors(self) -> int:
        """Number of different optical frequencies requried to produce this
        pathway."""
        return self.leaf.color_tier()

    def _trans_label(self, i: int, unique: bool=True) -> str:
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

    def _phi_angles(self, theta_angles: AnglesT) -> List[float]:
        """Order pulse/detection angles to evaluate R-factor."""
        if isinstance(theta_angles, abc.Sequence):
            theta_angles = dict(zip(('omg1', 'omg2', 'omg3', 'mu'),
                                    theta_angles))
        ints = self.leaf.interactions()

        return [theta_angles[ints[i].name] for i in self.light_inds]

    def geometric_factor(self, relative: bool=False,
                         angles: Optional[AnglesT]=None) -> float:
        r"""Geometric R-factor for pathway intensity for isotropic initial density
        matrix.

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
        r"""Geometric factors for `k=0,1,2`.

        .. math::

            G(J_i, J_j, J_k, J_l; k) = (2k+1)\begin{Bmatrix} k & k & 0\\ J_i & J_i & J_k \end{Bmatrix} \begin{Bmatrix} 1 & 1 & k\\ J_k & J_i & J_j \end{Bmatrix} \begin{Bmatrix} 1 & 1 & k\\ J_k & J_i & J_l \end{Bmatrix}.
        """
        js = list(self.js)
        return (cp.G(*(js + [0])), cp.G(*(js + [1])), cp.G(*(js + [2])))

    def T00s(self, angles: AnglesT):
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

    def custom_str(self, fields: Optional[Sequence[str]]=None) -> str:
        """String representation including only attributes in ``fields``."""
        fields = fields if fields is not None else self.fields
        s = ', '.join(["{:s}={!s}".format(f, getattr(self, f)) for f in fields])

        return "Pathway({:s})".format(s)

    @staticmethod
    def _ketbra_symbols(kb: KetBra) -> Tuple[str, str]:
        l, r = '  ', '  '
        if kb.parent:
            arr = '->' if kb.parent.sign == KSign.POS else '<-'
            if kb.parent.side == Side.KET:
                l = arr
            else:
                r = arr
        return l, r

    def print_diagram(self, abstract: bool=False, print: Callable=print):
        """Pretty print double-sided Feynmann diagram.

        Parameters
        ----------
        abstract
            Use spectroscopic transition notation (P, Q, R) instead of absolute
            J values.
        print
            Use this callable instead of :func:`print` builtin.
        """
        for kb in self.leaf.ketbras()[::-1]:
            l, r = Pathway._ketbra_symbols(kb)
            if abstract:
                print("{l:s}|{ket:2s}><{bra:>2s}|{r:s}".format(
                    l=l, r=r, ket=abstract_state_label(kb.ket, kb.root.ket),
                    bra=abstract_state_label(kb.bra, kb.root.ket)))
            else:
                print(f"{l}|{kb.ket.name}><{kb.bra.name}|{r}")

    def _tw_pprint(self, end: str='\n', print: Callable=print):
        if self.tw_coherence:
            import rotsim2d.symbolic.functions as sym
            tw_coherence = sym.rcs_expression(
                self.coherences[1], self.leaf.root.ket.j)
        else:
            tw_coherence = False
        print(f"Coherence during waiting time: {tw_coherence!r}", end=end)

    def pprint(self, abstract: bool=False,
               angles: Optional[AnglesT]=None,
               print: Callable=print):
        """Pretty print this pathway.

        Parameters
        ----------
        abstract
            Use spectroscopic transition notation (P, Q, R) instead of absolute
            J values.
        angles
            Evaluate polarization components and R-factor for these polarizations.
        print
            Use this callable instead of :func:`print` built-in.
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
    """Interface expected by functions in :module:`rotsim2d.propagate`."""
    @abstractmethod
    def nu(self, i: int) -> float:
        """Frequency of `i` resonance."""
        pass

    @abstractmethod
    def gamma(self, i: int) -> float:
        """Decay rate of `i` resonance."""
        pass

    @abstractmethod
    def amplitude(self, tw: Optional[float]=None, angles=None) -> complex:
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
    """
    def __init__(self, leaf: KetBra, vib_mode: mol.VibrationalMode, T: float):
        Pathway.__init__(self, leaf)
        self.vib_mode: mol.VibrationalMode = vib_mode
        """Vibrational mode associated with this pathway."""
        self.T: float = T
        """Temperature in Kelvin."""
        sides = [li.side for li in leaf.interactions()]
        for pair, side in zip(self.transitions, sides):
            if side == Side.BRA:
                pair = pair[::-1]
            self.const *= vib_mode.mu(pair)
        self.const: complex = self.const*vib_mode.equilibrium_pop(self.leaf.root.ket, T)
        r"""The `const` factor of :class:`Pathway` multiplied by the fractional
        population of the initial state of the pathway in thermal equilibrium,
        :meth:`molspecutils.molecule.VibrationalMode.equilibrium_pop`, and by the
        four-fold reduced matrix element, :math:`\langle \nu_i J_i\|T^{(0)}(\tilde{\mu})\|\nu_i J_i\rangle`:

        .. math::

            \langle \nu_i J_i\|T^{(0)}(\tilde{\mu})\|\nu_i J_i\rangle = \langle \nu_i J_i\|\mu_1\|\nu_1 J_1\rangle \langle \nu_1 J_1\|\mu_2\|\nu_2 J_2\rangle \langle \nu_2 J_2\|\mu_3\|\nu_3 J_3\rangle \langle \nu_3 J_3\|\mu_4\|\nu_i J_i\rangle"""

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
        E12
           Electric field integral for the pump pulses (same for both pulses).
        """
        rmu2 = self.vib_mode.mu(self.transitions[0])*\
            self.vib_mode.mu(self.transitions[1])
        deg_fac = np.sqrt(2*self.leaf.root.ket.j+1)

        return 3/4/C.hbar**2*deg_fac*rmu2*E12**2*\
            self.geometric_factor()

    def gamma(self, i: int) -> float:
        """Pressure-broadening coefficient of `i` coherence."""
        return self.vib_mode.gamma(self.coherences[i])

    def amplitude(self, tw: Optional[float]=None,
                  angles: Optional[AnglesT]=None) -> complex:
        r"""Amplitude of the pathway, given by the product of isotropic coefficient of
        the initial density matrix, :math:`\langle T^{(0)}_0(\eta_i
        J_i)^\dagger\rangle=(2J_i+1)^{-1/2}`, :attr:`const` and :meth:`Pathway.geometric_factor`.
        """
        ret = self.isotropy*self.const*self.geometric_factor(angles=angles)
        if tw is not None:
            ret *= np.exp(-2.0*np.pi*tw*(1.0j*self.nu(1)))

        return ret

    def _tw_pprint(self, print: Callable=print):
        Pathway._tw_pprint(self, end='', print=print)
        if self.tw_coherence:
            print(", {:.2f} cm-1".format(u.nu2wn(self.nu(1))))
        else:
            print()

    def pprint(self, abstract: bool=False,
               angles: Optional[AnglesT]=None,
               print: Callable=print):
        """Pretty print this pathway.

        Parameters
        ----------
        abstract
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


AbstractPathway = Union[Pathway, DressedPathway]

def split_by_js(kbl: Iterable[Pathway]) -> Dict[Tuple[int, ...], List[Pathway]]:
    """Collect pathways with the same `js`."""
    ret: Dict[Tuple[int, ...], List[Pathway]] = {}
    for dl in kbl:
        ret.setdefault(dl.js, []).append(dl)

    return ret


def perm_pols(pols: Tuple) -> Tuple:
    """Rearrangement of polarizations allowed by symmetry."""
    return (pols[2], pols[3], pols[0], pols[1])


def perm_js(js: Tuple) -> Tuple:
    """Rearrangement of J values allowed by symmetry."""
    return (js[0], js[3], js[2], js[1])


def undegenerate_js(ljs: Sequence[Tuple]) -> List[Tuple]:
    """Remove degenerate tuples of Js."""
    nondeg = []
    for jset in ljs:
        if jset not in nondeg and perm_js(jset) not in nondeg:
            nondeg.append(jset)

    return nondeg


def split_by_peaks(kbl: Iterable[AbstractPathway], abstract: bool=False,
                   collapse_directions: bool=True)\
    -> Dict[Tuple[str, str], List[AbstractPathway]]:
    r"""Collect pathways with the same 2D resonance.

    Parameters
    ----------
    kbl
        Iterable of pathways.
    abstract
        Splits by subbranches instead of individual peaks. Looks at peaks
        independent of :math:`J_i`.
    collapse_directions
        Assigns positive- and negative-frequency coherences to the same peak.

    Returns
    -------
    Dict[Tuple[str, str], List[AbstractPathway]]
        Dictionary from :attr:`Pathway.peak` or :attr:`Pathway.abstract_peak` to
        lists of pathways.
    """
    ret: Dict[Tuple[str, str], List[Pathway]] = {}
    attr: str = 'abstract_peak' if abstract else 'peak'
    for dl in kbl:
        peak = getattr(dl, attr)
        if collapse_directions and dl.leaf.is_SI():
            peak = peak_conj(peak)
        ret.setdefault(peak, []).append(dl)

    return ret


def pprint_dllist(dllist: Sequence[AbstractPathway],
                  abstract: bool=False,
                  angles: Optional[AnglesT]=None,
                  print: Callable=print):
    """Pretty print a list of :class:`Pathway` or :class:`DressedPathway`.

    Parameters
    ----------
    dllist
        Sequence of :class:`Pathway` or :class:`DressedPathway`, will print more
        information for the latter.
    abstract
        Use spectroscopic state labels relative to the ground state, instead of
        actual J values.
    angles
        Evaluate polarization tensor components and R-factor for these
        polarizations.
    print
        Use provided callable instead of :func:`print` built-in.
    """
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


def print_dl_dict(dldict: Mapping[Any, Sequence[Pathway]],
                  fields: Sequence[str]=None):
    """Pretty print a dict of pathways."""
    for k in dldict:
        print(k)
        for dl in dldict[k]:
            print("   ", dl.custom_str(fields=fields))


def print_dl_tuple_dict(dldict: Mapping[Any, Sequence[Tuple[Pathway, float]]],
                        fields: Optional[Sequence[str]]=None):
    """Pretty print a dict of (pathway, float) tuples."""
    for k in dldict:
        print(k)
        for dl, a in dldict[k]:
            print("   {:s}, {:f}".format(dl.custom_str(fields=fields), a))


# * Peaks without line shapes
@define
class Peak2D:
    """Represent 2D IR resonance produced by multiple pathways."""
    pump_wl: float
    "Pump wavenumber."
    probe_wl: float
    "Probe wavenumber."
    peak: Tuple[str, str]
    """Peak identifier.

    See also
    --------
    Pathway.peak, Pathway.abstract_peak
    """
    amplitude: complex
    """Sum of pathway amplitudes producing the peak.

    Rephasing pathway amplitudes have their sign flipped.

    Notes
    -----
    The Maxwell wave equation flips the sign of the negative frequency material
    polarization component, such that the positive and negative components add
    up to a real cosine wave. At the level of TD perturbation theory the signs
    are oppposite and summing corresponding rephasing and non-rephasing pathway
    amplitudes (or a pathway with its complex conjugate) gives zero. This means
    that if you want to estimate some signal amplitude/intensity based on
    amplitudes, you might get a zero unless you normalize the signs.

    See also
    --------
    Pathway.amplitude
    """
    intensity: complex
    """Pathway intensity."""
    max_intensity: complex
    """Pathway intensity at the peak of the line, assuming Lorentizan profile."""
    dp_list: Optional[List[DressedPathway]] = field(eq=False, default=None)
    """List of pathways backing this 2D resonance."""
    params: Optional[Dict[str, Any]] = None
    """`p`, pressure in atm; `tw`, waiting time in s; `angles`, polarization."""

    def max_abs_coeff(self, E12: float, conc: float) -> complex:
        """Return max abs. coeff from :attr:`max_intensity`.

        Assumes self-heterodyne probe detection.

        Parameters
        ----------
        E12 : float
            Product of field integrals of pump pulses.
        conc : float
            Concentration in 1/m**3.
        """
        return self.max_intensity/np.pi*conc*E12

    @classmethod
    def from_dp_list(cls, dp_list: List[DressedPathway],
                     tw: Optional[float]=0.0,
                     angles: Optional[AnglesT]=None,
                     p: Optional[float]=1.0) -> "Peak2D":
        r"""Make :class:`Peak2D` from a list of :class:`DressedPathway` associated
        with the same 2D resonance.

        Parameters
        ----------
        dp_list
            List of dressed pathway, can mix rephasing and nonrephasing pathways.
        tw
            Delay between second and third pulse (in s).
        angles
            Sequence of polarization angles.
        p
            Pressure (in atm), only matters for :attr:`max_intensity`.

        Returns
        -------
        Peak2D
        """
        peak = dp_list[0].peak
        pu, pr = u.nu2wn(dp_list[0].nu(0)), abs(u.nu2wn(dp_list[0].nu(2)))
        amplitude, intensity, max_intensity = 0.0, 0.0, 0.0
        for dp in dp_list:
            pre_amp = np.imag(dp.amplitude(tw=tw, angles=angles))
            amplitude += pre_amp*np.sign(dp.nu(2))
            intensity += pre_amp*np.pi*2*np.pi*dp.nu(2)/4/C.epsilon_0/C.c
            max_intensity += pre_amp*np.pi*2*np.pi*dp.nu(2)/4/C.epsilon_0/C.c\
                /dp.gamma(2)/p

        return cls(pu, pr, peak, amplitude, intensity, max_intensity, dp_list,
                   params=dict(p=p, tw=tw, angles=angles))


class Peak2DList(list):
    """List of 2D peaks with easy access to pump, probe frequencies, peak
    intensities and peak identifiers.

    Calling :meth:`copy` or slicing returns a regular list.
    """
    def __init__(self, iterable=()):
        super().__init__(iterable)
        self.normalized = False

    @classmethod
    def from_dp_list(cls, dpl: List[DressedPathway],
                     tw: Optional[float]=0.0,
                     angles: Optional[AnglesT]=None,
                     p: float=1.0) -> "Peak2DList":
        split_dls = split_by_peaks(dpl)
        pl = Peak2DList()
        for dp in split_dls.values():
            pl.append(Peak2D.from_dp_list(dp, tw, angles, p))
        pl.sort_by_amplitudes()

        return pl

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> "Peak2DList":
        """Read peak list from HDF5 file."""
        with h5py.File(path, mode='r') as f:
            pl = cls()
            params = json.loads(f.attrs['params']) # type: ignore
            for pu, pr, amp, inte, minte, peak in zip(
                    f['pumps'], f['probes'], f['amplitudes'], # type: ignore
                    f['intensities'], f['max_intensities'], f['peaks']): # type: ignore
                pl.append(
                    Peak2D(pu, pr, tuple(json.loads(peak)), amp, inte, minte,
                           params=params)) # type: ignore

        return pl

    @classmethod
    def from_params_dict(cls, params: Mapping) -> "Peak2DList":
        """Calculate list of peaks based on toml input data."""
        if params['spectrum']['type'] != 'peaks':
            raise ValueError("Wrong spectrum type requested.")

        dpws = DressedPathway.from_params_dict(params['pathways'])
        pl = Peak2DList.from_dp_list(
            dpws, tw=params['spectrum']['tw']*1e-12,
            angles=params['spectrum']['angles'])

        return pl

    @property
    def pumps(self) -> List[float]:
        """List of pump frequencies."""
        return [peak.pump_wl for peak in self]

    @property
    def probes(self) -> List[float]:
        """List of probe frequencies."""
        return [peak.probe_wl for peak in self]

    @property
    def peaks(self) -> List[Tuple[str, str]]:
        """Peak strings."""
        return [peak.peak for peak in self]

    @property
    def amplitudes(self) -> List[complex]:
        """Peak amplitude--sum of :meth:`DressedPathway.amplitude` over all pathways
        contributing to a 2D peak.
        """
        return [peak.amplitude for peak in self]

    @property
    def intensities(self) -> List[complex]:
        """Peak pathway intensities."""
        return [peak.intensity for peak in self]

    @property
    def max_intensities(self) -> List[complex]:
        """Max peak intensity assuming Lorentzian profile."""
        return [peak.max_intensity for peak in self]

    @staticmethod
    def _sort_func(peak):
        return abs(peak.amplitude)

    def sort_by_amplitudes(self):
        """Sort peaks by amplitude.

        Ensures that strong peaks are not covered by weak ones in scatter plot.
        """
        self.sort(key=self._sort_func)

    def get_by_peak(self, peak: Tuple[str, str]) -> Peak2D:
        """Return peak with :attr:`Peak2D.peak` equal to ``peak``."""
        for p in self:
            if p.peak == peak:
                return p
        else:
            raise IndexError("'{!s}' not found in the list".format(peak))

    def to_file(self, path: Union[str, Path], metadata: Optional[Dict]=None):
        """Save peak list to HDF5 file."""
        with h5py.File(path, mode='w') as f:
            f.create_dataset("pumps", data=self.pumps)
            f.create_dataset("probes", data=self.probes)
            f.create_dataset("amplitudes", data=self.amplitudes)
            f.create_dataset("intensities", data=self.intensities)
            f.create_dataset("max_intensities", data=self.max_intensities)
            f.create_dataset(
                "peaks", data=[json.dumps(peak) for peak in self.peaks])
            f.attrs['params'] = json.dumps(self[0].params)
            if metadata:
                f.attrs['metadata'] = json.dumps(metadata)


def equiv_peaks(pw: Union[Pathway, List[Tuple[RotState]]],
                pl: Peak2DList) -> Peak2DList:
    """Return peaks from ``pl`` based on pathways equivalent to ``pw``."""
    new_pl = Peak2DList()
    for peak in pl:
        if any(dp.leaf.is_equiv_pathway(pw) for dp in peak.dp_list):
            new_pl.append(peak)
    new_pl.sort_by_amplitudes()

    return new_pl
