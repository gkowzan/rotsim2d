r"""Associate polarizations, angular momentum factors and pathways amplitudes
with pathways, i.e. everything but line shapes.

TODO:
- change sys_params and pops
"""
from __future__ import annotations
from typing import Optional, Iterable, Tuple, List, Sequence
from collections import namedtuple
from math import isclose
import numpy as np
import knickknacks.units as u
import molspecutils.molecule as mol
import scipy.constants as C
from rotsim2d.pathways import KetBra, Side, KSign
from rotsim2d.couple import four_couple
import rotsim2d.couple as cp

#: Spectroscopic notation for transitions
dj_to_letter = {-2: "O", -1: "P", 0: "Q", 1: "R", 2: "S"}
def abstract_format(dnu: int, dj: int):
    if dnu==0 and dj==0:
        return "0"
    return str(dnu)+dj_to_letter[dj]


def abstract_state_label(state: mol.RotState, ref_state: mol.RotState) -> str:
    return abstract_format(state.nu-ref_state.nu, state.j-ref_state.j)


def abstract_pair_label(pair: Tuple[mol.RotState], ref_state: mol.RotState) -> str:
    return "|{:s}><{:s}|".format(
        abstract_state_label(pair[0], ref_state),
        abstract_state_label(pair[1], ref_state)
    )

class Pathway:
    """Collect information on a pathway based on KetBra tree leaf without
    specializing it to any specific vibrational mode of a molecule."""
    fields = ['leaf', 'coherences', 'transitions', 'js', 'angles', 'const',
              'geo_label', 'trans_label', 'tw_coherence', 'peak', 'abstract_peak',
              'experimental_label', 'colors']

    def __init__(self, leaf: KetBra):
        self.leaf = leaf
        self.isotropy = 1/np.sqrt(2*leaf.root.ket.j+1)
        self.coherences, self.transitions, wbras = leaf._pathway_info()
        self.js = tuple(x[0] for x in wbras)
        self.angles = tuple(x[1] for x in wbras)
        self.const = (1.0j/C.hbar)**(len(self.transitions)-1)*self.leaf.total_side()
        self.tw_coherence = not KetBra(*self.coherences[1]).is_diagonal()
        self.peak = ("|{:s}><{:s}|".format(self.coherences[0][0].name, self.coherences[0][1].name),
                     "|{:s}><{:s}|".format(self.coherences[2][0].name, self.coherences[2][1].name))
        self.abstract_peak = (abstract_pair_label(self.coherences[0], self.leaf.root.ket),
                              abstract_pair_label(self.coherences[2], self.leaf.root.ket))

    def __eq__(self, o):
        if not isinstance(o, Pathway):
            return NotImplemented
        return tuple(self.leaf.to_statelist()) == tuple(o.leaf.to_statelist())

    def __hash__(self):
        return hash(tuple(self.leaf.to_statelist()))

    @property
    def geo_label(self):
        """G-factor label."""
        return geometric_labels[tuple(j-self.js[0] for j in self.js)]

    @property
    def trans_label(self):
        """Three-fold transition label (Murdock style)."""
        return ''.join((self._trans_label(i) for i in (0, 1, 2)))

    @property
    def experimental_label(self):
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
        return self.leaf.color_tier()

    def _trans_label(self, i):
        labels = {-1: 'P', 0: 'Q', 1: 'R'}
        trans = sorted(self.transitions[i], key=lambda x: x.nu)

        return labels[trans[1].j-trans[0].j]

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
        ret = four_couple(self.js, self.angles)
        if relative:
            ret /= four_couple(self.js, [0.0]*4)

        return ret

    def gfactors(self):
        js = list(self.js)
        return (cp.G(*(js + [0])), cp.G(*(js + [1])), cp.G(*(js + [2])))

    def T00s(self):
        angles = list(self.angles)
        return (cp.T00(*(angles + [0])), cp.T00(*(angles + [1])), cp.T00(*(angles + [2])))

    @classmethod
    def from_kb_tree(cls, kb_tree: KetBra) -> List[Pathway]:
        """Make a list of Pathway's from KetBra tree."""
        return [cls(leaf) for leaf in kb_tree.root.leaves]

    @classmethod
    def from_kb_list(cls, kb_list: KetBra) -> List[Pathway]:
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

    def print_diagram(self, abstract=False):
        for kb in self.leaf.ketbras()[::-1]:
            l, r = Pathway._ketbra_symbols(kb)
            if abstract:
                print("{l:s}|{ket:2s}><{bra:>2s}|{r:s}".format(
                    l=l, r=r, ket=abstract_state_label(kb.ket, kb.root.ket),
                    bra=abstract_state_label(kb.bra, kb.root.ket)))
            else:
                print(f"{l}|{kb.ket.name}><{kb.bra.name}|{r}")

    def _tw_pprint(self):
        print(f"Coherence during waiting time: {self.tw_coherence!r}")

    def pprint(self, abstract=False):
        print("diagram:")
        self.print_diagram(abstract=abstract)
        print(f"G-factor label: {self.geo_label}")
        print(f"Transition chain label: {self.trans_label}")
        print(f"Experimental label: {self.experimental_label}")
        print("Colors: {:d}".format(self.leaf.color_tier()))
        pol = self.angles[0]
        if not all((isclose(pol, x) for x in self.angles)):
            print("Intensity relative to XXXX polarization: {}".format(self.geometric_factor(True)))
        self._tw_pprint()

    def __repr__(self):
        return f"Pathway(leaf={self.leaf!r})"


class DressedPathway(Pathway):
    """Excitation pathway specialized to a vibrational mode."""
    def __init__(self, leaf: KetBra, vib_mode: mol.VibrationalMode, T: float):
        Pathway.__init__(self, leaf)
        self.vib_mode = vib_mode
        self.T = T
        for pair in self.transitions:
            self.const *= vib_mode.mu(pair)
        self.const *= vib_mode.equilibrium_pop(self.leaf.root.ket, T)

    def __eq__(self, o):
        if not isinstance(o, DressedPathway):
            return NotImplemented
        return Pathway.__eq__(self, o) and isclose(self.T, o.T) and self.vib_mode == o.vib_mode
        
    def nu(self, i: int) -> float:
        return self.vib_mode.nu(self.coherences[i])

    def gamma(self, i: int) -> float:
        return self.vib_mode.gamma(self.coherences[i])

    def intensity(self, tw: Optional[float]=None) -> float:
        ret = self.isotropy*self.const*self.geometric_factor()
        if tw is not None:
            ret *= np.exp(-2.0*np.pi*tw*(1.0j*self.nu(1) + self.gamma(1)))

        return ret

    def _tw_pprint(self):
        print(f"Coherence during waiting time: {self.tw_coherence!r}", end='')
        if self.tw_coherence:
            print(", {:.2f} cm-1".format(u.nu2wn(self.nu(1))))
        else:
            print()


    @classmethod
    def from_kb_tree(cls, kb_tree: KetBra, vib_mode: mol.VibrationalMode,
                     T: float) -> List[DressedPathway]:
        """Make a list of DressedPathway's from KetBra tree."""
        return [cls(leaf, vib_mode, T) for leaf in kb_tree.root.leaves]

    @classmethod
    def from_kb_list(cls, kb_list: KetBra, vib_mode: mol.VibrationalMode,
                     T: float) -> List[DressedPathway]:
        """Make a list of DressedPathway's from KetBra list."""
        return sum((cls.from_kb_tree(kb_tree, vib_mode, T) for kb_tree in kb_list), [])

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
        if tw is not None:
            ret *= np.exp(-2.0*np.pi*tw*(1.0j*self.nu(1) + self.gamma(1)))

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


def split_by_peaks(kbl: Iterable[DressedLeaf]):
    ret = {}
    for dl in kbl:
        ret.setdefault(dl.peak, []).append(dl)
    return ret


def pprint_dllist(dllist, abstract=False):
    for i, dl in enumerate(dllist):
        if i == 0:
            print('-'*10)
            if isinstance(dl, DressedPathway):
                print('pump = {:.2f} cm-1, probe = {:.2f} cm-1'.format(
                    u.nu2wn(dl.nu(0)), u.nu2wn(dl.nu(2))))
            dl.pprint(abstract=abstract)
        else:
            print()
            dl.pprint(abstract=abstract)
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
    intensities."""
    @property
    def pumps(self):
        return [peak.pump_wl for peak in self]

    @property
    def probes(self):
        return [peak.probe_wl for peak in self]

    @property
    def sigs(self):
        return [peak.sig for peak in self]

    @staticmethod
    def _sort_func(peak):
        return abs(peak.sig)

    def sort_by_sigs(self):
        self.sort(key=self._sort_func)


def peak_list(ll: List[DressedPathway], tw: Optional[float]=0.0, return_dls: bool=False) -> Peak2DList:
    """Create a list of 2D peaks from a list of DressedLeaf's.

    Optionally return sorted list of DressedLeaf's corresponding to peaks.
    """
    ll = split_by_peaks(ll)
    pl = Peak2DList()
    dls = []
    for peak, dll in ll.items():
        pu, pr = u.nu2wn(dll[0].nu(0)), u.nu2wn(dll[0].nu(2))
        sig = 0.0
        for dl in dll:
            # with wigxjpf(300, 6):
            sig += np.imag(dl.intensity(tw=tw))
        pl.append(Peak2D(pu, pr, sig, peak))
        if return_dls:
            dls.append(dll)
    if return_dls:
        pairs = sorted(zip(pl, dls), key=lambda x: abs(x[0].sig))
        return Peak2DList([x[0] for x in pairs]), [x[1] for x in pairs]

    return pl


def equiv_peaks(pw: Pathway, pl: Peak2DList, dll: Sequence[Pathway]) -> Peak2DList:
    """Return peaks from `pl` which are polarization-equivalent to `pw`."""
    new_pl = Peak2DList()
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
        Map between detectio angle and a list of pathways.
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
