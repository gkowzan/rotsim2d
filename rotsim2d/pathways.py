r"""Generate all Liouville pathways for nth order rovibrational excitation."""
# * Imports, constants and enums
import enum
from copy import deepcopy
import operator as op
from functools import reduce, lru_cache
from typing import List, Union, Tuple, Iterable, Sequence, Optional, Callable
import numpy as np

import anytree as at
from anytree.exporter import UniqueDotExporter

from molspecutils.molecule import DiatomState, SymTopState, RotState

#: Right-circular polarized light
right_pol = (5/4*np.pi, -np.pi/2)
#: Left-circular polarized light
left_pol = (np.pi/4, np.pi/2)

def nodeattrfunc(node):
    if isinstance(node, LightInteraction):
        return "shape=box,label=\"{:s}\"".format(node.fullname)
    else:
        return "shape=ellipse,label=\"{:s}\"".format(node.name)


class Side(enum.IntEnum):
    """Ket- or Bra-side excitation by :class:`LightInteraction`"""
    KET = 1
    BRA = -1

class KSign(enum.IntEnum):
    """Positive or negative part light spectrum."""
    POS = 1
    NEG = -1


ks = {
    (KSign.NEG, KSign.POS, KSign.POS): 'SI',
    (KSign.POS, KSign.NEG, KSign.NEG): 'SI',
    (KSign.POS, KSign.NEG, KSign.POS): 'SII',
    (KSign.NEG, KSign.POS, KSign.NEG): 'SII',
    (KSign.POS, KSign.POS, KSign.NEG): 'SIII',
    (KSign.NEG, KSign.NEG, KSign.POS): 'SIII'
}

# * LightInteraction
class LightInteraction(at.NodeMixin):
    """Represents (dipole) interaction between the system and a light beam.

    Parameters
    ----------
    name : str
        Identifier for the light beam.
    side : Side
        Ket or bra excitation.
    sign : KSign
        Wavevector sign of created polarization contribution.
    readout : bool
        Readout or actual light interaction.
    """
    def __init__(self, name: str, side: Side, sign: KSign, readout: bool=False,
                 angle: Union[float, Tuple[float]]=0.0, parent=None, children=None):
        super(LightInteraction, self).__init__()
        self.separator = "->"
        self.name = name
        self.side = side
        self.sign = sign
        if isinstance(angle, tuple):
            self.angle = (angle[0], angle[1]*sign)
        else:
            self.angle = angle
        self.readout = readout
        self.fullname = "{:s}(side={:d}, sign={:d})".format(self.name, self.side, self.sign)
        self.parent = parent
        if children:
            self.children = children

    @property
    def tier(self):
        return sum((1 for x in self.ancestors if isinstance(x, LightInteraction)))

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return "LightInteraction({:s}, side={:d}, sign={:d})".format(self.name, self.side, self.sign)


# * KetBra
class KetBra(at.NodeMixin):
    """Node in a tree of excited states.

    Each KetBra is described by :attr:`ket` and :attr:`bra` which describes the
    density matrix element.
    """
    def __init__(self, ket: RotState, bra: RotState, parent=None, children=None):
        super(KetBra, self).__init__()
        self.ket = ket
        self.bra = bra
        self.separator = "->"
        self.name = "|{:s}><{:s}|".format(self.ket.name, self.bra.name)
        self.parent = parent
        if children:
            self.children = children
        self.to_statelist = lru_cache(None)(self.to_statelist)
        self._pathway_info_cache = None

    def get(self, side: Side) -> RotState:
        """Index KetBra by Side enum."""
        if side == Side.KET:
            return self.ket
        elif side == Side.BRA:
            return self.bra

    def print_tree(self):
        """Pretty print excitation tree."""
        for pre, _, node in at.RenderTree(self):
            treestr = "{:s}{:s}".format(pre, node.name)
            print(treestr)

    def __repr__(self):
        return "Ketbra({!r}, {!r})".format(self.ket, self.bra)

    def __str__(self):
        return self.name

    def __eq__(self, o):
        if not isinstance(o, KetBra):
            return NotImplemented
        return self.ket == o.ket and self.bra == o.bra

    def _pathway_info(self):
        """Return transitions, coherences and 4-fold dipole for pathway.

        The key task of this function is to prepare arguments for evaluation of
        G-factors and polarization tensors. At the end, `wbras` contains j
        values (`wbras[i][0]`) and depths of LightInteraction objects
        (`wbras[i][1]). Depths are counted by ignoring KetBra nodes and give
        information which of the four LightInteraction's acting on the system in
        sequence produced the j-state on the left (if the whole sequence is read
        from right to left). `wbras` can then be used to construct the four-fold
        dipole interaction operator::

            <wbras[0][0]|wbras[0][1]|wbras[1][0]><wbras[1][0]|wbras[1][1]|wbras[2][0]>\
            <wbras[2][0]|wbras[2][1]|wbras[3][0]><wbras[3][0]|wbras[3][1]|wbras[0][0]>

        The edge j-states are the same by definition of the process. Having
        defined this four-fold operator, we can then use spherical tensor
        operator decomposition and calculate polarization tensors and G-factors.
        """
        if self._pathway_info_cache is None:
            kb_series = self.ketbras()
            coherences, transitions = [], []
            wkets, wbras = [], []
            for i in range(1, len(kb_series)):
                kb, kbp = kb_series[i], kb_series[i-1]
                if kb.parent.side is Side.KET:
                    pair = (kbp.ket, kb.ket)
                    wkets.insert(0, (kb.ket.j,  kb.parent.tier))
                else:
                    pair = (kbp.bra, kb.bra)
                    wbras.append((kb.parent.parent.bra.j, kb.parent.tier))
                transitions.append(pair)

                if kb.parent.readout:
                    break
                coherences.append((kb.bra, kb.ket))

            wbras.extend(wkets)
            self._pathway_info_cache = (coherences, transitions, wbras)

        return self._pathway_info_cache

    def transitions(self):
        """Return list of transitions as a list of state pair."""
        return self._pathway_info()[1]

    def color_tier(self):
        """Number of colors neede to create this pathway."""
        if self.parent.readout:
            return len(set((tuple(sorted(x)) for x in self.transitions()[:-1])))
        return len(set((tuple(sorted(x)) for x in self.transitions())))

    def to_statelist(self, diatom=False, normalize=False) -> List[Tuple[RotState]]:
        """KetBras leading to this one as a list of state pairs.

        Drops `k` quantum number if `diatom` is True. Subtracts the initial `j`
        value from all states if `normalize` is True.
        """
        statelist = [(kb.ket, kb.bra) for kb in self.ketbras()]
        if diatom:
            statelist = [(DiatomState.from_symtop(state1), DiatomState.from_symtop(state2))
                         for state1, state2 in statelist]
        if normalize:
            startj = statelist[0][0].j
            statelist = [(state1._replace(j=state1.j-startj), state2._replace(j=state2.j-startj))
                         for state1, state2 in statelist]

        return statelist

    def copy(self):
        return deepcopy(self)

    def savepng(self, path):
        """Save excitation tree as an image."""
        UniqueDotExporter(self, nodeattrfunc=nodeattrfunc).to_picture(str(path))

        return path

    def conj(self):
        """Return conjugate of this KetBra."""
        return KetBra(self.bra, self.ket)

    def normalized(self):
        """Return copy of self with ket being the lower nu, j level."""
        if self.ket.nu > self.bra.nu or (self.ket.nu == self.bra.nu and self.ket.j > self.bra.j):
            return self.conj()
        return self

    def kb_ancestor(self, ancestor=None):
        """Return first KetBra ancestor."""
        if ancestor is None:
            ancestor = self.ancestors[0]
        if isinstance(ancestor, KetBra):
            return ancestor
        else:
            return self.kb_ancestor(ancestor.ancestors[0])

    def diagonals(self, sort=False):
        """Collect diagonal leaves."""
        pops = [x for x in self.leaves if x.ket == x.bra]
        if sort:
            pops.sort(key=lambda x: x.ket.j)
            pops.sort(key=lambda x: x.ket.nu)

        return pops

    def is_diagonal(self):
        return self.ket == self.bra

    def is_rephasing(self):
        ints = self.interactions()

        return ints[0].sign != ints[2].sign

    def is_SI(self, order=None):
        r"""Check if :math:`\vec{k}_s = -\vec{k}_1+\vec{k}_2+\vec{k}_3` (rephasing).

        `order` is a list of :class:`LightInteraction` names specifying which
        interaction corresponds to which wavevector.  First name is
        :math:`\vec{k}_1`, second one is :math:`\vec{k}_2`, etc.
        """
        if order is None:
            order = ('omg1', 'omg2', 'omg3')
        return ks.get(tuple(self.interaction(name).sign for name in order)) == 'SI'

    def is_SII(self, order=None):
        r"""Check if :math:`\vec{k}_s=\vec{k}_1-\vec{k}_2+\vec{k}_3` (non-rephasing).

        See also
        --------
        is_SI
        """
        if order is None:
            order = ('omg1', 'omg2', 'omg3')
        return ks.get(tuple(self.interaction(name).sign for name in order)) == 'SII'

    def is_SIII(self, order=None):
        r"""Check if :math:`\vec{k}_s=\vec{k}_1+\vec{k}_2-\vec{k}_3` (double quantum).

        See also
        --------
        is_SI
        """
        if order is None:
            order = ('omg1', 'omg2', 'omg3')
        return ks.get(tuple(self.interaction(name).sign for name in order)) == 'SIII'

    def is_Pinitial(self):
        r"""Check if initial excitation is P-branch."""
        sl = self.to_statelist()
        return any(x.j==sl[0][0].j-1 for x in sl[1])

    def is_Rinitial(self):
        r"""Check if initial excitation is R-branch."""
        sl = self.to_statelist()
        return any(x.j==sl[0][0].j+1 for x in sl[1])

    def is_Qinitial(self):
        r"""Check if initial excitation is Q-branch."""
        sl = self.to_statelist()
        return all(x.j==sl[0][0].j for x in sl[1])

    def is_esa(self):
        """Check if pathway corresponds to excited-state absorption.

        Applies to two-color and three-color excitations and pathways not
        starting from the ground state.
        """
        sl = self.to_statelist()
        coh3 = sl[3]
        rnu = sl[0][0].nu
        return any((x.nu == rnu+1 for x in coh3)) and any((x.nu == rnu+2 for x in coh3))

    def is_sep(self):
        """Check if pathway corresponds to stimulated emission pumping.

        Second coherence is either excited population of rotational coherence
        state. Third coherence is between first excited state and ground state.
        """
        sl = self.to_statelist()
        coh2, coh3 = sl[2], sl[3]
        rnu = sl[0][0].nu
        return coh2[0].nu == rnu+1 and coh2[1].nu == rnu+1\
            and any(x.nu == rnu+1 for x in coh3)\
            and any(x.nu == rnu for x in coh3)

    def is_gshb(self):
        """Check if pathway corresponds to ground-state hole-burning.

        Second coherence is either ground population of rotational coherence
        state. Third coherence is between first excited state and ground state.
        """
        sl = self.to_statelist()
        coh2, coh3 = sl[2], sl[3]
        rnu = sl[0][0].nu
        return coh2[0].nu == rnu and coh2[1].nu == rnu\
            and any(x.nu == rnu+1 for x in coh3)\
            and any(x.nu == rnu for x in coh3)

    def is_doublequantum(self):
        sl = self.to_statelist()
        return any(x.nu==sl[0][0].nu for x in sl[2])

    def is_interstate(self):
        """Check for coherent state after second interaction."""
        return not self.ketbras()[2].is_diagonal()

    def is_overtone(self):
        return abs(self.ket.nu-self.bra.nu)>1

    def has_overtone(self):
        return len([kb for kb in self.ketbras() if kb.is_overtone()]) > 0

    def is_between(self, pump, probe):
        """Check if this pathway produces cross-peak between `kb1` and `kb2`.

        Called with self being a leaf.
        """
        pump_kb, probe_kb = self.ketbras()[1], self.ketbras()[3] 
        return (pump == pump_kb or pump == pump_kb.conj()) and (probe == self or probe == self.conj())

    def is_dfwm(self):
        """Check if this pathway contains only coherences corresponding to a single dipole transition."""
        return self.color_tier() == 1

    is_onecolor = is_dfwm

    def is_twocolor(self):
        """Check if pathway is two-color."""
        return self.color_tier() == 2

    def is_threecolor(self):
        """Check if pathway is three-color."""
        return self.color_tier() == 3

    def is_equiv_pathway(self, o):
        """Check if other pathway is R-factor-equivalent.

        `o` should either be :class:`KetBra` or a result of :meth:`to_statelist`
        with `normalize` set to True.

        This won't work if initial state is vibrationally excited.
        """
        s = self.to_statelist(diatom=True, normalize=True)
        if isinstance(o, KetBra):
            o = o.to_statelist(diatom=True, normalize=True)
        for sp, op in zip(s, o):
            if not (sp[0].j == op[0].j and sp[1].j == op[1].j and
                    sp[0].nu == op[0].nu and sp[1].nu == op[1].nu):
                return False
        return True

    def is_pathway(self, *kbs):
        """Match self to pathway consisting of `kbs`."""
        return tuple(self.ketbras()) == kbs

    def is_some_pathway(self, kbs):
        """Match self to one of pathways in `kbs`."""
        return any(self.is_pathway(*kb) for kb in kbs)

    def ketbras(self):
        """Ancestor KetBras and self."""
        return [x for x in self.ancestors if isinstance(x, KetBra)] + [self]

    def interactions(self) -> List[LightInteraction]:
        """Interactions which generated this KetBra."""
        return [x for x in self.ancestors if isinstance(x, LightInteraction)]

    def interaction(self, name: str) -> Optional[LightInteraction]:
        """Return :class:`LightInteraction` with given name."""
        for x in self.ancestors:
            if isinstance(x, LightInteraction) and x.name == name:
                return x

    def ksigns(self) -> Tuple[KSign, ...]:
        """Return signs of wavevectors of interactions."""
        return tuple(i.sign for i in self.interactions())

    def sides(self) -> Tuple[Side, ...]:
        """Return sides of DM on which interactions acted."""
        return tuple(i.side for i in self.interactions())

    def total_ksign(self) -> int:
        """Cumulative sign of the term."""
        return reduce(op.mul, self.ksigns(), 1)

    def total_side(self) -> int:
        """Cumulative side of the term."""
        return reduce(op.mul, self.sides(), 1)

# * Tree filtering functions
def make_remove(func: Callable) -> Callable:
    """Tree filtering function factory.

    Parameters
    ----------
    func
        Callable taking :class:`KetBra` instance and returning True of False.

    Returns
    -------
    Function taking a `KetBra` excitation tree and removing all branches for
    which `func` returns True.
    """
    def remove_func(ketbra):
        maxdepth = max(leaf.depth for leaf in ketbra.leaves)
        for l in ketbra.leaves:
            if func(l):
                l.parent = None

        return prune(ketbra, depth=maxdepth)
    return remove_func


def make_only(func: Callable) -> Callable:
    """Tree filtering function factory.

    Parameters
    ----------
    func
        Callable taking :class:`KetBra` instance and returning True or False.

    Returns
    -------
    Function taking a `KetBra` excitation tree and leaving only branches for
    which `func` returns True.
    """
    def only_func(ketbra):
        maxdepth = max(leaf.depth for leaf in ketbra.leaves)
        for l in ketbra.leaves:
            if not func(l):
                l.parent = None

        return prune(ketbra, depth=maxdepth)
    return only_func


only_SI = make_only(lambda kb: kb.is_SI())
only_SII = make_only(lambda kb: kb.is_SII())
only_SIII = make_only(lambda kb: kb.is_SIII())
only_nonrephasing = make_remove(lambda kb: kb.is_rephasing())
only_rephasing = make_only(lambda kb: kb.is_rephasing())
remove_nondiagonal = make_only(lambda kb: kb.is_diagonal())
remove_overtones = make_remove(lambda kb: kb.has_overtone())
only_esa = make_only(lambda kb: kb.is_esa())
only_sep = make_only(lambda kb: kb.is_sep())
only_gshb = make_only(lambda kb: kb.is_gshb())
only_dfwm = make_only(lambda kb: kb.is_dfwm())
only_twocolor = make_only(lambda kb: kb.is_twocolor())
only_threecolor = make_only(lambda kb: kb.is_threecolor())
remove_threecolor = make_remove(lambda kb: kb.is_threecolor())
remove_interstates = make_remove(lambda kb: kb.is_interstate())
only_interstates = make_only(lambda kb: kb.is_interstate())
only_ketside = make_only(lambda kb: all(s == Side.KET for s in kb.sides()[:3]))
remove_ketside = make_remove(lambda kb: all(s == Side.KET for s in kb.sides()[:3]))
only_Pinitial = make_only(lambda kb: kb.is_Pinitial())
only_Rinitial = make_only(lambda kb: kb.is_Rinitial())
only_Qinitial = make_only(lambda kb: kb.is_Qinitial())


def only_between(ketbra: KetBra, pump: KetBra, probe: KetBra) -> KetBra:
    """Limit tree to pathways bewteen `kb1` and `kb2`."""
    maxdepth = max(leaf.depth for leaf in ketbra.leaves)
    for l in ketbra.leaves:
        if not l.is_between(pump, probe):
            l.parent = None

    return prune(ketbra, depth=maxdepth)


def only_pathway(ketbra: KetBra, pathway: KetBra) -> KetBra:
    maxdepth = max(leaf.depth for leaf in ketbra.leaves)
    for l in ketbra.leaves:
        if not l.is_pathway(pathway):
            l.parent = None

    return prune(ketbra, depth=maxdepth)


def only_some_pathway(ketbra: KetBra, pathways: List[KetBra]) -> KetBra:
    maxdepth = max(leaf.depth for leaf in ketbra.leaves)
    for l in ketbra.leaves:
        if not l.is_some_pathway(pathways):
            l.parent = None

    return prune(ketbra, depth=maxdepth)


def prune(ketbra: KetBra, depth: int) -> KetBra:
    """Remove leaves whose depth is less than `depth`."""
    found = True
    while found:
        found = False
        for leaf in ketbra.leaves:
            if leaf.is_root:
                return leaf
            if leaf.depth < depth:
                found = True
                leaf.parent = None

    return ketbra


# * Tree-modifying functions
def readout(ketbra: KetBra) -> KetBra:
    """Generate populations from excitations."""
    for kb in ketbra.leaves:
        excite(kb, 'mu', 'ket', True)
    maxdepth = max(leaf.depth for leaf in ketbra.leaves)

    return prune(remove_nondiagonal(ketbra), depth=maxdepth)


def excited_states_symtop(state: SymTopState, dnu: int) -> List[int]:
    """Return states reachable from `state` by dipole interaction."""
    djs = (-1, 0, 1)
    states = []

    if state.nu+dnu>=0:
        for dj in djs:
            if (state.k==0 or state.j==0) and dj==0:
                continue
            if state.j+dj < state.k:
                continue
            if state.j+dj>=0:
                states.append(SymTopState(state.nu+dnu, state.j+dj, state.k))

    return states


def excited_states_diatom(state: DiatomState, dnu: int) -> List[int]:
    """Return states reachable from `state` by dipole interaction."""
    djs = (-1, 1)
    states = []

    if state.nu+dnu >= 0:
        for dj in djs:
            if state.j+dj >=0:
                states.append(DiatomState(state.nu+dnu, state.j+dj))

    return states


#: Poor man's polymorphism
excited_states = {
    DiatomState: excited_states_diatom,
    SymTopState: excited_states_symtop
}

ksigns = {
    (Side.BRA, 1): KSign.NEG,
    (Side.BRA, -1): KSign.POS,
    (Side.KET, 1): KSign.POS,
    (Side.KET, -1): KSign.NEG
}

def _excite(ketbra: KetBra, light_name: str, side: Side=Side.KET,
            readout: bool=False) -> KetBra:
    for dnu in (-1, 1):
        states = excited_states[type(ketbra.get(side))](ketbra.get(side), dnu)
        if states:
            li = LightInteraction(light_name, side, ksigns[(side, dnu)],
                                  readout, parent=ketbra)
            for state in states:
                if side==Side.KET:
                    KetBra(state, ketbra.bra, parent=li)
                elif side==Side.BRA:
                    KetBra(ketbra.ket, state, parent=li)

    return ketbra


def excite(ketbra: KetBra, light_name: str, part: str='ket',
           readout: bool=False) -> KetBra:
    """Generate all excitations of `ketbra`.

    Modifies `ketbra` in place. Only parallel transitions.

    Parameters
    ----------
    ketbra : KetBra
        State to excite.
    light_name : str
        Identifier for the EM field doing the excitation.
    part : str
        'ket', 'bra' or 'both', consider ket, bra or double-sided excitations.
    readout : bool, optional
        Readout or actual light interaction.

    Returns
    -------
    KetBra
    """
    # ket excitation
    if part not in ('ket', 'bra', 'both'):
        raise ValueError("`part` has to be either 'ket', 'bra' or 'both'")
    if part in ('ket', 'both'):
        ketbra = _excite(ketbra, light_name, Side.KET, readout)
    if part in ('bra', 'both'):
        ketbra = _excite(ketbra, light_name, Side.BRA, readout)

    return ketbra


def multi_excite(ketbra: KetBra, light_names: List[str], parts: Optional[List]=None) -> KetBra:
    """Generate multiple excitations of `ketbra`.

    Parameters
    ----------
    ketbra : KetBra
        State to excite.
    light_names : list of str
        Names of EM fields, length sets the number of excitations.
    parts : None or list
        None for first ket excitation and rest 'both' excitations.
    """
    if parts is not None and len(parts) != len(light_names):
        raise ValueError("len(parts) != len(light_names)")
    if parts is None:
        parts = ['ket']
        parts.extend(['both']*(len(light_names)-1))

    if light_names:
        excite(ketbra, light_names.pop(0), parts.pop(0), False)
        for kb in ketbra.leaves:
            multi_excite(kb, light_names[:], parts[:])

    return ketbra.root


def gen_roots(jiter: Iterable, rotor: str='linear', kiter_func: Callable=None) -> List[KetBra]:
    roots = []
    for j in jiter:
        if rotor == 'linear':
            roots.append(KetBra(DiatomState(0, j), DiatomState(0, j)))
        elif rotor == 'symmetric':
            if kiter_func is None:
                kiter = range(0, j+1)
            else:
                kiter = kiter_func(j)
            for k in kiter:
                roots.append(KetBra(SymTopState(0, j, k), SymTopState(0, j, k)))
    return roots


def gen_excitations(root, light_names, parts, meths=None) -> KetBra:
    root = multi_excite(root, light_names, parts=parts)
    if meths is not None:
        for meth in meths:
            root = meth(root)
    root = readout(root)

    return root


def gen_pathways(jiter: Iterable, meths: Optional[Sequence[Callable]]=None,
                 rotor: str='linear', kiter_func: Callable=None,
                 pump_overlap: bool=False) -> List[KetBra]:
    roots = gen_roots(jiter, rotor, kiter_func)
    pws = [gen_excitations(root, ['omg1', 'omg2', 'omg3'],
                           ['ket', 'both', 'both'], meths)
           for root in roots]
    if pump_overlap:
        roots = gen_roots(jiter, rotor, kiter_func)
        pws.extend(
            [gen_excitations(root, ['omg2', 'omg1', 'omg3'],
                             ['ket', 'both', 'both'],
                             meths)
             for root in roots])
    pws = [p for p in pws if len(p.children)]

    return pws


def geometric_factor(leaf: KetBra):
    """Return polarization and angular momentum sequence for a pathway.

    Same as in leaf_response."""
    kbs = leaf.ketbras()
    wkets, wbras = [], []
    for i in range(1, len(kbs)):
        kb, kbp = kbs[i], kbs[i-1]
        if kb.parent.side is Side.KET:
            wkets.insert(0, (kb.kj, kb.parent.tier))
        else:
            wbras.append((kb.parent.parent.bj, kb.parent.tier))
        if kb.parent.readout:
            wbras.extend(wkets)
            break
    return tuple([x[0] for x in wbras]), tuple([x[1] for x in wbras])


def tensor_analysis(kb: KetBra) -> List[Tuple]:
    """Assign cross peak and geometric factor to each pathway."""
    ret = []
    for l in kb.leaves:
        js, pols = geometric_factor(l)
        kbs = l.ketbras()
        peak = (kbs[1].normalized().name, kbs[3].normalized().name)
        ret.append((peak, js, pols))

    return ret
