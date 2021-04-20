r"""Generate all Liouville pathways for nth order rovibrational excitation.
"""
# * Imports, constants and enums
import enum
from copy import deepcopy
import operator as op
from functools import reduce
from typing import List, Union, Tuple, Dict, Iterable, Sequence, Optional, Callable
import numpy as np

import anytree as at
from anytree.exporter import UniqueDotExporter

from spectroscopy.molecule import DiatomState, SymTopState, RotState 

dnus = [-1, +1] #: possible changes of vibrational state

right_pol = (5/4*np.pi, -np.pi/2)
left_pol = (np.pi/4, np.pi/2)

def nodeattrfunc(node):
    if isinstance(node, LightInteraction):
        return "shape=box,label=\"{:s}\"".format(node.fullname)
    else:
        return "shape=ellipse,label=\"{:s}\"".format(node.name)


class Side(enum.IntEnum):
    KET = 1
    BRA = -1

class KSign(enum.IntEnum):
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

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return "LightInteraction({:s}, side={:d}, sign={:d})".format(self.name, self.side, self.sign)


# * KetBra
class KetBra(at.NodeMixin):
    def __init__(self, ket, bra, pop: float=1.0, parent=None, children=None):
        super(KetBra, self).__init__()
        self.ket = ket
        self.bra = bra
        self.separator = "->"
        self.name = "|{:s}><{:s}|".format(self.ket.name, self.bra.name)
        self.parent = parent
        if children:
            self.children = children

        self.pop = pop                # fractional occupancy

    def print_tree(self):
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
        else:
            return self.ket == o.ket and self.bra == o.bra

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
        UniqueDotExporter(self, nodeattrfunc=nodeattrfunc).to_picture(str(path))

        return path

    def conj(self):
        """Return conjugate of this KetBra."""
        return KetBra(self.bra, self.ket)

    def normalized(self):
        """Return copy with ket being the lower nu, j level."""
        if self.ket.nu > self.bra.nu or (self.ket.nu == self.bra.nu and self.ket.j > self.bra.j):
            return self.conj()
        else:
            return self

    def kb_ancestor(self, ancestor=None):
        """Return first KetBra ancestor."""
        if ancestor is None:
            ancestor = self.ancestors[0]
        if isinstance(ancestor, KetBra):
            return ancestor
        else:
            return self.kb_ancestor(ancestor.ancestors[0])

    def diagonals(self, sort=False, rot=True):
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

    def is_esa(self):
        return (self.ket.nu != self.root.ket.nu) and (self.bra.nu != self.root.ket.nu)

    def is_sep(self):
        return all((x.ket.nu-x.root.ket.nu)<2 and (x.bra.nu-x.root.ket.nu)<2 for x in self.ketbras())\
            and not self.is_dfwm() and self.is_twocolor()

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
        kbs = self.ketbras()
        return all((kb.is_diagonal() or kb == kbs[1] or kb == kbs[1].conj() for kb in kbs))

    def is_twocolor(self, same=None):
        """Check if pathway is at most two-color.

        By default interactions 'omg1' and 'omg2' are checked. If another pair
        should be checked then the names should be listed as the argument.
        """
        chain = self.ancestors + (self,)
        # find indices of interactions 
        if same is None:
            same = ('omg1', 'omg2')
        i1 = [x[0] for x in enumerate(chain) if isinstance(x[1], LightInteraction) and x[1].name == same[0]][0]
        i2 = [x[0] for x in enumerate(chain) if isinstance(x[1], LightInteraction) and x[1].name == same[1]][0]

        # get transitions
        if chain[i1].side == Side.BRA:
            trans1 = (chain[i1-1].bra, chain[i1+1].bra)
        else:
            trans1 = (chain[i1-1].ket, chain[i1+1].ket)
        if chain[i2].side == Side.BRA:
            trans2 = (chain[i2-1].bra, chain[i2+1].bra)
        else:
            trans2 = (chain[i2-1].ket, chain[i2+1].ket)

        # check if transitions are between the same states
        if trans1 == trans2 or trans1 == trans2[::-1]:
            return True
        return False

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

    def interaction(self, name: str) -> LightInteraction:
        """Return :class:`LightInteraction` with given name."""
        for x in self.ancestors:
            if isinstance(x, LightInteraction) and x.name == name:
                return x

    def ksigns(self) -> Tuple[KSign]:
        """Return signs of wavevectors of interactions."""
        return tuple(i.sign for i in self.interactions())

    def sides(self) -> Tuple[Side]:
        """Return sides of DM on which interactions acted."""
        return tuple(i.side for i in self.interactions())

    def total_ksign(self) -> int:
        """Cumulative sign of the term."""
        return reduce(op.mul, self.ksigns(), 1)

    def total_side(self) -> int:
        """Cumulative side of the term."""
        return reduce(op.mul, self.sides(), 1)

# * Tree filtering functions
def remove_nonrephasing(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_rephasing():
            l.parent = None

    return prune(ketbra)


def remove_rephasing(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if l.is_rephasing():
            l.parent = None

    return prune(ketbra)


def remove_nondiagonal(ketbra: KetBra) -> KetBra:
    """Remove branches ending with non-diagonal leaves."""
    for l in ketbra.leaves:
        if not l.is_diagonal():
            l.parent = None

    return prune(ketbra)


def remove_esa(ketbra: KetBra) -> KetBra:
    """Remove excited state absorption."""
    for l in ketbra.leaves:
        if l.is_esa():
            l.parent = None

    return prune(ketbra)


def remove_overtones(ketbra: KetBra) -> KetBra:
    """Remove leaves with pathways containing overtone coherences."""
    for l in ketbra.leaves:
        if l.has_overtone():
            l.parent = None

    return prune(ketbra)


def only_between(ketbra: KetBra, pump: KetBra, probe: KetBra) -> KetBra:
    """Limit tree to pathways bewteen `kb1` and `kb2`."""
    for l in ketbra.leaves:
        if not l.is_between(pump, probe):
            l.parent = None

    return prune(ketbra)


def only_dfwm(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_dfwm():
            l.parent = None

    return prune(ketbra)


def only_sep(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_sep():
            l.parent = None

    return prune(ketbra)


def only_pathway(ketbra: KetBra, pathway: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_pathway(pathway):
            l.parent = None

    return prune(ketbra)


def only_SII(ketbra: KetBra) -> KetBra:
    """SII are non-rephasing."""
    for l in ketbra.leaves:
        if not l.is_SII():
            l.parent = None

    return prune(ketbra)


def only_SI(ketbra: KetBra) -> KetBra:
    """SI are rephasing without overtone coherences."""
    for l in ketbra.leaves:
        if not l.is_SI():
            l.parent = None

    return prune(ketbra)


def only_SIII(ketbra: KetBra) -> KetBra:
    """SIII are rephasing with overtone coherences."""
    for l in ketbra.leaves:
        if not l.is_SIII():
            l.parent = None

    return prune(ketbra)


def only_some_pathway(ketbra: KetBra, pathways: List[KetBra]) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_some_pathway(pathways):
            l.parent = None

    return prune(ketbra)


def remove_interstates(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.ketbras()[2].is_diagonal():
            l.parent = None

    return prune(ketbra)


def only_interstates(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if l.ketbras()[2].is_diagonal():
            l.parent = None

    return prune(ketbra)


def only_twocolor(ketbra: KetBra, same=None) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_twocolor(same):
            l.parent = None

    return prune(ketbra)


def prune(ketbra: KetBra) -> KetBra:
    """Remove leaves whose depth is less than maximum."""
    maxdepth = max(leaf.depth for leaf in ketbra.leaves)

    found = True
    while found:
        found = False
        for leaf in ketbra.leaves:
            if leaf.depth < maxdepth:
                found = True
                leaf.parent = None

    return ketbra


# * Tree-modifying functions
def readout(ketbra: KetBra, angle: Union[float, Tuple[float]]=0.0) -> KetBra:
    """Generate populations from excitations."""
    for kb in ketbra.leaves:
        excite(kb, 'mu', 'ket', True, angle)

    return prune(remove_nondiagonal(ketbra))


def excite(ketbra: KetBra, light_name: str, part: str='ket', readout: bool=False,
           angle: Union[float, Tuple[float]]=0.0) -> KetBra:
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
    angle : float or tuple of float, optional
        Either linear angle or a tuple of linear angle and circularity.

    Returns
    -------
    KetBra
    """
    # ket excitation
    if part not in ('ket', 'bra', 'both'):
        raise ValueError("`part` has to be either 'ket', 'bra' or 'both'")
    if isinstance(ketbra.ket, DiatomState):
        djs = [-1, 1]
    elif isinstance(ketbra.ket, SymTopState):
        djs = [-1, 0, 1]
    else:
        raise TypeError("`ketbra.ket` is neither DiatomState nor SymtopState")
    if part == 'ket' or part == 'both':
        for dnu in dnus:
            children = []
            nnu = ketbra.ket.nu+dnu
            if nnu < 0:
                continue
            for dj in djs:
                if ketbra.ket.j == 0 and dj == 0:
                    continue
                if isinstance(ketbra.ket, SymTopState) and ketbra.ket.k == 0 and dj == 0:
                    continue
                nj = ketbra.ket.j+dj
                if nj < 0:
                    continue
                if isinstance(ketbra.ket, SymTopState) and ketbra.ket.k > nj:
                    continue
                children.append(KetBra(ketbra.ket._replace(nu=nnu, j=nj), ketbra.bra))
            if dnu > 0:             # ket absorption
                LightInteraction(light_name, Side.KET, KSign.POS, readout, angle, parent=ketbra, children=children)
            elif dnu < 0:           # ket emission
                LightInteraction(light_name, Side.KET, KSign.NEG, readout, angle, parent=ketbra, children=children)

    if part == 'bra' or part == 'both':
        for dnu in dnus:
            children = []
            nnu = ketbra.bra.nu+dnu
            if nnu < 0:
                continue
            for dj in djs:
                if ketbra.bra.j == 0 and dj == 0:
                    continue
                if isinstance(ketbra.bra, SymTopState) and ketbra.bra.k == 0 and dj == 0:
                    continue
                nj = ketbra.bra.j+dj
                if nj < 0:
                    continue
                if isinstance(ketbra.bra, SymTopState) and ketbra.bra.k > nj:
                    continue
                children.append(KetBra(ketbra.ket, ketbra.bra._replace(nu=nnu, j=nj)))
            if dnu > 0:         # bra absorption
                LightInteraction(light_name, Side.BRA, KSign.NEG, readout, angle, parent=ketbra, children=children)
            elif dnu < 0:       # bra emission
                LightInteraction(light_name, Side.BRA, KSign.POS, readout, angle, parent=ketbra, children=children)

    return ketbra


def multi_excite(ketbra: KetBra, light_names: List[str], parts: Optional[List]=None,
                 light_angles: Optional[Union[List[float], List[Tuple[float]]]]=None) -> KetBra:
    """Generate multiple excitations of `ketbra`.

    Parameters
    ----------
    ketbra : KetBra
        State to excite.
    light_name : list of str
        Names of EM fields, length sets the number of excitations.
    parts : None or list
        None for first ket excitation and rest 'both' excitations.
    """
    if parts is not None and len(parts) != len(light_names):
        raise ValueError("len(parts) != len(light_names)")
    if parts is None:
        parts = ['ket']
        parts.extend(['both']*(len(light_names)-1))

    if light_angles is not None and len(light_angles) != len(light_names):
        raise ValueError("len(light_angles) != len(light_names)")
    if light_angles is None:
        light_angles = [0.0]*len(light_names)

    if light_names:
        excite(ketbra, light_names.pop(0), parts.pop(0), False, light_angles.pop(0))
        for kb in ketbra.leaves:
            multi_excite(kb, light_names[:], parts[:], light_angles[:])

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


def gen_excitations(root, light_names, parts, pols, meths=None) -> List[KetBra]:
    root = multi_excite(root, light_names, parts=parts, light_angles=pols[:3])
    if meths is not None:
        for meth in meths:
            root = meth(root)
    root = readout(root, pols[3])

    return root


def gen_pathways(jiter: Iterable, pols: Sequence,
                 meths: Optional[Sequence[Callable]]=None,
                 rotor: str='linear', kiter_func: Callable=None,
                 pump_overlap: bool=False) -> List[KetBra]:
    roots = gen_roots(jiter, rotor, kiter_func)
    pws = [gen_excitations(root, ['omg1', 'omg2', 'omg3'],
                           ['ket', 'both', 'both'], pols, meths)
           for root in roots]
    if pump_overlap:
        roots = gen_roots(jiter, rotor, kiter_func)
        pws.extend(
            [gen_excitations(root, ['omg2', 'omg1', 'omg3'],
                             ['ket', 'both', 'both'],
                             [pols[1], pols[0], pols[2], pols[3]],
                             meths)
             for root in roots])

    return pws


def geometric_factor(leaf: KetBra):
    """Return polarization and angular momentum sequence for a pathway.

    Same as in leaf_response."""
    kbs = leaf.ketbras()
    wkets, wbras = [], []
    for i in range(1, len(kbs)):
        kb, kbp = kbs[i], kbs[i-1]
        if kb.parent.side is Side.KET:
            wkets.insert(0, (kb.kj, kb.parent.angle))
        else:
            wbras.append((kb.parent.parent.bj, kb.parent.angle))
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
