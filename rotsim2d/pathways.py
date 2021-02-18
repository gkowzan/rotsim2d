r"""Generate all Liouville pathways for nth order rovibrational excitation.
"""
# * Imports, constants and enums
import enum
from copy import deepcopy
import operator as op
from functools import reduce
from typing import List, Union, Tuple, Dict, Iterable, Sequence, Optional
import numpy as np

import anytree as at
from anytree.exporter import UniqueDotExporter

dnus = [-1, +1] #: possible changes of vibrational state
djs = [-1, +1] #: possible changes of rotational state

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
    def __init__(self, name: str, side: Side, sign: KSign, readout: bool=False, angle: Union[float, Tuple[float]]=0.0,
                 parent=None, children=None):
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
    def __init__(self, knu: int, kj: int, bnu: int, bj: int, pop: float=1.0, parent=None, children=None):
        super(KetBra, self).__init__()
        self.separator = "->"
        self.name = "|{:d},{:d}><{:d},{:d}|".format(knu, kj, bnu, bj)
        self.knu, self.kj, self.bnu, self.bj = knu, kj, bnu, bj
        self.parent = parent
        if children:
            self.children = children

        self.pop = pop                # fractional occupancy

    def print_tree(self):
        for pre, _, node in at.RenderTree(self):
            treestr = "{:s}{:s}".format(pre, node.name)
            print(treestr)

    def __repr__(self):
        return "Ketbra({:d}, {:d}, {:d}, {:d})".format(self.knu, self.kj, self.bnu, self.bj)

    def __str__(self):
        return self.name

    def __eq__(self, o):
        if not isinstance(o, KetBra):
            return NotImplemented
        else:
            return self.knu == o.knu and self.kj == o.kj and self.bnu == o.bnu and self.bj == o.bj

    def copy(self):
        return deepcopy(self)

    def savepng(self, path):
        UniqueDotExporter(self, nodeattrfunc=nodeattrfunc).to_picture(str(path))

        return path
    
    def conj(self):
        """Return conjugate of this KetBra."""
        return KetBra(self.bnu, self.bj, self.knu, self.kj)

    def kb_ancestor(self, ancestor=None):
        """Return first KetBra ancestor."""
        if ancestor is None:
            ancestor = self.ancestors[0]
        if isinstance(ancestor, KetBra):
            return ancestor
        else:
            self.kb_ancestor(ancestor.ancestors[0])

    def diagonals(self, sort=False, rot=True):
        """Collect diagonal leaves."""
        pops = [x for x in self.leaves if x.knu == x.bnu and x.kj == x.bj]
        if sort:
            pops.sort(key=lambda x: x.kj)
            pops.sort(key=lambda x: x.knu)

        return pops

    def is_diagonal(self):
        return self.knu == self.bnu and self.kj == self.bj

    def is_rephasing(self):
        ints = self.interactions()

        return ints[0].sign != ints[2].sign

    def is_esa(self):
        return (self.knu != self.root.knu) and (self.bnu != self.root.knu)

    def is_overtone(self):
        return abs(self.knu-self.bnu)>1

    def has_overtone(self):
        return len([kb for kb in self.ketbras() if kb.is_overtone()]) > 0

    def is_between(self, kb1, kb2):
        """Check if this pathway produces cross-peak between `kb1` and `kb2`."""
        pump_kb = self.ketbras()[1]
        return (kb1 == pump_kb or kb1 == pump_kb.conj()) and (kb2 == self or kb2 == self.conj())

    def is_dfwm(self):
        """Check is this pathway contains only coherences corresponding to a single dipole transition."""
        kbs = self.ketbras()
        return all([kb.is_diagonal() or kb == kbs[1] or kb == kbs[1].conj() for kb in kbs])

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

    def ksigns(self) -> List[KSign]:
        return [i.sign for i in self.interactions()]

    def sides(self) -> List[Side]:
        return [i.side for i in self.interactions()]

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


def only_between(ketbra: KetBra, kb1: KetBra, kb2: KetBra) -> KetBra:
    """Limit tree to pathways bewteen `kb1` and `kb2`."""
    for l in ketbra.leaves:
        if not l.is_between(kb1, kb2):
            l.parent = None

    return prune(ketbra)


def only_dfwm(ketbra: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_dfwm():
            l.parent = None

    return prune(ketbra)


def only_pathway(ketbra: KetBra, pathway: KetBra) -> KetBra:
    for l in ketbra.leaves:
        if not l.is_pathway(pathway):
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

    return ketbra


def excite(ketbra: KetBra, light_name: str, part: str='ket', readout: bool=False,
           angle: Union[float, Tuple[float]]=0.0) -> KetBra:
    """Generate all excitations of `ketbra`.

    Modifies `ketbra` in place.

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
    if part == 'ket' or part == 'both':
        for dnu in dnus:
            children = []
            nnu = ketbra.knu+dnu
            if nnu < 0:
                continue
            for dj in djs:
                nj = ketbra.kj+dj
                if nj < 0:
                    continue
                children.append(KetBra(nnu, nj, ketbra.bnu, ketbra.bj))
            if dnu > 0:             # ket absorption
                LightInteraction(light_name, Side.KET, KSign.POS, readout, angle, parent=ketbra, children=children)
            elif dnu < 0:           # ket emission
                LightInteraction(light_name, Side.KET, KSign.NEG, readout, angle, parent=ketbra, children=children)

    if part == 'bra' or part == 'both':
        for dnu in dnus:
            children = []
            nnu = ketbra.bnu+dnu
            if nnu < 0:
                continue
            for dj in djs:
                nj = ketbra.bj+dj
                if nj < 0:
                    continue
                children.append(KetBra(ketbra.knu, ketbra.kj, nnu, nj))
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
