r"""Associate polarizations, angular momentum factors and pathways amplitudes
with pathways, i.e. everything but line shapes."""
from typing import Optional, Iterable, Tuple, List
from collections import namedtuple
import numpy as np
import scipy.constants as C
from rotsim2d.pathways import KetBra, Side
from rotsim2d.couple import four_couple
from pywigxjpf import wigxjpf
import shed.units as u


class DressedLeaf:
    """Assign polarizations in proper order, line intensities, coherence frequencies
    and allow for calculating the intensity of the pathway. Do everything except
    calculating the associated line shape/response time dependence."""
    fields = ['peak', 'const', 'nus', 'gams', 'js', 'pols', 'tw_coherence', 'diagonal', 'geo_label']

    def __init__(self, leaf: KetBra, sys_params):
        kb_series = [x for x in leaf.ancestors if isinstance(x, KetBra)] + [leaf]
        self.peak = (kb_series[1].name, kb_series[3].name)
        wkets, wbras = [], []
        self.nus, self.gams = [], []
        self.const = np.complex(1.0)
        for i in range(1, len(kb_series)):
            kb, kbp = kb_series[i], kb_series[i-1]

            # dipole interaction reduced matrix element
            if kb.parent.side is Side.KET:
                pair = ((kbp.knu, kbp.kj), (kb.knu, kb.kj))
                wkets.insert(0, (kb.kj,  kb.parent.angle))
            else:
                pair = ((kbp.bnu, kbp.bj), (kb.bnu, kb.bj))
                wbras.append((kb.parent.parent.bj, kb.parent.angle))
            try:
                mu = sys_params['line_params'][pair]['mu']
            except KeyError:
                mu = sys_params['line_params'][(pair[1], pair[0])]['mu']

            if kb.parent.readout:
                self.const *= mu*kb.root.pop
                wbras.extend(wkets)
                self.js = tuple(x[0] for x in wbras)
                self.pols = tuple(x[1] for x in wbras)
                break
            self.const *= 1.0j/C.hbar*kb.parent.side*mu

            pair = ((kb.bnu, kb.bj), (kb.knu, kb.kj))
            if pair in sys_params['line_params']:
                self.nus.append(sys_params['line_params'][pair]['nu'])
                self.gams.append(sys_params['line_params'][pair]['gam'])
            elif pair[::-1] in sys_params['line_params']:
                self.nus.append(-sys_params['line_params'][pair[::-1]]['nu'])
                self.gams.append(sys_params['line_params'][pair[::-1]]['gam'])
            else:
                raise KeyError("Missing `nu` and `gam` for {!s}".format(pair))
        self.tw_coherence = not (self.nus[1] == 0.0)
        self.diagonal = self.peak[0] == self.peak[1]
        self.geo_label = geometric_label(self)

    def intensity(self, tw: Optional[float]=None) -> float:
        """Intensity of the pathway."""
        ret = self.const*four_couple(self.js, self.pols)
        if tw is not None:
            ret *= np.exp(-2.0*np.pi*tw*(1.0j*self.nus[1] + self.gams[1]))

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

    def __str__(self):
        return self.custom_str()

    def __repr__(self):
        return self.custom_str()


def split_by_js(kbl: Iterable[DressedLeaf]):
    ret = {}
    for dl in kbl:
        ret.setdefault(dl.js, []).append(dl)
    return ret


def perm_pols(pols: Tuple):
    return (pols[2], pols[3], pols[0], pols[1])


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


def geometric_label(dl: DressedLeaf):
    ji = dl.js[0]
    if dl.js[1] == ji-1 and dl.js[2] == ji-2 and dl.js[3] == ji-1:
        return "O"
    if dl.js[1] == ji+1 and dl.js[2] == ji+2 and dl.js[3] == ji+1:
        return "S"
    if dl.js[1] == ji-1 and dl.js[2] == ji and dl.js[3] == ji-1:
        return "PP"
    if dl.js[1] == ji+1 and dl.js[2] == ji and dl.js[3] == ji+1:
        return "RR"
    if dl.js[1] == ji-1 and dl.js[2] == ji and dl.js[3] == ji+1:
        return "PR"
    if dl.js[1] == ji+1 and dl.js[2] == ji and dl.js[3] == ji-1:
        return "RP"


def split_by_pols_highjs(kbl: Iterable[DressedLeaf]):
    sopr = {'S', 'O', 'PR', 'RP'}
    rrpp = {'RR', 'PP'}
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


def dress_pws(pws, sys_params):
    return sum(([DressedLeaf(l, sys_params) for l in root.leaves] for root in pws),
               start=[])


# * Peaks without line shapes
Peak2D = namedtuple("Peak2D", "pump_wl probe_wl sig")

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


def peak_list(ll: List[DressedLeaf], tw: Optional[float]=0.0) -> Peak2DList:
    """Create a list of 2D peaks from a list of DressedLeaf's."""
    ll = split_by_peaks(ll)
    pl = Peak2DList()
    for peak, dll in ll.items():
        pu, pr = u.nu2wn(dll[0].nus[0]), u.nu2wn(dll[0].nus[2])
        sig = 0.0
        for dl in dll:
            with wigxjpf(300, 6):
                sig += np.imag(dl.intensity(tw=tw))
        pl.append(Peak2D(pu, pr, sig))
    pl.sort(key=lambda x: abs(x.sig))

    return pl
