"""Helper class and functions to more easily plot RC pathway interference."""
import itertools as it
from typing import Iterator, Any, Tuple, Union

import numpy as np
import scipy.constants as C
from molspecutils.molecule import (CH3ClAlchemyMode, COAlchemyMode,
                                   VibrationalMode)

import rotsim2d.dressedleaf as dl
import rotsim2d.pathways as pw
import rotsim2d.propagate as prop

TBs = {'CH3Cl': 1/(0.44340*100.0*C.c),
       'CO': 1/(1.93128087*100.0*C.c)}


def nth(it1: Iterator[Any], index: int) -> Any:
    return next(it.islice(it1, index, None), None)


class RCPeaks:
    def __init__(self, molecule: str, direction: str, j: int, k: int):
        if molecule == "CH3Cl":
            pws = pw.gen_pathways(
                [j], meths=[getattr(pw, 'only_'+direction),
                            pw.only_interstates],
                rotor='symmetric', kiter_func="[{:d}]".format(k))
            mode: VibrationalMode = CH3ClAlchemyMode()
        elif molecule == "CO":
            pws = pw.gen_pathways(
                [j], meths=[getattr(pw, 'only_'+direction),
                            pw.only_interstates],
                rotor='linear')
            mode: VibrationalMode = COAlchemyMode()
        else:
            raise ValueError("Unknown molecule")
        self.rc_peaks = dl.split_by_peaks(
            dl.DressedPathway.from_kb_list(pws, mode, 296.0),
            abstract=True)
        self.peaks = self.rc_peaks.keys()
        self.dps = self.rc_peaks.values()

    def __len__(self) -> int:
        return len(self.rc_peaks)

    def render(self, index: int) -> str:
        key = nth(iter(self.peaks), index)
        return "{branch:s}, {peak:s}, ({num:d})".format(
            branch=self.rc_peaks[key][0].peak_label,
            peak=str(key),
            num=len(self.rc_peaks[key]))

    def response(self, index: Union[int, Tuple[str, str]],
                 tws: np.ndarray) -> np.ndarray:
        if isinstance(index, tuple):
            pws = self.rc_peaks[index]
        else:
            pws = nth(iter(self.dps), index)
        resp = np.zeros((len(pws), tws.size), dtype=np.complex128)
        for i in range(len(pws)):
            resp[i] = prop.dressed_leaf_response(
                pws[i], [None, tws, None],
                ['t', 't', 't'], p=1e-4)

        return resp
