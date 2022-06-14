from typing import Sequence, Tuple

import pytest
gkpywigxjpf = pytest.importorskip("gkpywigxjpf")
wig = gkpywigxjpf.gkpywigxjpf
import rotsim2d.dressedleaf as dl
import rotsim2d.pathways as pw
import rotsim2d.symbolic.functions as sym


def Gxjpf(ji: int, jj: int, jk: int, jl: int, k: int):
    """Wigner-6j part of four-fold reduced matrix element."""
    return (2*k+1)*wig.wig6j(k, k, 0, ji, ji, jk)*\
        wig.wig6j(1, 1, k, jk, ji, jj)*\
        wig.wig6j(1, 1, k, ji, jk, jl)


def gfactors_xjpf(js: Sequence[int]) -> Tuple[float]:
    """G-factors calculated with wigxjpf."""
    js = list(js)
    return (Gxjpf(*(js + [0])), Gxjpf(*(js + [1])), Gxjpf(*(js + [2])))


@pytest.fixture
def pws():
    """Generate a list of pathways for further testing."""
    kbs = pw.gen_pathways([0, 1, 2], rotor='symmetric',
                          kiter_func="range(j+1)", pump_overlap=False)
    pws = dl.Pathway.from_kb_list(kbs)

    return pws


def test_wigner_numeric_symbolic(pws):
    """Compare G-factor `evalf`ed and numpy'ed."""
    for pw in pws:
        try:
            analytical = sym.dl_gfactors(pw)
            numeric = pw.gfactors()
            for n, a in zip(numeric, analytical):
                assert n == pytest.approx(a)
        except ZeroDivisionError:
            pw.pprint()
            raise


def test_wigner_wigxjpf_numeric(pws):
    """Compare G-factor based on wigxjp and analytical expressions."""
    wig.table_init(1000, 9)
    wig.temp_init(1000)
    for pw in pws:
        try:
            wigxjpf = gfactors_xjpf(pw.js)
            numeric = pw.gfactors()
            for w, n in zip(wigxjpf, numeric):
                assert n == pytest.approx(w)
        except ZeroDivisionError:
            pw.pprint()
            print(wigxjpf)
            raise
    wig.temp_free()
    wig.table_free()


def test_polarization_numeric_symbolic(pws):
    for pw in pws:
        numeric = pw.T00s([0, 1, 2, 3])
        analytical = sym.dl_T00s(pw, [0, 1, 2, 3])
        for n, a in zip(numeric, analytical):
            assert n == pytest.approx(a)
