import pytest
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.symbolic.functions as sym


@pytest.fixture
def pws():
    """Generate a list of pathways for further testing."""
    kbs = pw.gen_pathways([5], [0, 1, 2, 3], rotor='symmetric',
                          kiter_func=lambda x: [1], pump_overlap=False)
    pws = dl.Pathway.from_kb_list(kbs)

    return pws


def test_wigner_numeric_analytical(pws):
    for pw in pws:
        numeric = pw.gfactors()
        analytical = sym.dl_gfactors(pw)
        for n, a in zip(numeric, analytical):
            assert n == pytest.approx(a) 


def test_polarization_numeric_analytical(pws):
    for pw in pws:
        numeric = pw.T00s()
        analytical = sym.dl_T00s(pw)
        for n, a in zip(numeric, analytical):
            assert n == pytest.approx(a)
