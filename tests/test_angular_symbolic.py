import pytest
from rotsim2d.dressedleaf import geometric_labels
import rotsim2d.symbolic.results as symr
import rotsim2d.symbolic.functions as sym
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl

def test_gfactors():
    for js, label in geometric_labels.items():
        args = [symr.J_i+j for j in js]
        gfacs = tuple(
            (symr.factor(symr.powdenest(sym.gfactor_expr(*(args + [i])), force=True), deep=True)
             for i in range(3)))
        for i in range(3):
            assert symr.factor(symr.powdenest(symr.nsimplify(
                gfacs[i]-symr.gfactors[label][i]), force=True), deep=True) == symr.S(0)


def test_gfactors_highj():
    gfactors_highj = {k: tuple(symr.nsimplify(symr.limit(x*(2*symr.J_i+1)**(3/2), symr.J_i, symr.oo)) for x in v)
                      for k, v in symr.gfactors.items()}
    for k in gfactors_highj.keys():
        for i in range(3):
            assert symr.factor(symr.nsimplify(
                symr.gfactors_highj[k][i]-gfactors_highj[k][i]), deep=True) == symr.S(0)

@pytest.fixture
def pws():
    """Generate a list of pathways for further testing."""
    kbs = pw.gen_pathways([5], rotor='symmetric',
                          kiter_func="[1]", pump_overlap=False)
    pws = dl.Pathway.from_kb_list(kbs)

    return pws


def test_rfactors(pws):
    for p in pws:
        print(p.trans_label)
        rfactor = sym.RFactor.from_pathway(p)
        assert rfactor == sym.RFactor(symr.rfactors_dict[p.trans_label],
                                      angles='experimental')

def test_rfactors_highj(pws):
    for p in pws:
        rfactor = sym.RFactor.from_pathway(p, highj=True)
        assert rfactor == sym.RFactor(symr.rfactors_highj_dict[p.trans_label],
                                      angles='experimental')


