"""Some sanity checks for ketbras."""
import pytest
from anytree.iterators.preorderiter import PreOrderIter 
import rotsim2d.pathways as pw
import molspecutils.molecule as mol

@pytest.fixture
def kbs_tree():
    return pw.gen_pathways([5], [0, 1, 2, 3], rotor='symmetric',
                           kiter_func=lambda x: [1], pump_overlap=False)[0]


def test_leafs_are_ketbras(kbs_tree):
    for l in kbs_tree.leaves:
        assert isinstance(l, pw.KetBra)


def test_four_interactions(kbs_tree):
    for l in kbs_tree.leaves:
        assert len(l.interactions()) == 4


def test_five_ketbras(kbs_tree):
    for l in kbs_tree.leaves:
        assert len(l.ketbras()) == 5


def test_rot_states(kbs_tree):
    for kb in PreOrderIter(kbs_tree, filter_=lambda x: isinstance(x, pw.KetBra)):
        assert isinstance(kb.ket, mol.RotState)
        assert isinstance(kb.bra, mol.RotState)


def test_leaves_are_diagonal(kbs_tree):
    for kb in kbs_tree.leaves:
        assert kb.is_diagonal()


def test_statelist_test_structure(kbs_tree):
    for kb in kbs_tree.leaves:
        sl = kb.to_statelist()
        assert len(sl) == 5
        for pair in sl:
            assert len(pair) == 2
            assert isinstance(pair[0], mol.RotState)
            assert isinstance(pair[1], mol.RotState)
