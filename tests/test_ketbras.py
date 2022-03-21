"""Some sanity checks for ketbras."""
import pytest
from anytree.iterators.preorderiter import PreOrderIter
from anytree.walker import WalkError, Walker
import rotsim2d.pathways as pw
import molspecutils.molecule as mol

@pytest.fixture
def kbs_tree():
    return pw.gen_pathways([5], rotor='symmetric',
                           kiter_func="[1]", pump_overlap=False)[0]


def test_leafs_are_ketbras(kbs_tree):
    for l in kbs_tree.leaves:
        assert isinstance(l, pw.KetBra)


def test_first_interaction_is_absorption(kbs_tree):
    for li in kbs_tree.root.children:
        assert li.is_absorption()


def test_four_interactions(kbs_tree):
    for l in kbs_tree.leaves:
        assert len(l.interactions()) == 4


def test_flip_readout(kbs_tree):
    for l in kbs_tree.leaves:
        side = l.parent.side
        absorption = l.parent.is_absorption()
        pw.flip_readout(l)
        assert l.parent.side != side
        assert l.parent.is_absorption() != absorption


def test_conjugate_chain(kbs_tree):
    conj = pw.conjugate_chain(kbs_tree.leaves[0])
    orig_chain = [kbs_tree.leaves[0]] + list(kbs_tree.ancestors)
    conj_chain = [conj] + list(conj.ancestors)
    for o, c in zip(orig_chain, conj_chain):
        assert o == c.conj()
    with pytest.raises(WalkError):
        walker = Walker()
        walker.walk(conj, kbs_tree.leaves[0])


def test_ketbra_evolve(kbs_tree):
    for l in kbs_tree.leaves:
        assert l.evolve() == l


def test_five_ketbras(kbs_tree):
    for l in kbs_tree.leaves:
        assert len(l.ketbras()) == 5


def test_rot_states(kbs_tree):
    for kb in PreOrderIter(kbs_tree, filter_=lambda x: isinstance(x, pw.KetBra)):
        assert isinstance(kb.ket, (mol.DiatomState, mol.SymTopState))
        assert isinstance(kb.bra, (mol.DiatomState, mol.SymTopState))


def test_leaves_are_diagonal(kbs_tree):
    for kb in kbs_tree.leaves:
        assert kb.is_diagonal()


def test_statelist_test_structure(kbs_tree):
    for kb in kbs_tree.leaves:
        sl = kb.to_statelist()
        assert len(sl) == 5
        for pair in sl:
            assert len(pair) == 2
            assert isinstance(pair[0], (mol.DiatomState, mol.SymTopState))
            assert isinstance(pair[1], (mol.DiatomState, mol.SymTopState))
