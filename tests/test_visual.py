from pathlib import Path
import pytest
import rotsim2d.pathways as pw
import rotsim2d.visual as vis
import molspecutils.molecule as mol

@pytest.fixture
def kb2():
    kb2 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
    kb2 = pw.multi_excite(kb2, light_names=['omg1', 'omg2', 'omg3'],
                          parts=['ket', 'both', 'both'],
                          light_angles=[0]*3)
    return kb2


def test_latex_tikz(kb2, tmp_path):
    latex_code = vis.tikz_diagrams(kb2)
    latex_path = str(Path(tmp_path) / 'kb2_tikz.tex')
    vis.latex_compile(latex_path, vis.LATEX_PRE+latex_code+vis.LATEX_POST)
