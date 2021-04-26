import rotsim2d.pathways as pw
import molspecutils.molecule as mol
kb1 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
print(kb1)

pw.excite(kb1, light_name='omg1')
kb1.print_tree()

pw.readout(kb1)
kb1.print_tree()

kb2 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
kb2 = pw.multi_excite(kb2, light_names=['omg1', 'omg2', 'omg3'],
                      parts=['ket', 'both', 'both'],
                      light_angles=[0]*3)
print()
kb2.print_tree()
kb2.savepng('../images/kb2.png')

import rotsim2d.visual as vis
latex_code = vis.tikz_diagrams(kb2)
vis.latex_compile('../images/kb2_tikz.tex', vis.LATEX_PRE + latex_code + vis.LATEX_POST)

kb2 = pw.only_SII(kb2)
kb2 = pw.readout(kb2)
kb2.savepng('../images/kb2_SII.png')

# kb2 = pw.only_SI(kb2)
# kb2 = pw.readout(kb2)
# kb2.savepng('../images/kb2_SI.png')

from pprint import pprint
js = (0, 1, 2, 3)
kbs = pw.gen_pathways(js, pols=[0]*4, meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric', kiter_func=lambda j: range(j+1))
pprint(kbs)
