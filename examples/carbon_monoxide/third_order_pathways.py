from pathlib import Path
import rotsim2d.pathways as pw
import rotsim2d.visual as vis
import rotsim2d.couple as cp
import numpy as np
from pywigxjpf import wigxjpf

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide')
j = 5
ma = 54.7356103172453*np.pi/180.0 # magic angle

# * Check G factor
with wigxjpf(300, 6):
    # print(cp.G(5,6,5,4,2))
    print(cp.four_couple_linear([5,6,7,6], [np.pi+ma, ma, 0.0, 0.0]))

# * Investigate graphs
# ** All
root_all = pw.KetBra(0, j, 0, j)
root_all = pw.multi_excite(root_all, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[1, 2, 3])
root_all = pw.readout(root_all, angle=4)
root_all.savepng(str(OUTPUT / '3rd_noesa.png'))

factors = pw.tensor_analysis(root_all)
factors_by_js = pw.split_by_js(factors)
js_list = list(factors_by_js.keys())

# ** S_III
root_nonreph = pw.KetBra(0, j, 0, j)
root_nonreph = pw.multi_excite(root_nonreph, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'],
                               light_angles=[1, 2, 3])
root_nonreph = pw.only_SIII(root_nonreph)
root_nonreph = pw.readout(root_nonreph, 4)
vis.latex_compile(OUTPUT / 'diagrams/SIII_tikz.tex',
                  vis.LATEX_PRE + vis.tikz_diagrams(
                      root_nonreph, direction=False, abstract=(0, j), hspace="1.2cm")
                  + vis.LATEX_POST)

# ** S_II
root_nonreph = pw.KetBra(0, j, 0, j)
root_nonreph = pw.multi_excite(root_nonreph, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'],
                               light_angles=[1, 2, 3])
root_nonreph = pw.only_SII(root_nonreph)
root_nonreph = pw.readout(root_nonreph, 4)
vis.latex_compile(OUTPUT / 'diagrams/SII_tikz.tex',
                  vis.LATEX_PRE + vis.tikz_diagrams(
                      root_nonreph, direction=False, abstract=(0, j), hspace="1.2cm")
                  + vis.LATEX_POST)

# ** S_I
root_nonreph = pw.KetBra(0, j, 0, j)
root_nonreph = pw.multi_excite(root_nonreph, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'],
                               light_angles=[1, 2, 3])
root_nonreph = pw.only_SI(root_nonreph)
root_nonreph = pw.readout(root_nonreph, 4)
vis.latex_compile(OUTPUT / 'diagrams/SI_tikz.tex',
                  vis.LATEX_PRE + vis.tikz_diagrams(
                      root_nonreph, direction=False, abstract=(0, j), hspace="1.2cm")
                  + vis.LATEX_POST)

# ** Non-rephasing
root_nonreph = pw.KetBra(0, j, 0, j)
root_nonreph = pw.multi_excite(root_nonreph, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'],
                               light_angles=[1, 2, 3])
root_nonreph = pw.remove_rephasing(root_nonreph)
root_nonreph = pw.readout(root_nonreph, 4)
factors = pw.tensor_analysis(root_nonreph)
# factors_by_js = pw.split_by_js(factors)
# factors_by_crosspeaks = pw.split_by_crosspeaks(factors)
root_nonreph.savepng(str(OUTPUT / '3rd_nonrephasing.png'))
vis.latex_compile(OUTPUT / 'diagrams/3rd_nonrephasing_tikz.tex',
                  vis.LATEX_PRE + vis.tikz_diagrams(
                      root_nonreph, direction=False, abstract=(0, j), hspace="1.2cm")
                  + vis.LATEX_POST)

root_nonreph = pw.KetBra(0, j, 0, j)
root_nonreph = pw.multi_excite(root_nonreph, ['omg1', 'omg2', 'omg3'],
                            parts=['ket', 'both', 'both'])
root_nonreph = pw.remove_nonrephasing(pw.remove_esa(root_nonreph))
root_nonreph = pw.remove_overtones(root_nonreph)
root_nonreph = pw.remove_nondiagonal(pw.readout(root_nonreph))
root_nonreph.savepng(str(OUTPUT / "3rd_noesa_nonrephasing_noovertones.png"))

# ** Rephasing
root_reph = pw.KetBra(0, j, 0, j)
root_reph = pw.multi_excite(root_reph, ['omg1', 'omg2', 'omg3'],
                            parts=['ket', 'both', 'both'])
root_reph = pw.remove_rephasing(pw.remove_esa(root_reph))
root_reph = pw.remove_nondiagonal(pw.readout(root_reph))
root_reph.savepng(str(OUTPUT / "3rd_noesa_rephasing.png"))

root_reph = pw.KetBra(0, j, 0, j)
root_reph = pw.multi_excite(root_reph, ['omg1', 'omg2', 'omg3'],
                            parts=['ket', 'both', 'both'])
root_reph = pw.remove_rephasing(pw.remove_esa(root_reph))
root_reph = pw.remove_overtones(root_reph)
root_reph = pw.remove_nondiagonal(pw.readout(root_reph))
root_reph.savepng(str(OUTPUT / "3rd_noesa_rephasing_noovertones.png"))

# ** P5 P5
# *** All
j = 5
root_all = pw.KetBra(0, j, 0, j)
root_all = pw.multi_excite(root_all, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
root_all = pw.only_between(root_all, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
root_all = pw.remove_nondiagonal(pw.readout(root_all))
root_all.savepng(str(OUTPUT / '3rd_P5_P5.png'))
vis.latex_compile(OUTPUT / 'diagrams/3rd_P5_P5_tikz.tex',
                  vis.LATEX_PRE + vis.tikz_diagrams(root_all, direction=True) + vis.LATEX_POST)


# *** Rephasing
root_all = pw.KetBra(0, j, 0, j)
root_all = pw.multi_excite(root_all, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
root_all = pw.only_between(root_all, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
root_all = pw.remove_nonrephasing(root_all)
root_all = pw.remove_nondiagonal(pw.readout(root_all))
root_all.savepng(str(OUTPUT / '3rd_P5_P5_rephasing.png'))

# *** Non-rephasing
root_all = pw.KetBra(0, j, 0, j)
root_all = pw.multi_excite(root_all, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
root_all = pw.only_between(root_all, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
root_all = pw.remove_rephasing(root_all)
root_all = pw.remove_nondiagonal(pw.readout(root_all))
root_all.savepng(str(OUTPUT / '3rd_P5_P5_nonrephasing.png'))

# *** No interstates
j = 5
root_all = pw.KetBra(0, j, 0, j)
root_all = pw.multi_excite(root_all, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
root_all = pw.only_between(root_all, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
root_all = pw.remove_interstates(root_all)
root_all = pw.remove_nondiagonal(pw.readout(root_all))
root_all.savepng(str(OUTPUT / '3rd_P5_P5_no_interstates.png'))
