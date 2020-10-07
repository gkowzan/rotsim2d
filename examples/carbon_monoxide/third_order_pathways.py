from pathlib import Path
import rotsim2d.pathways as pw

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide')
j = 5

# * Investigate graphs
# ** No ESA
root_all = pw.KetBra(0, j, 0, j)
root_all = pw.multi_excite(root_all, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
root_all = pw.remove_esa(root_all)
root_all = pw.remove_nondiagonal(pw.readout(root_all))
root_all.savepng(str(OUTPUT / '3rd_noesa.png'))

# ** Non-rephasing
root_nonreph = pw.KetBra(0, j, 0, j)
root_nonreph = pw.multi_excite(root_nonreph, ['omg1', 'omg2', 'omg3'],
                            parts=['ket', 'both', 'both'])
root_nonreph = pw.remove_nonrephasing(pw.remove_esa(root_nonreph))
root_nonreph = pw.remove_nondiagonal(pw.readout(root_nonreph))
root_nonreph.savepng(str(OUTPUT / '3rd_noesa_nonrephasing.png'))

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
