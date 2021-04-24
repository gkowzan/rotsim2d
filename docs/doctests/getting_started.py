import rotsim2d.pathways as pw
import molspecutils.molecule as mol
kb1 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
print(kb1)

pw.excite(kb1, light_name='omg1')
kb1.print_tree()

pw.readout(kb1)
kb1.print_tree()
