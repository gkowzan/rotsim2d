# * Imports
import re
from rotsim2d.symbolic.functions import *
import rotsim2d.symbolic.results as symr
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
from rotsim2d.dressedleaf import geometric_labels
from sympy.printing.python import PythonPrinter
from pprint import pprint as pprint2
from yapf.yapflib.yapf_api import FormatCode
init_printing(pretty_print=False)
print("Script used to generate contents of rotsim2d.symbolic.results module.")

# * Printer
pyprinter = PythonPrinter()
filter_phis = re.compile(r"\\phi_([ijkl])")
filter_thetas = re.compile(r"\\theta_([ijkl])")

def pformat(name):
    obj = globals()[name]
    exprp = pyprinter._str(pyprinter.doprint(obj))

    return FormatCode("{:s} = {:s}".format(
        name, re.sub(filter_thetas, r"theta_\1",
                     re.sub(filter_phis, r"ph\1", exprp))))[0]


def pprint(name):
    print(pformat(name))

# * Generate all G-factors
gfactors = {}
for js, label in geometric_labels.items():
    args = [J_i+j for j in js]
    gfactors[label] = tuple((factor(powdenest(gfactor_expr(*(args + [i])),
                                              force=True), deep=True)
                             for i in range(3)))
pprint('gfactors')

gfactors_highj = {k: tuple(nsimplify(limit(x*(2*J_i+1)**(3/2), J_i, oo)) for x in v)
                  for k, v in gfactors.items()}
pprint('gfactors_highj')

gfactors_highj_numeric = {k: np.array(tuple(float(x) for x in v))
                          for k, v in gfactors_highj.items()}
with np.printoptions(sign=' '):
    pprint2(gfactors_highj_numeric)

# * Polarization expressions
phis = [phi, phj, phk, phl]
T00_exprs = [nsimplify(simplify(T00_phis(k, phis).rewrite(cos)))
             for k in (0, 1, 2)]
pprint('T00_exprs')

T00_theta_exprs = [nsimplify(simplify(T00_phis(k, thetas).rewrite(cos)))
                   for k in (0, 1, 2)]
pprint('T00_theta_exprs')

# * Generate all R-factors
kbs = pw.gen_pathways([5], rotor='symmetric',
                      kiter_func=lambda x: [1])
pws = dl.Pathway.from_kb_list(kbs)
rfactors = [RFactor.from_pathway(pw) for pw in pws]
rfactors_dict = {p.trans_label: rfactor.dict for p, rfactor
                 in zip(pws, rfactors)}
pprint('rfactors_dict')

rfactors_highj = [RFactor.from_pathway(pw, highj=True) for pw in pws]
rfactors_highj_dict = {p.trans_label: rfactor.dict for p, rfactor
                       in zip(pws, rfactors_highj)}
pprint('rfactors_highj_dict')
