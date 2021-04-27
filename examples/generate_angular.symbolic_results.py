# * Imports
from pprint import pprint, pformat
from rotsim2d.symbolic.functions import *
import rotsim2d.symbolic.results as symr
from rotsim2d.dressedleaf import geometric_labels
init_printing(pretty_print=False)
print("Script used to generate contents of rotsim2d.symbolic.results module.")
print("Does not print the results, source code demonstrates how to use")
print("rotsim2d.symbolic.funtions for e.g. custom derivations.")

# * Generate all G-factors
# ** General G-factors
gfactors = {}
for js, label in geometric_labels.items():
    args = [J_i+j for j in js]
    gfactors[label] = tuple((factor(powdenest(gfactor_expr(*(args + [i])), force=True), deep=True)
                             for i in range(3)))

# ** High-J G-factors
gfactors_highj = {k: tuple(nsimplify(limit(x*(2*J_i+1)**(3/2), J_i, oo)) for x in v)
                  for k, v in gfactors.items()}

# * Generate all R-factors
rfactors = {k: rfactorize(gfactors[k], symr.T00_exprs) for k in gfactors.keys()}
rfactors_dict = {k: rfactorize(gfactors[k], symr.T00_exprs, cfac=False, coeffs=True) for k in gfactors.keys()}
rfactors_xxxx = {k: factor(powdenest(v.subs({phi: S(0), phj: S(0), phk: S(0), phl: S(0)}), force=True), deep=True)
                 for k, v in rfactors.items()}
rfactors_relative = {k: rfactorize(gfactors[k], symr.T00_exprs, cfac=False, relative=True) for k in gfactors.keys()}
rfactors_highj = {k: rfactorize(gfactors_highj[k], symr.T00_exprs, cfac=False) for k in gfactors_highj.keys()}
rfactors_highj_dict = {k: rfactorize(gfactors_highj[k], symr.T00_exprs, cfac=False, coeffs=True)
                       for k in gfactors_highj.keys()}
