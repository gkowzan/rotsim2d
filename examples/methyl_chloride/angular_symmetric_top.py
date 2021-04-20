# * Imports
from pathlib import Path
import numpy as np
from sqlalchemy import create_engine
import matplotlib.pyplot as plt
import spectroscopy.happier as h
from spectroscopy.molecule import CH3ClAlchemyMode
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.angular.symbolic as sym
import rotsim2d.angular.symbolic_results as symr
plt.ion()

# * Vibrational mode
sql_path = Path(h.hitran_cache) / 'CH3Cl_nu3.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))
ch3cl_mode = CH3ClAlchemyMode(engine)
T = 296.0

# * Distinct R-factors
# Check how many distinct R-factors are there.
rfactors_list = [(v['c0'], v['c12'], v['c13'], v['c14']) for v in symr.rfactors_dict.values()]
rfactors_pol_list = [(v['c12'], v['c13'], v['c14']) for v in symr.rfactors_dict.values()]
rfactors_pol_dict = {k: (v['c12'], v['c13'], v['c14']) for k, v in symr.rfactors_dict.items()}
rfactors_pol_classify = {}
for k, v in rfactors_pol_dict.items():
    rfactors_pol_classify.setdefault(v, []).append(k)

rfactors_highj_dict = {k: (v['c12'], v['c13'], v['c14']) for k, v in symr.rfactors_highj_dict.items()}
rfactors_highj_classify = {}
for k, v in rfactors_highj_dict.items():
    rfactors_highj_classify.setdefault(v, []).append(k)

# * Classify pathways
kbs = pw.gen_pathways([5], [0, 1, 2, 3], meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric', kiter_func=lambda x: [1])
pws = dl.Pathway.from_kb_list(kbs)
rfactors = [sym.dl_to_rfactor(pw, symr.rfactors_highj) for pw in pws]
classified = sym.classify_dls(pws, rfactors)
classified_states = {k: [pw.leaf.to_statelist(diatom=True, normalize=True) for pw in v]
                     for k, v in classified.items()}

# * Find zeroing angles
zeroing_angles_gk = sym.suppression_angles(classified.keys(), [0, sym.pi/4, sym.pi/2])
zeroing_angles_vaccaro = sym.suppression_angles(classified.keys(), [0, sym.pi/4, 0])

# * Classify Pathway's with respect to suppression angles
# Make a map between zeroing angles and the pathways they suppress.
angles_pws_gk = sym.classify_suppression(classified_states, zeroing_angles_gk)
angles_pws_vaccaro = sym.classify_suppression(classified_states, zeroing_angles_vaccaro)

dl.print_dl_dict(classified, fields=['peak', 'pols', 'geo_label', 'tw_coherence'])

# * Collect coefficients from `classified`
phis = [sym.phi, sym.phj, sym.phk, sym.phl]
T00_trigs = [trig.subs(dict(zip(phis, sym.thetas))) for trig in sym.T00_trigs]
classified_coeffs = [sym.collect(trig.subs(zip(T00_trigs, [sym.x1, sym.x2, sym.x3])),
                                 [sym.x1, sym.x2, sym.x3], evaluate=False)
                     for trig in classified.keys()]
classified_coeffs = [(d[sym.x1], d[sym.x2], d[sym.x3]) for d in classified_coeffs]
