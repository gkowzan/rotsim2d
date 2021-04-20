# * Imports
from pathlib import Path
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

# * Classify pathways
# - add/remove pw.only_twocolor to switch between two-color and three-color
# - pump_overlap=True gives additional time ordering (pump pulses overlapping in
#   time, probe separate). The only effect on 2D spectra is that now we also have
#   negative frequencies
kbs = pw.gen_pathways([5], [0, 1, 2, 3], meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric', kiter_func=lambda x: [1], pump_overlap=False)
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

# this will print out expressions for polarization dependence and pathways which
# correspond to each expression
dl.print_dl_dict(classified, fields=['peak', 'angles', 'geo_label', 'tw_coherence'])
