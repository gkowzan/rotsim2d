# * Imports
import matplotlib.pyplot as plt
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.symbolic.functions as sym
import rotsim2d.symbolic.results as symr

# * Classify pathways
# - add/remove pw.only_twocolor to switch between two-color and three-color
# - pump_overlap=True gives additional time ordering (pump pulses overlapping in
#   time, probe separate). The only effect on 2D spectra is that now we also have
#   negative frequencies
kbs = pw.gen_pathways([5], meths=[pw.only_SII],
                      rotor='symmetric', kiter_func=lambda x: [1], pump_overlap=False)
pws = dl.Pathway.from_kb_list(kbs)
rfactors = [sym.RFactor.from_pathway(pw, True, True) for pw in pws]
classified = sym.classify_dls(pws, rfactors)
classified_states = {k: [pw.leaf.to_statelist(diatom=True, normalize=True) for pw in v]
                     for k, v in classified.items()}

print("Polarization response expressions")
dl.print_dl_dict(classified, fields=['peak', 'light_inds', 'geo_label', 'trans_label', 'tw_coherence'])
print()

# * Find zeroing angles
zeroing_angles_gk = sym.suppression_angles(classified.keys(), [0, sym.pi/4, sym.pi/2])
zeroing_angles_vaccaro = sym.suppression_angles(classified.keys(), [sym.pi/2, sym.pi/4, sym.pi/2])

# * Classify Pathway's with respect to suppression angles
# Make a map between zeroing angles and the pathways they suppress.
angles_pws_gk = sym.classify_suppression(classified, zeroing_angles_gk)
angles_pws_vaccaro = sym.classify_suppression(classified, zeroing_angles_vaccaro)

# this will print out expressions for polarization dependence and pathways which
# correspond to each expression
print("Suppression angles:")
dl.print_dl_dict(angles_pws_vaccaro, fields=['peak', 'light_inds', 'geo_label', 'trans_label', 'tw_coherence'])
