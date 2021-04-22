"""Show polarization dependence of 2D symmetric top spectra."""
# * Imports
from pathlib import Path
import numpy as np
from sqlalchemy import create_engine
import matplotlib.pyplot as plt
import spectroscopy.happier as h
from spectroscopy.molecule import CH3ClAlchemyMode
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import rotsim2d.angular.symbolic_results as symr
plt.ion()

# * Vibrational mode
sql_path = Path(h.hitran_cache) / 'CH3Cl.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))
ch3cl_mode = CH3ClAlchemyMode(engine)
T = 296.0

# * Vaccaro scheme
# Two angles are the same
# ** Angles
#: Zeroing angles in Vaccaro scheme
angles = list(symr.angles_pws_vaccaro.keys())

# ** XXXX
# - add/remove pw.only_twocolor to switch between two-color and three-color
# - pump_overlap=True gives additional time ordering
#   the only effect on 2D spectra is that now we also have negative frequencies
# HITRAN does not contain spectroscopic data for all k values, so I'm limiting
# them with `kiter_func`
kbs = pw.gen_pathways(range(1, 10), [0]*4, meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric', kiter_func=lambda x: range(x if x<10 else 10),
                      pump_overlap=True)
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)
fig_dict = vis.plot2d_scatter(peaks)

# *** Peaks zeroed by angle
# Show which angles will be zeroed by one of angles in `angles`.
# peaks_by_angles[angle] contains a list of peaks that will be zeroed. This is
# useful for three-color scheme in which many pathways overalp and it is not
# always clear that some pathways were in fact supressed.
peaks_by_angles = {}
for angle in symr.angles_pws_vaccaro:
    peaks_by_angles[angle] = dl.Peak2DList()
    for statelist in symr.angles_pws_vaccaro[angle]:
        peaks_by_angles[angle].extend(dl.equiv_peaks(statelist, peaks, dls))
    peaks_by_angles[angle].sort_by_sigs()


# Mark the zeroed pathways by black crosses
ang = angles[1]
fig_dict['ax'].scatter(peaks_by_angles[ang].probes, peaks_by_angles[ang].pumps,
                       s=10.0, c='black', marker='x')
 
# ** Some other detection angle
# change last element of `pols` argument to some value from `angles` list to
# zero some other class of pathways
kbs = pw.gen_pathways(range(1, 10), pols=[0.0, np.pi/4, 0.0, np.arctan(2)],
                      meths=[pw.only_SII, pw.only_twocolor], rotor='symmetric',
                      kiter_func=lambda x: range(x if x<10 else 10),
                      pump_overlap=False)
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)
fig_dict = vis.plot2d_scatter(peaks)
fig_dict['ax'].set(xlim=(714, 742), ylim=(722, 742))

# * GK scheme
# All angles different
# ** Angles
angles = list(symr.angles_pws_gk.keys())

# ** XXXX
kbs = pw.gen_pathways(range(1, 10), [0]*4, meths=[pw.only_SII], rotor='symmetric',
                      kiter_func=lambda x: range(x if x<10 else 10))
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)

peaks_by_angles = {}
for angle in symr.angles_pws_gk:
    peaks_by_angles[angle] = dl.Peak2DList()
    for statelist in symr.angles_pws_gk[angle]:
        peaks_by_angles[angle].extend(dl.equiv_peaks(statelist, peaks, dls))
    peaks_by_angles[angle].sort_by_sigs()

fig_dict = vis.plot2d_scatter(peaks)
fig_dict['ax'].set(title='XXXX', xlim=(714, 742), ylim=(722, 742))

ang = angles[3]
fig_dict['ax'].scatter(peaks_by_angles[ang].probes, peaks_by_angles[ang].pumps,
                       s=10.0, c='black', marker='x')

# ** Some other detection angle
# change last element of `pols` argument to some value from `angles` list to
# zero some other class of pathways
kbs = pw.gen_pathways(range(1, 10), pols=[0.0, np.pi/4, np.pi/2, np.arctan(2)],
                      meths=[pw.only_SII, pw.only_twocolor], rotor='symmetric',
                      kiter_func=lambda x: range(x if x<10 else 10),
                      pump_overlap=False)
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)
fig_dict = vis.plot2d_scatter(peaks)
fig_dict['ax'].set(xlim=(714, 742), ylim=(722, 742))
