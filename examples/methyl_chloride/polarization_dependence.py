"""Show polarization dependence of 2D symmetric top spectra."""
# * Imports
import numpy as np
import matplotlib.pyplot as plt
from molspecutils.molecule import CH3ClAlchemyMode
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import rotsim2d.symbolic.functions as s

# * Vibrational mode
ch3cl_mode = CH3ClAlchemyMode()
T = 296.0

# * Vaccaro scheme
# First and third pulse polarizations are the same
vac_angles_map = s.detection_angles([s.pi/2, s.pi/4, s.pi/2])
vac_angles = list(vac_angles_map.keys())

# ** XXXX
# - add/remove pw.only_twocolor to switch between two-color and three-color
# - pump_overlap=True gives additional time ordering
#   the only effect on 2D spectra is that now we also have negative frequencies
# HITRAN does not contain spectroscopic data for all k values, so I'm limiting
# them with `kiter_func`
kbs = pw.gen_pathways(
    range(1, 10),
    meths=[pw.only_SII, pw.only_twocolor],
    rotor="symmetric",
    kiter_func=lambda x: range(x if x < 10 else 10),
    pump_overlap=False,
)
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)
fig_dict = vis.plot2d_scatter(peaks)

# *** Peaks zeroed by angle
# Show which angles will be zeroed by one of angles in `vac_angles`.
# vac_peaks_by_angles[angle] contains a list of peaks that will be zeroed. This is
# useful for three-color scheme in which many pathways overlap and it is not
# always clear that some pathways were in fact supressed.
vac_peaks_by_angles = dl.split_by_equiv_peaks(vac_angles_map, peaks, dls)

# Mark the zeroed pathways by black crosses
ang = vac_angles[1]
fig_dict['ax'].scatter(vac_peaks_by_angles[ang].probes, vac_peaks_by_angles[ang].pumps,
                       s=10.0, c='black', marker='x')

# ** Some other detection angle
# change last element of `pols` argument to some value from `angles` list to
# zero some other class of pathways
kbs = pw.gen_pathways(
    range(1, 10),
    meths=[pw.only_SII, pw.only_twocolor],
    rotor="symmetric",
    kiter_func=lambda x: range(x if x < 10 else 10),
    pump_overlap=False,
)
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(
    dressed_pws,
    angles=[np.pi / 2, np.pi / 4, np.pi / 2, np.float64(ang)],
    return_dls=True,
)
fig_dict = vis.plot2d_scatter(peaks)
fig_dict['ax'].set(xlim=(714, 742), ylim=(722, 742))

# * GK scheme
# All angles different
# ** Angles
gk_angles_map = s.detection_angles([0, s.pi/4, s.pi/2])
gk_angles = list(gk_angles_map.keys())

# ** XXXX
kbs = pw.gen_pathways(range(1, 10), meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric',
                      kiter_func=lambda x: range(x if x<10 else 10))
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)
gk_peaks_by_angles = dl.split_by_equiv_peaks(gk_angles_map, peaks, dls)

fig_dict = vis.plot2d_scatter(peaks)
fig_dict['ax'].set(title='XXXX', xlim=(714, 742), ylim=(722, 742))

ang = gk_angles[4]
fig_dict['ax'].scatter(gk_peaks_by_angles[ang].probes, gk_peaks_by_angles[ang].pumps,
                       s=10.0, c='black', marker='x')

# ** Some other detection angle
# change last element of `pols` argument to some value from `angles` list to
# zero some other class of pathways
kbs = pw.gen_pathways(
    range(1, 10),
    meths=[pw.only_SII, pw.only_twocolor],
    rotor="symmetric",
    kiter_func=lambda x: range(x if x < 10 else 10),
    pump_overlap=False,
)
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, T)
peaks, dls = dl.peak_list(
    dressed_pws, angles=[0.0, np.pi / 4, np.pi / 2, np.float64(ang)], return_dls=True
)
fig_dict = vis.plot2d_scatter(peaks)
fig_dict['ax'].set(xlim=(714, 742), ylim=(722, 742))
plt.show()
