"""Test DressedLeaf and PeakList with new SQLAlchemy classes."""
# * Imports
import matplotlib.pyplot as plt
import knickknacks.units as u
from molspecutils.molecule import CH3ClAlchemyMode
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis

# * Vibrational mode
print('Initializing vibrational mode')
ch3cl_mode = CH3ClAlchemyMode()
T = 296.0

# * Pathways
print('Calculating peak list')
pws = pw.gen_pathways(range(1, 37), meths=[pw.only_SII], rotor='symmetric',
                      kiter_func=lambda x: range(x if x<10 else 10))
dressed_pws = dl.DressedPathway.from_kb_list(pws, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)

# * Visualize
fig_dict = vis.plot2d_scatter(peaks, scatter_kwargs=dict(picker=True))


def scatter_onpick(event):
    """Show information about the peak pathway."""
    if event.artist != fig_dict['sc']:
        return
    dl.pprint_dllist(dls[event.ind[0]])


fig_dict['fig'].canvas.mpl_connect('pick_event', scatter_onpick)
plt.show()
