"""Test DressedLeaf and PeakList with new SQLAlchemy classes."""
# * Imports
from pathlib import Path
from sqlalchemy import create_engine
import matplotlib.pyplot as plt
import molspecutils.happier as h
from molspecutils.molecule import CH3ClAlchemyMode
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
plt.ion()

# * Vibrational mode
sql_path = Path(h.hitran_cache) / 'CH3Cl.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))
ch3cl_mode = CH3ClAlchemyMode(engine)
T = 296.0

# * Pathways
pws = pw.gen_pathways(range(1, 37), [0]*4, meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric',
                      kiter_func=lambda x: range(x if x<10 else 10))
dressed_pws = dl.DressedPathway.from_kb_list(pws, ch3cl_mode, T)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)

# * Visualize
fig_dict = vis.plot2d_scatter(peaks, scatter_kwargs=dict(picker=True))

def scatter_onpick(event):
    """Show information about the peak pathway."""
    if event.artist != fig_dict['sc']:
        return
    for i, dl in enumerate(dls[event.ind[0]]):
        dl.pprint()
        if i == len(dls[event.ind[0]])-1:
            print()

fig_dict['fig'].canvas.mpl_connect('pick_event', scatter_onpick)
