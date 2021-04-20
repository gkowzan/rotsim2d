"""Test DressedLeaf and PeakList with new SQLAlchemy classes."""
# * Imports
from pathlib import Path
from sqlalchemy import create_engine
import matplotlib.pyplot as plt
import spectroscopy.happier as h
import spectroscopy.molecule
from spectroscopy.molecule import COAlchemyMode, DiatomState
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import rotsim2d.propagate as prop
import pywigxjpf as wig
plt.ion()

# * Vibrational mode
sql_path = Path(h.hitran_cache) / 'CO.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))
co_mode = COAlchemyMode(engine)
T = 296.0

# * Pathways
pws = pw.gen_pathways(range(10), [0]*4, meths=[pw.only_SII], rotor='linear')
dressed_pws = dl.dress_pws(pws, co_mode, T)
peaks = dl.peak_list(dressed_pws)

# * Visualize
vis.plot2d_scatter(peaks)

# * Response function
import shed.units as u
import numpy as np
import time

df = 10e9
F = u.wn2nu(700.0)
N = np.round(F/df).astype(np.int64)
fs = np.arange(-N, N+1)*df
fs_cm = u.nu2wn(fs)

t = time.time()
resp = np.zeros((fs.size, fs.size), dtype=np.complex128)
with wig.wigxjpf(300, 6):
    for dressed_leaf in dressed_pws:
        resp += prop.dressed_leaf_response(
            dressed_leaf, coords=[fs[:, None], None, fs[None, :]],
            domains=['f', 't', 'f'],
            freq_shifts=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
print(time.time()-t)
fig_dict = vis.plot2d_im([fs_cm, fs_cm], np.imag(resp))
