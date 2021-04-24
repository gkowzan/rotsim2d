# * Imports and functions
from typing import Tuple
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
from molspecutils import happier
import knickknacks.units as u
import rotsim2d.propagate as prop
prop.ignore_missing = False
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import gkpywigxjpf as wig
plt.ion()
PAPER_OUTPUT = Path('/mnt/d/DFCS/ARTYKULY/2021 rotsim PRL/figures')

def co_props(co_env, I):
    co_levels = happier.energy_levels([0,1,2], 5, I)
    co_pops = happier.equilibrium_pops(co_levels, co_env['T'], 5, I)
    co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2), (0, 0), (1, 1), (2, 2)], 5, 1)
    co_params_suppl = happier.generate_line_params(
        [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
    co_params_suppl2 = happier.generate_line_params([(0, 1, 0), (1, 2, 0)], co_levels, co_params)
    co_params.update(co_params_suppl)
    co_params.update(co_params_suppl2)
    co_params = happier.apply_env(co_params, co_env, to_nu=True)

    return {'elevels': co_levels, 'populations': co_pops, 'line_params': co_params}


def response_calc(pws, props, resp_xs, scale=1.0):
    root_prop = prop.Spectrum(props, filter=lambda kb: True,
                              freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
    multi_prop = prop.MultiPropagator(pws, root_prop)
    with wig.wigxjpf(300, 6):
        resp_xs += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])*scale

    return resp_xs


def response_calc_t2(pws, props, resp_xs, scale=1.0):
    root_prop = prop.Spectrum(props, filter=lambda kb: True,
                              freq_shift=[0.0, 0.0, u.wn2nu(1800.0)])
    multi_prop = prop.MultiPropagator(pws, root_prop)
    with wig.wigxjpf(300, 6):
        resp_xs += multi_prop.response([None, fs2[:, np.newaxis], fs[np.newaxis, :]])*scale

    return resp_xs


def load_npz(fpath, name):
    temp = np.load(fpath)
    arr = temp[name]
    temp.close()
    
    return arr

iso1_props = co_props({'p': 15.0,  'T': 296.0}, 1)

# * Response function parameters
df = 10e9
F = u.wn2nu(700.0)
N = np.round(F/df).astype(np.int)
fs = np.arange(-N, N+1)*df
fs_cm = u.nu2wn(fs)

F2 = u.wn2nu(100.0)
N2 = np.round(F2/df).astype(np.int)
fs2 = np.arange(-N2, N2+1)*df
fs2_cm = u.nu2wn(fs2)

jmax = 30

# * Base case
# ** XXXX (SI)
pws = pw.gen_pathways(range(jmax), [0, np.pi/4, np.pi/2, np.arctan(2/3)], iso1_props['populations'], [pw.only_SI])
resp_xs_xxxx_SI = np.zeros((fs2.size, fs.size), dtype=np.complex)
resp_xs_xxxx_SI = response_calc_t2(pws, iso1_props, resp_xs_xxxx_SI)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(fs2_cm, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_xxxx_SI),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
# fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
# fig_dict['fig'].savefig(NIST_OUTPUT / 'full_absorptive_spectrum.png', dpi=300.0)
# None
