"""Compare rotationally-resolved signals for different pathways and directions."""
# * Imports and functions
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyfftw.interfaces.scipy_fftpack as fftp
import rotsim2d.pathways as pw
import rotsim2d.propagate as prop
import rotsim2d.visual as vis
import shed.units as u
from shed.experiment import find_index
from spectroscopy import happier
plt.ion()
OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide')

# * System properties
pressure, Tgas, length = 20.0, 296.0, 3.0
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, Tgas, 5, 1)
co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2)], 5, 1)
co_params_suppl = happier.generate_line_params(
    [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
co_params.update(co_params_suppl)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env, to_nu=True)

# * Response function parameters
dt = 1.0/(2.0*u.wn2nu(343.0))/2.0
T = 1.0/co_params[((0,0), (1,1))]['gam']*4.0/5.0         # relaxation time constant
N = np.round(T/dt).astype(np.int)
t2 = 0.0
ts = np.arange(N)*dt

# * Calculate molecular coherence evolution
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0)])
resp_xs_freqs = root_prop.cross_section_freqs([ts])
resp_xs_freqs_wn = u.nu2wn(resp_xs_freqs)

# ** First order
resp_xs_spec_o1 = np.zeros(ts.size, dtype=np.complex)
for j in range(5, 6):
    # first order
    root = pw.KetBra(0, j, 0, j)
    root = pw.excite(root, 'omg1', 'both')
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    resp_xs_spec_o1 += pop_factor*root_prop.cross_section_spectrum(root, [ts])

# ** Third order
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

resp_xs = np.zeros((ts.size, ts.size), dtype=np.complex)
for j in range(0, 16):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    # root = pw.excite(root, 'omg1', 'both')
    # root = pw.remove_interstate(root)
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    for i, t1 in enumerate(ts):
        resp_xs[i, :] += pop_factor*root_prop.response(root, [t1, t2, ts])

# *** Fourier transform in pump axis
resp_xs_freqs_wn = fftp.fftshift(u.nu2wn(root_prop.cross_section_freqs([ts])))
resp_xs_spec_fft = fftp.fftshift(fftp.fft2(resp_xs, planner_effort='FFTW_ESTIMATE'))
if resp_xs_freqs_wn.size % 2 == 0:
    # Make positive and negative frequencies symmetric around DC.
    resp_xs_freqs_wn = resp_xs_freqs_wn[1:]
    resp_xs_spec_fft = resp_xs_spec_fft[1:, 1:]
# Invert the order elements to agree with Mukamel's FT convention.
resp_xs_spec_fft = resp_xs_spec_fft[::-1, ::-1]

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn+1800.0, resp_xs_freqs_wn+1800.0),
                         spec2d=prop.absorptive(resp_xs_spec_fft),
                         scale='linear')
fig_dict['fig'].canvas.set_window_title('Time domain -- absorptive -- all pathways')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'mandal_full_linear.png'))

# ** Third order -- nointerstates
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

resp_xs = np.zeros((ts.size, ts.size), dtype=np.complex)
for j in range(0, 16):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    # root = pw.excite(root, 'omg1', 'both')
    root = pw.remove_interstates(root)
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    for i, t1 in enumerate(ts):
        resp_xs[i, :] += pop_factor*root_prop.response(root, [t1, t2, ts])

# *** Fourier transform in pump axis
resp_xs_freqs_wn = fftp.fftshift(u.nu2wn(root_prop.cross_section_freqs([ts])))
resp_xs_spec_fft_noint = fftp.fftshift(fftp.fft2(resp_xs, planner_effort='FFTW_ESTIMATE'))
if resp_xs_freqs_wn.size % 2 == 0:
    # Make positive and negative frequencies symmetric around DC.
    resp_xs_freqs_wn = resp_xs_freqs_wn[1:]
    resp_xs_spec_fft_noint = resp_xs_spec_fft_noint[1:, 1:]
# Invert the order elements to agree with Mukamel's FT convention.
resp_xs_spec_fft_noint = resp_xs_spec_fft_noint[::-1, ::-1]

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn+1800.0, resp_xs_freqs_wn+1800.0),
                         spec2d=prop.absorptive(resp_xs_spec_fft_noint),
                         scale='linear')
fig_dict['fig'].canvas.set_window_title('Time domain -- absorptive -- no interstates')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'mandal_nointerstates_linear.png'))

# ** Third order -- no overtones
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

resp_xs = np.zeros((ts.size, ts.size), dtype=np.complex)
for j in range(0, 16):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    # root = pw.excite(root, 'omg1', 'both')
    root = pw.remove_overtones(root)
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    for i, t1 in enumerate(ts):
        resp_xs[i, :] += pop_factor*root_prop.response(root, [t1, t2, ts])

# *** Fourier transform in pump axis
resp_xs_freqs_wn = fftp.fftshift(u.nu2wn(root_prop.cross_section_freqs([ts])))
resp_xs_spec_fft_noover = fftp.fftshift(fftp.fft2(resp_xs, planner_effort='FFTW_ESTIMATE'))
if resp_xs_freqs_wn.size % 2 == 0:
    # Make positive and negative frequencies symmetric around DC.
    resp_xs_freqs_wn = resp_xs_freqs_wn[1:]
    resp_xs_spec_fft_noover = resp_xs_spec_fft_noover[1:, 1:]
# Invert the order elements to agree with Mukamel's FT convention.
resp_xs_spec_fft_noover = resp_xs_spec_fft_noover[::-1, ::-1]

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn+1800.0, resp_xs_freqs_wn+1800.0),
                         spec2d=prop.absorptive(resp_xs_spec_fft_noover),
                         scale='linear')
fig_dict['fig'].canvas.set_window_title('Time domain -- absorptive -- no overtones')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'mandal_noovertones_linear.png'))

# ** Compare w/ and wo/ overtones
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn+1800.0, resp_xs_freqs_wn+1800.0),
                         spec2d=np.imag(resp_xs_spec_fft-resp_xs_spec_fft_noover))
fig_dict['fig'].suptitle('purely absorptive')
fig_dict['fig'].canvas.set_window_title('Time domain -- absorptive')
