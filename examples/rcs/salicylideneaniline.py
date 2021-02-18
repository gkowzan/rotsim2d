# * Imports
from typing import Union
from pathlib import Path
import xarray as xr
import numpy as np
import scipy.constants as C
import scipy.optimize as opt
import matplotlib.pyplot as plt
import rotsim2d.rcs as rcs
plt.ion()
OUTPUT = Path('~/notes/.attach/sb/:8c55d411-1a6f-4227-bf33-653ce4043c67').expanduser()
INPUT = Path('/mnt/d/fizyka/DFCS/POMIARY/CETAS/January2021/Jan272021')

# * Functions
def read_pump_probe(fpath: Union[str, Path], anisotropy: bool=True, magic: bool=False, not_anisotropy: bool=False) -> xr.DataArray:
    fpath = Path(fpath)
    np_data = np.genfromtxt(fpath)
    xr_data = xr.DataArray(np_data[:, 1:],
                           coords=[('time', np_data.T[0]),
                                   ('param', ['par', 'par_err',
                                              'perp', 'perp_err'])],
                           dims=('time', 'param'),)
    if anisotropy:
        xr_r = r(xr_data)
        xr_data = xr.concat((xr_data, xr_r), dim='param')
    if magic:
        xr_ma = ma(xr_data)
        xr_data = xr.concat((xr_data, xr_ma), dim='param')
    if not_anisotropy:
        xr_notr = notr(xr_data)
        xr_data = xr.concat((xr_data, xr_notr), dim='param')
    return xr_data


def r(xr_decay: xr.DataArray) -> xr.DataArray:
    xr_data = (xr_decay.loc[:, 'par']-xr_decay.loc[:, 'perp'])/(xr_decay.loc[:, 'par']+2*xr_decay.loc[:, 'perp'])
    return xr_data.expand_dims({'param': ['r']}, 0)


def notr(xr_decay: xr.DataArray) -> xr.DataArray:
    xr_data = -(xr_decay.loc[:, 'par']-xr_decay.loc[:, 'perp'])
    return xr_data.expand_dims({'param': ['notr']}, 0)


def ma(xr_decay: xr.DataArray) -> xr.DataArray:
    xr_data = (xr_decay.loc[:, 'par']+2*xr_decay.loc[:, 'perp'])
    return xr_data.expand_dims({'param': ['ma']}, 0)

# * Load data
SA_Ar = {3.5: read_pump_probe(INPUT / 'SA_Ar_3_5PSI.txt', magic=True, not_anisotropy=True),
         12.3: read_pump_probe(INPUT / 'SA_Ar_12_3PSI.txt', magic=True, not_anisotropy=True),
         31.3: read_pump_probe(INPUT / 'SA_Ar_31_3PSI.txt', magic=True, not_anisotropy=True)}
SA_He = {3.6: read_pump_probe(INPUT / 'SA_He_3_6PSI.txt', magic=True, not_anisotropy=True),
         12.9: read_pump_probe(INPUT / 'SA_He_12_9PSI.txt', magic=True, not_anisotropy=True)}

# * Plot
# ** Parallel
fig, ax = plt.subplots()
for k, v in SA_Ar.items():
    ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'par'],
            label='Ar, {:.1f} psi'.format(k))
for k, v in SA_He.items():
    ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'par'],
            label='He, {:.1f} psi'.format(k))
ax.set(ylim=(-0.1, 3.0), title='SA-Ar/He', xlabel='Time (ps)', ylabel='Parallel')
ax.legend(loc='best')

# ** Perpendicular
fig, ax = plt.subplots()
for k, v in SA_Ar.items():
    ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'perp'],
            label='Ar, {:.1f} psi'.format(k))
for k, v in SA_He.items():
    ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'perp'],
            label='He, {:.1f} psi'.format(k))
ax.set(ylim=(-0.1, 3.0), title='SA-Ar/He', xlabel='Time (ps)', ylabel='Parallel')
ax.legend(loc='best')

# ** Anisotropy
fig, ax = plt.subplots()
for k, v in SA_Ar.items():
    ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'r'],
            label='Ar, {:.1f} psi'.format(k))
# for k, v in SA_He.items():
#     ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'r'],
#             label='He, {:.1f} psi'.format(k))
ax.set(ylim=(-0.5, 1.2), title='SA-Ar/He', xlabel='Time (ps)', ylabel='Anisotropy')
ax.legend(loc='best')

# ** Magic angle
fig, ax = plt.subplots()
for k, v in SA_Ar.items():
    ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'ma'],
            label='Ar, {:.1f} psi'.format(k))
# for k, v in SA_He.items():
#     ax.plot(v.loc[0.0:].coords['time'], v.loc[0.0:, 'r'],
#             label='He, {:.1f} psi'.format(k))
ax.set(title='SA-Ar/He', xlabel='Time (ps)', ylabel='MA')
ax.legend(loc='best')

# * Rotational coherence
A = 0.10621666
B = (0.04030298+0.02921689)/2
T = 5.0

ts = SA_Ar[3.5].loc[0.0:].coords['time'].data*1e-12
fig, ax = plt.subplots()
ax.plot(ts*1e12, rcs.decay_total_anisotropy(rcs.est_jmax(A, T), A, B, T, ts))
ax.plot(ts*1e12, SA_Ar[3.5].loc[0.0:, 'r'].data)
ax.plot(SA_Ar[12.3].loc[0.0:].coords['time'].data, SA_Ar[12.3].loc[0.0:, 'r'].data)
ax.plot(SA_Ar[31.3].loc[0.0:].coords['time'].data, SA_Ar[31.3].loc[0.0:, 'r'].data)
ax.set(ylim=(-0.1, 1.2))
