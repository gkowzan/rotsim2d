# * Imports and constants
from pathlib import Path
import numpy as np
import scipy.constants as C
import matplotlib.pyplot as plt
import rotsim2d.rcs as rcs
plt.ion()
OUTPUT = Path('~/notes/.attach/sb/:8c55d411-1a6f-4227-bf33-653ce4043c67').expanduser()

# trans-stilbene
A = 2.678e9/C.c*1e-2
B = 0.256e9/C.c*1e-2

# CO
# B = 1.9224966626745492                         # CO, cm-1
# A = 0.0

# * Simulate and plot
T = 50.0
Jmax = rcs.est_jmax(A, T)
ts = np.linspace(0, 5e-9, 10000)
# parallel = rcs.decay_total(Jmax, A, B, T, ts, False)

fig, ax = plt.subplots()
ax.plot(ts*1e12, rcs.normalized(rcs.decay_total(Jmax, A, B, T, ts, False)), label='B')
ax.plot(ts*1e12, rcs.normalized(rcs.decay_total(Jmax, A, B/5, T, ts, False)), label='B/5')
ax.plot(ts*1e12, rcs.normalized(rcs.decay_total(rcs.est_jmax(A, 100.0), A, B, 100.0, ts, False)),
        label='B, T=100 K')
ax.set(xlabel='Time (ps)', ylabel='RCS', xlim=(0, 70),
       title='trans-stilbene @ T={:.0f} K'.format(T))
ax.legend(loc='best')
fig.savefig(OUTPUT / 'trans-stilbene.png')
