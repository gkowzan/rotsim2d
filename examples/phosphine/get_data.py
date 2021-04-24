# * Imports and constants
from pathlib import Path
import re
import appdirs
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as C
import molspecutils.foreign.hapi3 as h3
from molspecutils.data.pytips import (Tdat, TIPS_ISO_HASH, TIPS_GSI_HASH,
                                      TIPS_NPT)
import knickknacks.units as u

plt.ion()

dirs = appdirs.AppDirs('happier', 'gkowzan')
hitran_cache = str(Path(dirs.user_cache_dir) / 'db')
#: Boltzmann constant in erg/K
cBolts = 1.380648813E-16
Tref = 296.0

M, I = 28, 1
table_name = "{:d}_{:d}".format(M, I)
h3.db_begin(hitran_cache)
if table_name not in h3.getTableList():
    h3.fetch(table_name, M, I, 0, 99999, ParameterGroups=['160-char', 'Labels'])

def label_to_dict(s: str) -> dict:
    """Convert HITRAN state label to dict."""
    r = {k: v for k, v in [x.strip().split('=') for x in s.split(';')]}

    return r

def group4_llq_to_dict(s: str) -> dict:
    """Extract local quantum numbers from HITRAN string.

    For group 4: symmetric rotors."""
    return {'J': int(s[:3]),
            'K': int(s[3:6]),
            'l': int(s[6:8]),
            'C': s[8:10].strip(),
            'Sym': s[10].strip(),
            'F': s[11:].strip()}

def group4_glq_to_tuple(s: str) -> tuple:
    """Extract global quantum numbers from HITRAN string.

    For group 4: symmetric rotors."""
    return (int(s[5:7]), int(s[7:9]), int(s[9:11]), int(s[11:13]), int(s[13:15]) if s[13:15].strip() else None)


# * Get and visualize fundamental \nu_1 transitions
total_list = h3.getColumns('28_1', ['nu', 'a', 'gamma_air', 'delta_air', 'n_air', 'statep', 'statepp'])
total_list = [x for x in zip(*total_list)]

def v1_func(row):
    statepp = label_to_dict(row[-1])
    statep = label_to_dict(row[-2])
    cond = statep['ElecStateLabel'] == statepp['ElecStateLabel'] and\
        statep['v4'] == str(int(statepp['v4'])+1) and\
        statep['v2'] == statepp['v2'] and statep['v1'] == statepp['v1'] and\
        statep['v3'] == statepp['v3'] and statep['K'] == '0' and statepp['K'] == '0' and\
        statep['l'] == statepp['l'] and statep['vibSym'] == statepp['vibSym']

    return cond

v1_list = [x for x in total_list if v1_func(x)]
v1_sticks = ([x[0] for x in v1_list], [x[1] for x in v1_list])

fig, ax = plt.subplots()
ax.stem(v1_sticks[0], v1_sticks[1], markerfmt='C0.')
