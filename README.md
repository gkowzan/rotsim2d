Package to simulate 2D rovibrational spectra of gas-phase samples.

See examples in `examples` directory:
- `methyl_chloride/classify_pathways.py` prints out a mapping between
  polarization-dependence expressions and pathways they correspond to. Modify
  `gen_pathways` call to obtain classifications for different experimental
  schemes.
- `methyl_chloride/peak_picker_three_color.py` and
  `methyl_chloride/peak_picker_three_color.py` allow one to inspect 2D peak
  scatter plots. Clicking on a peak prints out double-sided Feynmann diagrams
  that contribute to the selected peak. See comments in the code to modify what
  gets plotted.
- `methyl_chloride/polarization_dependence.py` allows one to inspect
  polarization dependece in the case wher two beam polarizations are the same
  and in the case when all polarizations are different. See comments in the code
  to modify the suppressed pathways and excitation scheme.

# Dependencies
- numpy, scipy, anytree, sympy
- [pywigxjpf](https://gitlab.com/allisonlab/mdcs/pywigxjpf)
- [shed](https://gitlab.com/allisonlab/mdcs/shed)
- [spectroscopy](https://gitlab.com/allisonlab/mdcs/spectroscopy)

# Installation
Install the package by downloading or cloning the repo and calling the following inside main repository directory (containing `setup.py`):

``` sh
python -m pip install -e .
```

or by installing directly from the repo with pip

``` sh
python -m pip install git+ssh://git@gitlab.com/allisonlab/mdcs/rotsim2d.git@symmetrictop
```

[spectroscopy](https://gitlab.com/allisonlab/mdcs/spectroscopy) package needs to
be initialized by executing `spectroscopy/scripts/initialize_dbs.py` script to
download spectroscopic data and fill the local database.

# TODO
+ Refactor angular momentum code.

## Speed up calculations
+ Use explicit formulas for G-factors instead of calculating Wigner coefficients numerically.
+ Transform the problem into n-dimensional array calculations.
+ Parallelize the code. Orders of magnitude improvement is needed, this won't be fast enough on its own.
+ Translate the code to C++. *Profile first.*
+ ~~Calculate only for initial ket excitation.~~
+ ~~Calculate directly in the frequency domain. This is not a general solution, can't be done if pulse-shape effects are included or dephasing function does not have analytical Fourier transform.~~

## Functionality

+ Make time-domain and frequency-domain results agree.
+ Constants.
+ Add symmetric top and physical model.
+ Add line-mixing (proper dephasing) and Doppler.
+ Generate real spectra (relative to linear absorption).
+ ~~Extract coherences and positions.~~
+ ~~Polarization and rotational coherence.~~
+ ~~Look at signals for different T_w times.~~ 
+ ~~Compare w/ overtones and wo/ overtones.~~
+ ~~Improve graph export.~~
+ ~~Add interstate coherences.~~
+ ~~Add filtering of pathways: k-vectors.~~
