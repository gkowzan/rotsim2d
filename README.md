Package to simulate 2D rovibrational spectra of gas-phase samples.

See examples in `examples/carbon_monoxide`.

# Dependencies
- numpy, scipy, anytree
- [pywigxjpf](https://gitlab.com/allisonlab/mdcs/pywigxjpf)
- [shed](https://gitlab.com/allisonlab/mdcs/shed)
- [spectroscopy](https://gitlab.com/allisonlab/mdcs/spectroscopy)

# Installation
Install the package by downloading or cloning the repo and calling the following inside main repository directory (containing `setup.py`):

``` sh
python -m pip install --no-deps -e .
```

or by installing directly from the repo with pip

``` sh
python -m pip git+ssh://git@gitlab.com:allisonlab/mdcs/rotsim2d.git@master
```

# TODO

## Speed up calculations

+ Transform the problem into n-dimensional array calculations.
+ Parallelize the code. Orders of magnitude improvement is needed, this won't be fast enough on its own.
+ Translate the code to C++. *Profile first.*
+ ~~Calculate only for initial ket excitation.~~
+ ~~Calculate directly in the frequency domain. This is not a general solution, can't be done if pulse-shape effects are included or dephasing function does not have analytical Fourier transform.~~

## Functionality

+ Make time-domain and frequency-domain results agree.
+ Constants.
+ Polarization and rotational coherence.
+ Add symmetric top and physical model.
+ Add line-mixing (proper dephasing) and Doppler.
+ Extract coherences and positions.
+ Generate real spectra (relative to linear absorption).
+ ~~Look at signals for different T_w times.~~ 
+ ~~Compare w/ overtones and wo/ overtones.~~
+ ~~Improve graph export.~~
+ ~~Add interstate coherences.~~
+ ~~Add filtering of pathways: k-vectors.~~
