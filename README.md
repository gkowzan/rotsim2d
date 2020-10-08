Package to simulate 2D rovibrational spectra of gas-phase samples.

# TODO

## Speed up calculations

+ Transform the problem into n-dimensional array calculations.
+ Parallelize the code. Orders of magnitude improvement is needed, this won't be fast enough on its own.
+ Translate the code to C++. *Profile first.*
+ ~~Calculate only for initial ket excitation.~~
+ ~~Calculate directly in the frequency domain. This is not a general solution, can't be done if pulse-shape effects are included or dephasing function does not have analytical Fourier transform.~~

## Functionality

+ Add symmetric top and physical model.
+ Polarization and rotational coherence.
+ Add line-mixing (proper dephasing) and Doppler.
+ Improve graph export.
+ Extract coherences and positions.
+ Generate real spectra (relative to linear absorption).
+ Constants.
+ ~~Add interstate coherences.~~
+ ~~Add filtering of pathways: k-vectors.~~
