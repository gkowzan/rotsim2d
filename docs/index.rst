.. _rotsim2d:

###################################################
rotsim2d: rotationally-resolved 2D infrared spectra
###################################################

**rotsim2d** is a library to model two-dimensional infrared (2DIR) rovibrational spectra.
It uses Mukamelian density-matrix perturbation theory treatment of third-order interaction.
It assumes explicit time-ordering of interactions, *i.e.* it is suitable for time-domain spectroscopy with time-separated optical pulses.

**rotsim2d** generates complete third-order excitation trees of rovibrational eigenstates, which can be manipulated and filtered in various ways.
Each pathway is also represented by an object, which provides its amplitude and polarization dependence and other information about the pathway.
This includes symbolic expressions for polarization and angular momentum dependence, which can be manipulated with SymPy library.
The user can generate and plot 2D resonance maps, which do not include information on the line shapes, but make it easire to visualize resonance amplitudes.
One can also simulate mixed time-frequency experimental signals which do include line shapes and which are directly comparable with measured signals.

**rotsim2d** also includes simple GUI applications that can be used to investigate branch structure of RR2DIR spectra, polarization dependence and evolution of rotational coherences during the waiting time between the second and third interaction.

.. image:: /images/cartoon_rcs_docs.png
   :width: 50%
   :align: center

.. toctree::
   :maxdepth: 2
   :hidden:

   getting-started/index
   user-guide/index
   api

Citations
=========
Please cite the following articles when publishing the results of using this library:

1. |LETTER|
2. |THEORY|
3. |HAPI|

License
=======

rotsim2d is available under the open source `MIT license <https://opensource.org/licenses/MIT>`_.
