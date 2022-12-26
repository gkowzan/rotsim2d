.. |LETTER| replace:: G.\ Kowzan and T. K. Allison, “Controlling Rotationally Resolved Two-Dimensional Infrared Spectra with Polarization,” J. Phys. Chem. Lett., vol. 13, no. 50, pp. 11650–11654, Dec. 2022, doi: 10.1021/acs.jpclett.2c03331.

.. |THEORY| replace:: G.\ Kowzan and T. K. Allison, “Theory of rotationally resolved two-dimensional infrared spectroscopy including polarization dependence and rotational coherence dynamics,” Phys. Rev. A, vol. 106, no. 4, p. 042819, Oct. 2022, doi: 10.1103/PhysRevA.106.042819.

.. |HAPI| replace:: R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcisło, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016).

.. |HITRAN| replace:: I.E. Gordon et al., The HITRAN2016 molecular spectroscopic database, J. Quant. Spectrosc. Radiat. Transfer 203, 3-69 (2017).

.. |MUKAMEL| replace:: S.\ Mukamel, "Principles of nonlinear spectroscopy". Oxford University Press, New York, 1995.

.. |HAMM| replace:: P.\ Hamm, M. Zanni,"Concepts and methods of 2D infrared spectroscopy". Cambridge University Press, 2011.

###################################################
rotsim2d: rotationally-resolved 2D infrared spectra
###################################################

.. image:: https://zenodo.org/badge/503824832.svg
   :target: https://zenodo.org/badge/latestdoi/503824832

**rotsim2d** is a library to model two-dimensional infrared (2DIR) rovibrational spectra.
It uses density-matrix perturbation theory treatment of third-order interaction.
It assumes explicit time-ordering of interactions, *i.e.* it is suitable for time-domain spectroscopy with time-separated optical pulses.

**rotsim2d** generates complete third-order excitation trees of rovibrational eigenstates, which can be manipulated and filtered in various ways.
Each pathway is also represented by an object, which provides its amplitude and polarization dependence and other information about the pathway.
This includes symbolic expressions for polarization and angular momentum dependence, which can be manipulated with SymPy library.
The user can generate and plot 2D resonance maps, which do not include information on the line shapes, but make it easier to visualize resonance amplitudes.
One can also simulate mixed time-frequency experimental signals which do include line shapes and which are directly comparable with measured signals.

**rotsim2d_apps** includes simple GUI applications that can be used to investigate branch structure of RR2DIR spectra, polarization dependence and evolution of rotational coherences during the waiting time between the second and third interaction.

.. image:: docs/images/cartoon_rcs_docs.png
   :width: 50%
   :align: center

Install
=======
The package is available on PyPI and can be most easily installed with `pip`::

  pip install rotsim2d

This will install the library.
To install GUI applications install the following package::

  pip install rotsim2d_apps

Documentation
=============
Documentation is available at `<https://rotsim2d.readthedocs.io/>`_.

Citations
=========
Please cite the following articles when publishing the results of using this library:

1. |LETTER|
2. |THEORY|
3. |HAPI|

License
=======
rotsim2d is available under the open source `MIT license <https://opensource.org/licenses/MIT>`_.

Funding
=======
.. list-table::
   :widths: auto
   :header-rows: 0

   * - .. image:: docs/images/flag_yellow_low.jpg
          :width: 200px
     - This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 101028278.
