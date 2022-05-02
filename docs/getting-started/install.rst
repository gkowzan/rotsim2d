Installation
============

**Dependencies**

- numpy, scipy, anytree, sympy
- `gkpywigxjpf <https://gitlab.com/allisonlab/mdcs/pywigxjpf>`_
- `knickknacks <https://gitlab.com/allisonlab/mdcs/shed>`_
- `molspecutils <https://gitlab.com/allisonlab/mdcs/spectroscopy>`_

**Install**

.. highlight:: sh

Install the package from our private GitLab repository by executing::

  pip install rotsim2d --extra-index-url https://<token_name>:<token>@gitlab.com/api/v4/projects/26140156/packages/pypi

where `<token_name>` and `<token>` are obtained by creating a personal token
with `read_api`` scope. See `Creating personal access tokens
<https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#creating-a-personal-access-token>`_
for details.

Molecular data
++++++++++++++

`rotsim2d` needs data on rovibrational states and transitions to produce 2D peak
plots, simulate 2D lineshapes or time-dependent signals. For simulations of
transitions within an isolated vibrational mode,
:class:`rotsim2d.dressedleaf.DressedPathway` needs to be provided with an object
implementing `VibrationalMode` interface defined in
:mod:`molspecutils.molecule`. The `molspecutils
<https://gitlab.com/allisonlab/mdcs/spectroscopy>`_ package currently contains
two classes implementing the interface: `COAlchemyMode` and `CH3ClAlchemyMode`,
where the latter is for the :math:`\nu_3` mode. The molecular data is not
bundled with the package and during the first run :mod:`molspecutils` will use
HAPI [1]_ to fetch transitions line lists from HITRAN [2]_ and cache them for
future use.

.. [1] R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, "HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data", J. Quant. Spectrosc. Radiat. Transfer *177*, 15-30 (2016).
.. [2] I.E. Gordon *et al.* "The HITRAN2016 molecular spectroscopic database", J. Quant. Spectrosc. Radiat. Transfer *203*, 3-69 (2017).

