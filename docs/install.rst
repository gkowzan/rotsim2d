Installation
============

**Dependencies**

- numpy, scipy, anytree, sympy
- `gkpywigxjpf <https://gitlab.com/allisonlab/mdcs/gkpywigxjpf>`_
- `knickknacks <https://gitlab.com/allisonlab/mdcs/knickknacks>`_
- `molspecutils <https://gitlab.com/allisonlab/mdcs/molspecutils>`_

**Install**

.. highlight:: sh

Install the package by downloading or cloning the repo and calling the following
inside main repository directory (containing `setup.py`)::

  python -m pip install -e .

or by installing directly from the repo with pip::

  python -m pip install git+ssh://git@gitlab.com/allisonlab/mdcs/rotsim2d.git@master

`molspecutils <https://gitlab.com/allisonlab/mdcs/molspecutils>`_ package needs to
be initialized by executing `molspecutils_init` script to
download spectroscopic data and fill the local database.
