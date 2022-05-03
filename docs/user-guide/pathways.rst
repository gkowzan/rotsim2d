.. currentmodule:: rotsim2d.pathways

#########################################
Creating and manipulating excitation tree
#########################################

First-order excitation
======================

Create :class:`KetBra` in ground vibrational state

.. ipython:: python
   
   import rotsim2d.pathways as pw
   import molspecutils.molecule as mol
   kb1 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
   print(kb1)

and excite it,

.. ipython:: python

   pw.excite(kb1, light_name='omg1')
   kb1.print_tree()
   
`kb1` is modified in-place. Use :func:`readout` to simulate measurement and see which
branches will contribute to macroscopic polarization, i.e. would survive
another one-sided application of the dipole operator and taking the trace of the
density matrix.

.. ipython:: python

   pw.readout(kb1)
   kb1.print_tree()


**********************
Third-order excitation
**********************

Create :class:`KetBra` in ground vibrational state and perform three excitations
with :func:`multi_excite`:

.. ipython:: python

   kb2 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
   kb2 = pw.multi_excite(kb2, light_names=['omg1', 'omg2', 'omg3'],
                         parts=['ket', 'both', 'both'])
   kb2.print_tree()

The excitation tree in this case is very long.  We can save it's graphical
representation with :meth:`KetBra.savepng`. This uses
:class:`anytree.exporter.DotExporter` and requires a working installation of
`GraphViz <https://graphviz.org/>`_.

.. ipython:: python

   kb2.savepng('docs/images/kb2.png')

.. image:: /images/kb2.png
   :target: ../_images/kb2.png

We can also generate double-sided Feynmann diagrams corresponding to the tree with LaTeX/Tikz:

.. ipython:: python

   import rotsim2d.visual as vis
   latex_code = vis.tikz_diagrams(kb2)
   #vis.latex_compile('images/kb2_tikz.tex', vis.latex_document(latex_code))

See the result: :download:`kb2_tikz.pdf <../images/kb2_tikz.pdf>`. See
:func:`rotsim2d.visual.tikz_diagrams` for more details.

******************
Filtering the tree
******************

`kb2` contains all rephasing and non-rephasing pathways, emitted in all
directions.  The tree can be filtered with one of :mod:`rotsim2d.pathways`
functions starting with `remove_` or only `only_`, see
:ref:`tree-filtering-functions`.  For example, to obtain only pathways emitted
in :math:`\vec{k}_s = \vec{k}_1 - \vec{k_2} + \vec{k}_3` (non-rephasing
pathways) we can use:

.. ipython:: python

   kb2_SII = pw.only_SII(kb2.copy())
   kb2_SII = pw.readout(kb2_SII)
   kb2_SII.savepng('docs/images/kb2_SII.png')

.. image:: /images/kb2_SII.png
   :target: ../_images/kb2_SII.png

For :math:`\vec{k}_s = \vec{k}_1 + \vec{k_2} - \vec{k}_3` direction (rephasing
without double-quantum coherences), we can use:

.. ipython:: python

   kb2_SI = pw.only_SI(kb2.copy())
   kb2_SI = pw.readout(kb2_SI)
   kb2_SI.savepng('docs/images/kb2_SI.png')

.. image:: /images/kb2_SI.png
   :target: ../_images/kb2_SI.png
   
*************************
Multiple excitations tree
*************************

Excitation trees for many initial `j` values can be generated with :func:`gen_pathways`:

.. ipython:: python

   from pprint import pprint
   js = (0, 1, 2, 3)
   kbs = pw.gen_pathways(js, meths=[pw.only_SII, pw.only_twocolor],
                         rotor='symmetric', kiter_func="range(j+1)")
   pprint(kbs)

This call created a list of third-order excitation trees for `j` values in
`js`. Each tree was filtered by :func:`only_SII` and
:func:`only_twocolor`. `meths` argument can contain any callable that takes
:class:`KetBra` instance as its sole argument and returns it, presumably after
filtering it in some way, see :ref:`tree-filtering-functions`.
`rotor='symmetric'` causes the functions to generate
tree for different `k` quantum number values. The range of these values is
determined by `kiter_func` [#f1]_. `kiter_func` is a Python expression evaluated
with `ASTEVAL <https://newville.github.io/asteval/>`_, which should return an
iterable over integers. The expression is evaluated for each `j` value from `js`
argument. The default expression is ``"range(j+1)"``, which iterates from 0 to `j`.

.. rubric:: Footnotes

.. [#f1] This argument could be called more aptly `kiter_expr`.
