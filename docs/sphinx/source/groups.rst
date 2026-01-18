Groups and selections
=====================

PILOTS supports two selection inputs:

* **Atom groups**: ``[groups]`` (based on dump fields, regions, id/type/mol, ...)
* **Topology groups**: ``[topo_groups]`` (based on topology bonds, components, ...)

Both are converted into a per-frame **atom set**, and measures can combine them
using boolean logic.

Atom groups
-----------

Atom groups are defined under ``[groups]``. Example:

.. code-block:: ini

   [groups]
   backbone = type==1
   slab = region:slab zmin=0 zmax=10

Topology groups
--------------

Topology groups are defined under ``[topo_groups]`` and always produce an atom set.
These require that topology bonds are loaded.

.. code-block:: ini

   [topo_groups]
   bonded = atoms_in_any_bond
   big_components = atoms_in_component_size>=100

Measure-level combination
-------------------------

Each measure may specify:

* ``group``: atom-group expression (default: ``all``)
* ``topo_group``: topo-group expression (default: ``all``)
* ``combine``: boolean expression over ``A`` and ``T`` (default: ``A&T``)

The complement operators ``!A`` and ``!T`` are defined with respect to the current
frame's atom universe.
