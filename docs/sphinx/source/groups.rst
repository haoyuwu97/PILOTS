Selections
==========

PILOTS supports two selection sources:

* **Atom groups** from ``[groups]`` based on dump fields, regions, id/type/mol, and similar predicates.
* **Topology groups** from ``[topo_groups]`` based on bonds, connectivity, components, and related topology-derived queries.

Both are converted into a per-frame atom set.
Measures can combine them with boolean logic at measure scope.

Atom groups
-----------

Atom groups are defined under ``[groups]``. Example:

.. code-block:: ini

   [groups]
   backbone = type==1
   slab = region:slab zmin=0 zmax=10

Topology groups
---------------

Topology groups are defined under ``[topo_groups]`` and always produce an atom set.
They require topology bonds to be loaded.

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

The complement operators ``!A`` and ``!T`` are defined with respect to the current frame's atom universe.
Many time-correlation measures additionally require the resulting selection to be **static** across frames.
