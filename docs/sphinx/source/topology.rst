Topology
========

PILOTS can optionally read LAMMPS ``data`` topology information and inject it
into measures via the :c++:type:`SystemContext`.

Loaded topology sections are **audited** in ``results.json``.

Inputs
------

* ``topology.format = lammps_data``
* ``topology.path = system.data``
* ``topology.load = masses,bonds,angles,...`` (mask; only parse what is needed)

Sections
--------

The topology container can hold:

* Masses (by type)
* Bonds
* Angles
* Dihedrals
* Impropers

Mol-first
---------

Measures that aggregate by molecule should prefer ``mol`` from the dump.
Optionally, the runner can derive molecule IDs from bond connectivity
(disabled by default to avoid implicit semantics).
