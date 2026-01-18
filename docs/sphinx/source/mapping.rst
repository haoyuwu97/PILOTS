Mapping and backmapping
=======================

Mapping converts atom-level trajectories into **bead-level** trajectories
for AA/UA/complex CG workflows.

The mapping subsystem is part of the reusable algorithm layer (K6) and is
audited and checkpointed via a stable *spec hash*.

Model scale
-----------

In the INI file:

.. code-block:: ini

   [model]
   scale = aa    ; aa|ua|cg|auto

Mapping modes
-------------

.. code-block:: ini

   [mapping]
   mode = identity   ; none|identity|by_mol|file
   position = com_mass  ; com_mass|com_geom|representative_atom

The meaning of each mode:

``identity``
   Beads == atoms.

``by_mol``
   One bead per molecule id (requires ``mol`` in dump).

``file``
   External mapping file defines atom_id -> bead_id (and optional weights).

Backmapping
-----------

Backmapping is provided as a utility to scatter bead-level results back to
per-atom arrays for output or visualization.
