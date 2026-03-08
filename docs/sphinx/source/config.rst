Configuration
=============

PILOTS is configured with an INI file.
A typical layout is:

.. code-block:: ini

   [general]
   input = dump.lammpstrj
   dt = 0.01
   output_dir = ./out

   [groups]
   ...

   [measure.msd]
   type = msd
   enabled = true
   group = all

[general]
---------

``input``
   Path to the LAMMPS text dump file.

``dt``
   Time step in simulation units. Used to convert frame lags into physical time.

``output_dir``
   Directory for text outputs, ``results.json``, and checkpoints.

``follow``
   If true, PILOTS blocks at EOF and waits for new frames appended to the dump file.

``flush_every_frames`` / ``flush_every_seconds``
   Periodic partial output and ``results.json`` updates during long runs.

``checkpoint_out`` / ``resume_from``
   Binary checkpointing and resume.

Measures
--------

Each measure is configured in a section named ``[measure.<instance>]``.
At minimum, define:

``type``
   The registered measure type, for example ``msd``.

``enabled``
   Whether the measure is active.

Common measure-level keys include ``group``, ``topo_group``, ``combine``, ``output``,
``frame_start``, and ``frame_end``.

See also: :doc:`groups`, :doc:`topology`, :doc:`mapping`, and :doc:`measures-reference`.
