Configuration file
==================

PILOTS is configured via an INI file. The typical structure is:

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
   Time step in your simulation units. Used to convert frame lags into physical time.

``output_dir``
   Directory for outputs (text datasets, results.json, checkpoints).

``follow``
   If true, PILOTS will block at EOF and wait for new frames appended to the dump file.

``flush_every_frames`` / ``flush_every_seconds``
   Periodic partial output and results.json updates during long runs.

``checkpoint_out`` / ``resume_from``
   Binary state checkpointing and resume.

Measures
--------

Each measure is configured in a section named ``[measure.<instance>]``.
At minimum:

``type``
   The registered type, e.g. ``msd``.

``enabled``
   Whether to run the measure.

See also: :doc:`groups`, :doc:`topology`, :doc:`mapping`.
