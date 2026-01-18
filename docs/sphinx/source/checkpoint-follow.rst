Follow, flush, checkpoint, resume
=================================

PILOTS is designed to run robustly on long HPC jobs.

Follow mode
-----------

If ``general.follow=true``, the LAMMPS dump reader waits for additional data
instead of failing at EOF. This supports real-time analysis while LAMMPS is still writing.

Flush
-----

Results are periodically flushed to output files and ``results.json`` using atomic write patterns.

Checkpoint/resume
-----------------

Use ``general.checkpoint_out`` to write state and ``general.resume_from`` to continue a killed job.
