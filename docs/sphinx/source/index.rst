PILOTS Manual
=============

PILOTS is a trajectory-analysis runner for large-scale molecular simulations with three priorities:

* **Scientific correctness** — fail fast on missing fields, keep selections auditable, and preserve identity consistency when a measure requires it.
* **Single-pass analysis** — run multiple measures in one trajectory pass.
* **Engineering reliability** — support follow mode, periodic flushing, and checkpoint/resume for long jobs.

For most users, the shortest path through the manual is:

1. :doc:`getting-started`
2. :doc:`config`
3. :doc:`groups`
4. :doc:`measures-reference`

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   getting-started
   cli
   config
   groups
   topology
   mapping
   checkpoint-follow
   results-json
   python

.. toctree::
   :maxdepth: 2
   :caption: Measure Reference

   measures-reference

.. toctree::
   :maxdepth: 2
   :caption: Internals and Extension

   architecture
   algorithms
   sdk
   adding-measures
   developer
