PILOTS Manual
=============

PILOTS is a trajectory analysis runner designed for large-scale molecular simulations
with a focus on:

* **Scientific correctness** (fail-fast on missing fields, auditable selections)
* **Single-pass pipeline** (multiple measures in one trajectory pass)
* **Engineering reliability** (follow + flush + checkpoint/resume)

The manual is intended to read like an LAMMPS-style reference: you should be able
to navigate from CLI/configuration to algorithmic primitives and then to measure
authoring.

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

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide

   architecture
   algorithms
   sdk
   developer
