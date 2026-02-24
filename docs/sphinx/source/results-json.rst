Results Index
=============

Alongside measure-specific text outputs, PILOTS writes a machine-readable
``results.json`` index that captures:

* input fingerprints (config/input/topology)
* model scale and mapping spec hash
* group and topo_group audits (including dynamic size stats)
* topology sections loaded
* measure metadata and profiling

This is intended to make experiments **reproducible and auditable**.

Schema version
--------------

The schema is versioned; see ``schema_version`` in the file.

Dynamic selections
------------------

When a selection depends on a dynamic group/topo_group, PILOTS records
``size_min/size_max/size_mean`` for the selection spec.

Audits and quality control (M0)
-------------------------------

PILOTS always writes a small set of *run-level audits* that make it easy to
sanity-check a run and keep results reproducible.

These audits are produced even if you only compute a single measure.

Run-level frame stats
~~~~~~~~~~~~~~~~~~~~~

``run.frame_stats`` is computed over frames that were actually processed by the
Runner (i.e. after ``resume_from``, it refers to the resumed segment).

Fields:

* ``frames``: number of processed frames.
* ``timestep_first/timestep_last``: first/last processed timestep.
* ``timestep_delta_mode``: the most common timestep stride between adjacent
  processed frames (often the dump interval).
* ``irregular_stride_count``: how many adjacent timestep deltas differ from the
  mode.
  This can indicate non-uniform output sampling or gaps, but it is *not* a
  guarantee of missing data.
* ``natoms_min/natoms_max``: a quick check for atom loss/growth.
* ``box_changes``: how many processed frames have a box different from the
  first processed frame (common for NPT).

Selections and topology audits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``groups.atom_groups`` and ``groups.topo_groups`` record each group expression
  and, for dynamic groups, ``size_min/size_max/size_mean``.
* ``groups.combined_selections`` records audit stats for combined selections
  (AtomGroup + TopoGroup + boolean combine).
* ``run.topology`` includes ``loaded_sections`` and the ``bond_graph`` definition
  used to build adjacency.

Mapping audit
~~~~~~~~~~~~~

When mapping is enabled, ``mapping.spec_hash_fnv1a64`` is recorded so that a
checkpoint/resume cannot silently change bead identity.
