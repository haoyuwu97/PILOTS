results.json
============

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
