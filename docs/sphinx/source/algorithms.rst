Reusable algorithm layer (K)
============================

The K-layer contains reusable, dependency-light algorithms that are shared between:

* measures
* topo_groups
* polymer classification
* mapping/backmapping
* future analysis CLI tools

Key modules (v1)
---------------

* Graph views: edge lists and adjacency
* Components / clusters / percolation primitives
* Graph shortest paths
* Graph structural descriptors
* Contact graph builder (placeholder)
* Stats/numerics primitives
* Mapping/backmapping
* Polymer classifier
* Chain ordering utilities

Design contract
---------------

Algorithms depend only on:

* :c++:type:`GraphView`
* :c++:type:`GeometryView`
* indices (MoleculeIndex/ChainIndex)

They do not depend on Runner, Reader, or IMeasure lifecycles.
