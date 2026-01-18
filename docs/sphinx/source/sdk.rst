Measure SDK (L)
===============

The SDK layer wraps the K-layer algorithms and runner plumbing with:

* strongly-typed handles (fields, selections, topology, mapping)
* output builders that automatically populate ``results.json`` descriptors
* stable state IO helpers (checkpoint/resume safety)

This keeps the physical definition of a measure close to the formula, while the platform enforces correctness.

Typical structure
-----------------

PILOTS recommends a *one-translation-unit* measure layout:

.. code-block:: c++

   // src/measures/MyMeasure.cpp
   // 1) caps()
   // 2) create()
   // 3) class MyMeasure : IMeasure
   // 4) static MeasureRegistrar ...

See ``src/measures/template_measure.cpp``.
