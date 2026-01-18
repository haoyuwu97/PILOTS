Architecture overview
=====================

PILOTS is structured as:

* **IO layer**: readers produce canonical :c++:type:`Frame` objects
* **Selection layer**: AtomGroup + TopoGroup -> combined atom selections
* **Topology layer**: optional; parsed once and injected via :c++:type:`SystemContext`
* **Algorithm layer (K)**: reusable graph/index/geometry + stats + mapping + polymer classifier
* **SDK layer (L)**: handles/views + output+state helpers for measures
* **Measures**: physics definitions that plug into the Runner pipeline

The key goal is: adding a measure should not require changing Runner or main.