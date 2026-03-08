Quickstart
==========

This page is the shortest path from a fresh checkout to a validated run.
If you are new to PILOTS, build the executable, validate a small config, then move to
:doc:`config`, :doc:`groups`, and :doc:`measures-reference`.

Build
-----

PILOTS uses CMake.

.. code-block:: bash

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j

The executable will be available as ``build/pilots``.

Run
---

PILOTS is configured with an INI file.

.. code-block:: bash

   ./build/pilots --config path/to/run.ini

Useful first-run commands:

* List the currently registered measure types:

  .. code-block:: bash

     ./build/pilots --list-measures

* Validate a configuration before processing the trajectory:

  .. code-block:: bash

     ./build/pilots --config run.ini --validate-config

* Run with an explicit thread count:

  .. code-block:: bash

     ./build/pilots --config run.ini --threads 8

View the manual
---------------

This repository contains a prebuilt manual under ``docs/manual/``.
Open ``docs/manual/index.html`` in a browser.

Rebuild the manual locally
--------------------------

To rebuild the HTML manual with Sphinx:

.. code-block:: bash

   python -m pip install -r docs/sphinx/requirements.txt
   sphinx-build -b html docs/sphinx/source docs/manual

Python integration
------------------

If you prefer a Python-first workflow for config generation and result loading, see :doc:`python`.
