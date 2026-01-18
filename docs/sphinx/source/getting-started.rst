Getting started
===============

Build
-----

PILOTS uses CMake.

.. code-block:: bash

   cmake -S . -B build
   cmake --build build -j

The executable will be available as ``build/pilots``.

Run
---

PILOTS is configured via an INI file:

.. code-block:: bash

   ./build/pilots --config path/to/run.ini

Common workflows:

* List available measure types:

  .. code-block:: bash

     ./build/pilots --list-measures

* Validate a configuration without processing the trajectory:

  .. code-block:: bash

     ./build/pilots --config run.ini --validate-config

View the manual
---------------

This repository contains a prebuilt manual under ``docs/manual/``.
Open ``docs/manual/index.html`` in a browser.

Long-term manual build (Sphinx)
-------------------------------

To rebuild the HTML manual with Sphinx:

.. code-block:: bash

   python -m pip install -r docs/sphinx/requirements.txt
   sphinx-build -b html docs/sphinx/source docs/manual
