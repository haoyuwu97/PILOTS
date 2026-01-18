Command line interface
======================

The ``pilots`` executable supports a small, stable CLI.

.. code-block:: bash

   ./build/pilots --help

Options
-------

``--config <path>``
   Path to the INI configuration file.

``--threads <N>``
   Set the OpenMP thread count (if PILOTS is built with OpenMP).

``--list-measures``
   Print all registered measure types (one per line) and exit.

``--validate-config``
   Run the full platform-level validation (read first frame, build groups,
   build topology, instantiate measures) and exit without processing the full
   trajectory.

Exit codes
----------

* ``0`` success
* non-zero indicates a validation or runtime error
