Python integration (pilotsio)
============================

Overview
--------

PILOTS is implemented as a high-performance C++ runner (the ``pilots`` executable).
The optional **Python wrapper** (``pilotsio``) provides:

- a small configuration DSL to build PILOTS INI files from Python,
- a launcher that calls ``pilots`` via ``subprocess`` (no ABI coupling),
- a results helper that loads ``results.json`` and reads datasets.

This design keeps the scientific and auditing contract stable:
``results.json`` is the interface between C++ and Python.

Install / build the ``pilots`` executable
----------------------------------------

Linux
^^^^^

Install build dependencies and build from source:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install -y cmake g++ make

   # in the repository root
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ctest --test-dir build --output-on-failure

The executable will be at ``build/pilots``.

macOS
^^^^^

On macOS, OpenMP typically requires ``libomp``. Using Homebrew LLVM is often the simplest:

.. code-block:: bash

   brew install cmake llvm libomp

   export CC="$(brew --prefix llvm)/bin/clang"
   export CXX="$(brew --prefix llvm)/bin/clang++"

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ctest --test-dir build --output-on-failure

Windows
^^^^^^^

Recommended for scientific workflows: **WSL2** (Ubuntu) and follow the Linux steps.

Native Windows builds are also possible with Visual Studio (MSVC) and CMake:

.. code-block:: powershell

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --config Release
   ctest --test-dir build -C Release --output-on-failure

The executable is typically ``build\Release\pilots.exe``.

Install the Python wrapper (pilotsio)
-------------------------------------

``pilotsio`` is a pure-Python package shipped in this repository under ``python/``.

Create an environment (recommended):

.. code-block:: bash

   python --help  # ensure Python is available

   python -m venv .venv

   # Linux/macOS
   source .venv/bin/activate

   # Windows PowerShell
   # .venv\Scripts\Activate.ps1

Install in editable mode:

.. code-block:: bash

   python -m pip install -U pip
   python -m pip install -e python

Optional (dataset loading as pandas DataFrames):

.. code-block:: bash

   python -m pip install -e "python[pandas]"

Make ``pilots`` discoverable
----------------------------

``pilotsio`` locates the C++ executable using the following resolution order:

1. an explicit ``pilots_bin=...`` argument in Python,
2. the environment variable ``PILOTS_BIN``,
3. ``pilots`` / ``pilots.exe`` on ``PATH``,
4. local build fallbacks (``./build/pilots``).

Recommended:

.. code-block:: bash

   export PILOTS_BIN=/abs/path/to/build/pilots

Example: run from Python
------------------------

The Python DSL creates an INI file and stores it as ``<output_dir>/config_used.ini``
for reproducibility.

.. code-block:: python

   from pilotsio import PilotsConfig, run

   cfg = (
       PilotsConfig()
       .general(dt=1.0, output_dir="out")
       .group("all", "type:*")
       .measure(
           "msd",
           type="msd",
           correlator="exact",
           group="all",
           remove_drift=True,
           drift_group="all",
           output="msd_all.dat",
       )
   )

   res = run(cfg, dump="traj.dump", topology=None)
   df = res.load_dataset_dataframe("msd")
   print(df.head())

Validate-only (no side effects)
------------------------------

You can validate a configuration without producing outputs:

.. code-block:: python

   from pilotsio import PilotsConfig, validate

   cfg = PilotsConfig().general(dt=1.0, output_dir="out").measure("msd", type="msd")
   validate(cfg, dump="traj.dump")

This runs ``pilots --validate-config`` and should not create output files.

OpenMP threads
--------------

PILOTS uses OpenMP. You can control threads via:

- environment variable ``OMP_NUM_THREADS``
- or ``run(..., omp_threads=N)``

.. code-block:: python

   res = run(cfg, dump="traj.dump", omp_threads=16)
