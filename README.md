# PILOTS: Physical Integrated Library for Observables & Trajectories in Simulations

## Overview
PILOTS is a high-performance, reproducible trajectory-analysis runner for large-scale molecular simulations. The current v1 scope targets **LAMMPS** trajectories (text dump, with optional LAMMPS data topology).
The core design goal is simple: when a new physical observable is needed, it should usually only require adding a new *measure* implementation. Follow/flush/checkpoint, selection, topology/graph primitives, auditing, and `results.json` indexing are handled by the platform.

---
## Documentation

- [Online Manual](https://haoyuwu97.github.io/PILOTS/)
- [Quickstart](https://haoyuwu97.github.io/PILOTS/getting-started.html)
- [Configuration](https://haoyuwu97.github.io/PILOTS/config.html)
- [Measures Reference](https://haoyuwu97.github.io/PILOTS/measures-reference.html)
- [Adding a measure](https://haoyuwu97.github.io/PILOTS/adding-measures.html)

---
## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

---
## Run

```bash
./build/pilots --config path/to/config.ini --threads 8
```

**Useful CLI helpers**

- `--list-measures` prints the registered measure types.
- `--validate-config` parses and validates a configuration without processing the trajectory.

---
## Python Wrapper (`pilotsio`)

This repository also ships an optional pure-Python wrapper under `python/`. It provides a small config DSL, a `subprocess` launcher for the `pilots` binary, and a helper for reading `results.json` and text datasets.

From the repository root:

```bash
python -m pip install -e python
```

Optional DataFrame support:

```bash
python -m pip install -e "python[pandas]"
```

---
## Capabilities (v1)

- **Reader:** LAMMPS text dump (`xu/yu/zu` or `x/y/z + ix/iy/iz`), correct triclinic/PBC unwrapped reconstruction, robust follow mode.
- **Engineering loop:** `flush_every_frames/seconds` streaming output plus `checkpoint_out/resume_from` for restartable long jobs.
- **Topology:** optional loading of LAMMPS data `masses/bonds/angles/dihedrals/impropers`, audited in `results.json`.
- **Selection:** `[groups]` + `[topo_groups]` + boolean combination (`& | !`) with dynamic-selection auditing.
- **Algorithms (K layer):** reusable graph/cluster primitives, shortest paths, descriptors, statistics, mapping/backmapping, polymer classifier.
- **Measure framework:** runner + registry, multiple measures in a single trajectory pass.
- **Outputs:** text datasets plus `results.json` (`schema_version=1.0`) with atomic updates.

---
## File Descriptions

- **`src/`**  
  Measure implementations, registration units, and runner internals.

- **`include/pilots/`**  
  Public headers for core, SDK, algorithms, selection, topology, and correlators.

- **`python/`**  
  Optional pure-Python wrapper and convenience helpers.

- **`docs/`**  
  Long-form documentation sources for the online manual.

- **`tests/`**  
  Smoke-test inputs and regression fixtures.

- **`bench/`**  
  Benchmark-related inputs or performance utilities.

- **`CMakeLists.txt`**  
  CMake build definition for the project.

- **`LICENSE`**  
  The license text for the repository.

---
## Adding a Measure

Copy `src/measures/template_measure.cpp` to a new file under `src/measures/`, implement the physics, and keep the registration line in the same translation unit. No changes to `main` or `Runner` should be needed.

For implementation guidance, see the manual sections on the measure reference, SDK, and developer notes.

---
