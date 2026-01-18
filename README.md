# PILOTS

PILOTS is a high-performance, reproducible trajectory analysis runner for large-scale molecular simulations.
The current v1 scope targets **LAMMPS** (text dump + optional LAMMPS data topology).

**Design goal:** when adding a new physical observable, you primarily write a new *measure* implementation.
Follow/flush/checkpoint, selection, topology/graph primitives, auditing, and `results.json` indexing are provided by the platform.

## Documentation

The manual is written in Sphinx and can be built into a static HTML site.

- **Sphinx sources:** `docs/sphinx/source/`
- **Build locally:** `./docs/build_docs.sh`
- **HTML output directory:** `docs/manual/` (generated)

### Online manual (GitHub Pages)

This repository includes an automated GitHub Pages deployment workflow:

- `.github/workflows/docs-pages.yml`

One-time setup:

1. Go to **Settings -> Pages**.
2. Under **Build and deployment**, set **Source** to **GitHub Actions**.

After that, editing any `docs/sphinx/source/*.rst` file directly on GitHub (web UI) and committing to `main`
will automatically rebuild and deploy the manual.

The published URL will be:

- https://OWNER.github.io/REPO/

Replace `OWNER` and `REPO` with your GitHub username/organization and repository name.

## Build

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

## Run

```bash
./build/pilots --config path/to/config.ini --threads 8
```

Useful CLI helpers:

- `--list-measures` (print available measure types)
- `--validate-config` (parse + validate configuration and exit)

## Python wrapper (pilotsio)

This repository also ships an optional pure-Python wrapper under `python/`.
It provides a small config DSL, a `subprocess` launcher for the `pilots` binary,
and a helper for reading `results.json` + datasets.

From the repository root:

```bash
python -m pip install -e python
```

Optional (DataFrame loading):

```bash
python -m pip install -e "python[pandas]"
```

See the manual for cross-platform build/install instructions.

## Capabilities (v1)

- **Reader:** LAMMPS text dump (`xu/yu/zu` or `x/y/z + ix/iy/iz`), correct triclinic/PBC unwrapped reconstruction, robust follow mode.
- **Engineering loop:** `flush_every_frames/seconds` streaming output + `checkpoint_out/resume_from` for restartable long jobs.
- **Topology:** LAMMPS data optional loading of `masses/bonds/angles/dihedrals/impropers`, audited in `results.json`.
- **Selection:** `[groups]` + `[topo_groups]` + boolean combination (`& | !`) with dynamic-selection auditing.
- **Algorithms (K layer):** reusable graph/cluster primitives, shortest paths, descriptors, statistics, mapping/backmapping, polymer classifier.
- **Measure framework:** Runner + registry, multiple measures in a single trajectory pass.
- **Outputs:** text outputs (kept) + `results.json` (schema_version=1.0) atomic updates.

## Developing: add a new measure

Copy `src/measures/template_measure.cpp` to a new file under `src/measures/`, implement the physics, and enable the registration line
at the bottom (same translation unit). No changes to `main` or `Runner` are needed.

See the manual sections:

- "Measure SDK" and "Developer notes".
