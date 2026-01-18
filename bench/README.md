# PILOTS benchmarking

This folder contains a small, maintainable Python benchmark driver.

## Build first

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

## Run benchmark

From the project root:

```bash
python3 bench/run_benchmark.py --pilots ./build/pilots --config tests/configs/msd.ini --threads 1,2,4,8 --runs 3
```

It writes:
- `benchmark.csv` (default)
- per-run output directories under `bench/_runs/`

The script reads each run's `results.json` to extract wall time and breakdowns.
