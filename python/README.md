# pilotsio

`pilotsio` is a pure-Python wrapper for the **PILOTS** trajectory analysis runner.

It provides:

- A small **configuration DSL** to build PILOTS INI configs from Python.
- A `run()` helper that executes the `pilots` binary via `subprocess`.
- A results helper that loads `results.json` and reads output datasets.

## Important

`pilotsio` does **not** compile or bundle the `pilots` C++ binary.
You must install/build `pilots` separately and make it discoverable via:

- `PILOTS_BIN=/path/to/pilots` (recommended), or
- putting `pilots` (or `pilots.exe`) on your `PATH`.

## Installation (from the PILOTS repository)

From the repository root:

```bash
python -m pip install -e python
```

Optional (for DataFrame loading):

```bash
python -m pip install -e "python[pandas]"
```

See the HTML manual for end-to-end installation and usage.
