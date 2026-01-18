"""pilotsio: Pure-Python wrapper for the PILOTS runner.

This package intentionally keeps a strict separation:

- `pilots` (C++): compute engine (CLI executable)
- `pilotsio` (Python): configuration DSL + subprocess launcher + results loader

The stable machine-readable contract is `results.json`.
"""

from .config import PilotsConfig
from .runner import run, validate, find_pilots_binary
from .results import PilotsResults

__all__ = [
    "PilotsConfig",
    "PilotsResults",
    "run",
    "validate",
    "find_pilots_binary",
]
