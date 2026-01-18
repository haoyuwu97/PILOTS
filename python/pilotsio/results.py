from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union


def _resolve_output_path(output_dir: Path, p: str) -> Path:
    pp = Path(p)
    return pp if pp.is_absolute() else (output_dir / pp).resolve(strict=False)


@dataclass
class PilotsResults:
    """Helper for reading PILOTS outputs.

    The stable contract for downstream tooling is `results.json`.
    This helper loads the index and provides convenience accessors.
    """

    output_dir: Path
    results_json_name: str = "results.json"

    _index: Optional[Dict[str, Any]] = None

    def index_path(self) -> Path:
        p = Path(self.results_json_name)
        return p if p.is_absolute() else (self.output_dir / p)

    def load_index(self, *, force: bool = False) -> Dict[str, Any]:
        if self._index is None or force:
            path = self.index_path()
            self._index = json.loads(path.read_text(encoding="utf-8"))
        return self._index

    def schema_version(self) -> str:
        idx = self.load_index()
        return str(idx.get("schema_version", ""))

    def measures(self) -> List[Dict[str, Any]]:
        idx = self.load_index()
        ms = idx.get("measures", [])
        if not isinstance(ms, list):
            return []
        return ms

    def list_measure_instances(self) -> List[str]:
        return [str(m.get("instance", "")) for m in self.measures()]

    def find_measure(self, name_or_type: str) -> Dict[str, Any]:
        matches: List[Dict[str, Any]] = []
        for m in self.measures():
            if m.get("instance") == name_or_type or m.get("type") == name_or_type:
                matches.append(m)

        if not matches:
            raise KeyError(f"Measure not found: {name_or_type}")
        if len(matches) > 1:
            names = [f"{mm.get('instance')} (type={mm.get('type')})" for mm in matches]
            raise KeyError(
                "Ambiguous measure selector. Provide an instance name. Candidates: " + ", ".join(names)
            )
        return matches[0]

    def dataset_path(self, measure: Union[str, Dict[str, Any]], output_index: int = 0) -> Path:
        m = self.find_measure(measure) if isinstance(measure, str) else measure
        outputs = m.get("outputs", [])
        if not outputs:
            raise KeyError(f"Measure '{m.get('instance')}' has no outputs")
        if output_index < 0 or output_index >= len(outputs):
            raise IndexError(f"output_index out of range: {output_index}")
        p = str(outputs[output_index].get("path", ""))
        if not p:
            raise KeyError(f"Measure '{m.get('instance')}' output has empty path")
        return _resolve_output_path(self.output_dir, p)

    def dataset_columns(self, measure: Union[str, Dict[str, Any]], output_index: int = 0) -> List[str]:
        m = self.find_measure(measure) if isinstance(measure, str) else measure
        outputs = m.get("outputs", [])
        cols = outputs[output_index].get("columns", []) if outputs else []
        return [str(c) for c in cols]

    def load_dataset_table(self, measure: Union[str, Dict[str, Any]], output_index: int = 0) -> List[List[float]]:
        """Load a dataset as a numeric table (list of rows).

        This uses a simple whitespace split and ignores comment lines starting with '#'.
        """

        path = self.dataset_path(measure, output_index=output_index)
        rows: List[List[float]] = []
        with path.open("r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                rows.append([float(x) for x in s.split()])
        return rows

    def load_dataset_dataframe(self, measure: Union[str, Dict[str, Any]], output_index: int = 0):
        """Load a dataset into a pandas DataFrame.

        Requires `pandas`.
        """

        try:
            import pandas as pd  # type: ignore
        except Exception as e:
            raise ImportError(
                "pandas is required for load_dataset_dataframe(). Install with: pip install 'pilotsio[pandas]'"
            ) from e

        path = self.dataset_path(measure, output_index=output_index)
        cols = self.dataset_columns(measure, output_index=output_index)

        df = pd.read_csv(
            path,
            delim_whitespace=True,
            comment="#",
            header=None,
            names=cols if cols else None,
        )
        return df
