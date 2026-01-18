from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Optional


def _as_str(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, Path):
        return v.as_posix()
    return str(v)


def _csv(items: Iterable[Any]) -> str:
    return ",".join(_as_str(x) for x in items)


@dataclass
class PilotsConfig:
    """Build a PILOTS INI configuration from Python.

    This is a thin, explicit builder. It does not attempt to validate every key,
    but it provides convenience methods for the most common sections.

    The underlying representation is an ordered mapping:

        section -> (key -> value)

    Values are rendered as raw strings in INI output.
    """

    _sections: "OrderedDict[str, OrderedDict[str, str]]" = None  # type: ignore

    def __post_init__(self) -> None:
        if self._sections is None:
            self._sections = OrderedDict()

    # ---------- Low-level API ----------
    def set(self, section: str, key: str, value: Any) -> "PilotsConfig":
        sec = self._sections.setdefault(section, OrderedDict())
        sec[key] = _as_str(value)
        return self

    def set_bool(self, section: str, key: str, value: bool) -> "PilotsConfig":
        return self.set(section, key, value)

    def set_list(self, section: str, key: str, values: Iterable[Any]) -> "PilotsConfig":
        return self.set(section, key, _csv(values))

    def update(self, section: str, mapping: Mapping[str, Any]) -> "PilotsConfig":
        for k, v in mapping.items():
            self.set(section, k, v)
        return self

    def has(self, section: str, key: Optional[str] = None) -> bool:
        if section not in self._sections:
            return False
        if key is None:
            return True
        return key in self._sections[section]

    def get(self, section: str, key: str, default: Optional[str] = None) -> Optional[str]:
        if section not in self._sections:
            return default
        return self._sections[section].get(key, default)

    def section(self, section: str) -> "OrderedDict[str, str]":
        return self._sections.setdefault(section, OrderedDict())

    def copy(self) -> "PilotsConfig":
        out = PilotsConfig()
        for sec, kv in self._sections.items():
            out._sections[sec] = OrderedDict(kv)
        return out

    # ---------- Convenience builders ----------
    def general(self, **kwargs: Any) -> "PilotsConfig":
        return self.update("general", kwargs)

    def model(self, **kwargs: Any) -> "PilotsConfig":
        return self.update("model", kwargs)

    def mapping(self, **kwargs: Any) -> "PilotsConfig":
        return self.update("mapping", kwargs)

    def mapping_file(self, **kwargs: Any) -> "PilotsConfig":
        return self.update("mapping.file", kwargs)

    def topology(self, **kwargs: Any) -> "PilotsConfig":
        return self.update("topology", kwargs)

    def group(self, name: str, expr: str) -> "PilotsConfig":
        return self.set("groups", name, expr)

    def topo_group(self, name: str, expr: str) -> "PilotsConfig":
        return self.set("topo_groups", name, expr)

    def measure(self, instance: str, *, type: Optional[str] = None, enabled: bool = True, **kwargs: Any) -> "PilotsConfig":
        sec = f"measure.{instance}"
        self.set(sec, "enabled", enabled)
        if type is not None:
            self.set(sec, "type", type)
        for k, v in kwargs.items():
            self.set(sec, k, v)
        return self

    # ---------- Rendering ----------
    def to_ini(self) -> str:
        """Render the configuration as an INI string.

        Ordering policy:

        - Core sections first (general/model/mapping/mapping.file/topology)
        - Selection sections (groups/topo_groups)
        - measure.* sections sorted by instance name
        - Other sections in insertion order
        """

        def write_section(lines: list[str], sec: str, kv: MutableMapping[str, str]) -> None:
            lines.append(f"[{sec}]")
            for k, v in kv.items():
                lines.append(f"{k} = {v}")
            lines.append("")

        lines: list[str] = []

        priority = [
            "general",
            "model",
            "mapping",
            "mapping.file",
            "topology",
            "groups",
            "topo_groups",
        ]

        seen: set[str] = set()

        for sec in priority:
            if sec in self._sections:
                write_section(lines, sec, self._sections[sec])
                seen.add(sec)

        measure_secs = sorted([s for s in self._sections.keys() if s.startswith("measure.")])
        for sec in measure_secs:
            write_section(lines, sec, self._sections[sec])
            seen.add(sec)

        for sec, kv in self._sections.items():
            if sec in seen:
                continue
            write_section(lines, sec, kv)

        return "\n".join(lines).rstrip() + "\n"

    def write(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        text = self.to_ini()
        path.write_text(text, encoding="utf-8")
        return path
