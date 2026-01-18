from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from .config import PilotsConfig
from .results import PilotsResults


class PilotsRunError(RuntimeError):
    """Raised when the PILOTS binary exits with a non-zero return code."""

    def __init__(self, message: str, *, returncode: int, cmd: Sequence[str]) -> None:
        super().__init__(message)
        self.returncode = returncode
        self.cmd = list(cmd)


def _resolve_path(p: Union[str, Path], *, base: Optional[Path] = None) -> Path:
    pp = Path(p).expanduser()
    if not pp.is_absolute() and base is not None:
        pp = base / pp
    # resolve(strict=False) is available on Python>=3.9
    return pp.resolve(strict=False)


def find_pilots_binary(pilots_bin: Optional[Union[str, Path]] = None) -> Path:
    """Locate the `pilots` executable.

    Resolution order:

    1) Explicit `pilots_bin` argument
    2) Environment variable `PILOTS_BIN`
    3) `pilots` (or `pilots.exe`) on PATH
    4) Common local build locations relative to the current working directory

    Returns an absolute Path.
    """

    candidates: List[Union[str, Path]] = []
    if pilots_bin is not None:
        candidates.append(pilots_bin)

    env_bin = os.environ.get("PILOTS_BIN", "").strip()
    if env_bin:
        candidates.append(env_bin)

    # PATH search
    which = shutil.which("pilots") or shutil.which("pilots.exe")
    if which:
        candidates.append(which)

    # Local build fallbacks
    candidates.extend([
        Path("build") / "pilots",
        Path("build") / "Release" / "pilots.exe",
        Path("build") / "Debug" / "pilots.exe",
    ])

    for c in candidates:
        p = _resolve_path(c)
        if p.exists() and p.is_file():
            return p

    msg = (
        "Could not find the PILOTS executable.\n"
        "- Build it from source: cmake -S . -B build && cmake --build build\n"
        "- Then either add build/ to PATH or set PILOTS_BIN=/path/to/pilots\n"
    )
    raise FileNotFoundError(msg)


def _materialize_config(
    cfg: PilotsConfig,
    *,
    dump: Union[str, Path],
    output_dir: Union[str, Path],
    topology: Optional[Union[str, Path]] = None,
    base_dir: Optional[Path] = None,
) -> PilotsConfig:
    out = cfg.copy()

    # Required inputs
    dump_abs = _resolve_path(dump, base=base_dir).as_posix()
    out.set("general", "input", dump_abs)

    out_dir_abs = _resolve_path(output_dir, base=base_dir).as_posix()
    out.set("general", "output_dir", out_dir_abs)

    if topology is not None:
        topo_abs = _resolve_path(topology, base=base_dir).as_posix()
        out.set("topology", "file", topo_abs)

    # If mapping.file.path exists, normalize it to an absolute path for robustness.
    mp = out.get("mapping.file", "path")
    if mp:
        mp_abs = _resolve_path(mp, base=base_dir).as_posix()
        out.set("mapping.file", "path", mp_abs)

    # sanity: dt must be present for scientific correctness
    if not out.has("general", "dt"):
        raise ValueError("PILOTS requires [general] dt to be set (simulation time per frame).")

    return out


def validate(
    cfg: PilotsConfig,
    *,
    dump: Union[str, Path],
    output_dir: Union[str, Path] = "out",
    topology: Optional[Union[str, Path]] = None,
    pilots_bin: Optional[Union[str, Path]] = None,
    workdir: Optional[Union[str, Path]] = None,
    extra_args: Optional[Sequence[str]] = None,
    env: Optional[Dict[str, str]] = None,
) -> None:
    """Validate a configuration via `pilots --validate-config`.

    This should have no side effects (no output files, no directory creation),
    but it performs the same dependency checks as a real run.
    """

    _run_impl(
        cfg,
        dump=dump,
        output_dir=output_dir,
        topology=topology,
        pilots_bin=pilots_bin,
        workdir=workdir,
        validate_only=True,
        extra_args=extra_args,
        env=env,
    )


def run(
    cfg: PilotsConfig,
    *,
    dump: Union[str, Path],
    output_dir: Union[str, Path] = "out",
    topology: Optional[Union[str, Path]] = None,
    pilots_bin: Optional[Union[str, Path]] = None,
    workdir: Optional[Union[str, Path]] = None,
    validate_config: bool = True,
    omp_threads: Optional[int] = None,
    extra_args: Optional[Sequence[str]] = None,
    env: Optional[Dict[str, str]] = None,
) -> PilotsResults:
    """Run PILOTS and return a `PilotsResults` helper.

    Parameters
    ----------
    cfg:
        PilotsConfig built in Python.
    dump:
        Path to a LAMMPS text dump.
    output_dir:
        Where outputs (including results.json) will be written.
    topology:
        Optional topology file (e.g. LAMMPS data).
    validate_config:
        If True, run a fast `--validate-config` pass first.
    omp_threads:
        If provided, sets OMP_NUM_THREADS for the subprocess.
    """

    if validate_config:
        validate(
            cfg,
            dump=dump,
            output_dir=output_dir,
            topology=topology,
            pilots_bin=pilots_bin,
            workdir=workdir,
            extra_args=extra_args,
            env=env,
        )

    res = _run_impl(
        cfg,
        dump=dump,
        output_dir=output_dir,
        topology=topology,
        pilots_bin=pilots_bin,
        workdir=workdir,
        validate_only=False,
        extra_args=extra_args,
        env=env,
        omp_threads=omp_threads,
    )

    # output_dir is written to config as absolute; resolve the same way here.
    base = Path(workdir).resolve() if workdir else Path.cwd().resolve()
    out_dir_abs = _resolve_path(output_dir, base=base)
    return PilotsResults(out_dir_abs)


def _run_impl(
    cfg: PilotsConfig,
    *,
    dump: Union[str, Path],
    output_dir: Union[str, Path],
    topology: Optional[Union[str, Path]],
    pilots_bin: Optional[Union[str, Path]],
    workdir: Optional[Union[str, Path]],
    validate_only: bool,
    extra_args: Optional[Sequence[str]],
    env: Optional[Dict[str, str]],
    omp_threads: Optional[int] = None,
) -> subprocess.CompletedProcess[str]:
    pilots_path = find_pilots_binary(pilots_bin)

    base = Path(workdir).resolve() if workdir else Path.cwd().resolve()

    # Materialize a fully-specified config. We normalize key paths to absolute.
    cfg_used = _materialize_config(cfg, dump=dump, output_dir=output_dir, topology=topology, base_dir=base)

    # Write config:
    # - validate: temp file (avoid side effects in output_dir)
    # - run: store config in output_dir for reproducibility
    if validate_only:
        with tempfile.TemporaryDirectory(prefix="pilotsio-") as td:
            cfg_path = Path(td) / "config.ini"
            cfg_used.write(cfg_path)
            return _invoke_pilots(
                pilots_path,
                cfg_path,
                workdir=base,
                validate_only=True,
                extra_args=extra_args,
                env=env,
                omp_threads=omp_threads,
            )

    out_dir_abs = _resolve_path(output_dir, base=base)
    out_dir_abs.mkdir(parents=True, exist_ok=True)
    cfg_path = out_dir_abs / "config_used.ini"
    cfg_used.write(cfg_path)

    return _invoke_pilots(
        pilots_path,
        cfg_path,
        workdir=base,
        validate_only=False,
        extra_args=extra_args,
        env=env,
        omp_threads=omp_threads,
    )


def _invoke_pilots(
    pilots_path: Path,
    cfg_path: Path,
    *,
    workdir: Path,
    validate_only: bool,
    extra_args: Optional[Sequence[str]],
    env: Optional[Dict[str, str]],
    omp_threads: Optional[int],
) -> subprocess.CompletedProcess[str]:
    cmd: List[str] = [str(pilots_path), "--config", str(cfg_path)]
    if validate_only:
        cmd.append("--validate-config")
    if extra_args:
        cmd.extend(list(extra_args))

    proc_env = dict(os.environ)
    if env:
        proc_env.update(env)
    if omp_threads is not None:
        proc_env["OMP_NUM_THREADS"] = str(int(omp_threads))

    cp = subprocess.run(cmd, cwd=str(workdir), env=proc_env, text=True)
    if cp.returncode != 0:
        raise PilotsRunError(
            f"PILOTS failed with return code {cp.returncode} (config={cfg_path}).",
            returncode=cp.returncode,
            cmd=cmd,
        )
    return cp
