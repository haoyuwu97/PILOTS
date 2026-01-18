#!/usr/bin/env python3

import argparse
import configparser
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def parse_threads(s: str):
    vals = []
    for part in s.split(','):
        part = part.strip()
        if not part:
            continue
        v = int(part)
        if v <= 0:
            raise ValueError("thread count must be positive")
        vals.append(v)
    if not vals:
        raise ValueError("no threads specified")
    return vals


def load_ini(path: Path) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser()
    cfg.optionxform = str  # preserve case
    with path.open('r', encoding='utf-8') as f:
        cfg.read_file(f)
    return cfg


def write_ini(cfg: configparser.ConfigParser, path: Path) -> None:
    with path.open('w', encoding='utf-8') as f:
        cfg.write(f)


def main() -> int:
    ap = argparse.ArgumentParser(description="Benchmark PILOTS using results.json")
    ap.add_argument("--pilots", required=True, help="Path to pilots executable")
    ap.add_argument("--config", required=True, help="Path to base INI config")
    ap.add_argument("--threads", default="1,2,4,8", help="Comma-separated thread counts")
    ap.add_argument("--runs", type=int, default=3, help="Repeats per thread count")
    ap.add_argument("--out", default="benchmark.csv", help="CSV output file")
    ap.add_argument("--keep", action="store_true", help="Keep per-run output dirs")
    args = ap.parse_args()

    pilots = Path(args.pilots)
    base_cfg_path = Path(args.config)
    out_csv = Path(args.out)

    if not pilots.exists():
        print(f"ERROR: pilots executable not found: {pilots}", file=sys.stderr)
        return 2
    if not base_cfg_path.exists():
        print(f"ERROR: config not found: {base_cfg_path}", file=sys.stderr)
        return 2

    threads_list = parse_threads(args.threads)
    if args.runs <= 0:
        print("ERROR: --runs must be positive", file=sys.stderr)
        return 2

    runs_root = Path("bench/_runs")
    runs_root.mkdir(parents=True, exist_ok=True)

    rows = []
    for th in threads_list:
        for r in range(args.runs):
            run_id = f"t{th}_r{r}_{int(time.time()*1000)}"
            run_dir = runs_root / run_id
            run_dir.mkdir(parents=True, exist_ok=True)

            cfg = load_ini(base_cfg_path)
            if "general" not in cfg:
                cfg["general"] = {}

            # Ensure per-run outputs do not clobber each other.
            cfg["general"]["output_dir"] = str(run_dir)
            cfg["general"]["results_json"] = "results.json"
            cfg["general"]["profile"] = "false"  # keep profiling in results.json, but avoid noisy stderr

            tmp_cfg_path = run_dir / "config.ini"
            write_ini(cfg, tmp_cfg_path)

            cmd = [str(pilots), "--config", str(tmp_cfg_path), "--threads", str(th)]
            t0 = time.time()
            p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            wall = time.time() - t0
            if p.returncode != 0:
                print("\n".join(["ERROR: pilots failed", f"  cmd: {' '.join(cmd)}", f"  rc: {p.returncode}", "  stderr:", p.stderr]), file=sys.stderr)
                return 1

            res_path = run_dir / "results.json"
            if not res_path.exists():
                print(f"ERROR: results.json not found for run: {run_dir}", file=sys.stderr)
                return 1
            with res_path.open('r', encoding='utf-8') as f:
                res = json.load(f)

            run = res.get("run", {})
            prof = res.get("profiling", {})
            rows.append({
                "threads": th,
                "run": r,
                "wall_seconds_external": wall,
                "wall_seconds_pilots": run.get("wall_seconds", None),
                "reader_seconds": prof.get("reader_seconds", None),
                "group_setup_seconds": prof.get("group_setup_seconds", None),
                "run_dir": str(run_dir),
            })

            if not args.keep:
                # Keep results.json and outputs if desired, otherwise delete run directory.
                shutil.rmtree(run_dir, ignore_errors=True)

    # Write CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["threads", "run", "wall_seconds_external", "wall_seconds_pilots", "reader_seconds", "group_setup_seconds", "run_dir"]
    with out_csv.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"Wrote {out_csv} with {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
