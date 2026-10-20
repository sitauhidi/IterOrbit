#!/usr/bin/env python3
"""Builds analysis/data/: every input of analysis/analysis.py as a CSV, so the analysis never reads a log.

Needs the experiments' workspace (their logs and the sweeps' own CSVs), which is not part of this repository: by
default the folder that contains this IterOrbit checkout, or the folder named by $ITERORBIT_EXPERIMENTS. Parsing
only, no experiments are run:
    python analysis/build_data.py

Every absolute path in the data is written relative to the workspace (dataset-imsm/..., SOTA/..., ...).

Two kinds of input:
  COPIED   CSVs the sweeps' own parsers already wrote (imsm-analysis/csv, mccreesh-analysis/csv, SOTA/,
           the instance selections), copied unchanged under a name that says what they hold.
  PARSED   sweeps the analyses used to read from their logs, parsed here into one CSV each:
           iterorbit_runs.csv           every IterOrbit run of the extra sweeps (no accelerations, the
                                        one-at-a-time acceleration ablation, the OrbitSI-filter variants,
                                        the 2026-10-17 re-sweep), one row per run, timers as logged
           orbitsi_published_runs.csv   the OrbitSI package's published configuration, end to end and
                                        filter only, one row per run with its stage timers
           glasgow_count.csv            the Glasgow solver's counting runs on the unlabelled benchmark
           index_build_logs.csv         the index-build and OrbitSI load-split measurements
           labelled_graphs.csv          size and labels of the eight labelled data graphs

data/MANIFEST.csv lists every file with its source and row count.
"""
import csv
import glob
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ROOT = os.environ.get("ITERORBIT_EXPERIMENTS") or os.path.dirname(os.path.dirname(HERE))
WORKSPACE_DIRS = ("dataset-imsm", "dataset-mccreesh-graph", "imsm-analysis", "mccreesh-analysis", "SOTA", "IterOrbit", "OrbitSI")
ABSOLUTE = re.compile(r"(?:/[^/,\"\s]+)+?/(?=(?:%s)/)" % "|".join(WORKSPACE_DIRS))


def relative(text):
    """text with every absolute path into the workspace made relative to it."""
    return ABSOLUTE.sub("", text)


def write_csv(name, text):
    with open(os.path.join(DATA, name), "w") as f:
        f.write(relative(text))

COPIED = {
    "labelled_patterns.csv": "dataset-imsm/selected_profiling_patterns.csv",
    "unlabelled_pairs.csv": "dataset-mccreesh-graph/selected_profiling_pairs.csv",
    "labelled_filter_singlepass.csv": "imsm-analysis/csv/filter_efficacy_results.csv",
    "labelled_filter_fixpoint.csv": "imsm-analysis/csv/filter_efficacy_fixpoint_results.csv",
    "labelled_filter_fixpoint_rounds.csv": "imsm-analysis/csv/filter_efficacy_fixpoint_rounds.csv",
    "labelled_filter_no_orbit.csv": "imsm-analysis/csv/filter_efficacy_no_orbit_results.csv",
    "unlabelled_filter_singlepass.csv": "mccreesh-analysis/csv/filter_efficacy_results.csv",
    "unlabelled_filter_fixpoint.csv": "mccreesh-analysis/csv/filter_efficacy_fixpoint_results.csv",
    "unlabelled_filter_fixpoint_rounds.csv": "mccreesh-analysis/csv/filter_efficacy_fixpoint_rounds.csv",
    "unlabelled_filter_no_orbit.csv": "mccreesh-analysis/csv/filter_efficacy_no_orbit_results.csv",
    "labelled_iterorbit_e2e.csv": "SOTA/_orbitex_accelerated_sweep_results.csv",
    "labelled_sota_runs.csv": "SOTA/_sota_sweep_results.csv",
    "labelled_e2e_per_pattern.csv": "SOTA/_orbitex_vs_sota_per_pattern.csv",
    "unlabelled_e2e_phases.csv": "mccreesh-analysis/csv/e2e_phases.csv",
    "unlabelled_e2e_results.csv": "mccreesh-analysis/csv/e2e_results.csv",
    "unlabelled_glasgow_ground_truth.csv": "mccreesh-analysis/csv/glasgow_ground_truth.csv",
    "unlabelled_sota_filter_candidates.csv": "mccreesh-analysis/csv/sota_filter_candidates.csv",
    "unlabelled_induced_e2e.csv": "mccreesh-analysis/csv/induced_e2e.csv",
    "index_build_times.csv": "imsm-analysis/csv/index_build_times.csv",
    "enum_e2e.csv": "imsm-analysis/csv/enum_e2e.csv",
    "enum_mem.csv": "imsm-analysis/csv/enum_mem.csv",
    "accel_ablation.csv": "imsm-analysis/csv/accel_ablation.csv",
    "single_pass_e2e.csv": "imsm-analysis/csv/orbitsi_mode_e2e.csv",
    "accel_ablation_instances.csv": "imsm-analysis/logs/accel_ablation/instances.csv",
}

ITERORBIT_SWEEPS = [
    ("imsm", "no_accel", "imsm-analysis/logs/no_accel_e2e"),
    ("mccreesh", "no_accel", "mccreesh-analysis/logs/no_accel_e2e"),
    ("imsm", "accel_ablation", "imsm-analysis/logs/accel_ablation"),
    ("imsm", "orbitsi_filter", "imsm-analysis/logs/orbitsi_nosj_e2e"),
    ("imsm", "orbitsi_filter_nosym", "imsm-analysis/logs/orbitsi_nosj_nosym_e2e"),
    ("imsm", "resweep_20261017", "imsm-analysis/logs/iterorbit_e2e_20261017"),
    ("mccreesh", "resweep_20261017", "mccreesh-analysis/logs/iterorbit_e2e_20261017"),
]


def num(pattern, text, cast=float):
    m = re.search(pattern, text, re.M)
    return cast(m.group(1)) if m else None


def job_task(path):
    m = re.search(r"_(\d+)(?:_(\d+))?\.(?:out|err)$", path)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2)) if m.group(2) else None


def read(path):
    with open(path, errors="replace") as f:
        return f.read()


def iterorbit_run(path):
    t = read(path)
    err = path[:-4] + ".err"
    stderr = read(err) if os.path.exists(err) else ""
    job, task = job_task(path)
    cmd = re.search(r"^Command: (.*)$", t, re.M)
    flags = ""
    if cmd:
        toks = cmd.group(1).split()
        keep, skip = [], False
        for i, tok in enumerate(toks):
            if skip:
                skip = False
                continue
            if tok in ("--data", "--pattern"):
                skip = True
                continue
            if tok.startswith("--") and tok not in ("--kill-after=10s",):
                keep.append(tok)
                if i + 1 < len(toks) and not toks[i + 1].startswith("--"):
                    keep.append(toks[i + 1])
                    skip = True
        flags = " ".join(keep)
    sym = re.search(r"^Symmetry breaking: (\d+) levels, multiplier (\d+)( \(truncated\))?", t, re.M)
    return {
        "job": job, "task": task, "file": os.path.basename(path),
        "header_bench": num(r"^\s*Bench:\s+(\S+)", t, str),
        "dataset": num(r"^\s*Dataset:\s+(\S+)", t, str), "pattern_name": num(r"^\s*Pattern:\s+(\S+)", t, str),
        "pattern_path": num(r"^Command: .* --pattern (\S+)", t, str),
        "flags": flags,
        "status": num(r"^STATUS:? (\w+)", t, str) or "MISSING",
        "elapsed_ms": num(r"^ELAPSED_MS: (\d+)", t),
        "load_ms": num(r"Total loading time: ([\d.]+) ms", t),
        "filter_ms": num(r"Total filtering time: ([\d.]+) ms", t),
        "order_us": num(r"Total ordering time: ([\d.]+) us", t),
        "symmetry_ms": num(r"Total symmetry time: ([\d.]+) ms", t),
        "symmetry_levels": int(sym.group(1)) if sym else None,
        "symmetry_multiplier": sym.group(2) if sym else None,
        "symmetry_truncated": bool(sym and sym.group(3)),
        "search_ms": num(r"Total searching time: ([\d.]+) ms", t),
        "exec_ms": num(r"Total execution time: ([\d.]+) ms", t),
        "count": num(r"^Matches found: (\d+)", t, str),
        "saturated": "exceeds 64 bits" in stderr or "exceeds 64 bits" in t,
    }


def iterorbit_runs():
    rows = []
    for bench, sweep, d in ITERORBIT_SWEEPS:
        for p in sorted(glob.glob(os.path.join(ROOT, d, "*.out"))):
            r = iterorbit_run(p)
            if r["dataset"] is None:
                continue
            hb = r.pop("header_bench")
            rows.append({"bench": hb or bench, "sweep": sweep, **r})
    return pd.DataFrame(rows)


def orbitsi_published_runs():
    rows = []
    for bench in ("imsm", "mccreesh"):
        for sweep, sub, prefix, line in (("e2e", "orbitsi_published_e2e", "osi_pub_e2e", "STAGES"),
                                         ("filter", "orbitsi_published_filter", "osi_pub_filter", "FILTER")):
            for p in sorted(glob.glob(os.path.join(ROOT, f"{bench}-analysis", "logs", sub, f"{prefix}_*_[0-9]*.out"))):
                t = read(p)
                job, task = job_task(p)
                r = {"bench": bench, "sweep": sweep, "job": job, "task": task, "file": os.path.basename(p),
                     "dataset": num(r"^\s*Dataset:\s+(\S+)", t, str), "pattern_name": num(r"^\s*Pattern:\s+(\S+)", t, str),
                     "status": num(r"^STATUS:? (\w+)", t, str) or "MISSING",
                     "load_ms": num(r"^LOAD \S+ ([\d.]+) ms", t), "has_line": False}
                m = re.search(rf"^{line} .*$", t, re.M)
                if m:
                    r["has_line"] = True
                    r["line"] = m.group(0)
                    for k, v in re.findall(r"(\S+)=(\S+)", m.group(0)):
                        r[f"kv:{k}"] = v
                rows.append(r)
    return pd.DataFrame(rows)


def glasgow_count():
    rows = []
    for p in sorted(glob.glob(os.path.join(ROOT, "mccreesh-analysis", "logs", "glasgow_count", "*.txt"))):
        kv = dict(re.findall(r"^(\w+) = (.*)$", read(p), re.M))
        for k in ("command", "commandline", "pattern_file", "target_file"):
            kv.pop(k, None)
        rows.append({"file": os.path.basename(p), **kv})
    return pd.DataFrame(rows)


def index_build_logs():
    rows = []
    for p in sorted(glob.glob(os.path.join(ROOT, "imsm-analysis", "logs", "index_build_times", "*.out"))):
        for m in re.finditer(r"^(INDEX|LOADSPLIT) (\S+) (.*)$", read(p), re.M):
            r = {"file": os.path.basename(p), "kind": m.group(1), "path": m.group(2)}
            r.update({k: float(v) for k, v in re.findall(r"(\w+)=([\d.]+)", m.group(3))})
            rows.append(r)
    return pd.DataFrame(rows)


def labelled_graphs():
    rows = []
    for p in sorted(glob.glob(os.path.join(ROOT, "dataset-imsm", "data", "*.graph"))):
        name = os.path.basename(p)[:-6]
        if name == "test":
            continue
        n = m = 0
        labels = set()
        with open(p) as f:
            for line in f:
                c = line[:2]
                if c == "e ":
                    m += 1
                elif c == "v ":
                    n += 1
                    labels.add(line.split()[2])
        rows.append({"dataset": name, "vertices": n, "edges": m, "mean_degree": 2 * m / n, "labels": len(labels)})
    return pd.DataFrame(rows)


PARSED = {
    "iterorbit_runs.csv": iterorbit_runs,
    "orbitsi_published_runs.csv": orbitsi_published_runs,
    "glasgow_count.csv": glasgow_count,
    "index_build_logs.csv": index_build_logs,
    "labelled_graphs.csv": labelled_graphs,
}


def main(only=None):
    if not os.path.isdir(os.path.join(ROOT, "imsm-analysis")):
        sys.exit(f"{ROOT} is not the experiments' workspace: set ITERORBIT_EXPERIMENTS")
    os.makedirs(DATA, exist_ok=True)
    manifest = []
    for name, src in COPIED.items():
        if only and name not in only:
            continue
        s = os.path.join(ROOT, src)
        if not os.path.exists(s):
            print(f"MISSING {src}", file=sys.stderr)
            continue
        with open(s) as f:
            text = f.read()
        write_csv(name, text)
        rows = sum(1 for _ in csv.reader(text.splitlines())) - 1
        manifest.append({"file": name, "kind": "copied", "source": src, "rows": rows})
        print(f"copied  {name:42s} {rows:7d} rows  <- {src}")
    for name, build in PARSED.items():
        if only and name not in only:
            continue
        d = build()
        write_csv(name, d.to_csv(index=False))
        manifest.append({"file": name, "kind": "parsed", "source": build.__name__, "rows": len(d)})
        print(f"parsed  {name:42s} {len(d):7d} rows")
    if not only:
        pd.DataFrame(manifest).to_csv(os.path.join(DATA, "MANIFEST.csv"), index=False)


if __name__ == "__main__":
    main(set(sys.argv[1:]) or None)
