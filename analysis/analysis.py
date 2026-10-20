# %% [markdown]
# # IterOrbit: data, figures and tables
#
# Every number, table and figure of the paper, and the exploratory analyses behind it, from the CSVs in
# `analysis/data/` (built once from the experiments' logs by `analysis/build_data.py`; this script reads no log).
#
# Run it as a script (`python analysis/analysis.py`) or open it as a notebook: the `# %%` lines are cells
# (VS Code and Jupyter with Jupytext read them; `jupytext --to notebook analysis/analysis.py` writes an
# `.ipynb`). Each section is one function that prints its results and returns its main tables. In a notebook the
# figures are shown under their cells and not saved; as a script they are written as PDFs to `analysis/figures/`
# (the paper's figures to `../paper/figures/` instead when the paper's checkout sits next to this repository;
# `ANALYSIS_SAVE_FIGURES=1`/`0` overrides either default; `ANALYSIS_FIG_DIR` writes them elsewhere).
#
# Conventions used throughout:
# - **End-to-end time** is each system's own filter + order + enumerate timers, graph loading excluded, plus the
#   build time of the indexes its filter reads (`with_index`), floored at 0.1 ms (`TIME_FLOOR_MS`).
# - **Speed-up** of IterOrbit over a system is the geometric mean of (system's time / IterOrbit's) over the
#   instances both completed; "fastest" counts ties for every tied system.
# - **Reduction** is measured against round 1's LDF + NLF + arc-consistency candidates.

# %%
import os
import re
import sys
import warnings
from pathlib import Path

import matplotlib

INTERACTIVE = "ipykernel" in sys.modules
if not INTERACTIVE:
    matplotlib.use("Agg")
import matplotlib._mathtext
import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
pd.set_option("display.width", 220)


def _find_root():
    here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    for p in [here, *here.parents]:
        if (p / "analysis" / "data" / "MANIFEST.csv").exists():
            return p
    raise FileNotFoundError("analysis/data not found: run analysis/build_data.py first")


ROOT = _find_root()
DATA = ROOT / "analysis" / "data"
_PAPER = ROOT.parent / "paper" / "figures"   # the paper's checkout next to this repository, if there is one
PAPER_FIG_DIR = Path(os.environ.get("ANALYSIS_FIG_DIR", _PAPER if _PAPER.is_dir() else ROOT / "analysis" / "figures"))
EXTRA_FIG_DIR = Path(os.environ.get("ANALYSIS_FIG_DIR", ROOT / "analysis" / "figures"))
# In a notebook the figures are only shown; as a script they are written as PDFs. ANALYSIS_SAVE_FIGURES=1 or 0 overrides.
SAVE_FIGURES = os.environ.get("ANALYSIS_SAVE_FIGURES", "0" if INTERACTIVE else "1") == "1"


def csv(name, **kw):
    """One of the data files (analysis/data/MANIFEST.csv lists them and their sources)."""
    return pd.read_csv(DATA / name, **kw)


def _rel(p):
    return os.path.relpath(p, ROOT)


print(f"data: {_rel(DATA)}\n" + (f"paper figures -> {_rel(PAPER_FIG_DIR)}\nother figures -> {_rel(EXTRA_FIG_DIR)}"
                                  if SAVE_FIGURES else "figures: shown only, not saved"))

# %% [markdown]
# ## Toolkit
# Benchmarks, tools, colours, the print-scaled figure style (`visualisation_style.md`), and the loaders every
# section shares.

# %%
IMSM = ["hprd", "human", "patents", "wordnet", "yeast"]   # the five labelled graphs of the earlier drafts
IMSM_ALL = ["dblp", "eu2005", "hprd", "human", "patents", "wordnet", "yeast", "youtube"]
PLOT_IMSM = IMSM_ALL
MCC = ["LV", "images-CVIU11", "meshes-CVIU11", "scalefree", "si"]
LABEL = {"dblp": "DBLP", "eu2005": "EU2005", "youtube": "YouTube",
         "hprd": "HPRD", "human": "Human", "patents": "Patents", "wordnet": "WordNet", "yeast": "Yeast",
         "LV": "LV", "images-CVIU11": "Images", "meshes-CVIU11": "Meshes", "scalefree": "Scalefree",
         "si": "SI"}
BENCH = {"imsm": "Labelled", "mccreesh": "Unlabelled"}

S = sns.color_palette("Set2")
TOOLS = ["IterOrbit", "RM", "DP-iso", "VEQ", "CaLiG"]
TOOL_COLOR = {"IterOrbit": S[1], "RM": S[0], "VEQ": S[2], "DP-iso": S[3], "CaLiG": S[4], "OrbitSI": S[5]}
TOOL_MARKER = {"IterOrbit": "D", "RM": "s", "CaLiG": "^", "VEQ": "o", "DP-iso": "v", "OrbitSI": "P"}
PLOT_TOOLS = ["IterOrbit", "OrbitSI", "RM", "DP-iso", "VEQ", "CaLiG"]
SHORT = {"IterOrbit": "IO", "OrbitSI": "OSI"}
MODES = ["Baseline", "OrbitSI", "IterOrbit"]
MODE_LABEL = {"Baseline": "Baseline", "OrbitSI": r"$\mathbf{I_{max}=1}$ (≈ OSI)",
              "IterOrbit": r"$\mathbf{I_{max}=0}$ (IO, to fixpoint)"}
MODE_LABEL_PLAIN = {"Baseline": "Baseline", "OrbitSI": r"$\mathbf{I_{max}=1}$", "IterOrbit": r"$\mathbf{I_{max}=0}$"}
MODE_COLOR = {"Baseline": S[2], "OrbitSI": S[0], "IterOrbit": S[1]}
MODE_HATCH = {"Baseline": "//", "OrbitSI": "\\\\", "IterOrbit": "xx"}
SIZES = [3, 4, 5]
SIZE_MARKER = {3: "o", 4: "s", 5: "^"}
SIZE_LINESTYLE = {3: "solid", 4: "dashed", 5: "dotted"}

COLUMN_IN = 238.96 / 72.27   # svjour3 two-column \columnwidth
TEXT_IN = 495.0 / 72.27      # \textwidth
PRINT_TICK_PT = 7.5          # printed size of a style-guide size-10 tick label
FONT = 1.0                   # multiplier on every style-guide font size, set by setup()


def fs(size):
    """A style-guide font size (ticks 10, axis labels 12, global labels 14, ...), scaled for print."""
    return size * FONT


def setup(fig_width=None, print_width=COLUMN_IN):
    """Start a figure drawn fig_width inches wide that prints print_width wide: sets the font scale."""
    global FONT
    FONT = PRINT_TICK_PT / 10 * fig_width / print_width if fig_width else 1.0
    matplotlib._mathtext.SHRINK_FACTOR = 0.8
    sns.set_context("paper", font_scale=1.2 * FONT)
    sns.set_style("whitegrid")


def gmean(x, floor=1e-3):
    x = np.asarray(x, dtype=float)
    return float(np.exp(np.log(np.clip(x, floor, None)).mean())) if len(x) else np.nan


def style_axis(ax, title=None, tick_size=10):
    if title:
        ax.set_title(title, fontweight="bold", fontsize=fs(11))
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontsize(fs(tick_size))
    ax.tick_params(labelbottom=True, labelleft=True)
    sns.despine(ax=ax)


def _extent(fig, axes):
    """Union of the axes' full extents (tick labels, titles, ...) in figure fractions: x0, y0, x1, y1."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    boxes = [inv.transform(ax.get_tightbbox(r)) for ax in axes]
    return (min(b[0][0] for b in boxes), min(b[0][1] for b in boxes),
            max(b[1][0] for b in boxes), max(b[1][1] for b in boxes))


def block_heading(fig, axes, text, pad=0.008, fontsize=13):
    """A benchmark's name centred above a block of panels."""
    _, _, _, y1 = _extent(fig, axes)
    pos = [ax.get_position() for ax in axes]
    xc = (min(p.x0 for p in pos) + max(p.x1 for p in pos)) / 2
    fig.text(xc, y1 + pad, text, ha="center", va="bottom", weight="bold", fontsize=fs(fontsize))


def shared_ylabel(fig, axes, text, pad=0.01, fontsize=12):
    """One y-axis label for the whole figure, left of every panel's tick labels."""
    x0, _, _, _ = _extent(fig, axes)
    pos = [ax.get_position() for ax in axes]
    yc = (min(p.y0 for p in pos) + max(p.y1 for p in pos)) / 2
    fig.text(x0 - pad, yc, text, rotation=90, ha="right", va="center", weight="bold", fontsize=fs(fontsize))


def shared_xlabel(fig, axes, text, pad=0.008, fontsize=12):
    """One x-axis label for the whole figure, below every panel's tick labels."""
    _, y0, _, _ = _extent(fig, axes)
    pos = [ax.get_position() for ax in axes]
    xc = (min(p.x0 for p in pos) + max(p.x1 for p in pos)) / 2
    fig.text(xc, y0 - pad, text, ha="center", va="top", weight="bold", fontsize=fs(fontsize))


def bench_grid(figsize, hspace=0.62, wspace=0.32, block_gap=0, **adjust):
    """Two rows of four labelled panels, then one row of five unlabelled ones, every row the full width.
    Returns the figure, the labelled axes (PLOT_IMSM order) and the unlabelled ones (MCC order)."""
    fig = plt.figure(figsize=figsize)
    if block_gap:
        adjust = {**adjust, "bottom": adjust.get("bottom", 0.1) + block_gap}
    gs = fig.add_gridspec(3, 20, hspace=hspace, wspace=wspace * 4, **adjust)
    lab = [fig.add_subplot(gs[i // 4, (i % 4) * 5:(i % 4) * 5 + 5]) for i in range(8)]
    unl = [fig.add_subplot(gs[2, j * 4:j * 4 + 4]) for j in range(5)]
    for ax in unl:
        p = ax.get_position()
        ax.set_position([p.x0, p.y0 - block_gap, p.width, p.height])
    return fig, lab, unl


def equal_gaps(fig, axes, left, right, gap):
    """Repositions one row of axes so the space between neighbouring panels' full extents is `gap`."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    pads = []
    for ax in axes:
        box = ax.get_position()
        (x0, _), (x1, _) = inv.transform(ax.get_tightbbox(r))
        pads.append((box.x0 - x0, x1 - box.x1))
    w = (right - left - sum(a + b for a, b in pads) - gap * (len(axes) - 1)) / len(axes)
    x = left
    for ax, (lpad, rpad) in zip(axes, pads):
        box = ax.get_position()
        ax.set_position([x + lpad, box.y0, w, box.height])
        x += lpad + w + rpad + gap


def legend_row(fig, title, handles, y, fontsize=10, **kw):
    """One legend row with its title as the first entry (an invisible handle), centred."""
    t = plt.Line2D([0], [0], linestyle="none", marker=None, label=title)
    leg = fig.legend(handles=[t] + handles, bbox_to_anchor=(0.5, y), loc="center", ncol=len(handles) + 1,
                     frameon=False, **kw)
    for text in leg.get_texts():
        text.set_fontweight("bold")
        text.set_fontsize(fs(fontsize))
    leg.get_texts()[0].set_fontsize(fs(fontsize + 1))
    return leg


def row_major(handles, ncol):
    """Reorders legend handles so a legend with ncol columns reads in the handles' order row by row."""
    nrow = -(-len(handles) // ncol)
    rows = [handles[i * ncol:(i + 1) * ncol] for i in range(nrow)]
    return [r[c] for c in range(ncol) for r in rows if c < len(r)]


def save(fig, name, paper=True):
    """Shows the figure in a notebook and, with SAVE_FIGURES, writes it as PDF (paper figures to PAPER_FIG_DIR,
    the others to EXTRA_FIG_DIR); then closes it."""
    boxes = [t.get_bbox_patch() for t in fig.findobj(matplotlib.text.Text) if t.get_bbox_patch() is not None]
    extra = fig.get_default_bbox_extra_artists() + fig.texts + fig.legends + boxes if boxes else None
    if SAVE_FIGURES:
        out = PAPER_FIG_DIR if paper else EXTRA_FIG_DIR
        out.mkdir(parents=True, exist_ok=True)
        path = out / f"{name}.pdf"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.03, bbox_extra_artists=extra)
        print(f"Saved {path}")
    if INTERACTIVE:
        from IPython.display import display
        display(fig)
    plt.close(fig)


def ms_tick(v, _):
    """0.1, 1, 10 as plain numbers, 100 and up as bold powers of ten."""
    if v >= 100:
        return rf"$\mathbf{{10^{{{int(round(np.log10(v)))}}}}}$"
    return f"{v:g}"


def common_patterns(df, by=("dataset", "graphlet_size")):
    """Restrict every mode to the patterns present in all modes, per (dataset, graphlet size)."""
    out = []
    for _, grp in df.groupby(list(by)):
        common = None
        for m in grp["mode"].unique():
            pats = set(grp.loc[grp["mode"] == m, "pattern_name"])
            common = pats if common is None else common & pats
        out.append(grp[grp.pattern_name.isin(common)])
    return pd.concat(out, ignore_index=True)


# ---- filtering --------------------------------------------------------------------------------------
def imsm_filtering(datasets=None, modes=None):
    """Labelled single pass (Baseline = full graph, OrbitSI = candidate subgraph) and fixpoint (IterOrbit)
    filter runs, with reduction_pct and cost_ms; patterns intersected over the modes. EU2005's g=5 full-graph
    Baseline is left out (it completes 6 of 100 patterns)."""
    s = csv("labelled_filter_singlepass.csv")
    s["mode"] = np.where(s.use_full_graph == 1, "Baseline", "OrbitSI")
    f = csv("labelled_filter_fixpoint.csv")
    f["mode"] = "IterOrbit"
    cols = ["dataset", "pattern_name", "graphlet_size", "mode", "total_filter_ms", "reduction_vs_semijoin"]
    d = pd.concat([s[cols], f[cols]], ignore_index=True)
    d = d[d.dataset.isin(IMSM if datasets is None else datasets)].copy()
    d = d[~((d.dataset == "eu2005") & (d.graphlet_size == 5) & (d["mode"] == "Baseline"))]
    d = d[d["mode"].isin(MODES if modes is None else modes)]
    d["reduction_pct"] = 100.0 * d.reduction_vs_semijoin
    d["cost_ms"] = d.total_filter_ms
    return common_patterns(d)


def mcc_filtering():
    """As imsm_filtering() for the unlabelled pairs. A refuted pair counts as 100% at its full filtering
    time; pairs LDF+NLF+arc consistency already refute in round 1 are flagged (cheap_refuted)."""
    s = csv("unlabelled_filter_singlepass.csv")
    s["mode"] = np.where(s.use_full_graph == 1, "Baseline", "OrbitSI")
    f = csv("unlabelled_filter_fixpoint.csv")
    f["mode"] = "IterOrbit"
    cols = ["dataset", "pattern_name", "graphlet_size", "mode", "outcome", "refuted_stage", "refuted_round",
            "total_filter_ms", "reduction_vs_semijoin"]
    d = pd.concat([s[cols], f[cols]], ignore_index=True)
    d = d[d.outcome.isin(["completed", "refuted"])].copy()
    d["cheap_refuted"] = ((d.outcome == "refuted") & (d.refuted_round.fillna(1) == 1)
                          & d.refuted_stage.isin(["LDF", "NLF", "Semi-Join"]))
    d["reduction_pct"] = 100.0 * d.reduction_vs_semijoin
    d["cost_ms"] = d.total_filter_ms
    return d


# ---- end to end ---------------------------------------------------------------------------------------
INCLUDE_LOADING = False   # end-to-end times leave graph loading out
INCLUDE_INDEX = True      # ... and charge each filter's index build (with_index)
TIME_FLOOR_MS = 0.1       # the baselines' framework prints seconds to four decimals


def imsm_phases(datasets=None):
    """Long frame tool, dataset, pattern_name, read/filter/order/enumerate (ms): every completed labelled
    counting run of IterOrbit and the four baselines on the given datasets (default IMSM)."""
    ox = csv("labelled_iterorbit_e2e.csv")
    ox = ox[ox.status == "COMPLETED"]
    a = pd.DataFrame({"tool": "IterOrbit", "dataset": ox.dataset, "pattern_name": ox.pattern_name,
                      "read": ox.load_time_ms, "filter": ox.filter_time_ms,
                      "order": ox.order_time_us.fillna(0) / 1000.0, "enumerate": ox.search_time_ms.fillna(0)})
    s = csv("labelled_sota_runs.csv")
    s = s[s.status == "COMPLETED"]
    b = pd.DataFrame({"tool": s.tool.map({"rapidmatch": "RM", "dpiso": "DP-iso", "veq": "VEQ", "calig": "CaLiG"}),
                      "dataset": s.dataset, "pattern_name": s.pattern_name,
                      "read": s.load_graphs_time_s * 1e3, "filter": s.filter_vertices_time_s * 1e3,
                      "order": s.generate_query_plan_time_s * 1e3,
                      "enumerate": (s.build_table_time_s.fillna(0) + s.enumerate_time_s.fillna(0)) * 1e3})
    d = pd.concat([a, b], ignore_index=True)
    return d[d.dataset.isin(IMSM if datasets is None else datasets)].dropna(subset=["read", "filter"])


def mcc_phases():
    """As imsm_phases() for the unlabelled pairs."""
    d = csv("unlabelled_e2e_phases.csv")
    d = d[d.status == "COMPLETED"].replace({"tool": {"OrbitEX": "IterOrbit", "RapidMatch": "RM"}})
    d = d[d.tool.isin(TOOLS)]
    return pd.DataFrame({"tool": d.tool, "dataset": d.dataset, "pattern_name": d.pattern_name,
                         "read": d.read_ms, "filter": d.filter_ms, "order": d.order_ms.fillna(0),
                         "enumerate": d.enumerate_ms.fillna(0)}).dropna(subset=["read", "filter"])


def phase_time(d):
    t = d["filter"] + d["order"] + d["enumerate"]
    if INCLUDE_INDEX and "index" in d:
        t = t + d["index"]
    return (t + d["read"] if INCLUDE_LOADING else t).clip(lower=TIME_FLOOR_MS)


def with_index(d, bench):
    """d (a phase frame of one benchmark) with an "index" column (ms): the build time of the indexes the
    tool's filter reads on that instance's data graph; 0 for IterOrbit; GuP's own index_ms."""
    ix = csv("index_build_times.csv")
    ix = ix[ix.bench == bench].set_index("graph")
    if bench == "imsm":
        graph = d.dataset
    else:
        pairs = csv("unlabelled_pairs.csv")
        target = dict(zip(zip(pairs.dataset, pairs.instance), pairs.target_path))
        graph = pd.Series([target[k] for k in zip(d.dataset, d.pattern_name)], index=d.index)
    col = np.where(d.tool == "OrbitSI", "orbitsi_used_ms", "framework_used_ms")
    used = np.array([ix.at[g, c] for g, c in zip(graph, col)])
    gup = d.index_ms.to_numpy() if "index_ms" in d else np.full(len(d), np.nan)
    return d.assign(index=np.select([d.tool == "IterOrbit", d.tool == "GuP"], [0.0, gup], used))


def all_completed(d, tools=TOOLS):
    """(dataset, pattern_name) keys every one of the tools completed."""
    k = d[d.tool.isin(tools)].groupby(["dataset", "pattern_name"]).tool.nunique()
    return set(k[k == len(tools)].index)


def instance_keys(bench):
    """Every instance of a benchmark: 800 labelled patterns or 480 unlabelled pairs, as (dataset, pattern_name)."""
    if bench == "imsm":
        return [tuple(x) for x in csv("labelled_patterns.csv")[["dataset", "pattern_name"]].values]
    return [tuple(x) for x in csv("unlabelled_pairs.csv")[["dataset", "instance"]].values]


# ---- the OrbitSI package (published configuration, C++ core) --------------------------------------------
def orbitsi_runs(bench, sweep):
    """The package's runs: sweep "e2e" (counting) or "filter" (--filter-only); `line` is the run's STAGES or
    FILTER line as logged."""
    o = csv("orbitsi_published_runs.csv", dtype={"line": str}, low_memory=False)
    return o[(o.bench == bench) & (o.sweep == sweep)]


def _stage_ms(line):
    return {k: float(v) for k, v in re.findall(r"(\S+)=([\d.]+)", line) if not k.startswith(("n_", "m_", "count"))}


def orbitsi_phases(bench):
    """The package's completed end-to-end runs in the frame of imsm_phases(): filter = pattern orbits +
    LDF/NLF + orbit stage + "front" (the Python side of the call); order; enumerate = search; read = loading."""
    rows = []
    for r in orbitsi_runs(bench, "e2e").itertuples():
        if not isinstance(r.line, str):
            continue
        kv = dict(re.findall(r"(\S+)=(\S+)", r.line))
        ms = {k: float(v) for k, v in kv.items() if not k.startswith(("n_", "m_", "count")) and "." not in k}
        rows.append({"tool": "OrbitSI", "dataset": r.dataset, "pattern_name": r.pattern_name, "read": r.load_ms,
                     "filter": sum(v for k, v in ms.items() if k not in ("order", "search")),
                     "order": ms.get("order", 0.0), "enumerate": ms.get("search", 0.0)})
    return pd.DataFrame(rows)


def orbitsi_filter(bench):
    """The package's filter on every instance it finished: refuted, cand (mean candidates per pattern vertex,
    0 when refuted) and cost_ms (the filter stages plus the Python side of the call)."""
    rows = []
    for r in orbitsi_runs(bench, "filter").itertuples():
        if not isinstance(r.line, str):
            continue
        kv = dict(re.findall(r"(\S+)=(\S+)", r.line))
        ms = {k: float(v) for k, v in kv.items()
              if k not in ("status", "n_pattern", "cand_sum", "cand_avg") and k not in ("order", "search")}
        rows.append({"dataset": r.dataset, "pattern_name": r.pattern_name, "refuted": kv["status"] == "refuted",
                     "cand": float(kv["cand_avg"]), "cost_ms": sum(ms.values())})
    return pd.DataFrame(rows)


def orbitsi_e2e_totals(bench):
    """Per run of the package's counting sweep: load (ms), t (sum of its top-level stage timers, completed runs
    only), count (exact integer) and wall = load + t; indexed by (dataset, pattern_name)."""
    rows = []
    for r in orbitsi_runs(bench, "e2e").itertuples():
        x = {"dataset": r.dataset, "pattern_name": r.pattern_name, "load": r.load_ms}
        if isinstance(r.line, str):
            x["t"] = sum(v for k, v in _stage_ms(r.line).items() if "." not in k)
            x["count"] = int(re.search(r" count=(\d+)", r.line).group(1))
        rows.append(x)
    d = pd.DataFrame(rows).set_index(["dataset", "pattern_name"])
    d["wall"] = d.load + d.t
    return d


# ---- IterOrbit's extra sweeps (no accelerations, ablations, OrbitSI-filter variants, re-sweep) -----------
def iterorbit_sweep(sweep, bench=None):
    """Runs of one of IterOrbit's extra sweeps (iterorbit_runs.csv), with filter/order/enumerate in ms (order
    converted from us) and count as a float, converted by Python's float() as the logs' parsers did (pandas'
    own string parser rounds 2^63-1, the saturated count, to a different float)."""
    r = csv("iterorbit_runs.csv", dtype={"count": str, "symmetry_multiplier": str})
    r = r[r.sweep == sweep]
    if bench is not None:
        r = r[r.bench == bench]
    return r.assign(filter=r.filter_ms, order=r.order_us / 1e3, enumerate=r.search_ms,
                    count=r["count"].map(lambda c: float(c) if isinstance(c, str) else np.nan))


def iterorbit_resweep(bench):
    """The IterOrbit end-to-end re-sweep of 2026-10-17: status, t (filter + order + search, ms), sym (symmetry
    time), t_sym, count, indexed by (dataset, pattern_name)."""
    r = iterorbit_sweep("resweep_20261017", bench)
    done = r.status == "COMPLETED"
    d = pd.DataFrame({"dataset": r.dataset, "pattern_name": r.pattern_name, "status": r.status,
                      "t": (r.filter_ms + r.order_us.fillna(0) / 1e3 + r.search_ms.fillna(0)).where(done),
                      "sym": r.symmetry_ms.fillna(0).where(done), "count": r["count"].where(done),
                      "elapsed_ms": r.elapsed_ms})
    d = d.set_index(["dataset", "pattern_name"])
    d["t_sym"] = d.t + d.sym
    return d


def glasgow_runs():
    """The Glasgow solver's counting runs on the unlabelled benchmark: decided (status true/false), runtime
    (ms, its own timer, after reading the graphs; whole ms), wall (ms), count; indexed by (dataset, pattern_name)."""
    g = csv("glasgow_count.csv", dtype={"status": str, "solution_count": str})
    g = g.rename(columns={"instance": "pattern_name"})
    g["decided"] = g.status.isin(["true", "false"])
    g["runtime"] = pd.to_numeric(g.runtime, errors="coerce").where(g.decided)
    g["wall"] = (pd.to_numeric(g.wall_s, errors="coerce") * 1000.0).where(g.decided)
    return g.set_index(["dataset", "pattern_name"])


# ---- tool order in the figures ------------------------------------------------------------------------
_ORDER = None


def plot_tool_order():
    """PLOT_TOOLS sorted by overall end-to-end time, slowest first: each tool's geometric mean over the
    instances all six completed, both benchmarks pooled."""
    global _ORDER
    if _ORDER is None:
        frames = []
        for bench, loader in (("imsm", lambda: imsm_phases(PLOT_IMSM)), ("mccreesh", mcc_phases)):
            d = with_index(pd.concat([loader(), orbitsi_phases(bench)], ignore_index=True), bench)
            keys = all_completed(d, PLOT_TOOLS)
            frames.append(d[[k in keys for k in zip(d.dataset, d.pattern_name)]])
        d = pd.concat(frames, ignore_index=True).assign(t=lambda f: phase_time(f))
        g = {t: gmean(d[d.tool == t].t, TIME_FLOOR_MS) for t in PLOT_TOOLS}
        _ORDER = sorted(PLOT_TOOLS, key=g.get, reverse=True)
        print("tool order (overall geo ms): " + ", ".join(f"{t} {g[t]:.3g}" for t in _ORDER))
    return _ORDER


# %% [markdown]
# # Paper
# ## Benchmarks (Tables `tab:datasets-overview`, `tab:mccreesh-overview`)

# %%
def benchmarks():
    """The labelled data graphs (vertices, edges, mean degree, labels) and the unlabelled benchmark's families
    (pairs, pattern and target vertex ranges, pairs without an embedding per the Glasgow solver)."""
    lab = csv("labelled_graphs.csv").assign(dataset=lambda d: d.dataset.map(LABEL)).set_index("dataset")
    print("Labelled data graphs:\n" + lab.round(3).to_string())
    pairs = csv("unlabelled_pairs.csv")
    gt = csv("unlabelled_glasgow_ground_truth.csv", dtype={"glasgow_status": str})
    nomatch = gt[gt.glasgow_status.str.lower() == "false"].groupby("dataset").size()
    unl = pairs.groupby("dataset").agg(pairs=("instance", "size"), pattern_min=("pattern_n", "min"),
                                       pattern_max=("pattern_n", "max"), target_min=("target_n", "min"),
                                       target_max=("target_n", "max")).reindex(MCC)
    unl["no_match"] = nomatch.reindex(MCC).fillna(0).astype(int)
    unl.index = [LABEL[d] for d in unl.index]
    print("\nUnlabelled families:\n" + unl.to_string())
    return lab, unl


benchmarks();

# %% [markdown]
# ## Filtering on labelled graphs and without labels
# The single pass ($I_{max}=1$, "OrbitSI" mode in the CSVs) and the fixpoint ($I_{max}=0$, "IterOrbit"), as
# plotted in `filter_tradeoff.pdf` (mean reduction, geometric-mean time) and `iteration_gain.pdf`.

# %%
def filtering_section():
    """Every number of the "Filtering on Labelled Graphs" and "Filtering Without Labels" subsections."""
    MODES2 = ["OrbitSI", "IterOrbit"]
    NAME = {"OrbitSI": "single", "IterOrbit": "fixpoint"}

    def table(d, title):
        m = d.groupby(["dataset", "graphlet_size", "mode"]).agg(red=("reduction_pct", "mean"), ms=("cost_ms", gmean),
                                                               n=("pattern_name", "size")).unstack("mode")
        print(f"\n== {title}: mean reduction % / geo ms per dataset, g (single | fixpoint | +pp, cost x)")
        for (ds, g), r in m.iterrows():
            print(f"  {LABEL[ds]:10s} g={g} n={int(r['n']['OrbitSI']):3d}  single {r['red']['OrbitSI']:5.1f}% {r['ms']['OrbitSI']:9.2f}ms"
                  f" | fix {r['red']['IterOrbit']:5.1f}% {r['ms']['IterOrbit']:9.2f}ms | +{r['red']['IterOrbit'] - r['red']['OrbitSI']:5.1f}pp"
                  f" {r['ms']['IterOrbit'] / r['ms']['OrbitSI']:.2f}x")
        pooled = d.groupby(["graphlet_size", "mode"]).agg(red=("reduction_pct", "mean"), ms=("cost_ms", gmean)).unstack("mode")
        print("  pooled:\n" + pooled.round(2).to_string())
        per = m["red"].groupby(level="graphlet_size").mean()
        print("  mean over datasets of per-dataset mean reduction:\n" + per.round(2).to_string())
        return m

    def rounds(fixpoint_csv, datasets, g=None):
        f = csv(fixpoint_csv)
        f = f[f.dataset.isin(datasets) & f.reduction_vs_semijoin.notna()]
        if g is not None:
            f = f[f.graphlet_size == g]
        return f

    imsm = imsm_filtering(PLOT_IMSM, MODES2)
    m = table(imsm, "Labelled")
    others = [d for d in PLOT_IMSM if d != "wordnet"]
    s = m["red"]["OrbitSI"]
    print("  single pass, mean over the seven datasets other than WordNet: "
          + " ".join(f"g={g} {s.loc[others].xs(g, level='graphlet_size').mean():.1f}%" for g in SIZES))
    print("  single pass WordNet: " + " ".join(f"g={g} {s.loc['wordnet', g]:.1f}%" for g in SIZES))
    gain = (m["red"]["IterOrbit"] - m["red"]["OrbitSI"]).groupby(level="graphlet_size").mean()
    ratio = (m["ms"]["IterOrbit"] / m["ms"]["OrbitSI"])
    print("  fixpoint - single, mean over datasets (pp): " + " ".join(f"g={g} {v:.1f}" for g, v in gain.items()))
    print("  fixpoint / single cost, geo mean over datasets: " + " ".join(f"g={g} {gmean(ratio.xs(g, level='graphlet_size')):.2f}x" for g in SIZES)
          + f"  range {ratio.min():.2f}-{ratio.max():.2f}x")
    f = rounds("labelled_filter_fixpoint.csv", PLOT_IMSM)
    print(f"  rounds to converge (all g, {len(f)} runs): mean {f.n_rounds.mean():.1f} median {f.n_rounds.median():.0f} max {f.n_rounds.max()}")
    for ds, x in f.groupby("dataset"):
        print(f"    {LABEL[ds]:10s} mean {x.n_rounds.mean():.1f} median {x.n_rounds.median():.0f} max {x.n_rounds.max()}")
    f4 = f[f.graphlet_size == 4]
    print(f"  g=4: mean {f4.n_rounds.mean():.1f} median {f4.n_rounds.median():.0f} max {f4.n_rounds.max()}")

    mcc = mcc_filtering()
    cheap = mcc[mcc.cheap_refuted]
    print(f"\nUnlabelled: pairs refuted in round 1 by LDF+NLF+AC: {cheap.drop_duplicates(['dataset', 'pattern_name']).shape[0]}")
    for mode in MODES2:
        for g in SIZES:
            x = mcc[(mcc["mode"] == mode) & (mcc.graphlet_size == g) & (mcc.outcome == "refuted")]
            print(f"  {NAME[mode]:8s} g={g}: refuted {x.pattern_name.size}")
    mm = table(common_patterns(mcc[~mcc.cheap_refuted & mcc["mode"].isin(MODES2)]), "Unlabelled (round-1 refutations left out)")
    full = mcc[mcc["mode"].isin(MODES2)]
    t = full.groupby(["graphlet_size", "mode"]).cost_ms.agg([gmean, "max"]).unstack("mode")
    print("  geo filtering ms over all pairs (incl. cheap refutations), and max:\n" + t.round(2).to_string())
    for g in SIZES:
        f = rounds("unlabelled_filter_fixpoint.csv", MCC, g)
        print(f"  g={g} rounds (all pairs): mean {f.n_rounds.mean():.1f} median {f.n_rounds.median():.0f} max {f.n_rounds.max()}")
        im = f[f.dataset == "images-CVIU11"]
        print(f"     Images: median rounds {im.n_rounds.median():.0f}, geo ms {gmean(im.total_filter_ms):.1f}")
    f5 = rounds("unlabelled_filter_fixpoint.csv", MCC, 5)
    for ds, x in f5.groupby("dataset"):
        print(f"    g=5 {LABEL[ds]:10s} mean {x.n_rounds.mean():.1f} median {x.n_rounds.median():.0f} max {x.n_rounds.max()}")
    return m, mm


filtering_section();

# %% [markdown]
# ## Choosing the iteration cap and graphlet size
# From the labelled fixpoint trajectories; a run that converged earlier is carried forward.

# %%
def params_section():
    """Every number of the "Choosing the Iteration Cap and Graphlet Size" subsection."""
    r = csv("labelled_filter_fixpoint_rounds.csv")
    r = r[r.dataset.isin(PLOT_IMSM)].sort_values("round")
    traj = []
    for (d, p, g), grp in r.groupby(["dataset", "pattern_name", "graphlet_size"]):
        pct = grp.set_index("round").pct_of_sj_baseline
        cost = grp.set_index("round").cumulative_cost_ms
        last = pct.index.max()
        row = dict(dataset=d, pattern_name=p, graphlet_size=g, n_rounds=last, final_red=100.0 - pct.loc[last], final_cost=cost.loc[last])
        for t in range(1, 6):
            tt = min(t, last)
            row[f"red_{t}"], row[f"cost_{t}"] = 100.0 - pct.loc[tt], cost.loc[tt]
        traj.append(row)
    traj = pd.DataFrame(traj)
    tr = traj[traj.final_red > 1e-9]
    print(f"trajectories {len(traj)}, with a gain beyond arc consistency {len(tr)}, none {len(traj) - len(tr)}")
    print(tr.groupby("dataset").size().to_string())
    rows = []
    for t in range(1, 6):
        cap = 100 * tr[f"red_{t}"] / tr.final_red
        cst = 100 * tr[f"cost_{t}"] / tr.final_cost
        rows.append(dict(I_max=t, cap_median=cap.median(), cap_mean=cap.mean(), cost_median=cst.median(), cost_mean=cst.mean()))
    caps = pd.DataFrame(rows)
    print(caps.round(1).to_string(index=False))
    c3 = tr.assign(c=100 * tr.red_3 / tr.final_red).groupby("dataset").c.agg(["mean", "min"])
    print("captured by round 3 per dataset (mean, min):\n" + c3.round(1).to_string())
    print("rounds per dataset (all g): max\n" + traj.groupby("dataset").n_rounds.max().to_string())
    print("rounds by dataset x g max:\n" + traj.pivot_table(index="dataset", columns="graphlet_size", values="n_rounds", aggfunc="max").to_string())

    k = traj[["dataset", "pattern_name", "graphlet_size", "red_3", "cost_3"]]
    n = k.groupby(["dataset", "pattern_name"]).graphlet_size.nunique()
    k = k.merge(n[n == 3].reset_index()[["dataset", "pattern_name"]])
    print("\npatterns with all three g:", k.groupby("dataset").pattern_name.nunique().to_dict())
    per = k.groupby(["dataset", "graphlet_size"]).agg(red_mean=("red_3", "mean"), cost_mean=("cost_3", "mean"), cost_geo=("cost_3", gmean))
    print(per.round(2).unstack().to_string())
    for label, sub in (("all eight", per), ("without EU2005", per.drop("eu2005", level="dataset"))):
        avg = sub.groupby("graphlet_size").mean()
        avg["roi_mean"] = avg.red_mean.diff() / avg.cost_mean.diff()
        avg["roi_geo"] = avg.red_mean.diff() / avg.cost_geo.diff()
        avg["cost_ratio_mean"] = avg.cost_mean / avg.cost_mean.shift()
        print(f"\naveraged across datasets ({label}):\n" + avg.to_string(float_format=lambda v: f"{v:.4g}"))
    w = per.unstack("graphlet_size")
    print("\nper dataset g4->g5: +pp, geo cost ratio, mean cost ratio; g3->g4 likewise")
    for ds in w.index:
        print(f"  {LABEL[ds]:9s} g4->5 +{w.red_mean[5][ds] - w.red_mean[4][ds]:.2f}pp ({w.red_mean[4][ds]:.2f}->{w.red_mean[5][ds]:.2f}%)"
              f" x{w.cost_geo[5][ds] / w.cost_geo[4][ds]:.1f} geo ({w.cost_geo[4][ds]:.1f}->{w.cost_geo[5][ds]:.1f} ms)"
              f" | g3->4 +{w.red_mean[4][ds] - w.red_mean[3][ds]:.2f}pp x{w.cost_geo[4][ds] / w.cost_geo[3][ds]:.2f}")
    return caps, per


params_section();

# %% [markdown]
# ### Figure `filter_tradeoff.pdf`
# Filtering time vs mean reduction, single pass vs fixpoint, four datasets per benchmark; markers are the
# geometric-mean time and arithmetic-mean reduction, lines their 95% confidence intervals.

# %%
def fig_filter_tradeoff():
    setup(6.6)
    SHOWN = {"imsm": ["wordnet", "patents", "eu2005", "human"],
             "mccreesh": ["LV", "images-CVIU11", "scalefree", "si"]}
    MODES2 = ["OrbitSI", "IterOrbit"]
    Z = 1.96

    def tick(v, _):
        return f"{v / 1000:g}k" if v >= 1000 else f"{v:g}"

    def panel(ax, d, title):
        m = d.groupby(["mode", "graphlet_size"]).agg(
            x=("cost_ms", lambda c: gmean(c)), y=("reduction_pct", "mean"),
            lse=("cost_ms", lambda c: Z * np.log(np.clip(c, 1e-3, None)).std(ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0.0),
            yse=("reduction_pct", lambda r: Z * r.std(ddof=1) / np.sqrt(len(r)) if len(r) > 1 else 0.0),
        ).reset_index()
        m["x1"], m["x3"] = m.x / np.exp(m.lse), m.x * np.exp(m.lse)
        m["y1"], m["y3"] = m.y - m.yse, m.y + m.yse
        for _, r in m.iterrows():
            col, ls = MODE_COLOR[r["mode"]], SIZE_LINESTYLE[r.graphlet_size]
            _, _, bars = ax.errorbar(r.x, r.y, xerr=[[r.x - r.x1], [r.x3 - r.x]], yerr=[[r.y - r.y1], [r.y3 - r.y]],
                                     fmt="none", ecolor=col, elinewidth=1.1, capsize=2, capthick=1.1, alpha=0.8, zorder=2)
            for b in bars:
                b.set_linestyle(ls)
            ax.plot(r.x, r.y, marker=SIZE_MARKER[r.graphlet_size], markersize=7, markerfacecolor="none",
                    markeredgecolor=col, markeredgewidth=1.5, linestyle="none", zorder=3)
        ax.set_xscale("log")
        lo, hi = min(0.0, m.y1.min()), max(m.y.max(), m.y3.max())
        span = max(hi - lo, 1.0)
        ax.set_ylim(lo - 0.08 * span, hi + 0.12 * span)
        if hi < 1.0:
            ax.set_ylim(-0.25, 5)
        xmin, xmax = m.x1.min(), m.x3.max()
        pad = max((xmax / xmin) ** 0.12, 1.25)
        ax.set_xlim(xmin / pad, xmax * pad)
        decades = np.log10(ax.get_xlim()[1] / ax.get_xlim()[0])
        if decades < 0.6:
            ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=3))
        else:
            ax.xaxis.set_major_locator(mticker.LogLocator(base=10, subs=(1, 2, 5) if decades < 2.2 else (1,), numticks=12))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(tick))
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
        style_axis(ax, title, tick_size=10)
        ax.title.set_fontsize(fs(11))

    fig = plt.figure(figsize=(6.6, 8.8))
    BLOCK_GAP = 0.04
    gs = fig.add_gridspec(4, 2, hspace=0.78, wspace=0.26, top=0.85, bottom=0.075 + BLOCK_GAP, left=0.12, right=0.975)
    axes = [fig.add_subplot(gs[r, c]) for r in range(4) for c in range(2)]
    for ax in axes[4:]:
        p = ax.get_position()
        ax.set_position([p.x0, p.y0 - BLOCK_GAP, p.width, p.height])
    imsm = imsm_filtering(SHOWN["imsm"], MODES2)
    for ax, ds in zip(axes[:4], SHOWN["imsm"]):
        panel(ax, imsm[imsm.dataset == ds], LABEL[ds])
    mcc = mcc_filtering()
    mcc = common_patterns(mcc[~mcc.cheap_refuted & mcc["mode"].isin(MODES2)])
    for ax, ds in zip(axes[4:], SHOWN["mccreesh"]):
        panel(ax, mcc[mcc.dataset == ds], LABEL[ds])
    block_heading(fig, axes[:4], "Labelled", fontsize=12)
    block_heading(fig, axes[4:], "Unlabelled", fontsize=12)
    shared_ylabel(fig, axes, "Mean reduction (%)")
    shared_xlabel(fig, axes, "Total filtering time (ms, log scale)")
    mode_h = [plt.Line2D([0], [0], marker="o", linestyle="none", markersize=8, markerfacecolor="none",
                         markeredgecolor=MODE_COLOR[m], markeredgewidth=1.6, label=MODE_LABEL[m]) for m in MODES2]
    size_h = [plt.Line2D([0], [0], marker=SIZE_MARKER[s], linestyle=SIZE_LINESTYLE[s], color="black", markersize=8,
                         markerfacecolor="none", markeredgecolor="black", markeredgewidth=1.6, label=f"K{s}") for s in SIZES]
    leg = fig.legend(handles=mode_h, bbox_to_anchor=(0.5, 0.975), loc="center", ncol=2, frameon=False,
                     columnspacing=1.0, handletextpad=0.2)
    for t in leg.get_texts():
        t.set_fontweight("bold")
        t.set_fontsize(fs(10))
    legend_row(fig, "Graphlet Size:", size_h, 0.94, fontsize=10, columnspacing=1.0, handletextpad=0.4, handlelength=2.4)
    save(fig, "filter_tradeoff")


fig_filter_tradeoff()

# %% [markdown]
# ### Figure `iteration_gain.pdf`
# Per-round cost and gain of iterating in the end-to-end configuration (labelled g=4, solid to $I_{max}=3$;
# unlabelled g=5 to the fixpoint).

# %%
def fig_iteration_gain():
    setup(6.6)
    Z = 1.96
    E2E_ROUND = 3

    def trajectories(rounds_csv, g, datasets, drop=frozenset()):
        r = csv(rounds_csv)
        r = r[(r.graphlet_size == g) & r.dataset.isin(datasets)]
        r = r[[k not in drop for k in zip(r.dataset, r.pattern_name)]]
        out = []
        for ds, d in r.groupby("dataset"):
            T = int(d["round"].max())
            red = d.pivot_table(index="pattern_name", columns="round", values="pct_of_sj_baseline").reindex(columns=range(1, T + 1)).ffill(axis=1)
            cost = d.pivot_table(index="pattern_name", columns="round", values="cumulative_cost_ms").reindex(columns=range(1, T + 1)).ffill(axis=1)
            gain = red[[1]].values - red
            lratio = np.log(cost.clip(lower=1e-3).div(cost[1].clip(lower=1e-3), axis=0))
            for t in range(1, T + 1):
                out.append(pd.DataFrame({"dataset": ds, "round": t, "gain": gain[t], "lratio": lratio[t]}))
        return pd.concat(out, ignore_index=True)

    def panel(ax, p, datasets, title, mark_round=None):
        m = p.groupby(["dataset", "round"]).agg(
            gain=("gain", "mean"), gse=("gain", lambda x: Z * x.std(ddof=1) / np.sqrt(len(x))),
            lr=("lratio", "mean"), lse=("lratio", lambda x: Z * x.std(ddof=1) / np.sqrt(len(x))), n=("gain", "size"))
        colors = dict(zip(datasets, sns.color_palette("Dark2", 8)))
        for ds in datasets:
            r = m.loc[ds]
            T = r.index.max()
            k = min(mark_round or T, T)
            x, y = np.exp(r.lr), r.gain
            col = colors[ds]
            ax.plot(x.loc[:k], y.loc[:k], color=col, linewidth=1.6, marker="o", markersize=4, markerfacecolor="white",
                    markeredgecolor=col, zorder=2)
            if T > k:
                ax.plot(x.loc[k:], y.loc[k:], color=col, linewidth=1.2, linestyle="dotted", marker="o", markersize=3,
                        markerfacecolor="white", markeredgecolor=col, alpha=0.8, zorder=2)
            e = r.loc[k]
            xc = np.exp(e.lr)
            ax.errorbar(xc, e.gain, xerr=[[xc - xc / np.exp(e.lse)], [xc * np.exp(e.lse) - xc]], yerr=e.gse,
                        fmt=TOOL_MARKER["IterOrbit"], markersize=8, markerfacecolor=col, markeredgecolor="black",
                        markeredgewidth=0.8, ecolor=col, elinewidth=1.2, capsize=3, capthick=1.2, zorder=4)
            last = r.loc[T]
            print(f"  {ds:14s} n={int(e.n):3d}  rounds 1..{T}: " + "  ".join(f"r{t} {np.exp(r.loc[t].lr):.2f}x/{r.loc[t].gain:.1f}pp" for t in r.index[: min(T, 6)])
                  + (f"  ... r{T} {np.exp(last.lr):.2f}x/{last.gain:.1f}pp" if T > 6 else ""))
        ax.set_xscale("log")
        ax.axhline(0, color="gray", linewidth=0.8, zorder=1)
        ax.xaxis.set_major_locator(mticker.FixedLocator([1, 1.5, 2, 3, 4, 5]))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}" + r"$\times$"))
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xlim(0.97, np.exp(m.lr + m.lse).max() * 1.08)
        ax.set_ylim(-2, (m.gain + m.gse).max() * 1.08)
        style_axis(ax, title)
        handles = [plt.Line2D([0], [0], color=colors[ds], linewidth=3, label=LABEL[ds]) for ds in datasets]
        leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=fs(10),
                        handlelength=1.4, borderaxespad=0)
        for t in leg.get_texts():
            t.set_fontweight("bold")

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 6.6))
    print("Labelled, g=4")
    panel(axes[0], trajectories("labelled_filter_fixpoint_rounds.csv", 4, PLOT_IMSM), PLOT_IMSM,
          r"Labelled ($\mathbf{g=4}$; solid to $\mathbf{I_{max}=3}$)", E2E_ROUND)
    mcc = mcc_filtering()
    cheap = mcc[(mcc["mode"] == "IterOrbit") & (mcc.graphlet_size == 5) & mcc.cheap_refuted]
    print("Unlabelled, g=5")
    panel(axes[1], trajectories("unlabelled_filter_fixpoint_rounds.csv", 5, MCC, frozenset(zip(cheap.dataset, cheap.pattern_name))), MCC,
          r"Unlabelled ($\mathbf{g=5}$, to the fixpoint)")
    fig.subplots_adjust(top=0.95, bottom=0.1, left=0.12, right=0.76, hspace=0.38)
    shared_ylabel(fig, axes, "Gain over round 1 (pp)")
    fig.canvas.draw()
    y0 = fig.transFigure.inverted().transform(axes[1].get_tightbbox(fig.canvas.get_renderer()))[0][1]
    fig.text(0.5, y0 - 0.008, "Cumulative filtering cost relative to round 1", ha="center", va="top",
             fontweight="bold", fontsize=fs(12))
    save(fig, "iteration_gain")


fig_iteration_gain()

# %% [markdown]
# ### Figure `refutation.pdf`
# Unsatisfiable unlabelled pairs proved so by filtering alone (LV, Images, Meshes; Scalefree and SI pairs all
# have embeddings).

# %%
def fig_refutation():
    setup(7.2)
    FAMILIES = ["LV", "images-CVIU11", "meshes-CVIU11"]
    d = mcc_filtering()
    gt = csv("unlabelled_glasgow_ground_truth.csv", dtype={"glasgow_status": str})
    unsat = gt[gt.glasgow_status.astype(str).str.lower() == "false"]
    n_unsat = unsat.groupby("dataset").size()
    unsat_keys = set(zip(unsat.dataset, unsat.pattern_name))
    ref = d[d.outcome == "refuted"]
    assert all(k in unsat_keys for k in zip(ref.dataset, ref.pattern_name)), "a refuted pair is satisfiable"
    osi = orbitsi_filter("mccreesh")
    osi_ref = osi[osi.refuted]
    assert all(k in unsat_keys for k in zip(osi_ref.dataset, osi_ref.pattern_name)), "OrbitSI refutes a satisfiable pair"
    print(f"OrbitSI filter finished {len(osi)}/480 pairs")

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.4), sharey=True)
    width = 0.36
    MODES2 = ["OrbitSI", "IterOrbit"]
    for ax, fam in zip(axes, FAMILIES):
        n = int(n_unsat[fam])
        r = ref[ref.dataset == fam]
        cheap = 100.0 * r[r.cheap_refuted].pattern_name.nunique() / n
        print(f"{fam}: {n} unsatisfiable; cheap filters refute {cheap:.0f}%")
        for j, mode in enumerate(MODES2):
            xs = np.arange(len(SIZES)) + (j - 0.5) * width
            vals = [100.0 * r[(r["mode"] == mode) & (r.graphlet_size == g)].pattern_name.nunique() / n for g in SIZES]
            print(f"  {mode:9s} " + " ".join(f"K{g} {v:5.1f}%" for g, v in zip(SIZES, vals)))
            ax.bar(xs, vals, width=width, color=MODE_COLOR[mode], edgecolor="black", linewidth=1.5,
                   hatch=MODE_HATCH[mode], zorder=2)
        ax.axhline(cheap, color="black", linestyle="--", linewidth=1.5, zorder=3)
        o = 100.0 * osi_ref[osi_ref.dataset == fam].pattern_name.nunique() / n
        print(f"  OrbitSI   {o:5.1f}%")
        ax.axhline(o, color=TOOL_COLOR["OrbitSI"], linestyle=":", linewidth=2.2, zorder=3)
        ax.set_xticks(range(len(SIZES)))
        ax.set_xticklabels([f"K{g}" for g in SIZES])
        ax.set_ylim(0, 108)
        style_axis(ax, f"{LABEL[fam]}\n({n} unsat.)")
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    axes[0].set_ylabel("Pairs refuted (%)", fontweight="bold", fontsize=fs(11))
    fig.subplots_adjust(top=0.7, bottom=0.08, left=0.14, right=0.995, wspace=0.12)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=MODE_COLOR[m], edgecolor="black", linewidth=1.5,
                             hatch=MODE_HATCH[m], label=MODE_LABEL_PLAIN[m]) for m in MODES2]
    lines = [plt.Line2D([0], [0], color="black", linestyle="--", linewidth=1.5, label="LDF+NLF+AC"),
             plt.Line2D([0], [0], color=TOOL_COLOR["OrbitSI"], linestyle=":", linewidth=2.2, label="OrbitSI")]
    fig.add_artist(legend_row(fig, "Mode:", handles, 0.965, columnspacing=1.0, handletextpad=0.4, handlelength=1.8))
    legend_row(fig, "Refuted by:", lines, 0.9, columnspacing=1.0, handletextpad=0.4, handlelength=1.8)
    save(fig, "refutation")


fig_refutation()


# %% [markdown]
# ## End-to-end performance (Table `tab:e2e`)
# All eight labelled graphs and the five unlabelled families: IterOrbit, IterOrbit without its search
# accelerations, the OrbitSI package's published configuration, RM, DP-iso, VEQ and CaLiG.

# %%
SIX = PLOT_TOOLS                     # IterOrbit, OrbitSI, RM, DP-iso, VEQ, CaLiG
OTHERS = SIX[1:]
NO_ACC = "IterOrbit (no acc.)"


def e2e_times(bench, datasets):
    """The counting runs of one benchmark: (long frame with t, instance x tool pivot of t), every system plus
    IterOrbit without its search accelerations."""
    loader = (lambda: imsm_phases(datasets)) if bench == "imsm" else mcc_phases
    d = with_index(pd.concat([loader(), orbitsi_phases(bench)], ignore_index=True), bench)
    d = d[d.dataset.isin(datasets)]
    raw = iterorbit_sweep("no_accel", bench)[["dataset", "pattern_name", "status", "filter", "order", "enumerate"]]
    na = raw[(raw.status == "COMPLETED") & raw.dataset.isin(datasets)].assign(tool=NO_ACC, read=0.0, index=0.0)
    na[["order", "enumerate"]] = na[["order", "enumerate"]].fillna(0)
    d = pd.concat([d, na], ignore_index=True)
    d["t"] = phase_time(d)
    return d, d.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")


def e2e_section():
    """Every number of the "End-to-End Performance" subsection: per benchmark, completed, fastest and
    IterOrbit's speed-up over each system; per dataset, the geometric means over the instances all six completed."""

    def speedup(t, a, ref="IterOrbit"):
        k = t[a].notna() & t[ref].notna()
        r = t[a][k] / t[ref][k]
        return f"{gmean(r):.2f}x ({(r > 1).sum()}/{(r < 1).sum()})"

    def report(bench, datasets, keys_total):
        d, t = e2e_times(bench, datasets)
        t = t.reindex(pd.MultiIndex.from_tuples(keys_total, names=["dataset", "pattern_name"]))
        six = t[SIX]
        fastest = six.eq(six.min(axis=1), axis=0).sum()
        print(f"\n==== {BENCH[bench]} ({len(t)} instances)")
        for s in SIX + [NO_ACC]:
            line = f"  {s:20s} done {int(t[s].notna().sum()):3d}"
            if s in SIX:
                line += f"  fastest {int(fastest[s]):3d}"
            if s != "IterOrbit":
                line += f"  speed-up of IterOrbit over it {speedup(t, s)}"
            if s in OTHERS:
                line += f"  | over it, IterOrbit without acc. {speedup(t, s, NO_ACC)}"
            print(line)
        only = t.IterOrbit.notna() & t[OTHERS].isna().all(axis=1)
        print(f"  completed by IterOrbit and none of the other five: {int(only.sum())}; by none of the six: {int(six.isna().all(axis=1).sum())}")
        keys = all_completed(d, SIX)
        print("  per dataset (geo ms over all-six-completed; IterOrbit filter share; fastest and runner-up):")
        for ds in datasets:
            x = d[(d.dataset == ds) & [k in keys for k in zip(d.dataset, d.pattern_name)]]
            g = {s: gmean(x[x.tool == s].t, TIME_FLOOR_MS) for s in SIX}
            order = sorted(SIX, key=g.get)
            io = x[x.tool == "IterOrbit"]
            share = (io["filter"] / phase_time(io)).mean()
            idx = {s: (x[x.tool == s]["index"] / x[x.tool == s].t).mean() for s in OTHERS}
            print(f"    {LABEL[ds]:10s} n={len(io):3d}  " + "  ".join(f"{s} {g[s]:.3g}" for s in SIX)
                  + f"  | fastest {order[0]}, {g[order[1]] / g[order[0]]:.2f}x ahead of {order[1]}"
                  + (f"; IterOrbit {g['IterOrbit'] / g[order[0]]:.2f}x behind" if order[0] != "IterOrbit" else "")
                  + f" | IO filter share {100 * share:.0f}% | index share " + " ".join(f"{s} {100 * v:.0f}%" for s, v in idx.items()))
            f = six.loc[ds]
            print(f"    {'':10s} fastest counts (all {len(f)}): " + " ".join(f"{s} {int(f.eq(f.min(axis=1), axis=0).sum()[s])}" for s in SIX))
        return t

    tabs = {"imsm": report("imsm", IMSM_ALL, instance_keys("imsm")),
            "mccreesh": report("mccreesh", MCC, instance_keys("mccreesh"))}
    s = csv("labelled_sota_runs.csv")
    g = s.groupby("dataset")[["data_vertices", "data_edges", "data_num_labels"]].first()
    g["avg_deg"] = 2 * g.data_edges / g.data_vertices
    print("\nlabelled data graphs:\n" + g.to_string())
    return tabs


e2e_section();

# %% [markdown]
# ### Sensitivity to the index charge (Table `tab:index-amortised`)
# IterOrbit's speed-up and fastest counts with each system's index build charged in full (as in the paper),
# shared by a data graph's 100 patterns, and fully amortised.

# %%
def index_amortised():
    SCALES = (("charged", 1.0), ("shared by 100 patterns", 0.01), ("fully amortised", 0.0))

    def pivot(d0, scale):
        d = d0.assign(index=d0["index"] * scale)
        d["t"] = phase_time(d)
        return d.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first").reindex(columns=SIX)

    def summary(t):
        sp = []
        for s in SIX[1:]:
            k = t[s].notna() & t["IterOrbit"].notna()
            sp.append(f"{SHORT.get(s, s)} {gmean(t[s][k] / t['IterOrbit'][k]):6.2f}x")
        f = t.eq(t.min(axis=1), axis=0).sum()
        return "  ".join(sp) + "  | fastest " + ", ".join(f"{SHORT.get(s, s)} {int(f[s])}" for s in SIX)

    out = {}
    for bench, datasets in (("imsm", IMSM_ALL), ("mccreesh", MCC)):
        d0, _ = e2e_times(bench, datasets)
        d0 = d0[d0.tool.isin(SIX)].copy()
        print(f"\n==== {BENCH[bench]}")
        tabs = {name: pivot(d0, sc) for name, sc in SCALES}
        out[bench] = tabs
        for name, _ in SCALES:
            print(f"  {name:24s} {summary(tabs[name])}")
        if bench != "imsm":
            continue
        rm = d0[d0.tool == "RM"].assign(t=lambda x: x["filter"] + x["order"] + x["enumerate"] + x["index"])
        share = (rm["index"] / rm["t"]).groupby(rm.dataset).median()
        print("  per dataset (index share of RM's time, median):")
        for ds in datasets:
            print(f"   {LABEL[ds]:9s} {100 * share[ds]:3.0f}%")
            for name in ("charged", "fully amortised"):
                print(f"     {name:16s} {summary(tabs[name].loc[ds])}")
    return out


index_amortised();

# %% [markdown]
# ### Figure `cactus.pdf`
# Instances completed within each time, per benchmark; the Glasgow solver on the unlabelled benchmark.

# %%
CACTUS_D = sns.color_palette("Dark2")
CACTUS_COLOR = {"IterOrbit": CACTUS_D[1], "RM": CACTUS_D[0], "VEQ": CACTUS_D[2], "DP-iso": CACTUS_D[3], "CaLiG": CACTUS_D[4],
                "OrbitSI": CACTUS_D[5], "Glasgow": CACTUS_D[7], "GuP": CACTUS_D[6]}
CACTUS_MARKER = {**TOOL_MARKER, "Glasgow": "X", "GuP": "h"}
CACTUS_LINE_WIDTH = {"IterOrbit": 3.4}
CACTUS_DEFAULT_WIDTH = 2.2
CACTUS_N_MARKERS = 8


def _cactus_line(ax, t, y):
    ax.plot(np.arange(1, len(y) + 1), y, color=CACTUS_COLOR[t], marker=CACTUS_MARKER[t],
            markevery=np.unique(np.linspace(0, len(y) - 1, CACTUS_N_MARKERS).astype(int)).tolist(),
            markerfacecolor="none", markeredgecolor=CACTUS_COLOR[t], markeredgewidth=1.5,
            markersize=8, linewidth=CACTUS_LINE_WIDTH.get(t, CACTUS_DEFAULT_WIDTH), linestyle="solid",
            zorder=3 if t == "IterOrbit" else 2, label=SHORT.get(t, t))


def _cactus_axis(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(mticker.FixedLocator([1, 1e2, 1e4, 1e6]))
    ax.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs=(1.0,), numticks=12))
    ax.grid(which="minor", axis="y", linewidth=0.5, alpha=0.6)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(ms_tick))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())


def _cactus_finish(fig, axes, order, name):
    fig.subplots_adjust(top=0.83, bottom=0.095, left=0.13, right=0.975, hspace=0.42)
    for ax in axes:
        style_axis(ax)
    axes[1].set_xlabel("Instances completed", fontweight="bold", fontsize=fs(12))
    handles = [plt.Line2D([0], [0], color=CACTUS_COLOR[t], marker=CACTUS_MARKER[t], markerfacecolor="none",
                          markeredgecolor=CACTUS_COLOR[t], markeredgewidth=1.5, markersize=8,
                          linewidth=CACTUS_LINE_WIDTH.get(t, CACTUS_DEFAULT_WIDTH), label=SHORT.get(t, t)) for t in order]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    left = min(inv.transform(ax.get_tightbbox(r))[0][0] for ax in axes)
    top, bottom = axes[0].get_position(), axes[1].get_position()
    fig.text(left - 0.012, (top.y1 + bottom.y0) / 2, "Time (ms)", rotation=90, ha="right", va="center",
             weight="bold", fontsize=fs(12))
    for row, y in ((handles[:4], 0.976), (handles[4:], 0.926)):
        leg = fig.legend(handles=row, bbox_to_anchor=((top.x0 + top.x1) / 2, y), loc="center",
                         ncol=len(row), frameon=False, borderpad=0, borderaxespad=0,
                         columnspacing=0.8, handletextpad=0.3, handlelength=1.8)
        for text in leg.get_texts():
            text.set_fontweight("bold")
            text.set_fontsize(fs(10))
    save(fig, name)


CACTUS_TITLE = {"imsm": f"Labelled ({len(PLOT_IMSM)} graphs, 800 instances)",
                "mccreesh": f"Unlabelled ({len(MCC)} families, 480 instances)"}


def fig_cactus():
    setup(7.4)
    GLASGOW, FLOOR = "Glasgow", 1.0
    g = glasgow_runs()
    g = g[g.decided & g.index.get_level_values(0).isin(MCC)]
    glasgow = g.runtime.clip(lower=FLOOR)

    def legend_order(mcc):
        order = plot_tool_order()
        t = mcc.assign(t=phase_time(mcc)).pivot_table(index=["dataset", "pattern_name"], columns="tool",
                                                       values="t", aggfunc="first")[order]
        t[GLASGOW] = glasgow.reindex(t.index)
        t = t.dropna().clip(lower=FLOOR)
        gm = {c: gmean(t[c]) for c in t.columns}
        i = next((k for k, c in enumerate(order) if gm[c] < gm[GLASGOW]), len(order))
        print(f"legend: Glasgow placed by {len(t)} pairs all seven completed: " + ", ".join(f"{c} {gm[c]:.3g}" for c in t.columns))
        return order[:i] + [GLASGOW] + order[i:]

    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.2))
    for ax, (bench, loader, datasets) in zip(axes, [("imsm", lambda: imsm_phases(PLOT_IMSM), PLOT_IMSM),
                                                     ("mccreesh", mcc_phases, MCC)]):
        d = with_index(pd.concat([loader(), orbitsi_phases(bench)], ignore_index=True), bench)
        sub = d[d.dataset.isin(datasets)]
        series = {t: np.sort(phase_time(sub[sub.tool == t]).to_numpy()) for t in plot_tool_order()}
        if bench == "mccreesh":
            series[GLASGOW] = np.sort(glasgow.to_numpy())
        done = {}
        for t, y in series.items():
            done[t] = len(y)
            if len(y):
                _cactus_line(ax, t, y)
        _cactus_axis(ax)
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6, integer=True))
        style_axis(ax, CACTUS_TITLE[bench])
        print(f"{BENCH[bench]}: instances completed  " + "  ".join(f"{t} {n}" for t, n in done.items()))
        if bench == "mccreesh":
            order = legend_order(sub)
        ax.set_xlim(0, 500)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(100))
    _cactus_finish(fig, axes, order, "cactus")


fig_cactus()

# %% [markdown]
# ### Figures `e2e_phases.pdf` (counting) and `e2e_phases_enum.pdf` (enumeration)
# Geometric-mean time per dataset and system over the instances all systems completed, with 95% bootstrap
# intervals, each bar divided by the system's mean share of time per phase. Above each panel: the instances on
# which each system is the fastest (leaders bold and underlined); in each x label: the instances it completed
# (the most completed in bold). `e2e_phases_enum.pdf` is not in the paper.

# %%
def fig_e2e_phases(enum=False):
    N_BOOT = 2000
    PHASES = ((["read"] if INCLUDE_LOADING else []) + (["index"] if INCLUDE_INDEX else []) + ["filter", "order", "enumerate"])
    PHASE_LABEL = {"read": "Read", "index": "Index", "filter": "Filter", "order": "Order", "enumerate": "Enumerate"}
    PHASE_COLOR = {"filter": S[0], "order": S[1], "enumerate": S[2], "index": S[6], "read": S[7]}
    PHASE_HATCH = {"filter": "//", "order": "xx", "enumerate": "\\\\", "index": "O", "read": "--"}
    setup(13.5, TEXT_IN)
    underlined = []
    COUNT_SIZE = 9
    COUNT_ROW_PT = fs(9)
    TITLE_PAD_PT = fs(19)
    rng = np.random.default_rng(0)
    GUP = "GuP"
    tools = PLOT_TOOLS + ([GUP] if enum else [])
    stagger = 0.0
    if enum:
        stagger = fs(10)
        TITLE_PAD_PT += stagger
        enum_data = csv("enum_e2e.csv")
        enum_data = enum_data[enum_data.status == "COMPLETED"]

    def bar_order():
        order = plot_tool_order()
        if not enum:
            return order
        x = pd.concat([with_index(enum_data[enum_data.bench == b], b) for b in ("imsm", "mccreesh")])
        x["t"] = phase_time(x)
        t = x.pivot_table(index=["bench", "dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")[order + [GUP]].dropna()
        g = {c: gmean(t[c]) for c in t.columns}
        i = next((k for k, c in enumerate(order) if g[c] < g[GUP]), len(order))
        return order[:i] + [GUP] + order[i:]

    def load(bench):
        if enum:
            return with_index(enum_data[enum_data.bench == bench].copy(), bench)
        loader = (lambda: imsm_phases(PLOT_IMSM)) if bench == "imsm" else mcc_phases
        return with_index(pd.concat([loader(), orbitsi_phases(bench)], ignore_index=True), bench)

    ORDER = bar_order()

    def summarise(d, keys, ds):
        sub = d[d.dataset == ds]
        sub = sub[[k in keys for k in zip(sub.dataset, sub.pattern_name)]]
        out = {}
        for t in tools:
            x = sub[sub.tool == t]
            total = phase_time(x).to_numpy()
            logs = np.log(total)
            boots = rng.choice(logs, size=(N_BOOT, len(logs)), replace=True).mean(axis=1)
            shares = (x[PHASES].clip(lower=0).to_numpy() / total[:, None]).mean(axis=0)
            out[t] = (float(np.exp(logs.mean())), float(np.exp(np.percentile(boots, 2.5))),
                      float(np.exp(np.percentile(boots, 97.5))), shares / shares.sum(), len(x))
        return out

    def instance_counts(d, ds, total):
        x = d[d.dataset == ds].assign(t=lambda f: phase_time(f))
        t = x.pivot_table(index="pattern_name", columns="tool", values="t", aggfunc="first").reindex(columns=tools)
        best = t.min(axis=1)
        fastest = t.eq(best, axis=0).sum()
        done = t.notna().sum()
        assert done.max() <= total
        return {tool: (int(done[tool]), int(fastest[tool])) for tool in tools}

    def leaders(values):
        best = max(values.values())
        return {t for t, v in values.items() if v == best}

    def panel(ax, summ, counts, title):
        most_done = leaders({t: c[0] for t, c in counts.items()})
        most_fastest = leaders({t: c[1] for t, c in counts.items()})
        ymin = 10 ** np.floor(np.log10(min(v[0] for v in summ.values()) / 2))
        for i, t in enumerate(ORDER):
            g, lo, hi, shares, _ = summ[t]
            span, bottom, cum = np.log(g / ymin), ymin, 0.0
            for p, s in zip(PHASES, shares):
                cum += s
                top = ymin * np.exp(span * cum)
                ax.bar(i, top - bottom, 0.8, bottom=bottom, color=PHASE_COLOR[p], edgecolor="black",
                       hatch=PHASE_HATCH[p], linewidth=0)
                if bottom > ymin:
                    ax.hlines(bottom, i - 0.4, i + 0.4, colors="black", linewidth=1.5, zorder=3)
                bottom = top
            ax.bar(i, g - ymin, 0.8, bottom=ymin, fill=False, edgecolor="black", linewidth=1.5)
            ax.errorbar(i, g, yerr=[[g - lo], [hi - g]], ecolor="black", elinewidth=1.5, capsize=3,
                        capthick=1.5, fmt="none")
            lead = t in most_fastest
            a = ax.annotate(str(counts[t][1]), xy=(i, 1), xycoords=("data", "axes fraction"),
                            xytext=(0, COUNT_ROW_PT + (i % 2) * stagger), textcoords="offset points", ha="center",
                            va="center", fontsize=fs(COUNT_SIZE), fontweight="bold" if lead else "normal", annotation_clip=False)
            if lead:
                underlined.append(a)
        ax.set_yscale("log")
        top = max(v[2] for v in summ.values())
        ax.set_ylim(bottom=ymin, top=top * (top / ymin) ** 0.06)
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(ms_tick))
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels([f"{SHORT.get(t, t)} ({counts[t][0]})" for t in ORDER], rotation=90,
                           rotation_mode="anchor", ha="right", va="center")
        ax.tick_params(axis="x", length=0, pad=4)
        style_axis(ax)
        ax.set_title(title, fontweight="bold", fontsize=fs(11), pad=TITLE_PAD_PT)
        for t, label in zip(ORDER, ax.get_xticklabels()):
            label.set_fontweight("bold" if t in most_done else "normal")

    fig, lab, unl = bench_grid((13.5, 12.4 if enum else 11.8), top=0.86 if enum else 0.875, bottom=0.11, left=0.075,
                               right=0.995, hspace=1.6 if enum else 1.42, wspace=0.32, block_gap=0.03)
    for bench, datasets, axes in [("imsm", PLOT_IMSM, lab), ("mccreesh", MCC, unl)]:
        d = load(bench)
        keys = all_completed(d, tools)
        sel = csv("labelled_patterns.csv") if bench == "imsm" else csv("unlabelled_pairs.csv")
        totals = sel.groupby("dataset").size()
        print(f"\n{BENCH[bench]}: geometric-mean time (ms) and shares ({'/'.join(PHASES)}) over instances all {len(tools)} completed")
        for ax, ds in zip(axes, datasets):
            summ = summarise(d, keys, ds)
            n = summ[tools[0]][4]
            print(f"  {ds:14s} n={n:3d}  " + "  ".join(f"{t} {summ[t][0]:.4g} [" + "/".join(f"{100 * s:.0f}" for s in summ[t][3]) + "]"
                                                   for t in tools))
            counts = instance_counts(d, ds, int(totals[ds]))
            print(f"  {'':14s} of {int(totals[ds])}: " + "  ".join(f"{SHORT.get(t, t)} done {c[0]} fastest {c[1]}" for t, c in counts.items()))
            panel(ax, summ, counts, f"{LABEL[ds]} (n={n})")

    for row in (lab[:4], lab[4:], unl):
        equal_gaps(fig, row, left=0.035, right=0.995, gap=0.016)
    block_heading(fig, lab, "Labelled")
    block_heading(fig, unl, "Unlabelled")
    shared_ylabel(fig, lab + unl, "Time (ms)")
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=PHASE_COLOR[p], edgecolor="black", linewidth=1.5,
                             hatch=PHASE_HATCH[p], label=PHASE_LABEL[p]) for p in PHASES]
    legend_row(fig, "Phase:", handles, 0.982, columnspacing=1.4, handletextpad=0.5, handlelength=2.0)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for a in underlined:
        (x0, y0), (x1, _) = inv.transform(a.get_window_extent(r))
        fig.add_artist(plt.Line2D([x0, x1], [y0 - 0.001] * 2, transform=fig.transFigure, color="black", linewidth=1.2))
    save(fig, "e2e_phases_enum" if enum else "e2e_phases")


fig_e2e_phases()

# %% [markdown]
# ### The Glasgow solver (unlabelled)
# Glasgow's runtime starts after both graphs are read and is reported in whole milliseconds, so both systems are
# floored at 1 ms here.

# %%
def glasgow_section():
    FLOOR = 1.0

    def gm(x):
        return float(np.exp(np.log(np.asarray(x, dtype=float)).mean()))

    g = glasgow_runs()
    g = pd.DataFrame({"g_done": g.decided, "g_ms": g.runtime, "g_wall": g.wall})
    ph = mcc_phases()
    ox = ph[ph.tool == "IterOrbit"].set_index(["dataset", "pattern_name"])
    e = csv("unlabelled_e2e_results.csv")
    od = e[e.tool == "OrbitEX"].set_index(["dataset", "pattern_name"])
    o_done = od.done.astype(bool)
    g_done = g.g_done.reindex(od.index).fillna(False).astype(bool)
    print(f"completed: IterOrbit {int(o_done.sum())}, Glasgow {int(g_done.sum())}; "
          f"only IterOrbit {int((o_done & ~g_done).sum())}, only Glasgow {int((~o_done & g_done).sum())}")
    j = pd.concat([g, (ox["filter"] + ox["order"] + ox["enumerate"]).rename("o_ms"),
                   od.wall_ms.where(o_done).rename("o_wall")], axis=1).dropna(subset=["g_ms", "o_ms"])
    r = j.g_ms.clip(lower=FLOOR) / j.o_ms.clip(lower=FLOOR)
    print(f"both completed {len(j)}; excl. loading (floor {FLOOR} ms): {gm(r):.2f}x, "
          f"IterOrbit wins {(r > 1).sum()}, losses {(r < 1).sum()}, ties {(r == 1).sum()}")
    t = r.groupby(level=0).agg(n="size", speedup=gm, wins=lambda x: (x > 1).sum(), losses=lambda x: (x < 1).sum(),
                               ties=lambda x: (x == 1).sum())
    t["glasgow_ms"] = j.g_ms.clip(lower=FLOOR).groupby(level=0).agg(gm)
    t["iterorbit_ms"] = j.o_ms.clip(lower=FLOOR).groupby(level=0).agg(gm)
    print(t.round(2).to_string())
    long_ = (j.g_ms >= 1000) | (j.o_ms >= 1000)
    rl = r[long_]
    print(f"pairs either takes >= 1 s: {int(long_.sum())}; {gm(rl):.2f}x, IterOrbit wins {(rl > 1).sum()}, "
          f"losses {(rl < 1).sum()}; by family " + ", ".join(f"{d} {len(x)}: {gm(x):.2f}x" for d, x in rl.groupby(level=0)))
    w = j.dropna(subset=["o_wall"])
    rw = w.g_wall / w.o_wall
    print(f"wall clock: {gm(rw):.2f}x, wins {(rw > 1).sum()}, losses {(rw < 1).sum()}")
    print("\nfamily         pairs  done_IO  done_G  both  excl   W/L/T        wall")
    for fam in MCC + ["All"]:
        sel = (lambda ix: ix.get_level_values(0) == fam) if fam != "All" else (lambda ix: np.ones(len(ix), bool))
        o, gd, rr, ww = o_done[sel(o_done.index)], g_done[sel(g_done.index)], r[sel(r.index)], rw[sel(rw.index)]
        print(f"{fam:14s} {len(o):5d} {int(o.sum()):8d} {int(gd.sum()):7d} {len(rr):5d} {gm(rr):5.2f}x "
              f"{(rr > 1).sum():3d}/{(rr < 1).sum():3d}/{(rr == 1).sum():3d} {gm(ww):6.2f}x")
    return t


glasgow_section();

# %% [markdown]
# ### Induced subgraph isomorphism (unlabelled)
# IterOrbit (g=5 fixpoint, `--count-only --induced`) against the Glasgow solver (`--count-solutions --induced`),
# all 480 pairs; times floored at 1 ms. IterOrbit's time is filter + order + search, and also with its symmetry
# timer added.

# %%
def induced_section():
    FLOOR = 1.0
    d = csv("unlabelled_induced_e2e.csv")
    io = d[d.tool == "IterOrbit"].set_index(["dataset", "pattern_name"])
    gl = d[d.tool == "Glasgow"].set_index(["dataset", "pattern_name"])
    print(f"runs: IterOrbit {len(io)}, Glasgow {len(gl)}")
    iod, gld = io[io.status == "COMPLETED"], gl[gl.status == "COMPLETED"]
    print(f"completed: IterOrbit {len(iod)}, Glasgow {len(gld)}; only IterOrbit {len(iod.index.difference(gld.index))}, "
          f"only Glasgow {len(gld.index.difference(iod.index))}, neither {480 - len(iod.index.union(gld.index))}")
    io_t = (iod.filter_ms + iod.order_ms.fillna(0) + iod.search_ms.fillna(0)).clip(lower=FLOOR)
    io_ts = (io_t + iod.symmetry_ms.fillna(0)).clip(lower=FLOOR)
    gl_t = gld.runtime_ms.clip(lower=FLOOR)
    j = pd.DataFrame({"io": iod["count"], "sat": iod.saturated, "gl": gld["count"]}).dropna(subset=["io", "gl"])
    bad = j[(j.io != j.gl) & ~j.sat.astype(bool)]
    print(f"counts compared on {len(j)} pairs: {len(bad)} differ (excluding {int(j.sat.astype(bool).sum())} saturated at 2^63)")
    if len(bad):
        print(bad.head(10).to_string())
    unsat = gld[gld.satisfiable == False].index   # noqa: E712
    ref = iod[iod.search_ms.isna()].index
    print(f"unsatisfiable (Glasgow): {len(unsat)}; IterOrbit refutes {len(ref.intersection(unsat))} of them by filtering alone"
          f" (and {len(ref.difference(unsat))} pairs Glasgow does not call unsatisfiable)")
    both = io_t.index.intersection(gl_t.index)
    r, rs = gl_t[both] / io_t[both], gl_t[both] / io_ts[both]
    print(f"\nboth complete: {len(both)}; Glasgow/IterOrbit geo {gmean(r):.2f}x (with IterOrbit's symmetry timer "
          f"{gmean(rs):.2f}x); IterOrbit faster on {(r > 1).sum()}, Glasgow on {(r < 1).sum()}, tied {(r == 1).sum()}")
    print(f"  IterOrbit symmetry timer: median {iod.symmetry_ms.median():.2g} ms, max {iod.symmetry_ms.max():.3g} ms")
    print("per family (n both, geo speed-up, IterOrbit done / Glasgow done of the family):")
    for ds in MCC:
        b = [k for k in both if k[0] == ds]
        n_all = (d[(d.tool == "IterOrbit") & (d.dataset == ds)]).shape[0]
        print(f"  {LABEL[ds]:10s} n={len(b):3d}  {gmean(r[b]):6.2f}x   IO {(iod.index.get_level_values(0) == ds).sum():3d} / "
              f"Glasgow {(gld.index.get_level_values(0) == ds).sum():3d} of {n_all}")
    ni = mcc_phases()
    ni = ni[ni.tool == "IterOrbit"]
    ni_t = pd.Series(phase_time(ni).clip(lower=FLOOR).to_numpy(), index=pd.MultiIndex.from_arrays([ni.dataset, ni.pattern_name]))
    k = ni_t.index.intersection(io_t.index)
    print(f"\nIterOrbit induced vs non-induced (paper's runs, different build), {len(k)} pairs both complete: "
          f"induced/non-induced time geo {gmean(io_t[k] / ni_t[k]):.2f}; non-induced completed {len(ni_t)}, induced {len(io_t)}")


induced_section()

# %% [markdown]
# ## What the orbit stage contributes (Table `tab:orbit-ablation`)
# 1. The orbit-free fixpoint against the fixpoint with orbits (labelled g=4, unlabelled g=5), both to convergence.
# 2. The search accelerations end to end (IterOrbit vs its no-acceleration rerun).
# 3. The one-at-a-time acceleration ablation (process wall-clock times, timeout = 1800 s; HPRD left out).

# %%
def orbit_section():
    def ablation(fixpoint_csv, no_orbit_csv, datasets, g):
        fx = csv(fixpoint_csv)
        fx = fx[(fx.graphlet_size == g) & fx.dataset.isin(datasets)]
        no = csv(no_orbit_csv)
        no = no[no.dataset.isin(datasets)]
        if "outcome" not in fx:
            fx = fx.assign(outcome="completed")
        m = fx.merge(no, on=["dataset", "pattern_name"], suffixes=("_o", "_n"))
        m = m[m.reduction_vs_semijoin_o.notna() & m.reduction_vs_semijoin_n.notna()]
        out = m.groupby("dataset").agg(
            n=("pattern_name", "size"),
            orb=("reduction_vs_semijoin_o", lambda r: 100 * r.mean()), none=("reduction_vs_semijoin_n", lambda r: 100 * r.mean()),
            ref_o=("outcome_o", lambda o: int((o == "refuted").sum())), ref_n=("outcome_n", lambda o: int((o == "refuted").sum())),
            t_o=("total_filter_ms_o", gmean), t_n=("total_filter_ms_n", gmean))
        out["gap"] = out.orb - out["none"]
        out["t_ratio"] = out.t_n / out.t_o
        tighter = (m.reduction_vs_semijoin_n > m.reduction_vs_semijoin_o + 1e-9).sum()
        print(f"\n== orbit ablation g={g}: none tighter than orbits on {tighter} of {len(m)}")
        print(out.round(2).rename(index=LABEL).to_string())
        return out

    lab = ablation("labelled_filter_fixpoint.csv", "labelled_filter_no_orbit.csv", PLOT_IMSM, 4)
    unl = ablation("unlabelled_filter_fixpoint.csv", "unlabelled_filter_no_orbit.csv", MCC, 5)
    rich = ["dblp", "eu2005", "hprd", "human"]
    print(f"label-rich {rich}: max gap {lab.loc[rich].gap.max():.2f} pp, none/orbits time {lab.loc[rich].t_ratio.min():.2f}-{lab.loc[rich].t_ratio.max():.2f}")

    print("\n== search accelerations end to end")
    for bench, datasets in (("imsm", IMSM_ALL), ("mccreesh", MCC)):
        d, t = e2e_times(bench, datasets)
        done = t.notna().sum()
        k = t.IterOrbit.notna() & t[NO_ACC].notna()
        r = t[NO_ACC][k] / t.IterOrbit[k]
        print(f"  {BENCH[bench]}: done IO {done['IterOrbit']} no-acc {done[NO_ACC]}; IO over no-acc {gmean(r):.2f}x on {k.sum()}")
        per = pd.Series({ds: gmean(r.xs(ds, level='dataset')) for ds in datasets if ds in r.index.get_level_values(0)})
        print("    per dataset: " + "  ".join(f"{LABEL[ds]} {v:.2f}x" for ds, v in per.sort_values(ascending=False).items()))
        best = t[OTHERS].min(axis=1)
        kk = t[NO_ACC].notna() & best.notna()
        print(f"    no-acc vs fastest other per instance: {gmean(best[kk] / t[NO_ACC][kk]):.2f}x on {kk.sum()};"
              f" most completed by another system: {t[OTHERS].notna().sum().idxmax()} {t[OTHERS].notna().sum().max()}")
        print("    no-acc fastest-of-six-with-no-acc-in-place-of-IO: " + str(int(t[[NO_ACC] + OTHERS].eq(t[[NO_ACC] + OTHERS].min(axis=1), axis=0).sum()[NO_ACC])))

    print("\n== one-at-a-time ablation (wall clock, timeout 1800 s), HPRD left out")
    a = csv("accel_ablation.csv")
    x = a[a.dataset != "hprd"]
    g = lambda v: float(np.exp(np.log(v).mean()))
    print(f"  {len(x)} instances: {x.groupby('bench').size().to_dict()}; per dataset {x.groupby('dataset').size().to_dict()}")
    print("  completed:", {c: int(x[f'{c}_done'].sum()) for c in ('old', 'edge', 'sym', 'new')})
    for name, s in [("all", x)] + list(x.groupby("bench")):
        print(f"  {name}: speed-up over both off  edge {g(s.old_eff / s.edge_eff):.2f}  sym {g(s.old_eff / s.sym_eff):.2f}  both {g(s.old_eff / s.new_eff):.2f}")
    sf = x.dataset == "scalefree"
    for label, s in (("scalefree", x[sf]), ("not scalefree", x[~sf])):
        print(f"  {label} ({len(s)}): edge-only/both-off {g(s.edge_eff / s.old_eff):.3f}; sym-only/both-off {g(s.sym_eff / s.old_eff):.3f}")
    print("  symmetric (multiplier > 1):", int((x.multiplier > 1).sum()), "of", len(x), "; median multiplier", float(x[x.multiplier > 1].multiplier.median()))
    return lab, unl


orbit_section();

# %% [markdown]
# ### Figure `sota_filter.pdf`
# Filtering (+ index) time vs candidates per pattern vertex: IterOrbit's filter in its end-to-end configuration
# against OrbitSI's and the four baselines'; markers are the geometric-mean time and mean candidates, lines their
# 95% confidence intervals.

# %%
def fig_sota_filter():
    import math
    setup(13.5, TEXT_IN)
    CI_ALPHA, CI_WIDTH, CI_CAP = 0.7, 1.3, 3
    Z = 1.96
    ORDER = plot_tool_order()

    def orbitsi(bench):
        o = orbitsi_filter(bench)
        f = csv("labelled_filter_singlepass.csv" if bench == "imsm" else "unlabelled_filter_singlepass.csv")
        f = f[(f.graphlet_size == 4) & (f.use_full_graph == 0)]
        cand = f.avg_cand_after_orbit.where(f.outcome != "refuted", 0.0) if "outcome" in f else f.avg_cand_after_orbit
        cand = dict(zip(zip(f.dataset, f.pattern_name), cand))
        return pd.DataFrame({"dataset": o.dataset, "pattern_name": o.pattern_name, "method": "OrbitSI",
                             "cost_ms": o.cost_ms, "cand": [cand.get(k, np.nan) for k in zip(o.dataset, o.pattern_name)]})

    def labelled():
        r = csv("labelled_filter_fixpoint_rounds.csv")
        r = r[(r.graphlet_size == 4) & (r["round"] <= 3)].sort_values("round")
        ox = r.groupby(["dataset", "pattern_name"], as_index=False).last()
        ox = pd.DataFrame({"dataset": ox.dataset, "pattern_name": ox.pattern_name, "method": "IterOrbit",
                           "cost_ms": ox.cumulative_cost_ms, "cand": ox.avg_candidates})
        s = csv("labelled_sota_runs.csv")
        s = s[s.status == "COMPLETED"]
        s = pd.DataFrame({"dataset": s.dataset, "pattern_name": s.pattern_name,
                          "method": s.tool.map({"rapidmatch": "RM", "calig": "CaLiG", "veq": "VEQ", "dpiso": "DP-iso"}),
                          "cost_ms": s.filter_vertices_time_s * 1000.0, "cand": s.num_candidates_avg})
        return pd.concat([ox, orbitsi("imsm"), s], ignore_index=True), PLOT_IMSM

    def unlabelled():
        f = csv("unlabelled_filter_fixpoint.csv")
        f = f[f.graphlet_size == 5]
        ox = pd.DataFrame({"dataset": f.dataset, "pattern_name": f.pattern_name, "method": "IterOrbit",
                           "cost_ms": f.total_filter_ms,
                           "cand": np.where(f.outcome == "refuted", 0.0, f.avg_cand_after_orbit_final)})
        s = csv("unlabelled_sota_filter_candidates.csv")
        s = s[s.status == "COMPLETED"]
        s = pd.DataFrame({"dataset": s.dataset, "pattern_name": s.pattern_name, "method": s.tool,
                          "cost_ms": s.filter_ms, "cand": s.avg_candidates})
        return pd.concat([ox, orbitsi("mccreesh"), s], ignore_index=True), MCC

    def restrict(d):
        d = d.dropna(subset=["cost_ms", "cand"])
        o = d[d.method == "OrbitSI"]
        osi_keys = set(zip(o.dataset, o.pattern_name))
        d = d[[k in osi_keys for k in zip(d.dataset, d.pattern_name)]]
        base = d[~d.method.isin(["IterOrbit", "OrbitSI"])]
        nb = base.method.nunique()
        k = base.groupby(["dataset", "pattern_name"]).method.nunique()
        keys = set(k[k == nb].index)
        return d[[x in keys for x in zip(d.dataset, d.pattern_name)]]

    def yfmt(v, _):
        return f"{v / 1000:g}k" if v >= 1000 else f"{v:g}"

    def xfmt(v, _):
        return rf"$\mathbf{{10^{{{int(round(math.log10(v)))}}}}}$" if v >= 100 else f"{v:g}"

    fig, lab, unl = bench_grid((13.5, 9.6), top=0.88, bottom=0.085, left=0.075, right=0.995, hspace=0.72, wspace=0.28, block_gap=0.045)
    for row, (bench, loader, axes) in enumerate((("imsm", labelled, lab), ("mccreesh", unlabelled, unl))):
        d, datasets = loader()
        d = restrict(d[d.dataset.isin(datasets)])
        if INCLUDE_INDEX:
            d = with_index(d.assign(tool=d.method), bench).drop(columns="tool")
            d["cost_ms"] = d.cost_ms + d["index"]
            print(f"\n{bench}: mean index share of the plotted time: " + "  ".join(
                f"{t} {100 * (g['index'] / g.cost_ms).mean():.0f}%" for t, g in d.groupby("method")))
        pooled = d.groupby("method").agg(n=("cand", "size"), cost=("cost_ms", lambda c: gmean(c, 1e-3)),
                                         cand=("cand", "mean")).reindex(ORDER).dropna()
        print(f"\n{['Labelled', 'Unlabelled'][row]}: pooled (n, geo ms, mean cand/vertex)\n{pooled.round(2).to_string()}")
        for ax, ds in zip(axes, datasets):
            x = d[d.dataset == ds]
            m = x.groupby("method").agg(
                cost=("cost_ms", lambda c: gmean(c, 1e-3)), cand=("cand", "mean"),
                lse=("cost_ms", lambda c: Z * np.log(np.maximum(c, 1e-3)).std(ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0.0),
                kse=("cand", lambda c: Z * c.std(ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0.0),
                refuted=("cand", lambda c: int((c == 0).sum()))).reindex(ORDER).dropna(subset=["cost"])
            m["c1"], m["c3"] = m.cost / np.exp(m.lse), m.cost * np.exp(m.lse)
            m["k1"], m["k3"] = (m.cand - m.kse).clip(lower=0), m.cand + m.kse
            n = x[x.method == "RM"].pattern_name.nunique()
            print(f"  {ds:14s} n={n:3d}  " + "  ".join(f"{t} {m.loc[t, 'cost']:.3g}ms/{m.loc[t, 'cand']:.3g}"
                                                       f"(ref {m.loc[t, 'refuted']:.0f})" for t in m.index))
            ax.set_xscale("log")
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(xfmt))
            xmin, xmax = max(x.cost_ms.min(), 0.05), x.cost_ms.max()
            x0, x1 = xmin / 1.3, xmax * (xmax / xmin) ** 0.15
            ax.set_xlim(x0, x1)
            top = max(m.cand.max(), m.k3.max())
            y1 = top * 1.15 if top > 0 else 1
            ax.set_ylim(0, y1)
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(yfmt))
            for t, r in m.iterrows():
                col = TOOL_COLOR[t]
                ax.errorbar(r.cost, r.cand, xerr=[[r.cost - r.c1], [r.c3 - r.cost]], yerr=[[r.cand - r.k1], [r.k3 - r.cand]],
                            fmt="none", ecolor=col, elinewidth=CI_WIDTH, capsize=CI_CAP, capthick=CI_WIDTH, alpha=CI_ALPHA, zorder=2)
                ax.plot(r.cost, r.cand, marker=TOOL_MARKER[t], markersize=10, markerfacecolor="none",
                        markeredgecolor=col, markeredgewidth=1.8, linestyle="none", zorder=4 if t == "IterOrbit" else 3)
            style_axis(ax, f"{LABEL[ds]} (n={n})")

    for row in (lab[:4], lab[4:], unl):
        equal_gaps(fig, row, left=0.045, right=0.995, gap=0.022)
    block_heading(fig, lab, "Labelled")
    block_heading(fig, unl, "Unlabelled")
    shared_ylabel(fig, lab + unl, "Candidates / vertex")
    shared_xlabel(fig, lab + unl, "Filtering + Index Time (ms, log scale)" if INCLUDE_INDEX else "Filtering Time (ms, log scale)")
    handles = [plt.Line2D([0], [0], marker=TOOL_MARKER[t], linestyle="none", markersize=9, markerfacecolor="none",
                          markeredgecolor=TOOL_COLOR[t], markeredgewidth=1.8, label=SHORT.get(t, t)) for t in ORDER]
    legend_row(fig, "Filter:", handles, 0.975, columnspacing=1.2, handletextpad=0.3)
    save(fig, "sota_filter")


fig_sota_filter()


# %% [markdown]
# ## Enumeration (first $10^5$ embeddings, one at a time)
# Every system visits the embeddings one at a time up to $10^5$ (Sun and Luo's protocol), without symmetry
# breaking or last-level counting; GuP joins as a seventh system (unlabelled patterns of at most 256 vertices).

# %%
SEVEN = SIX + ["GuP"]
ENUM_LIMIT = 100000


def enum_section():
    """Count checks against the counting runs; per benchmark, completed, fastest and IterOrbit's speed-up over
    each system, next to the same numbers when counting; per dataset, the speed-ups."""
    e = csv("enum_e2e.csv")
    ref_l = csv("labelled_e2e_per_pattern.csv")[["dataset", "pattern_name", "ox_emb"]]
    ref_m = csv("unlabelled_e2e_results.csv")
    ref_m = ref_m[ref_m.tool == "OrbitEX"][["dataset", "pattern_name", "count"]].rename(columns={"count": "ox_emb"})
    ref = pd.concat([ref_l, ref_m]).dropna()
    c = e[e.status == "COMPLETED"].merge(ref, on=["dataset", "pattern_name"], how="left")
    c["expect"] = np.minimum(c.ox_emb, ENUM_LIMIT)
    ok = np.where(c.ox_emb >= ENUM_LIMIT, c["count"] >= ENUM_LIMIT, c["count"] == c.ox_emb)
    bad = c[c.ox_emb.notna() & ~ok]
    over = c[c["count"] > ENUM_LIMIT]
    print(f"count checks: {int(c.ox_emb.notna().sum())} completed runs with a reference, {len(bad)} wrong; "
          f"over the limit: {over.groupby(['bench', 'tool']).size().to_dict()} (max {int(over['count'].max()) if len(over) else 0})")
    if len(bad):
        print(bad.groupby(["bench", "tool"]).size().to_string())
        print(bad[["bench", "tool", "dataset", "pattern_name", "count", "expect"]].head(20).to_string(index=False))
    nr = c[c.ox_emb.isna()].assign(capped=lambda f: np.minimum(f["count"], ENUM_LIMIT))
    spread = nr.groupby(["bench", "dataset", "pattern_name"]).capped.agg(["min", "max", "size"])
    dis = spread[spread["min"] != spread["max"]]
    print(f"without a reference: {len(nr)} completed runs on {len(spread)} instances "
          f"({int((spread['size'] > 1).sum())} completed by two or more systems), {len(dis)} with disagreeing capped counts")
    if len(dis):
        print(dis.head(20).to_string())
    print(f"runs reaching the limit: " + e[e.status == "COMPLETED"].assign(r=lambda f: f.limit_reached.astype(float)).groupby(["bench", "tool"]).r.mean().round(2).to_dict().__repr__())

    def speedup(t, a, ref="IterOrbit"):
        k = t[a].notna() & t[ref].notna()
        r = t[a][k] / t[ref][k]
        return (gmean(r) if k.any() else np.nan), int(k.sum())

    for bench, datasets in (("imsm", IMSM_ALL), ("mccreesh", MCC)):
        x = e[(e.bench == bench) & (e.status == "COMPLETED") & e.dataset.isin(datasets)].copy()
        x = with_index(x, bench)
        x["t"] = phase_time(x)
        t = x.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first").reindex(columns=SEVEN)
        _, tc = e2e_times(bench, datasets)
        tc = tc.reindex(columns=SEVEN)
        print(f"\n==== {BENCH[bench]}: enumeration (10^5) vs counting (tab:e2e)")
        keys = t.index.union(tc.index)
        t, tc = t.reindex(keys), tc.reindex(keys)
        fe = t.eq(t.min(axis=1), axis=0).sum()
        fe6 = t[SIX].eq(t[SIX].min(axis=1), axis=0).sum()
        fc = tc.eq(tc.min(axis=1), axis=0).sum()
        print(f"  fastest of the six without GuP: " + ", ".join(f"{q} {int(fe6[q])}" for q in SIX))
        for s in SEVEN:
            line = f"  {s:10s} done {int(t[s].notna().sum()):3d} (count {int(tc[s].notna().sum()):3d})  fastest {int(fe[s]):3d} (count {int(fc[s]):3d})"
            if s != "IterOrbit":
                ge, ne = speedup(t, s)
                gc, nc = speedup(tc, s)
                line += f"  IO speed-up {ge:.2f}x on {ne} (count {gc:.2f}x on {nc})"
            print(line)
        g = x[x.tool == "GuP"]
        if len(g):
            gi = x[x.tool.isin(["GuP", "IterOrbit"])].assign(t=lambda f: np.where(f.tool == "GuP", (f.t - f["index"]).clip(lower=TIME_FLOOR_MS), f.t))
            ts = gi.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")
            gs, ns = speedup(ts, "GuP")
            print(f"  GuP: index build median {g['index'].median():.1f} ms, {100 * (g['index'] / g.t).median():.0f}% of its time; "
                  f"IO speed-up over GuP without its index {gs:.2f}x on {ns}")
        io = x[x.tool == "IterOrbit"]
        print(f"  IterOrbit filter share of its time: median {100 * (io['filter'] / io.t).median():.0f}%")
        print("  per dataset, IO speed-up over each (enumeration | counting):")
        for ds in datasets:
            te, tcc = t.loc[ds], tc.loc[ds]
            print(f"    {LABEL[ds]:10s} " + "  ".join(f"{SHORT.get(s, s)} {speedup(te, s)[0]:.2f}|{speedup(tcc, s)[0]:.2f}" for s in SEVEN[1:])
                  + f"  | fastest IO {int(te.eq(te.min(axis=1), axis=0).sum()['IterOrbit'])}/{int(tcc.eq(tcc.min(axis=1), axis=0).sum()['IterOrbit'])}")


enum_section()

# %% [markdown]
# ### By pattern size and density (Table `tab:size-density` from the counting runs)
# Per bin: instances any system completed, how many each system completed and won (fastest of six), and
# IterOrbit's speed-up over each system. Dense = average degree >= 3 (Sun and Luo's split; the labelled corpus's
# own names are kept, which call every 4-vertex pattern dense).

# %%
def size_density(mode="count"):
    assert mode in ("enum", "count")
    SIZE_BINS = {"imsm": [0, 8, 16, 24, 32], "mccreesh": [0, 16, 32, 128, 10 ** 6]}
    e = csv("enum_e2e.csv")
    lab = csv("labelled_patterns.csv").assign(bench="imsm", nv=lambda d: d.vertex_count, nE=lambda d: d.edge_count)
    unl = csv("unlabelled_pairs.csv")
    unl = unl.assign(bench="mccreesh", pattern_name=unl.instance, nv=unl.pattern_n, nE=unl.pattern_m)
    cols = ["bench", "dataset", "pattern_name", "nv", "nE"]
    inst = pd.concat([lab[cols], unl[cols]], ignore_index=True)
    inst["dense"] = 2 * inst.nE / inst.nv >= 3
    ok = inst.bench == "imsm"
    inst.loc[ok, "dense"] = inst[ok].pattern_name.str.contains("dense")

    def speedup(t, a):
        k = t[a].notna() & t["IterOrbit"].notna()
        return (gmean(t[a][k] / t["IterOrbit"][k]) if k.any() else np.nan), int(k.sum())

    def table(t, groups):
        fast = t.eq(t.min(axis=1), axis=0)
        for name, keys in groups:
            tt, ff = t.loc[keys], fast.loc[keys]
            line = f"  {name:12s} n={len(tt):3d}  "
            line += "  ".join(f"{SHORT.get(s, s)} {int(tt[s].notna().sum()):3d}/{int(ff[s].sum()):3d}" for s in SIX)
            line += "  | IO speed-up " + "  ".join(
                f"{SHORT.get(s, s)} {g:5.2f}x({n})" for s in SIX[1:] for g, n in [speedup(tt, s)])
            print(line)

    for bench, datasets in (("imsm", IMSM_ALL), ("mccreesh", MCC)):
        if mode == "count":
            _, t = e2e_times(bench, datasets)
            t = t.reindex(columns=SIX).dropna(how="all")
        else:
            x = e[(e.bench == bench) & (e.status == "COMPLETED") & e.dataset.isin(datasets)].copy()
            x = with_index(x, bench)
            x["t"] = phase_time(x)
            t = x.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first").reindex(columns=SIX)
        meta = inst[inst.bench == bench].set_index(["dataset", "pattern_name"]).reindex(t.index)
        assert meta.nv.notna().all()
        size = pd.cut(meta.nv, SIZE_BINS[bench])
        print(f"\n==== {BENCH[bench]} ({mode}): completed/fastest per system, IterOrbit's speed-up (instances both completed)")
        print(" by |V(q)|:")
        table(t, [(f"{iv.left + 1}-{iv.right}" if iv.right < 10 ** 6 else f">{iv.left}", t.index[size == iv])
                  for iv in size.cat.categories])
        print(" by density:")
        table(t, [(n, t.index[meta.dense == v]) for n, v in (("sparse", False), ("dense", True))])
        print(" by |V(q)| x density:")
        table(t, [(f"{'D' if v else 'S'} {iv.left + 1}-{iv.right}" if iv.right < 10 ** 6 else f"{'D' if v else 'S'} >{iv.left}",
                   t.index[(size == iv) & (meta.dense == v)]) for iv in size.cat.categories for v in (False, True)])
        if mode == "count":
            continue
        lim = x.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="limit_reached", aggfunc="first")
        r = lim.reindex(t.index)["IterOrbit"].astype(float)
        print(" IterOrbit runs reaching 10^5, by |V(q)|: " + ", ".join(
            f"{iv}: {r[size == iv].mean():.2f}" for iv in size.cat.categories))


size_density("count")
size_density("enum")

# %% [markdown]
# ### Figure `enum_cactus.pdf`
# The cactus plot of the enumeration sweep (both panels to the full 800 / 480 instances; GuP placed in the
# legend by its geometric mean over the instances all seven completed).

# %%
def fig_enum_cactus():
    setup(7.4)
    GUP = "GuP"
    X = {"imsm": (800, 200), "mccreesh": (480, 100)}
    e = csv("enum_e2e.csv")

    def legend_order():
        order = plot_tool_order()
        x = e[(e.status == "COMPLETED")].copy()
        x = pd.concat([with_index(x[x.bench == b], b) for b in ("imsm", "mccreesh")])
        x["t"] = phase_time(x)
        t = x.pivot_table(index=["bench", "dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")[order + [GUP]].dropna()
        g = {c: gmean(t[c]) for c in t.columns}
        i = next((k for k, c in enumerate(order) if g[c] < g[GUP]), len(order))
        print(f"legend: GuP placed by {len(t)} instances all seven completed: " + ", ".join(f"{c} {g[c]:.3g}" for c in t.columns))
        return order[:i] + [GUP] + order[i:]

    order = legend_order()
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.2))
    for ax, (bench, datasets) in zip(axes, [("imsm", PLOT_IMSM), ("mccreesh", MCC)]):
        d = e[(e.bench == bench) & (e.status == "COMPLETED") & e.dataset.isin(datasets)].copy()
        d = with_index(d, bench)
        done = {}
        for t in order:
            y = np.sort(phase_time(d[d.tool == t]).to_numpy())
            done[t] = len(y)
            if len(y):
                _cactus_line(ax, t, y)
        print(f"{BENCH[bench]}: instances completed  " + "  ".join(f"{t} {n}" for t, n in done.items()))
        _cactus_axis(ax)
        ax.set_xlim(0, X[bench][0])
        ax.xaxis.set_major_locator(mticker.MultipleLocator(X[bench][1]))
        style_axis(ax, CACTUS_TITLE[bench])
    _cactus_finish(fig, axes, order, "enum_cactus")


fig_enum_cactus()

# %% [markdown]
# ### Figure `e2e_phases_enum.pdf` (not in the paper)
# `e2e_phases.pdf` for the enumeration sweep, with GuP as a seventh bar (its single timer counts as enumeration,
# plus its index build); the unlabelled bars cover only the pairs GuP accepts.

# %%
fig_e2e_phases(enum=True)


# %% [markdown]
# # Exploratory analyses
# Analyses behind decisions and checks that are not, or no longer, in the paper.
#
# ## IterOrbit without its search accelerations, and the one-at-a-time ablation
# The no-acceleration rerun (`--no-symmetry --no-edge-orbit-filter`, same frozen binary) against the default run
# and the four baselines; then, on 70 instances, each acceleration alone (a timeout counts as 1800 s).

# %%
def no_accel(all_imsm=False):
    BASE = ["RM", "DP-iso", "VEQ", "CaLiG"]
    NO, ON = "No accel.", "IterOrbit"
    T_MS = 1800000.0

    def load(bench):
        d = iterorbit_sweep("no_accel", bench)[["dataset", "pattern_name", "status", "filter", "order", "enumerate", "count"]]
        print(f"no_accel ({bench}): {len(d)} logs, {d.status.value_counts().to_dict()}")
        return d

    def report(bench, loader, datasets):
        raw = load(bench)
        raw = raw[raw.dataset.isin(datasets)]
        old = raw[raw.status == "COMPLETED"].assign(tool=NO, read=0.0)
        old[["order", "enumerate"]] = old[["order", "enumerate"]].fillna(0)
        d = pd.concat([loader(), old], ignore_index=True)
        d = d[d.dataset.isin(datasets)]
        d["t"] = phase_time(d)
        t = d.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")
        print(f"\n================ {BENCH[bench]}: {len(raw)} instances run ================")
        print("completed: " + ", ".join(f"{k} {int(t[k].notna().sum())}" for k in [NO, ON] + BASE))
        best = t[BASE].min(axis=1)
        out = []
        for me in (NO, ON):
            for o in BASE + ["fastest of four"]:
                other = best if o == "fastest of four" else t[o]
                both = t[me].notna() & other.notna()
                r = other[both] / t[me][both]
                out.append({"config": me, "vs": o, "both": int(both.sum()), "speedup": round(gmean(r), 2),
                            "wins": int((r > 1).sum()), "losses": int((r < 1).sum())})
        print(pd.DataFrame(out).to_string(index=False))
        both = t[NO].notna() & t[ON].notna()
        r = t[NO][both] / t[ON][both]
        print(f"accelerations' speed-up (no accel. / default), {int(both.sum())} both completed: {gmean(r):.2f}x")
        per = []
        for ds in datasets:
            s, b = t.loc[ds], best.loc[ds]
            row = {"dataset": ds, "done_no": int(s[NO].notna().sum()), "done_on": int(s[ON].notna().sum()),
                   "done_RM": int(s.RM.notna().sum()), "done_best": int(s[BASE].notna().sum().max())}
            for me, tag in ((NO, "no"), (ON, "on")):
                k = s[me].notna() & b.notna()
                row[f"{tag}_vs_fastest"] = round(gmean(b[k] / s[me][k]), 2)
                k = s[me].notna() & s.RM.notna()
                row[f"{tag}_vs_RM"] = round(gmean(s.RM[k] / s[me][k]), 2)
            k = s[NO].notna() & s[ON].notna()
            row["accel_speedup"] = round(gmean(s[NO][k] / s[ON][k]), 2)
            per.append(row)
        print(pd.DataFrame(per).to_string(index=False))
        anyb = t[BASE].notna().any(axis=1)
        print(f"only no-accel completes: {int((t[NO].notna() & ~anyb).sum())}; "
              f"no-accel times out but a baseline completes: {int((t[NO].isna() & anyb).sum())}")
        return raw

    def counts_check(raw, ref, name):
        j = raw.dropna(subset=["count"]).merge(ref, on=["dataset", "pattern_name"])
        bad = j[j["count"] != j.ref]
        print(f"{name}: counts compared {len(j)}, mismatches {len(bad)}")
        if len(bad):
            print(bad.to_string(index=False))

    imsm = report("imsm", imsm_phases, IMSM + (["dblp", "eu2005", "youtube"] if all_imsm else []))
    mcc = report("mccreesh", mcc_phases, MCC)
    print()
    ox = csv("labelled_iterorbit_e2e.csv")
    emb = next(c for c in ox.columns if "emb" in c.lower() or "match" in c.lower())
    counts_check(imsm, ox[ox.status == "COMPLETED"][["dataset", "pattern_name", emb]].rename(columns={emb: "ref"}), "Labelled")
    e = csv("unlabelled_e2e_results.csv")
    e = e[(e.tool == "OrbitEX") & (e.status == "COMPLETED")]
    counts_check(mcc, e[["dataset", "pattern_name", "count"]].rename(columns={"count": "ref"}), "Unlabelled")

    abl = csv("accel_ablation_instances.csv")
    for r in iterorbit_sweep("accel_ablation").itertuples():
        row, cfg = divmod(int(r.task), 2)
        name = ("edge", "sym")[cfg]
        done = r.status == "COMPLETED"
        ms = sum(np.nan_to_num(v) for v in (r.filter_ms, r.order_us / 1e3, r.search_ms))
        abl.loc[row, name] = max(ms, TIME_FLOOR_MS) if done else T_MS
        if name == "sym" and isinstance(r.symmetry_multiplier, str):
            abl.loc[row, "multiplier"] = int(r.symmetry_multiplier)

    def eff(raw, loader, datasets):
        off = raw.assign(t=np.where(raw.status == "COMPLETED",
                                    phase_time(raw.assign(read=0.0).fillna({"order": 0, "enumerate": 0})), T_MS))
        on = loader()
        on = on[(on.tool == "IterOrbit") & on.dataset.isin(datasets)].assign(t=lambda x: phase_time(x))
        return off.set_index(["dataset", "pattern_name"]).t, on.set_index(["dataset", "pattern_name"]).t

    off_i, on_i = eff(imsm, imsm_phases, IMSM)
    off_m, on_m = eff(mcc, mcc_phases, MCC)
    off_all, on_all = pd.concat([off_i, off_m]), pd.concat([on_i, on_m])
    k = list(zip(abl.dataset, abl.pattern_name))
    abl["off"] = off_all.reindex(k).to_numpy()
    abl["on"] = on_all.reindex(k).fillna(T_MS).to_numpy()
    x = abl[abl.dataset.isin(IMSM + MCC) & abl.off.notna()].copy()
    print(f"\n================ single-feature ablation: {len(x)} instances "
          f"({x.groupby('bench').size().to_dict()}) ================")
    for c in ("off", "edge", "sym", "on"):
        x[f"{c}_done"] = x[c] < T_MS
    print("completed:", {c: int(x[f"{c}_done"].sum()) for c in ("off", "edge", "sym", "on")})
    for name, g in [("all", x)] + list(x.groupby("bench")):
        print(f"{name:9s} speed-up over both off: edge only {gmean(g.off / g.edge):.2f}x, "
              f"symmetry only {gmean(g.off / g.sym):.2f}x, both {gmean(g.off / g.on):.2f}x")
    sf = x.dataset == "scalefree"
    for label, g in (("scalefree", x[sf]), ("not scalefree", x[~sf])):
        print(f"{label} ({len(g)}): edge-only time / both off {gmean(g.edge / g.off):.3f} "
              f"({gmean(g.off / g.edge):.1f}x); symmetry-only / both off {gmean(g.sym / g.off):.3f}")
    print(f"symmetric patterns (multiplier > 1): {int((x.multiplier > 1).sum())} of {len(x)}, "
          f"median multiplier {x[x.multiplier > 1].multiplier.median():.0f}")
    gain = ~x.off_done & x.on_done
    print(f"both off times out, both on completes: {int(gain.sum())}; symmetry only completes "
          f"{int((gain & x.sym_done).sum())}, edge only {int((gain & x.edge_done).sum())}")
    return x


no_accel();

# %% [markdown]
# ## The original OrbitSI filter with IterOrbit's search
# One round of LDF + NLF, the induced candidate subgraph and the orbit filter (`--iterate 1
# --disable-semi-join`, g=4, no edge-orbit index filter), with IterOrbit's order, symmetry breaking and search;
# all 13 datasets, against IterOrbit and the four baselines.

# %%
def orbitsi_filter_e2e():
    ME = "OrbitSI"
    BASE = ["RM", "DP-iso", "VEQ", "CaLiG"]
    key = ["dataset", "pattern_name"]
    r = iterorbit_sweep("orbitsi_filter")
    d = pd.DataFrame({"bench": r.bench, "dataset": r.dataset, "pattern_name": r.pattern_name, "pattern_path": r.pattern_path,
                      "status": r.status, "filter": r["filter"], "order": np.nan_to_num(r["order"]),
                      "enumerate": np.nan_to_num(r["enumerate"]), "count": r["count"],
                      "refuted": (r.status == "COMPLETED") & r.search_ms.isna()})
    print(f"logs {len(d)}, duplicates {int(d.duplicated(key).sum())}")
    print("status:", d.groupby("bench").status.value_counts().unstack(fill_value=0).to_dict("index"))

    ox = d[d.status == "COMPLETED"].assign(tool=ME, read=0.0)
    for bench, loader, datasets in (("imsm", lambda: imsm_phases(IMSM_ALL), IMSM_ALL), ("mccreesh", mcc_phases, MCC)):
        x = pd.concat([loader(), ox[ox.bench == bench]], ignore_index=True)
        x = x[x.dataset.isin(datasets)]
        x["t"] = phase_time(x)
        t = x.pivot_table(index=key, columns="tool", values="t", aggfunc="first")
        tools = [ME, "IterOrbit"] + BASE
        for label, dsets in ([("all 8", IMSM_ALL), ("paper's 5", IMSM)] if bench == "imsm" else [("all 5", MCC)]):
            s = t[t.index.get_level_values(0).isin(dsets)]
            print(f"\n================ {BENCH[bench]}, {label} ({len(s)} instances) ================")
            print("completed: " + ", ".join(f"{k} {int(s[k].notna().sum())}" for k in tools))
            best = s[BASE].min(axis=1)
            out = []
            for o in ["IterOrbit"] + BASE + ["fastest of four"]:
                other = best if o == "fastest of four" else s[o]
                both = s[ME].notna() & other.notna()
                rr = other[both] / s[ME][both]
                out.append({"vs": o, "both": int(both.sum()), "speedup (other/OrbitSI)": round(gmean(rr), 2),
                            "OrbitSI wins": int((rr > 1).sum()), "losses": int((rr < 1).sum())})
            print(pd.DataFrame(out).to_string(index=False))
        per = []
        for ds in datasets:
            s = t.loc[ds]
            best = s[BASE].min(axis=1)
            row = {"dataset": ds, **{f"done_{k}": int(s[k].notna().sum()) for k in [ME, "IterOrbit", "RM"]},
                   "done_best_base": int(s[BASE].notna().sum().max())}
            for o, other in (("IterOrbit", s.IterOrbit), ("RM", s.RM), ("fastest4", best)):
                k = s[ME].notna() & other.notna()
                row[f"{o}/OrbitSI"] = round(gmean(other[k] / s[ME][k]), 2)
            per.append(row)
        print("\nper dataset (speed-up columns: other's time / OrbitSI's, >1 = OrbitSI faster)")
        print(pd.DataFrame(per).to_string(index=False))

    print("\n================ filtering ================")
    ref_i = csv("labelled_iterorbit_e2e.csv")
    ref_m = csv("unlabelled_e2e_phases.csv")
    ref_m = ref_m[ref_m.tool == "OrbitEX"]
    f = pd.concat([ref_i[key + ["filter_time_ms"]].rename(columns={"filter_time_ms": "io_filter"}),
                   ref_m[key + ["filter_ms"]].rename(columns={"filter_ms": "io_filter"})])
    f = d.merge(f, on=key).dropna(subset=["filter", "io_filter"])
    for b, g in f.groupby("bench"):
        print(f"{BENCH[b]}: OrbitSI filter geo {gmean(g['filter']):.2f} ms vs IterOrbit {gmean(g.io_filter):.2f} ms "
              f"({gmean(g.io_filter / g['filter'].clip(lower=1e-3)):.2f}x); refuted by OrbitSI's filter: {int(g.refuted.sum())}")
    emb = ref_i[ref_i.status == "COMPLETED"][key + ["embeddings"]].rename(columns={"embeddings": "ref"})
    e = csv("unlabelled_e2e_results.csv")
    e = e[(e.tool == "OrbitEX") & (e.status == "COMPLETED")][key + ["count"]].rename(columns={"count": "ref"})
    j = d.dropna(subset=["count"]).merge(pd.concat([emb, e]), on=key)
    bad = j[j["count"] != j.ref]
    print(f"\ncounts vs IterOrbit: {len(j)} compared, {len(bad)} mismatches")
    if len(bad):
        print(bad[key + ["count", "ref"]].to_string(index=False))


orbitsi_filter_e2e()

# %% [markdown]
# ## Unlabelled: the OrbitSI filter with and without symmetry breaking
# Against IterOrbit with and without its search accelerations and the four baselines. OrbitSI-nosym and
# IterOrbit-noacc differ only in the filter.

# %%
def orbitsi_filter_nosym():
    KEY = ["dataset", "pattern_name"]
    BASE = ["RM", "DP-iso", "VEQ", "CaLiG"]

    def load(sweep, tool):
        r = iterorbit_sweep(sweep)
        r = r[r.dataset.isin(MCC)]
        d = pd.DataFrame({"tool": tool, "dataset": r.dataset, "pattern_name": r.pattern_name, "status": r.status, "read": 0.0,
                          "filter": r["filter"], "order": np.nan_to_num(r["order"]), "enumerate": np.nan_to_num(r["enumerate"]),
                          "count": r["count"], "searched": r.search_ms.notna()})
        assert not d.duplicated(KEY).any(), sweep
        print(f"{tool}: {len(d)} logs, {d.status.value_counts().to_dict()}")
        return d

    runs = pd.concat([load("orbitsi_filter", "OrbitSI"), load("orbitsi_filter_nosym", "OrbitSI-nosym"),
                      load("no_accel", "IterOrbit-noacc")], ignore_index=True)
    done = runs[runs.status == "COMPLETED"]
    d = pd.concat([mcc_phases(), done], ignore_index=True)
    d["t"] = phase_time(d)
    t = d.pivot_table(index=KEY, columns="tool", values="t", aggfunc="first")
    CONF = ["OrbitSI", "OrbitSI-nosym", "IterOrbit", "IterOrbit-noacc"]
    print("\ncompleted (of 480): " + ", ".join(f"{k} {int(t[k].notna().sum())}" for k in CONF + BASE))
    print("decided by the filter alone (completed, no search): "
          + ", ".join(f"{k} {int((~g.searched & (g.status == 'COMPLETED')).sum())}" for k, g in runs.groupby("tool")))
    best = t[BASE].min(axis=1)
    rows = []
    for me in CONF:
        r = {"config": me}
        for o in CONF + BASE + ["fastest of four"]:
            if o == me:
                continue
            other = best if o == "fastest of four" else t[o]
            both = t[me].notna() & other.notna()
            r[o] = round(gmean(other[both] / t[me][both]), 2)
        rows.append(r)
    print("\nspeed-up of each configuration (row) over each system (column): column's time / row's, >1 = row faster")
    print(pd.DataFrame(rows).set_index("config").to_string())
    print("\nper family: speed-up of OrbitSI-nosym over ...")
    per = []
    for ds in MCC:
        s = t.loc[ds]
        row = {"family": ds, **{f"done {k}": int(s[k].notna().sum()) for k in CONF}}
        for o in ["OrbitSI", "IterOrbit-noacc", "IterOrbit", "RM"]:
            k = s["OrbitSI-nosym"].notna() & s[o].notna()
            row[o] = round(gmean(s[o][k] / s["OrbitSI-nosym"][k]), 2)
        k = s["OrbitSI-nosym"].notna() & best.loc[ds].notna()
        row["fastest4"] = round(gmean(best.loc[ds][k] / s["OrbitSI-nosym"][k]), 2)
        per.append(row)
    print(pd.DataFrame(per).to_string(index=False))
    lost = t["OrbitSI"].notna() & t["OrbitSI-nosym"].isna()
    print(f"\ncompleted with symmetry, not without: {int(lost.sum())} "
          f"{lost[lost].groupby(level=0).size().to_dict()}")
    e = csv("unlabelled_e2e_results.csv")
    ref = e[(e.tool == "OrbitEX") & (e.status == "COMPLETED")][KEY + ["count"]].rename(columns={"count": "ref"})
    for k, g in runs.groupby("tool"):
        j = g.dropna(subset=["count"]).merge(ref, on=KEY)
        print(f"counts {k} vs IterOrbit: {len(j)} compared, {int((j['count'] != j.ref).sum())} mismatches")


orbitsi_filter_nosym()

# %% [markdown]
# ## IterOrbit with a single filtering round
# IterOrbit at g=4 with one round (`--iterate 1`, otherwise the end-to-end build and options) against the paper's
# configuration and the baselines (labelled: HPRD, Human, Patents, WordNet, Yeast; unlabelled: all 480 pairs).

# %%
def single_pass():
    d = csv("single_pass_e2e.csv")
    print(f"finished tasks: {len(d)} {d.groupby(['bench', 'status']).size().to_dict()}")

    def gm(x):
        return float(np.exp(np.log(np.asarray(x, dtype=float)).mean())) if len(x) else np.nan

    BASE = ["RM", "DP-iso", "VEQ", "CaLiG"]
    for bench, ph in (("imsm", imsm_phases()), ("mccreesh", mcc_phases())):
        o = d[d.bench == bench].set_index(["dataset", "pattern_name"])
        if o.empty:
            continue
        nol = ph.assign(t=phase_time(ph)).pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t")
        if bench == "imsm":
            pp = csv("labelled_e2e_per_pattern.csv").set_index(["dataset", "pattern_name"])
            cols = {"IterOrbit": "ox", "RM": "RapidMatch", "DP-iso": "DP-iso", "VEQ": "VEQ", "CaLiG": "CaLiG"}
            done = pd.DataFrame({k: pp[f"{v}_done"].astype(bool) for k, v in cols.items()})
            emb = pp.ox_emb
        else:
            e = csv("unlabelled_e2e_results.csv").replace({"tool": {"OrbitEX": "IterOrbit", "RapidMatch": "RM"}})
            done = e.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="done", aggfunc="first").astype(bool)
            emb = e[e.tool == "IterOrbit"].set_index(["dataset", "pattern_name"])["count"]
        idx = o.index
        done = done.reindex(idx)
        t = nol.reindex(idx)[TOOLS].where(done[TOOLS])
        t["OrbitSI"] = o.nol_ms.where(o.done)
        print(f"\n===== {BENCH[bench]}: {len(idx)} instances finished so far "
              f"({o.groupby(level=0).size().to_dict()}) =====")
        comp = {k: int(t[k].notna().sum()) for k in ["OrbitSI"] + TOOLS}
        comp["fastest of four"] = int(t[BASE].notna().any(axis=1).sum())
        print("completed:", comp)
        both = o.done & done.IterOrbit
        mism = (o["count"][both] != emb.reindex(idx)[both]).sum()
        print(f"count mismatches vs IterOrbit (both completed): {int(mism)} of {int(both.sum())}")
        print("speed-up of the single pass over X (X's time / the single pass's, excluding loading; >1 = single pass faster):")
        for x in TOOLS + ["fastest of four"]:
            other = t[BASE].min(axis=1) if x == "fastest of four" else t[x]
            ok = t.OrbitSI.notna() & other.notna()
            r = other[ok] / t.OrbitSI[ok]
            print(f"  vs {x:15s} n={int(ok.sum()):4d}  {gm(r):6.2f}x  single pass wins {(r > 1).sum():4d}, losses {(r < 1).sum():4d}")
        ok = t.OrbitSI.notna() & t.IterOrbit.notna()
        r = (t.IterOrbit[ok] / t.OrbitSI[ok]).rename("r").reset_index()
        print("per dataset, IterOrbit's time / the single pass's (<1 = iterating pays off):")
        print(r.groupby("dataset").r.agg(n="size", ratio=gm, single_faster=lambda x: (x > 1).sum(),
                                         iterorbit_faster=lambda x: (x < 1).sum()).round(2).to_string())
        print("only the single pass completes:", int((t.OrbitSI.notna() & t.IterOrbit.isna()).sum()),
              "; only IterOrbit completes:", int((t.OrbitSI.isna() & t.IterOrbit.notna()).sum()))


single_pass()

# %% [markdown]
# ## Oracle: the faster of the paper's configuration and a single pass, per instance
# How much a per-instance choice of filter effort could gain, from runs that already exist (labelled: the five
# datasets the single pass ran on; unlabelled: all pairs). OSI is the OrbitSI package, for reference.

# %%
def oracle():
    FIVE = ["hprd", "human", "patents", "wordnet", "yeast"]

    def io_times(bench):
        d = imsm_phases(PLOT_IMSM) if bench == "imsm" else mcc_phases()
        d = d[d.tool == "IterOrbit"]
        return pd.Series(phase_time(d).to_numpy(), index=pd.MultiIndex.from_arrays([d.dataset, d.pattern_name]))

    def pass_times(bench):
        o = csv("single_pass_e2e.csv")
        o = o[(o.bench == bench) & (o.status == "COMPLETED")]
        t = (o.filter_ms + o.order_ms + o.search_ms).clip(lower=TIME_FLOOR_MS)
        return pd.Series(t.to_numpy(), index=pd.MultiIndex.from_arrays([o.dataset, o.pattern_name])), \
            pd.Series(o["count"].to_numpy(), index=pd.MultiIndex.from_arrays([o.dataset, o.pattern_name]))

    def osi_times(bench):
        d = with_index(orbitsi_phases(bench), bench)
        return pd.Series(phase_time(d).to_numpy(), index=pd.MultiIndex.from_arrays([d.dataset, d.pattern_name]))

    def report(bench, datasets, universe):
        io, ps, osi = io_times(bench), pass_times(bench)[0], osi_times(bench)
        keep = lambda s: s[s.index.get_level_values(0).isin(datasets)]
        io, ps, osi = keep(io), keep(ps), keep(osi)
        t = pd.DataFrame({"IO": io, "pass": ps, "OSI": osi})
        print(f"\n=== {BENCH[bench]} ({', '.join(LABEL[d] for d in datasets)}; {universe} instances) ===")
        print(f"completed: IO {t.IO.notna().sum()}, single pass {t['pass'].notna().sum()}, "
              f"either (oracle) {(t.IO.notna() | t['pass'].notna()).sum()}, OSI {t.OSI.notna().sum()}")
        both = t.dropna(subset=["IO", "pass"])
        faster = (both["pass"] < both.IO).sum()
        print(f"both IO and single pass complete: {len(both)}; single pass faster on {faster} ({100 * faster / len(both):.0f}%), "
              f"IO faster on {(both.IO < both['pass']).sum()}")
        orc = both[["IO", "pass"]].min(axis=1)
        print(f"  geo time: IO {gmean(both.IO):.3g} ms, single pass {gmean(both['pass']):.3g} ms, oracle {gmean(orc):.3g} ms "
              f"-> oracle {gmean(both.IO) / gmean(orc):.2f}x faster than always-IO, {gmean(both['pass']) / gmean(orc):.2f}x than always-pass")
        save_ = (both.IO - orc)
        print(f"  total time: IO {both.IO.sum() / 1e3:.1f} s, oracle {orc.sum() / 1e3:.1f} s "
              f"({100 * save_.sum() / both.IO.sum():.1f}% saved); where single pass wins, it saves a median "
              f"{np.median(save_[save_ > 0]):.2g} ms (IO median there {np.median(both.IO[save_ > 0]):.2g} ms)")
        big = both[(both.IO >= 1000) | (both["pass"] >= 1000)]
        if len(big):
            bo = big[["IO", "pass"]].min(axis=1)
            print(f"  instances taking >= 1 s in either: {len(big)}; single pass faster on {(big['pass'] < big.IO).sum()}; "
                  f"oracle {gmean(big.IO) / gmean(bo):.2f}x over IO there")
        only_io = t.IO.notna() & t["pass"].isna()
        only_ps = t["pass"].notna() & t.IO.isna()
        print(f"  completed only by IO: {only_io.sum()}, only by the single pass: {only_ps.sum()}")
        print("  per dataset (n both, single pass faster, oracle speed-up over IO):")
        for ds in datasets:
            b = both[both.index.get_level_values(0) == ds]
            if len(b):
                o = b[["IO", "pass"]].min(axis=1)
                print(f"    {LABEL[ds]:10s} n={len(b):3d}  pass faster {(b['pass'] < b.IO).sum():3d}  "
                      f"oracle {gmean(b.IO) / gmean(o):.2f}x  (IO/pass {gmean(b.IO) / gmean(b['pass']):.2f})")
        tri = t.dropna()
        if len(tri):
            o3 = tri[["IO", "pass"]].min(axis=1)
            print(f"  vs OSI on the {len(tri)} instances all three complete: IO/OSI {gmean(tri.IO) / gmean(tri.OSI):.2f}, "
                  f"oracle/OSI {gmean(o3) / gmean(tri.OSI):.2f} (below 1: faster than OSI); OSI fastest of the three on "
                  f"{((tri.OSI < tri.IO) & (tri.OSI < tri['pass'])).sum()}")
        ioc = None
        if bench == "imsm":
            ox = csv("labelled_iterorbit_e2e.csv")
            ox = ox[ox.status == "COMPLETED"]
            ioc = pd.Series(ox.embeddings.to_numpy(), index=pd.MultiIndex.from_arrays([ox.dataset, ox.pattern_name]))
        else:
            e = csv("unlabelled_e2e_phases.csv")
            e = e[(e.status == "COMPLETED") & (e.tool == "OrbitEX")]
            if "count" in e:
                ioc = pd.Series(e["count"].to_numpy(), index=pd.MultiIndex.from_arrays([e.dataset, e.pattern_name]))
        if ioc is not None:
            pc = pass_times(bench)[1]
            j = pd.DataFrame({"io": ioc, "ps": pc}).dropna()
            print(f"  counts: {len(j)} compared, {(j.io != j.ps).sum()} differ")

    report("imsm", FIVE, 500)
    report("mccreesh", MCC, 480)


oracle()

# %% [markdown]
# ## The IterOrbit re-sweep of 2026-10-17
# A later build of IterOrbit (marked IterOrbit* elsewhere) against the paper's sweep and the other systems;
# "+sym" adds the binary's separate symmetry timer.

# %%
def resweep():
    ALL8 = IMSM_ALL
    SOTA = ["RM", "DP-iso", "VEQ", "CaLiG"]

    def sp(a, b):
        k = a.notna() & b.notna()
        r = a[k].clip(lower=TIME_FLOOR_MS) / b[k].clip(lower=TIME_FLOOR_MS)
        return f"{gmean(r):.2f}x ({(r > 1).sum()}/{(r < 1).sum()})" if k.any() else "-"

    def report(name, new, old_t, old_done, old_count, others):
        done = new.status == "COMPLETED"
        od = old_done.reindex(new.index).fillna(False).astype(bool)
        print(f"\n== {name}: new completed {int(done.sum())}/{len(new)} (old {int(od.sum())}); "
              f"newly completed {int((done & ~od).sum())}, newly timed out {int((~done & od).sum())}; "
              + ", ".join(f"{c} {int(others[c].notna().sum())}" for c in others.columns))
        oc = old_count.reindex(new.index)
        k = done & od & oc.notna()
        sat = k & ((oc >= 2**63 - 1) | (new["count"] >= 2**63 - 1))
        print(f"  counts vs old sweep: {int(k.sum())} compared, mismatches {int(((new['count'] != oc) & k & ~sat).sum())}, saturated {int(sat.sum())}")
        nt = new.t.where(done)
        print(f"  new vs old IterOrbit (paper convention): {sp(old_t.reindex(new.index).where(od), nt)}")
        for conv, col in (("paper", "t"), ("+sym", "t_sym")):
            x = new[col].where(done)
            fastest = others[SOTA].min(axis=1)
            print(f"  [{conv}] " + "; ".join(f"{c} {sp(others[c], x)}" for c in others.columns) + f"; fastest-SOTA {sp(fastest, x)}")

    new = iterorbit_resweep("imsm")
    ox = csv("labelled_iterorbit_e2e.csv").set_index(["dataset", "pattern_name"])
    ox_done = ox.status == "COMPLETED"
    ox_t = (ox.filter_time_ms + ox.order_time_us.fillna(0) / 1e3 + ox.search_time_ms.fillna(0)).where(ox_done)
    s = csv("labelled_sota_runs.csv")
    s = s[s.status == "COMPLETED"]
    s["t"] = (s.filter_vertices_time_s + s.generate_query_plan_time_s.fillna(0) + s.build_table_time_s.fillna(0)
              + s.enumerate_time_s.fillna(0)) * 1e3
    oth = s.assign(tool=s.tool.map({"rapidmatch": "RM", "dpiso": "DP-iso", "veq": "VEQ", "calig": "CaLiG"})) \
        .pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first").reindex(new.index)[SOTA]
    oth["OrbitSI"] = orbitsi_e2e_totals("imsm").t.reindex(new.index)
    for name, ds in (("LABELLED, 5 graphs", IMSM), ("LABELLED, 8 graphs", ALL8)):
        m = new.index.get_level_values(0).isin(ds)
        report(name, new[m], ox_t, ox_done, ox.embeddings, oth[m])
    print("  per dataset, new vs old IterOrbit / vs RM (paper convention):")
    for d in ALL8:
        m = new.index.get_level_values(0) == d
        nt = new.t.where(new.status == "COMPLETED")[m]
        print(f"    {d:8s} {sp(ox_t.reindex(nt.index), nt):18s} {sp(oth.RM[m], nt)}")

    new = iterorbit_resweep("mccreesh")
    e = csv("unlabelled_e2e_results.csv").replace({"tool": {"OrbitEX": "IterOrbit", "RapidMatch": "RM"}}) \
        .set_index(["dataset", "pattern_name"])
    old = e[e.tool == "IterOrbit"]
    ph = mcc_phases()
    ph["t"] = ph["filter"] + ph["order"] + ph["enumerate"]
    pv = ph.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")
    oth = pv.reindex(new.index)[SOTA]
    oth["OrbitSI"] = orbitsi_e2e_totals("mccreesh").t.reindex(new.index)
    g = glasgow_runs()
    oth["Glasgow (1 ms floor)"] = g.runtime[g.decided].reindex(new.index).clip(lower=1.0)
    report("UNLABELLED, 480", new, pv.IterOrbit, old.done.astype(bool), e[e.tool == "Glasgow"].ref_count.fillna(old["count"]), oth)
    print("  per family, new vs old IterOrbit / vs RM / vs OrbitSI (paper convention):")
    for d, grp in new.groupby(level=0):
        nt = grp.t.where(grp.status == "COMPLETED")
        print(f"    {d:14s} {sp(pv.IterOrbit.reindex(grp.index), nt):18s} {sp(oth.RM.loc[grp.index], nt):18s} {sp(oth.OrbitSI.loc[grp.index], nt)}")


resweep()

# %% [markdown]
# ## The OrbitSI package end to end
# The package's published configuration (C++ core) on both benchmarks: completions, counts against the other
# systems' and Glasgow's, and speed-ups without loading (its stage timers) and wall clock (its load + stages).

# %%
def orbitsi_package_e2e():
    ALL8 = IMSM_ALL
    BASE = ["IterOrbit", "RM", "DP-iso", "VEQ", "CaLiG"]

    def gm(x):
        return gmean(np.asarray(x, dtype=float))

    def compare(o, t_no, t_wall, label):
        print(f"  {label}: OrbitSI completed {int(o.t.notna().sum())}/{len(o)}; "
              + ", ".join(f"{b} {int(t_no[b].notna().sum())}" for b in t_no.columns))
        for name, t, oc in (("without loading", t_no, o.t), ("wall clock", t_wall, o.wall.where(o.t.notna()))):
            out = []
            for b in list(t.columns) + ["fastest-SOTA"]:
                other = t[["RM", "DP-iso", "VEQ", "CaLiG"]].min(axis=1) if b == "fastest-SOTA" else t[b]
                k = oc.notna() & other.notna()
                r = other[k].clip(lower=TIME_FLOOR_MS) / oc[k].clip(lower=TIME_FLOOR_MS)
                out.append(f"{b} {gm(r):.2f}x ({(r > 1).sum()}/{(r < 1).sum()})")
            print(f"    {name:16s} " + "; ".join(out))

    o = orbitsi_e2e_totals("imsm")
    pp = csv("labelled_e2e_per_pattern.csv").set_index(["dataset", "pattern_name"])
    col = {"IterOrbit": "ox", "RM": "RapidMatch", "DP-iso": "DP-iso", "VEQ": "VEQ", "CaLiG": "CaLiG"}
    done = pd.DataFrame({k: pp[f"{v}_done"].astype(bool) for k, v in col.items()})
    emb = pd.DataFrame({k: pp[f"{v}_emb"] for k, v in col.items()}).where(done)
    wall = pd.DataFrame({k: pp[f"{v}_ms"] for k, v in col.items()}).where(done)
    ox = csv("labelled_iterorbit_e2e.csv")
    ox = ox[ox.status == "COMPLETED"].set_index(["dataset", "pattern_name"])
    s = csv("labelled_sota_runs.csv")
    s = s[s.status == "COMPLETED"]
    s["t"] = (s.filter_vertices_time_s + s.generate_query_plan_time_s.fillna(0) + s.build_table_time_s.fillna(0)
              + s.enumerate_time_s.fillna(0)) * 1e3
    nol = s.assign(tool=s.tool.map({"rapidmatch": "RM", "dpiso": "DP-iso", "veq": "VEQ", "calig": "CaLiG"})) \
        .pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")
    nol["IterOrbit"] = ox.filter_time_ms + ox.order_time_us.fillna(0) / 1e3 + ox.search_time_ms.fillna(0)
    nol = nol.reindex(o.index)[BASE].where(done.reindex(o.index))
    wall = wall.reindex(o.index)
    ref = emb.IterOrbit.copy()
    for b in ["RM", "DP-iso", "CaLiG", "VEQ"]:
        e = emb[b].where(emb[b] < 2**32) if b == "VEQ" else emb[b]
        ref = ref.fillna(e)
    ref = ref.reindex(o.index)
    oc = o["count"].astype(float)
    k = oc.notna() & ref.notna()
    sat = k & (ref >= 2**63 - 1)
    print("== LABELLED")
    print(f"  counts: {int(k.sum())} checkable of {int(oc.notna().sum())} completed; mismatches "
          f"{int(((oc != ref) & k & ~sat).sum())}; IterOrbit saturated at 2^63-1 (OrbitSI's u128 count shown): "
          + ", ".join(f"{i[0]}/{i[1]} {int(o['count'][i])}" for i in oc[sat].index)
          + f"; completed only by OrbitSI (no reference): {int((oc.notna() & ref.isna()).sum())}")
    for name, ds in (("5 graphs", IMSM), ("8 graphs", ALL8)):
        m = o.index.get_level_values(0).isin(ds)
        compare(o[m], nol[m], wall[m], name)
    print("  per dataset (without loading, vs RM / DP-iso / IterOrbit; completed OrbitSI/RM/DP-iso/IterOrbit):")
    for d in ALL8:
        m = o.index.get_level_values(0) == d
        oo, tt = o[m], nol[m]
        sp = []
        for b in ["RM", "DP-iso", "IterOrbit"]:
            kk = oo.t.notna() & tt[b].notna()
            sp.append(f"{b} {gm(tt[b][kk].clip(lower=TIME_FLOOR_MS) / oo.t[kk].clip(lower=TIME_FLOOR_MS)):.2f}x")
        print(f"    {d:8s} {int(oo.t.notna().sum()):3d}/{int(tt.RM.notna().sum()):3d}/{int(tt['DP-iso'].notna().sum()):3d}/"
              f"{int(tt.IterOrbit.notna().sum()):3d}  " + ", ".join(sp) + f"  load {oo.load.median():.0f} ms (median)")

    o = orbitsi_e2e_totals("mccreesh")
    e = csv("unlabelled_e2e_results.csv").replace({"tool": {"OrbitEX": "IterOrbit", "RapidMatch": "RM"}})
    e = e.set_index(["dataset", "pattern_name"])
    ref = e[e.tool == "Glasgow"].ref_count.reindex(o.index)
    ph = mcc_phases()
    nol = ph.assign(t=ph["filter"] + ph["order"] + ph["enumerate"]) \
        .pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first").reindex(o.index)[BASE]
    wall = e[e.tool.isin(BASE) & e.done.astype(bool)].reset_index() \
        .pivot_table(index=["dataset", "pattern_name"], columns="tool", values="wall_ms", aggfunc="first").reindex(o.index)[BASE]
    oc = o["count"].astype(float)
    k = oc.notna() & ref.notna()
    print("\n== UNLABELLED")
    print(f"  counts: {int(k.sum())} checkable against the reference (ref_count) of {int(oc.notna().sum())} completed; "
          f"mismatches {int(((oc != ref) & k & (ref < 2**63 - 1)).sum())}; reference saturated: {int((k & (ref >= 2**63 - 1)).sum())}")
    compare(o, nol, wall, "all 480")
    for d, grp in o.groupby(level=0):
        tt = nol.loc[grp.index]
        sp = []
        for b in ["RM", "DP-iso", "IterOrbit"]:
            kk = grp.t.notna() & tt[b].notna()
            sp.append(f"{b} {gm(tt[b][kk].clip(lower=TIME_FLOOR_MS) / grp.t[kk].clip(lower=TIME_FLOOR_MS)):.2f}x")
        print(f"    {d:14s} {int(grp.t.notna().sum()):3d}/{int(tt.RM.notna().sum()):3d}/{int(tt['DP-iso'].notna().sum()):3d}/"
              f"{int(tt.IterOrbit.notna().sum()):3d}  " + ", ".join(sp))
    g = glasgow_runs()
    g = g.runtime[g.decided].reindex(o.index)
    kk = g.notna() & o.t.notna()
    r = g[kk].clip(lower=1.0) / o.t[kk].clip(lower=1.0)
    print(f"  vs Glasgow (completed {int(g.notna().sum())}; both floored at 1 ms, without loading): {gm(r):.2f}x "
          f"({(r > 1).sum()}/{(r < 1).sum()})")


orbitsi_package_e2e()

# %% [markdown]
# ## PAR-5
# The mean time over every instance, a run that did not complete within the 1800 s limit counting as 5 x 1800 s;
# without loading and wall clock. IterOrbit* is the 2026-10-17 re-sweep.

# %%
def par5():
    LIMIT_S = 1800.0
    PENALTY_S = 5 * LIMIT_S

    def score(t_ms, n):
        t = np.asarray(t_ms, dtype=float) / 1e3
        return (t.sum() + (n - len(t)) * PENALTY_S) / n if n else np.nan

    def table(times, keys, groups, title):
        rows = []
        for name, ds in groups:
            k = [x for x in keys if x[0] in ds]
            r = {"set": name, "n": len(k)}
            for tool, s in times.items():
                s = s.reindex(pd.MultiIndex.from_tuples(k)).dropna()
                r[tool] = score(s.clip(lower=TIME_FLOOR_MS), len(k))
                r[f"{tool} done"] = len(s)
            rows.append(r)
        d = pd.DataFrame(rows).set_index("set")
        tools = list(times)
        print(f"\n{title} (PAR-5, s; completed in brackets)")
        print(d.apply(lambda r: pd.Series({"n": int(r.n), **{t: f"{r[t]:8.1f} ({int(r[f'{t} done'])})" for t in tools}}),
                      axis=1).to_string())
        return d

    def series(d, col):
        return d.set_index(["dataset", "pattern_name"])[col]

    keys = instance_keys("imsm")
    ph = imsm_phases(IMSM_ALL)
    ph["t"] = ph["filter"] + ph["order"] + ph["enumerate"]
    osi = orbitsi_phases("imsm")
    osi["t"] = osi["filter"] + osi["order"] + osi["enumerate"]
    new = iterorbit_resweep("imsm")
    nol = {t: series(ph[ph.tool == t], "t") for t in TOOLS}
    nol["IterOrbit*"] = new.t.where(new.status == "COMPLETED").dropna()
    nol["OrbitSI"] = series(osi, "t")
    pp = csv("labelled_e2e_per_pattern.csv").set_index(["dataset", "pattern_name"])
    col = {"IterOrbit": "ox", "RM": "RapidMatch", "DP-iso": "DP-iso", "VEQ": "VEQ", "CaLiG": "CaLiG"}
    wall = {k: pp[f"{v}_ms"].where(pp[f"{v}_done"].astype(bool)).dropna() for k, v in col.items()}
    wall["IterOrbit*"] = new.elapsed_ms.where(new.status == "COMPLETED").dropna()
    wall["OrbitSI"] = series(osi.assign(w=osi.read + osi.t), "w")
    order = ["IterOrbit", "IterOrbit*", "OrbitSI", "RM", "DP-iso", "VEQ", "CaLiG"]
    groups = [(LABEL[d], [d]) for d in IMSM_ALL] + [("all 8", IMSM_ALL), ("paper's 5", IMSM)]
    table({t: nol[t] for t in order}, keys, groups, "LABELLED, without loading")
    table({t: wall[t] for t in order}, keys, groups, "LABELLED, wall clock")

    e = csv("unlabelled_e2e_results.csv").replace({"tool": {"OrbitEX": "IterOrbit", "RapidMatch": "RM"}})
    keys = [tuple(x) for x in e[e.tool == "IterOrbit"][["dataset", "pattern_name"]].values]
    ph = mcc_phases()
    ph["t"] = ph["filter"] + ph["order"] + ph["enumerate"]
    osi = orbitsi_phases("mccreesh")
    osi["t"] = osi["filter"] + osi["order"] + osi["enumerate"]
    new = iterorbit_resweep("mccreesh")
    nol = {t: series(ph[ph.tool == t], "t") for t in TOOLS}
    nol["IterOrbit*"] = new.t.where(new.status == "COMPLETED").dropna()
    nol["OrbitSI"] = series(osi, "t")
    g = glasgow_runs()
    g = g[g.decided]
    nol["Glasgow"] = g.runtime
    done = e[e.done.astype(bool)]
    wall = {t: series(done[done.tool == t], "wall_ms") for t in TOOLS}
    wall["IterOrbit*"] = new.elapsed_ms.where(new.status == "COMPLETED").dropna()
    wall["OrbitSI"] = series(osi.assign(w=osi.read + osi.t), "w")
    wall["Glasgow"] = g.wall
    order += ["Glasgow"]
    groups = [(LABEL[d], [d]) for d in MCC] + [("all 480", MCC)]
    table({t: nol[t] for t in order}, keys, groups, "UNLABELLED, without loading")
    table({t: wall[t] for t in order}, keys, groups, "UNLABELLED, wall clock")


par5()

# %% [markdown]
# ## The earlier five-graph end-to-end summary
# The end-to-end numbers as the drafts of 2026-10-04 to 10-17 reported them: labelled restricted to HPRD, Human,
# Patents, WordNet and Yeast, no index charge, IterOrbit against the four baselines, wall clock and without loading.

# %%
def five_graph_e2e():
    BASE = ["RM", "DP-iso", "VEQ", "CaLiG"]

    def head_to_head(t, done, label):
        print(f"\n--- {label} ---")
        rows = []
        for o in BASE + ["fastest of four"]:
            other = t[BASE].min(axis=1) if o == "fastest of four" else t[o]
            both = t["IterOrbit"].notna() & other.notna()
            r = other[both] / t["IterOrbit"][both]
            rows.append({"vs": o, "both": int(both.sum()), "speedup": round(gmean(r), 2),
                         "wins": int((r > 1).sum()), "losses": int((r < 1).sum())})
        print(pd.DataFrame(rows).to_string(index=False))
        best = t[BASE].min(axis=1)
        both = t["IterOrbit"].notna() & best.notna()
        r = (best[both] / t["IterOrbit"][both]).rename("r").reset_index()
        print("per dataset vs fastest of four:")
        print(r.groupby("dataset").r.agg(n="size", speedup=gmean, wins=lambda x: int((x > 1).sum()),
                                         losses=lambda x: int((x < 1).sum())).round(2).to_string())
        for o in BASE:
            b = t["IterOrbit"].notna() & t[o].notna()
            rr = (t[o][b] / t["IterOrbit"][b]).rename("r").reset_index()
            print(f"  per dataset vs {o}: " + ", ".join(f"{d} {gmean(g.r):.2f}" for d, g in rr.groupby("dataset")))
        return r

    def report(name, done, wall, ph):
        print(f"\n================ {name} ================")
        n = len(done)
        comp = done.sum().astype(int)
        anyb = done[BASE].any(axis=1)
        print(f"instances {n}; completed: " + ", ".join(f"{k} {comp[k]}" for k in TOOLS)
              + f", fastest of four {int(anyb.sum())}")
        print("completed per dataset:")
        print(done.groupby(level="dataset").sum().astype(int).to_string())
        only_ox = done.IterOrbit & ~anyb
        lost = ~done.IterOrbit & anyb
        print(f"only IterOrbit completes: {int(only_ox.sum())}; IterOrbit times out but a baseline completes: "
              f"{int(lost.sum())} {lost[lost].groupby(level='dataset').size().to_dict()}")
        for d in done.index.get_level_values(0).unique():
            sub = done.loc[d]
            print(f"  {d}: IterOrbit {int(sub.IterOrbit.sum())}, best baseline {int(sub[BASE].sum().max())}, "
                  f"only IterOrbit {int((sub.IterOrbit & ~sub[BASE].any(axis=1)).sum())}")
        head_to_head(wall.where(done), done, "wall clock (loading included)")
        nol = ph.assign(t=(ph["filter"] + ph["order"] + ph["enumerate"]).clip(lower=TIME_FLOOR_MS)) \
            .pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t")
        nol = nol.reindex(done.index)[TOOLS].where(done)
        r = head_to_head(nol, done, "without loading (filter + order + enumerate)")
        ox = ph[ph.tool == "IterOrbit"]
        tot = ox.read + ox["filter"] + ox.order + ox.enumerate
        print(f"\nIterOrbit: median loading share {100 * (ox.read / tot).median():.0f}%; "
              f"filter share of non-loading time, median {100 * (ox['filter'] / (tot - ox.read).clip(lower=1e-3)).median():.0f}%")
        lossy = r[r.r < 1]
        print(f"losses to fastest of four (without loading): {len(lossy)}; of those under 10 ms (IterOrbit): "
              f"{int((nol.IterOrbit.reindex(pd.MultiIndex.from_frame(lossy[['dataset', 'pattern_name']])) < 10).sum())}")

    pp = csv("labelled_e2e_per_pattern.csv")
    pp = pp[pp.dataset.isin(IMSM)].set_index(["dataset", "pattern_name"])
    col = {"IterOrbit": "ox", "RM": "RapidMatch", "DP-iso": "DP-iso", "VEQ": "VEQ", "CaLiG": "CaLiG"}
    done = pd.DataFrame({k: pp[f"{v}_done"].astype(bool) for k, v in col.items()})
    wall = pd.DataFrame({k: pp[f"{v}_ms"] for k, v in col.items()})
    report("LABELLED (5 datasets)", done, wall, imsm_phases())
    emb = pd.DataFrame({k: pp[f"{v}_emb"] for k, v in col.items()})
    print(f"\nVEQ counts above 2^32 differing from IterOrbit's: "
          f"{int(((emb.VEQ != emb.IterOrbit) & done.VEQ & done.IterOrbit & (emb.IterOrbit >= 2**32)).sum())}; "
          f"other count mismatches (both completed, below 2^32): "
          f"{int(sum(((emb[k] != emb.IterOrbit) & done[k] & done.IterOrbit & (emb.IterOrbit < 2**32)).sum() for k in BASE))}; "
          f"IterOrbit saturated (2^63-1): {int((emb.IterOrbit >= 2**63 - 1).sum())}")
    e = csv("unlabelled_e2e_results.csv").replace({"tool": {"OrbitEX": "IterOrbit", "RapidMatch": "RM"}})
    e = e[e.tool.isin(TOOLS)]
    done = e.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="done", aggfunc="first").astype(bool)[TOOLS]
    wall = e.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="wall_ms", aggfunc="first")[TOOLS]
    report("UNLABELLED", done, wall, mcc_phases())


five_graph_e2e()

# %% [markdown]
# ## Charging the index builds
# What charging each system's loader-built indexes does to IterOrbit's speed-ups: none, the indexes the filter
# reads ("used", as in the paper), or every index the loader builds ("all"); with the build times themselves.

# %%
def index_charge():
    ib = csv("index_build_logs.csv")
    ix = ib[ib.kind == "INDEX"].groupby("path")[["load", "reverse", "edge", "nlf"]].median()
    ix["used"] = ix.reverse + ix.nlf
    ix["all"] = ix.used + ix.edge
    BASE = ["OrbitSI", "RM", "DP-iso", "VEQ", "CaLiG"]
    osi_ix = ib[ib.kind == "LOADSPLIT"].dropna(axis=1, how="all").groupby("path").median(numeric_only=True)
    osi_ix["used"] = osi_ix.label_index + osi_ix.nlf
    osi_ix["all"] = osi_ix.used + osi_ix.nlf_table

    def compare(name, ph, charge, datasets, osi_charge):
        ph = ph.assign(t=ph["filter"] + ph["order"] + ph["enumerate"])
        t = ph.pivot_table(index=["dataset", "pattern_name"], columns="tool", values="t", aggfunc="first")
        print(f"\n== {name}")
        for label, ds in [(LABEL.get(d, d), [d]) for d in datasets] + [("pooled", datasets)]:
            x = t[t.index.get_level_values(0).isin(ds)]
            out = []
            for b in BASE:
                k = x.IterOrbit.notna() & x[b].notna()
                cells = []
                for extra in ("none", "used", "all"):
                    src = osi_charge if b == "OrbitSI" else charge
                    add = 0 if extra == "none" else src[extra].reindex(x.index[k]).fillna(0)
                    r = (x[b][k] + add).clip(lower=TIME_FLOOR_MS) / x.IterOrbit[k].clip(lower=TIME_FLOOR_MS)
                    cells.append(f"{gmean(r):.2f}x {(r > 1).sum()}/{(r < 1).sum()}")
                out.append(f"{b}: " + " -> ".join(cells))
            print(f"  {label:12s} " + " | ".join(out))

    lab = ix[ix.index.str.contains("dataset-imsm")].copy()
    lab.index = [p.split("/")[-1].replace(".graph", "") for p in lab.index]
    print("Labelled index build times (ms, median of 3):")
    print(lab.round(1).to_string())
    ol = osi_ix[osi_ix.index.str.contains("dataset-imsm")].copy()
    ol.index = [p.split("/")[-1].replace(".graph", "") for p in ol.index]
    print("\nOrbitSI loading steps (ms, median of 3):")
    print(ol[["parse", "graph_index", "nlf_table", "copy", "label_index", "nlf", "used", "all"]].round(1).to_string())
    print(f"one-time warm-up: {osi_ix.warm_up.max():.1f} ms")
    ph = pd.concat([imsm_phases(IMSM_ALL), orbitsi_phases("imsm")], ignore_index=True)
    keys = ph[["dataset", "pattern_name"]].drop_duplicates().set_index(["dataset", "pattern_name"])
    charge = pd.DataFrame({c: [lab.loc[d, c] for d, _ in keys.index] for c in ("used", "all")}, index=keys.index)
    osi_charge = pd.DataFrame({c: [ol.loc[d, c] for d, _ in keys.index] for c in ("used", "all")}, index=keys.index)
    compare("LABELLED (speed-up of IterOrbit: no index charge -> used indexes -> all loader indexes; wins/losses)",
            ph, charge, IMSM_ALL, osi_charge)

    pairs = csv("unlabelled_pairs.csv").set_index(["dataset", "instance"])
    unl = ix[ix.index.str.contains("dataset-mccreesh")]
    tgt = pairs.target_path
    charge = pd.DataFrame({c: tgt.map(unl[c]) for c in ("used", "all")})
    charge.index.names = ["dataset", "pattern_name"]
    print("\nUnlabelled index build times per family (ms, median over target graphs):")
    print(charge.groupby(level=0).median().round(2).to_string())
    ou = osi_ix[osi_ix.index.str.contains("dataset-mccreesh")]
    osi_charge = pd.DataFrame({c: tgt.map(ou[c]) for c in ("used", "all")})
    osi_charge.index.names = ["dataset", "pattern_name"]
    print("\nOrbitSI index build times per family (ms, median over target graphs):")
    print(osi_charge.groupby(level=0).median().round(2).to_string())
    ph = pd.concat([mcc_phases(), orbitsi_phases("mccreesh")], ignore_index=True)
    compare("UNLABELLED", ph, charge, MCC, osi_charge)


index_charge()

# %% [markdown]
# ## Convergence of the fixpoint per dataset (the former Table `tab:convergence`)
# Labelled: reduction after rounds 1-3 and at the fixpoint as 1 - geometric mean of the surviving fractions
# (five datasets, all graphlet sizes); unlabelled: arithmetic means. Replaced in the paper by `iteration_gain.pdf`.

# %%
def convergence_table():
    def geo_red(pct):
        surv = np.clip(1.0 - np.asarray(pct, dtype=float) / 100.0, 1e-9, None)
        return 100.0 * (1.0 - float(np.exp(np.log(surv).mean())))

    def rows(prefix, datasets, avg):
        fix = csv(f"{prefix}_filter_fixpoint.csv")
        fix = fix[fix.dataset.isin(datasets) & fix.reduction_vs_semijoin.notna()]
        r = csv(f"{prefix}_filter_fixpoint_rounds.csv")
        r = r[r.dataset.isin(datasets)].assign(pct=lambda x: 100.0 - x.pct_of_sj_baseline)
        keys = ["dataset", "pattern_name", "graphlet_size"]
        r = r.merge(fix[keys])
        wide = r.pivot_table(index=keys, columns="round", values="pct").reindex(columns=range(1, 61)).ffill(axis=1)
        final = fix.set_index(keys).reduction_vs_semijoin * 100.0
        out = []
        for ds in datasets:
            w = wide.loc[ds]
            f = final.loc[ds]
            gain = f > 1e-9
            cap3 = (100.0 * w[3][gain.values] / f[gain]).mean()
            out.append({"dataset": LABEL[ds], "round 1": avg(w[1]), "round 2": avg(w[2]), "round 3": avg(w[3]), "fixpoint": avg(f),
                        "mean rounds": fix[fix.dataset == ds].n_rounds.mean(), "median rounds": fix[fix.dataset == ds].n_rounds.median(),
                        "max rounds": fix[fix.dataset == ds].n_rounds.max(), "with gain": f"{int(gain.sum())}/{len(f)}",
                        "captured by round 3 (%)": cap3})
        return pd.DataFrame(out).set_index("dataset")

    for bench, prefix, datasets, avg in (("imsm", "labelled", IMSM, geo_red), ("mccreesh", "unlabelled", MCC, np.mean)):
        print(f"\n{BENCH[bench]}:\n" + rows(prefix, datasets, avg).to_string(float_format=lambda v: f"{v:.1f}"))


convergence_table()
