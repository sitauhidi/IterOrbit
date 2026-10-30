# IterOrbit

**IterOrbit** is a subgraph isomorphism engine: given a small
*pattern* graph and a large *data* graph, both vertex-labelled and undirected, it finds or counts
the embeddings of the pattern in the data graph. It extends OrbitSI's graphlet-orbit filter by
iterating it: the candidates that survive a round induce a smaller data graph, and the filters
are evaluated again on that graph, until nothing more is removed.

The executable is called `iterorbit`.

## Requirements

- A C++17 compiler with OpenMP (developed with GCC 14.1) and `make`.
- Python 3 with `pytest`, for the tests.

The orbit-counting library Escape/EVOKE is included in `extern/Escape`.

## Building

```bash
make                      # -> build/iterorbit
make GCC_DIR=/path/to/gcc # use another GCC installation (the default path is site-specific)
make BUILD_DIR=out        # build into another directory
```

The build uses `-march=native`, so run the binary on the kind of machine it was built on.
`make clean` removes the whole build directory.

## Input format

One graph per file, as text:

```
t <number of vertices> <number of edges>
v <id> <label> <degree>
e <u> <v>
```

Only the `v` lines (vertex ID and label) and `e` lines (one per undirected edge) are read.
Vertex IDs are `0 .. n-1`.

## Usage

```bash
./build/iterorbit --data <data graph> --pattern <pattern graph> [options]
```

Typical runs:

```bash
# Count the embeddings, iterating the filter to a fixpoint (the configuration of the paper's
# end-to-end runs)
./build/iterorbit --data test/data_graph/HPRD.graph \
                --pattern test/query_graph/query_dense_16_104.graph --iterate 0 --count-only

# List every embedding, with a single filtering round
./build/iterorbit --data test/data_graph/HPRD.graph \
                --pattern test/query_graph/query_dense_16_104.graph --verbose

# Filtering only: candidate statistics and timings per stage and round
./build/iterorbit --data test/data_graph/HPRD.graph \
                --pattern test/query_graph/query_dense_16_104.graph --iterate 0 --filter-only
```

### Options

`./build/iterorbit` without arguments prints this list.

Filtering:

| Option | Meaning |
|---|---|
| `--graphlet-size <3\|4\|5>` | Graphlet size for orbit counting (default 4). |
| `--iterate <N>` | Maximum number of filtering rounds, `I_max` (default 1). `0` iterates until a round removes no candidate. |
| `--use-full-graph` | Count orbits on the whole data graph instead of the candidate-induced subgraph. One round only; cannot be combined with `--iterate`. |
| `--filter-only` | Stop after filtering and print candidate statistics. |
| `--no-arc-consistency` | Run the semi-join once per round instead of to a fixpoint. |
| `--disable-semi-join` | Skip the semi-join stage. |
| `--no-orbit-filter` | Skip the orbit stage in every round; the data graph is still contracted. |
| `--edge-thinning` | Count orbits only on the edges of the candidate subgraph that some pattern edge can map to (off by default; `--no-edge-thinning` is accepted and keeps it off). |
| `--label-index` | Build the data graph's neighbour-label index at load, for LDF and NLF (off by default: they count labels from the adjacency; `--no-label-index` is accepted and keeps it off). |
| `--single-induction` | Induce the data graph once per round (before the orbit stage) instead of twice (also after it). |
| `--orbit-throttle`, `--orbit-throttle-threshold <f>` | Skip the orbit stage in rounds where the candidate sets have barely shrunk since it last ran (off by default). |

Ordering and search:

| Option | Meaning |
|---|---|
| `--induced` | Induced subgraph isomorphism (default: non-induced). |
| `--count-only` | Count the embeddings without storing them. Cannot be combined with `--verbose`. |
| `--max-embeddings <N>` | Enumerate the embeddings one by one and stop after `N` of them (no symmetry breaking or last-level counting). Cannot be combined with `--count-only`. |
| `--cbj` | Conflict-directed backjumping. |
| `--order-mode <core-orbit\|core\|orbit>` | Matching order (default `core-orbit`); see `include/OrderEngine.h`. |
| `--order <u,u,...>` | Use this matching order (a permutation of the pattern's vertex IDs). |
| `--no-edge-orbit-filter` | Do not prune the candidate-space index by per-edge triangle counts. |
| `--no-symmetry` | With `--count-only`: count every embedding one by one instead of one per class of pattern automorphisms. |
| `--search-stats` | Print the matching order and per-depth search-tree statistics. |
| `--verbose` | Print every embedding and per-stage filtering statistics. |

The defaults are the configuration the paper describes: the data graph induced twice per round, orbits
counted on the whole candidate subgraph (no edge thinning), no label index, and the `core-orbit` order.
The options that switch a component off exist for ablation; `--edge-thinning` and `--label-index` switch
on optional speed-ups. None of them changes the result.

### Output

The run prints its configuration, then one line per phase and the result:

```
Edge thinning: off; label index: off
...
Total loading time: ... ms.
Total filtering time: ... ms.
Total ordering time: ... us.
Total symmetry time: ... ms.        (with --count-only)
Total searching time: ... ms.
Matches found: <number of embeddings>
Total execution time: ... ms.
```

If filtering empties a candidate set, the pattern has no embedding: the run reports
`Matches found: 0` without searching. A count that exceeds 64 bits is reported as `INT64_MAX`,
with a warning.

## Code map

| Paper | Code |
|---|---|
| Candidate sets, filter operators, Algorithm 1 (IterOrbitFilter) | `include/FilterEngine.h`, `src/FilterEngine.cpp` |
| Graphlet orbit counts | `include/OrbitCounter.h`, `src/OrbitCounter.cpp`, `extern/Escape` |
| Candidate-induced subgraph | `AppGraph::createSubgraph` in `src/AppGraph.cpp` |
| Ordering strategy | `include/OrderEngine.h`, `src/OrderEngine.cpp` |
| Enumeration, edge-orbit index filter | `include/SearchEngine.h`, `src/SearchEngine.cpp` |
| Counting with symmetry breaking | `include/Symmetry.h`, `src/Symmetry.cpp` |
| Command line, phase timers | `src/main.cpp` |

The headers document each component's interface. `src/OrbitCounter.cpp` reimplements several
Escape routines; each says which routine it replaces and what differs.

## API documentation

The headers and file summaries are written for Doxygen. With Doxygen installed, run

```bash
doxygen        # -> docs/html/index.html
```

in this directory (`Doxyfile` holds the configuration).

## Testing

`test/test_orbitsi.py` runs every pattern in `test/query_graph` against `test/data_graph/HPRD.graph`
and compares the number of embeddings with `test/expected_output.res`.

```bash
make test                         # graphlet sizes 3, 4 and 5, default configuration
GRAPHLET_SIZE=5 ITERATE=0 CBJ=1 COUNT_ONLY=1 pytest test/test_orbitsi.py
```

The configuration is set through environment variables: `GRAPHLET_SIZE`, `ITERATE`,
`USE_FULL_GRAPH`, `CBJ`, `COUNT_ONLY`, `ARC_CONSISTENCY` (`0` for `--no-arc-consistency`),
`ORDER_MODE`, `EXTRA_FLAGS` (any other options) and `ORBITSI_BIN` (the binary to test, default
`build/iterorbit`).

`test/orbit_equivalence.cpp` checks the reimplemented orbit-counting routines against Escape's
originals on subgraphs of a given graph; its header comment gives the build line and usage.

## Reproducing the paper's results

`analysis/` holds every number, table and figure of the paper, and the exploratory analyses behind it:

- `analysis/analysis.ipynb`: the notebook that computes everything from the experiments' results, one section
  per figure, table or subsection, with its outputs. `analysis/analysis.py` is the same code as a script; the
  notebook is generated from it. Both read the results as CSV files from `analysis/data/`, which is not part of
  this repository.
- `slurm/`: the SLURM scripts that run IterOrbit's experiments and write those logs (see `slurm/README.md`).
- `baseline/`: the patches the paper's comparison applied to the baselines' upstream sources (`subgraphmatching-mccreesh` and `-filteronly` for SubgraphMatchingSurvey `eb3ab43`, whose own candidate printing is left as upstream has it; `gup-timed` for GuP `fbaf103`), each applied with `patch -p1` to a copy of the checkout.

```bash
pip install -r analysis/requirements.txt
jupyter notebook analysis/analysis.ipynb      # figures shown inline
python analysis/analysis.py                   # or as a script: figures written as PDFs to analysis/figures/
```
