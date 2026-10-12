# OrbitEX (IterOrbit)

OrbitEX is the implementation of **IterOrbit**, a subgraph isomorphism engine: given a small
*pattern* graph and a large *data* graph, both vertex-labelled and undirected, it finds or counts
the embeddings of the pattern in the data graph. It extends OrbitSI's graphlet-orbit filter by
iterating it: the candidates that survive a round induce a smaller data graph, and the filters
are evaluated again on that graph, until nothing more is removed.

The executable is called `orbitsi`.

## Requirements

- A C++17 compiler with OpenMP (developed with GCC 14.1) and `make`.
- Python 3 with `pytest`, for the tests.

The orbit-counting library Escape/EVOKE is included in `extern/Escape`.

## Building

```bash
make                      # -> build/orbitsi
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
./build/orbitsi --data <data graph> --pattern <pattern graph> [options]
```

Typical runs:

```bash
# Count the embeddings, iterating the filter to a fixpoint (the configuration of the paper's
# end-to-end runs)
./build/orbitsi --data test/data_graph/HPRD.graph \
                --pattern test/query_graph/query_dense_16_104.graph --iterate 0 --count-only

# List every embedding, with a single filtering round
./build/orbitsi --data test/data_graph/HPRD.graph \
                --pattern test/query_graph/query_dense_16_104.graph --verbose

# Filtering only: candidate statistics and timings per stage and round
./build/orbitsi --data test/data_graph/HPRD.graph \
                --pattern test/query_graph/query_dense_16_104.graph --iterate 0 --filter-only
```

### Options

`./build/orbitsi` without arguments prints this list.

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
| `--single-induction` | Induce the data graph once per round (before the orbit stage) instead of twice (also after it). |
| `--orbit-throttle`, `--orbit-throttle-threshold <f>` | Skip the orbit stage in rounds where the candidate sets have barely shrunk since it last ran (off by default). |

Ordering and search:

| Option | Meaning |
|---|---|
| `--induced` | Induced subgraph isomorphism (default: non-induced). |
| `--count-only` | Count the embeddings without storing them. Cannot be combined with `--verbose`. |
| `--cbj` | Conflict-directed backjumping. |
| `--order-mode <core\|orbit\|core-orbit>` | Matching order (default `core`); see `include/OrderEngine.h`. |
| `--order <u,u,...>` | Use this matching order (a permutation of the pattern's vertex IDs). |
| `--no-edge-orbit-filter` | Do not prune the candidate-space index by per-edge triangle counts. |
| `--no-symmetry` | With `--count-only`: count every embedding one by one instead of one per class of pattern automorphisms. |
| `--search-stats` | Print the matching order and per-depth search-tree statistics. |
| `--verbose` | Print every embedding and per-stage filtering statistics. |

The options that switch a component off exist for ablation; none of them changes the result.

### Output

The run prints its configuration, then one line per phase and the result:

```
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
`build/orbitsi`).

`test/orbit_equivalence.cpp` checks the reimplemented orbit-counting routines against Escape's
originals on subgraphs of a given graph; its header comment gives the build line and usage.
