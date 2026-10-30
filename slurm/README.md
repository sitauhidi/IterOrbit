# Running the experiments

Four SLURM scripts run IterOrbit's experiments in the paper. Each script is standalone: its settings (paths, the time limit and IterOrbit's flags) are at its top, and it takes no environment variables. Each run writes one log per instance.

| Script | What it runs |
|---|---|
| `build.sbatch` | Compiles IterOrbit into `IterOrbit/build/iterorbit`, the binary the experiments run. One job. |
| `filter.sbatch <bench>` | Filtering only: seven runs per instance, g = 3, 4, 5 as a single pass and to the fixpoint, plus the fixpoint without orbits. |
| `count.sbatch <bench> [<variant>]` | Counting end to end. Variants: `default` (as in the paper), `nosym`, `noedge` and `noacc` (the search-acceleration ablation), and `induced` (unlabelled only). |
| `enum.sbatch <bench>` | Enumeration of the first 10^5 embeddings. |

`<bench>` is `imsm` (800 labelled patterns) or `mccreesh` (480 unlabelled pairs). Every array task runs one instance. The arrays are `0-799`, and the unlabelled benchmark needs `--array=0-479`; a task past the last instance exits at once.

## Setup

Edit the settings block at the top of each script. `WORKSPACE` is the folder that holds this repository and the two benchmarks:

```
<workspace>/
  IterOrbit/                this repository
  dataset-imsm/             8 labelled data graphs, patterns, selected_profiling_patterns.csv (800)
  dataset-mccreesh-graph/   unlabelled pairs, selected_profiling_pairs.csv (480)
  results/paper/            output
```

## Running

Submit from this folder. The cluster limits how many tasks one user can have queued, so submit a few arrays at a time.

```bash
b=$(sbatch --parsable build.sbatch)
D=--dependency=afterok:$b
for bench in imsm mccreesh; do
    A=--array=0-799; [ $bench = mccreesh ] && A=--array=0-479
    sbatch $D $A filter.sbatch $bench
    for v in default nosym noedge noacc; do
        sbatch $D $A count.sbatch $bench $v
    done
    sbatch $D $A enum.sbatch $bench
done
sbatch $D --array=0-479 count.sbatch mccreesh induced
```

Logs go to:
- `results/paper/filter/<bench>/<configuration>/<dataset>/<pattern>.out`, where the configuration is `gs4_singlepass`, `gs4_fixpoint`, `gs4_fixpoint_no_orbit`, and so on;
- `results/paper/count/<bench>/<variant>/<dataset>/<pattern>.out`;
- `results/paper/enum/<bench>/<dataset>/<pattern>.out`.

Each log records the command. It ends with `STATUS` (`COMPLETED`, `TIMEOUT` or `ERROR`) and `ELAPSED_MS`. A timeout is a recorded outcome, so every task ends as COMPLETED in SLURM. `build.sbatch` writes `build_<job>.out` in this folder.

## Conventions

- **Resources:** one NUMA domain per task (16 cores, 120 GB). `OMP_NUM_THREADS=1` caps IterOrbit's only OpenMP section, in g=5 orbit counting.
- **Time limit:** 30 minutes per run.
- **IterOrbit's configuration:** labelled g=4 and I_max=3; unlabelled g=5 to the fixpoint. The scripts pass only the flags that differ from the binary's defaults. So every run induces the candidate subgraph after every stage, as in Algorithm 1, counts orbits on the whole candidate subgraph (no edge thinning), builds no label index, and orders with `core-orbit`.
