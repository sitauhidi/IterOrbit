#!/bin/bash -l
#SBATCH --job-name=filter
#SBATCH --partition=k2-living-labs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=16
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --array=0-799
#SBATCH --output=/dev/null
#
# Filtering only (Sections "Filtering Strength and Cost" and "What the Orbit Stage Contributes"): IterOrbit's
# filter on one instance per array task, stopped after filtering. Seven runs per task, one after another:
# graphlet size g = 3, 4, 5, each as a single pass (I_max = 1) and iterated to the fixpoint, and the
# fixpoint without the orbit stage at the end-to-end graphlet size (g = 4 labelled, 5 unlabelled).
#     sbatch filter.sbatch <bench>
#     sbatch --array=0-479 filter.sbatch mccreesh     (the unlabelled benchmark has 480 pairs)
#   <bench>   imsm (800 labelled patterns) | mccreesh (480 unlabelled pairs)
# Log: one per configuration, RESULTS/filter/<bench>/<configuration>/<dataset>/<pattern>.out
# (gs<g>_singlepass, gs<g>_fixpoint, gs<g>_fixpoint_no_orbit), each ending with STATUS (COMPLETED, TIMEOUT
# or ERROR) and ELAPSED_MS. A task past the benchmark's last instance exits at once.

# ---- settings ----------------------------------------------------------------------------------------
WORKSPACE=/mnt/scratch2/users/40334799/oex
RESULTS=${WORKSPACE}/results/paper
TIMEOUT_S=1800
export OMP_NUM_THREADS=1
IO_FILTER="--filter-only"
# ------------------------------------------------------------------------------------------------------
BENCH=$1
IO_BIN=${WORKSPACE}/IterOrbit/build/iterorbit     # compiled by build.sbatch
T=${SLURM_ARRAY_TASK_ID:?run this with sbatch}

# ---- the instance: row T of the benchmark's selection --------------------------------------------------
case ${BENCH} in
    imsm)
        ROW=$(sed -n "$((T + 2))p" "${WORKSPACE}/dataset-imsm/selected_profiling_patterns.csv")
        [ -n "${ROW}" ] || exit 0
        IFS=, read -r DATASET PATTERN _ <<< "${ROW}"
        PATTERN=${PATTERN%.graph}
        DATA=${WORKSPACE}/dataset-imsm/data/${DATASET}.graph
        PAT=${WORKSPACE}/dataset-imsm/patterns/${DATASET}/${PATTERN}.graph
        G_E2E=4 ;;
    mccreesh)
        MCC=${WORKSPACE}/dataset-mccreesh-graph
        ROW=$(sed -n "$((T + 2))p" "${MCC}/selected_profiling_pairs.csv")
        [ -n "${ROW}" ] || exit 0
        IFS=, read -r DATASET PATTERN _ _ _ _ PAT_PATH DATA_PATH <<< "${ROW}"
        PAT=${MCC}/${PAT_PATH#*dataset-mccreesh-graph/}
        DATA=${MCC}/${DATA_PATH#*dataset-mccreesh-graph/}
        G_E2E=5 ;;
    *) echo "usage: sbatch filter.sbatch imsm|mccreesh" >&2; exit 1 ;;
esac

run() {   # run <configuration> <IterOrbit flags ...>: one filtering run, logged in that configuration's folder
    local config=$1
    local log=${RESULTS}/filter/${BENCH}/${config}/${DATASET}/${PATTERN}.out
    shift
    mkdir -p "$(dirname "${log}")"
    {
        echo "  Configuration: ${config}"
        echo "  Bench:    ${BENCH}"
        echo "  Dataset:  ${DATASET}"
        echo "  Pattern:  ${PATTERN}"
        echo "  Host:     $(hostname)"
        echo "  Date:     $(date)"
        echo "  Timeout:  ${TIMEOUT_S} s"
        echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
        echo "Command: ${IO_BIN} --data ${DATA} --pattern ${PAT} ${IO_FILTER} $*"
        echo ""
        local start=$(date +%s%N)
        timeout --kill-after=10s "${TIMEOUT_S}s" "${IO_BIN}" --data "${DATA}" --pattern "${PAT}" ${IO_FILTER} "$@"
        local rc=$?
        echo ""
        case ${rc} in
            0)       echo "STATUS: COMPLETED" ;;
            124|137) echo "STATUS: TIMEOUT" ;;
            *)       echo "STATUS: ERROR (exit code ${rc})" ;;
        esac
        echo "ELAPSED_MS: $(( ($(date +%s%N) - start) / 1000000 ))"
    } > "${log}" 2>&1
}

for g in 3 4 5; do
    run gs${g}_singlepass --graphlet-size ${g}
    run gs${g}_fixpoint --graphlet-size ${g} --iterate 0
done
run gs${G_E2E}_fixpoint_no_orbit --graphlet-size ${G_E2E} --iterate 0 --no-orbit-filter
exit 0