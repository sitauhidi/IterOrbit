#!/bin/bash -l
#SBATCH --job-name=count
#SBATCH --partition=k2-living-labs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=16
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=00:31:00
#SBATCH --array=0-799
#SBATCH --output=/dev/null
#
# Counting end to end (Section "Counting and Induced Matching"): IterOrbit counts every embedding of one
# instance per array task.
#     sbatch count.sbatch <bench> [<variant>]
#     sbatch --array=0-479 count.sbatch mccreesh [<variant>]     (the unlabelled benchmark has 480 pairs)
#   <bench>    imsm (800 labelled patterns) | mccreesh (480 unlabelled pairs)
#   <variant>  default   IterOrbit as in the paper (the default)
#              nosym     without symmetry breaking
#              noedge    without the edge-triangle filter on the search index
#              noacc     without both ("no acc.")
#              induced   induced counting (mccreesh only)
# Log: RESULTS/count/<bench>/<variant>/<dataset>/<pattern>.out, ending with STATUS (COMPLETED, TIMEOUT or
# ERROR) and ELAPSED_MS (wall clock, reading the graphs included). A task past the benchmark's last
# instance exits at once.

# ---- settings ----------------------------------------------------------------------------------------
WORKSPACE=/mnt/scratch2/users/40334799/oex
RESULTS=${WORKSPACE}/results/paper
TIMEOUT_S=1800
export OMP_NUM_THREADS=1
IO_IMSM="--iterate 3"
IO_MCC="--graphlet-size 5 --iterate 0"
IO_COUNT="--count-only"

BENCH=$1
VARIANT=${2:-default}
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
        IO_CONFIG=${IO_IMSM} ;;
    mccreesh)
        MCC=${WORKSPACE}/dataset-mccreesh-graph
        ROW=$(sed -n "$((T + 2))p" "${MCC}/selected_profiling_pairs.csv")
        [ -n "${ROW}" ] || exit 0
        IFS=, read -r DATASET PATTERN _ _ _ _ PAT_PATH DATA_PATH <<< "${ROW}"
        PAT=${MCC}/${PAT_PATH#*dataset-mccreesh-graph/}
        DATA=${MCC}/${DATA_PATH#*dataset-mccreesh-graph/}
        IO_CONFIG=${IO_MCC} ;;
    *) echo "usage: sbatch count.sbatch imsm|mccreesh [<variant>]" >&2; exit 1 ;;
esac
case ${VARIANT} in
    default) EXTRA=() ;;
    nosym)   EXTRA=(--no-symmetry) ;;
    noedge)  EXTRA=(--no-edge-orbit-filter) ;;
    noacc)   EXTRA=(--no-symmetry --no-edge-orbit-filter) ;;
    induced) [ "${BENCH}" = mccreesh ] || { echo "induced runs on mccreesh only" >&2; exit 1; }
             EXTRA=(--induced) ;;
    *) echo "unknown variant ${VARIANT}" >&2; exit 1 ;;
esac

LOG=${RESULTS}/count/${BENCH}/${VARIANT}/${DATASET}/${PATTERN}.out
mkdir -p "$(dirname "${LOG}")"
exec > "${LOG}" 2>&1
echo "  Variant:  ${VARIANT}"
echo "  Bench:    ${BENCH}"
echo "  Dataset:  ${DATASET}"
echo "  Pattern:  ${PATTERN}"
echo "  Host:     $(hostname)"
echo "  Date:     $(date)"
echo "  Timeout:  ${TIMEOUT_S} s"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "Command: ${IO_BIN} --data ${DATA} --pattern ${PAT} ${IO_CONFIG} ${IO_COUNT} ${EXTRA[*]}"
echo ""
START=$(date +%s%N)
timeout --kill-after=10s "${TIMEOUT_S}s" "${IO_BIN}" --data "${DATA}" --pattern "${PAT}" ${IO_CONFIG} ${IO_COUNT} "${EXTRA[@]}"
RC=$?
echo ""
case ${RC} in
    0)       echo "STATUS: COMPLETED" ;;
    124|137) echo "STATUS: TIMEOUT" ;;
    *)       echo "STATUS: ERROR (exit code ${RC})" ;;
esac
echo "ELAPSED_MS: $(( ($(date +%s%N) - START) / 1000000 ))"
exit 0