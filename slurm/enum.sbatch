#!/bin/bash -l
#SBATCH --job-name=enum
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
# Enumeration end to end (Section "Enumeration End to End"): IterOrbit lists the embeddings of one instance
# one at a time and stops after the first 10^5 (Sun and Luo's protocol), one instance per array task, with
# no symmetry breaking or last-level counting.
#     sbatch enum.sbatch <bench>
#     sbatch --array=0-479 enum.sbatch mccreesh     (the unlabelled benchmark has 480 pairs)
#   <bench>   imsm (800 labelled patterns) | mccreesh (480 unlabelled pairs)
# Log: RESULTS/enum/<bench>/<dataset>/<pattern>.out, ending with STATUS (COMPLETED, TIMEOUT or ERROR) and
# ELAPSED_MS (wall clock, reading the graphs included). A task past the benchmark's last instance exits
# at once.

# ---- settings ----------------------------------------------------------------------------------------
WORKSPACE=/mnt/scratch2/users/40334799/oex
RESULTS=${WORKSPACE}/results/paper
TIMEOUT_S=1800
LIMIT=100000                  # embeddings listed per instance
export OMP_NUM_THREADS=1
IO_IMSM="--graphlet-size 4 --iterate 3"
IO_MCC="--graphlet-size 5 --iterate 0"
IO_ENUM="--max-embeddings ${LIMIT} --order-mode core-orbit --no-label-index --no-edge-thinning"
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
        IO_CONFIG=${IO_IMSM} ;;
    mccreesh)
        MCC=${WORKSPACE}/dataset-mccreesh-graph
        ROW=$(sed -n "$((T + 2))p" "${MCC}/selected_profiling_pairs.csv")
        [ -n "${ROW}" ] || exit 0
        IFS=, read -r DATASET PATTERN _ _ _ _ PAT_PATH DATA_PATH <<< "${ROW}"
        PAT=${MCC}/${PAT_PATH#*dataset-mccreesh-graph/}
        DATA=${MCC}/${DATA_PATH#*dataset-mccreesh-graph/}
        IO_CONFIG=${IO_MCC} ;;
    *) echo "usage: sbatch enum.sbatch imsm|mccreesh" >&2; exit 1 ;;
esac

LOG=${RESULTS}/enum/${BENCH}/${DATASET}/${PATTERN}.out
mkdir -p "$(dirname "${LOG}")"
exec > "${LOG}" 2>&1
echo "  Bench:    ${BENCH}"
echo "  Dataset:  ${DATASET}"
echo "  Pattern:  ${PATTERN}"
echo "  Host:     $(hostname)"
echo "  Date:     $(date)"
echo "  Timeout:  ${TIMEOUT_S} s"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "Command: ${IO_BIN} --data ${DATA} --pattern ${PAT} ${IO_CONFIG} ${IO_ENUM}"
echo ""
START=$(date +%s%N)
timeout --kill-after=10s "${TIMEOUT_S}s" "${IO_BIN}" --data "${DATA}" --pattern "${PAT}" ${IO_CONFIG} ${IO_ENUM}
RC=$?
echo ""
case ${RC} in
    0)       echo "STATUS: COMPLETED" ;;
    124|137) echo "STATUS: TIMEOUT" ;;
    *)       echo "STATUS: ERROR (exit code ${RC})" ;;
esac
echo "ELAPSED_MS: $(( ($(date +%s%N) - START) / 1000000 ))"
exit 0
