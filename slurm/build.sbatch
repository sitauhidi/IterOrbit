#!/bin/bash -l
#SBATCH --job-name=build
#SBATCH --partition=k2-living-labs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=16
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --output=build_%j.out
#
# Compiles IterOrbit on a compute node (-march=native: build where it runs) into IterOrbit/build/iterorbit,
# the binary the experiments run. Submit it first, from this folder:
#     sbatch build.sbatch
# and make the experiment arrays wait for it (--dependency=afterok:<job id>; see README.md). The log is
# build_<job id>.out in this folder.

# ---- settings ----------------------------------------------------------------------------------------
WORKSPACE=/mnt/scratch2/users/40334799/oex     # holds IterOrbit/, dataset-imsm/, dataset-mccreesh-graph/
GCC_DIR=/opt/gridware/depots/54e7fb3c/el8/pkg/compilers/gcc/14.1.0         # as IterOrbit/Makefile
# ------------------------------------------------------------------------------------------------------
set -e
echo "host $(hostname)  $(date)"
cd "${WORKSPACE}/IterOrbit"
make -j16 GCC_DIR="${GCC_DIR}"
echo "BUILD DONE: ${WORKSPACE}/IterOrbit/build/iterorbit"
