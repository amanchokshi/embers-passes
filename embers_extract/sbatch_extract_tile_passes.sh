#!/bin/bash

#SBATCH --job-name='extract-passes'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0:16:00
#SBATCH --account=oz048
#SBATCH --mem=4gb
#SBATCH --output=extract_passes-%A.out
#SBATCH --error=extract_passes-%A.err

set -euo pipefail

CONTAINER=/fred/oz048/achokshi/embers/embers_1.0.1_py313.sif
WORKDIR=/fred/oz048/achokshi/embers

module load apptainer

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd "${WORKDIR}"

echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Container: ${CONTAINER}"

apptainer exec \
    -B /fred:/fred \
    "${CONTAINER}" \
    python extract_tile_passes.py

echo "Job finished at: $(date)"
