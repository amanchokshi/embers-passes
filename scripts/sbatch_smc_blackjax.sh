#!/bin/bash
#SBATCH --job-name=smc-blackjax
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --time=10:00:00
#SBATCH --output=logs/smc_blackjax_%j.out
#SBATCH --error=logs/smc_blackjax_%j.err

set -uo pipefail

mkdir -p logs
mkdir -p ../data/smc_nuts_blackjax

module purge || true
module use /project/rrg-sievers/achokshi/software/modulefiles
module load embers-passes


# -------------------------------------------------------------------------
# Shared settings
# -------------------------------------------------------------------------

PASS_FILE="../passes/monotonic_concave_5hz_pval/S06XX_rf0XX_passes.h5"
OUTPUT_DIR="../data/smc_nuts_blackjax"

N_PARTICLES=3072
N_DEVICES=48
SEED=314159
CHUNK_DAYS=20
MAX_DOUBLINGS=10


# Avoid hidden CPU oversubscription from BLAS/OpenMP libraries.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1


# -------------------------------------------------------------------------
# Launch four independent 48-core JAX processes.
#
# Whole-node CPU layout:
#
#   chunk 0 -> CPUs   0-47
#   chunk 1 -> CPUs  48-95
#   chunk 2 -> CPUs  96-143
#   chunk 3 -> CPUs 144-191
#
# Each process exposes 48 JAX CPU devices and handles
#
#     3072 / 48 = 64 particles per device.
# -------------------------------------------------------------------------

pids=()
chunks=()

for CHUNK_INDEX in 0 1 2 3
do
    CPU_START=$(( CHUNK_INDEX * N_DEVICES ))
    CPU_END=$(( CPU_START + N_DEVICES - 1 ))

    TAG="N3072_S${SEED}_C20_I${CHUNK_INDEX}_AMP_COMPLEX"

    STDOUT_LOG="logs/${TAG}.out"
    STDERR_LOG="logs/${TAG}.err"

    CACHE_ROOT="${SLURM_TMPDIR}/embers_smc_chunk_${CHUNK_INDEX}"

    mkdir -p \
        "${CACHE_ROOT}/matplotlib" \
        "${CACHE_ROOT}/numba"

    echo "Launching chunk ${CHUNK_INDEX}"
    echo "  CPUs:   ${CPU_START}-${CPU_END}"
    echo "  tag:    ${TAG}"
    echo "  stdout: ${STDOUT_LOG}"
    echo "  stderr: ${STDERR_LOG}"

    taskset -c "${CPU_START}-${CPU_END}" \
        env \
            XDG_CACHE_HOME="${CACHE_ROOT}" \
            MPLCONFIGDIR="${CACHE_ROOT}/matplotlib" \
            NUMBA_CACHE_DIR="${CACHE_ROOT}/numba" \
            XLA_FLAGS="--xla_force_host_platform_device_count=${N_DEVICES}" \
            python -u run_smc_blackjax.py \
                "${PASS_FILE}" \
                --kind both \
                --particles "${N_PARTICLES}" \
                --devices "${N_DEVICES}" \
                --seed "${SEED}" \
                --max-doublings "${MAX_DOUBLINGS}" \
                --chunk-days "${CHUNK_DAYS}" \
                --chunk-index "${CHUNK_INDEX}" \
		--single-pass-z-cut -3 \
                --output-dir "${OUTPUT_DIR}" \
                --tag "${TAG}" \
        > "${STDOUT_LOG}" \
        2> "${STDERR_LOG}" &

    pids+=("$!")
    chunks+=("${CHUNK_INDEX}")
done


# -------------------------------------------------------------------------
# Wait for all four processes and report failures independently.
# -------------------------------------------------------------------------

echo
echo "Waiting for chunk runs..."

failed=0

for i in "${!pids[@]}"
do
    pid="${pids[$i]}"
    chunk="${chunks[$i]}"

    if wait "${pid}"; then
        echo "Chunk ${chunk}: SUCCESS"
    else
        status=$?
        echo "Chunk ${chunk}: FAILED (exit code ${status})" >&2
        failed=1
    fi
done

echo
echo "All chunk processes have exited."
echo "Check individual logs in ./logs/"

exit "${failed}"
