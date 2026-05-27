#!/bin/bash
#SBATCH --job-name=embers-smc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --time=09:00:00
#SBATCH --output=logs/hyperbeam_smc_%j.out
#SBATCH --error=logs/hyperbeam_smc_%j.err

set -euo pipefail

mkdir -p logs

module purge
module use /project/rrg-sievers/achokshi/software/modulefiles
module load embers-passes

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1


BASE_TMP="${SLURM_TMPDIR:-/scratch/achokshi}/pytensor_${SLURM_JOB_ID}"
mkdir -p "${BASE_TMP}"

PASS_FILE="../passes/rf0/S06XX_rf0XX_passes.h5"
BEAM_FILE="${MWA_BEAM_FILE}"

DRAWS=2048
NSEEDS=64
RAYON_THREADS=3

OUTDIR="/scratch/achokshi/embers-passes/data/SMC/N${NSEEDS}_D${DRAWS}_T${RAYON_THREADS}_${SLURM_JOB_ID}"
mkdir -p "${OUTDIR}"

NSIDE=32
FREQ_HZ=137000000
POINTING=0
POL=xx

echo "Job ID: ${SLURM_JOB_ID}"
echo "Output: ${OUTDIR}"
echo "Draws: ${DRAWS}"
echo "Seeds: ${NSEEDS}"
echo "RAYON_NUM_THREADS per process: ${RAYON_THREADS}"

for i in $(seq 0 $((NSEEDS - 1))); do
    seed=$((42 + i))

    (
        export RAYON_NUM_THREADS="${RAYON_THREADS}"
        export PYTENSOR_FLAGS="base_compiledir=${BASE_TMP}/seed_${seed}"
        mkdir -p "${BASE_TMP}/seed_${seed}"

        echo "[seed=${seed}] starting $(date)"

        /usr/bin/time -v python dipole_beam_fit.py \
            --seed "${seed}" \
            --draws "${DRAWS}" \
            --outdir "${OUTDIR}" \
            --pass-file "${PASS_FILE}" \
            --beam-file "${BEAM_FILE}" \
            --nside "${NSIDE}" \
            --freq-hz "${FREQ_HZ}" \
            --pointing "${POINTING}" \
            --pol "${POL}"

        echo "[seed=${seed}] finished $(date)"
    ) > "logs/seed_${seed}_${SLURM_JOB_ID}.out" \
      2> "logs/seed_${seed}_${SLURM_JOB_ID}.err" &
done

wait

echo "All runs complete."
echo "Outputs written to ${OUTDIR}"
