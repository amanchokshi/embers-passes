#!/bin/bash
#SBATCH --job-name=emb-smc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=768
#SBATCH --time=12:00:00
#SBATCH --tmp=256
#SBATCH --output=logs/dipole_beam_fit_%j.out
#SBATCH --error=logs/dipole_beam_fit_%j.err

module purge
module use /fred/oz048/achokshi/software/modulefiles
module load embers-passes/0.1.0.ac

# Prevent NumPy and related libraries from independently spawning threads.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export RAYON_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "RAYON_NUM_THREADS: ${RAYON_NUM_THREADS}"
echo "Python: $(command -v python)"
python --version

DRAWS=1024
SEED=44

NSIDE=32
FREQ_HZ=137000000
POINTING=0

cp /fred/oz048/achokshi/software/mwa_beams/mwa_full_embedded_element_pattern.h5 ${JOBFS}
BEAM_FILE=${JOBFS}/mwa_full_embedded_element_pattern.h5

PASS_DIR="../passes/monotonic_concave_cal/rf0"
OUTPUT_ROOT="/fred/oz048/achokshi/embers/embers-passes/data/SMC/rf0_gaussian"

TILE=06
POL=XX
POL_LOWER=xx

PASS_FILE="${PASS_DIR}/S${TILE}${POL}_rf0${POL}_passes.h5"
OUTDIR="${OUTPUT_ROOT}/S${TILE}${POL}"

BASE_TMP="${SLURM_TMPDIR:-/tmp/${USER}}/pytensor_${SLURM_JOB_ID}"
PYTENSOR_CACHE="${BASE_TMP}/S${TILE}${POL}_seed_${SEED}"

mkdir -p "${OUTDIR}" "${PYTENSOR_CACHE}"

if [[ ! -f "${PASS_FILE}" ]]; then
    echo "ERROR: Missing pass file: ${PASS_FILE}" >&2
    exit 1
fi

export PYTENSOR_FLAGS="base_compiledir=${PYTENSOR_CACHE}"

echo "Started:             $(date)"
echo "Host:                $(hostname)"
echo "Job ID:              ${SLURM_JOB_ID}"
echo "Tile:                S${TILE}${POL}"
echo "Seed:                ${SEED}"
echo "Draws:               ${DRAWS}"
echo "Beam file:           ${BEAM_FILE}"
echo "Pass file:           ${PASS_FILE}"
echo "Output:              ${OUTDIR}"
echo "PyTensor cache:      ${PYTENSOR_CACHE}"
echo "Python:              $(which python)"
python --version

/usr/bin/time -v python dipole_beam_fit.py \
    --seed "${SEED}" \
    --draws "${DRAWS}" \
    --outdir "${OUTDIR}" \
    --pass-file "${PASS_FILE}" \
    --beam-file "${BEAM_FILE}" \
    --nside "${NSIDE}" \
    --freq-hz "${FREQ_HZ}" \
    --pointing "${POINTING}" \
    --pol "${POL_LOWER}"

touch "${OUTDIR}/seed_${SEED}.done"

echo "Finished: $(date)"
