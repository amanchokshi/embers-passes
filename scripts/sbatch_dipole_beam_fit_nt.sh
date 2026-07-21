#!/bin/bash
#SBATCH --job-name=embers-smc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --time=3:00:00
#SBATCH --tmp=1G
#SBATCH --output=../data/SMC/rf0_gaussian/logs/dipole_beam_fit_%A_%a.out
#SBATCH --error=../data/SMC/rf0_gaussian/logs/dipole_beam_fit_%A_%a.err

set -euo pipefail
trap 'echo "FAILED at line ${LINENO}" >&2' ERR

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 TILEPOL" >&2
    echo "Example: sbatch $0 S06XX" >&2
    exit 2
fi

TILEPOL="$1"

if [[ ! "${TILEPOL}" =~ ^S[0-9]+(XX|YY)$ ]]; then
    echo "Invalid tile/polarization: ${TILEPOL}" >&2
    echo "Expected a value such as S06XX or S36YY." >&2
    exit 2
fi

POL="${TILEPOL: -2}"
POL_LOWER="${POL,,}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "This script must be run as a Slurm array job." >&2
    echo "Example: sbatch --array=0-15 $0 S06XX" >&2
    exit 1
fi

# Array tasks 0--15 correspond to seeds 42--57.
SEED=$((42 + SLURM_ARRAY_TASK_ID))

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------

set +e
module purge
module use /fred/oz048/achokshi/software/modulefiles
module load embers-passes/0.1.0.ac
set -e

# Prevent NumPy and related libraries from creating extra thread pools.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Hyperbeam/Rayon may use all CPUs allocated to this array task.
export RAYON_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DRAWS=1024
NSIDE=32
FREQ_HZ=137000000
POINTING=0

PASS_DIR="../passes/monotonic_concave_cal/rf0"
PASS_FILE="${PASS_DIR}/${TILEPOL}_rf0${POL}_passes.h5"

MWA_BEAM_FILE="/fred/oz048/achokshi/software/mwa_beams/mwa_full_embedded_element_pattern.h5"

OUTPUT_ROOT="/fred/oz048/achokshi/embers/embers-passes/data/SMC/rf0_gaussian"
DATASET_OUT="${OUTPUT_ROOT}/${TILEPOL}"
DONE_FILE="${DATASET_OUT}/seed_${SEED}.done"

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------

if [[ ! -f "${PASS_FILE}" ]]; then
    echo "Missing pass file: ${PASS_FILE}" >&2
    exit 1
fi

if [[ ! -f "${MWA_BEAM_FILE}" ]]; then
    echo "Missing beam file: ${MWA_BEAM_FILE}" >&2
    exit 1
fi

if [[ -z "${JOBFS:-}" ]]; then
    echo "JOBFS is not defined." >&2
    exit 1
fi

mkdir -p "${DATASET_OUT}"

# This makes it safe to resubmit the full array after partial completion.
if [[ -f "${DONE_FILE}" ]]; then
    echo "${TILEPOL}, seed ${SEED} is already complete."
    exit 0
fi

# ---------------------------------------------------------------------------
# Node-local files
# ---------------------------------------------------------------------------

BEAM_BASENAME="$(basename "${MWA_BEAM_FILE}")"
LOCAL_BEAM_FILE="${JOBFS}/${BEAM_BASENAME}"

echo "Copying beam model to node-local JOBFS:"
echo "  ${MWA_BEAM_FILE}"
echo "  -> ${LOCAL_BEAM_FILE}"

cp "${MWA_BEAM_FILE}" "${LOCAL_BEAM_FILE}"

PYTENSOR_DIR="${JOBFS}/pytensor_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${PYTENSOR_DIR}"

export PYTENSOR_FLAGS="base_compiledir=${PYTENSOR_DIR}"

# ---------------------------------------------------------------------------
# Job information
# ---------------------------------------------------------------------------

echo "Started:            $(date)"
echo "Host:               $(hostname)"
echo "Slurm job ID:       ${SLURM_JOB_ID}"
echo "Array job ID:       ${SLURM_ARRAY_JOB_ID}"
echo "Array task ID:      ${SLURM_ARRAY_TASK_ID}"
echo "Tile/polarization:  ${TILEPOL}"
echo "Polarization:       ${POL_LOWER}"
echo "Seed:               ${SEED}"
echo "Draws:              ${DRAWS}"
echo "Allocated CPUs:     ${SLURM_CPUS_PER_TASK}"
echo "Rayon threads:      ${RAYON_NUM_THREADS}"
echo "Pass file:          ${PASS_FILE}"
echo "Beam file:          ${LOCAL_BEAM_FILE}"
echo "Output directory:   ${DATASET_OUT}"
echo "PyTensor directory: ${PYTENSOR_DIR}"
echo "Python:             $(which python)"
python --version

# ---------------------------------------------------------------------------
# Run one tile/polarization and one seed
# ---------------------------------------------------------------------------

/usr/bin/time -v python dipole_beam_fit.py \
    --seed "${SEED}" \
    --draws "${DRAWS}" \
    --outdir "${DATASET_OUT}" \
    --pass-file "${PASS_FILE}" \
    --beam-file "${LOCAL_BEAM_FILE}" \
    --nside "${NSIDE}" \
    --freq-hz "${FREQ_HZ}" \
    --pointing "${POINTING}" \
    --pol "${POL_LOWER}"

touch "${DONE_FILE}"

echo "Finished: $(date)"
