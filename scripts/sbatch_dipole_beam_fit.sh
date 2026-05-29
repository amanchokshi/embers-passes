#!/bin/bash
#SBATCH --job-name=embers-smc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=logs/hyperbeam_smc_%A_%a.out
#SBATCH --error=logs/hyperbeam_smc_%A_%a.err

set -euo pipefail
trap 'echo "FAILED at line ${LINENO}" >&2' ERR

mkdir -p logs

module purge || true
module use /project/rrg-sievers/achokshi/software/modulefiles
module load embers-passes

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

DRAWS=1024
NSEEDS=16
RAYON_THREADS=3

NSIDE=32
FREQ_HZ=137000000
POINTING=0

PASS_DIR="../passes/rf0"
BEAM_FILE="${MWA_BEAM_FILE}"
GROUP_FILE="tile_groups.txt"

mapfile -t TILE_GROUPS < "${GROUP_FILE}" || true

if (( SLURM_ARRAY_TASK_ID >= ${#TILE_GROUPS[@]} )); then
    echo "Invalid array index ${SLURM_ARRAY_TASK_ID}; only ${#TILE_GROUPS[@]} tile groups."
    exit 1
fi

read -ra TILES <<< "${TILE_GROUPS[$SLURM_ARRAY_TASK_ID]}"

BASE_TMP="${SLURM_TMPDIR:-/scratch/achokshi}/pytensor_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${BASE_TMP}"

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Array job: ${SLURM_ARRAY_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Tiles: ${TILES[*]}"
echo "Draws: ${DRAWS}"
echo "Seeds per dataset: ${NSEEDS}"
echo "RAYON_NUM_THREADS per process: ${RAYON_THREADS}"
echo "PyTensor cache: ${BASE_TMP}"
echo "Python: $(which python)"
python --version

for tile in "${TILES[@]}"; do
    for POL in XX YY; do
        pol_lower=$(echo "${POL}" | tr '[:upper:]' '[:lower:]')

        pass_file="${PASS_DIR}/S${tile}${POL}_rf0${POL}_passes.h5"
        dataset_out="/scratch/achokshi/embers-passes/data/SMC/rf0/S${tile}${POL}"

        if [[ ! -f "${pass_file}" ]]; then
            echo "Missing pass file: ${pass_file}" >&2
            exit 1
        fi

        mkdir -p "${dataset_out}"

        echo "Launching S${tile}${POL}"
        echo "  Pass file: ${pass_file}"
        echo "  Output:    ${dataset_out}"

        for i in $(seq 0 $((NSEEDS - 1))); do
            seed=$((42 + i))

            (
                export RAYON_NUM_THREADS="${RAYON_THREADS}"
                export PYTENSOR_FLAGS="base_compiledir=${BASE_TMP}/S${tile}${POL}_seed_${seed}"

                mkdir -p "${BASE_TMP}/S${tile}${POL}_seed_${seed}"

                echo "[S${tile}${POL} seed=${seed}] starting $(date)"

                /usr/bin/time -v python dipole_beam_fit.py \
                    --seed "${seed}" \
                    --draws "${DRAWS}" \
                    --outdir "${dataset_out}" \
                    --pass-file "${pass_file}" \
                    --beam-file "${BEAM_FILE}" \
                    --nside "${NSIDE}" \
                    --freq-hz "${FREQ_HZ}" \
                    --pointing "${POINTING}" \
                    --pol "${pol_lower}"

                touch "${dataset_out}/seed_${seed}.done"

                echo "[S${tile}${POL} seed=${seed}] finished $(date)"
            ) > "logs/S${tile}${POL}_seed_${seed}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
              2> "logs/S${tile}${POL}_seed_${seed}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" &
        done
    done
done

wait

echo "All runs complete: $(date)"
