#!/usr/bin/env bash
set -u -o pipefail

# Usage:
#   bash run_weak_scaling.sh [NX] [REPEAT]
#
# Example:
#   bash run_weak_scaling.sh 500 5

NX=${1:-500}
REPEAT=${2:-5}

OMP_LIST="${OMP_LIST:-1 2 4 8 16}"
MPI_LIST="${MPI_LIST:-1 2 4 8 16}"

SERIAL_EXE="${SERIAL_EXE:-./serial_scaling.exe}"
OMP_EXE="${OMP_EXE:-./omp_scaling.exe}"
MPI_EXE="${MPI_EXE:-./mpi_scaling.exe}"

OUTDIR="${OUTDIR:-results/weak_scaling}"
mkdir -p "${OUTDIR}"

RAW_CSV="${OUTDIR}/weak_scaling_raw.csv"
LOGFILE="${OUTDIR}/weak_scaling.log"
FAILFILE="${OUTDIR}/weak_scaling_failed.log"

echo "mode,p,run_id,nx,copies_y,ny,wall_seconds,status" > "${RAW_CSV}"
: > "${LOGFILE}"
: > "${FAILFILE}"

extract_field() {
    local key="$1"
    awk -v k="$key" '{
        for(i=1;i<=NF;i++){
            split($i,a,"=");
            if(a[1]==k){ print a[2]; exit }
        }
    }'
}

run_and_record() {
    local mode="$1"
    local p="$2"
    local run_id="$3"
    local nx="$4"
    local copies_y="$5"
    shift 5

    local ny=$((197 * copies_y))

    echo "Running mode=${mode} p=${p} run=${run_id} nx=${nx} copies_y=${copies_y} ny=${ny}" | tee -a "${LOGFILE}"

    local output
    local rc=0

    if output=$("$@" 2>&1); then
        echo "${output}" >> "${LOGFILE}"

        local wall
        wall=$(echo "${output}" | extract_field wall_seconds)

        if [[ -z "${wall}" ]]; then
            echo "FAILED: missing wall_seconds | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
            echo "${mode},${p},${run_id},${nx},${copies_y},${ny},,missing_wall_seconds" >> "${RAW_CSV}"
        else
            echo "${mode},${p},${run_id},${nx},${copies_y},${ny},${wall},ok" >> "${RAW_CSV}"
        fi
    else
        rc=$?
        echo "${output}" >> "${LOGFILE}"
        echo "FAILED: exit_code=${rc} | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},${nx},${copies_y},${ny},,failed_exit_${rc}" >> "${RAW_CSV}"
    fi
}

echo "[checkpoint] start serial" | tee -a "${LOGFILE}"
for r in $(seq 1 "${REPEAT}"); do
    run_and_record "serial" 1 "${r}" "${NX}" 1 \
        "${SERIAL_EXE}" "${NX}" 1
done
echo "[checkpoint] finished serial" | tee -a "${LOGFILE}"

echo "[checkpoint] start omp" | tee -a "${LOGFILE}"
for p in ${OMP_LIST}; do
    for r in $(seq 1 "${REPEAT}"); do
        run_and_record "omp" "${p}" "${r}" "${NX}" "${p}" \
            env OMP_NUM_THREADS="${p}" "${OMP_EXE}" "${NX}" "${p}"
    done
done
echo "[checkpoint] finished omp" | tee -a "${LOGFILE}"

echo "[checkpoint] start mpi" | tee -a "${LOGFILE}"
for p in ${MPI_LIST}; do
    for r in $(seq 1 "${REPEAT}"); do
        run_and_record "mpi" "${p}" "${r}" "${NX}" "${p}" \
            mpirun -np "${p}" "${MPI_EXE}" "${NX}" "${p}"
    done
done
echo "[checkpoint] finished mpi" | tee -a "${LOGFILE}"

echo "Done."
echo "Raw data:   ${RAW_CSV}"
echo "Main log:   ${LOGFILE}"
echo "Fail log:   ${FAILFILE}"