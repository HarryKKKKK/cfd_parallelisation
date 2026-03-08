#!/usr/bin/env bash
set -u -o pipefail

# Usage:
#   bash run_strong_scaling.sh
#
# Optional environment overrides:
#   OMP_LIST="1 2 4 8 16" MPI_LIST="1 2 4 8 16" OUTDIR="results/strong_scaling" bash run_strong_scaling.sh

OMP_LIST="${OMP_LIST:-1 2 4 8 16}"
MPI_LIST="${MPI_LIST:-1 2 4 8 16}"

SERIAL_EXE="${SERIAL_EXE:-./serial_base.exe}"
OMP_EXE="${OMP_EXE:-./omp_base.exe}"
MPI_EXE="${MPI_EXE:-./mpi_base.exe}"

OUTDIR="${OUTDIR:-results/strong_scaling}"
mkdir -p "${OUTDIR}"

RAW_CSV="${OUTDIR}/strong_scaling_raw.csv"
LOGFILE="${OUTDIR}/strong_scaling.log"
FAILFILE="${OUTDIR}/strong_scaling_failed.log"

echo "mode,p,wall_seconds,steps,status" > "${RAW_CSV}"
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
    shift 2

    echo "Running mode=${mode} p=${p}" | tee -a "${LOGFILE}"

    local output
    local rc=0

    if output=$("$@" 2>&1); then
        echo "${output}" >> "${LOGFILE}"

        local wall
        local steps
        wall=$(echo "${output}" | extract_field wall_seconds)
        steps=$(echo "${output}" | extract_field steps)

        if [[ -z "${wall}" ]]; then
            echo "FAILED: missing wall_seconds | mode=${mode} p=${p}" | tee -a "${LOGFILE}" "${FAILFILE}"
            echo "${mode},${p},,${steps:-},missing_wall_seconds" >> "${RAW_CSV}"
        else
            echo "${mode},${p},${wall},${steps:-},ok" >> "${RAW_CSV}"
        fi
    else
        rc=$?
        echo "${output}" >> "${LOGFILE}"
        echo "FAILED: exit_code=${rc} | mode=${mode} p=${p}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},,,failed_exit_${rc}" >> "${RAW_CSV}"
    fi
}

echo "[checkpoint] start serial" | tee -a "${LOGFILE}"
run_and_record "serial" 1 "${SERIAL_EXE}"
echo "[checkpoint] finished serial" | tee -a "${LOGFILE}"

echo "[checkpoint] start omp" | tee -a "${LOGFILE}"
for p in ${OMP_LIST}; do
    run_and_record "omp" "${p}" \
        env OMP_NUM_THREADS="${p}" "${OMP_EXE}"
done
echo "[checkpoint] finished omp" | tee -a "${LOGFILE}"

echo "[checkpoint] start mpi" | tee -a "${LOGFILE}"
for p in ${MPI_LIST}; do
    run_and_record "mpi" "${p}" \
        mpirun -np "${p}" "${MPI_EXE}"
done
echo "[checkpoint] finished mpi" | tee -a "${LOGFILE}"

echo "Done."
echo "Raw data:   ${RAW_CSV}"
echo "Main log:   ${LOGFILE}"
echo "Fail log:   ${FAILFILE}"