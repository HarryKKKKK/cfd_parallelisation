#!/usr/bin/env bash
set -euo pipefail

# =========================
# Config
# =========================
NX=${1:-500}
NY=${2:-197}
TEND=${3:-0.0011741}
REPEAT=${4:-5}

OMP_LIST="${OMP_LIST:-1 2 4 8 16}"
MPI_LIST="${MPI_LIST:-1 2 4 8 16}"

SERIAL_EXE="${SERIAL_EXE:-./serial_scaling.exe}"
OMP_EXE="${OMP_EXE:-./omp_scaling.exe}"
MPI_EXE="${MPI_EXE:-./mpi_scaling.exe}"

OUTDIR="${OUTDIR:-results}"
mkdir -p "${OUTDIR}"

RAW_CSV="${OUTDIR}/strong_scaling_raw.csv"
LOGFILE="${OUTDIR}/strong_scaling.log"

# =========================
# CSV header
# =========================
echo "mode,p,run_id,nx,ny,t_end,wall_seconds" > "${RAW_CSV}"
: > "${LOGFILE}"

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
    shift 3

    echo "Running mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}"

    local output
    output=$("$@")

    echo "${output}" >> "${LOGFILE}"

    local wall
    wall=$(echo "${output}" | extract_field wall_seconds)

    if [[ -z "${wall}" ]]; then
        echo "ERROR: wall_seconds not found in output" | tee -a "${LOGFILE}"
        echo "Output was: ${output}" | tee -a "${LOGFILE}"
        exit 1
    fi

    echo "${mode},${p},${run_id},${NX},${NY},${TEND},${wall}" >> "${RAW_CSV}"
}

# =========================
# Serial
# =========================
for r in $(seq 1 "${REPEAT}"); do
    run_and_record "serial" 1 "${r}" \
        "${SERIAL_EXE}" "${NX}" "${NY}" "${TEND}"
done

# =========================
# OpenMP
# =========================
for p in ${OMP_LIST}; do
    for r in $(seq 1 "${REPEAT}"); do
        run_and_record "omp" "${p}" "${r}" \
            env OMP_NUM_THREADS="${p}" "${OMP_EXE}" "${NX}" "${NY}" "${TEND}"
    done
done

# =========================
# MPI
# =========================
for p in ${MPI_LIST}; do
    for r in $(seq 1 "${REPEAT}"); do
        run_and_record "mpi" "${p}" "${r}" \
            mpirun -np "${p}" "${MPI_EXE}" "${NX}" "${NY}" "${TEND}"
    done
done

echo "Done."
echo "Raw data: ${RAW_CSV}"
echo "Log file: ${LOGFILE}"