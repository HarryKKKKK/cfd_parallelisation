#!/bin/bash
#SBATCH --job-name=res_scaling
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=hk597
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6:00:00 
#SBATCH --output=res_scaling_%j.out

set -u -o pipefail

TEST_NX=${1:-1500}
TEST_NY=${2:-591}

OMP_LIST="${OMP_LIST:-1 2 4 8 16}"
MPI_LIST="${MPI_LIST:-1 2 4 8 16}"
REPEATS="${REPEATS:-1}"

SERIAL_EXE="${SERIAL_EXE:-./serial_base.exe}"
OMP_EXE="${OMP_EXE:-./omp_base.exe}"
MPI_EXE="${MPI_EXE:-./mpi_base.exe}"

OUTDIR="${OUTDIR:-results/resolution}"
mkdir -p "${OUTDIR}"

RAW_CSV="${OUTDIR}/res_${TEST_NX}_raw.csv"
LOGFILE="${OUTDIR}/res_.log"
FAILFILE="${OUTDIR}/res_failed_.log"

echo "===== Environment ====="
date
hostname
echo "Target Resolution: nx=${TEST_NX}, ny=${TEST_NY}"
echo "Output CSV: ${RAW_CSV}"
echo "======================="

# 初始化 CSV 表头
echo "mode,p,run,nx,ny,wall_seconds,steps,status" > "${RAW_CSV}"
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
    shift 3

    echo "Running mode=${mode} p=${p} run=${run_id} (nx=${TEST_NX}, ny=${TEST_NY})" | tee -a "${LOGFILE}"

    local output
    local rc=0

    output=$("$@" "${TEST_NX}" "${TEST_NY}" 2>&1) || rc=$?
    echo "${output}" >> "${LOGFILE}"

    local wall steps nx ny
    wall=$(echo "${output}" | extract_field wall_seconds)
    steps=$(echo "${output}" | extract_field steps)
    nx=$(echo "${output}" | extract_field nx)
    ny=$(echo "${output}" | extract_field ny)

    if [[ "${rc}" -ne 0 ]]; then
        echo "FAILED: exit_code=${rc} | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},${TEST_NX},${TEST_NY},,${steps:-},failed_exit_${rc}" >> "${RAW_CSV}"
        return
    fi

    if [[ -z "${wall}" ]]; then
        echo "FAILED: missing wall_seconds" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},${TEST_NX},${TEST_NY},,${steps:-},missing_wall_seconds" >> "${RAW_CSV}"
        return
    fi

    echo "${mode},${p},${run_id},${nx},${ny},${wall},${steps},ok" >> "${RAW_CSV}"
}

echo "[checkpoint] start serial" | tee -a "${LOGFILE}"
for ((r=1; r<=REPEATS; r++)); do
    run_and_record "serial" 1 "${r}" "${SERIAL_EXE}"
done

echo "[checkpoint] start omp" | tee -a "${LOGFILE}"
for p in ${OMP_LIST}; do
    for ((r=1; r<=REPEATS; r++)); do
        run_and_record "omp" "${p}" "${r}" env OMP_NUM_THREADS="${p}" "${OMP_EXE}"
    done
done

echo "[checkpoint] start mpi" | tee -a "${LOGFILE}"
for p in ${MPI_LIST}; do
    for ((r=1; r<=REPEATS; r++)); do
        run_and_record "mpi" "${p}" "${r}" mpiexec -n "${p}" "${MPI_EXE}"
    done
done

echo "Done. Results saved to ${RAW_CSV}"