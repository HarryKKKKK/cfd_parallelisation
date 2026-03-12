#!/bin/bash
#SBATCH --job-name=strong_scaling
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=hk597
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=strong_scaling_%j.out

set -u -o pipefail

# Optional environment overrides
OMP_LIST="${OMP_LIST:-1 2 4 8 16}"
MPI_LIST="${MPI_LIST:-1 2 4 8 16}"
REPEATS="${REPEATS:-2}"

SERIAL_EXE="${SERIAL_EXE:-./serial_base.exe}"
OMP_EXE="${OMP_EXE:-./omp_base.exe}"
MPI_EXE="${MPI_EXE:-./mpi_base.exe}"

OUTDIR="${OUTDIR:-results/strong_scaling/hllc}"
mkdir -p "${OUTDIR}"

RAW_CSV="${OUTDIR}/strong_scaling_raw.csv"
LOGFILE="${OUTDIR}/strong_scaling.log"
FAILFILE="${OUTDIR}/strong_scaling_failed.log"

echo "===== Environment ====="
date
hostname
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-}"
echo "PWD=$PWD"
echo "PATH=$PATH"
which srun || true
which mpirun || true
which mpiexec || true
echo "======================="

echo "mode,p,run,wall_seconds,steps,status" > "${RAW_CSV}"
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

    echo "Running mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}"

    local output
    local rc=0

    output=$("$@" 2>&1) || rc=$?
    echo "${output}" >> "${LOGFILE}"

    local wall
    local steps
    wall=$(echo "${output}" | extract_field wall_seconds)
    steps=$(echo "${output}" | extract_field steps)

    if [[ "${rc}" -ne 0 ]]; then
        echo "FAILED: exit_code=${rc} | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},,${steps:-},failed_exit_${rc}" >> "${RAW_CSV}"
        return
    fi

    if [[ -z "${wall}" ]]; then
        echo "FAILED: missing wall_seconds | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},,${steps:-},missing_wall_seconds" >> "${RAW_CSV}"
        return
    fi

    echo "${mode},${p},${run_id},${wall},${steps:-},ok" >> "${RAW_CSV}"
}

summarise_mode() {
    local mode="$1"
    awk -F, -v target="${mode}" '
        NR==1 {next}
        $1==target && $6=="ok" {
            cnt[$2]++
            sum[$2]+=$4
            if (!($2 in min) || $4<min[$2]) min[$2]=$4
        }
        END {
            for (p in cnt) {
                printf "%s,p=%s,avg=%.6f,min=%.6f,n=%d\n", target, p, sum[p]/cnt[p], min[p], cnt[p]
            }
        }
    ' "${RAW_CSV}" | sort -t= -k2,2n
}

echo "[checkpoint] start serial" | tee -a "${LOGFILE}"
for ((r=1; r<=REPEATS; r++)); do
    run_and_record "serial" 1 "${r}" "${SERIAL_EXE}"
done
echo "[checkpoint] finished serial" | tee -a "${LOGFILE}"

echo "[checkpoint] start omp" | tee -a "${LOGFILE}"
for p in ${OMP_LIST}; do
    for ((r=1; r<=REPEATS; r++)); do
        run_and_record "omp" "${p}" "${r}" \
            env OMP_NUM_THREADS="${p}" "${OMP_EXE}"
    done
done
echo "[checkpoint] finished omp" | tee -a "${LOGFILE}"

echo "[checkpoint] start mpi" | tee -a "${LOGFILE}"
for p in ${MPI_LIST}; do
    for ((r=1; r<=REPEATS; r++)); do
        run_and_record "mpi" "${p}" "${r}" \
            mpiexec -n "${p}" "${MPI_EXE}"
    done
done
echo "[checkpoint] finished mpi" | tee -a "${LOGFILE}"

echo "" | tee -a "${LOGFILE}"
echo "===== SUMMARY =====" | tee -a "${LOGFILE}"
summarise_mode "serial" | tee -a "${LOGFILE}"
summarise_mode "omp" | tee -a "${LOGFILE}"
summarise_mode "mpi" | tee -a "${LOGFILE}"
echo "===================" | tee -a "${LOGFILE}"

echo "Done."
echo "Raw data:   ${RAW_CSV}"
echo "Main log:   ${LOGFILE}"
echo "Fail log:   ${FAILFILE}"