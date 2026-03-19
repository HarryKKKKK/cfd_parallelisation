#!/bin/bash
#SBATCH --job-name=weak_scaling
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=hk597
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=weak_scaling_%j.out

set -u -o pipefail

# Optional environment overrides
OMP_LIST="${OMP_LIST:-1 2 4 8 16}"
MPI_LIST="${MPI_LIST:-1 2 4 8 16}"
REPEATS="${REPEATS:-2}"

# Make sure these point to your weak scaling executables
SERIAL_EXE="${SERIAL_EXE:-./serial_scaling.exe}"
OMP_EXE="${OMP_EXE:-./omp_scaling.exe}"
MPI_EXE="${MPI_EXE:-./mpi_scaling.exe}"

OUTDIR="${OUTDIR:-results/weak_scaling}"
mkdir -p "${OUTDIR}"

RAW_CSV="${OUTDIR}/weak_scaling_raw.csv"
LOGFILE="${OUTDIR}/weak_scaling.log"
FAILFILE="${OUTDIR}/weak_scaling_failed.log"

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

# Added weak scaling specific columns: nx, copies_y, ny
echo "mode,p,run,nx,copies_y,ny,wall_seconds,steps,status" > "${RAW_CSV}"
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

    local wall steps nx ny copies_y
    wall=$(echo "${output}" | extract_field wall_seconds)
    steps=$(echo "${output}" | extract_field steps)
    nx=$(echo "${output}" | extract_field nx)
    ny=$(echo "${output}" | extract_field ny)
    copies_y=$(echo "${output}" | extract_field copies_y)

    # Fallbacks in case execution failed before printing INIT
    nx=${nx:-}
    ny=${ny:-}
    copies_y=${copies_y:-$p}
    steps=${steps:-}

    if [[ "${rc}" -ne 0 ]]; then
        echo "FAILED: exit_code=${rc} | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},${nx},${copies_y},${ny},,${steps},failed_exit_${rc}" >> "${RAW_CSV}"
        return
    fi

    if [[ -z "${wall}" ]]; then
        echo "FAILED: missing wall_seconds | mode=${mode} p=${p} run=${run_id}" | tee -a "${LOGFILE}" "${FAILFILE}"
        echo "${mode},${p},${run_id},${nx},${copies_y},${ny},,${steps},missing_wall_seconds" >> "${RAW_CSV}"
        return
    fi

    echo "${mode},${p},${run_id},${nx},${copies_y},${ny},${wall},${steps},ok" >> "${RAW_CSV}"
}

summarise_mode() {
    local mode="$1"
    # Note: Column indices shifted because of nx, copies_y, ny
    # $1=mode, $2=p, $7=wall_seconds, $9=status
    awk -F, -v target="${mode}" '
        NR==1 {next}
        $1==target && $9=="ok" {
            cnt[$2]++
            sum[$2]+=$7
            if (!($2 in min) || $7<min[$2]) min[$2]=$7
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
    # Pass 1 as copies_y for serial
    run_and_record "serial" 1 "${r}" "${SERIAL_EXE}" 1
done
echo "[checkpoint] finished serial" | tee -a "${LOGFILE}"

echo "[checkpoint] start omp" | tee -a "${LOGFILE}"
for p in ${OMP_LIST}; do
    for ((r=1; r<=REPEATS; r++)); do
        # Pass $p as copies_y for OpenMP
        run_and_record "omp" "${p}" "${r}" \
            env OMP_NUM_THREADS="${p}" "${OMP_EXE}" "${p}"
    done
done
echo "[checkpoint] finished omp" | tee -a "${LOGFILE}"

echo "[checkpoint] start mpi" | tee -a "${LOGFILE}"
for p in ${MPI_LIST}; do
    for ((r=1; r<=REPEATS; r++)); do
        # Pass $p as copies_y for MPI
        run_and_record "mpi" "${p}" "${r}" \
            mpiexec -n "${p}" "${MPI_EXE}" "${p}"
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