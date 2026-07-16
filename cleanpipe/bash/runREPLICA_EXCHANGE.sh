#!/bin/bash

# Build and launch a temperature replica-exchange MD (T-REMD) simulation.
#
# Usage:
#   ./runREPLICA_EXCHANGE.sh NAME GRO TOP NS_PER_REPLICA TMIN TMAX FF ARCH NTOMP NTMPI_PER_REPLICA [NREP]
#
# Required arguments:
#   1  NAME                 Name used in runREPLICA_EXCHANGE_NAME and the scheduler job name
#   2  GRO                  Starting coordinates; must already be suitable for MD (normally minimized/equilibrated)
#   3  TOP                  GROMACS topology
#   4  NS_PER_REPLICA       Length of every replica, in ns
#   5  TMIN                 Lowest temperature, in K
#   6  TMAX                 Highest temperature, in K
#   7  FF                   charmm36 or martini3
#   8  ARCH                 pc, slurm, or rome
#   9  NTOMP                OpenMP threads per MPI rank
#   10 NTMPI_PER_REPLICA    MPI ranks assigned to each replica
#
# Optional argument:
#   11 NREP                 Number of replicas. If omitted or "auto", estimate it from the atom count
#                           using the GROMACS temperature-spacing heuristic.
#
# Example for 32 replicas on one 128-core Rome node:
#   ./runREPLICA_EXCHANGE.sh CHYO500 chyo500.gro topol.top 500 298.15 340 charmm36 rome 4 1 32
#
# Optional environment variables:
#   PCOUPLTYPE=semiisotropic   isotropic or semiisotropic; default: semiisotropic
#   REPLEX_PS=1.0              exchange-attempt interval in ps; default: 1.0
#   GMX_MPI_BIN=gmx_mpi        external-MPI GROMACS command used by the generated job
#   SLURM_PARTITION=calcul     Slurm partition; default preserves the example script
#   SLURM_GPUS=1               GPUs requested by Slurm; set to 0 for CPU-only
#   SLURM_EXCLUDE=node-15      Slurm excluded node; empty value disables exclusion
#
# Important:
#   Replica exchange in GROMACS requires an external-MPI build. A thread-MPI-only `gmx`
#   executable is not sufficient. The generated job therefore expects gmx_mpi by default.

set -euo pipefail

usage() {
    sed -n '3,42p' "$0"
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_number() {
    awk -v value="$1" 'BEGIN { exit !(value + 0 > 0) }'
}

absolute_path() {
    local path=$1
    if command -v realpath >/dev/null 2>&1; then
        realpath -- "$path"
    else
        local directory basename_value
        directory=$(dirname -- "$path")
        basename_value=$(basename -- "$path")
        (cd "$directory" && printf '%s/%s\n' "$PWD" "$basename_value")
    fi
}

if [[ $# -lt 10 || $# -gt 11 ]]; then
    usage
    fail "10 required arguments and at most 1 optional argument are accepted."
fi

NAME=$1
GRO=$2
TOP=$3
NS_PER_REPLICA=$4
TMIN=$5
TMAX=$6
FF=$7
ARCHITECTURE=$8
NTOMP=$9
NTMPI_PER_REPLICA=${10}
NREP_REQUESTED=${11:-auto}

PCOUPLTYPE=${PCOUPLTYPE:-semiisotropic}
REPLEX_PS=${REPLEX_PS:-1.0}
GMX_MPI_BIN_SETTING=${GMX_MPI_BIN:-gmx_mpi}
SLURM_PARTITION_SETTING=${SLURM_PARTITION:-calcul}
SLURM_GPUS_SETTING=${SLURM_GPUS:-1}
SLURM_EXCLUDE_SETTING=${SLURM_EXCLUDE:-node-15}

[[ "$NAME" =~ ^[A-Za-z0-9._-]+$ ]] || fail "NAME may contain only letters, numbers, dot, underscore, and hyphen."
[[ -f "$GRO" ]] || fail "GRO file not found: $GRO"
[[ -f "$TOP" ]] || fail "TOP file not found: $TOP"
[[ "$FF" == "charmm36" || "$FF" == "martini3" ]] || fail "FF must be charmm36 or martini3."
[[ "$ARCHITECTURE" == "pc" || "$ARCHITECTURE" == "slurm" || "$ARCHITECTURE" == "rome" ]] || fail "ARCH must be pc, slurm, or rome."
[[ "$NTOMP" =~ ^[1-9][0-9]*$ ]] || fail "NTOMP must be a positive integer."
[[ "$NTMPI_PER_REPLICA" =~ ^[1-9][0-9]*$ ]] || fail "NTMPI_PER_REPLICA must be a positive integer."
[[ "$NREP_REQUESTED" == "auto" || "$NREP_REQUESTED" =~ ^[0-9]+$ ]] || fail "NREP must be auto or an integer of at least 2."
if [[ "$NREP_REQUESTED" != "auto" ]]; then
    (( NREP_REQUESTED >= 2 )) || fail "NREP must be at least 2."
fi
[[ "$PCOUPLTYPE" == "isotropic" || "$PCOUPLTYPE" == "semiisotropic" ]] || fail "PCOUPLTYPE must be isotropic or semiisotropic."
[[ "$SLURM_GPUS_SETTING" =~ ^[0-9]+$ ]] || fail "SLURM_GPUS must be a non-negative integer."
is_positive_number "$NS_PER_REPLICA" || fail "NS_PER_REPLICA must be positive."
is_positive_number "$TMIN" || fail "TMIN must be positive."
is_positive_number "$TMAX" || fail "TMAX must be positive."
is_positive_number "$REPLEX_PS" || fail "REPLEX_PS must be positive."
awk -v lo="$TMIN" -v hi="$TMAX" 'BEGIN { exit !(hi > lo) }' || fail "TMAX must be greater than TMIN."

GRO_ABS=$(absolute_path "$GRO")
TOP_ABS=$(absolute_path "$TOP")
INPUT_DIR=$(dirname "$TOP_ABS")

NATOMS=$(awk 'NR == 2 { gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0); if ($0 ~ /^[0-9]+$/) print $0 }' "$GRO_ABS")
[[ "$NATOMS" =~ ^[1-9][0-9]*$ ]] || fail "Could not read a valid atom count from line 2 of $GRO_ABS."

if [[ "$FF" == "charmm36" ]]; then
    DT=0.002
    NSTLIST=10
    RLIST=1.2
    RCOULOMB=1.2
    COULOMBTYPE=PME
    EPSILON_R=1
    RVDW=1.2
    VDW_MODIFIER=force-switch
    RVDW_SWITCH=1.0
    DISPCORR=no
    CONSTRAINTS=h-bonds
    TAU_P=5
    COMPRESSIBILITY=4.5e-5
else
    DT=0.02
    NSTLIST=20
    RLIST=1.1
    RCOULOMB=1.1
    COULOMBTYPE=reaction-field
    EPSILON_R=15
    RVDW=1.1
    VDW_MODIFIER=potential-shift-verlet
    RVDW_SWITCH=0
    DISPCORR=no
    CONSTRAINTS=none
    TAU_P=12
    COMPRESSIBILITY=3e-4
fi

STEPS=$(awk -v ns="$NS_PER_REPLICA" -v dt="$DT" 'BEGIN { printf "%.0f", ns * 1000.0 / dt }')
REPLEX_STEPS=$(awk -v ps="$REPLEX_PS" -v dt="$DT" 'BEGIN { value = int(ps / dt + 0.5); if (value < 1) value = 1; print value }')
NSTENERGY=$(awk -v dt="$DT" 'BEGIN { value = int(10.0 / dt + 0.5); if (value < 1) value = 1; print value }')
NSTXTC=$(awk -v dt="$DT" 'BEGIN { value = int(100.0 / dt + 0.5); if (value < 1) value = 1; print value }')

# GROMACS manual estimate: epsilon ~= 1/sqrt(Natoms), then use a geometric ladder.
AUTO_NREP=$(awk -v n="$NATOMS" -v lo="$TMIN" -v hi="$TMAX" 'BEGIN {
    epsilon = 1.0 / sqrt(n)
    intervals = log(hi / lo) / log(1.0 + epsilon)
    replicas = int(intervals)
    if (replicas < intervals) replicas++
    replicas++
    if (replicas < 4) replicas = 4
    print replicas
}')

if [[ "$NREP_REQUESTED" == "auto" ]]; then
    NREP=$AUTO_NREP
else
    NREP=$NREP_REQUESTED
fi

TOTAL_MPI=$((NREP * NTMPI_PER_REPLICA))
TOTAL_CORES=$((TOTAL_MPI * NTOMP))

if (( NREP < AUTO_NREP )); then
    echo "WARNING: You requested $NREP replicas, while the atom-count heuristic suggests about $AUTO_NREP."
    echo "         The run can work, but neighboring energy distributions may overlap poorly."
fi

if [[ "$PCOUPLTYPE" == "semiisotropic" ]]; then
    REF_P="1.0 1.0"
    COMPRESSIBILITY_LINE="$COMPRESSIBILITY $COMPRESSIBILITY"
else
    REF_P="1.0"
    COMPRESSIBILITY_LINE="$COMPRESSIBILITY"
fi

RUN_DIR="runREPLICA_EXCHANGE_${NAME}"
RUN_DIR_ABS=$(absolute_path "$(pwd)")/${RUN_DIR}
SCRIPT_NAME="script.${NAME}.replica_exchange.sh"

CONFIG_DESCRIPTION=$(cat <<EOF
NAME=$NAME
GRO_ABS=$GRO_ABS
TOP_ABS=$TOP_ABS
NS_PER_REPLICA=$NS_PER_REPLICA
TMIN=$TMIN
TMAX=$TMAX
FF=$FF
ARCHITECTURE=$ARCHITECTURE
NTOMP=$NTOMP
NTMPI_PER_REPLICA=$NTMPI_PER_REPLICA
NREP=$NREP
PCOUPLTYPE=$PCOUPLTYPE
REPLEX_PS=$REPLEX_PS
GMX_MPI_BIN=$GMX_MPI_BIN_SETTING
SLURM_PARTITION=$SLURM_PARTITION_SETTING
SLURM_GPUS=$SLURM_GPUS_SETTING
SLURM_EXCLUDE=$SLURM_EXCLUDE_SETTING
EOF
)

launch_existing_or_new_job() {
    cd "$RUN_DIR_ABS"
    case "$ARCHITECTURE" in
        pc)
            "./$SCRIPT_NAME"
            ;;
        slurm)
            sbatch "$SCRIPT_NAME" && echo "Job was submitted."
            ;;
        rome)
            ccc_msub "$SCRIPT_NAME" && echo "Job was submitted."
            ;;
    esac
}

if [[ -d "$RUN_DIR_ABS" ]]; then
    [[ -f "$RUN_DIR_ABS/00_SETUP/config_description.txt" ]] || fail "$RUN_DIR_ABS exists but is not a recognized run created by this script."
    if [[ "$(cat "$RUN_DIR_ABS/00_SETUP/config_description.txt")" != "$CONFIG_DESCRIPTION" ]]; then
        fail "$RUN_DIR_ABS already exists with different parameters. Use another NAME or remove the old folder deliberately."
    fi
    [[ -x "$RUN_DIR_ABS/$SCRIPT_NAME" ]] || fail "Existing generated job script is missing or not executable: $RUN_DIR_ABS/$SCRIPT_NAME"
    echo "Existing matching run found. Relaunching its saved job script."
    launch_existing_or_new_job
    exit 0
fi

mkdir -p "$RUN_DIR_ABS/00_SETUP"
printf '%s\n' "$CONFIG_DESCRIPTION" > "$RUN_DIR_ABS/00_SETUP/config_description.txt"

# Generate a geometric temperature ladder and ordered replica-directory list.
awk -v n="$NREP" -v lo="$TMIN" -v hi="$TMAX" 'BEGIN {
    for (i = 0; i < n; i++) {
        if (n == 1) temperature = lo
        else temperature = lo * exp(log(hi / lo) * i / (n - 1))
        printf "%03d %.6f replica_%03d\n", i, temperature, i
    }
}' > "$RUN_DIR_ABS/00_SETUP/temperature_ladder.dat"

awk '{ print $3 }' "$RUN_DIR_ABS/00_SETUP/temperature_ladder.dat" > "$RUN_DIR_ABS/00_SETUP/replica_directories.txt"

# Save configuration as shell-safe assignments for the generated job script.
{
    printf 'NAME=%q\n' "$NAME"
    printf 'GRO_ABS=%q\n' "$GRO_ABS"
    printf 'TOP_ABS=%q\n' "$TOP_ABS"
    printf 'INPUT_DIR=%q\n' "$INPUT_DIR"
    printf 'RUN_DIR_ABS=%q\n' "$RUN_DIR_ABS"
    printf 'SCRIPT_NAME=%q\n' "$SCRIPT_NAME"
    printf 'FF=%q\n' "$FF"
    printf 'ARCHITECTURE=%q\n' "$ARCHITECTURE"
    printf 'NATOMS=%q\n' "$NATOMS"
    printf 'NS_PER_REPLICA=%q\n' "$NS_PER_REPLICA"
    printf 'TMIN=%q\n' "$TMIN"
    printf 'TMAX=%q\n' "$TMAX"
    printf 'NREP=%q\n' "$NREP"
    printf 'AUTO_NREP=%q\n' "$AUTO_NREP"
    printf 'NTOMP=%q\n' "$NTOMP"
    printf 'NTMPI_PER_REPLICA=%q\n' "$NTMPI_PER_REPLICA"
    printf 'TOTAL_MPI=%q\n' "$TOTAL_MPI"
    printf 'TOTAL_CORES=%q\n' "$TOTAL_CORES"
    printf 'DT=%q\n' "$DT"
    printf 'STEPS=%q\n' "$STEPS"
    printf 'REPLEX_STEPS=%q\n' "$REPLEX_STEPS"
    printf 'GMX_MPI_BIN_SETTING=%q\n' "$GMX_MPI_BIN_SETTING"
} > "$RUN_DIR_ABS/00_SETUP/config.sh"

cat > "$RUN_DIR_ABS/00_SETUP/run_summary.txt" <<EOF
Replica-exchange run: $NAME
Atoms per replica: $NATOMS
Force field preset: $FF
Length per replica: $NS_PER_REPLICA ns
Temperature range: $TMIN to $TMAX K
Replicas used: $NREP
Atom-count heuristic: $AUTO_NREP replicas
Temperature ladder: geometric
Exchange attempts: every $REPLEX_PS ps = $REPLEX_STEPS steps
MPI ranks per replica: $NTMPI_PER_REPLICA
OpenMP threads per MPI rank: $NTOMP
Total MPI ranks: $TOTAL_MPI
Total requested CPU cores: $TOTAL_CORES
Pressure coupling type: $PCOUPLTYPE
EOF

while read -r replica_index replica_temperature replica_directory; do
    replica_path="$RUN_DIR_ABS/$replica_directory"
    mkdir -p "$replica_path"
    printf '%s\n' "$replica_temperature" > "$replica_path/temperature.K.txt"
    GEN_SEED=$((100000 + 1009 * 10#$replica_index))

    {
        cat <<EOF
; Temperature replica-exchange MD
; Replica index: $replica_index
; Replica temperature: $replica_temperature K
; Generated by runREPLICA_EXCHANGE.sh

; Run control
integrator               = md
tinit                    = 0
dt                       = $DT
nsteps                   = $STEPS
continuation             = no

; Center-of-mass motion
comm-mode                = Linear
nstcomm                  = 100
comm-grps                = System

; Output control
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = $NSTENERGY
nstenergy                = $NSTENERGY
nstcalcenergy            = 100
nstxout-compressed       = $NSTXTC
compressed-x-precision   = 1000

; Periodic box and neighbor searching
pbc                      = xyz
periodic-molecules       = no
cutoff-scheme            = Verlet
nstlist                  = $NSTLIST
verlet-buffer-tolerance  = 0.005
rlist                    = $RLIST

; Electrostatics
coulombtype              = $COULOMBTYPE
rcoulomb                 = $RCOULOMB
epsilon_r                = $EPSILON_R
epsilon_rf               = 0
EOF
        if [[ "$FF" == "charmm36" ]]; then
            cat <<EOF
pme_order                = 4
fourierspacing           = 0.12
ewald_rtol               = 1e-5
EOF
        fi
        cat <<EOF

; Van der Waals
vdwtype                  = cutoff
vdw-modifier             = $VDW_MODIFIER
rvdw-switch              = $RVDW_SWITCH
rvdw                     = $RVDW
DispCorr                 = $DISPCORR

; Initial velocities for this replica
continuation             = no
gen-vel                  = yes
gen-temp                 = $replica_temperature
gen-seed                 = $GEN_SEED

; Temperature coupling
Tcoupl                   = V-rescale
tc-grps                  = System
ref-t                    = $replica_temperature
tau-t                    = 1.0

; Pressure coupling
Pcoupl                   = C-rescale
Pcoupltype               = $PCOUPLTYPE
ref-p                    = $REF_P
tau-p                    = $TAU_P
compressibility          = $COMPRESSIBILITY_LINE
refcoord-scaling         = com

; Constraints
constraints              = $CONSTRAINTS
constraint-algorithm     = lincs
lincs-order              = 4
lincs-iter               = 1
EOF
    } > "$replica_path/remd.mdp"
done < "$RUN_DIR_ABS/00_SETUP/temperature_ladder.dat"

# Build scheduler-specific header. The executable body is common to all architectures.
{
    echo '#!/bin/bash'
    echo
    if [[ "$ARCHITECTURE" == "rome" ]]; then
        cat <<EOF
#MSUB   -r ${NAME}.remd
#MSUB   -n ${TOTAL_MPI}
#MSUB   -c ${NTOMP}
#MSUB   -W yes
#MSUB   -o out.scheduler.%I.${NAME}
#MSUB   -e err.scheduler.%I.${NAME}
#MSUB   -q rome
#MSUB   -A gen13458
#MSUB   -m scratch,work,store
#MSUB   -Q normal
#MSUB   -T 86400
EOF
    elif [[ "$ARCHITECTURE" == "slurm" ]]; then
        cat <<EOF
#SBATCH --partition=${SLURM_PARTITION_SETTING}
#SBATCH --ntasks=${TOTAL_MPI}
#SBATCH --cpus-per-task=${NTOMP}
#SBATCH --job-name=${NAME}.remd
#SBATCH --output=outanderr.slurm.%j.${NAME}
EOF
        if (( SLURM_GPUS_SETTING > 0 )); then
            echo "#SBATCH --gres=gpu:${SLURM_GPUS_SETTING}"
        fi
        if [[ -n "$SLURM_EXCLUDE_SETTING" ]]; then
            echo "#SBATCH --exclude=${SLURM_EXCLUDE_SETTING}"
        fi
    fi
    cat <<'JOB_BODY'

set -euo pipefail

RUN_DIR=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=/dev/null
source "$RUN_DIR/00_SETUP/config.sh"
cd "$RUN_DIR"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

if command -v module >/dev/null 2>&1; then
    if [[ "$ARCHITECTURE" == "rome" ]]; then
        module purge
        module load gnu/11
        module load nvhpc/24.3
        module load mpi/openmpi/4
        module load gromacs/2025.0
    elif [[ "$ARCHITECTURE" == "slurm" ]]; then
        module purge
        module load cuda/11.8
        module load gromacs/2024.5
    elif [[ "$ARCHITECTURE" == "pc" ]]; then
        # Preserve the module choices from the example script when a module system exists.
        module purge
        module load cuda/11.8
        module load gromacs/2024.5
    fi
fi

GMX_MPI_BIN=${GMX_MPI_BIN:-$GMX_MPI_BIN_SETTING}
command -v "$GMX_MPI_BIN" >/dev/null 2>&1 || fail "External-MPI GROMACS executable not found: $GMX_MPI_BIN. Set GMX_MPI_BIN if it has another name."

export OMP_NUM_THREADS=$NTOMP
export OMP_DYNAMIC=FALSE

if [[ "$ARCHITECTURE" == "rome" ]]; then
    export GMX_DISABLE_GPU_DETECTION=1
    export I_MPI_PIN_CELL=core
    export I_MPI_PIN_DOMAIN=auto
    MDRUN_LAUNCHER=(ccc_mprun "$GMX_MPI_BIN")
elif [[ "$ARCHITECTURE" == "slurm" ]]; then
    command -v srun >/dev/null 2>&1 || fail "srun was not found."
    MDRUN_LAUNCHER=(srun --ntasks="$TOTAL_MPI" --cpus-per-task="$NTOMP" "$GMX_MPI_BIN")
else
    command -v mpirun >/dev/null 2>&1 || fail "mpirun was not found; replica exchange requires external MPI."
    MDRUN_LAUNCHER=(mpirun -np "$TOTAL_MPI" "$GMX_MPI_BIN")
fi

mapfile -t REPLICA_DIRS < "$RUN_DIR/00_SETUP/replica_directories.txt"
(( ${#REPLICA_DIRS[@]} == NREP )) || fail "Replica-directory count does not match NREP."
(( TOTAL_MPI % NREP == 0 )) || fail "TOTAL_MPI must be a multiple of NREP."
(( TOTAL_MPI / NREP == NTMPI_PER_REPLICA )) || fail "MPI-ranks-per-replica consistency check failed."

all_replicas_complete() {
    local replica_directory checkpoint checkpoint_step
    for replica_directory in "${REPLICA_DIRS[@]}"; do
        checkpoint="$RUN_DIR/$replica_directory/remd.cpt"
        [[ -f "$checkpoint" ]] || return 1
        checkpoint_step=$(
            "$GMX_MPI_BIN" dump -cp "$checkpoint" 2>/dev/null |
                awk '/^[[:space:]]*step[[:space:]]*=/ { print $3; exit }'
        )
        [[ "$checkpoint_step" =~ ^[0-9]+$ ]] || return 1
        (( checkpoint_step >= STEPS )) || return 1
    done
    return 0
}

if [[ -f "$RUN_DIR/REMD_COMPLETE" ]] || all_replicas_complete; then
    touch "$RUN_DIR/REMD_COMPLETE"
    echo "Replica-exchange simulation is already complete."
    exit 0
fi

# Build each .tpr only once. Run grompp from the topology directory so relative #include paths work.
while read -r replica_index replica_temperature replica_directory; do
    replica_path="$RUN_DIR/$replica_directory"
    if [[ ! -f "$replica_path/remd.tpr" ]]; then
        echo "#################################################################"
        echo "Creating TPR for $replica_directory at $replica_temperature K"
        echo "#################################################################"
        (
            cd "$INPUT_DIR"
            "$GMX_MPI_BIN" grompp \
                -f "$replica_path/remd.mdp" \
                -c "$GRO_ABS" \
                -r "$GRO_ABS" \
                -p "$TOP_ABS" \
                -o "$replica_path/remd.tpr" \
                -po "$replica_path/mdout.mdp"
        ) 2>&1 | tee "$replica_path/outanderr.grompp"
    fi
    [[ -s "$replica_path/remd.tpr" ]] || fail "Missing or empty TPR: $replica_path/remd.tpr"
done < "$RUN_DIR/00_SETUP/temperature_ladder.dat"

# Queue exactly one dependent continuation on TGCC Rome. It will start only if this job ends successfully.
if [[ "$ARCHITECTURE" == "rome" ]]; then
    [[ -n "${BRIDGE_MSUB_JOBID:-}" ]] || fail "BRIDGE_MSUB_JOBID is not defined by the Rome scheduler."
    echo "Submitting the dependent continuation job before starting mdrun."
    ccc_msub -E "--dependency=afterok:${BRIDGE_MSUB_JOBID}" "$RUN_DIR/$SCRIPT_NAME"
fi

MDRUN_ARGUMENTS=(
    mdrun
    -v
    -multidir "${REPLICA_DIRS[@]}"
    -deffnm remd
    -replex "$REPLEX_STEPS"
    -reseed -1
    -cpi remd.cpt
    -append
    -cpt 30
    -ntomp "$NTOMP"
    -pin on
)

if [[ "$ARCHITECTURE" == "rome" ]]; then
    MDRUN_ARGUMENTS+=( -maxh 23 )
fi

echo "#################################################################"
echo "Starting T-REMD: $NREP replicas, $TOTAL_MPI total MPI ranks, $NTOMP OpenMP threads/rank"
echo "#################################################################"
printf 'Command:'
printf ' %q' "${MDRUN_LAUNCHER[@]}" "${MDRUN_ARGUMENTS[@]}"
printf '\n'

"${MDRUN_LAUNCHER[@]}" "${MDRUN_ARGUMENTS[@]}" 2>&1 | tee -a "$RUN_DIR/outanderr.mdrun"

if all_replicas_complete; then
    touch "$RUN_DIR/REMD_COMPLETE"
    echo "Replica-exchange simulation completed all requested steps."
else
    echo "The run stopped cleanly before all requested steps were completed."
    echo "Relaunch this same saved job script; GROMACS will continue from each remd.cpt."
fi
JOB_BODY
} > "$RUN_DIR_ABS/$SCRIPT_NAME"

chmod +x "$RUN_DIR_ABS/$SCRIPT_NAME"

echo
cat "$RUN_DIR_ABS/00_SETUP/run_summary.txt"
echo
echo "Folder structure created in: $RUN_DIR_ABS"
echo "Saved job script: $RUN_DIR_ABS/$SCRIPT_NAME"
echo

launch_existing_or_new_job
