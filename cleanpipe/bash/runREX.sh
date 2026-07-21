#!/bin/bash

# usage example:  ./runREX.sh bulk bulk_at_310.gro bulk.top 1 313 333 charmm36 rome 7 18
#                 ./runREX.sh bulk bulk_at_310.gro bulk.top 1 313 333 charmm36 genoa 8 18
#                 ./runREX.sh bulk bulk_at_310.gro bulk.top 1 313 333 charmm36 MI300 8 18

# ARGUMENTS:
# 1-name that goes on the runREX_<...> 
# 2-gro
# 3-top
# 4-nanoseconds of production
# 5-temperature min
# 6 temperature max
# 7-forcefield to be used in mdp construction (must be "charmm36" or "martini3")
# 8-architecture (pc, slurm, rome, genoa, MI300) # rome at tgcc. genoa and MI300 at adastra
# 9-ntOMP
#10-ntMPI





if [ $# -lt 10 ]; then
    echo "10 arguments needed : name filename.gro filename.top numberOfNanoseconds temperatureMin temperatureMax forceFieldName architecture ntOMP ntMPI"
    exit 1
fi




######## obtain the arguments defined by the user when he called the function ########
NAME=$1 #name of the thing to be modeled
echo " "
echo "Name of thing to be modeled: ${NAME}"
echo " "

GRO=$2 #name of gro
TOP=$3 #name of top
echo "GRO and TOP file names: ${GRO} ${TOP} (the main inputs)"
echo " "

if [[ ! -f "$GRO" ]]; then
    echo "Error: $GRO not found!"
    exit 1
fi

[[ -r "$GRO" ]] || { echo "ERROR: Cannot read GRO file: $GRO" >&2; exit 1; }

N_ATOMS=$(awk 'NR == 2 { gsub(/[[:space:]\r]/, ""); print; exit }' "$GRO")
[[ "$N_ATOMS" =~ ^[1-9][0-9]*$ ]] ||
    { echo "ERROR: Invalid atom count on line 2 of $GRO: '$N_ATOMS'" >&2; exit 1; }


if [[ ! -f "$TOP" ]]; then
    echo "Error: $TOP not found!"
    exit 1
fi


T_MIN=$5
echo "Temperature (min): ${T_MIN}"
echo " "

T_MAX=$6
echo "Temperature (max): ${T_MAX}"
echo " "

FF=$7
echo "Force field (this will be used to setup the correct mdp parameters. the only possible values are charmm36 or martini3: ${FF}"

if [[ $FF == "charmm36" || $FF == "martini3" ]]; then
    echo "FF is valid"
else
    echo "Error: the force field must be either 'charmm36' or 'martini3'"
    exit 1
fi
echo " "


PRODUCTION_DURATION=$4 #how many nanoseconds
if [[ $FF == "charmm36" ]]; then
  STEPS=$(echo "scale=0; ($PRODUCTION_DURATION / 0.002) * 1000" | bc)
fi

if [[ $FF == "martini3" ]]; then
  STEPS=$(echo "scale=0; ($PRODUCTION_DURATION / 0.02) * 1000" | bc)
fi


echo "Production Duration: ${PRODUCTION_DURATION} ns (implemented in the mdp by setting: ${STEPS} steps x dt )"
echo " "


ARCHITECTURE=$8
echo "Architecture:${ARCHITECTURE}"

if [[ $ARCHITECTURE == "pc" || $ARCHITECTURE == "slurm" || $ARCHITECTURE == "rome" || $ARCHITECTURE == "genoa" || $ARCHITECTURE == "MI300" ]]; then
    echo "Architecture is valid"
else
    echo "Error: architecture must be 'pc' 'slurm' 'rome' 'genoa' 'MI300' "
    exit 1
fi
echo " "


NTOMP=$9
if [[ "$NTOMP" =~ ^[0-9]+$ ]]; then
    echo "NTOMP OK (is an integer)"
else
    echo "NTOMP is NOT an integer"
    exit 1
fi



NTMPI=${10}
if [[ "$NTMPI" =~ ^[0-9]+$ ]]; then
    echo "NTMPI OK (is an integer)"
else
    echo "NTMPI is NOT an integer"
    exit 1
fi
echo " "


########### check top file for [ distance_restraints ] or [ dihedral_restraints ]

# Check if .top contains [ distance_restraints ]
if grep -q "\[ distance_restraints \]" "$TOP"; then
    # If the text is found
    distance_restraints_option="simple"
else
    # If the text is not found
    distance_restraints_option="no"
fi

# Check if .top contains [ dihedral_restraints ] 
if grep -q "\[ dihedral_restraints \]" "$TOP"; then
    # If the text is found
    dihedral_restraints_option="yes"
else
    # If the text is not found
    dihedral_restraints_option="no"
fi











######################## create folder structure ########################

REPLICA_ROOT="runREX_${NAME}"

awk -v min="$T_MIN" -v max="$T_MAX" '
    BEGIN {
        number = "^[0-9]+([.][0-9]+)?$"
        exit !(min ~ number && max ~ number && min > 0 && max > min)
    }
' || { echo "ERROR: T_MIN and T_MAX must satisfy 0 < T_MIN < T_MAX." >&2; exit 1; }

# Calculate the GROMACS estimate and generate a geometric ladder that
# includes T_MIN and T_MAX exactly.
mapfile -t REMD_DATA < <(
    awk -v tmin="$T_MIN" -v tmax="$T_MAX" -v n="$N_ATOMS" '
        function ceiling(x) {
            return (x == int(x)) ? int(x) : int(x) + 1
        }

        BEGIN {
            epsilon = 1 / sqrt(n)

            n_replicas = ceiling(
                log(tmax / tmin) / log(1 + epsilon)
            ) + 1

            ratio = exp(
                log(tmax / tmin) / (n_replicas - 1)
            )

            printf "%.10g %d %.10g\n", epsilon, n_replicas, ratio

            for (i = 0; i < n_replicas; i++)
                printf "%.3f\n", tmin * ratio^i
        }
    '
)

read -r EPSILON N_REPLICAS TEMP_RATIO <<< "${REMD_DATA[0]}"
TEMPERATURES=("${REMD_DATA[@]:1}")

mkdir -p "$REPLICA_ROOT"
printf "%s\n" "${TEMPERATURES[@]}" > "$REPLICA_ROOT/temperatures.dat"

REPLICA_DIRS=()

for i in "${!TEMPERATURES[@]}"; do
    TEMPERATURE="${TEMPERATURES[$i]}"
    printf -v DIR "%s/replica_%03d_T%sK" "$REPLICA_ROOT" "$i" "$TEMPERATURE"

    mkdir -p "$DIR"
    REPLICA_DIRS+=("$DIR")
done
	
echo "folder structure created"



######################## generate mdp files #############################
#this is a funcion so I can read clearly some mdp file options
options() {
  declare -A map
  for pair in "$@"; do
    key="${pair%%=*}"
    val="${pair#*=}"
    map["$key"]="$val"
  done
  echo "${map[$FF]}"
}




# Definition of the name of the mdp files of the current lambda


file_prod_mdp="prod.mdp"
REPLICA_DIRS=()

for i in "${!TEMPERATURES[@]}"; do
    current_t="${TEMPERATURES[$i]}"

    printf -v DIR "%s/replica_%03d_T%sK" \
        "$REPLICA_ROOT" "$i" "$current_t"

    mkdir -p "$DIR"
    REPLICA_DIRS+=("$DIR")

    echo "creating ${DIR}/${file_prod_mdp}"

    cat > "${DIR}/${file_prod_mdp}" <<EOT


; Run control
integrator               = md 
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)
nsteps                   = ${STEPS}
nstcomm                  = 100

; Output control
nstxout                  = 50000
nstvout                  = 50000
nstfout                  = 50000
nstlog                   = 50000
nstenergy                = 50000
nstxout-compressed       = 50000


; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10  martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw




; non-bonded electrostatic forces
rcoulomb             	 =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype         	 =     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
;PME demands for 3 extra parameters
pme_order                =     4                 
fourierspacing           =     0.12              
ewald_rtol               =     1e-05
;dieletrical constant for short and long ranges
epsilon_r                =     $(options charmm36=1 martini3=15) ;for martini3 polarizable water this should be 2.5
epsilon_rf               =     0 ; zero means infinity




; non-bonded Van Der Waals forces
rvdw                     =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
vdw_type                 =     cutoff             ; what to do after that radius
;cutoff demands for 3 extra parameters to deal with the discontinuity
vdw-modifier             =     $(options charmm36=force-switch martini3=potential-shift-verlet)
rvdw-switch              =     1.0
DispCorr                 =     $(options charmm36=EnerPres martini3=no)





; Naive velocities
continuation             = yes 
gen_vel                  = no 

; Temperature
Tcoupl                   = V-rescale
tc_grps                  = system
ref_t                    = ${current_t} ;K 
tau_t                    = 1

; Pressure
Pcoupl                   = C-rescale
ref_p                    = 1.0  ;bar
tau_p                    = $(options charmm36=5      martini3=12)
compressibility          = $(options charmm36=4.5e-5 martini3=3e-4)
refcoord_scaling         = com


; making bonds stiff       (to avoid the need of calculating fast vibrations. something that would require dividing ts by 4!)
constraints              = $(options charmm36=h-bonds martini3=none)
constraint-algorithm     = lincs
lincs-order              = 4

; restraints configuration
disre = ${distance_restraints_option}
disre_fc = 1000

dihre = ${dihedral_restraints_option}
dihre_fc = 1000


; how to restraint protein position, and how to make water flexible
; define = -DPOSRES -DFLEXIBLE

EOT

done # end of loop that creates mdps for all temperatures, inside their respective folders
echo "all mdp files created"



GRO_ABS=$(realpath -- "$GRO")                    # Preserve the coordinate-file location before changing directory.
TOP_ABS=$(realpath -- "$TOP")                    # Preserve the topology-file location before changing directory.
REPLICA_ROOT_ABS=$(realpath -- "$REPLICA_ROOT")  # Obtain the absolute replica-exchange directory.
cd -- "$REPLICA_ROOT_ABS" || exit 1              # Enter the replica-exchange directory. e.g. runREX_${NAME}



REPLEX=500         # steps between exchanges 
RESEED=123




##########################################################################################################
cat <<'EOT' > "script.${NAME}.sh"
#!/bin/bash

cd -- "$(dirname -- "${BASH_SOURCE[0]}")" || exit 1  # Always run from the folder containing this script.

EOT


##########################################################################################################
if [[ $ARCHITECTURE == "rome" ]]; then # insert the rome header, if the user chose this architecture

cat <<EOT >> "script.${NAME}.sh"

#MSUB   -r ${NAME}.repl       # Job name
#MSUB   -n ${NTMPI}                # Number of tasks in parallel mode
#MSUB   -c ${NTOMP}                       # Number of cores per parallel task
#MSUB   -W yes                     # Let multiple jobs sharing same name & user run simultaneously
#MSUB   -o out.scheduler.%I.${NAME}            # Output file
#MSUB   -e err.scheduler.%I.${NAME}            # Output file for errors
#MSUB   -q rome                    # Partition:    rome        
#MSUB   -A gen13458                # Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store      # File system:  scratch,work,store
#MSUB   -Q normal                  # Quality of Service (test,normal,long) (ccc_mqinfo)
#MSUB   -T 86400                   # Maximum walltime in seconds

set -x # echo commands

module purge  # retire tous les modules déchargeables de l'environnement
module load gnu/11 # charge gnu/11 et définit gnu/11 comme compilateur dans votre environnement
module load nvhpc/24.3 # besoin de mettre avant OpenMPI comme ce dernier charge un cuda qui n'est pas compatible avec nvhpc/24.3
module load mpi/openmpi/4 # charge la souche OpenMPI
module load gromacs/2025.0 # charge le produit

export GMX_DISABLE_GPU_DETECTION=1 # prevent GROMACS from using GPUs
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto

export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads
export OMP_DYNAMIC=FALSE


# ---- 24h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain. I'll know this checking for a done.txt file, that is created in the post processing
if [[ -f "done.txt" ]]; then
    echo "Simulation already complete. Exiting."
    exit 0
fi
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
ccc_msub -E "--dependency=afterok:\${BRIDGE_MSUB_JOBID}" script.${NAME}.sh
# --------------------------------------------------------------



EOT

GMX="gmx_mpi"                                           # Serial GROMACS tools such as grompp and dump.
MDRUN="ccc_mprun gmx_mpi mdrun"                         # External-MPI REMD launch using the allocated TGCC ranks.
#MDRUN_OPTIONS="-dd 3 3 3 -npme 13 -dlb yes" #ATENTION: this must be coherent with #MSUB -n 40 #domain decomposition dont work with steep 
#MDRUN_OPTIONS="-nt 1" #nt cant be used in rome, just set OMP_NUM_THREADS and #MSUB -n
MDRUN_OPTIONS="-maxh 23 -cpi"   # stop cleanly at ~23h, auto-continue from checkpoint




#########################################################################################################
elif [[ $ARCHITECTURE == "genoa" ]]; then # insert the adastra-genoa header, if the user chose this architecture

cat <<EOT >>  "script.${NAME}.sh"

#SBATCH --account=c1613458
#SBATCH -J ${NAME}.repl
#SBATCH --constraint=GENOA         # GENOA(192)(CPU) or MI250(64)(GPU)
##SBATCH --nodes=
#SBATCH --ntasks-per-node=${NTMPI} 
#SBATCH --cpus-per-task=${NTOMP}
##SBATCH --exclusive
#SBATCH -o ${NAME}.repl.scheduler.out
#SBATCH -e ${NAME}.repl.scheduler.err 


module purge
#module load CCE-CPU-4.0.0
develop CCE-CPU-5.0.0
module spider gromacs/2025.2-omp-mpi
#module load gromacs/2024.3-omp-mpi

module list

# Let OpenMP + gmx handle pinning (see -pin off below)
export OMP_PROC_BIND=CLOSE
export OMP_PLACES=THREADS
export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)




# ---- 24h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain (done.txt is made in post-processing).
if [[ -f "done.txt" ]]; then
    echo "Simulation already complete. Exiting."
    exit 0
fi
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
sbatch "--dependency=afterok:\${SLURM_JOB_ID}" script.${NAME}.sh
# --------------------------------------------------------------




EOT

GMX="gmx_mpi"                                                                   # External-MPI GROMACS (grompp, dump).
MDRUN="srun --cpus-per-task=${NTOMP} --threads-per-core=1 gmx_mpi mdrun"        # srun spreads the NTMPI ranks.
MDRUN_OPTIONS="-ntomp ${NTOMP} -pin off -maxh 23 -cpi"   


#########################################################################################################
elif [[ $ARCHITECTURE == "MI300" ]]; then # insert the MI300-genoa header, if the user chose this architecture

cat <<EOT >>  "script.${NAME}.sh"

#SBATCH --account=c1613458 #cad17773
#SBATCH --job-name=${NAME}.repl
#SBATCH --constraint=MI300
#SBATCH --ntasks-per-node=${NTMPI} 
#SBATCH --cpus-per-task=${NTOMP}
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -o ${NAME}.repl.scheduler.out
#SBATCH -e ${NAME}.repl.scheduler.err 




module purge

module load develop
module use /lus/work/CT7/cad17773/SHARED/Configuration.spack-user-develop/modules/tcl/linux-rhel9-x86_64
module use /lus/work/CT7/cad17773/SHARED/Configuration.spack-user-develop/modules/tcl/linux-rhel9-zen3
module use /lus/work/CT7/cad17773/SHARED/Configuration.spack-user-develop/modules/tcl/linux-rhel9-zen4
module load cce/20.0.0/zen4/gromacs/2026.1-aux5
module list


export OMP_PROC_BIND=CLOSE
export OMP_PLACES=THREADS
export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
export MPICH_GPU_SUPPORT_ENABLED=1

# ---- 24h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain (done.txt is made in post-processing).
if [[ -f "done.txt" ]]; then
    echo "Simulation already complete. Exiting."
    exit 0
fi
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
sbatch "--dependency=afterok:\${SLURM_JOB_ID}" script.${NAME}.sh
# --------------------------------------------------------------




EOT

GMX="gmx_mpi"                                                                   # External-MPI GROMACS (grompp, dump).
MDRUN="srun --cpus-per-task=${NTOMP} --threads-per-core=1 gmx_mpi mdrun"        # srun spreads the NTMPI ranks.
MDRUN_OPTIONS="-ntomp ${NTOMP} -pin off -maxh 23 -cpi"   



##########################################################################################################
elif [[ $ARCHITECTURE == "slurm" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >> "script.${NAME}.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=${NTOMP}
#SBATCH --ntasks=${NTMPI}
#SBATCH --gres=gpu:1
##SBATCH --mem-per-cpu=1GB
##SBATCH --nodes=1
#SBATCH --job-name=${NAME}.realistic
#SBATCH --output=outanderr.slurm.${NAME}
#SBATCH --exclude=node-15

module purge
module load cuda/11.8
module load gromacs/2024.5


# ---- 48h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain. I'll know this checking for a done.txt file, that is created in the post processing
if [[ -f "done.txt" ]]; then
    echo "Simulation already complete. Exiting."
    exit 0
fi
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
sbatch "--dependency=afterok:\${SLURM_JOB_ID}" script.${NAME}.sh
# --------------------------------------------------------------



EOT

GMX="gmx_mpi"                                                   # External-MPI GROMACS executable.
MDRUN="srun gmx_mpi mdrun"                                     # Slurm launches the allocated MPI ranks.
MDRUN_OPTIONS="-ntomp ${NTOMP} -maxh 47 -cpi" #ATENTION: -ntomp and -ntmpi must be coherent with #SBATCH --cpus-per-task







##########################################################################################################
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >> "script.${NAME}.sh"

module purge
module load cuda/11.8
module load gromacs/2024.5

EOT

GMX="gmx_mpi"                                                   # External-MPI GROMACS executable.
MDRUN="mpirun -np ${NTMPI} gmx_mpi mdrun"                      # Launch the requested total number of MPI ranks.
MDRUN_OPTIONS="-ntomp ${NTOMP} -cpi"           # OpenMP threads per MPI rank and checkpoint continuation.





fi # end of if that inserts preparations before the gromacs commands
##########################################################################################################


cat <<EOT >> "script.${NAME}.sh"

set -o pipefail  

GRO_ABS=$(printf '%q' "$GRO_ABS")  # Starting coordinates, safely embedded as an absolute path.
TOP_ABS=$(printf '%q' "$TOP_ABS")  # Topology, safely embedded as an absolute path.

EXPECTED_REPLICAS=${N_REPLICAS}    # Number of replica folders created by the outer script.
TOTAL_MPI_RANKS=${NTMPI}           # Total external MPI ranks requested from the scheduler or mpirun.

##### function to check if a simulation reached the planned number of steps #####
planned_steps_reached() {
    local TPR="\$1" CPT="\$2"                                                        # the inputs are the TPR filename, and checkpoint filename.

    [[ -s "\$TPR" && -s "\$CPT" ]] || return 1                                       # Return false if one of the input files is missing or empty.

    local CURRENT_STEP PLANNED_STEPS
    CURRENT_STEP=\$(${GMX} dump -cp "\$CPT" 2>/dev/null | awk -F= '/^[[:space:]]*step[[:space:]]*=/{gsub(/[[:space:]]/,"",\$2); print \$2; exit}')
    PLANNED_STEPS=\$(${GMX} dump -s "\$TPR" 2>/dev/null | awk -F= '/^[[:space:]]*nsteps[[:space:]]*=/{gsub(/[[:space:]]/,"",\$2); print \$2; exit}')

    [[ "\$CURRENT_STEP" =~ ^[0-9]+$ && "\$PLANNED_STEPS" =~ ^[0-9]+$ ]] || return 1  # Return false if one of the value is empty or negative 
    (( CURRENT_STEP >= PLANNED_STEPS ))                                              # Return true if the planned number of steps has been reached.
}

##### function to inspect gmx outanderr file, looking for failure messages #####
gmx_failed() {
    local LABEL="\$1"
    local OUTANDERR_FILE="\$2"
    local PATTERN

    PATTERN='fatal[[:space:]]+error|error[[:space:]]+in[[:space:]]+user[[:space:]]+input|ERROR[[:space:]]+[1-9][0-9]*|there (was|were) [1-9][0-9]* errors?|inconsistency[[:space:]]+in[[:space:]]+user[[:space:]]+input|too[[:space:]]+many[[:space:]]+warnings|failed|failure|assertion[[:space:]]+failed|segmentation[[:space:]]+fault|floating[[:space:]]+point[[:space:]]+exception|bus[[:space:]]+error|core[[:space:]]+dumped|aborted|killed|out[[:space:]]+of[[:space:]]+memory|cannot[[:space:]]+allocate[[:space:]]+memory|permission[[:space:]]+denied|no[[:space:]]+such[[:space:]]+file|cannot[[:space:]]+open|could[[:space:]]+not[[:space:]]+be[[:space:]]+opened|command[[:space:]]+not[[:space:]]+found'


    if LC_ALL=C grep -Eiq "\$PATTERN" "\$OUTANDERR_FILE"; then
        echo "error during \$LABEL. this is reported in \$OUTANDERR_FILE — stopping script."
        return 0
    fi

    return 1
}







echo "#############################################################"
echo "######################### gromacs commands ##################"
echo "#############################################################"

# Verify that the external-MPI GROMACS executable is available.
command -v ${GMX} >/dev/null 2>&1 || { echo "ERROR: ${GMX} was not found; replica exchange requires an external-MPI GROMACS build."; exit 1; }

# Obtain the replica directories in their zero-padded numerical order, which is also the ascending-temperature order.
shopt -s nullglob; REPLICA_DIRS=(replica_[0-9][0-9][0-9]_T*K); shopt -u nullglob

# Confirm that all expected replica folders were found.
(( \${#REPLICA_DIRS[@]} == EXPECTED_REPLICAS )) || { echo "ERROR: expected \$EXPECTED_REPLICAS replica folders, but found \${#REPLICA_DIRS[@]}."; exit 1; }

# Confirm that the total MPI-rank count can be divided equally among the replicas.
(( TOTAL_MPI_RANKS >= EXPECTED_REPLICAS && TOTAL_MPI_RANKS % EXPECTED_REPLICAS == 0 )) || { echo "ERROR: NTMPI=\$TOTAL_MPI_RANKS must be a multiple of the \$EXPECTED_REPLICAS replicas."; exit 1; }

echo "Replicas: \${#REPLICA_DIRS[@]}"
echo "Total MPI ranks: \$TOTAL_MPI_RANKS"
echo "MPI ranks per replica: \$((TOTAL_MPI_RANKS / EXPECTED_REPLICAS))"
printf '  %s\n' "\${REPLICA_DIRS[@]}"

# Stop immediately when every replica has already reached its planned number of steps.
all_done=yes; for replica_dir in "\${REPLICA_DIRS[@]}"; do planned_steps_reached "\$replica_dir/prod.tpr" "\$replica_dir/prod.cpt" || { all_done=no; break; }; done
[[ "\$all_done" == yes ]] && { touch done.txt; echo "All replicas have reached their planned number of steps."; exit 0; }

# Refuse to continue when only some replicas have checkpoints, because all REMD replicas must remain synchronized.
checkpoint_count=0; for replica_dir in "\${REPLICA_DIRS[@]}"; do [[ -s "\$replica_dir/prod.cpt" ]] && ((checkpoint_count += 1)); done
(( checkpoint_count == 0 || checkpoint_count == EXPECTED_REPLICAS )) || { echo "ERROR: only \$checkpoint_count/\$EXPECTED_REPLICAS replicas have prod.cpt; refusing an inconsistent restart."; exit 1; }

# Resolve the topology directory so that relative #include paths in the topology continue to work.
TOP_DIR=\$(dirname -- "\$TOP_ABS")
TOP_FILE=\$(basename -- "\$TOP_ABS")

# Generate each missing prod.tpr, while preserving any TPR already associated with a checkpoint.
for replica_dir in "\${REPLICA_DIRS[@]}"; do
    [[ -s "\$replica_dir/prod.mdp" ]] || { echo "ERROR: missing \$replica_dir/prod.mdp."; exit 1; }  # Every replica needs its temperature-specific MDP.
    [[ ! -s "\$replica_dir/prod.cpt" || -s "\$replica_dir/prod.tpr" ]] || { echo "ERROR: \$replica_dir has prod.cpt but no prod.tpr."; exit 1; }  # Never rebuild a TPR behind an existing checkpoint.
    [[ -s "\$replica_dir/prod.tpr" ]] && { echo "keeping existing \$replica_dir/prod.tpr"; continue; }  # Relaunches preserve the original TPR.

    replica_abs="\$PWD/\$replica_dir"
    echo "creating \$replica_dir/prod.tpr"

    # Run grompp from the topology directory so relative topology includes remain resolvable.
    if ! (cd -- "\$TOP_DIR" && ${GMX} grompp -f "\$replica_abs/prod.mdp" -c "\$GRO_ABS" -r "\$GRO_ABS" -p "\$TOP_FILE" -o "\$replica_abs/prod.tpr" -po "\$replica_abs/mdout.mdp") 2>&1 | tee "\$replica_abs/outanderr.grompp"; then
        echo "ERROR: grompp returned a failure for \$replica_dir."
        exit 1
    fi

    # Inspect grompp output and verify that its TPR was actually created.
    if gmx_failed "\$replica_dir grompp" "\$replica_abs/outanderr.grompp"; then exit 1; fi
    [[ -s "\$replica_abs/prod.tpr" ]] || { echo "ERROR: grompp finished without creating \$replica_dir/prod.tpr."; exit 1; }
done

#################################################################################################
# Launch all replicas in one external-MPI multi-simulation; each prod.cpt is interpreted relative to its replica directory.
if ! ${MDRUN} -multidir "\${REPLICA_DIRS[@]}" -deffnm prod ${MDRUN_OPTIONS} -replex ${REPLEX} -reseed ${RESEED} 2>&1 | tee outanderr.mdrun; then
    echo "ERROR: replica-exchange mdrun returned a nonzero status."
    exit 1
fi
#################################################################################################

# Inspect only the current mdrun launch; tee overwrites this diagnostic log while GROMACS appends its actual simulation outputs.
if gmx_failed "replica-exchange mdrun" "outanderr.mdrun"; then exit 1; fi

# A clean maxh stop or completed run must leave a checkpoint for every replica.
for replica_dir in "\${REPLICA_DIRS[@]}"; do [[ -s "\$replica_dir/prod.cpt" ]] || { echo "ERROR: \$replica_dir/prod.cpt was not created."; exit 1; }; done

# Mark the complete REMD calculation as done only when every replica reached the TPR target step.
all_done=yes; for replica_dir in "\${REPLICA_DIRS[@]}"; do planned_steps_reached "\$replica_dir/prod.tpr" "\$replica_dir/prod.cpt" || { all_done=no; break; }; done
[[ "\$all_done" == yes ]] && { touch done.txt; echo "All replicas completed successfully."; } || echo "Replica exchange stopped cleanly and will continue in the chained job."

echo "###################################################################"





EOT
	
chmod +x script.${NAME}.sh


if [[ $ARCHITECTURE == "slurm" ]]; then
    sbatch script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "rome" ]]; then
    ccc_msub script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "genoa" ]]; then
    sbatch script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "MI300" ]]; then
    sbatch script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "pc" ]]; then
    ./script.${NAME}.sh && echo "script finished"
    

fi












