#!/bin/bash

# usage example:  ./runREALISTIC.sh a12HW a12HW.gro a12HW.top 1 298 charmm36 pc 2 8

# ARGUMENTS:
# 1-name that goes onthe runREALISTIC to be created and the job name
# 2-gro
# 3-top
# 4-nanoseconds of production
# 5-temperature (remeber that martini3 was parametrized at 310)
# 6-forcefield to be used in mdp construction (must be "charmm36" or "martini3")
# 7-architecture (pc, oxygen, rome)
# 8-ntOMP
# 9-ntMPI

# hardcoded option to use Soft Core or Slow Groth before EM. You can set JUST ONE to "yes"
SC="no" # Soft Core
SG="no" # Slow Groth

if [ $# -lt 9 ]; then
    echo "9 arguments needed : name filename.gro filename.top numberOfNanoseconds temperatureList forceFieldName architecture ntOMP ntMPI"
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

if [[ ! -f "$TOP" ]]; then
    echo "Error: $TOP not found!"
    exit 1
fi


t=$5
echo "Temperature: ${t}"
echo " "

FF=$6
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


ARCHITECTURE=$7
echo "Architecture:${ARCHITECTURE}"

if [[ $ARCHITECTURE == "pc" || $ARCHITECTURE == "oxygen" || $ARCHITECTURE == "rome" ]]; then
    echo "Architecture is valid"
else
    echo "Error: architecture must be 'pc' or 'oxygen' or 'rome' "
    exit 1
fi
echo " "


NTOMP=$8
if [[ "$NTOMP" =~ ^[0-9]+$ ]]; then
    echo "NTOMP OK (is an integer)"
else
    echo "NTOMP is NOT an integer"
    exit 1
fi



NTMPI=$9
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

if [[ $SC == "yes" ]]; then # Soft Core option
 mkdir -p runREALISTIC_${NAME}/0_SC 
fi
if [[ $SG == "yes" ]]; then # Slow groth option
 mkdir -p runREALISTIC_${NAME}/0_SG
fi
mkdir -p runREALISTIC_${NAME}/1_EM
mkdir -p runREALISTIC_${NAME}/2_NVT
mkdir -p runREALISTIC_${NAME}/3_NPT
mkdir -p runREALISTIC_${NAME}/4_PROD
	
echo "folder structure created"

cd runREALISTIC_${NAME} || exit

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
if [[ $SC == "yes" ]]; then # Soft Core option
 file_sc_mdp="sc.mdp"
fi
if [[ $SG == "yes" ]]; then # Slow Groth option
 file_sg_mdp="sg.mdp"
fi
file_em_mdp="em.mdp"
file_nvt_mdp="nvt.mdp"
file_npt_mdp="npt.mdp"
file_prod_mdp="prod.mdp"



if [[ $SC == "yes" ]]; then #Soft Core will be used before EM
  
echo "creating ${file_sc_mdp}"
cat <<EOT > "0_SC/${file_sc_mdp}"



integrator              = steep
emtol                   = 500 
nsteps                  = 5000

; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10 martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw


; non-bonded electrostatic forces
rcoulomb                =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype             =     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
;PME demands for 3 extra parameters
pme_order               =     4                 
fourierspacing          =     0.12              
ewald_rtol              =     1e-05
;dieletrical constant for short and long ranges
epsilon_r               =     $(options charmm36=1 martini3=15) ;for martini3 polarizable water this should be 2.5
epsilon_rf              =     0 ; zero means infinity

; non-bonded Van Der Waals forces
rvdw                    =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
vdw_type                =     cutoff             ; what to do after that radius
;cutoff demands for 3 extra parameters to deal with the discontinuity
vdw-modifier            =     $(options charmm36=force-switch martini3=potential-shift-verlet)
rvdw-switch             =     1.0
DispCorr                =     $(options charmm36=EnerPres martini3=no)


; CONSTRAINTS
constraints             = h-bonds
constraint_algorithm    = LINCS

; FREE ENERGY
free-energy         = yes
couple-moltype      = System  
couple-lambda0      = vdw-q         
couple-lambda1      = none    
init-lambda         = 0.10   ;theoretical 0.50 dont work     
nstdhdl             = 0
couple-intramol     = yes

; SOFT CORE
sc-alpha            = 4  ;standard that dont work: $(options charmm36=0.5 martini3=1.3)
sc-power            = 2  ;standard that dont work: 1
sc-coul             = yes          
sc-sigma            = $(options charmm36=0.3 martini3=0.47)


EOT

fi #end of condition defining that SC will be used before EM



if [[ $SG == "yes" ]]; then #SlowGroth will be used before EM
  
echo "creating ${file_sg_mdp}"
cat <<EOT > "0_SG/${file_sg_mdp}"


integrator              = sd
dt                      = 0.002
nsteps                  = 5000
nstxtcout               = 5000
nstvout                 = 5000
nstfout                 = 5000
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10 martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw



; non-bonded electrostatic forces
rcoulomb                 =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype              =     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
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
continuation             = no 
gen_vel                  = yes
gen_temp                 = ${t}
gen_seed                 = -1

; Temperature
Tcoupl                   = V-rescale
tc_grps                  = system
ref_t                    = ${t} ;K
tau_t                    = 1

; Pressure
Pcoupl                   = C-rescale
ref_p                    = 1.0 ;bar
tau_p                    = $(options charmm36=5      martini3=12)
compressibility          = $(options charmm36=4.5e-5 martini3=3e-4)
refcoord_scaling         = com
nstcomm                  = 100
comm_mode                = linear
comm_grps                = 

; making bonds stiff       (to avoid the need of calculating fast vibrations. something that would require dividing ts by 4!)
constraints              = $(options charmm36=h-bonds martini3=none)
constraint-algorithm     = lincs
lincs-order              = 4


; Free energy variables
free-energy              = yes
couple-moltype           = System ; original OOOTG
couple-lambda0           = none
couple-lambda1           = vdw-q
couple-intramol          = yes
init-lambda              = 0
delta-lambda             = 0.0002
nstdhdl                  = 60
fep-lambdas              = 
mass-lambdas             = 
bonded-lambdas           = 
restraint-lambdas        = 
temperature-lambdas      = 
calc-lambda-neighbors    = 1
init-lambda-weights      = 
dhdl-print-energy        = no
sc-alpha                 = 4 ; $(options charmm36=0.5 martini3=1.3)
sc-power                 = 2 ; 1
;sc-r-power               = 6     ; this value came from a martini example, should this be different for charmm36?
sc-sigma                 = $(options charmm36=0.3 martini3=0.47)
sc-coul                  = no
separate-dhdl-file       = yes
dhdl-derivatives         = yes
dh_hist_size             = 0    ; this value came from a martini example, should this be different for charmm36?
dh_hist_spacing          = 0.1  ; this value came from a martini example, should this be different for charmm36?



EOT

fi #end of condition defining that SG will be used before EM



echo "creating ${file_em_mdp}"
cat <<EOT > "1_EM/${file_em_mdp}"

; Run control
integrator               = steep 
nsteps                   = 10000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.005
;niter                    = 20
;nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1

; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10 martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw


; non-bonded electrostatic forces
rcoulomb                =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype             =     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
;PME demands for 3 extra parameters
pme_order               =     4                 
fourierspacing          =     0.12              
ewald_rtol              =     1e-05
;dieletrical constant for short and long ranges
epsilon_r               =     $(options charmm36=1 martini3=15) ;for martini3 polarizable water this should be 2.5
epsilon_rf              =     0 ; zero means infinity


; non-bonded Van Der Waals forces
rvdw                	=     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
vdw_type             	=     cutoff             ; what to do after that radius
;cutoff demands for 3 extra parameters to deal with the discontinuity
vdw-modifier         	=     $(options charmm36=force-switch martini3=potential-shift-verlet)
rvdw-switch          	=     1.0
DispCorr            	=     $(options charmm36=EnerPres martini3=no)



; Velocities, Temperature and Pressure dont apply during EM
gen_vel                  = no 
continuation             = no
Tcoupl                   = no
Pcoupl                   = no




; making bonds stiff
constraints              = none 




EOT




echo "creating ${file_nvt_mdp}"
cat <<EOT > "2_NVT/${file_nvt_mdp}"

; Run control
integrator               = md
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)

nsteps                   = 50000
nstcomm                  = 100

; Output control
nstxout                  = 5000
nstvout                  = 5000
nstfout                  = 5000
nstlog                   = 5000
nstenergy                = 5000
nstxout-compressed       = 0


; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10 martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw



; non-bonded electrostatic forces
rcoulomb                 =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype              =     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
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
continuation             = no 
gen_vel                  = yes
gen_temp                 = ${t}
gen_seed                 = -1

; Temperature
Tcoupl                   = V-rescale
tc_grps                  = system
ref_t                    = ${t} ;K
tau_t                    = 1

; Pressure
Pcoupl                   = No


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
define = -DPOSRES ;-DFLEXIBLE
EOT




echo "creating ${file_npt_mdp}"
cat <<EOT > "3_NPT/${file_npt_mdp}"


; Run control
integrator               = md
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)
nsteps                   = 3000000
nstcomm                  = 100

; Output control
nstxout                  = 5000
nstvout                  = 5000
nstfout                  = 5000
nstlog                   = 5000
nstenergy                = 5000
nstxout-compressed       = 0


; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10  martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw



; non-bonded electrostatic forces
rcoulomb                 =     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype              =     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
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
DispCorr                 =     $(options charmm36=EnerPres     martini3=no)

; Naive velocities
continuation             = yes 
gen_vel                  = no 

; Temperature
Tcoupl                   = V-rescale
tc_grps                  = system
ref_t                    = ${t} ;K
tau_t                    = 1
 
; Pressure
Pcoupl                   = C-rescale
ref_p                    = 1.0 ;bar
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
define = -DPOSRES ;-DFLEXIBLE

EOT





echo "creating ${file_prod_mdp}"
cat <<EOT > "4_PROD/${file_prod_mdp}"


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
ref_t                    = ${t} ;K 
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








echo "all mdp files created"


##########################################################################################################
cat <<EOT > "script.${NAME}.sh"
#!/bin/bash

EOT


##########################################################################################################
if [[ $ARCHITECTURE == "rome" ]]; then # insert the rome header, if the user chose this architecture

cat <<EOT >> "script.${NAME}.sh"

#MSUB   -r ${NAME}.realistic       # Job name
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
if [[ -f "4_PROD/done.txt" ]]; then
    echo "Simulation already complete. No need for follow-up job"
else
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
ccc_msub -E "--dependency=afterany:\${BRIDGE_MSUB_JOBID}" script.${NAME}.sh
fi
# --------------------------------------------------------------



EOT

GMX="ccc_mprun gmx_mpi"
#MDRUN_OPTIONS="-dd 3 3 3 -npme 13 -dlb yes" #ATENTION: this must be coherent with #MSUB -n 40 #domain decomposition dont work with steep 
#MDRUN_OPTIONS="-nt 1" #nt cant be used in rome, just set OMP_NUM_THREADS and #MSUB -n
MDRUN_OPTIONS="-maxh 23 -cpi"   # stop cleanly at ~23h, auto-continue from checkpoint



##########################################################################################################
elif [[ $ARCHITECTURE == "oxygen" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >> "script.${NAME}.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=${NTOMP}
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
if [[ -f "4_PROD/done.txt" ]]; then
    echo "Simulation already complete. No need for follow-up job"
 
else
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
THIS_SCRIPT_PATH="\$(readlink -f "\$0")"
sbatch --dependency=afterany:\$SLURM_JOB_ID "\$THIS_SCRIPT_PATH"

fi
# --------------------------------------------------------------



EOT

GMX="gmx"
MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI} -maxh 47 -cpi" #ATENTION: -ntomp and -ntmpi must be coherent with #SBATCH --cpus-per-task







##########################################################################################################
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >> "script.${NAME}.sh"

module purge
module load cuda/11.8
module load gromacs/2024.5

EOT

GMX="gmx"
MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI}"





fi # end of if that inserts preparations before the gromacs commands
##########################################################################################################


cat <<EOT >> "script.${NAME}.sh"

set -o pipefail  



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

if [[ ${SC} == "yes" ]]; then  #use SC before EM
   echo "#############################################################"
   echo "######################### SC grompp #########################"
   echo "#############################################################"
   module unload gromacs
   module load gromacs/2023 # this is an ugly fix because in newer versions the SC run will give an arror related to decupling and coulomb. But here that error, whatever it means, is irrelevant. Im just trying to get rid of superpositions.
   
   cd 0_SC || exit
          
   if [[ -s "sc.gro" ]]; then
     echo "skipping sc (sc.gro already there)"
   else
   
   ${GMX} grompp -f ${file_sc_mdp} -c "../../${GRO}" -p "../../${TOP}" -o "sc.tpr" -maxwarn 2 2>&1 | tee "outanderr.grompp"
   if gmx_failed "sc" "outanderr.grompp"; then exit 1; fi
   if [[ ! -f "sc.tpr" ]]; then echo "sc finished without saving a TPR file"; exit 1; fi
   
   
   echo "########################## SC mdrun #############################"
   
   
   ${GMX} mdrun -v -deffnm "sc" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
   if [[ ! -f "sc.gro" ]]; then echo "sc finished without saving a GRO file"; exit 1; fi
   
   
   
   module unload gromacs
   module load gromacs/2024.5
   fi

fi #end of condition to use SC before EM



if [[ ${SG} == "yes" ]]; then  #use SG before EM
   echo "#############################################################"
   echo "######################### SG grompp #########################"
   echo "#############################################################"

   
   cd 0_SG || exit
          
   if [[ -s "sg.gro" ]]; then
     echo "skipping sg (sg.gro already there)"
   else
   
   ${GMX} grompp -f ${file_sg_mdp} -c "../../${GRO}" -p "../../${TOP}" -o "sg.tpr" -maxwarn 2 2>&1 | tee "outanderr.grompp"
   if gmx_failed "sg" "outanderr.grompp"; then exit 1; fi
   if [[ ! -f "sg.tpr" ]]; then echo "sg finished without saving a TPR file"; exit 1; fi
   
   
   echo "########################## SG mdrun #############################"
   
   
   ${GMX} mdrun -v -deffnm "sg" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
   if [[ ! -f "sg.gro" ]]; then echo "sg finished without saving a GRO file"; exit 1; fi
   
   

   fi

fi #end of condition to use SG before EM

echo "#############################################################"
echo "######################### EM grompp #########################"
echo "#############################################################"

if [[ ${SC} == "yes" ]]; then #SC was used before EM
    cd ../1_EM || exit #come from 0_SC
elif [[ ${SG} == "yes" ]]; then #SG was used before EM
    cd ../1_EM || exit #come from 0_SG
else
    cd 1_EM || exit    #begin at 1_EM
fi

if [[ -s "em.gro" ]]; then
    echo "skipping em (em.gro already there)"
else



if [[ ${SC} == "yes" ]]; then #SC was used before EM
  ${GMX} grompp -f ${file_em_mdp} -c "../0_SC/sc.gro" -p "../../${TOP}" -o "em.tpr" 2>&1 | tee "outanderr.grompp"
elif [[ ${SG} == "yes" ]]; then #SG was used before EM
  ${GMX} grompp -f ${file_em_mdp} -c "../0_SG/sg.gro" -p "../../${TOP}" -o "em.tpr" 2>&1 | tee "outanderr.grompp"
else
  ${GMX} grompp -f ${file_em_mdp} -c "../../${GRO}" -p "../../${TOP}" -o "em.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp" 
fi

  if gmx_failed "em" "outanderr.grompp"; then exit 1; fi
  if [[ ! -f "em.tpr" ]]; then echo "em finished without saving a TPR file"; exit 1; fi

  
  echo "######################### EM mdrun ##############################"
  ${GMX} mdrun -deffnm "em" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
  if gmx_failed "em" "outanderr.mdrun"; then exit 1; fi
  if [[ ! -f "em.gro" ]]; then echo "em finished without saving a GRO file"; exit 1; fi



fi
echo "#############################################################"
echo "######################### NVT grompp ########################"
echo "#############################################################"
cd ../2_NVT || exit
if planned_steps_reached "nvt.tpr" "nvt.cpt"; then
    echo "skipping nvt (steps reached)"
else




  ${GMX} grompp -f ${file_nvt_mdp} -c "../1_EM/em.gro" -r "../1_EM/em.gro" -p "../../${TOP}" -o "nvt.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"
  if gmx_failed "nvt" "outanderr.grompp"; then exit 1; fi
  if [[ ! -f "nvt.tpr" ]]; then echo "nvt finished without saving a TPR file"; exit 1; fi

  
  echo "######################### NVT mdrun ##############################"
  ${GMX} mdrun -v -deffnm "nvt" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
  if gmx_failed "nvt" "outanderr.mdrun"; then exit 1; fi
  if [[ ! -f "nvt.gro" ]]; then echo "nvt finished without saving a GRO file"; exit 1; fi



fi
echo "##############################################################"
echo "######################### NPT grompp #########################"
echo "##############################################################"
cd ../3_NPT || exit
if planned_steps_reached "npt.tpr" "npt.cpt"; then
    echo "skipping npt (steps reached)"
else




  ${GMX} grompp -f ${file_npt_mdp} -c "../2_NVT/nvt.gro" -r "../2_NVT/nvt.gro" -p "../../${TOP}" -o "npt.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"
  if gmx_failed "npt" "outanderr.grompp"; then exit 1; fi
  if [[ ! -f "npt.tpr" ]]; then echo "npt finished without saving a TPR file"; exit 1; fi

  
  echo "######################### NPT mdrun ##############################"
  ${GMX} mdrun -v -deffnm "npt" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
  if gmx_failed "npt" "outanderr.mdrun"; then exit 1; fi
  if [[ ! -f "npt.gro" ]]; then echo "npt finished without saving a GRO file"; exit 1; fi
  
  
  
fi
echo "#############################################################"
echo "#################### PRODUCTION grompp ######################"
echo "#############################################################"
cd ../4_PROD || exit
if planned_steps_reached "prod.tpr" "prod.cpt"; then
    echo "skipping prod (steps reached)"
else




  ${GMX} grompp -f ${file_prod_mdp} -c "../3_NPT/npt.gro" -p "../../${TOP}" -o "prod.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"
  if gmx_failed "prod" "outanderr.grompp"; then exit 1; fi
  if [[ ! -f "prod.tpr" ]]; then echo "prod finished without saving a TPR file"; exit 1; fi

  
  echo "################### PRODUCTION mdrun ############################"
  ${GMX} mdrun -v -deffnm "prod" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
  if gmx_failed "prod" "outanderr.mdrun"; then exit 1; fi
  if [[ ! -f "prod.gro" ]]; then echo "prod finished without saving a GRO file"; exit 1; fi



fi
echo "#################### CENTER AND FIT ###############################"
if planned_steps_reached "prod.tpr" "prod.cpt"; then # only post-process once PROD has truly finished
    touch "done.txt"



    printf '1\n0' | ${GMX} trjconv -s "prod.tpr" -f "prod.xtc" -o "prod.centered.xtc" -center -pbc mol 2>&1 | tee "outanderr.center"
    if gmx_failed "center" "outanderr.center"; then exit 1; fi


    printf '1\n0' | ${GMX} trjconv -s "prod.tpr" -f "prod.centered.xtc" -o "prod.fitted.xtc" -fit progressive 2>&1 | tee "outanderr.fit"
    if gmx_failed "fit" "outanderr.fit"; then exit 1; fi
    


    printf '1\n0' | ${GMX} trjconv -f prod.gro -s prod.tpr -pbc mol -o prod.whole.gro


    printf 'Density\n\n' | ${GMX} energy -f prod.edr -o density20toEND.xvg -b 20000 2>&1 | tee "outanderr.density"





    
fi
echo "###################################################################"

cd ../.. # get out of runREALISTIC





EOT
	
chmod +x script.${NAME}.sh


if [[ $ARCHITECTURE == "oxygen" ]]; then
    sbatch script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "rome" ]]; then
    ccc_msub script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "pc" ]]; then
    ./script.${NAME}.sh && echo "script finished"
    

fi












