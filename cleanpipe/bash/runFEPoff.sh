#!/bin/bash

# usage example: ./runFEP.sh alaHO alaHO.gro alaHO.top 100 "283 298 313" charmm36 pc molecule_0 8 2


# ARGUMENTS:
# 1-name that goes on the runFEP the folder to be created, and job name
# 2-gro
# 3-top
# 4-nanoseconds of production
# 5-temperatures (remeber that martini3 was parametrized at 310)
# 6-forcefield to use in mdp construction (must be "charmm36" or "martini3")
# 7-architecture (pc, slurm, rome)
# 8-molecule to be decoupled
# 9-ntOMP
# 10-ntMPI

if [ $# -lt 10 ]; then
    echo "10 arguments needed : name filename.gro filename.top numberOfNanoseconds temperatureList forceFieldName architecture moleculeName ntOMP ntMPI"
    exit 1
fi






######## obtain the arguments defined by the user when he called the function ########
NAME=$1 #name of the system to be simulated
echo " "
echo "Name of system: ${NAME}"
echo " "

GRO=$2 #name of gro
TOP=$3 #name of top
echo "GRO: ${GRO}"
echo " "
echo "TOP: ${TOP}"
echo " "

if [[ ! -f "$GRO" ]]; then
    echo "Error: $GRO not found!"
    exit 1
fi

if [[ ! -f "$TOP" ]]; then
    echo "Error: $TOP not found!"
    exit 1
fi








TEMPERATURE_LIST=$5
echo "Temperature list: ${TEMPERATURE_LIST}"
echo " "



FF=$6
echo "Force field: ${FF}"
echo " "

if [[ $FF == "charmm36" || $FF == "martini3" ]]; then
    echo "    FF is valid"
else
    echo "Error: the force field must be either 'charmm36' or 'martini3'"
    exit 1
fi
echo " "



PRODUCTION_DURATION=$4 #how many nanoseconds
if [[ $FF == "charmm36" ]]; then
  STEPS=$(echo "scale=0; ($PRODUCTION_DURATION / 0.002) * 10000" | bc)
fi

if [[ $FF == "martini3" ]]; then
  STEPS=$(echo "scale=0; ($PRODUCTION_DURATION / 0.02) * 10000" | bc)
fi


echo "Production Duration: ${PRODUCTION_DURATION} ns (implemented in the mdp by setting: ${STEPS} steps x dt )"
echo " "


ARCHITECTURE=$7
echo "Architecture: ${ARCHITECTURE}"
echo " "
if [[ $ARCHITECTURE == "pc" || $ARCHITECTURE == "slurm" || $ARCHITECTURE == "rome" ]]; then
    echo "    Architecture is valid"
else
    echo "Error: architecture must be 'pc' or 'slurm' or 'rome' "
    exit 1
fi
echo " "


MOL_TO_DECOUPLE=$8
echo "Molecule to decouple: ${MOL_TO_DECOUPLE}"
echo " "

NTOMP=$9
if [[ "$NTOMP" =~ ^[0-9]+$ ]]; then
    echo "NTOMP OK (is an integer)"
else
    echo "NTOMP is NOT an integer"
fi



NTMPI=$10
if [[ "$NTMPI" =~ ^[0-9]+$ ]]; then
    echo "NTMPI OK (is an integer)"
else
    echo "NTMPI is NOT an integer"
fi
echo " "







########## check top file for [ distance_restraints ] or [ dihedral_restraints ]

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











echo "####################### create folder structure ########################"


mkdir "${NAME}__runFEP" || exit 1
cd "${NAME}__runFEP" || exit 1

for t in $TEMPERATURE_LIST; do

    mkdir "t${t}"
    cd "t${t}" || exit 1



    for i in $(seq -w 0 20); do
	mkdir "Lambda_${i}"
	cd "Lambda_${i}"  || exit 1

	mkdir 1_EM
	mkdir 2_NVT
	mkdir 3_NPT
	mkdir 4_PROD
	
	cd .. #back to T${t} folder
    done # lambda loop

    cd .. #back to __runFEP folder
done # temperature loop

echo "folder structure created"
echo " "

echo "######################## generate mdp files #############################"
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




for t in $TEMPERATURE_LIST; do
for i in $(seq -w 0 20); do # from 00 to 20. this impacts the mdp pararamenter init_lambda_state


# Definition of the name of the mdp files for the current temperature and lambda
file_em_mdp="em_FEP_t${t}_lambda${i}.mdp"
file_nvt_mdp="nvt_FEP_t${t}_lambda${i}.mdp"
file_npt_mdp="npt_FEP_t${t}_lambda${i}.mdp"
file_prod_mdp="prod_FEP_t${t}_lambda${i}.mdp"



echo "creating ${file_em_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/1_EM/${file_em_mdp}"

; Run control
integrator               = steep 
nsteps                   = 100000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
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


; Decoupling control
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1                ; only immediate neighboring windows
couple-moltype           = ${MOL_TO_DECOUPLE}  ; name of molecule to decouple
couple-lambda0           = vdw-q            ; at lambda 0: interactions are ON  
couple-lambda1           = none             ; at lambda 1: interactions are OFF
couple-intramol          = yes              ; for big molecules, it’s necessary to turn off the internal interactions, not just the molecule-surroundings interactions.
nstdhdl                  = 10

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; decouple gently         (to avoid singularities during decoupling)
sc-alpha                 = $(options charmm36=0.5 martini3=1.3)
sc-coul                  = no          
sc-power                 = 1
sc-sigma                 = $(options charmm36=0.3 martini3=0.47)



; Velocities, Temperature and Pressure dont apply during EM
gen_vel                  = no 
continuation             = no
Tcoupl                   = no
Pcoupl                   = no




; making bonds stiff
constraints              = none 


EOT




echo "creating ${file_nvt_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/2_NVT/${file_nvt_mdp}"

; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)

nsteps                   = 50000
nstcomm                  = 50000

; Output control
nstxout                  = 50000
nstvout                  = 50000
nstfout                  = 0
nstlog                   = 50000
nstenergy                = 50000
nstxout-compressed       = 0


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



; Decoupling control
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1                ; only immediate neighboring windows
couple-moltype           = ${MOL_TO_DECOUPLE}  ; name of molecule to decouple
couple-lambda0           = vdw-q            ; at lambda 0: interactions are ON  
couple-lambda1           = none             ; at lambda 1: interactions are OFF
couple-intramol          = yes              ; for big molecules, it’s necessary to turn off the internal interactions, not just the molecule-surroundings interactions.
nstdhdl                  = 10

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; decouple gently         (to avoid singularities during decoupling)
sc-alpha                 = $(options charmm36=0.5 martini3=1.3)
sc-coul                  = no          
sc-power                 = 1
sc-sigma                 = $(options charmm36=0.3 martini3=0.47)


; Naive velocities
continuation             = no 
gen_vel                  = yes
gen_temp                 = ${t}
gen_seed                 = -1

; Temperature
;Tcoupl is implicitly handled by the sd integrator
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
; define = -DPOSRES -DFLEXIBLE
EOT




echo "creating ${file_npt_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/3_NPT/${file_npt_mdp}"


; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)
nsteps                   = 100000
nstcomm                  = 50000

; Output control
nstxout                  = 50000
nstvout                  = 50000
nstfout                  = 0
nstlog                   = 50000
nstenergy                = 50000
nstxout-compressed       = 0


; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10  martini3=20)
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
DispCorr                =     $(options charmm36=EnerPres     martini3=no)


; Decoupling control
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1                ; only immediate neighboring windows
couple-moltype           = ${MOL_TO_DECOUPLE}  ; name of molecule to decouple
couple-lambda0           = vdw-q            ; at lambda 0: interactions are ON  
couple-lambda1           = none             ; at lambda 1: interactions are OFF
couple-intramol          = yes              ; for big molecules, it’s necessary to turn off the internal interactions, not just the molecule-surroundings interactions.
nstdhdl                  = 10

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; decouple gently         (to avoid singularities during decoupling)
sc-alpha                 = $(options charmm36=0.5 martini3=1.3)
sc-coul                  = no          
sc-power                 = 1
sc-sigma                 = $(options charmm36=0.3 martini3=0.47)






; Naive velocities
continuation             = yes 
gen_vel                  = no 

; Temperature
;Tcoupl is implicitly handled by the sd integrator
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
; define = -DPOSRES -DFLEXIBLE

EOT





echo "creating ${file_prod_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/4_PROD/${file_prod_mdp}"


; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)
nsteps                   = ${STEPS}
nstcomm                  = 100

; Output control
nstxout                  = 500000
nstvout                  = 500000
nstfout                  = 0
nstlog                   = 500000
nstenergy                = 5000
nstxout-compressed       = 500000


; box config
pbc                      = xyz


; neighborsearch algorith
cutoff-scheme            = verlet
nstlist                  = $(options charmm36=10  martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1) ;this has to match the rcoulomb and rvdw




; non-bonded electrostatic forces
rcoulomb            	=     $(options charmm36=1.2 martini3=1.1)                ; ideal potential up to this radius
coulombtype         	=     $(options charmm36=PME martini3=reaction-field)     ; what to do after that radius
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


; Decoupling control
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1                ; only immediate neighboring windows
couple-moltype           = ${MOL_TO_DECOUPLE}  ; name of molecule to decouple
couple-lambda0           = vdw-q            ; at lambda 0: interactions are ON  
couple-lambda1           = none             ; at lambda 1: interactions are OFF
couple-intramol          = yes              ; for big molecules, it’s necessary to turn off the internal interactions, not just the molecule-surroundings interactions.
nstdhdl                  = 10

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; decouple gently         (to avoid singularities during decoupling)
sc-alpha                 = $(options charmm36=0.5 martini3=1.3)
sc-coul                  = no          
sc-power                 = 1
sc-sigma                 = $(options charmm36=0.3 martini3=0.47)




; Naive velocities
continuation             = yes 
gen_vel                  = no 

; Temperature
;Tcoupl is implicitly handled by the sd integrator
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





done # lambda loop
done #temperature loop

echo "all mdp files created"
echo " "




for t in $TEMPERATURE_LIST; do
for i in $(seq -w 0 20); do 

# reset the name of the mdp files for the current temperature and lambda
file_em_mdp="em_FEP_t${t}_lambda${i}.mdp"
file_nvt_mdp="nvt_FEP_t${t}_lambda${i}.mdp"
file_npt_mdp="npt_FEP_t${t}_lambda${i}.mdp"
file_prod_mdp="prod_FEP_t${t}_lambda${i}.mdp"


pwd
cd t${t}/Lambda_${i}

cat <<EOT > "t${t}.l${i}.sh"
#!/bin/bash

EOT




########################################################################################################
if [[ $ARCHITECTURE == "rome" ]]; then # insert the rome header, if the user chose this architecture


cat <<EOT >> "t${t}.l${i}.sh"

#MSUB   -r ${NAME}.${t}.${i}.fep             # Job name
#MSUB   -n ${NTMPI}                          # Number of tasks in parallel mode (ntmpi)
#MSUB   -c 1                                 # Number of cores per parallel task
#MSUB   -W yes                               # Let multiple jobs sharing same name & user run simultaneously
#MSUB   -o t${t}.l${i}.%I.scheduler.out      # Output file
#MSUB   -e t${t}.l${i}.%I.scheduler.err      # Output file for errors
#MSUB   -q rome                              # Partition:    rome        
#MSUB   -A gen13458                          # Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store                # File system:  scratch,work,store
#MSUB   -Q normal                            # Quality of Service (test,normal,long) (ccc_mqinfo)
#MSUB   -T 86400                             # Maximum walltime in seconds

set -x # echo commands

module purge  # retire tous les modules déchargeables de l'environnement
module load gnu/11 # charge gnu/11 et définit gnu/11 comme compilateur dans votre environnement
module load nvhpc/24.3 # besoin de mettre avant OpenMPI comme ce dernier charge un cuda qui n'est pas compatible avec nvhpc/24.3
module load mpi/openmpi/4 # charge la souche OpenMPI
module load gromacs/2025.0 # charge le produit

export GMX_DISABLE_GPU_DETECTION=1 # prevent GROMACS from using GPUs
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto

OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)

EOT

#GMX ENGINE 
GMX="ccc_mprun gmx_mpi"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #MSUB -n
#MDRUN_OPTIONS="-dd 3 3 3 -npme 13 -dlb yes" #this requires #MSUB -n 40 , but domain decomposition dont work with steep 
#MDRUN_OPTIONS="-nt 1" #doesbt work, because -nt, -ntomp, -ntmpi, cant be used in rome, you have to set OMP_NUM_THREADS and #MSUB -n instead
MDRUN_OPTIONS=""

#########################################################################################################
elif [[ $ARCHITECTURE == "slurm" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >> "t${t}.l${i}.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=${NTOMP}
##SBATCH --gres=gpu:1
##SBATCH --nodes=1
#SBATCH --job-name=${NAME}.${t}.${i}.fep
#SBATCH --output=t${t}.l${i}.scheduler.outanderr
#SBATCH --exclude=node-15

module purge
module load cuda/11.8
module load gromacs/2024.5

#alternative:
#module purge
#module load cuda/12.2
#module load gromacs/2025.0

EOT

#GMX ENGINE
GMX="gmx"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #SBATCH --cpus-per-task
MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI}" #this requires  #SBATCH --cpus-per-task=1

##########################################################################################################3
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >> "t${t}.l${i}.sh"

module purge
module load cuda/11.8
module load gromacs/2024.5

#alternative:
#module purge
#module load cuda/12.2
#module load gromacs/2025.0

EOT

#GMX ENGINE
GMX="gmx"

#GMX MDRUN ADITIONAL OPTIONS
#MDRUN_OPTIONS="-ntomp 16 -ntmpi 1"
#MDRUN_OPTIONS="-ntmpi 1 -ntomp 16 -gpu_id 0"
#MDRUN_OPTIONS="-nt 1 -pin on"
MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI}" 


########################################################################################################
fi # end of if that inserts script headers and module loading before the gromacs commands








#now the gromacs commands will be appended to the headers and moldule loading
cat <<EOT >> "t${t}.l${i}.sh"



set -o pipefail  # stop if any part of a pipeline fails




echo "############## EM - Temperature ${t} Lambda ${i} #########################################"
cd 1_EM || exit 1
	

${GMX} grompp -f ${file_em_mdp} -c "../../../../${GRO}" -p "../../../../${TOP}" -o "${NAME}_em.tpr" 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "##########################################################################################"

${GMX} mdrun -deffnm "${NAME}_em" ${MDRUN_OPTIONS} 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi





echo "############# NVT - Temperature ${t} Lambda ${i}  ########################################"
cd ../2_NVT || exit 1


${GMX} grompp -f ${file_nvt_mdp} -c "../1_EM/${NAME}_em.gro" -r "../1_EM/${NAME}_em.gro" -p "../../../../${TOP}" -o "${NAME}_nvt.tpr" -maxwarn 1 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "##########################################################################################"

${GMX} mdrun -v -deffnm "${NAME}_nvt" ${MDRUN_OPTIONS} 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi





echo "##################### NPT - Temperature ${t} Lambda ${i}  ###############################"
cd ../3_NPT || exit 1

${GMX} grompp -f ${file_npt_mdp} -c "../2_NVT/${NAME}_nvt.gro" -r "../2_NVT/${NAME}_nvt.gro" -p "../../../../${TOP}" -o "${NAME}_npt.tpr" -maxwarn 1 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "#########################################################################################"

${GMX} mdrun -v -deffnm "${NAME}_npt" ${MDRUN_OPTIONS} 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi




echo "#################### PRODUCTION - Temperature ${t} Lambda ${i}  #######################"
cd ../4_PROD || exit 1	

${GMX} grompp -f ${file_prod_mdp} -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${NAME}_prod_${t}_${i}.tpr" -maxwarn 1 2>&1 | tee "log"
#gmx grompp -f xxx I must create a mdp file that is the same as prod but with shorter lenght xxx  -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${NAME}_quick_${t}_${i}.tpr" -maxwarn 1

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



echo "#######################################################################################"

${GMX} mdrun -v -deffnm "${NAME}_prod_${t}_${i}" ${MDRUN_OPTIONS} 2>&1 | tee "log" 

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi




echo "############################## CENTER AND FIT - Temperature ${t} Lambda ${i}  ###################################"


printf '1\n0' | ${GMX} trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.xtc" -o "${NAME}_prod_${t}_${i}.centered.xtc" -center -pbc mol 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



printf '1\n0' | ${GMX} trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.centered.xtc" -o "${NAME}_prod_${t}_${i}.fitted.xtc" -fit progressive 2>&1 | tee "log"

if grep -q "Error" "log"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi






cd .. # go back to the current lambda
EOT
	


# xxx in the future, I should put an analysis script here, that procced with the analysis. this has to be a different script, because for slurm and rome, I should run just after all the jobs are finished



chmod +x t${t}.l${i}.sh


if [[ $ARCHITECTURE == "slurm" ]]; then
    sbatch t${t}.l${i}.sh && echo "job was sent (t: ${t} Lambda: ${i})"

elif [[ $ARCHITECTURE == "rome" ]]; then
    ccc_msub t${t}.l${i}.sh && echo "job was sent (t: ${t} Lambda: ${i})"

elif [[ $ARCHITECTURE == "pc" ]]; then
    nohup ./t${t}.l${i}.sh > t${t}.l{l}.redirected.out.and.err 2>&1 &
    echo "script was lounched (t: ${t} Lambda: ${i})"

fi











cd ../.. #go back to runFEP so we can go to the next temperature and lambda
done # lambda loop
done # temperature loop




