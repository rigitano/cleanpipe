#!/bin/bash


# usage example: 
# martini
#./runFEPoff.sh alaHO alaHO.gro alaHO.top 100 "283 298 313" martini3 pc molecule_0 2 8
# aa
#./runFEPoff.sh alaHO alaHO.gro alaHO.top 100 "283 298 313" charmm36 pc molecule_0 ???? ?????

# ARGUMENTS:
# 1-name that goes on the runFEP the folder to be created, and job name
# 2-gro
# 3-top
# 4-nanoseconds of production
# 5-temperatures (remeber that martini3 was parametrized at 310)
# 6-forcefield to use in mdp construction (must be "charmm36" or "martini3")
# 7-architecture (pc, slurm, rome, adastra)
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
    echo "Error: the force field must be either 'charmm36' or 'martini3' !"
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
echo "Architecture: ${ARCHITECTURE}"
echo " "
if [[ $ARCHITECTURE == "pc" || $ARCHITECTURE == "slurm" || $ARCHITECTURE == "rome" ]]; then
    echo "    Architecture is valid"
else
    echo "Error: architecture must be 'pc' or 'slurm' or 'rome' !"
    exit 1
fi
echo " "


MOL_TO_DECOUPLE=$8
echo "Molecule to decouple: ${MOL_TO_DECOUPLE}"
echo " "

NTOMP=$9
echo "NTOMP: ${NTOMP}"
if [[ "$NTOMP" =~ ^[0-9]+$ ]]; then
    echo "NTOMP OK (is an integer)"
else
    echo "Error: NTOMP is NOT an integer!"
    exit 1
fi



NTMPI=${10}
echo "NTMPI: ${NTMPI}"
if [[ "$NTMPI" =~ ^[0-9]+$ ]]; then
    echo "NTMPI OK (is an integer)"
else
    echo "Error: NTMPI is NOT an integer!"
    exit 1
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


cat <<EOT >>  "t${t}.l${i}.sh"

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

export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)

EOT

#GMX ENGINE 
GMX="ccc_mprun gmx_mpi"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #MSUB -n
#MDRUN_OPTIONS="-dd 3 3 3 -npme 13 -dlb yes" #this requires #MSUB -n 40 , but domain decomposition dont work with steep 
#MDRUN_OPTIONS="-nt 1" #doesbt work, because -nt, -ntomp, -ntmpi, cant be used in rome, you have to set OMP_NUM_THREADS and #MSUB -n instead
MDRUN_OPTIONS=""

#########################################################################################################
elif [[ $ARCHITECTURE == "slurm" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >>  "t${t}.l${i}.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=${NTOMP}
#SBATCH --job-name=${NAME}.${t}.${i}.fep
#SBATCH --output=t${t}.l${i}.scheduler.outanderr
#SBATCH --gres=gpu:1
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


#########################################################################################################
elif [[ $ARCHITECTURE == "adastra" ]]; then # insert the adastra header, if the user chose this architecture

cat <<EOT >>  "t${t}.l${i}.sh"

#SBATCH --account=c1613458
#SBATCH -J ${NAME}.${t}.${i}.fep
#SBATCH --constraint=GENOA         # GENOA(192)(CPU) or MI250(64)(GPU)
##SBATCH --nodes=
#SBATCH --ntasks-per-node=${NTMPI} 
#SBATCH --cpus-per-task=${NTOMP}
##SBATCH --exclusive
#SBATCH -o t${t}.l${i}.scheduler.out
#SBATCH -e t${t}.l${i}.scheduler.err 


module purge
#module load CCE-CPU-4.0.0
develop CCE-CPU-5.0.0
module spider gromacs/2025.2-omp-mpi
#module load gromacs/2024.3-omp-mpi

module list


export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)

EOT

#GMX ENGINE 
GMX="gmx_mpi"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #SBATCH
MDRUN_OPTIONS=""




##########################################################################################################3
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >>  "t${t}.l${i}.sh"

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
cat <<EOT >>  "t${t}.l${i}.sh"



set -o pipefail  # stop if any part of a pipeline fails




echo "############## EM - Temperature ${t} Lambda ${i} #########################################"
cd 1_EM || exit 1
	

${GMX} grompp -f ${file_em_mdp} -c "../../../../${GRO}" -p "../../../../${TOP}" -o "${NAME}_em.tpr" 2>&1 | tee "log.grompp"

if grep -q "Error" "log.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "##########################################################################################"

${GMX} mdrun -deffnm "${NAME}_em" 2>&1 | tee "log.mdrun"

if grep -q "Error" "log.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi





echo "############# NVT - Temperature ${t} Lambda ${i}  ########################################"
cd ../2_NVT || exit 1


${GMX} grompp -f ${file_nvt_mdp} -c "../1_EM/${NAME}_em.gro" -r "../1_EM/${NAME}_em.gro" -p "../../../../${TOP}" -o "${NAME}_nvt.tpr" -maxwarn 1 2>&1 | tee "log.grompp"

if grep -q "Error" "log.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "##########################################################################################"

${GMX} mdrun -v -deffnm "${NAME}_nvt" ${MDRUN_OPTIONS} 2>&1 | tee "log.mdrun"

if grep -q "Error" "log.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi





echo "##################### NPT - Temperature ${t} Lambda ${i}  ###############################"
cd ../3_NPT || exit 1

${GMX} grompp -f ${file_npt_mdp} -c "../2_NVT/${NAME}_nvt.gro" -r "../2_NVT/${NAME}_nvt.gro" -p "../../../../${TOP}" -o "${NAME}_npt.tpr" -maxwarn 1 2>&1 | tee "log.grompp"

if grep -q "Error" "log.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "#########################################################################################"

${GMX} mdrun -v -deffnm "${NAME}_npt" ${MDRUN_OPTIONS} 2>&1 | tee "log.mdrun"

if grep -q "Error" "log.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi




echo "#################### PRODUCTION - Temperature ${t} Lambda ${i}  #######################"
cd ../4_PROD || exit 1	

${GMX} grompp -f ${file_prod_mdp} -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${NAME}_prod_${t}_${i}.tpr" -maxwarn 1 2>&1 | tee "log.grompp"
#gmx grompp -f xxx I must create a mdp file that is the same as prod but with shorter lenght xxx  -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${NAME}_quick_${t}_${i}.tpr" -maxwarn 1

if grep -q "Error" "log.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



echo "#######################################################################################"

${GMX} mdrun -v -deffnm "${NAME}_prod_${t}_${i}" ${MDRUN_OPTIONS} 2>&1 | tee "log.mdrun" 

if grep -q "Error" "log.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi




echo "############################## CENTER AND FIT - Temperature ${t} Lambda ${i}  ###################################"


printf '1\n0' | ${GMX} trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.xtc" -o "${NAME}_prod_${t}_${i}.centered.xtc" -center -pbc mol 2>&1 | tee "log.center"

if grep -q "Error" "log.center"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



printf '1\n0' | ${GMX} trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.centered.xtc" -o "${NAME}_prod_${t}_${i}.fitted.xtc" -fit progressive 2>&1 | tee "log.fit"

if grep -q "Error" "log.fit"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi






cd .. # go back to the current lambda
EOT
	








chmod +x t${t}.l${i}.sh


if [[ $ARCHITECTURE == "slurm" ]]; then
    jid=$(sbatch t${t}.l${i}.sh | awk '{print $4}')
    #build list of ids sent, to use to set the dependency of the analysis script 
    slurm_ids+=($jid)
    DEPENDENCY_STRING=$(printf "afterok:%s:" "${slurm_ids[@]}")
    DEPENDENCY_STRING=${DEPENDENCY_STRING%:}    # remove final colon

    echo "job was sent (t: ${t} Lambda: ${i}) - id ${jid}"


elif [[ $ARCHITECTURE == "rome" ]]; then
    jid=$(ccc_msub t${t}.l${i}.sh | grep -Eo '[0-9]+' | tail -n1)
    #build list of ids sent, to use to set the dependency of the analysis script 
    tgcc_ids+=($jid)
    DEPENDENCY_STRING=$(printf "%s," "${tgcc_ids[@]}")
    DEPENDENCY_STRING=${DEPENDENCY_STRING%,}      # remove final comma

    echo "job was sent (t: ${t} Lambda: ${i}) - id ${jid}"


elif [[ $ARCHITECTURE == "pc" ]]; then
    nohup ./t${t}.l${i}.sh > t${t}.l{l}.redirected.out.and.err 2>&1 &
    pid=$!
    DEPENDENCY_STRING+=($pid)
    echo "script was lounched (t: ${t} Lambda: ${i})"

fi











cd ../.. #go back to runFEP so we can go to the next temperature and lambda
done # lambda loop
done # temperature loop




echo "Jobs were sent for all lambdas and temperatures"
echo " "
echo "This is the dependency string:"
echo "$DEPENDENCY_STRING"
echo " "








# analysis script that will be run just after all the jobs are finished. each temperature will have one
# this analysis script will concatenate outputs for each lambda, if there is *part* in the name, and then calculate the bar and baring
for t in $TEMPERATURE_LIST; do




cat <<EOT > "${NAME}.${t}.ConcatAndBar.sh"
#!/bin/bash

EOT




########################################################################################################
if [[ $ARCHITECTURE == "rome" ]]; then # insert the rome header, if the user chose this architecture


cat <<EOT >>  "${NAME}.${t}.ConcatAndBar.sh"

#MSUB   -r ${NAME}.${t}.ConcatAndBar         # Job name
#MSUB   -n 1                                 # Number of tasks in parallel mode (ntmpi)
#MSUB   -c 8                                 # Number of cores per parallel task
#MSUB   -W yes                               # Let multiple jobs sharing same name & user run simultaneously
#MSUB   -o final.analysis.%I.scheduler.out   # Output file
#MSUB   -e final.analysis.%I.scheduler.err   # Output file for errors
#MSUB   -q rome                              # Partition:    rome        
#MSUB   -A gen13458                          # Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store                # File system:  scratch,work,store
#MSUB   -Q normal                            # Quality of Service (test,normal,long) (ccc_mqinfo)
#MSUB   -T 86400                             # Maximum walltime in seconds
#MSUB   -a ${DEPENDENCY_STRING}
#MSUB   -@ henrique.rigitano@ibcp.fr:end

set -x # echo commands

module purge  # retire tous les modules déchargeables de l'environnement
module load gnu/11 # charge gnu/11 et définit gnu/11 comme compilateur dans votre environnement
module load nvhpc/24.3 # besoin de mettre avant OpenMPI comme ce dernier charge un cuda qui n'est pas compatible avec nvhpc/24.3
module load mpi/openmpi/4 # charge la souche OpenMPI
module load gromacs/2025.0 # charge le produit

export GMX_DISABLE_GPU_DETECTION=1 # prevent GROMACS from using GPUs
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto

OMP_NUM_THREADS=1      # number of OpenMP threads (ntomp)

EOT

#GMX ENGINE 
GMX="ccc_mprun gmx_mpi"



#########################################################################################################
elif [[ $ARCHITECTURE == "slurm" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >>  "${NAME}.${t}.ConcatAndBar.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=1
#SBATCH --job-name=${NAME}.${t}.ConcatAndBar
#SBATCH --output=final.analysis.scheduler.outanderr
#SBATCH --exclude=node-15
#SBATCH --dependency=${DEPENDENCY_STRING}

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


##########################################################################################################3
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >>  "${NAME}.${t}.ConcatAndBar.sh"

sleep 14400

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



########################################################################################################
fi # end of if that inserts script headers and module loading before the gromacs commands



#now the gromacs commands will be appended to the headers and moldule loading
cat <<EOT >> "${NAME}.${t}.ConcatAndBar.sh"

set -o pipefail  # stop if any part of a pipeline fails




echo '##################################################################'
echo '########## concatenation of mdrun outputs, if necessary ##########'
echo '################ (for all temperatures and lamdas) ###############'
echo '##################################################################'





for i in {00..20}; do  # lambda loop

cd t${t}/Lambda_\${i}/4_PROD || exit 1
pwd



##################### xvg concatenation, if necessary #####################


# Count part files
count=\$(find . -name "*part*.xvg" -type f | wc -l)

if [[ "\$count" -gt 0 ]]; then

    echo "Found \$count *part*.xvg files. Sorting…"

    #I will create concatenated files with name all. so lets delete them, I case I run this script before
    rm -f -- *.all.xvg


    # Sort files like in your xtc script
    sorted_files=\$(find . -name "*part*.xvg" -type f \
        | sed -E 's#.*_([0-9]+)\.([0-9]+)\.part0*([0-9]+)\.(.*)#\2 \3 & #' \
        | sort -k1,1n -k2,2n \
        | awk '{print \$3}')

    nonempty_files=()

    echo "Checking which XVG files contain data…"

    for f in \$sorted_files; do
        echo -n " → Checking \$f ... "

        # Extract last non-comment line
        lastline=\$(grep -v '^[#@]' "\$f" | tail -n 1)

        if [[ -z "\$lastline" ]]; then
            echo "EMPTY — skipping"
        else
            echo "OK"
            nonempty_files+=("\$f")
        fi
    done

    if [[ \${#nonempty_files[@]} -eq 0 ]]; then
        echo "All xvg part files are empty at t${t} L\${i} — nothing to concatenate."
        exit 1
    fi

    echo "Non-empty files to be concatenated:"
    printf '   %s\n' "\${nonempty_files[@]}"

    # Define output name based on prefix of first file
    first="\${nonempty_files[0]}"
    prefix="\${first%%.part*}"
    out="\${prefix}.all.xvg"

    echo "Writing concatenated file to: \$out"

    # Write header from first file only
    grep '^[#@]' "\${nonempty_files[0]}" > "\$out"

    echo "@    legend \"Concatenated XVG\"" >> "\$out"
    echo >> "\$out"

    # Append only data (ignore comments) from each file
    for f in "\${nonempty_files[@]}"; do
        grep -v '^[#@]' "\$f" >> "\$out"
    done

    echo "Done."

else
    echo "No part-files found. Renaming *.xvg to *.all.xvg…"
    find . -maxdepth 1 -name "*.xvg" -type f | while read -r f; do
        base="\${f%.xvg}"
        new="\${base}.all.xvg"
        echo "Renaming: \$f → \$new"
        mv "\$f" "\$new"
    done
fi






##################### xtc concatenation, if necessary #####################


# Count part files
count=\$(find . -name "*part*.xtc" -type f | wc -l)

if [[ "\$count" -gt 0 ]]; then
    echo "Found \$count *part*.xtc files. Sorting…"

    #I will create concatenated files with name all. so lets delete them, I case I run this script before
    rm -f -- *.all.xvg


    # Sort files using your existing logic
    sorted_files=\$(find . -name "*part*.xtc" -type f \
        | sed -E 's#.*_([0-9]+)\.([0-9]+)\.part0*([0-9]+)\.(.*)#\2 \3 & #' \
        | sort -k1,1n -k2,2n \
        | awk '{print \$3}')

    echo "Checking which files are empty…"

    nonempty_files=()

    for f in \$sorted_files; do
        echo -n " → Checking \$f ... "

        # Run gmx check and detect emptiness
        # An empty xtc usually shows somthing like: "Last frame read 0" or "Read 0 frames"
        frames=\$(${GMX} check -f "\$f" 2>/dev/null \
                | awk '
                    /[Ff]rame/ {
                        # Extract the last numeric value in the line
                        for (i = NF; i > 0; i--) {
                            if (\$i ~ /^[0-9]+$/) { print \$i; exit }
                        }
                    }')

        # Default to 0 if empty (means no frame count found)
        frames=\${frames:-0}

        if (( frames > 0 )); then
            echo "OK (\$frames frames)"
            nonempty_files+=("\$f")
        else
            echo "EMPTY — skipping"
        fi


    done

    if [[ \${#nonempty_files[@]} -eq 0 ]]; then
        echo "All xtc part files are empty at t${t} L\${i} — nothing to concatenate."
        exit 1
    fi

    echo "Non-empty files to be concatenated:"
    printf '   %s\n' "\${nonempty_files[@]}"

    # Extract prefix from first file (safer than hard-coding)
    first="\${nonempty_files[0]}"
    prefix="\${first%%.part*}"   # removes .part0001.xtc etc.

    out="\${prefix}.all.xtc"

    echo "Concatenating into: \$out"

    # Run gmx trjcat with only valid files
    ${GMX} trjcat -f "\${nonempty_files[@]}" -o "\$out" -settime <<\EOF
0
EOF

    echo "Done."

else
    echo "No part-files found. Renaming *.xtc to *.all.xtc…"

    find . -maxdepth 1 -name "*.xtc" -type f | while read -r f; do
        base="\${f%.xtc}"
        new="\${base}.all.xtc"
        echo "Renaming: \$f → \$new"
        mv "\$f" "\$new"
    done
fi



cd ../../..  #go back to runFEP, thext iteration will jump to the correct 4_PROD

done # lambda loop





echo '######################################################################'
echo '##################### bar and barint calculation #####################'
echo '################## (for all temperatures and lamdas) #################'
echo '######################################################################'






pwd
cd t${t} || exit 1

if [[ -d bar ]]; then
    rm -rf bar
fi

mkdir bar


${GMX} bar -f Lambda_*/4_PROD/*.all.xvg -o bar/${NAME}_${t}_bar.xvg -oi bar/${NAME}_${t}_barint.xvg 2>&1 | tee "log.bar"

if grep -q "Error" "log.bar"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi




cd .. #go back to runFEP so we can go to the next temperature


EOT



chmod +x ${NAME}.${t}.ConcatAndBar.sh


if [[ $ARCHITECTURE == "slurm" ]]; then
    sbatch ${NAME}.${t}.ConcatAndBar.sh
    echo "analysis job was sent. it will wait until dependencies finish"


elif [[ $ARCHITECTURE == "rome" ]]; then
    ccc_msub ${NAME}.${t}.ConcatAndBar.sh
    echo "analysis job was sent. it will wait until dependencies finish"


elif [[ $ARCHITECTURE == "pc" ]]; then
    nohup ${NAME}.${t}.ConcatAndBar.sh 2>&1 &

    echo "analysis job was sent. it will wait 4 hours before starting"

fi





done # temperature loop
