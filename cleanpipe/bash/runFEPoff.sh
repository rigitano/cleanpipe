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
# 7-architecture (pc, oxygen, rome, genoa) # rome is a tgcc partition, genoa is an adastra partiion
# 8-molecule to be decoupled
# 9-ntOMP     #igonred in rome, because -dd works great!
# 10-ntMPI    #igonred in rome, because -dd works great!

if [ $# -lt 10 ]; then
    echo "10 arguments needed : name filename.gro filename.top numberOfNanoseconds temperatureList forceFieldName architecture moleculeName ntOMP ntMPI"
    exit 1
fi






######## obtain the arguments defined by the user when he called the function ########
NAME=$1 #name of the system to be simulated
echo " "
echo "Name of system: ${NAME}"
echo " "

#GRO=$2 #name of gro
GRO=$(readlink -f "$2")
#TOP=$3 #name of top
TOP=$(readlink -f "$3")
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
if [[ $ARCHITECTURE == "pc" || $ARCHITECTURE == "oxygen" || $ARCHITECTURE == "rome" || $ARCHITECTURE == "genoa" ]]; then
    echo "Architecture is valid"
else
    echo "Error: architecture must be 'pc' 'oxygen' 'rome' 'genoa' 'MI300' "
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











echo "creating folder structure..."


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

echo "generating mdp files..."
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



echo "${file_em_mdp}"
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
nstlog                   = 5000
nstenergy                = 5000

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




echo "${file_nvt_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/2_NVT/${file_nvt_mdp}"

; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)

nsteps                   = 50000
nstcomm                  = 100
comm_mode                = linear
comm_grps                = 

; Output control
nstxout                  = 5000
nstvout                  = 5000
nstfout                  = 0
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



; how to restraint protein position, and how to make water flexible
; define = -DPOSRES -DFLEXIBLE
EOT




echo "${file_npt_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/3_NPT/${file_npt_mdp}"


; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)
nsteps                   = 100000
nstcomm                  = 100
comm_mode                = linear
comm_grps                = 

; Output control
nstxout                  = 5000
nstvout                  = 5000
nstfout                  = 0
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



; how to restraint protein position, and how to make water flexible
; define = -DPOSRES -DFLEXIBLE

EOT





echo "${file_prod_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/4_PROD/${file_prod_mdp}"


; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = $(options charmm36=0.002 martini3=0.02)
nsteps                   = ${STEPS}
nstcomm                  = 100
comm_mode                = linear
comm_grps                = 

; Output control
nstxout                  = 50000
nstvout                  = 50000
nstfout                  = 50000
nstlog                   = 50000
nstenergy                = 5000
nstxout-compressed       = 50000
nstdhdl                  = 100
nstcalcenergy            = 100

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




; how to restraint protein position, and how to make water flexible
; define = -DPOSRES -DFLEXIBLE

EOT





done #lambda loop
done #temperature loop

echo "all mdp files created"
echo " "



echo "creating simulation scripts for all lambdas. repeating for each temperature..."
for t in $TEMPERATURE_LIST; do
for i in $(seq -w 0 20); do 

# get the name of the mdp files for the current temperature and lambda
file_em_mdp="em_FEP_t${t}_lambda${i}.mdp"
file_nvt_mdp="nvt_FEP_t${t}_lambda${i}.mdp"
file_npt_mdp="npt_FEP_t${t}_lambda${i}.mdp"
file_prod_mdp="prod_FEP_t${t}_lambda${i}.mdp"



cd t${t}/Lambda_${i}
pwd

cat <<EOT > "t${t}.l${i}.sh"
#!/bin/bash

EOT




########################################################################################################
if [[ $ARCHITECTURE == "rome" ]]; then # insert the tgcc-rome header, if the user chose this architecture


cat <<EOT >>  "t${t}.l${i}.sh"

#MSUB   -r ${NAME}.${t}                   # Job name
#MSUB   -n ${NTMPI}                       # Number of tasks in parallel mode
#MSUB   -c ${NTOMP}                       # Number of cores per parallel task
#MSUB   -W                                # Let multiple jobs sharing same name & user run simultaneously
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
module load gromacs/2025.4 # charge le produit

export GMX_DISABLE_GPU_DETECTION=1 # prevent GROMACS from using GPUs
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto

export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)


# ---- 24h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain. I'll know this checking for a done.txt file, that is created in the post processing
if [[ -f "4_PROD/done.txt" ]]; then
    echo "Simulation already complete. No need for follow-up job"
else
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
ccc_msub -E "--dependency=afterok:\${BRIDGE_MSUB_JOBID}" t${t}.l${i}.sh
fi
# --------------------------------------------------------------



EOT

#GMX ENGINE 
GMX="ccc_mprun gmx_mpi"
GMXS="ccc_mprun -n 1 gmx_mpi"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #MSUB -n
#MDRUN_OPTIONS="-dd 3 3 3 -npme 13 -dlb yes" #this requires #MSUB -n 40 , but domain decomposition dont work with steep 
#MDRUN_OPTIONS="-nt 1" #doesbt work, because -nt, -ntomp, -ntmpi, cant be used in rome, you have to set OMP_NUM_THREADS and #MSUB -n instead
#MDRUN_OPTIONS=""
MDRUN_OPTIONS="-maxh 20 -cpi" 


#########################################################################################################
elif [[ $ARCHITECTURE == "genoa" ]]; then # insert the adastra-genoa header, if the user chose this architecture

cat <<EOT >>  "t${t}.l${i}.sh"

#SBATCH --account=c1613458
#SBATCH -J ${NAME}.${t}
#SBATCH --constraint=GENOA         # GENOA(192)(CPU) or MI250(64)(GPU)
##SBATCH --nodes=
#SBATCH --ntasks-per-node=${NTMPI} 
#SBATCH --cpus-per-task=${NTOMP}
##SBATCH --exclusive
#SBATCH -o t${t}.l${i}.%j.scheduler.out
#SBATCH -e t${t}.l${i}.%j.scheduler.err 


module purge
#module load CCE-CPU-4.0.0
develop CCE-CPU-5.0.0
module spider gromacs/2025.2-omp-mpi
module load gromacs/2025.2-omp-mpi

module list


export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)


# ---- 48h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain. I'll know this checking for a done.txt file, that is created in the post processing
if [[ -f "4_PROD/done.txt" ]]; then
    echo "Simulation already complete. No need for follow-up job"
    exit 0
else
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
sbatch --dependency=afterok:\$SLURM_JOB_ID t${t}.l${i}.sh

fi
# --------------------------------------------------------------



EOT

#GMX ENGINE 
GMX="srun gmx_mpi"
GMXS="srun -n 1 gmx_mpi"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #SBATCH
MDRUN_OPTIONS="-maxh 20 -cpi"





#########################################################################################################
elif [[ $ARCHITECTURE == "oxygen" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >>  "t${t}.l${i}.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=${NTOMP}
#SBATCH --ntasks-per-node=${NTMPI}
#SBATCH --gres=gpu:1
#SBATCH --job-name=${NAME}.${t}
#SBATCH --output=t${t}.l${i}.%j.scheduler.outanderr
#SBATCH --exclude=node-15

module purge
module load gromacs/2025.4

# ---- 48h-wall self-chaining : queue the follow-up job now ----
# If production is already finished, stop the chain. I'll know this checking for a done.txt file, that is created in the post processing
if [[ -f "4_PROD/done.txt" ]]; then
    echo "Simulation already complete. No need for follow-up job"
    exit 0
else
# Otherwise queue the NEXT copy of this job, to start when THIS one ends OK.
sbatch --dependency=afterok:\$SLURM_JOB_ID t${t}.l${i}.sh

fi
# --------------------------------------------------------------

EOT

#GMX ENGINE
GMX="gmx"
GMXS="gmx"

#GMX MDRUN ADITIONAL OPTIONS. ATENTION: this must be coherent with #SBATCH --cpus-per-task
#MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI}" #this requires  #SBATCH --cpus-per-task=1
MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI} -maxh 44 -cpi" #ATENTION: -ntomp and -ntmpi must be coherent with #SBATCH --cpus-per-task




##########################################################################################################3
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >>  "t${t}.l${i}.sh"

module purge
module load gromacs/2025.4



EOT

#GMX ENGINE
GMX="gmx"
GMXS="gmx"

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

##### function to check if a simulation reached the planned number of steps #####
planned_steps_reached() {
    local TPR="\$1" CPT="\$2"                                                        # the inputs are the TPR filename, and checkpoint filename.

    [[ -s "\$TPR" && -s "\$CPT" ]] || return 1                                       # Return false if one of the input files is missing or empty.

    local CURRENT_STEP PLANNED_STEPS
    CURRENT_STEP=\$(${GMXS} dump -cp "\$CPT" 2>/dev/null | awk -F= '/^[[:space:]]*step[[:space:]]*=/{gsub(/[[:space:]]/,"",\$2); print \$2; exit}')
    PLANNED_STEPS=\$(${GMXS} dump -s "\$TPR" 2>/dev/null | awk -F= '/^[[:space:]]*nsteps[[:space:]]*=/{gsub(/[[:space:]]/,"",\$2); print \$2; exit}')

    [[ "\$CURRENT_STEP" =~ ^[0-9]+$ && "\$PLANNED_STEPS" =~ ^[0-9]+$ ]] || return 1  # Return false if one of the value is empty or negative 
    (( CURRENT_STEP >= PLANNED_STEPS ))                                              # Return true if the planned number of steps has been reached.
}

##### function to inspect gmx outanderr file, looking for failure messages #####
gmx_failed() {
    local LABEL="\$1"
    local OUTANDERR_FILE="\$2"
    local PATTERN

    PATTERN='fatal[[:space:]]+error|error[[:space:]]+in[[:space:]]+user[[:space:]]+input|ERROR[[:space:]]+[1-9][0-9]*|there (was|were) [1-9][0-9]* errors?|inconsistency[[:space:]]+in[[:space:]]+user[[:space:]]+input|too[[:space:]]+many[[:space:]]+warnings|failure|assertion[[:space:]]+failed|segmentation[[:space:]]+fault|floating[[:space:]]+point[[:space:]]+exception|bus[[:space:]]+error|core[[:space:]]+dumped|aborted|killed|out[[:space:]]+of[[:space:]]+memory|cannot[[:space:]]+allocate[[:space:]]+memory|permission[[:space:]]+denied|no[[:space:]]+such[[:space:]]+file|cannot[[:space:]]+open|could[[:space:]]+not[[:space:]]+be[[:space:]]+opened|command[[:space:]]+not[[:space:]]+found'


    if LC_ALL=C grep -Eiq "\$PATTERN" "\$OUTANDERR_FILE"; then
        echo "error during \$LABEL. this is reported in \$OUTANDERR_FILE — stopping script."
        return 0
    fi

    return 1
}



echo "############## EM - Temperature ${t} Lambda ${i} #########################################"
cd 1_EM || exit 1
	
if [[ -s "${NAME}_em.gro" ]]; then
    echo "skipping em (${NAME}_em.gro already there)"
else

    ${GMXS} grompp -f ${file_em_mdp} -c "${GRO}" -p "${TOP}" -o "${NAME}_em.tpr" 2>&1 | tee "outanderr.grompp"
    if gmx_failed "${NAME}_em" "outanderr.grompp"; then exit 1; fi
    if [[ ! -f "${NAME}_em.tpr" ]]; then echo "${NAME}_em finished without saving a TPR file"; exit 1; fi


    echo "##########################################################################################"

    ${GMX} mdrun -deffnm "${NAME}_em" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
    if gmx_failed "${NAME}_em" "outanderr.mdrun"; then exit 1; fi
    if [[ ! -f "${NAME}_em.gro" ]]; then echo "${NAME}_em finished without saving a GRO file"; exit 1; fi




fi
echo "############# NVT - Temperature ${t} Lambda ${i}  ########################################"
cd ../2_NVT || exit 1
if planned_steps_reached "${NAME}_nvt.tpr" "${NAME}_nvt.cpt"; then
    echo "skipping ${NAME}_nvt (steps reached)"
else


    ${GMXS} grompp -f ${file_nvt_mdp} -c "../1_EM/${NAME}_em.gro" -r "../1_EM/${NAME}_em.gro" -p "${TOP}" -o "${NAME}_nvt.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"
    if gmx_failed "${NAME}_nvt" "outanderr.grompp"; then exit 1; fi
    if [[ ! -f "${NAME}_nvt.tpr" ]]; then echo "${NAME}_nvt finished without saving a TPR file"; exit 1; fi


    echo "##########################################################################################"

    ${GMX} mdrun -v -deffnm "${NAME}_nvt" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
    if gmx_failed "${NAME}_nvt" "outanderr.mdrun"; then exit 1; fi
    if ! planned_steps_reached "${NAME}_nvt.tpr" "${NAME}_nvt.cpt"; then echo "NVT stopped before completion; job must resume it."; exit 0; fi
    if [[ ! -f "${NAME}_nvt.gro" ]]; then echo "${NAME}_nvt finished without saving a GRO file"; exit 1; fi



fi
echo "##################### NPT - Temperature ${t} Lambda ${i}  ###############################"
cd ../3_NPT || exit 1

if planned_steps_reached "${NAME}_npt.tpr" "${NAME}_npt.cpt"; then
    echo "skipping ${NAME}_npt (steps reached)"
else

    ${GMXS} grompp -f ${file_npt_mdp} -c "../2_NVT/${NAME}_nvt.gro" -r "../2_NVT/${NAME}_nvt.gro" -p "${TOP}" -o "${NAME}_npt.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"
    if gmx_failed "${NAME}_npt" "outanderr.grompp"; then exit 1; fi
    if [[ ! -f "${NAME}_npt.tpr" ]]; then echo "${NAME}_npt finished without saving a TPR file"; exit 1; fi


    echo "#########################################################################################"

    ${GMX} mdrun -v -deffnm "${NAME}_npt" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
    if gmx_failed "${NAME}_npt" "outanderr.mdrun"; then exit 1; fi
    if ! planned_steps_reached "${NAME}_npt.tpr" "${NAME}_npt.cpt"; then echo "NVT stopped before completion; job must resume it."; exit 0; fi
    if [[ ! -f "${NAME}_npt.gro" ]]; then echo "${NAME}_npt finished without saving a GRO file"; exit 1; fi



fi
echo "#################### PRODUCTION - Temperature ${t} Lambda ${i}  #######################"
cd ../4_PROD || exit 1	
if planned_steps_reached "${NAME}_prod_${t}_${i}.tpr" "${NAME}_prod_${t}_${i}.cpt"; then
    echo "skipping ${NAME}_prod_${t}_${i} (steps reached)"
else


    ${GMXS} grompp -f ${file_prod_mdp} -c "../3_NPT/${NAME}_npt.gro" -p "${TOP}" -o "${NAME}_prod_${t}_${i}.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"
    if gmx_failed "${NAME}_prod_${t}_${i}" "outanderr.grompp"; then exit 1; fi
    if [[ ! -f "${NAME}_prod_${t}_${i}.tpr" ]]; then echo "${NAME}_prod_${t}_${i} finished without saving a TPR file"; exit 1; fi



    echo "#######################################################################################"

    ${GMX} mdrun -v -deffnm "${NAME}_prod_${t}_${i}" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun" 
    if gmx_failed "${NAME}_prod_${t}_${i}" "outanderr.mdrun"; then exit 1; fi
    if ! planned_steps_reached "${NAME}_prod_${t}_${i}.tpr" "${NAME}_prod_${t}_${i}.cpt"; then echo "NVT stopped before completion; job must resume it."; exit 0; fi
    if [[ ! -f "${NAME}_prod_${t}_${i}.gro" ]]; then echo "${NAME}_prod_${t}_${i} finished without saving a GRO file"; exit 1; fi



fi
echo "############################## LAUNCH BAR CALCULATION  ###################################"
if [[ ! -f done.txt ]] && planned_steps_reached "${NAME}_prod_${t}_${i}.tpr" "${NAME}_prod_${t}_${i}.cpt"; then # only post-process once PROD has truly finished
    touch "done.txt"


    
    # this lambda is done. lets check if all the other lambdas are also done. if so, lets louch the bar calculation script.sh
    # notice that this test will be done for all lambdas, but only the last one to finish will enter the condition
    
    n=\$(find ../.. -path '*/4_PROD/done.txt' | wc -l)
    if (( n == 21 )); then
     
        if [[ ${ARCHITECTURE} == "rome" ]]; then
            cd ../../
            ccc_msub ${NAME}.${t}.BarForAllLambdas.sh
        elif [[ ${ARCHITECTURE} == "genoa" ]]; then
            cd ../../
            sbatch ${NAME}.${t}.BarForAllLambdas.sh
        elif [[ ${ARCHITECTURE} == "oxygen" ]]; then
            cd ../../
            sbatch ${NAME}.${t}.BarForAllLambdas.sh
        elif [[ ${ARCHITECTURE} == "pc" ]]; then
            cd ../../
            nohup ./${NAME}.${t}.BarForAllLambdas.sh > outanderror.BarForAllLambdas 2>&1 &
        fi
    fi


    echo "############################## CENTER AND FIT - Temperature ${t} Lambda ${i}  ###################################"


    printf '1\n0\n' | ${GMXS} trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.xtc" -o "${NAME}_prod_${t}_${i}.centered.xtc" -center -pbc mol 2>&1 | tee "outanderr.center"
    if gmx_failed "${NAME}_center" "outanderr.center"; then exit 1; fi



    printf '1\n0\n' | ${GMXS} trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.centered.xtc" -o "${NAME}_prod_${t}_${i}.fitted.xtc" -fit progressive 2>&1 | tee "outanderr.fit"
    if gmx_failed "${NAME}_fit" "outanderr.fit"; then exit 1; fi





fi








EOT
chmod +x t${t}.l${i}.sh





cd ../.. #go back to runFEP so we can go to the next temperature and lambda in the iteration
done # lambda loop
done # temperature loop

echo "simulation scripts were created for all lambdas and for all temperatures"
echo " "










echo "creating bar analysis scripts. 1 for each temperature..."
# analysis script that will be run just after all the jobs are finished. each temperature will have one

for t in $TEMPERATURE_LIST; do

cd t${t}
pwd

cat <<EOT > "${NAME}.${t}.BarForAllLambdas.sh"
#!/bin/bash

EOT




########################################################################################################
if [[ $ARCHITECTURE == "rome" ]]; then # insert the rome header, if the user chose this architecture


cat <<EOT >>  "${NAME}.${t}.BarForAllLambdas.sh"

#MSUB   -r ${NAME}.${t}                      # Job name
#MSUB   -n 1                                 # Number of tasks in parallel mode (ntmpi)
#MSUB   -c 1                                 # Number of cores per parallel task
#MSUB   -w                                   # PREVENT multiple jobs sharing same name run simultaneously. so this job will run only when the others are done
#MSUB   -o ${t}.barForAllLambdas.%I.scheduler.out        # Output file
#MSUB   -e ${t}.barForAllLambdas.%I.scheduler.err        # Output file for errors
#MSUB   -q rome                              # Partition:    rome        
#MSUB   -A gen13458                          # Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store                # File system:  scratch,work,store
#MSUB   -Q normal                            # Quality of Service (test,normal,long) (ccc_mqinfo)
#MSUB   -T 3600                              # Maximum walltime in seconds
#MSUB   -E "--dependency=singleton"          # probably the jobname will be enought to set the  
#MSUB   -@ henrique.rigitano@ibcp.fr:end

set -x # echo commands

module purge  # retire tous les modules déchargeables de l'environnement
module load gnu/11 # charge gnu/11 et définit gnu/11 comme compilateur dans votre environnement
module load nvhpc/24.3 # besoin de mettre avant OpenMPI comme ce dernier charge un cuda qui n'est pas compatible avec nvhpc/24.3
module load mpi/openmpi/4 # charge la souche OpenMPI
module load gromacs/2025.4 # charge le produit

export GMX_DISABLE_GPU_DETECTION=1 # prevent GROMACS from using GPUs
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto

export OMP_NUM_THREADS=1      # number of OpenMP threads (ntomp)

EOT

#GMX ENGINE 
GMX="ccc_mprun gmx_mpi"

#########################################################################################################
elif [[ $ARCHITECTURE == "genoa" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >>  "${NAME}.${t}.BarForAllLambdas.sh"

#SBATCH --account=c1613458
#SBATCH -J ${NAME}.${t}
#SBATCH --constraint=GENOA         # GENOA(192)(CPU) or MI250(64)(GPU)
##SBATCH --nodes=
#SBATCH --ntasks-per-node=${NTMPI} 
#SBATCH --cpus-per-task=${NTOMP}
##SBATCH --exclusive
#SBATCH -o ${t}.barForAllLambdas.%j.scheduler.out
#SBATCH -e ${t}.barForAllLambdas.%j.scheduler.err 


module purge
#module load CCE-CPU-4.0.0
develop CCE-CPU-5.0.0
module spider gromacs/2025.2-omp-mpi
module load gromacs/2025.2-omp-mpi

module list


export OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads (ntomp)



EOT

#GMX ENGINE
GMX="srun gmx_mpi"

#########################################################################################################
elif [[ $ARCHITECTURE == "oxygen" ]]; then # insert the slurm header, if the user chose this architecture

cat <<EOT >>  "${NAME}.${t}.BarForAllLambdas.sh"

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=1
#SBATCH --job-name=${NAME}.${t}
#SBATCH --output=final.analysis.%j.scheduler.outanderr
#SBATCH --dependency=singleton

module purge
module load gromacs/2025.4



EOT

#GMX ENGINE
GMX="gmx"


##########################################################################################################3
elif [[ $ARCHITECTURE == "pc" ]]; then # insert what should be the gromacs commands in my local pc

cat <<EOT >>  "${NAME}.${t}.BarForAllLambdas.sh"


module purge
module load gromacs/2025.4



EOT

#GMX ENGINE
GMX="gmx"



########################################################################################################
fi # end of if that inserts script headers and module loading before the gromacs commands



#now the gromacs commands will be appended to the headers and moldule loading
cat <<EOT >> "${NAME}.${t}.BarForAllLambdas.sh"

set -o pipefail  # stop if any part of a pipeline fails

##### function to inspect gmx outanderr file, looking for failure messages #####
gmx_failed() {
    local LABEL="\$1"
    local OUTANDERR_FILE="\$2"
    local PATTERN

    PATTERN='fatal[[:space:]]+error|error[[:space:]]+in[[:space:]]+user[[:space:]]+input|ERROR[[:space:]]+[1-9][0-9]*|there (was|were) [1-9][0-9]* errors?|inconsistency[[:space:]]+in[[:space:]]+user[[:space:]]+input|too[[:space:]]+many[[:space:]]+warnings|failure|assertion[[:space:]]+failed|segmentation[[:space:]]+fault|floating[[:space:]]+point[[:space:]]+exception|bus[[:space:]]+error|core[[:space:]]+dumped|aborted|killed|out[[:space:]]+of[[:space:]]+memory|cannot[[:space:]]+allocate[[:space:]]+memory|permission[[:space:]]+denied|no[[:space:]]+such[[:space:]]+file|cannot[[:space:]]+open|could[[:space:]]+not[[:space:]]+be[[:space:]]+opened|command[[:space:]]+not[[:space:]]+found'


    if LC_ALL=C grep -Eiq "\$PATTERN" "\$OUTANDERR_FILE"; then
        echo "error during \$LABEL. this is reported in \$OUTANDERR_FILE — stopping script."
        return 0
    fi

    return 1
}





for i in {00..20}; do  # lambda loop to see if all



    cd Lambda_\${i}/4_PROD || exit 1

    #check if done file is there for corrent lambda. if its not, it will exit the script
    if [[ -f "done.txt" ]]; then
        cd ../..  #go back to runFEP, the next iteration will jump to the correct 4_PROD
    else
        echo "Missing done.txt t\${t} Lambda_\${i}"
        exit 1
    fi




done # lambda loop





echo '######################################################################'
echo '##################### bar and barint calculation #####################'
echo '############ (using all lambdas of a given temperature) ##############'
echo '######################################################################'


echo "${NAME}.${t}.BarForAllLambdas.sh script initiated"


#create directory if it doesnt exist. if exists exit, becaus another cript already did the calculation
mkdir bar 2>/dev/null || exit 0


${GMX} bar -f Lambda_*/4_PROD/${NAME}_prod_${t}_*.xvg -o bar/${NAME}_${t}_bar.xvg -oi bar/${NAME}_${t}_barint.xvg 2>&1 | tee "outanderr.bar"
if gmx_failed "bar" "outanderr.bar"; then exit 1; fi







EOT
chmod +x ${NAME}.${t}.BarForAllLambdas.sh





cd .. #go back to runFEP so we can go to the next temperature
done # temperature loop
echo "all bar calculation scripts were created"
echo " "


echo "louching the simulation scripts for all lambdas. repeating for each temperature..." #. Remember that they contain a chunk of code that check when all are finished for a given temperature, so that the last one will lounch the bar calculation for that temperature"
for t in $TEMPERATURE_LIST; do
for i in $(seq -w 0 20); do 


cd t${t}/Lambda_${i}
pwd



if [[ $ARCHITECTURE == "oxygen" ]]; then
    jid=$(sbatch --parsable "t${t}.l${i}.sh")
    jid=${jid%%;*}          # federated clusters return jobid;cluster

    echo "job was sent - id ${jid}"

elif [[ $ARCHITECTURE == "genoa" ]]; then
    jid=$(sbatch --parsable "t${t}.l${i}.sh")
    jid=${jid%%;*}          # federated clusters return jobid;cluster

    echo "job was sent - id ${jid}"

elif [[ $ARCHITECTURE == "rome" ]]; then
    jid=$(ccc_msub t${t}.l${i}.sh | grep -Eo '[0-9]+' | tail -n1)
    #build list of ids sent, to use to set the dependency of the analysis script 

    echo "job was sent - id ${jid}"


elif [[ $ARCHITECTURE == "pc" ]]; then
    nohup ./t${t}.l${i}.sh > t${t}.l${i}.outanderr 2>&1 &
    pid=$!

    echo "script was lounched - pid ${pid}"

fi

cd ../.. 
done # lambda loop over
done # temperature loop over
echo "scripts were louched for all lambdas and temperatures"
echo " "
