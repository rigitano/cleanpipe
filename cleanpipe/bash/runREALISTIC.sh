#!/bin/bash

# usage example:  ./runREALISTIC.sh a12HW a12HW.gro a12HW.top 1 298 charmm36 pc 2 8

# ARGUMENTS:
# 1-name that goes onthe runREALISTIC to be created and the job name
# 2-gro
# 3-top
# 4-nanoseconds of production
# 5-temperature (remeber that martini3 was parametrized at 310)
# 6-forcefield to be used in mdp construction (must be "charmm36" or "martini3")
# 7-architecture (pc, slurm, rome)
# 8-ntOMP
# 9-ntMPI

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

if [[ $ARCHITECTURE == "pc" || $ARCHITECTURE == "slurm" || $ARCHITECTURE == "rome" ]]; then
    echo "Architecture is valid"
else
    echo "Error: architecture must be 'pc' or 'slurm' or 'rome' "
    exit 1
fi
echo " "


NTOMP=$8
if [[ "$NTOMP" =~ ^[0-9]+$ ]]; then
    echo "NTOMP OK (is an integer)"
else
    echo "NTOMP is NOT an integer"
fi



NTMPI=$9
if [[ "$NTMPI" =~ ^[0-9]+$ ]]; then
    echo "NTMPI OK (is an integer)"
else
    echo "NTMPI is NOT an integer"
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

#mkdir -p runREALISTIC_${NAME}/0_SC
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
file_sc_mdp="sc.mdp"
file_em_mdp="em.mdp"
file_nvt_mdp="nvt.mdp"
file_npt_mdp="npt.mdp"
file_prod_mdp="prod.mdp"

echo "creating ${file_sc_mdp}"
cat <<EOT > "0_SC/${file_sc_mdp}"

; the original  came from
; https://github.com/jacksoncrowley/TS2CG-Setup-Pipeline/blob/main/mdp/em1.mdp
; but the important part is of the free energy variables. mainly setting the forces as 1% and putting some sc potential. jackson has no strong feeling about the rest of the parameters.




define			 = -DFLEXIBLE

integrator               = steep
nsteps                   = 500
nstxout                  = 0
nstfout                  = 0
nstlog                   = 100 

; NEIGHBORSEARCHING PARAMETERS
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
periodic-molecules       = no
verlet-buffer-tolerance  = 0.005
rlist                    = 1

; OPTIONS FOR ELECTROSTATICS AND VDW
coulombtype              = cut-off
coulomb-modifier         = Potential-shift-Verlet
rcoulomb-switch          = 0
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw-switch              = 0
rvdw                     = 1.1

; Free energy variables
free-energy = yes
init-lambda              = 0.01
sc-alpha                 = 4
sc-power                 = 2
sc-coul                  = yes
nstdhdl                  = 0 
couple-moltype           = system
; we are changing both the vdw and the charge. In the initial state, both are on
couple-lambda0           = vdw-q
; in the final state, both are off.
couple-lambda1           = none
couple-intramol          = yes



EOT




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
nstcomm                  = 50000

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
dt                       = $(options charmm36=0.001 martini3=0.02)
nsteps                   = 10000000
nstcomm                  = 50000

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
#MSUB   -c 1                       # Number of cores per parallel task
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

OMP_NUM_THREADS=${NTOMP}      # number of OpenMP threads

EOT

GMX="ccc_mprun gmx_mpi"
#MDRUN_OPTIONS="-dd 3 3 3 -npme 13 -dlb yes" #ATENTION: this must be coherent with #MSUB -n 40 #domain decomposition dont work with steep 
#MDRUN_OPTIONS="-nt 1" #nt cant be used in rome, just set OMP_NUM_THREADS and #MSUB -n
MDRUN_OPTIONS=""




##########################################################################################################
elif [[ $ARCHITECTURE == "slurm" ]]; then # insert the slurm header, if the user chose this architecture

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

EOT

GMX="gmx"
MDRUN_OPTIONS="-ntomp ${NTOMP} -ntmpi ${NTMPI}" #ATENTION: this must be coherent with #SBATCH --cpus-per-task







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

set -o pipefail  # stop if any part of a pipeline fails



#echo "#############################################################"
#echo "######################### SC grompp #########################"
#echo "#############################################################"
#module unload gromacs
#module load gromacs/2023 # this is an ugly fix because in newer versions the SC run will give an arror related to decupling and coulomb. But here that error, whatever it means, is irrelevant. Im just trying to get rid of superpositions.
#
#cd 0_SC || exit
#       
#${GMX} grompp -f ${file_sc_mdp} -c "../../${GRO}" -p "../../${TOP}" -o "sc.tpr" 2>&1 | tee "outanderr.grompp"
#
#if grep -q "Error" "outanderr.grompp"; then
#    echo "GROMACS reported an error — stopping script."
#    exit 1
#fi
#
#
#
#echo "########################## SC mdrun #############################"
#
#
#${GMX} mdrun -v -deffnm "sc" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"
#
#if grep -q "Error" "outanderr.mdrun"; then
#    echo "GROMACS reported an error — stopping script."
#    exit 1
#fi
#
#
#
#module unload gromacs
#module load gromacs/2024.5


echo "#############################################################"
echo "######################### EM grompp #########################"
echo "#############################################################"

#cd ../1_EM || exit
cd 1_EM || exit

#${GMX} grompp -f ${file_em_mdp} -c "../0_SC/sc.gro" -p "../../${TOP}" -o "em.tpr" 2>&1 | tee "outanderr.grompp"	
${GMX} grompp -f ${file_em_mdp} -c "../../${GRO}" -p "../../${TOP}" -o "em.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"

if grep -q "Error" "outanderr.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



echo "######################### EM mdrun #############################"

${GMX} mdrun -deffnm "em" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"

if grep -q "Error" "outanderr.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi

echo "#############################################################"
echo "######################### NVT grompp #########################"
echo "#############################################################"

cd ../2_NVT || exit

${GMX} grompp -f ${file_nvt_mdp} -c "../1_EM/em.gro" -r "../1_EM/em.gro" -p "../../${TOP}" -o "nvt.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"

if grep -q "Error" "outanderr.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



echo "######################### NVT mdrun ##############################"

${GMX} mdrun -v -deffnm "nvt" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"

if grep -q "Error" "outanderr.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "##############################################################"
echo "######################### NPT grompp #########################"
echo "##############################################################"

cd ../3_NPT || exit

${GMX} grompp -f ${file_npt_mdp} -c "../2_NVT/nvt.gro" -r "../2_NVT/nvt.gro" -p "../../${TOP}" -o "npt.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"

if grep -q "Error" "outanderr.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



echo "######################### NPT mdrun ##############################"

${GMX} mdrun -v -deffnm "npt" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"

if grep -q "Error" "outanderr.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if the intent of the script is realy just make REALISTIC, you could quit here
# thats why I usually put time 1 ns. by doing this the next step is as short as possible

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





echo "#############################################################"
echo "#################### PRODUCTION grompp #######################"
echo "#############################################################"

cd ../4_PROD || exit
	
${GMX} grompp -f ${file_prod_mdp} -c "../3_NPT/npt.gro" -p "../../${TOP}" -o "prod.tpr" -maxwarn 1 2>&1 | tee "outanderr.grompp"

if grep -q "Error" "outanderr.grompp"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


echo "###################### PRODUCTION mdrun #################################"

${GMX} mdrun -v -deffnm "prod" ${MDRUN_OPTIONS} 2>&1 | tee "outanderr.mdrun"

if grep -q "Error" "outanderr.mdrun"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi




echo "############################## CENTER AND FIT ###################################"


printf '1\n0' | ${GMX} trjconv -s "prod.tpr" -f "prod.xtc" -o "prod.centered.xtc" -center -pbc mol 2>&1 | tee "log.center"

if grep -q "Error" "log.center"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi



printf '1\n0' | ${GMX} trjconv -s "prod.tpr" -f "prod.centered.xtc" -o "prod.fitted.xtc" -fit progressive 2>&1 | tee "log.fit"

if grep -q "Error" "log.fit"; then
    echo "GROMACS reported an error — stopping script."
    exit 1
fi


cd ../.. # get out of runREALISTIC





EOT
	
chmod +x script.${NAME}.sh


if [[ $ARCHITECTURE == "slurm" ]]; then
    sbatch script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "rome" ]]; then
    ccc_msub script.${NAME}.sh && echo "job was sent"

elif [[ $ARCHITECTURE == "pc" ]]; then
    ./script.${NAME}.sh && echo "script finished"
    

fi










