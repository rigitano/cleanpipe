#!/bin/bash

# usage example:  ./runFEPoff-rome.sh insulin insulin.gro insulin.top 100

if [ $# -lt 4 ]; then
    echo "4 arguments needed :          name       filename.gro   filename.top   numberOfNanoseconds"
    exit 1
fi





#the simulation pipeline (EM,NVT, NPT, PROD) will have to be repeated 40 times for the same system,
#so fist a huge file structure is created. what I mean is that 40 folders are created inside the _runFEB folder, each one with the 4 steps of the pipeline
#then all the mdp files are created (what changes between the lambdas is just the paramenter 'init_lambda_state' )





######## obtain the arguments defined by the user when he called the function ########
NAME=$1 #name of the thing to be modeled
echo " "
echo "Name of thing to be modeled: ${NAME} (That name will be used to create a folder with that name folowed by '_run' and will also be the prefix of the individual files to be created )"
echo " "

GRO=$2 #name of gro
TOP=$3 #name of top
echo "GRO and TOP file names: ${GRO} ${TOP} (Those are inputs for the gromacs modelisation pipeline)"
echo " "

PRODUCTION_DURATION=$4 #how many nanoseconds
STEPS=$(( (PRODUCTION_DURATION * 1000000) / 2 ))
echo "Duration of production MD : ${PRODUCTION_DURATION} ns"
echo "Steps for that duration : ${STEPS} (time x 1000 / 0.002)"



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


mkdir "${NAME}__runFEP"
cd "${NAME}__runFEP" || exit

for t in 283 298 313; do

    mkdir "t${t}"
    cd "t${t}" || exit



    for i in $(seq -w 0 20); do
	mkdir "Lambda_${i}"
	cd "Lambda_${i}"  || exit

	mkdir 1_EM
	mkdir 2_NVT
	mkdir 3_NPT
	mkdir 4_PROD
	
	cd .. #back to T${t} folder
    done # lambda loop

    cd .. #back to __runFEP folder
done # temperature loop

echo "folder structure created"


######################## generate mdp files #############################



for t in 283 298 313; do
for i in $(seq -w 0 20); do # from 00 to 20. this impacts the mdp pararamenter init_lambda_state


# Definition of the name of the mdp files of the current lambda
file_em_mdp="em_FEP_t${t}_lambda${i}.mdp"
file_nvt_mdp="nvt_FEP_t${t}_lambda${i}.mdp"
file_npt_mdp="npt_FEP_t${t}_lambda${i}.mdp"
file_prod_mdp="prod_FEP_t${t}_lambda${i}.mdp"
file_bench_mdp="bench_FEP_t${t}_lambda${i}.mdp"

echo "creating ${file_em_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/1_EM/${file_em_mdp}"

; Run control
integrator               = steep 
nsteps                   = 5000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
;niter                    = 20
;nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2


; non-bonded electrostatic forces
rcoulomb            	=     1.2                ; ideal potential up to this radius
coulombtype         	=     PME                ; what to do after that radius
pme_order               =     4                  ; PME PARAMETER: higher values improve the accuracy
fourierspacing          =     0.12               ; PME PARAMETER: lower values improve the accuracy
ewald_rtol              =     1e-05              ; PME PARAMETER: lower values improve the accuracy


; non-bonded Van Der Waals forces
rvdw                	=     1.2                ; ideal potential up to this radius
vdw_type             	=     cutoff             ; what to do after that radius
vdw-modifier         	=     force-switch       ; deal with that simple cutoff after the radius
rvdw-switch          	=     1.0                ; deal with that simple cutoff after the radius
DispCorr            	=     no                 ; deal with that simple cutoff after the radius


; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no

; Free energy control stuff
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1                ; only immediate neighboring windows
couple-moltype           = Protein_chain_A  ; name of moleculetype to decouple
couple-lambda0           = vdw-q            ; Van der Waals and Coulomb  
couple-lambda1           = none             ; turn off everything
couple-intramol          = yes              ; big molecule

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no     
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10


; basic setup
gen_vel                  = no 
continuation             = no

; bond restrictions
constraints              = none 


EOT




echo "creating ${file_nvt_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/2_NVT/${file_nvt_mdp}"

; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000
nstcomm                  = 50000

; Output control
nstxout                  = 50000
nstvout                  = 50000
nstfout                  = 0
nstlog                   = 50000
nstenergy                = 50000
nstxout-compressed       = 0

; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20 
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2

; non-bonded electrostatic forces
rcoulomb            	=     1.2                ; ideal potential up to this radius
coulombtype         	=     PME                ; what to do after that radius
pme_order               =     4                  ; PME PARAMETER: higher values improve the accuracy
fourierspacing          =     0.12               ; PME PARAMETER: lower values improve the accuracy
ewald_rtol              =     1e-05              ; PME PARAMETER: lower values improve the accuracy


; non-bonded Van Der Waals forces
rvdw                	=     1.2                ; ideal potential up to this radius
vdw_type             	=     cutoff             ; what to do after that radius
vdw-modifier         	=     force-switch       ; deal with that simple cutoff after the radius
rvdw-switch          	=     1.0                ; deal with that simple cutoff after the radius
DispCorr            	=     no                 ; deal with that simple cutoff after the radius


; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = ${t}

; Pressure coupling is off for NVT
Pcoupl                   = No


; Free energy control stuff
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = Protein_chain_A  ; name of moleculetype to decouple
couple-lambda0           = vdw-q    ; altered so to include Coulombic interactions  
couple-lambda1           = none     ; turn off everything
couple-intramol          = yes      ; big molecule

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no          
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10

; options for bonds
constraint-algorithm     = lincs
constraints              = h-bonds  ; we only have C-H bonds here
lincs-order              = 4

; general dynamic setup
continuation             = no 
gen_vel                  = yes
gen_temp                 = ${t}
gen_seed                 = -1



; restraints configuration
disre = ${distance_restraints_option}
disre_fc = 1000

dihre = ${dihedral_restraints_option}
dihre_fc = 1000
EOT




echo "creating ${file_npt_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/3_NPT/${file_npt_mdp}"


; Run control
define = -DFLEXIBLE
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 100000
nstcomm                  = 50000

; Output control
nstxout                  = 50000
nstvout                  = 50000
nstfout                  = 0
nstlog                   = 50000
nstenergy                = 50000
nstxout-compressed       = 0

; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2

; non-bonded electrostatic forces
rcoulomb            	=     1.2                ; ideal potential up to this radius
coulombtype         	=     PME                ; what to do after that radius
pme_order               =     4                  ; PME PARAMETER: higher values improve the accuracy
fourierspacing          =     0.12               ; PME PARAMETER: lower values improve the accuracy
ewald_rtol              =     1e-05              ; PME PARAMETER: lower values improve the accuracy


; non-bonded Van Der Waals forces
rvdw                	=     1.2                ; ideal potential up to this radius
vdw_type             	=     cutoff             ; what to do after that radius
vdw-modifier         	=     force-switch       ; deal with that simple cutoff after the radius
rvdw-switch          	=     1.0                ; deal with that simple cutoff after the radius
DispCorr            	=     no                 ; deal with that simple cutoff after the radius

; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 2.0
ref_t                    = ${t} 
; Pressure coupling is on for NPT
Pcoupl                   = C-rescale
tau_p                    = 5.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
refcoord_scaling         = com

; Free energy control stuff
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = Protein_chain_A  ; name of moleculetype to decouple
couple-lambda0           = vdw-q    ; altered so to include Coulombic interactions  
couple-lambda1           = none     ; turn off everything
couple-intramol          = yes      ; big molecule

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no        
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10



; options for bonds
constraint-algorithm     = lincs
constraints              = h-bonds  ; we only have C-H bonds here
lincs-order              = 4 

; general dynamic setup
continuation             = yes 
gen_vel                  = no 



; restraints configuration
disre = ${distance_restraints_option}
disre_fc = 1000

dihre = ${dihedral_restraints_option}
dihre_fc = 1000
EOT





echo "creating ${file_prod_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/4_PROD/${file_prod_mdp}"


; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = ${STEPS}
nstcomm                  = 100

; Output control
nstxout                  = 500000
nstvout                  = 500000
nstfout                  = 0
nstlog                   = 500000
nstenergy                = 5000
nstxout-compressed       = 500000

; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2

; non-bonded electrostatic forces
rcoulomb            	=     1.2                ; ideal potential up to this radius
coulombtype         	=     PME                ; what to do after that radius
pme_order               =     4                  ; PME PARAMETER: higher values improve the accuracy
fourierspacing          =     0.12               ; PME PARAMETER: lower values improve the accuracy
ewald_rtol              =     1e-05              ; PME PARAMETER: lower values improve the accuracy


; non-bonded Van Der Waals forces
rvdw                	=     1.2                ; ideal potential up to this radius
vdw_type             	=     cutoff             ; what to do after that radius
vdw-modifier         	=     force-switch       ; deal with that simple cutoff after the radius
rvdw-switch          	=     1.0                ; deal with that simple cutoff after the radius
DispCorr            	=     no                 ; deal with that simple cutoff after the radius



; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 2.0
ref_t                    = ${t} 
; Pressure coupling is on for NPT
Pcoupl                   = c-rescale 
tau_p                    = 5.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
refcoord_scaling         = com

; Free energy control stuff
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = Protein_chain_A  ; name of moleculetype to decouple
couple-lambda0           = vdw-q    ; altered so to include Coulombic interactions  
couple-lambda1           = none     ; turn off everything
couple-intramol          = yes      ; big molecule

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00


; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no     
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10

; options for bonds
constraint-algorithm     = lincs
constraints              = h-bonds  ; we only have C-H bonds here
lincs-order              = 4 


; general dynamic setup
continuation             = yes 
gen_vel                  = no 



; restraints configuration
disre = ${distance_restraints_option}
disre_fc = 1000

dihre = ${dihedral_restraints_option}
dihre_fc = 1000
EOT



echo "creating ${file_bench_mdp}"
cat <<EOT > "t${t}/Lambda_${i}/4_PROD/${file_bench_mdp}"


; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000 ; 1 ns, because this is just to perform benchmark
nstcomm                  = 100

; Output control
nstxout                  = 500000
nstvout                  = 500000
nstfout                  = 0
nstlog                   = 500000
nstenergy                = 5000
nstxout-compressed       = 500000

; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2

; non-bonded electrostatic forces
rcoulomb            	=     1.2                ; ideal potential up to this radius
coulombtype         	=     PME                ; what to do after that radius
pme_order               =     4                  ; PME PARAMETER: higher values improve the accuracy
fourierspacing          =     0.12               ; PME PARAMETER: lower values improve the accuracy
ewald_rtol              =     1e-05              ; PME PARAMETER: lower values improve the accuracy


; non-bonded Van Der Waals forces
rvdw                	=     1.2                ; ideal potential up to this radius
vdw_type             	=     cutoff             ; what to do after that radius
vdw-modifier         	=     force-switch       ; deal with that simple cutoff after the radius
rvdw-switch          	=     1.0                ; deal with that simple cutoff after the radius
DispCorr            	=     no                 ; deal with that simple cutoff after the radius



; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 2.0
ref_t                    = ${t} 
; Pressure coupling is on for NPT
Pcoupl                   = c-rescale 
tau_p                    = 5.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
refcoord_scaling         = com

; Free energy control stuff
free_energy              = yes
init_lambda_state        = ${i}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = Protein_chain_A  ; name of moleculetype to decouple
couple-lambda0           = vdw-q    ; altered so to include Coulombic interactions  
couple-lambda1           = none     ; turn off everything
couple-intramol          = yes      ; big molecule

; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20  
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00

; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no     
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10

; options for bonds
constraint-algorithm     = lincs
constraints              = h-bonds  ; we only have C-H bonds here
lincs-order              = 4 


; general dynamic setup
continuation             = yes 
gen_vel                  = no 


; restraints configuration
disre = ${distance_restraints_option}
disre_fc = 1000

dihre = ${dihedral_restraints_option}
dihre_fc = 1000
EOT



done # lambda loop
done #temperature loop

echo "all mdp files created"





for t in 283 298 313; do
cd "t${t}"
for i in $(seq -w 0 20); do 

TPR_WITHOUT_EXTENSION="${NAME}_prod_${t}_${i}"
                       

cat <<EOT > "job.moab"
#!/bin/bash

#  @: Editable Variable
#     Please specify the input according to your needs
#     The rest should be fine
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CLUSTER SETTINGS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#MSUB   -r ${TPR_WITHOUT_EXTENSION}-mdrun            # Job name
#MSUB   -n 40                       # Number of tasks in parallel mode
#MSUB   -c 1                       # Number of cores per parallel task
# #MSUB -N 1                       # Number of nodes to allocate (Inferred)
#MSUB   -W yes                     # Let multiple jobs sharing same name & user run simultaneously
#MSUB   -o GMX.job.IR.output.%I    # Output file
#MSUB   -e GMX.job.IR.outerr.%I    # Output file for errors
#MSUB   -q rome                    # Partition:    rome        
#MSUB   -A gen13458                # Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store      # File system:  scratch,work,store
#MSUB   -Q normal                  # Quality of Service (test,normal,long) (ccc_mqinfo)
#MSUB   -T 86400                    # Maximum walltime in seconds

#@@@@@@@@@@@@@@@@@@@@ Set variables for naming the files @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
GMX=gmx_mpi
# Rootname for gromacs files
ROOTNAME=$TPR_WITHOUT_EXTENSION
# Additional options to use in mdrun
MDRUN_OPT="-dd 3 3 3 -npme 13 -dlb yes"
# The file name of the submition file (this file)
BATCH_FNAME=job.moab #[IreneRome]
# Date of run
DATE=\$(date | awk '{print \$1 \$2 \$3 \$6}')
# Time to count for the post-MD checks (in seconds)
CHECK_DURATION=60
set -x #[IreneRome] echo launched commands
##################### load gromacs and set variables  ############################
############# in principle, no changes needed beyond this point ##################

ml purge
module load gnu/11 mpi/openmpi/4 gromacs/2023.2
# Note that sometimes it is advisable especially for small simulations to use more open MP
# threads and less MPI ranks. However, for large systems this appears to be the most
# efficent and reasonably fast setting.
# Note that in any case the number of OpenMP threads times the MPI ranks has to be equal
# the number of available CPUS or cores unless hyperthreading is used, which GROMACS can
# do.

### How many CPU were asked?
ncores=\${SLURM_NTASKS} # number of MPI ranks = number of cores
OMP_NUM_THREADS=1      # number of OpenMP threads = 1
echo "-------------CORES ----------------"
echo \${ncores}
echo "-----------------------------------"

export I_MPI_PIN_CELL=core
### export SLURM_CPU_BIND=none #[IreneRome]
export I_MPI_PIN_DOMAIN=auto
export WORKDIR=\$(readlink -f \$SLURM_SUBMIT_DIR)

MDRUN="ccc_mprun gmx_mpi mdrun" #[IreneRome]
GROMPP="ccc_mprun gmx_mpi grompp"


######################################################################################

### Read the walltime and convert it in maxh value
# We need to keep some time after the MD run for the checks,
# so we substract CHECK_DURATION from the walltime to get maxh.
echo "jobid is \$SLURM_JOB_ID "
echo \$(scontrol show jobid \$SLURM_JOB_ID)
walltime=\$(scontrol show jobid \$SLURM_JOB_ID | \\
           tr ' ' '\\n' | \\
           grep "TimeLimit" | \\
           cut -f 2 -d = | \\
           awk -F '-' -v check_duration=\$CHECK_DURATION \\
               'BEGIN{seconds=0} \\
                { \\
                    if (NF==1) {days=0; hms=\$1} \\
                    else {days=\$1; hms=\$2}; \\
                    seconds=seconds + (days * 24 * 3600); \\
                    time_split_length=split(hms, time_split, ":"); \\
                    for (i=time_split_length; i>0; i--) { \\
                        seconds=seconds + \\
                                (time_split[i] \\
                                     * 60**(time_split_length - i)) \\
                    } \\
                } \\
                END{print (seconds - check_duration)/3600}')
echo "walltime is \$walltime"

### Move in the working directory
# Make sure any symbolic links are resolved to absolute path
export WORKDIR=\$(readlink -f \$SLURM_SUBMIT_DIR)
echo "Work directory is \$WORKDIR"
  
# Change to the direcotry that the job was submitted from
cd \$WORKDIR

### Identify the run
# The information on the run number is stored in the last_cycle file.
# This file needs to be created if it does not already exist.
if [[ ! -e last_cycle ]]
then
    echo 0 > last_cycle
fi

# At this point, the last_cycle file exists. We read the index of the last run
# and we increment the number
prev_cycle=\$(tail -n1 last_cycle)
cycle=\$(( \$prev_cycle + 1 ))

echo "The current cycle is \$cycle"

# Update the last_cycle file for the next round
echo \$cycle >> last_cycle

### Identify the previous run
# By default, the previous run is the one we read in the last_cycle file. If a
# run crashed before it writes a checkpoint, then we cannot start from it. So
# we need to find the last checkpoint available. If there is no checkpoint,
# then prev_cycle is 0 and we need to start from the beginning.
while [[ (! -e "\$ROOTNAME.\${prev_cycle}.cpt") && (\${prev_cycle} -gt 0) ]]
do
    let prev_cycle--
done
checkpoint=\$ROOTNAME.\$prev_cycle.cpt

# Sometime the simulation started with an other script and the checkpoint is
# not numbered. We still want to continue from this checkpoint.
if [[ (\$prev_cycle -eq 0) && (-e \$ROOTNAME.cpt) ]]
then 
    checkpoint=\$ROOTNAME.cpt
fi

echo "We will use this checkpoint: \$checkpoint"

### Is the simulation finished already?
# If the simulation already reached the number of steps requested in the TPR,
# then it is useless to start a new run.
# We first read the TRP file to know how many steps were requested.
tpr_nsteps=\$(mpirun -np 1 \${GMX} dump -s \${ROOTNAME}.tpr 2> /dev/null | grep nsteps | cut -f 2 -d = | sed 's/[^0-9]//')
# tpr_nsteps=\$(\${GMX} dump -s \${ROOTNAME}.tpr 2> /dev/null | \\            grep nsteps | cut -f 2 -d = | sed 's/[^0-9]//')

# We need to find what is the last step that has been simulated. We read it
# from the log file of the previous cycle.
last_step=\$(grep "Writing checkpoint" \$ROOTNAME.\$prev_cycle.part*.log | tail -n1 | cut -f 4 -d ' ')
echo "Requested number of steps: \$tpr_nsteps"
echo "Last step of the previous run: \$last_step"

### Run the simulation if appropriate
# We do the run if we are not done or if there is no previous log file. This
# last case can happen if (1) it is the first run or (2) we start from a run
# done with an other script. If we start from a run that used an other script,
# then we try anyway because it is painful to detect and it will not cost much
# anyway.
if [[ (\$(ls \$ROOTNAME.\$prev_cycle.part*.log | wc -l) -eq 0) || \\
      (\$last_step -lt \$tpr_nsteps) ]]
    then

    ######################################################################################
    ################# this is the job that is paralelizble and relauchable ###############
    ######################################################################################





######################### EM  - Lambda i #########################
cd Lambda_${i}/1_EM || exit
	
\$GROMPP -f em_FEP_t${t}_lambda${i}.mdp -c "../../../../${GRO}" -p "../../../../${TOP}" -o "${NAME}_em.tpr"
	
	




cd ../../Lambda_${i}/1_EM || exit
	
\$MDRUN -deffnm "${NAME}_em" -ntomp 16 -ntmpi 1



######################### NVT  - Lambda i  #########################
cd ../../Lambda_${i}/2_NVT || exit

\$GROMPP -f nvt_FEP_t${t}_lambda${i}.mdp -c "../1_EM/${NAME}_em.gro" -r "../1_EM/${NAME}_em.gro" -p "../../../../${TOP}" -o "${NAME}_nvt.tpr" -maxwarn 1






cd ../../Lambda_${i}/2_NVT || exit
	
\$MDRUN -v -deffnm "${NAME}_nvt" -ntomp 16 -ntmpi 1





######################### NPT  - Lambda i  #########################
cd ../../Lambda_${i}/3_NPT || exit

\$GROMPP -f npt_FEP_t${t}_lambda${i}.mdp -c "../2_NVT/${NAME}_nvt.gro" -r "../2_NVT/${NAME}_nvt.gro" -p "../../../../${TOP}" -o "${NAME}_npt.tpr" -maxwarn 1






cd ../../Lambda_${i}/3_NPT || exit

\$MDRUN -v -deffnm "${NAME}_npt" -ntomp 16 -ntmpi 1






#################### PRODUCTION  - Lambda i  #######################
cd ../../Lambda_${i}/4_PROD || exit
	
\$GROMPP -f prod_FEP_t${t}_lambda${i}.mdp -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${ROOTNAME}.tpr" -maxwarn 1




cd ../../Lambda_${i}/4_PROD || exit

#gmx mdrun -v -deffnm "${NAME}_prod_${t}_${i}" -ntomp 16 -ntmpi 1

\$MDRUN -nice 0 \\
        -s \$ROOTNAME \\
        -deffnm \$ROOTNAME.\$cycle \\
        -v \\
        -stepout 1000 \\
        -maxh \$walltime \\
        -cpi \$checkpoint \\
        -noappend \\
        \$MDRUN_OPT >& \$ROOTNAME.\$cycle.runout 









    ###########################################################################################################################
    ############# Check if the trajectory and the energy file are not corrupted.                                     ##########
    ############# Relounch the job if they are not done yet, and there is no corruption                              ##########
    ############# to simplify the code, EM, NVT, NPT will be unecessarily redone. PROD will catchup where it stopped ##########
    ###########################################################################################################################


    mpirun -np 1 \${GMX} check -f \$ROOTNAME.\$cycle.part*.xtc
    #\${GMX} check -f \$ROOTNAME.\$cycle.part*.xtc
    check_xtc=\$?
    echo "XTC check output code is \$check_xtc"
    mpirun -np 1 \${GMX} check -e \$ROOTNAME.\$cycle.part*.edr
    #\${GMX} check -e \$ROOTNAME.\$cycle.part*.edr
    check_edr=\$?
    echo "EDR check output code is \$check_edr"
    integrity=\$(( \$check_xtc + \$check_edr ))
    if [[ \$integrity -ne 0 ]]
    then
        echo "Error with XTC or EDR. Check for corruption."
    else
        echo "Integrity OK"
    fi
 
    # If we are not done, then we need to requeue a job
    last_step=\$(grep "Writing checkpoint" \$ROOTNAME.\$cycle.part*.log | \\
            tail -n1 | cut -f 4 -d ' ')
    echo "Last simulated step is \$last_step"
    if [[ \$last_step -lt \$tpr_nsteps && \$integrity -eq 0 ]]
    then
        ccc_msub \$BATCH_FNAME
    fi
fi

echo "Run \$cycle is done"









EOT
	

ccc_msub job.moab
echo "job sent: t ${t} Lambda ${i}"

done # lambda loop
cd ..
done # temperature loop




