#!/bin/bash

# usage example:  ./job.sh insulin insulin.gro insulin.top 0.1

if [ $# -lt 4 ]; then
    echo "4 arguments needed :          name       filename.gro   filename.top   numberOfNanoseconds"
    echo "for example :        ./job.sh hemoglobin hemoglobin.gro hemoglobin.top 10"
    exit 1
fi




####### this is what changes because of FEP #######
#the simulation pipeline (EM,NVT, NPT, PROD) will have to be repeated 40 times for the same system,
#so fist a huge file structure is created. what I mean is that 40 folders are created inside the _runFEB folder, each one with the 4 steps of the pipeline
#then all the mdp files are created (what changes between the lambdas is just the paramenter 'init_lambda_state' )


###### this is what changes because of oxygen #####
#in every mdrun, "-nt 4" becomes "-dd 3 3 3 -ntmpi 40 -npme 13"
#the grompp and mdrun commands will be sent using sbatch, which requires the creation of a temporary file with a weird header





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
STEPS=$(echo "scale=0; ($PRODUCTION_DURATION / 0.02) * 10000" | bc)
echo "Duration of production MD: ${PRODUCTION_DURATION} ns (that duration is implemented in the mdp file by setting nsteps=${STEPS} and dt=0.002)"
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

cat <<EOT > "temp.sh"
#!/bin/bash

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=40
#SBATCH --gres=gpu:1
##SBATCH --mem-per-cpu=1GB
#SBATCH --nodes=1
#SBATCH --job-name=${NAME}${t}${i}
#SBATCH --output=log.out_${i}
#SBATCH --exclude=node-15


module load gromacs/2023




echo "######################### EM prep - Lambda ${i} #########################"
cd Lambda_${i}/1_EM || exit
	
gmx grompp -f em_FEP_t${t}_lambda${i}.mdp -c "../../../../${GRO}" -p "../../../../${TOP}" -o "${NAME}_em.tpr"
	
	



echo "######################### EM run - Lambda ${i}  ##############################"
cd ../../Lambda_${i}/1_EM || exit
	
gmx mdrun -deffnm "${NAME}_em" -ntomp 16 -ntmpi 1



echo "######################### NVT prep - Lambda ${i}  #########################"
cd ../../Lambda_${i}/2_NVT || exit

gmx grompp -f nvt_FEP_t${t}_lambda${i}.mdp -c "../1_EM/${NAME}_em.gro" -r "../1_EM/${NAME}_em.gro" -p "../../../../${TOP}" -o "${NAME}_nvt.tpr" -maxwarn 1





echo "######################### NVT run - Lambda ${i}  ##############################"
cd ../../Lambda_${i}/2_NVT || exit
	
gmx mdrun -v -deffnm "${NAME}_nvt" -ntomp 16 -ntmpi 1





echo "######################### NPT prep - Lambda ${i}  #########################"
cd ../../Lambda_${i}/3_NPT || exit

gmx grompp -f npt_FEP_t${t}_lambda${i}.mdp -c "../2_NVT/${NAME}_nvt.gro" -r "../2_NVT/${NAME}_nvt.gro" -p "../../../../${TOP}" -o "${NAME}_npt.tpr" -maxwarn 1





echo "######################### NPT run - Lambda ${i}  ##############################"
cd ../../Lambda_${i}/3_NPT || exit

gmx mdrun -v -deffnm "${NAME}_npt" -ntomp 16 -ntmpi 1






echo "#################### PRODUCTION prep - Lambda ${i}  #######################"
cd ../../Lambda_${i}/4_PROD || exit
	
gmx grompp -f prod_FEP_t${t}_lambda${i}.mdp -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${NAME}_prod_${t}_${i}.tpr" -maxwarn 1

#gmx grompp -f bench_FEP_t${t}_lambda${i}.mdp -c "../3_NPT/${NAME}_npt.gro" -p "../../../../${TOP}" -o "${NAME}_quick_${t}_${i}.tpr" -maxwarn 1





echo "#################### PRODUCTION run - Lambda ${i}  ###################################"
cd ../../Lambda_${i}/4_PROD || exit

#gmx mdrun -v -deffnm "${NAME}_prod_${t}_${i}" -ntomp 16 -ntmpi 1





echo "#################### PRODUCTION recenter ###################################"
cd ../../Lambda_${i}/4_PROD || exit

#printf '1\n0' | gmx trjconv -s "${NAME}_prod_${t}_${i}.tpr" -f "${NAME}_prod_${t}_${i}.xtc" -o "${NAME}_prod_${t}_${i}.centered.xtc" -center -pbc mol

cd ../..



EOT
	
chmod +x temp.sh
sbatch temp.sh
rm temp.sh
echo "job sent: t ${t} Lambda ${i}"

done # lambda loop
cd ..
done # temperature loop




