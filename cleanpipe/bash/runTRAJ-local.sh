#!/bin/bash

# usage example:  ./runTRAJ-local.sh 1234 dna_in_water.gro dna_in_water.top 100 "System" 298
# you have to be inside the system folder, where the top and gro are
# the folders will be created where the script is

# these are the arguments:
#procedures that must be done. 4 or 1234. (4 mean only PRODTRAJ, 1234 mean EM NVT NPT PRODTRAJ)
#name of gro
#name of top
#production time (nanoseconds)
#groups to monitor. ex: "Protein Non-Protein" or "System" or "Protein Water" or whatever you want
#temperature. ex: 298

#attention: for now this script assumes the existence of a Protein. this is used in temperature and pressure coupling and in the center and fit commands


if [ "$#" -ne 6 ]; then
    echo "Wrong number or arguments. Usage: $0 <procedures> <gro_file> <top_file> <production_duration_ns> <groups_to_monitor> <temperature>"
    exit 1
fi

PROCEDURE=$1
echo "Procedures to be done (4 mean only PRODTRAJ, 1234 mean EM NVT NPT PRODTRAJ): ${PROCEDURES}"

GRO=$2
echo "Input gro file: ${GRO}"

TOP=$6
echo "Input top file: ${TOP}"

PRODUCTION_DURATION=$4
STEPS=$(( (PRODUCTION_DURATION * 1000000) / 2 ))
echo "Duration of production MD : ${PRODUCTION_DURATION} ns"
echo "Steps for that duration : ${STEPS} (time x 1000 / 0.002)"


GROUPS_TO_MONITOR=$5
echo "Groups to monitor : ${GROUPS_TO_MONITOR}"

TEMPERATURE=$6
CLEAN_STRING=$(echo "$GROUPS_TO_MONITOR" | tr '\t' ' ' | xargs)
WORD_COUNT=$(echo "$CLEAN_STRING" | wc -w)
TEMPERATURES=$(yes $TEMPERATURE | head -n "$WORD_COUNT" | paste -sd ' ' -)
ONES=$(yes 1 | head -n "$WORD_COUNT" | paste -sd ' ' -)
echo "Temperature for each group: ${TEMPERATURES}"
echo "Tau_t for each group: ${ONES}"

echo "     "

if [[ ! -f "${GRO}" ]]; then
    echo "Error: GRO file '${GRO}' does not exist."
    exit 1
fi



if [[ ! -f "${TOP}" ]]; then
    echo "Error: TOP file '${TOP}' does not exist."
    exit 1
fi


###################### check top and itps for [ distance_restraints ] or [ dihedral_restraints ] ######################
echo "looking for restraints in the topology files..."




# === Config ===
DIHRE_OPTION="no"
DISRE_OPTION="no"

# === Ensure TOP is defined ===
if [[ -z "$TOP" ]]; then
    echo "Error: TOP variable is not set."
    exit 1
fi

# === Gather files ===
ITP_FILES=( ./*.itp )
ALL_FILES=( "${ITP_FILES[@]}" "$TOP" )

# === Check restraints in all files ===
for file in "${ALL_FILES[@]}"; do
    [[ -f "$file" ]] || continue  # Skip if not a real file

    if grep -q '\[ *dihedral_restraints *\]' "$file"; then
        DIHRE_OPTION="yes"
    fi
    if grep -q '\[ *distance_restraints *\]' "$file"; then
        DISRE_OPTION="simple"
    fi
done

# === Output ===
echo "Included .itp files:"
for itp in "${ITP_FILES[@]}"; do
    echo "  $itp"
done

echo ""
echo "DIHRE_OPTION=$DIHRE_OPTION"
echo "DISRE_OPTION=$DISRE_OPTION"





echo ""


######################### create folder structure ########################
echo "creating folders..."

if [[ "${PROCEDURES}" == *"1"* ]]; then
    if [[ ! -d "1_EM" ]]; then
        mkdir 1_EM
        echo "  folder 1_EM created"
    else
        echo "  folder 1_EM already exists, skipping creation."
    fi
fi


if [[ "${PROCEDURES}" == *"2"* ]]; then
    if [[ ! -d "2_NVT" ]]; then
        mkdir 2_NVT
        echo "  folder 2_NVT created"
    else
        echo "  folder 2_NVT already exists, skipping creation."
    fi
fi


if [[ "${PROCEDURES}" == *"3"* ]]; then
    if [[ ! -d "3_NPT" ]]; then
        mkdir 3_NPT
        echo "  folder 3_NPT created"
    else
        echo "  folder 3_NPT already exists, skipping creation."
    fi
fi


if [[ "${PROCEDURES}" == *"4"* ]]; then
    if [[ ! -d "4_PRODTRAJ" ]]; then
        mkdir 4_PRODTRAJ
        echo "  folder 4_PRODTRAJ created"
    else
        echo "  folder 4_PRODTRAJ already exists, skipping creation."
    fi
fi






#################### generate mdp files #####################
echo "creating mdp files..."

# The name of the file you want to write to
file_em_mdp="em.mdp"
file_nvt_mdp="nvt.mdp"
file_npt_mdp="npt.mdp"
file_prod_mdp="prod.mdp"




######### mdp for em #########
if [[ "${PROCEDURES}" == *"1"* ]]; then
cat <<EOT > "1_EM/${file_em_mdp}"



Integrator =	steep
emtol      =	1000 ;100
emstep     =	0.01
nsteps     =    100000 ; this is the max value to be used just if emtol is never reached

;box configuration	
pbc            = xyz


; fix distances
disre = ${DISRE_OPTION}
disre_fc = 1000
; fix dihedrals
dihre = ${DIHRE_OPTION}
dihre_fc = 1000



;LONG RANGE ESTIMATION	
rcoulomb       =	1.2
coulombtype    =	PME
;pme options	
pme_order      =	4
fourierspacing =	0.12
ewald_rtol     =	1.00E-05
	
rvdw           =	1.2
vdw_type       =	cutoff
;cutoff workarounds	
vdw-modifier   =	force-switch
rvdw-switch    =	1
DispCorr       =	no



;this is a bug workaround. this should be usefull just for FEP. this is not the case, but if I dont set it, the program gives an error
sc-r-power = 6


EOT
echo "  1_EM/${file_em_mdp} created"
fi




######### mdp for nvt #########
if [[ "${PROCEDURES}" == *"2"* ]]; then
cat <<EOT > "2_NVT/${file_nvt_mdp}"


Integrator =	md	
dt         =	0.002
nsteps     =	50 ;50000 ; (100 ps)

;box configuration	
pbc                     = 	xyz


; activate pinning of proteins or water flexibility	
define               =	-DPOSRES
refcoord_scaling     = 
; stiffen bonds	
constraints          =	h-bonds
constraint_algorithm =	lincs
lincs_iter           =	1
lincs_order          =	4
; fix distances
disre = ${DISRE_OPTION}
disre_fc = 1000
; fix dihedrals
dihre = ${DIHRE_OPTION}
dihre_fc = 1000

;LONG RANGE ESTIMATION	

rcoulomb       =	1.2
coulombtype    =	PME
;pme options	
pme_order      =	4
fourierspacing =	0.12
ewald_rtol     =	1.00E-05
	
rvdw           =	1.2
vdw_type       =	cutoff
;cutoff workarounds	
vdw-modifier   =	force-switch
rvdw-switch    =	1
DispCorr       =	EnerPres
	
;NEIGHBOUR TRACKING	
cutoff-scheme  =	Verlet
ns_type        = 	grid
rlist          =	1.2



; velocity assingment	
continuation =	no
gen_vel      =	yes
gen_temp     =	310


; Temperature coupling	
tcoupl    =	V-rescale
tc-grps   =	${GROUPS_TO_MONITOR}
ref_t     = ${TEMPERATURES}
tau_t     =	${ONES}

; Pressure coupling	
pcoupl                  = 	no



; output control	
; TRR	
nstxout    =	50000
nstvout    =	50000
nstfout    =	
; EDR	
nstenergy  =	50000
energygrps =	
; LOG	
nstlog     =	50000
; XTC instead of TRR	
nstxout-compressed =	
compressed-x-grps =	



EOT
echo "  2_NVT/${file_nvt_mdp} created"
fi




######### mdp for npt #########
if [[ "${PROCEDURES}" == *"3"* ]]; then
cat <<EOT > "3_NPT/${file_npt_mdp}"

Integrator =	md	
dt         =	0.002
nsteps     =	50 ;50000 ; (100 ps)


;box configuration	
pbc                     = 	xyz


; activate pinning of proteins or water flexibility	
define               =  -DPOSRES
refcoord_scaling     =  com
; stiffen bonds	
constraints          =	h-bonds
constraint_algorithm =	lincs
lincs_iter           =	1
lincs_order          =	4
; fix distances
disre = ${DISRE_OPTION}
disre_fc = 1000
; fix dihedrals
dihre = ${DIHRE_OPTION}
dihre_fc = 1000

;LONG RANGE ESTIMATION	

rcoulomb       =	1.2
coulombtype    =	PME
;pme options	
pme_order      =	4
fourierspacing =	0.12
ewald_rtol     =	1.00E-05
	
rvdw           =	1.2
vdw_type       =	cutoff
;cutoff workarounds	
vdw-modifier   =	force-switch
rvdw-switch    =	1
DispCorr       =	EnerPres
	
;NEIGHBOUR TRACKING	
cutoff-scheme  =	Verlet
ns_type        = 	grid
rlist          =	1.2



; velocity assingment	
continuation =	yes
gen_vel      =	no
gen_temp     =	


; Temperature coupling	
tcoupl    =	V-rescale
tc-grps   =	${GROUPS_TO_MONITOR}
ref_t     = ${TEMPERATURES}
tau_t     =	${ONES}

; Pressure coupling	
pcoupl          =	C-rescale
pcoupltype      =	isotropic
ref_p           =	1
tau_p           =	5
compressibility =	4.50E-05




; output control	
; TRR	
nstxout            =	50000
nstvout            =	50000
nstfout            =	
; EDR	
nstenergy          =	50000
energygrps         =	
; LOG	
nstlog             =	50000
; XTC instead of TRR	
nstxout-compressed =	
compressed-x-grps  =	


EOT
echo "  3_NPT/${file_npt_mdp} created"
fi




######### mdp for prod #########
if [[ "${PROCEDURES}" == *"4"* ]]; then
cat <<EOT > "4_PRODTRAJ/${file_prod_mdp}"


Integrator =	md	
dt         =	0.002 ; 2 femtoseconds. without lincs this would have to be 0.0001
nsteps     =	${STEPS}



;box configuration	
pbc                     = 	xyz


; activate pinning of proteins or water flexibility	
define               =  
refcoord_scaling     =  
; stiffen bonds	
constraints          =	h-bonds
constraint_algorithm =	lincs
lincs_iter           =	1
lincs_order          =	4
; fix distances
disre = ${DISRE_OPTION}
disre_fc = 1000
; fix dihedrals
dihre = ${DIHRE_OPTION}
dihre_fc = 1000

;LONG RANGE ESTIMATION	

rcoulomb       =	1.2
coulombtype    =	PME
;pme options	
pme_order      =	4
fourierspacing =	0.12
ewald_rtol     =	1.00E-05
	
rvdw           =	1.2
vdw_type       =	cutoff
;cutoff workarounds	
vdw-modifier   =	force-switch
rvdw-switch    =	1
DispCorr       =	EnerPres
	
;NEIGHBOUR TRACKING	
cutoff-scheme  =	Verlet
ns_type        = 	grid
rlist          =	1.2



; velocity assingment	
continuation =	yes
gen_vel      =	no
gen_temp     =	


; Temperature coupling	
tcoupl    =	V-rescale
tc-grps   =	${GROUPS_TO_MONITOR}
ref_t     = ${TEMPERATURES}
tau_t     =	${ONES}

; Pressure coupling	
pcoupl          =	C-rescale
pcoupltype      =	isotropic
ref_p           =	1
tau_p           =	5
compressibility =	4.50E-05




; output control	
; TRR	
nstxout            =	50000
nstvout            =	50000
nstfout            =	50000
; EDR	
nstenergy          =	50000
energygrps         =	
; LOG	
nstlog             =	50000
; XTC instead of TRR	
nstxout-compressed =	5000
compressed-x-grps  =	System



EOT
echo "  4_PRODTRAJ/${file_prod_mdp} created"
fi








if [[ "${PROCEDURES}" == "1"* ]]; then
    echo "first procedure is EM. using inputed GRO"
fi


######################### EM prep #########################"
if [[ "${PROCEDURES}" == *"1"* ]]; then
    echo "executing grompp for EM..."
    # Check if the mdp file exist
    if [[ ! -f "1_EM/${file_em_mdp}"  ]]; then
        echo "  Error: 1_EM/${file_em_mdp} do not exist."
        exit 1
    else
        
        cd 1_EM || exit
        gmx grompp -f $file_em_mdp -c "../${GRO}" -p "../${TOP}" -o "em.tpr" > grompp.out 2> grompp.err # HERE!!!!!!!!!
        # Check if there were any errors during grompp
        if [ $? -ne 0 ]; then
            echo "  grompp failed"
            cat grompp.err
            exit 1
        else
            echo "  grompp succeeded"
        fi

        cd ..

    fi
fi


######################### EM run #########################"
if [[ "${PROCEDURES}" == *"1"* ]]; then
    echo "executing mdrun for EM..."
    # Check if the tpr file exist
    if [[ ! -f "1_EM/em.tpr"  ]]; then
        echo "  Error: 1_EM/em.tpr do not exist."
        exit 1
    else
        
        cd 1_EM || exit
        gmx mdrun -deffnm "em" > mdrun.out 2> mdrun.err # HERE!!!!!!!!!
        # Check if there were any errors during mdrun
        if [ $? -ne 0 ]; then
            echo "  mdrun failed"
            cat mdrun.err
            exit 1
        else
            echo "  mdrun succeeded"
        fi
        cd ..
    fi
fi




if [[ "${PROCEDURES}" == "2"* ]]; then
    echo "first procedure is NVT. using inputed GRO"
else
    GRO="1_EM/em.gro"
fi


######################### NVT prep #########################
if [[ "${PROCEDURES}" == *"2"* ]]; then
    echo "executing grompp for NVT..."
    # Check if the mdp file exist
    if [[ ! -f "2_NVT/${file_nvt_mdp}"  ]]; then
        echo "  Error: 2_NVT/${file_nvt_mdp} do not exist."
        exit 1
    else
        
        cd 2_NVT || exit
        gmx grompp -f $file_nvt_mdp -c "../${GRO}" -r "../1_EM/em.gro" -p "../${TOP}" -o "nvt.tpr" -maxwarn 1 > grompp.out 2> grompp.err # HERE!!!!!!!!!
        # Check if there were any errors during grompp
        if [ $? -ne 0 ]; then
            echo "  grompp failed"
            cat grompp.err
            exit 1
        else
            echo "  grompp succeeded"
        fi
        cd ..

    fi
fi


######################### NVT run #########################
if [[ "${PROCEDURES}" == *"2"* ]]; then
    echo "executing mdrun for NVT..."
    # Check if the tpr file exist
    if [[ ! -f "2_NVT/nvt.tpr"  ]]; then
        echo "  Error: 2_NVT/nvt.tpr do not exist."
        exit 1
    else
        
        cd 2_NVT || exit
        gmx mdrun -v -deffnm "nvt" > mdrun.out 2> mdrun.err # HERE!!!!!!!!!
        # Check if there were any errors during mdrun
        if [ $? -ne 0 ]; then
            echo "  mdrun failed"
            cat mdrun.err
            exit 1
        else
            echo "  mdrun succeeded"
        fi
        cd ..

    fi
fi






if [[ "${PROCEDURES}" == "3"* ]]; then
    echo "first procedure is NPT. using inputed GRO"
else
    GRO="2_NVT/nvt.gro"
fi

######################### NPT prep #########################
if [[ "${PROCEDURES}" == *"3"* ]]; then
    echo "executing grompp for NPT..."
    # Check if the mdp file exist
    if [[ ! -f "3_NPT/${file_npt_mdp}"  ]]; then
        echo "  Error: 3_NPT/${file_npt_mdp} do not exist."
        exit 1
    else
        
        cd 3_NPT || exit
        gmx grompp -f $file_npt_mdp -c "../${GRO}" -r "../2_NVT/nvt.gro" -p "../${TOP}" -o "npt.tpr" -maxwarn 1 > grompp.out 2> grompp.err # HERE!!!!!!!!!
        # Check if there were any errors during grompp
        if [ $? -ne 0 ]; then
            echo "  grompp failed"
            cat grompp.err
            exit 1
        else
            echo "  grompp succeeded"
        fi
        cd ..

    fi
fi


######################### NPT run #########################
if [[ "${PROCEDURES}" == *"3"* ]]; then
    echo "executing mdrun for NPT..."
    # Check if the tpr file exist
    if [[ ! -f "3_NPT/npt.tpr"  ]]; then
        echo "  Error: 3_NPT/npt.tpr do not exist."
        exit 1
    else
        
        cd 3_NPT || exit
        gmx mdrun -v -deffnm "npt" > mdrun.out 2> mdrun.err # HERE!!!!!!!!!
        # Check if there were any errors during mdrun
        if [ $? -ne 0 ]; then
            echo "  mdrun failed"
            cat mdrun.err
            exit 1
        else
            echo "  mdrun succeeded"
        fi
        cd ..

    fi
fi




if [[ "${PROCEDURES}" == "4"* ]]; then
    echo "first procedure is PRODTRAJ. using inputed GRO"
else
    GRO="3_NPT/npt.gro"
fi

######################### PROD prep #########################
if [[ "${PROCEDURES}" == *"4"* ]]; then
    echo "executing grompp for PRODTRAJ..."
    # Check if the mdp file exist
    if [[ ! -f "4_PRODTRAJ/${file_prod_mdp}"  ]]; then
        echo "  Error: 4_PRODTRAJ/${file_prod_mdp} do not exist."
        exit 1
    else
        
        cd 4_PRODTRAJ || exit
        gmx grompp -f $file_prod_mdp -c "../${GRO}" -p "../${TOP}" -o "prod.tpr" -maxwarn 1 > grompp.out 2> grompp.err # HERE!!!!!!!!!
        # Check if there were any errors during grompp
        if [ $? -ne 0 ]; then
            echo "  grompp failed"
            cat grompp.err
            exit 1
        else
            echo "  grompp succeeded"
        fi
        cd ..

    fi
fi


######################### PROD run #########################
if [[ "${PROCEDURES}" == *"4"* ]]; then
    echo "executing mdrun for PRODTRAJ..."
    # Check if the tpr file exist
    if [[ ! -f "4_PRODTRAJ/prod.tpr"  ]]; then
        echo "  Error: 4_PRODTRAJ/prod.tpr do not exist."
        exit 1
    else
        
        cd 4_PRODTRAJ || exit
        gmx mdrun -v -deffnm "prod" > mdrun.out 2> mdrun.err # HERE!!!!!!!!!
        # Check if there were any errors during mdrun
        if [ $? -ne 0 ]; then
            echo "  mdrun failed"
            cat mdrun.err
            exit 1
        else
            echo "  mdrun succeeded"
        fi
        cd ..

    fi
fi


######################### PROD center and fit #########################
if [[ "${PROCEDURES}" == *"4"* ]]; then
    echo "executing center and fit for PROD..."
    # Check if the xtc file exist
    if [[ ! -f "4_PRODTRAJ/prod.xtc"  ]]; then
        echo "  Error: 4_PRODTRAJ/prod.xtc do not exist."
        exit 1
    else
        
        cd 4_PRODTRAJ || exit         
        printf '1\n0' | gmx trjconv -s "prod.tpr" -f "prod.xtc" -o "prod_centered.xtc" -center -pbc mol > center.out 2> center.err # HERE!!!!!!!!!
        if [ $? -ne 0 ]; then
            echo "  center failed"
            cat center.err
            exit 1
        else
            echo "  center succeeded"
        fi

        printf '1\n0' | gmx trjconv -s "prod.tpr" -f "prod_centered.xtc" -o "prod_fitted.xtc" -fit progressive > fit.out 2> fit.err   # HERE!!!!!!!!!
        if [ $? -ne 0 ]; then
            echo "  fit failed"
            cat fit.err
            exit 1
        else
            echo "  fit succeeded"
        fi

        cd ..

    fi
fi



