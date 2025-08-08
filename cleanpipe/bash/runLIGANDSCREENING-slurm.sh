#!/bin/bash
#SBATCH --job-name=full_md_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-N   #################### xxx  Replace N with number of systems


# this script has no arguments. you have to modify the two places with "xxx"
# its just the number of systems, the folder where the system are, and the output folder names you want
# the goal here is to choose a folder that contains several protein+ligant systems, 
# where the only difference is the ligant. we are doing a ligant screening after all. the
# goal of this type of screening is to find a good ligant among several candidates.
# this script is usefull because it will run all those system in paralel.
# these is a cool logic here, where after all are run. a analysis scritp will be automatically louched
# the goal of the anaysis is to show a table telling the binding energy of all ligand canditates

# =====================================================================
# Configuration - xxx MODIFY THESE VALUES, THEY ARE FOLDER NAMES (INPUT FOLDER AND OUTPUT FOLDERS)
# =====================================================================
INPUT_DIR="input_systems"      # Parent directory with system subfolders
# this is how the INPUT_DIR should look like:
#input_systems/
#├── system1/
#│   ├── system.gro        # Full system coordinates
#│   └── topol.top         # Topology (must contain [ protein ] and [ ligand ] definitions)
#├── system2/
#│   ├── system.gro
#│   └── topol.top
#└── ...


OUTPUT_BASE="outputs"          # Base directory for simulation outputs
MMGBSA_OUT="mmgbsa_analysis"   # Final MMGBSA results directory
# this is what the outputs folders will look like:
#outputs/
#├── system1/
#│   ├── em.*, nvt.*, ... 
#│   └── system1_results.tar.gz
#├── system2/
#│   ├── ... 
#│   └── system2_results.tar.gz
#└── ...
#mmgbsa_analysis/
#├── system1_mmgbsa.dat
#├── system2_mmgbsa.dat
#└── ...






# =====================================================================
# MDP Files content
# =====================================================================
EM_MDP=$(cat <<EOF
; Energy minimization
integrator               = steep
nsteps                   = 5000
emtol                    = 100.0
emstep                   = 0.01
nstxout                  = 100
cutoff-scheme            = Verlet
vdwtype                  = Cut-off
vdw-modifier             = Force-switch
rvdw                     = 1.2
rvdw-switch              = 1.0
coulombtype              = PME
rcoulomb                 = 1.2
constraints              = none
EOF
)

NVT_MDP=$(cat <<EOF
; NVT Equilibration
integrator               = md
nsteps                   = 50000
dt                       = 0.002
nstxout                  = 1000
nstvout                  = 1000
nstenergy                = 1000
nstlog                   = 1000
cutoff-scheme            = Verlet
vdwtype                  = Cut-off
vdw-modifier             = Force-switch
rvdw                     = 1.2
rvdw-switch              = 1.0
coulombtype              = PME
rcoulomb                 = 1.2
constraints              = h-bonds
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 298
EOF
)

NPT_MDP=$(cat <<EOF
; NPT Equilibration
integrator               = md
nsteps                   = 100000
dt                       = 0.002
nstxout                  = 1000
nstvout                  = 1000
nstenergy                = 1000
nstlog                   = 1000
cutoff-scheme            = Verlet
vdwtype                  = Cut-off
vdw-modifier             = Force-switch
rvdw                     = 1.2
rvdw-switch              = 1.0
coulombtype              = PME
rcoulomb                 = 1.2
constraints              = h-bonds
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 298
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 1.0
ref-p                    = 1.0
compressibility          = 4.5e-5
EOF
)

PROD_MDP=$(cat <<EOF
; Production MD
integrator               = md
nsteps                   = 5000000 ; 10 ns
dt                       = 0.002
nstxout                  = 0
nstvout                  = 0
nstenergy                = 5000
nstlog                   = 5000
nstxout-compressed       = 5000
cutoff-scheme            = Verlet
vdwtype                  = Cut-off
vdw-modifier             = Force-switch
rvdw                     = 1.2
rvdw-switch              = 1.0
coulombtype              = PME
rcoulomb                 = 1.2
constraints              = h-bonds
continuation             = yes
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 298
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 1.0
ref-p                    = 1.0
compressibility          = 4.5e-5
EOF
)

# =====================================================================
# Folders and files construction
# =====================================================================
module load gromacs/2023.3

# Get system directory for this array task
SYSTEM_LIST=($(ls -d ${INPUT_DIR}/*/))
SYSTEM_DIR=${SYSTEM_LIST[$SLURM_ARRAY_TASK_ID - 1]}
SYSTEM_NAME=$(basename ${SYSTEM_DIR})

# Create output directories
WORK_DIR="${OUTPUT_BASE}/${SYSTEM_NAME}"
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p ${WORK_DIR} ${LOG_DIR}

# Generate MDP files
echo "${EM_MDP}" > ${WORK_DIR}/em.mdp
echo "${NVT_MDP}" > ${WORK_DIR}/nvt.mdp
echo "${NPT_MDP}" > ${WORK_DIR}/npt.mdp
echo "${PROD_MDP}" > ${WORK_DIR}/prod.mdp

# Convert to absolute paths for reliability
SYSTEM_DIR=$(realpath ${SYSTEM_DIR})
WORK_DIR=$(realpath ${WORK_DIR})

# =====================================================================
# Simulations
# =====================================================================
cd ${WORK_DIR}

# 1. Energy Minimization
echo "Starting EM for ${SYSTEM_NAME}"
gmx grompp -f em.mdp \
           -c ${SYSTEM_DIR}/system.gro \
           -p ${SYSTEM_DIR}/topol.top \
           -o em.tpr \
           -po em_out.mdp \
           -maxwarn 1 > em_grompp.log 2>&1
gmx mdrun -v -deffnm em \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > em_mdrun.log 2>&1

# 2. NVT Equilibration
echo "Starting NVT for ${SYSTEM_NAME}"
gmx grompp -f nvt.mdp \
           -c em.gro \
           -r em.gro \
           -p ${SYSTEM_DIR}/topol.top \
           -o nvt.tpr \
           -po nvt_out.mdp \
           -maxwarn 1 > nvt_grompp.log 2>&1
gmx mdrun -v -deffnm nvt \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > nvt_mdrun.log 2>&1

# 3. NPT Equilibration
echo "Starting NPT for ${SYSTEM_NAME}"
gmx grompp -f npt.mdp \
           -c nvt.gro \
           -t nvt.cpt \
           -p ${SYSTEM_DIR}/topol.top \
           -o npt.tpr \
           -po npt_out.mdp \
           -maxwarn 1 > npt_grompp.log 2>&1
gmx mdrun -v -deffnm npt \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > npt_mdrun.log 2>&1

# 4. Production MD
echo "Starting Production for ${SYSTEM_NAME}"
gmx grompp -f prod.mdp \
           -c npt.gro \
           -t npt.cpt \
           -p ${SYSTEM_DIR}/topol.top \
           -o prod.tpr \
           -po prod_out.mdp \
           -maxwarn 1 > prod_grompp.log 2>&1
gmx mdrun -v -deffnm prod \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > prod_mdrun.log 2>&1

# =====================================================================
# Generate Analysis Script 
# (the script will be generated and run just once, because just the first job will match the 'if' below)
# (the analysis script will wait for all jobs to end, before the --dependency flag will be used)
# =====================================================================
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    cat > analysis_after_paralel_jobs-slurm.sh << 'ANALYSIS_EOF'
#!/bin/bash
#SBATCH --job-name=mmgbsa_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/mmgbsa_%j.out

# Inherit environment variables from parent
source $SLURM_SUBMIT_DIR/$0

# =====================================================================
# MMGBSA Analysis (will run just if all the paralel simulations are finished)
# =====================================================================
module load gromacs/2023.3

MMGBSA_MDP=$(cat << 'MMGBSA_EOF'
; MMGBSA Parameters
gb-method                = OBC
alpb                     = yes
sa-method                = Ace-approximation
sa-surface-tension       = 0.0226778 ; kJ/mol/nm²
temperature              = 298        ; K
nstgbradii               = 1
soft-core                = no
MMGBSA_EOF
)

mkdir -p ${MMGBSA_OUT}

# Process all systems
echo "Starting MMGBSA analysis"
for system_dir in ${OUTPUT_BASE}/*/; do
    SYSTEM_NAME=$(basename ${system_dir})
    cd "${system_dir}"
    
    echo "Processing ${SYSTEM_NAME}"
    
    # Create index groups
    echo -e "Protein\nLigand\nq" | gmx make_ndx -f prod.tpr -o mmgbsa_groups.ndx > make_ndx.log 2>&1
    
    # Get group numbers
    PROTEIN_GRP=$(grep -A1 "Protein" mmgbsa_groups.ndx | head -1 | cut -d' ' -f1 | tr -d [ | tr -d ])
    LIGAND_GRP=$(grep -A1 "Ligand" mmgbsa_groups.ndx | head -1 | cut -d' ' -f1 | tr -d [ | tr -d ])
    
    # Run MMGBSA
    echo "${MMGBSA_MDP}" > mmgbsa.mdp
    gmx mmgbsa -s prod.tpr \
               -i mmgbsa.mdp \
               -ci mmgbsa_groups.ndx \
               -cg ${PROTEIN_GRP} ${LIGAND_GRP} \
               -ct prod.trr \
               -o ${MMGBSA_OUT}/${SYSTEM_NAME}_mmgbsa.dat > mmgbsa.log 2>&1
               
    # Cleanup
    mkdir -p ${SYSTEM_NAME}_results
    mv em.* nvt.* npt.* prod.* mmgbsa* ${SYSTEM_NAME}_results/
    tar -czf ${SYSTEM_NAME}_results.tar.gz ${SYSTEM_NAME}_results
    rm -rf ${SYSTEM_NAME}_results
done

echo "Completed MMGBSA analysis for all systems"
ANALYSIS_EOF

    # Make executable and submit with dependency
    chmod +x analysis_after_paralel_jobs-slurm.sh
    ANALYSIS_JOBID=$(sbatch --parsable --dependency=afterok:$SLURM_ARRAY_JOB_ID analysis_after_paralel_jobs-slurm.sh)
    echo "Submitted analysis job ID: $ANALYSIS_JOBID with dependency on array: $SLURM_ARRAY_JOB_ID"
fi

# =====================================================================
# Final Message, indicating completion of the simulation for one of the paralel jobs
# =====================================================================
echo "Completed simulation for ${SYSTEM_NAME}"