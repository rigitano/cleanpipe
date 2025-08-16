#!/bin/bash
#SBATCH --job-name=fep_md_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=0-62   # 21 lambdas × 3 temperatures = 63 tasks

# =====================================================================
# Configuration - MODIFY THESE VALUES
# =====================================================================
INPUT_DIR="input_systems"      # Directory with lambda subfolders (L_00, L_01, etc.)
FEP_OUT_BASE="FEP"             # Base directory for FEP simulations
ANALYSIS_OUT="FEP_ANALYSIS"    # FEP analysis results directory

# Force field settings - MUST match your system preparation
FORCEFIELD="charmm36-mar2019"  # charmm36, amber99sb, etc.
WATER_MODEL="tip3p"             # tip3p, spce, tip4p, etc.

# Define temperatures and lambda range
TEMPERATURES=(290 300 310)
LAMBDAS=($(seq 0 20))
TOTAL_TASKS=${#TEMPERATURES[@]}
TASKS_PER_TEMP=$(( (${#LAMBDAS[@]} + TOTAL_TASKS - 1) / TOTAL_TASKS ))

# Calculate current temperature and lambda
TEMP_INDEX=$((SLURM_ARRAY_TASK_ID / ${#LAMBDAS[@]}))
LAMBDA_INDEX=$((SLURM_ARRAY_TASK_ID % ${#LAMBDAS[@]}))
CURRENT_TEMP=${TEMPERATURES[$TEMP_INDEX]}
CURRENT_LAMBDA=${LAMBDAS[$LAMBDA_INDEX]}

# Format lambda with leading zeros
LAMBDA_DIR=$(printf "L_%02d" $CURRENT_LAMBDA)

# =====================================================================
# MDP File Templates - EDIT PARAMETERS AS NEEDED
# =====================================================================
# Helper function to set temperature and lambda in MDP
set_params() {
    local mdp="$1"
    mdp=${mdp//REF_TEMP/$CURRENT_TEMP}
    mdp=${mdp//INIT_LAMBDA/$CURRENT_LAMBDA}
    echo "$mdp"
}

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

NVT_MDP=$(set_params "$(cat <<EOF
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
ref-t                    = REF_TEMP  ; Will be replaced
; Free energy parameters
free_energy              = yes
init_lambda_state        = INIT_LAMBDA  ; Will be replaced
calc_lambda_neighbors    = 1
fep_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
sc_alpha                 = 0.5
sc_power                 = 1
sc_sigma                 = 0.3
couple-moltype           = Ligand
couple-intramol          = no
EOF
)")

NPT_MDP=$(set_params "$(cat <<EOF
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
ref-t                    = REF_TEMP  ; Will be replaced
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 1.0
ref-p                    = 1.0
compressibility          = 4.5e-5
; Free energy parameters
free_energy              = yes
init_lambda_state        = INIT_LAMBDA  ; Will be replaced
calc_lambda_neighbors    = 1
fep_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
sc_alpha                 = 0.5
sc_power                 = 1
sc_sigma                 = 0.3
couple-moltype           = Ligand
couple-intramol          = no
EOF
)")

PROD_MDP=$(set_params "$(cat <<EOF
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
ref-t                    = REF_TEMP  ; Will be replaced
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 1.0
ref-p                    = 1.0
compressibility          = 4.5e-5
; Free energy parameters
free_energy              = yes
init_lambda_state        = INIT_LAMBDA  ; Will be replaced
calc_lambda_neighbors    = 1
fep_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
sc_alpha                 = 0.5
sc_power                 = 1
sc_sigma                 = 0.3
couple-moltype           = Ligand
couple-intramol          = no
EOF
)")

# =====================================================================
# Main Pipeline
# =====================================================================
module load gromacs/2023.3

# Create directory structure
TEMP_DIR="${FEP_OUT_BASE}/${CURRENT_TEMP}K"
LAMBDA_PATH="${TEMP_DIR}/${LAMBDA_DIR}"
WORK_DIR="${LAMBDA_PATH}/4_PROD"  # Production is final output
STAGE_DIRS=("1_EM" "2_NVT" "3_NPT" "4_PROD")

# Create all required directories
mkdir -p "${ANALYSIS_OUT}"
for stage in "${STAGE_DIRS[@]}"; do
    mkdir -p "${LAMBDA_PATH}/${stage}"
done

# Generate MDP files in respective directories
echo "${EM_MDP}" > "${LAMBDA_PATH}/1_EM/em.mdp"
echo "${NVT_MDP}" > "${LAMBDA_PATH}/2_NVT/nvt.mdp"
echo "${NPT_MDP}" > "${LAMBDA_PATH}/3_NPT/npt.mdp"
echo "${PROD_MDP}" > "${WORK_DIR}/prod.mdp"

# Convert to absolute paths
INPUT_DIR=$(realpath "${INPUT_DIR}/${LAMBDA_DIR}")
WORK_DIR=$(realpath "${WORK_DIR}")
LAMBDA_PATH=$(realpath "${LAMBDA_PATH}")

# =====================================================================
# Simulation Pipeline
# =====================================================================
# 1. Energy Minimization
cd "${LAMBDA_PATH}/1_EM"
gmx grompp -f em.mdp \
           -c "${INPUT_DIR}/system.gro" \
           -p "${INPUT_DIR}/topol.top" \
           -o em.tpr \
           -po em_out.mdp \
           -maxwarn 1 > em_grompp.log 2>&1

gmx mdrun -v -deffnm em \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > em_mdrun.log 2>&1

# 2. NVT Equilibration
cd "${LAMBDA_PATH}/2_NVT"
gmx grompp -f nvt.mdp \
           -c "../1_EM/em.gro" \
           -r "../1_EM/em.gro" \
           -p "${INPUT_DIR}/topol.top" \
           -o nvt.tpr \
           -po nvt_out.mdp \
           -maxwarn 1 > nvt_grompp.log 2>&1

gmx mdrun -v -deffnm nvt \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > nvt_mdrun.log 2>&1

# 3. NPT Equilibration
cd "${LAMBDA_PATH}/3_NPT"
gmx grompp -f npt.mdp \
           -c "../2_NVT/nvt.gro" \
           -t "../2_NVT/nvt.cpt" \
           -p "${INPUT_DIR}/topol.top" \
           -o npt.tpr \
           -po npt_out.mdp \
           -maxwarn 1 > npt_grompp.log 2>&1

gmx mdrun -v -deffnm npt \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} > npt_mdrun.log 2>&1

# 4. Production MD
cd "${WORK_DIR}"
gmx grompp -f prod.mdp \
           -c "../3_NPT/npt.gro" \
           -t "../3_NPT/npt.cpt" \
           -p "${INPUT_DIR}/topol.top" \
           -o prod.tpr \
           -po prod_out.mdp \
           -maxwarn 1 > prod_grompp.log 2>&1

gmx mdrun -v -deffnm prod \
          -ntomp ${SLURM_CPUS_PER_TASK} \
          -ntmpi ${SLURM_NTASKS} \
          -dhdl dhdl > prod_mdrun.log 2>&1

# =====================================================================
# Create Analysis Script (Only first task)
# =====================================================================
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    cat > fep_analysis_slurm.sh << 'ANALYSIS_EOF'
#!/bin/bash
#SBATCH --job-name=fep_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/fep_analysis_%j.out

# Inherit environment variables
source $SLURM_SUBMIT_DIR/$0

# Load module
module load gromacs/2023.3

# Create analysis directory
mkdir -p "${ANALYSIS_OUT}"

echo "Starting FEP analysis"

# Process each temperature separately
for TEMP in 290 300 310; do
    TEMP_DIR="${FEP_OUT_BASE}/${TEMP}K"
    OUTPUT_PREFIX="${ANALYSIS_OUT}/fep_${TEMP}K"
    
    # Collect all dhdl.xvg files for this temperature
    DHDL_FILES=()
    for LAMBDA in {00..20}; do
        DHDL_FILE="${TEMP_DIR}/L_${LAMBDA}/4_PROD/dhdl.xvg"
        if [ -f "$DHDL_FILE" ]; then
            DHDL_FILES+=("$DHDL_FILE")
        else
            echo "WARNING: Missing $DHDL_FILE"
        fi
    done

    # Run BAR analysis if we have files
    if [ ${#DHDL_FILES[@]} -gt 0 ]; then
        echo "Analyzing ${#DHDL_FILES[@]} lambda windows for ${TEMP}K"
        gmx bar -f "${DHDL_FILES[@]}" \
                -o "${OUTPUT_PREFIX}.xvg" \
                -oi "${OUTPUT_PREFIX}_int.xvg" \
                -oh "${OUTPUT_PREFIX}_hist.xvg" > "${OUTPUT_PREFIX}.log" 2>&1
    else
        echo "ERROR: No dhdl files found for ${TEMP}K"
    fi
done

# Create summary report
echo "FEP Analysis Summary" > "${ANALYSIS_OUT}/summary.txt"
for TEMP in 290 300 310; do
    OUTPUT_PREFIX="${ANALYSIS_OUT}/fep_${TEMP}K"
    if [ -f "${OUTPUT_PREFIX}.xvg" ]; then
        echo -e "\n===== ${TEMP}K =====" >> "${ANALYSIS_OUT}/summary.txt"
        # Extract free energy summary
        grep "Free energy" "${OUTPUT_PREFIX}.log" >> "${ANALYSIS_OUT}/summary.txt"
        # Extract final ΔG value
        tail -n 1 "${OUTPUT_PREFIX}.xvg" | awk '{print "Final ΔG = " $2 " kJ/mol"}' >> "${ANALYSIS_OUT}/summary.txt"
    fi
done

echo "FEP analysis completed. Results in ${ANALYSIS_OUT}"
ANALYSIS_EOF

    # Make executable and submit with dependency
    chmod +x fep_analysis_slurm.sh
    ANALYSIS_JOBID=$(sbatch --parsable --dependency=afterany:$SLURM_ARRAY_JOB_ID fep_analysis_slurm.sh)
    echo "Submitted FEP analysis job ID: $ANALYSIS_JOBID with dependency on array: $SLURM_ARRAY_JOB_ID"
fi

echo "Completed simulation for Temp=${CURRENT_TEMP}K Lambda=${LAMBDA_DIR}"