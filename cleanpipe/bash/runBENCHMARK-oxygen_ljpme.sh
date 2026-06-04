# usage example
# ./runBENCHMARKrome Helix_bench.tpr

# then a folder will be created containing all the tries


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



# xxx tell if there is LJ-PME on the system

#HAS_LJPME=false
HAS_LJPME=true


#based on that detection, define the proper set of keywords:


if [ "$HAS_LJPME" = true ]; then

    # ── LJ-PME system: PME must stay on CPU ──────────────────────────────
    keywords1=(
        "1"  "2"  "4"  "8"  "16" "20" "40"   # ntmpi 1, vary ntomp
        "2"  "4"  "8"  "16" "20"              # ntmpi 2, vary ntomp
        "4"  "8"  "16" "20"                   # ntmpi 4, vary ntomp
        "8"  "16" "40"                        # ntmpi 8, vary ntomp
    )
    keywords3=(
        # ntmpi 1, vary ntomp — nb on GPU, PME forced to CPU
        "-ntmpi 1 -ntomp 1  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 1 -ntomp 2  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 1 -ntomp 4  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 1 -ntomp 8  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 1 -ntomp 16 -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 1 -ntomp 20 -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 1 -ntomp 40 -gpu_id 0 -nb gpu -pme cpu -update gpu"
        # ntmpi 2, vary ntomp
        "-ntmpi 2 -ntomp 1  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 2 -ntomp 2  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 2 -ntomp 4  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 2 -ntomp 8  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 2 -ntomp 10 -gpu_id 0 -nb gpu -pme cpu -update gpu"
        # ntmpi 4, vary ntomp
        "-ntmpi 4 -ntomp 1  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 4 -ntomp 2  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 4 -ntomp 4  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 4 -ntomp 5  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        # ntmpi 8, vary ntomp
        "-ntmpi 8 -ntomp 1  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 8 -ntomp 2  -gpu_id 0 -nb gpu -pme cpu -update gpu"
        "-ntmpi 8 -ntomp 5  -gpu_id 0 -nb gpu -pme cpu -update gpu"
    )

else

    # ── Standard system: full GPU offload available ───────────────────────
    keywords1=(
        "1"  "2"  "4"  "8"  "16" "20" "40"   # ntmpi 1, vary ntomp
        "2"  "4"  "8"  "16" "20"              # ntmpi 2, vary ntomp
        "4"  "8"  "16" "20"                   # ntmpi 4, vary ntomp
        "8"  "16" "40"                        # ntmpi 8, vary ntomp
        # GPU offload strategy comparison (fixed threading)
        "20" "20" "20"
    )
    keywords3=(
        # ntmpi 1, vary ntomp — full GPU offload
        "-ntmpi 1 -ntomp 1  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 1 -ntomp 2  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 1 -ntomp 4  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 1 -ntomp 8  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 1 -ntomp 16 -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 1 -ntomp 20 -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 1 -ntomp 40 -gpu_id 0 -nb gpu -pme gpu -update gpu"
        # ntmpi 2, vary ntomp
        "-ntmpi 2 -ntomp 1  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 2 -ntomp 2  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 2 -ntomp 4  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 2 -ntomp 8  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 2 -ntomp 10 -gpu_id 0 -nb gpu -pme gpu -update gpu"
        # ntmpi 4, vary ntomp
        "-ntmpi 4 -ntomp 1  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 4 -ntomp 2  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 4 -ntomp 4  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 4 -ntomp 5  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        # ntmpi 8, vary ntomp
        "-ntmpi 8 -ntomp 1  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 8 -ntomp 2  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        "-ntmpi 8 -ntomp 5  -gpu_id 0 -nb gpu -pme gpu -update gpu"
        # GPU offload strategy comparison at fixed ntmpi 2 ntomp 10
        "-ntmpi 2 -ntomp 10 -gpu_id 0 -nb gpu -pme cpu -update cpu"
        "-ntmpi 2 -ntomp 10 -gpu_id 0 -nb gpu -pme gpu -update cpu"
        "-ntmpi 2 -ntomp 10 -gpu_id 0 -nb gpu -pme gpu -update gpu"
    )

fi




# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# obtain the size of the keywords, and exit if they are not the same
len1=${#keywords1[@]}
#len2=${#keywords2[@]}
len3=${#keywords3[@]}

if [[ $len1 -eq $len3 ]]; then
	echo "keywords lenght is ok"
else
	echo "the length of the keywords lists are not the same!"
	exit
fi



# obtain the paramenter that contains the full tpr name
TPR=$1

# Extract the file name without the extension
TPR_WITHOUT_EXTENSION="${TPR%.*}"

# Create the target directory that will contain all tries of the current benchmark
TARGET_DIR="${TPR_WITHOUT_EXTENSION}__runBENCHMARK"
mkdir -p "${TARGET_DIR}"
cd ${TARGET_DIR}




# Go throught the pairs of keywords and save a try for each pair
for ((i = 0; i < len1; i++)); do

    # Create a new directory for current iteration
    mkdir -p "try$((i+1))"

    cp -f "../${TPR}" "try$((i+1))"

    
    # Create the file that will later be executed:
    cat > "try$((i+1))/job.sh" <<EOF
#!/bin/bash

#  @: Editable Variable
#     Please specify the input according to your needs
#     The rest should be fine
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CLUSTER SETTINGS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#SBATCH --partition=calcul
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${keywords1[i]}
#SBATCH --gres=gpu:1            
#SBATCH --job-name=${TPR_WITHOUT_EXTENSION}-try$((i+1))
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=01:00:00                  # short walltime — benchmarks don't need long



module purge  # retire tous les modules déchargeables de l'environnement
module load gromacs/2024.5 # charge le produit

gmx mdrun -deffnm ${TPR_WITHOUT_EXTENSION} ${keywords3[i]} -nsteps 50000 > mdrun.out 2> mdrun.err
#####################################################################################################################



    # Check for .log files
    log_file=\$(ls *.log 2> /dev/null)

    if [ ! -z "\$log_file" ]; then
      # Extract the line containing "Performance:" and get the number after it
      performance_line=\$(grep "Performance:" "\$log_file" | awk '{print \$2}')
    else
      performance_line="" # no log file
    fi

    # Check if performance_line is empty
    if [ -n "\$performance_line" ]; then
        # Extract the integer part of the number before the decimal point
        integer_part=\${performance_line%%.*}
    else
        integer_part="fail"
    fi


# Define the suffix to be added to the directory name
SUFFIX="-\${integer_part}"


# Store the path where this script is in
CURRENT_PATH="\$(pwd)"
CURRENT_FOLDER_NAME="try$((i+1))"
cd ..
PARENT_DIR=\$(pwd)
cd "\$CURRENT_FOLDER_NAME"

# Determine the new folder name with suffix
NEW_FOLDER_NAME="\${CURRENT_FOLDER_NAME}\${SUFFIX}"

# Rename the folder
mv -f "\$CURRENT_PATH" "\$PARENT_DIR/\$NEW_FOLDER_NAME"

# Update the working directory
cd "\$PARENT_DIR/\$NEW_FOLDER_NAME"


EOF


    # Change directory to the new folder
    cd "try$((i+1))"

    # Submit the job
    sbatch "job.sh"

    # Go back to the original directory
    cd .. # now outside try



done

