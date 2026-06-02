# usage example
# ./runBENCHMARKrome Helix_bench.tpr

# then a folder will be created containing all the tries


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Define the pair of keywords to be replaced in the (.moab) file. each pair will be a specific try. 
#           1 | 2   3   4   5   6 | 7   8   9   10   11 | 12  13   14 | 15  16  17 | 18 19 |
keywords1=("1" "2" "4" "8" "16" "32" "1" "2" "4" "8" "16" "32" "2" "4" "8" "16" "32" "16" "32" "64" "9" "10" "12" "14" "15" "16" "25" "26" "27" "37" "40") # mpi cores
#keywords2=("1" "2" "4" "8" "16" "32" "1" "1" "1" "1" "1" "2" "2" "2" "8" "16" "32" "4" "8" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1") # ompi





keywords3=(\
	"-nt 1" \
	"-nt 2" \
	"-nt 4" \
	"-nt 8" \
	"-nt 16" \
	"-nt 32" \
	"-ntmpi 1 -ntomp 1" \
	"-ntmpi 1 -ntomp 2" \
	"-ntmpi 1 -ntomp 4" \
	"-ntmpi 1 -ntomp 8" \
	"-ntmpi 1 -ntomp 16" \
	"-ntmpi 1 -ntomp 32" \
	"-ntmpi 2 -ntomp 1" \
	"-ntmpi 4 -ntomp 1" \
	"-ntmpi 8 -ntomp 1" \
	"-ntmpi 16 -ntomp 1" \
	"-ntmpi 32 -ntomp 1" \
	"-ntmpi 8 -ntomp 2" \
	"-ntmpi 16 -ntomp 2" \
	"-ntmpi 32 -ntomp 2" \
	"-dd 2 2 2 -npme 1 -dlb yes" \
	"-dd 2 2 2 -npme 2 -dlb yes" \
	"-dd 2 2 2 -npme 4 -dlb yes" \
	"-dd 2 2 3 -npme 2 -dlb yes" \
	"-dd 2 2 3 -npme 3 -dlb yes" \
	"-dd 2 2 3 -npme 4 -dlb yes" \
	"-dd 2 3 3 -npme 7 -dlb yes" \
	"-dd 2 3 3 -npme 8 -dlb yes" \
	"-dd 2 3 3 -npme 9 -dlb yes" \
	"-dd 3 3 3 -npme 10 -dlb yes" \
	"-dd 3 3 3 -npme 13 -dlb yes")




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
#SBATCH --cpus-per-task=${keywords1[i]}
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --job-name=${TPR_WITHOUT_EXTENSION}-try$((i+1))
#SBATCH --output=scheduler.out.and.err
##SBATCH --exclude=node-15


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

