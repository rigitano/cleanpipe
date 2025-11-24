# usage example
# ./runBENCHMARKrome Helix_bench.tpr

# then a folder will be created containing all the tries


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Define the pair of keywords to be replaced in the (.moab) file. each pair will be a specific try. 
#           1 | 2   3   4   5   6 | 7   8   9   10   11 | 12  13   14 | 15  16  17 | 18 19 |
keywords1=("1" "1" "1" "1" "1" "1" "2" "4" "8" "16" "32" "8" "16" "32" "2" "2" "2" "8" "4" "9" "10" "12" "14" "15" "16" "25" "26" "27" "37" "40") # mpi cores
keywords2=("1" "2" "4" "8" "16" "32" "1" "1" "1" "1" "1" "2" "2" "2" "8" "16" "32" "4" "8" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1") # ompi





keywords3=(\
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
	"" \
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
len2=${#keywords2[@]}
len3=${#keywords3[@]}

if [[ $len1 -eq $len2 && $len2 -eq $len3 ]]; then
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
    cat > "try$((i+1))/job.moab" <<EOF
#!/bin/bash

#  @: Editable Variable
#     Please specify the input according to your needs
#     The rest should be fine
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CLUSTER SETTINGS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#MSUB   -r ${TPR_WITHOUT_EXTENSION}-try$((i+1))            # Job name
#MSUB   -n ${keywords1[i]}         # Number of tasks in parallel mode
#MSUB   -c 1                       # Number of cores per parallel task
# #MSUB -N 1                       # Number of nodes to allocate (Inferred)
#MSUB   -W yes                     # Let multiple jobs sharing same name & user run simultaneously
#MSUB   -o GMX.job.IR.output.%I    # Output file
#MSUB   -e GMX.job.IR.outerr.%I    # Output file for errors
#MSUB   -q rome                    # Partition:    rome        
#MSUB   -A gen13458                # Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store      # File system:  scratch,work,store
#MSUB   -Q normal                  # Quality of Service (test,normal,long) (ccc_mqinfo)
#MSUB   -T 3540                    # Maximum walltime in seconds

#@@@@@@@@@@@@@@@@@@@@ Set variables for naming the files @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
GMX=gmx_mpi
# Rootname for gromacs files
ROOTNAME=$TPR_WITHOUT_EXTENSION
# Additional options to use in mdrun
MDRUN_OPT="${keywords3[i]}"
# The file name of the submition file (this file)
BATCH_FNAME=job.moab #[IreneRome]
# Date of run
DATE=\$(date | awk '{print \$1 \$2 \$3 \$6}')
# Time to count for the post-MD checks (in seconds)
CHECK_DURATION=60
set -x #[IreneRome] echo launched commands
##################### load gromacs and set variables  ############################
############# in principle, no changes needed beyond this point ##################

module purge  # retire tous les modules déchargeables de l'environnement
module load gnu/11 # charge gnu/11 et définit gnu/11 comme compilateur dans votre environnement
module load nvhpc/24.3 # besoin de mettre avant OpenMPI comme ce dernier charge un cuda qui n'est pas compatible avec nvhpc/24.3
module load mpi/openmpi/4 # charge la souche OpenMPI
module load gromacs/2025.0 # charge le produit



# Note that sometimes it is advisable especially for small simulations to use more open MP
# threads and less MPI ranks. However, for large systems this appears to be the most
# efficent and reasonably fast setting.
# Note that in any case the number of OpenMP threads times the MPI ranks has to be equal
# the number of available CPUS or cores unless hyperthreading is used, which GROMACS can
# do.

### How many CPU were asked?
ncores=\${SLURM_NTASKS} # number of MPI ranks = number of cores
OMP_NUM_THREADS=${keywords2[i]}      # number of OpenMP threads = 1
echo "-------------CORES ----------------"
echo \${ncores}
echo "-----------------------------------"

export I_MPI_PIN_CELL=core
### export SLURM_CPU_BIND=none #[IreneRome]
export I_MPI_PIN_DOMAIN=auto
export WORKDIR=\$(readlink -f \$SLURM_SUBMIT_DIR)

MDRUN="ccc_mprun gmx_mpi mdrun" #[IreneRome]

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
    # Launch the parallel job
    \$MDRUN -nice 0 -s \$ROOTNAME -deffnm \$ROOTNAME.\$cycle -v \\
           -stepout 1000 -maxh \$walltime -cpi \$checkpoint -noappend \\
           \$MDRUN_OPT -nsteps 50000 >& \$ROOTNAME.\$cycle.runout 

    # Check if the trajectory and the energy file are not corrupted. We will
    # not relaunch a run if a corruption occured.
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




######################################################################################################################
# the run script is over. now its time to pick the performance results in the log file and put it on the folder name #
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
    ccc_msub "job.moab"

    # Go back to the original directory
    cd .. # now outside try



done

