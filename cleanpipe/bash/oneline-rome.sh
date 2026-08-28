#MSUB   -r 1line            # Job name
#MSUB   -n 4                          # Number of tasks in parallel mode (ntmpi)
#MSUB   -c 1                                 # Number of cores per parallel task
#MSUB   -W yes                               # Let multiple jobs sharing same name & user run simultaneously
#MSUB   -o 1line.%I.scheduler.out      # Output file
#MSUB   -e 1line.%I.scheduler.err      # Output file for errors
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

export OMP_NUM_THREADS=1      # number of OpenMP threads (ntomp)

ccc_mprun gmx_mpi mdrun -v -deffnm a1_inO_prod_298_01
