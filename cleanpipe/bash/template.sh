#!/bin/bash
#SBATCH --job-name=myBeautifulRun

#SBATCH --nodes=1  
#SBATCH --exclude=node-15

#SBATCH --cpus-per-task=40          #cpus/cores per task xxx is it the same?
#SBATCH --gres=gpu:1
##SBATCH --mem-per-cpu=1GB

#SBATCH   -o myBeautifulRun.sdout    # Output file
#SBATCH   -e myBeautifulRun.sderr    # Output file for errors
#SBATCH --output=log.out_${i}  #xxx I should remove that and split into out and err as above

#SBATCH --partition=calcul
xxxSBATCH -A whatsthenameofmyaccount?  #account name

xxxSBATCH -qos=expedite
xxxSBATCH -t 86400

#SBATCH -d 100,101,102   # list of job dependencies



module load gromacs/2023
gmx mdrun -v -deffnm $ROOTNAME.$cycle -ntomp 16 -ntmpi 1


#!/bin/bash
#MSUB   -r myBeautifulRun          # Job name
#MSUB   -W yes                     # Let multiple jobs sharing same name & user run simultaneously


# #MSUB -N 1                       # Number of nodes to allocate (Inferred) xxx the site says -N is for job name (not -r), and here it should be #MSUB -l nodes=1
#MSUB   -n 40                      # Number of tasks in parallel mode
#MSUB   -c 1                       # Number of cores per parallel task



#MSUB   -o myBeautifulRun.sdout    # Output file
#MSUB   -e myBeautifulRun.sderr    # Output file for errors

#MSUB   -q rome                    # Partition:    rome    
#MSUB   -A gen13458                # Account/Project code: gen10138 or spe00017
#MSUB   -m scratch,work,store      # File system:  scratch,work,store   xxx site says its a mail message???

xxxMSUB -q queue   #site says this is the wat moab has this instead of partition

#MSUB   -Q normal                  # Quality of Service (test,normal,long) (ccc_mqinfo). xxx site says #MSUB -l qos=expedite
#MSUB   -T 86400                   # Maximum walltime in seconds. xxx site says  #MSUB -l walltime=86400

#MSUB -l depend=100,101,102   # list of job dependencies

ml purge
module load gnu/11 mpi/openmpi/4 gromacs/2023.2
ccc_mprun gmx_mpi mdrun -nice 0 -s $ROOTNAME -deffnm $ROOTNAME.$cycle -v -stepout 1000 -maxh $walltime -cpi $checkpoint -noappend -dd 3 3 3 -npme 13 -dlb yes >& $ROOTNAME.$cycle.runout 
#ccc_msub job.moab


#cluster
# partition (means a specific set of the nodes in the cluster, with different users and queues. this word don have anything to do with memory partition!)
#  nodes
#   cpus
#     cores (made to recieve just 1 paralelizale chunk of my workload. so 1 palalelizable chunk= 1 core.  the program do that, maybe obeying totals set by the user. here is a concrete example: gromacs, if I set -ntmpi 2, groacs will probably do domain decomposition and send two halves as a paralelizable chunk. those papalel chuncks created by the program using the mpi protocol are called mpi ranks)
#            (but then someone invented a way to do fake paralelization called openMP. this is a algorithm that will split the workload and distribuite into the available cores, but being able to send several to the same core, looking like paralelization, but its just time switching. this is effictient inside a single node)
#           (if you are creative anought to want to use both. a good plan is to let the program split the total workload among the nodes (using mpi), and after that, openMP wil further split the chunk in a node into more subchunks. just make sure to match the total fake cores. but you should also test letting mpi splitting into less cores, and let openmp do some intercore work. this might be better depending on the case )
#   gpus (the program, ex: groamcs, will use cuda code to paralilize the workload)

#so to summarize. to paralelize is a good idea to leverage the spread out nature of clusters. but who does the paralilization? if its amongst gpus, its a cuda code somewere. if its amongst cpu1-the program, by itself, using a comunication protocol called mpi, or 2-openMP algorithm, or 3-the program by itself using cuda code
# here are a possible example on how to device total values for the program to then decide how to divide the workload: -nmpi will set the part done by gromacs (that used mpi to comunicate) & -tomp will set the part done by openmp & /-gputasks/-nb gpu/-bonded gpu/-pme gpu -npme 1/-update gpu will set the part done by gromacs using cuda code to talk to the hardware (this should be able to simulate the entire viral capsid)

#a scheduler, like slurm, assign the cores to me, so that roshan will not use them at the same time.  of course the number of assigned cores (gpus and cpus) have to match the paralelization previsouly decided 