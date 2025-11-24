#!/bin/bash

#SBATCH --partition=calcul
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --job-name=oneline
#SBATCH --output=oneline.scheduler.out.and.err
#SBATCH --exclude=node-15

module purge
module load gromacs/2024.5


gmx mdrun -v -s *prod*.tpr -ntomp 32 -ntmpi 1 -deffnm again
