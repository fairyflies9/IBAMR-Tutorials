#!/bin/bash
#SBATCH --job-name=CIB-ex
#SBATCH --ntasks=16
#SBATCH --time=2:00:00
#SBATCH --partition=debug_queue
#SBATCH --output ibamr2d.out

mpirun ./main2d input2d
