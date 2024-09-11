#!/bin/bash
#SBATCH --job-name=528_queue
#SBATCH --ntasks-per-node=32
#SBATCH -N 2
#SBATCH --time=72:00:00
#SBATCH --partition=528_queue
#SBATCH --output IBFE3D.out

mpirun ./main3d input3d
