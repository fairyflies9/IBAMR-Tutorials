#!/bin/bash
#SBATCH --job-name=3D_Cylinder
#SBATCH --ntasks-per-node=16
#SBATCH --time=2:00:00
#SBATCH --partition=debug_queue
#SBATCH --output IBAMR_3D.out

mpirun ./main3d input3d
