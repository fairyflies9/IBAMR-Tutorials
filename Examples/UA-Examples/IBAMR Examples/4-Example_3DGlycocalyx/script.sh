#!/bin/bash
#SBATCH --job-name=IBAMR-test
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=lauram9
#SBATCH --output IBAMR2D.out

mpirun ./main3d input3d
