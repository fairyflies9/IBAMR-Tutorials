#!/bin/bash
#SBATCH --job-name=AT3
#SBATCH --ntasks-per-node=94
#SBATCH -N 1
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=lauram9
#SBATCH --output output1

srun ./main2d input2d
