#!/bin/bash
#SBATCH --job-name=CIB_Examples
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=lauram9
#SBATCH --output=IBAMR2D.output

mpirun ./main2d input2d