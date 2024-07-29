#!/bin/bash
#SBATCH --job-name=IBAMR-test
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --output=IBAMR.out
#SBATCH --account=lauram9

srun ./main2d input2d
