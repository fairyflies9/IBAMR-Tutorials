#!/bin/bash
#SBATCH --job-name=IBAMR-test
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=windfall
#SBATCG --output IBAMR2D.out

srun ./main2d input2d
