#!/bin/bash
#SBATCH --job-name=first_ibfe_job
#SBATCH --ntasks=16
#SBATCH --time=2:00:00
#SBATCH --partition=debug_queue
#SBATCH --output IBFE2D.out

./main2d input2d

