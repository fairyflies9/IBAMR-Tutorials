#!/bin/bash
#SBATCH --job-name=ibamr-rubberband
#SBATCH --ntasks-per-node=16
#SBATCH --time=2:00:00
#SBATCH --mem=1000
#SBATCH --partition=debug_queue
#SBATCH --output IBFE2D.out

mpirun ./main2d input2d
