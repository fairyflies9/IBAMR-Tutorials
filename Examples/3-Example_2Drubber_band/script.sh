#!/bin/bash
#SBATCH --job-name=ibamr-rubberband
#SBATCH --ntasks-per-node=32
#SBATCH -N 2
#SBATCH --time=12:00:00
#SBATCH --mem=1000
#SBATCH --partition=528_queue
#SBATCH --output ibamr2d.out

mpirun ./main2d input2d
