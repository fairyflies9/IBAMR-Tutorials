#!/bin/bash
#SBATCH --job-name=ib2d
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=04:00:00
#SBATCH --partition=windfall

module load matlab
matlab -nodisplay -nosplash < main2d.m
