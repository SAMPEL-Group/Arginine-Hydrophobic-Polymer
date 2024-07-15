#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -t 01:59:58
#SBATCH --job-name=hewl-157arg-equil
#SBATCH -e stderr-name
#SBATCH -o stdout-name
#SBATCH  --gres=gpu:k40:1
#SBATCH  -p k40
#SBATCH --mem=20g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=zajac028@umn.edu

source ~/software/load_scripts/load_gromacs-2023.2.sh

title=hewl-157arg
top=$title

source continue.sh 
