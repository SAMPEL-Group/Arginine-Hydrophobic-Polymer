#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 03:59:58
#SBATCH --job-name=AAAA-CCCCBBBB-DDDD
#SBATCH -e stderr-name
#SBATCH -o stdout-name
#SBATCH  --gres=gpu:v100:1
#SBATCH  -p v100
#SBATCH --mem=20g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=zajac028@umn.edu

source ~/software/load_scripts/load_gromacs-2023.2.sh

rep=DDDD
title=AAAA-CCCCBBBB
top=$title

source nvt-equil-exc.sh 
source npt-equil-exc.sh
