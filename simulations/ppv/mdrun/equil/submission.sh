#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH -t 11:59:58
#SBATCH --job-name=ppv-5f-1132arg-equil
#SBATCH -e stderr-name
#SBATCH -o stdout-name
#SBATCH -p a100-4,a100-8
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=0
#SBATCH --mail-type=all
#SBATCH --mail-user=zajac028@umn.edu

module purge
source ~/software/load_scripts/load_gromacs-2023.2-rocky8.sh
#source ~/software/load_scripts/load_gromacs-2021.4.sh
#source ~/software/load_scripts/load_plumed-2.8.0.sh

WORKING_DIR=/scratch.global/zajac028/viruses/ppv/equil

title=ppv-5fold-1132arg-outer
top=ppv-5fold-1132arg-outer

source nvt-equil.sh
#source nvt-equil-cpu.sh 
#source npt-equil.sh
#source continue.sh
