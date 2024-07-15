#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks=128
#SBATCH -t 23:59:59
#SBATCH --job-name=ppv-5fold-prod
#SBATCH -e stderr-name
#SBATCH -o stdout-name
#SBATCH -p a100-8
#SBATCH --gres=gpu:a100:8
#SBATCH --mem=0
#SBATCH --mail-type=all
#SBATCH --mail-user=zajac028@umn.edu

#module purge
source ~/software/load_scripts/load_gromacs-2023.2-rocky8.sh
#source ~/software/load_scripts/load_gromacs-2021.4.sh
#source ~/software/load_scripts/load_plumed-2.8.0.sh
WORKING_DIR=/scratch.global/zajac028/viruses/ppv/prod

title=ppv-5fold-396arg-outer

source nvt-prod.sh
#source nvt-prod-reduced-restraints.sh
#source nvt-prod.sh
#source nvt-prod-cpu.sh 
#source nvt-prod-cpu-from-restraints.sh
#source continue.sh
#source continue-reduced-restraints.sh
#source nvt-prod-wat-dyn.sh
