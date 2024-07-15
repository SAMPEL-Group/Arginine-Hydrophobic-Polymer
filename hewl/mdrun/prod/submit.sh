#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH -t 23:59:58
#SBATCH --job-name=hewl-repro
#SBATCH -e stderr-name
#SBATCH -o stdout-name
#SBATCH  --gres=gpu:a100:2
#SBATCH  -p a100-4
#SBATCH --mem=40g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=zajac028@umn.edu

module purge
source ~/software/load_scripts/load_gromacs-2023.2.sh

title=hewl-24arg
top=$title

#source npt-prod-exc.sh
source mdrun-continue.sh
