#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -t 23:59:58
#SBATCH --job-name=ppv-topol
#SBATCH -e stderr-name
#SBATCH -o stdout-name
#SBATCH  --gres=gpu:k40:1
#SBATCH  -p k40
#SBATCH --mem=20g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=zajac028@umn.edu

module purge
module use /home/sarupria/shared/software/ModuleFiles/modules/linux-centos7-haswell
module load gromacs/2021.3-gcc/8.2.0-nompi-openmp-cuda10_2

source pdb2gmx.sh 
