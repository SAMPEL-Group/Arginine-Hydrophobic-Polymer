WORKING_DIR=/scratch.global/zajac028/hewl-bsa/vagenende-2013

export OMP_NUM_THREADS=6
gmx mdrun \
 -cpi $WORKING_DIR/$title-npt-prod.cpt \
 -deffnm $WORKING_DIR/$title-npt-prod -v \
 -nb gpu -ntomp 6 -ntmpi 4
