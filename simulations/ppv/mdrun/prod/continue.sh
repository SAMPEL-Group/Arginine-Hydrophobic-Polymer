export OMP_NUM_THREADS=16
gmx mdrun \
 -s $WORKING_DIR/$title-nvt-prod.tpr -cpi $WORKING_DIR/$title-nvt-prod.cpt \
 -deffnm $WORKING_DIR/$title-nvt-prod -ntomp 16 -ntmpi 4
