export OMP_NUM_THREADS=8
gmx mdrun \
 -s $WORKING_DIR/$title-nvt-equil.tpr -cpi $WORKING_DIR/$title-nvt-equil.cpt \
 -deffnm $WORKING_DIR/$title-nvt-equil \
 -ntomp 8 -ntmpi 8
