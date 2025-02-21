gmx grompp \
 -f ../mdp/npt-equil.mdp \
 -c $WORKING_DIR/$title-nvt-equil.gro \
 -r $WORKING_DIR/$title-nvt-equil.gro \
 -p ../topol/$top.top \
 -maxwarn 2 \
 -o $WORKING_DIR/$title-npt-equil.tpr

export OMP_NUM_THREADS=6
gmx mdrun \
 -s $WORKING_DIR/$title-npt-equil.tpr \
 -deffnm $WORKING_DIR/$title-npt-equil \
 -ntomp 6 -ntmpi 4
