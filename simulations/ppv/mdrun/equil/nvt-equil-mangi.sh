gmx grompp \
 -f ../mdp/nvt-equil-walls.mdp \
 -c ../emin/$title-emin-wall.gro \
 -r ../emin/$title-emin-wall.gro \
 -n ../ndx/$title.ndx \
 -p ../topol/$title.top \
 -o $WORKING_DIR/$title-nvt-equil.tpr

export OMP_NUM_THREADS=6
gmx mdrun \
 -s $WORKING_DIR/$title-nvt-equil.tpr \
 -deffnm $WORKING_DIR/$title-nvt-equil -ntomp 6 -ntmpi 4
