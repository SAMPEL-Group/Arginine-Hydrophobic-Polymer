gmx grompp \
 -f ../mdp/nvt-prod-walls.mdp \
 -c ../equil/$title-nvt-equil.gro \
 -r ../equil/$title-nvt-equil.gro \
 -n ../ndx/$title.ndx \
 -p ../topol/$title.top \
 -o $WORKING_DIR/$title-nvt-prod.tpr

export OMP_NUM_THREADS=32
gmx mdrun \
 -s $WORKING_DIR/$title-nvt-prod.tpr \
 -deffnm $WORKING_DIR/$title-nvt-prod \
 -ntomp 32 -ntmpi 4
