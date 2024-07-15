mpirun -np 1 gmx_mpi grompp \
 -f ../mdp/nvt-prod.mdp \
 -c ../equil/$title-nvt-equil.gro \
 -r ../equil/$title-nvt-equil.gro \
 -n ../ndx/$title.ndx \
 -p ../topol/$top.top \
 -o $WORKING_DIR/$title-nvt-prod.tpr

mpirun gmx_mpi mdrun \
 -s $WORKING_DIR/$title-nvt-prod.tpr \
 -deffnm $WORKING_DIR/$title-nvt-prod
