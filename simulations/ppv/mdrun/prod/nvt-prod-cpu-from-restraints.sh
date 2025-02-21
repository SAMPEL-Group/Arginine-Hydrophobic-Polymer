mpirun -np 1 gmx_mpi grompp \
 -f ../mdp/nvt-prod-nores.mdp \
 -c ./surface-restraints/$title-nvt-prod-res.gro \
 -r ./surface-restraints/$title-nvt-prod-res.gro \
 -n ../ndx/$title.ndx \
 -p ../topol/$top.top \
 -o $WORKING_DIR/$title-nvt-prod.tpr

mpirun gmx_mpi mdrun \
 -s $WORKING_DIR/$title-nvt-prod.tpr \
 -deffnm $WORKING_DIR/$title-nvt-prod
