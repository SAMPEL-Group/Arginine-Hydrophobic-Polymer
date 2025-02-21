mpirun -np 1 gmx_mpi grompp \
 -f ../mdp/nvt-equil-exc.mdp \
 -c ../emin/$title-emin-wall.gro \
 -r ../emin/$title-emin-wall.gro \
 -n ../ndx/$title.ndx \
 -p ../topol/$top.top \
 -o $WORKING_DIR/$title-nvt-equil.tpr

mpirun gmx_mpi mdrun \
 -s $WORKING_DIR/$title-nvt-equil.tpr \
 -deffnm $WORKING_DIR/$title-nvt-equil
