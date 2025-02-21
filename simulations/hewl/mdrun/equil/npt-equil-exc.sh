gmx grompp \
 -f ../mdp/npt-equil-exc.mdp \
 -n ../ndx/$title.ndx \
 -c ./$title-nvt-equil.gro \
 -r ./$title-nvt-equil.gro \
 -p ../topol/$top.top \
 -o ./$title-$rep-npt-equil.tpr

export OMP_NUM_THREADS=4
gmx mdrun \
 -s ./$title-$rep-npt-equil.tpr \
 -v -deffnm ./$title-$rep-npt-equil \
 -nb gpu -ntomp 4 -ntmpi 2
