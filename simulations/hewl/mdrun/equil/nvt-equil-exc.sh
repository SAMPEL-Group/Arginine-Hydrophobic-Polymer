gmx grompp \
 -f ../mdp/nvt-equil-exc.mdp \
 -n ../ndx/$title.ndx \
 -c ../emin/$title-emin.gro \
 -r ../emin/$title-emin.gro \
 -p ../topol/$top.top \
 -o ./$title-$rep-nvt-equil.tpr

export OMP_NUM_THREADS=4
gmx mdrun \
 -s ./$title-$rep-nvt-equil.tpr \
 -v -deffnm ./$title-$rep-nvt-equil\
 -nb gpu -ntomp 4 -ntmpi 2
