gmx grompp \
 -f ../mdp/npt-equil.mdp \
 -c ./$title-nvt-equil.gro \
 -r ./$title-nvt-equil.gro \
 -p ../topol/$top.top \
 -o ./$title-npt-equil.tpr

export OMP_NUM_THREADS=6
gmx mdrun \
 -s ./$title-npt-equil.tpr \
 -v -deffnm ./$title-npt-equil \
 -nb gpu -ntomp 6 -ntmpi 4

gmx trjconv \
 -f ./$title-npt-equil.xtc \
 -s ./$title-npt-equil.tpr \
 -o ./$title-npt-equil.xtc -pbc whole << EOF
0
0
EOF

