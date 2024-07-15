gmx grompp \
 -f ../mdp/nvt-equil.mdp \
 -c ../emin/$title-emin.gro \
 -r ../emin/$title-emin.gro \
 -p ../topol/$top.top \
 -o ./$title-nvt-equil.tpr

gmx mdrun \
 -s ./$title-nvt-equil.tpr \
 -v -deffnm ./$title-nvt-equil

gmx trjconv \
 -f ./$title-nvt-equil.xtc \
 -s ./$title-nvt-equil.tpr \
 -o ./$title-nvt-equil.xtc -pbc whole << EOF
0
0
EOF
