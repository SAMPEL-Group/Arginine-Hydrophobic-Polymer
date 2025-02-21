gmx mdrun \
 -s ./$title-npt-equil.tpr -cpi $title-npt-equil.cpt \
 -v -deffnm ./$title-npt-equil

gmx trjconv \
 -f ./$title-npt-equil.xtc \
 -s ./$title-npt-equil.tpr \
 -o ./$title-npt-equil.xtc -pbc whole << EOF
0
0
EOF
