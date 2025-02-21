title=ppv-5fold
gmx trjconv \
 -f ./$title-nvt-equil.xtc \
 -s ./$title-nvt-equil.tpr \
 -b 0 -e 0 \
 -o ./$title-nvt-equil.gro -center -pbc whole

gmx trjconv \
 -f ./$title-nvt-equil.xtc \
 -s ./$title-nvt-equil.tpr \
 -o ./$title-nvt-equil.xtc -center -pbc cluster
