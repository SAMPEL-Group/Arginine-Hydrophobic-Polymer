WORKING_DIR=/scratch.global/zajac028/hewl-bsa/vagenende-2013

gmx grompp \
 -f ./npt-prod-exc.mdp \
 -n ../ndx/$title.ndx \
 -c ../emin/$title-emin.gro \
 -r ../emin/$title-emin.gro \
 -p ../topol/$top.top \
 -maxwarn 1 \
 -o $WORKING_DIR/$title-npt-prod.tpr

export OMP_NUM_THREADS=6
gmx mdrun \
 -s $WORKING_DIR/$title-npt-prod.tpr \
 -v -deffnm $WORKING_DIR/$title-npt-prod \
  -nb gpu -ntomp 6 -ntmpi 4
