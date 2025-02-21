title=ppv-full
w=0.0
gmx pdb2gmx -f ../conf/$title-sliced-w$w.pdb -ignh \
	    -o ./conf/$title.gro \
	    -i ./chain-itp/vanilla-posres/$title-posre.itp \
	    -p ./$title.top << EOF
1
2
EOF
