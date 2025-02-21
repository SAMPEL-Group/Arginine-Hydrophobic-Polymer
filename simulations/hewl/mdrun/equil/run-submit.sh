prot=hewl
for mol in arg lys glu ; do
for n in 79 157 236 314; do
for rep in 2 3 ; do

sed s/AAAA/$prot/g submit.sh > temp1-submit.sh
sed s/BBBB/$mol/g temp1-submit.sh > temp2-submit.sh
sed s/CCCC/$n/g temp2-submit.sh > temp3-submit.sh
sed s/DDDD/$rep/g temp3-submit.sh > temp4-submit.sh

sbatch temp4-submit.sh
rm temp1-submit.sh temp2-submit.sh temp3-submit.sh temp4-submit.sh

done
done
done
