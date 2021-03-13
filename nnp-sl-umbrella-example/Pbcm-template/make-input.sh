p=$1
t=$2
r=$3

prefix='Pbcm-8'

lx=$(grep "$p $t" boxsize-${prefix}.dat | awk -v p=$p '$1==p{print $3*1.5}')
lz=$(grep "$p $t" boxsize-${prefix}.dat | awk -v p=$p '$1==p{print $7*1.5}')

if [ ! -e ../${prefix}-P-$p-T-$t-R-$r/ ]; then


mkdir ../${prefix}-P-$p-T-$t-R-$r/
cp equal-Pbcm-8-sl-P-$p.data plumed.dat single-npt.sh ../${prefix}-P-$p-T-$t-R-$r

echo $p $t $lx $lz 

sed -e "s/RANDOM/$r/" -e "s/TEMPERATURE/$t/" -e "s/PRESSURE/$p/" -e "s/XBOXSIZEX/$lx/" -e "s/ZBOXSIZEZ/$lz/" continue.lmp > ../${prefix}-P-$p-T-$t-R-$r/continue.lmp

cd ../${prefix}-P-$p-T-$t-R-$r/

sbatch single-npt.sh

cd ..


fi
