p=$1
t=$2
r=$3

prefix='fcc'
if [ ! -e ../${prefix}-P-$p-T-$t-R-$r/ ]; then



# check if cubic
cubic=$(grep "$p $t" boxsize-${prefix}.dat | awk -v p=$p '$1==p{cub=(($3-$5)**2.+($3-$7)**2.+($5-$7)**2.)**(1./2.); if(cub > 0.6 ){print 1} else {print 0}}')

if [ $cubic == 0 ]; then

mkdir ../${prefix}-P-$p-T-$t-R-$r/
cp equal-sl-P-$p.data plumed.dat single-npt.sh ../${prefix}-P-$p-T-$t-R-$r
	#echo $p $t "# checked the box is cubic"

b=$(grep "$p $t" boxsize-${prefix}.dat | awk -v p=$p '$1==p{print ($3*$5*$7)**(1./3.)}')
echo $p $t "boxisize" $b

sed -e "s/TEMPERATURE/$t/" -e "s/PRESSURE/$p/" -e "s/BOXSIZE/$b/g" -e "s/RANDOM/$r/" continue.lmp > ../${prefix}-P-$p-T-$t-R-$r/continue.lmp

cd ../${prefix}-P-$p-T-$t-R-$r/

sbatch single-npt.sh

cd ..

else
	echo $p $t "# box is not cubic"
fi


fi
