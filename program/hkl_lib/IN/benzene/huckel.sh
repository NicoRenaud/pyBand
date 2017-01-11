#!/bin/bash

#export OMP_NUM_THREADS=2
#export MKL_NUM_THREADS=$OMP_NUM_THREADS



IN=./benzene.in
OUT=./

EXE=../../huckel

date1=$(date +"%s")
$EXE $IN $OUT
date2=$(date +"%s")
diff=$(($date2-$date1))
echo "Done in $(($diff / 60)) minutes and $(($diff % 60)) seconds"

#cd $OUT
#gnuplot plotDOS.gnu
#cd ../

