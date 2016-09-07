#!/bin/sh

rm eff.dat

#for fd in 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.2 2.4 2.6 2.8 3.0
#do
#  echo $fd
#  cd $fd
#   ./collect.sh
#  cd ..
#done


for fd in 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0
do
  cd $fd
   ../../efficiency/a.out flux-freq.dat >tmpfile
   iqe=`tail -n 1 tmpfile | cut -c48-`
   echo $fd   $iqe >> ../eff.dat
   rm tmpfile fort.60
  cd ..
done
