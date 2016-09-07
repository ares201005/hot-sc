#!/bin/sh

/bin/rm -rf win*

step=1
for i in {1..80}
do
  fd=win$i
  lower=$(echo "($i-1)*$step+1"|bc -l)
  upper=$(echo "$i*$step"|bc -l)
  
  mkdir $fd
  cp scan.mscoup $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.mscoup
  sed -i -e "s/upper/$upper/g" scan.mscoup
  
  qsub scan.mscoup
  #./scan.mscoup >err.out 2>err.out &

  cd ..
done
