#!/bin/sh

/bin/rm -rf win*

step=18
for i in {1..32}
do
  fd=win$i
  lower=$(echo "($i-1)*$step+1"|bc -l)
  upper=$(echo "$i*$step"|bc -l)
  
  mkdir $fd
  cp scan.freq $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.freq
  sed -i -e "s/upper/$upper/g" scan.freq
  
  qsub scan.freq
  #./scan.freq >err.out 2>err.out &

  cd ..
done
