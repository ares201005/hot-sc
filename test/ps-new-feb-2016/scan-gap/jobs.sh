#!/bin/sh

/bin/rm -rf win*

step=4
for i in {3..10}
do
  fd=win$i
  lower=$(echo "($i-1)*$step+1"|bc -l)
  upper=$(echo "$i*$step"|bc -l)
  
  mkdir $fd
  cp scan.gap $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.gap
  sed -i -e "s/upper/$upper/g" scan.gap
  
  ./scan.gap >err.out 2>err.out &

  cd ..
done
