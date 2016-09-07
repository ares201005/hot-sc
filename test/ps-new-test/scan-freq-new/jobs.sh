#!/bin/sh

/bin/rm -rf win*

step=70
for i in {1..8}
do
  fd=win$i
  lower=$(echo "($i-1)*$step+1"|bc -l)
  upper=$(echo "$i*$step"|bc -l)
  
  mkdir $fd
  cp scan.freq $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.freq
  sed -i -e "s/upper/$upper/g" scan.freq
  
  ./scan.freq >err.out 2>err.out &

  cd ..
done
