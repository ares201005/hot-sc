#!/bin/sh

#/bin/rm -rf win*

step=5
for i in {3..8}
do
  fd=win$i
  lower=$(echo "($i-3)*$step+1+80"|bc -l)
  upper=$(echo "($i-2)*$step+80"|bc -l)
  interval=$(echo "1"|bc -l)
  
  mkdir $fd
  cp scan.iv $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.iv
  sed -i -e "s/upper/$upper/g" scan.iv
  sed -i -e "s/interval/$interval/g" scan.iv
  
  ./scan.iv >err.out 2>err.out &

  cd ..
done
