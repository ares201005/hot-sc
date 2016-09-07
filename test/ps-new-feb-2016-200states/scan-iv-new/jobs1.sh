#!/bin/sh

#/bin/rm -rf win*

step=1
for i in {9..35}
do
  fd=win$i
  lower=$(echo "($i-9)*$step+1+80"|bc -l)
  upper=$(echo "($i-8)*$step+80"|bc -l)
  interval=$(echo "1"|bc -l)
  
  mkdir $fd
  cp scan.iv $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.iv
  sed -i -e "s/upper/$upper/g" scan.iv
  sed -i -e "s/interval/$interval/g" scan.iv
  
  qsub scan.iv
  #./scan.iv >err.out 2>err.out &

  cd ..
done
