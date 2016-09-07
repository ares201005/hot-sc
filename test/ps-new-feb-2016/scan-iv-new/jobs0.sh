#!/bin/sh

#/bin/rm -rf win*

step=40
for i in {1..2}
do
  fd=win$i
  lower=$(echo "($i-1)*$step+1"|bc -l)
  upper=$(echo "$i*$step"|bc -l)
  interval=$(echo "5"|bc -l)
  
  mkdir $fd
  cp scan.iv $fd
  cd $fd

  sed -i -e "s/lower/$lower/g" scan.iv
  sed -i -e "s/upper/$upper/g" scan.iv
  sed -i -e "s/interval/$interval/g" scan.iv
  
  ./scan.iv >err.out 2>err.out &

  cd ..
done
