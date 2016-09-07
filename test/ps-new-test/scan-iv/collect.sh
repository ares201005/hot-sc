#!/bin/sh

echo "freq   currL1     currL2        currR1      currR1" > curr-iv.dat
echo "freq   flux        maximum-curr  IQE     Re[eps]      Im[eps]" > flux-iv.dat
step=70
for i in {1..16}
do
  fd=win$i
  cat $fd/curr-iv.dat >>curr-iv.dat
  cat $fd/flux-iv.dat >>flux-iv.dat
done

