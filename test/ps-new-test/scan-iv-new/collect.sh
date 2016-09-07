#!/bin/sh

echo "voltage   currL1     currL2        currR1      currR1" > curr-iv.dat
echo "voltage   flux        maximum-curr  IQE     Re[eps]      Im[eps]" > flux-iv.dat
for i in {1..8}
do
  fd=win$i
  cat $fd/curr-iv.dat >>curr-iv.dat
  cat $fd/flux-iv.dat >>flux-iv.dat
done

