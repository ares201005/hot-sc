#!/bin/sh

echo "freq   currL1     currL2        currR1      currR1" > curr-freq.dat
echo "freq   flux        maximum-curr  IQE     Re[eps]      Im[eps]" > flux-freq.dat
for i in {1..80}
do
  fd=win$i
  cat $fd/curr-freq.dat >>curr-freq.dat
  cat $fd/flux.dat >>flux-freq.dat
done

