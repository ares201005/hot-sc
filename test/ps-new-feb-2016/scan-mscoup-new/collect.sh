#!/bin/sh

echo "V_ms  relax     curr1           currL2            currR1           currR2" > curr-relax-coup.dat
echo "V_ms relax       flux          maximum-curr       IQE       Re[eps]      Im[eps]" > flux-relax-coup.dat
for i in {1..8}
do
  fd=win$i
  cat $fd/curr-relax-coup.dat >>curr-relax-coup.dat
  cat $fd/flux-relax-coup.dat >>flux-relax-coup.dat
done

