#!/bin/sh

for i in {1..1000}
do
 echo $i
 pkill negf-pt
done
