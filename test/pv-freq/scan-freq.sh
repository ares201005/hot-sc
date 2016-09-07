#!/bin/sh

file=curr-freq.dat
rm $file

a=0.05
for i in {0..100}
do 
 c=$(echo "$i*$a"|bc -l)
 echo "phfreq=$c, \$end" > in.2
 cat in.freq in.2 > new.td

 ../../bin/negf-pt <new.td >out 

 curr1=$(grep ' current  through lead:  1 is' out |cut -c32-47)
 curr2=$(grep ' current2 through lead:  1 is' out |cut -c32-47)
 curr3=$(grep ' current  through lead:  2 is' out |cut -c32-47)
 curr4=$(grep ' current2 through lead:  2 is' out |cut -c32-47)

 echo "$c $curr1  $curr2  $curr3  $curr4" >> $file

done
