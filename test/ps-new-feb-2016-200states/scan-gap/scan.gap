#!/bin/sh


cat > in.relax <<FIN
\$general lphoton=.false., tplasmon=.true., lscba=.false. \$end
\$temperature telec=300.0d0 \$end
\$land ngrid=10, downrange=0.5d0, uprange=0.5d0 \$end
\$photon mode=1,phfreq=3.60, \$end
FIN


file=curr-gap.dat
file2=flux-gap.dat
rm $file $file2

echo "gap  volt       curr1           curr2             curr3            curr4" > $file
echo "gap  volt       flux            maxcurr           IQE        eps[Re]      eps[Im]" > $file2

a=0.02
b=0.02
for i in {lower..upper}
do
 for j in {0..55}
 do

 c=$(echo "$i*$a"|bc -l)
 gap=$(echo "2*$i*$a"|bc -l)

 d=$(echo "$j*$a"|bc -l)
 volt=$(echo "2*$j*$a"|bc -l)

 echo "\$voltage ampG=0.d0, ampl=-$d, ampr=$d \$end" >in.2
 echo "\$system nLead=2, norbs=2, onsite=$c, hop=1.0d0, elcoup=1.0d0 miu0=0.d0, half_w(1)=1.d3, half_w(2)=1.d3 \$end" >> in.2
 echo "\$plasmon t_read_file=.false., barrier_e=$c, barrier_h=$c,eps_coup=2.0d-3,mscoup=1.d0,nst=100,eig(1)=-4.0,eig(2)=4.d0,relax=0.01, \$end" >> in.2

 cat in.relax in.2 > new.td

 ../../../../bin/negf-pt <new.td >out 

 curr1=$(grep ' current  through lead:  1 is' out |cut -c32-47)
 curr2=$(grep ' current2 through lead:  1 is' out |cut -c32-47)
 curr3=$(grep ' current  through lead:  2 is' out |cut -c32-47)
 curr4=$(grep ' current2 through lead:  2 is' out |cut -c32-47)

 echo "$gap   $volt $curr1  $curr2  $curr3  $curr4" >> $file

 flux=$(grep ' absorped photon flux is' out |cut -c27-42)
 maxcurr=$(grep ' maximum photocurrent is' out |cut -c27-42)
 IQE=$(grep ' IQE is  ' out |cut -c29-36)
 eps=$(grep -A 1 'eps+2' out |tail -n 1 |cut -c62-)  # (eps-1)/(eps+2) - 1

 echo "$gap    $volt $flux  $maxcurr  $IQE $eps " >> $file2

 done
 echo "" >> $file
 echo "" >> $file2
done
