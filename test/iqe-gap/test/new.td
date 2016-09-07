$general lphoton=.false., tplasmon=.true., lscba=.false. $end
$system nLead=2, norbs=2, onsite=0.5d0, hop=1.0d0, elcoup=1.0d0 miu0=0.d0, half_w(1)=1.d3, half_w(2)=1.d3 $end
$temperature telec=300.0d0 $end
$land ngrid=10, downrange=0.5d0, uprange=0.5d0 $end
$plasmon t_read_file=.false., barrier_e=0.5d0, barrier_h=0.5d0,eps_coup=2.d-3,mscoup=1.0d0,nst=400,eig(1)=-4.0,eig(2)=4.d0,relax=0.01d0 $end
$voltage ampG=0.d0, ampl=-0.0d0, ampr=0.0d0 $end
$photon mode=1,phfreq=1.0, $end
