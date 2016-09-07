program model
use transmod
use transphonon
implicit none

include '../include/parameter'
include '../include/trans_data'

character*10 :: date,time

real*8, allocatable :: fock0(:,:)
real*8, allocatable :: miu(:)

integer :: ierror(10)=0

logical :: lscba
namelist /general/ tplasmon,lphoton,lscba
integer :: norbs, nLead
real*8 :: onsite, hop, elcoup,half_w(5)
namelist /system/ nLead, norbs, onsite, hop, elcoup,miu0,half_w
namelist /temperature/ telec

real*8 :: bonded, repulsion,diele, length
namelist /ppp/ bonded, repulsion, diele, length
!
namelist /photon/ mode, epcoup, phfreq, pt_telec
namelist /voltage/ ampl, ampr,ampG
!
integer :: icgrid
real*8  :: engyl, engyr
namelist / transcoef /  icgrid, engyl, engyr
!
logical  :: t_read_file
type(mps):: metal_ps
real*8   :: barrier_e,barrier_h,eps_coup,mscoup,relax,eig(2)
integer  :: nst
namelist /plasmon/ barrier_e,barrier_h,eps_coup,relax,nst,eig,mscoup,t_read_file


real*8  :: width,dtmp
integer :: dim07
integer :: i, j, k, l
!
type(dBlkSeq), allocatable :: Lambda(:)

!--------------------------------
!
!--------------------------------

call date_and_time(date,time)
write(6,*)
write(6,1000) date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)


lscba=.false.
lphoton = .false.
tplasmon= .false.
rewind(5)
read(5,nml=general)
!
nLead = 2
norbs = 2
onsite = 0.d0
hop = 1.d0
elcoup = 1.d0
miu0 = 0.d0
half_w =0.1d0
rewind(5)
read(5,nml=system)
rewind(5)

write(6,1003) ' norbs   ', norbs
write(6,1003) ' nLead   ', nLead
write(6,1001) ' elcoup: ', elcoup
write(6,1001) ' hop     ', hop
write(6,1001) ' miu0    ', miu0
write(6,1001) ' onsite  ', onsite
write(6,'(A20,5f12.4)') ' half_w  ', half_w(1:5)

!
telec = 300  
rewind(5)
read(5,nml=temperature)
rewind(5)

write(6,1001) ' temperature: ', telec

beta = ev2j / (boltz * telec)         ! 1/kT in 1/ev  
telec = (telec * boltz)/ev2j          ! Kelvin to ev.

write(6,1001) ' beta(ev^-1): ', beta
write(6,1001) ' telec(ev):   ', telec
!

allocate( miu(nLead), Lambda(nLead) )
do i=1,nLead
  allocate(Lambda(i)%g0(norbs,norbs))
  allocate(Lambda(i)%ge(norbs,norbs))
  Lambda(i)%g0=0.d0
  Lambda(i)%ge=0.d0
enddo
!
ampl = 0.d0
ampr = 0.d0
ampG = 0.d0
rewind(5)
read(5,nml=voltage)

write(6,1001) " left voltage:  ", ampl
write(6,1001) " right voltage: ", ampr
write(6,1001) " gate voltage:  ", ampG

miu(1)=miu0+ampl
miu(2)=miu0+ampr
if(nLead>2) miu(3)=miu0+ampG

!----------------------
! get photon namelist !
!----------------------
mode    = 1
epcoup  = 1.d0
phfreq  = 0.05d0
pt_telec=6.d3
rewind(5)
read(5,photon, end=105)
105 continue

pt_beta = ev2j / (boltz * pt_telec)   ! 1/kT in 1/ev  

write(6,'(A,L2)') ' Is photon considered?  ', lphoton
write(6,'(A,L2)') ' Is plasmon considered? ', tplasmon
write(6,1003) ' number of mode:', mode
write(6,1001) ' photon freq:   ', phfreq
write(6,1001) ' e-p coupling:  ', epcoup
write(6,1001) ' photon telec:  ', pt_telec

if(lphoton.and.tplasmon) then
  write(6,'(A)') "warning: only one of lphoton and tplasmon can be true!"
  stop
endif

if(tPlasmon) then
  barrier_e = 0.0
  barrier_h = 0.0
  eps_coup  = 1.d0
  mscoup    = 1.d0
  relax     = 1.d0
  nst       = 2
  eig(1)     = -1.d0
  eig(2)     = 1.d0
  t_read_file= .false.
  rewind(5)
  read(5,plasmon,end=106)
  106 continue
  metal_ps%barrier_e= barrier_e
  metal_ps%barrier_h= barrier_h
  metal_ps%mscoup   = mscoup
  metal_ps%relax    = relax
  metal_ps%nst      = nst
  write(6,1001) " barrier for electron:      ", metal_ps%barrier_e
  write(6,1001) " barrier for hole:          ", metal_ps%barrier_h
  write(6,1001) " metal-system coupling:     ", metal_ps%mscoup
  write(6,1001) " metal-plasmon coupling:    ", eps_coup
  write(6,1001) " relaxation time:           ", metal_ps%relax
  write(6,1003) " number of states in metal: ", metal_ps%nst
  !
  allocate(metal_ps%eig(nst))
  metal_ps%eig=0.d0
  if(t_read_file) then
    call read_eig(nst,metal_ps%eig,"eig.dat")
  else
    write(6,*) "states of metal"
    do i=1,nst
      if(nst>1) metal_ps%eig(i) = eig(1) + dble(i-1)*(eig(2)-eig(1))/dble(nst-1)
      write(6,'(I5,e15.7)') i,metal_ps%eig(i)
    enddo
  endif
endif
!
ierror = 0
allocate( fock0(norbs,norbs))
allocate(lamdaL(norbs,norbs), lamdaR(norbs,norbs) )
allocate(deltapt(norbs,norbs), deltah(norbs,norbs))
allocate(newhmat(norbs,norbs) )

do i=1,10
  if(ierror(i)/=0) then
    write(6,*)'allocate the', i,'th momoery failed...'
    call flush(6)
    stop
  endif
enddo

fock0=0.d0
fock0(1,1) = -onsite
fock0(2,2) = onsite

!* pin the LUMO level to the left chemical potential
!fock0(1,1) = fock0(1,1) + ampl
!fock0(2,2) = fock0(2,2) + ampr
!fock0(1,2) = hop
!fock0(2,1) = hop 
!*
!-------------------------
!  line-width function   !
!-------------------------
lamdaL=0.d0
lamdaR=0.d0
lamdaL(1,1) = elcoup;
lamdaR(2,2) = elcoup
if(ampr<onsite) then
  dtmp=exp(-(onsite-ampr)/telec)
else
  dtmp=1.d0
endif
lamdaL(2,2) = elcoup*dtmp
lamdaR(1,1) = elcoup*dtmp 

Lambda(1)%g0 = LamdaL
Lambda(2)%g0 = LamdaR

if(nLead>2) then
  Lambda(3)%g0(1,1)=elcoup
  Lambda(3)%g0(2,2)=elcoup
endif

do i=1,nLead
  Lambda(i)%width=half_w(i)
enddo


write(6,*)
write(6,'(A26,2e15.7)') "  fock0(1,1), fock0(2,2) :", fock0(1,1),fock0(2,2)
write(6,'(A26,2e15.7)') "  lamdaL(1,1),lamdaL(2,2):", lamdaL(1,1),lamdaL(2,2)
write(6,'(A26,2e15.7)') "  lamdaR(1,1),lamdaR(2,2):", lamdaR(1,1),lamdaR(2,2)

icgrid = 20
engyl = - 5.d0
engyr = 5.d0
rewind(5)
read(5,transcoef,end=107)
107 continue
!
if(dabs(engyl-engyr)<=1.0d-1) then
  dim07 = icgrid
else
  dim07 = int(dabs(engyl-engyr)/1.0d-1) * icgrid
endif
!
call landauer(nLead,norbs,miu,beta,fock0,Lambda)

if(lphoton) then
  allocate( coup(norbs,norbs,mode))
  allocate( pheigv(mode),bocc(mode) )
  !call epcoupling(norbs,mode,pt_beta,phfreq,epcoup, pheigv, bocc,coup)
  call ep_am15(norbs,mode,pt_beta,phfreq,epcoup,pheigv,bocc,coup,.false.)
endif
!
if(tplasmon) then
  allocate( coup(metal_ps%nst,metal_ps%nst,mode))
  allocate( pheigv(mode),bocc(mode) )
  !if(t_read_file) then
    !------------------------------
    !read dipole matrix from file !
    !------------------------------
  !  call read_coup(nst,mode,pt_beta,phfreq,pheigv,bocc,coup,'dipole-z.dat')
  !else
    call ep_am15(metal_ps%nst,mode,pt_beta,phfreq,eps_coup,pheigv,bocc,coup,.true.)
  !endif
endif

!------------------------------
! get steady state current    !
!------------------------------

if(lphoton.or.tplasmon) then
  if(.not.lscba) then
    !-------
    !  loe !
    !-------
    if(lphoton)  call landauerphlop(nLead,norbs,mode,fock0,Lambda,miu,telec,pheigv,coup)

    !if(tplasmon) call LandauerPsLop(nLead,norbs,mode,fock0,Lambda,metal_ps,miu,telec,pheigv,coup)
    if(tplasmon) call LandauerPsLop2(nLead,norbs,mode,fock0,Lambda,metal_ps,miu,telec,pheigv,coup)
  else
    !-------------------------------
    ! self-consistent calculation  !
    !-------------------------------
    call landauerphscba(nLead,norbs,mode,fock0,Lambda,miu,telec,pheigv,coup)
  endif
endif

call date_and_time(date,time)
write(6,1014) date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)

1000 format(/,' WBL-Heom calculation started on ',/,&
            ' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2,/)
1001 format(A,e15.7)
1002 format(A,f15.7)
1003 format(A,I5)

1014 format(/,' WBL-Heom Calculation Finished on ',/,&
            ' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2)
!---------------
! deallocation !
!---------------

deallocate(deltapt,deltah,newhmat,lamdaL,lamdaR)
deallocate(miu)
do i=1,nLead
  deallocate(Lambda(i)%g0,Lambda(i)%ge)
enddo
deallocate(Lambda)

if(lphoton) then
  deallocate( pheigv, coup,bocc )
endif

!
end program model
