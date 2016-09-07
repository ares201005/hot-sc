program main
implicit none

integer :: i,io,nfreq
character(len=32) :: arg
character*25 :: filename
real*8 :: iqe(1000), freq(1000)
real*8 :: flux, curr,pin,pout,nocc
real*8, parameter :: kt=0.5170968d0 ! 6000k
!
do i=1,iargc()
  call getarg(i,arg)
  if(i==1) filename=trim(arg)
enddo

open(20,file=filename,status='old')
read(20,*)

freq= 0.d0
iqe = 0.d0
i=1
do
  read(20,*,iostat=io) freq(i),flux,curr,iqe(i)
  if(io>0) then
    write(6,*) "input error!"
    exit
  else if(io<0) then
    exit
  endif
  i = i+1
enddo
nfreq=i-1


curr=0.d0
pin=0.d0
do i=1,nfreq
  nocc=1.d0/(dexp(freq(i)/kt)-1.d0)
  flux=freq(i)*freq(i)*nocc
  !if(freq(i)>0.8d0.and.freq(i)<2.d0) then 
  !   curr = curr + flux*iqe(i+150)/100.d0
  !else
     curr = curr + flux*iqe(i)/100.d0
  !endif
  pin = pin + flux*freq(i)
  write(60,*) freq(i),nocc, flux
enddo

write(6,*) curr,pin,curr/pin*100.d0

end
