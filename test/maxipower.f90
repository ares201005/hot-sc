program main
implicit none
!
integer :: i,io
character(len=32) :: arg
character*25 :: filename
real*8 :: volt, curr(6), power,maxp
do i=1,iargc()
  call getarg(i,arg)
  if(i==1) filename=trim(arg)
  !write(6,'(2A)') "filename=",filename
enddo

open(20,file=filename,status='old')

maxp=0.d0
do
  read(20,*,iostat=io) volt,curr(1:6)
  if(io>0) then
    write(6,*) "input error!"
    exit
  else if(io<0) then
    exit
  endif

  power = 2*volt * (curr(2)/2.d0-curr(4)/2.d0 + curr(1))
  maxp = max(maxp,power)
enddo
write(6,'(A,f12.6)') "max power is: ", maxp

end 
