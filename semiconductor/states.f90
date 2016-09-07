program test
implicit none
!
integer :: N,Nm,i,j,k,ni,nk
real*8  :: width
real*8  :: Ec,Ev,dE,Ef,ePt,dos,dos2
real*8, allocatable :: Es(:)
real*8, allocatable :: Em(:)

nk = 4
Ef = 4.0d0
Ec = 1.d0
Ev =-1.d0

!N=(2*nk+1)*2
N=(2*nk+1)*(2*nk+1)*(2*nk+1)*2

Nm=(nint(sqrt(dble(N)))-1)/2


allocate( Es(N),Em( (2*Nm+1)*(2*Nm+1)) )

open(20,file='dos.dat',status='replace')

! get the electronic states
Es=0.d0
Em=0.d0

ni=1
i=0
j=0
do i=-nk,nk
  do j=-nk,nk
    do k=-nk,nk
      Es(ni)   = Ec + (i*i+j*j+k*k)*5.0d-2
      Es(ni+1) = Ev - (i*i+j*j+k*k)*5.0d-2

      write(6,'(I6,2f12.4)') ni, Es(ni),Es(ni+1)
      ni=ni+2
    enddo
  enddo
enddo

!
Em=0.d0
ni=1
do i=-Nm,Nm
  do j=-Nm,Nm
    Em(ni)= (i*i+j*j)*1.5d-2-Ef
   ! write(6,'(I6,f12.4)') ni, Em(ni)
    ni=ni+1
  enddo
enddo

! calculate dos
dE=0.01d0
width=0.10d0

ePt=0.d0
do k=-601,601
  ePt=dble(k-1)*dE
  dos=0.d0
  dos2=0.d0
  do i=1,N
     dos=dos + exp(-(ePt-Es(i))*(ePt-Es(i))/(width*width))
     !dos=dos+width/((ePt-Es(i))*(ePt-Es(i)) + width*width)
  enddo

  do i=1,(2*Nm+1)*(2*Nm+1)
     dos2=dos2+width/((ePt-Em(i))*(ePt-Em(i)) + width*width)
  enddo

  write(20,'(f9.4,2e15.7)') ePt, dos,dos2
enddo

close(20)

deallocate( Es, Em)
end 
