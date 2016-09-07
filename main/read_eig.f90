subroutine read_eig(nst,eig,fname)
implicit none
integer,         intent(in) :: nst
real*8,          intent(out):: eig(nst)
character(len=*),intent(in) :: fname
!
integer :: i,j

open(333,file=fname,status='old')
do i=1,nst
  read(333,*) j, eig(i)
  write(6,*) j, eig(i)
enddo
close(333)

end subroutine read_eig


subroutine read_coup(nst,mode,beta,freq,pheigv,bocc,coup,fname)
implicit none
integer,         intent(in)  :: nst
integer,         intent(in)  :: mode
real*8,          intent(in)  :: beta, freq
real*8,          intent(out) :: pheigv(mode),bocc(mode)
complex*16,      intent(out) :: coup(nst,nst,mode)
character(len=*),intent(in)  :: fname
!
complex*16 :: epsilonw,ztmp
real*8     :: power,dE,tmp,flux,factor
integer    :: i,k
!
real*8, allocatable :: dipole(:,:)
!

! read in dipole matrix
write(6,*)
write(6,*) "read dipole matrix from the file"

allocate(dipole(nst,nst))

coup=dcmplx(0.d0,0.d0)
open(444,file=fname,status='old')
do i=1,nst
  read(444,*) dipole(1:nst,i)
enddo
close(444)
write(6,*) 'readfile finished!'
!
pheigv= 0.d0
bocc  = 0.d0
coup  = cmplx(0.d0,0.d0)

dE=freq/dble(mode)

if(mode==1) dE=1.0d0

factor=0.d0
do k=1,mode
  pheigv(k)= dble(k-0.5d0)*dE
  if(mode==1) pheigv(k)=freq
  flux=pheigv(k)*pheigv(k)/(dexp(pheigv(k)*beta)-1.d0)
  factor=factor+flux*pheigv(k)*dE
enddo

write(6,*)
write(6,*) ' node   energy        occ           flux      coupling-factor , (eps-1/(eps+2)-1'
power=0.d0
do k=1,mode
  flux=pheigv(k)*pheigv(k)/(dexp(pheigv(k)*beta)-1.d0)/factor
  power=power+flux*pheigv(k)*dE

  call dielectric(2,'Ag',pheigv(k),epsilonw)
  ztmp = (epsilonw-1.d0)/(epsilonw+2.d0) - 1.d0

  bocc(k) = 1.d0/(dexp( beta * pheigv(k)) - 1.d0)
  !dE is the weight of summation over photon mode
  tmp=dsqrt(flux/bocc(k)/pheigv(k)*dE)  

  coup(1:nst,1:nst,k)=dcmplx(dipole(1:nst,1:nst),0.d0)*tmp * ztmp
  write(6,'(I4,2X, f9.4, 5e15.5)') k, pheigv(k), bocc(k), flux, tmp, ztmp
  !
enddo

write(6,'(A,f9.4,A)') "  power density is:", power,'kW/m^2/ev'

deallocate(dipole)
end subroutine read_coup
