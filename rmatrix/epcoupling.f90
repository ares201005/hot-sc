subroutine epcoupling(norbs, mode, beta,freq, epcoup, pheigv, bocc,coup)
implicit none

include '../include/parameter'
!
! electron-photon coupling parameters
!
integer, intent(in) :: mode,norbs
real*8,  intent(in) :: beta,freq, epcoup
real*8, intent(out) :: pheigv(mode), bocc(mode), coup(norbs,norbs,mode)
!
integer :: i, j, k
!

coup = 0.d0

write(6,*)
write(6,*) '   occupation number of photon '
do k = 1, mode
  pheigv(k) = freq
  !coup(i,i,i) = epcoup
  do i = 1, norbs
    do j = i+1, norbs
      coup(j, i, k) = epcoup
      coup(i, j, k) = epcoup
      write(6,'(A20,2I5,2e15.7)') ' coup(i,j,k)',j,i, coup(i,j,k),coup(j,i,k)
    enddo
  enddo
  !bocc(k) = 1.d0 !/(dexp( beta * pheigv(k)) - 1.d0)
  bocc(k) = 1.d0 /(dexp( beta * pheigv(k)) - 1.d0)
  write(6,'(A15,I4,3e15.7)') 'the ith mode:', k, pheigv(k), bocc(k), 1.d0 /(dexp( beta * pheigv(k)) - 1.d0)
enddo
!
end subroutine epcoupling

!-------------------------------------
!
!-------------------------------------
subroutine ep_am15(norbs,mode,beta,freq,epcoup,pheigv,bocc,coup,tdielectric)
implicit none
include '../include/parameter'
integer,   intent(in)  :: norbs, mode
real*8,    intent(in)  :: beta, freq, epcoup
real*8,    intent(out) :: pheigv(mode),bocc(mode)
complex*16,intent(out) :: coup(norbs,norbs,mode)
logical,   intent(in)  :: tdielectric
!
integer :: i,j,k
real*8  :: dE,tmp,flux,factor,power
complex*16 :: epsilonw,ztmp

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
if(tdielectric) then
  write(6,*) ' node   energy        occ           flux      coupling-factor , (eps-1/(eps+2)-1'
else
  write(6,*) ' node   energy        occ           flux      coupling-factor '
endif
power=0.d0
do k=1,mode
  flux=pheigv(k)*pheigv(k)/(dexp(pheigv(k)*beta)-1.d0)/factor
  power=power+flux*pheigv(k)*dE

  bocc(k) = 1.d0/(dexp( beta * pheigv(k)) - 1.d0)
  !bocc(k) = 1.d0
  !dE is the weight of summation over photon mode
  !tmp=dsqrt(flux/bocc(k)/pheigv(k)*dE)  
  tmp=dsqrt(flux/bocc(k)*pheigv(k)*dE)  

  do i=1,norbs/2
    do j=norbs/2+1,norbs
    !do j=i+1,norbs
      coup(j,i,k)=dcmplx(epcoup*tmp,0.d0)
      coup(i,j,k)=dcmplx(epcoup*tmp,0.d0)
    enddo
  enddo

  if(tdielectric) then
    call dielectric(3,'Au',pheigv(k),epsilonw)
    ztmp = (epsilonw-1.d0)/(epsilonw+2.d0) - 1.d0
    coup(1:norbs,1:norbs,k) = coup(1:norbs,1:norbs,k) * ztmp
    !coup(1:norbs,1:norbs,k) = coup(1:norbs,1:norbs,k) * dcmplx(0.d0,dimag(ztmp))
  endif

  if(tdielectric) then
    write(6,'(I4,2X, f9.4, 5e15.5)') k, pheigv(k), bocc(k), flux, tmp, ztmp
  else
    write(6,'(I4,2X, f9.4, 3e15.5)') k, pheigv(k), bocc(k), flux, tmp
  endif
  !
enddo

write(6,'(A,f9.4,A)') "  power density is:", power,'kW/m^2/ev'

end subroutine ep_am15
