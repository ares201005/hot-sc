subroutine SSHmodel(nat, val, vec)
implicit none

real*8, parameter  :: Mautosi=9.10938291d0, evtosi=1.602176565d0


integer, intent(in) :: Nat
real*8, intent(out) :: val(nat), vec(nat,nat)

integer :: i, j, k, lwork,info
real*8  :: spring0, mass0, dtmp1, dtmp2

real*8,allocatable :: spring(:), mass(:)
real*8,allocatable :: hess(:, :)

lwork = nat * nat


spring0 = 2.d0    ! ev/(A*A)
mass0 = 197.d0    ! a.m.u

Write(6,*) ' Number of sites: ', Nat
write(6,*) ' Spring constant: ', spring0
write(6,*) ' Atomic mass:     ', mass0

allocate( spring(nat), mass(nat), hess(nat, nat) )

! construct the V matrix
spring = spring0
mass = mass0
hess = 0.d0
do i = 1, nat
  if(i==1) then
     hess(i, i) = spring(i)+spring(1)
     hess(i+1, i) = - spring(i)
  else if(i>1.and.i<nat) then
     hess(i,  i) = spring(i) + spring(i-1)
     hess(i+1,i) = - spring(i)
     hess(i-1,i) = - spring(i-1)
  else if(i==nat) then
     hess(i,i) = spring(i) + spring(i-1)
     hess(i-1,i) = -spring(i-1)
  endif
  do j = max(1,i-1), min(i+1,nat)
    hess(j,i) = hess(j,i) / dsqrt(mass(j)*mass(i))/1836.d0
    write(6,'(2I5,f15.7)') j, i, hess(j,i)
  enddo
enddo 
    
!===============

call dsyev('v', 'u', nat, hess, nat, val, vec, lwork, info)
! unit transport
dtmp2 = dsqrt(1.d0/mautosi/evtosi) * 1.054571726d1
write(6,*) ' dtmp2: ', dtmp2

vec = hess
do i = 1, nat
  val(i) = dsqrt(val(i))*dtmp2
  write(6,*) ' ith Eigenvalue: ', i, val(i)
  dtmp1 = 0.d0
  write(6,*) ' displacement for ith mode! '
  do j = 1, nat
    write(6,*) j, vec(j,i)
    dtmp1 = dtmp1 + vec(j,i) * vec(j,i)
  enddo
  write(6,*) ' vec(*,i) * vec(*,i): ', dtmp1

enddo

deallocate( spring, mass, hess )


end subroutine SSHmodel
