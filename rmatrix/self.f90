subroutine self1(ePt,norbs,nst,miu,elcoup,eImag,eig,lambda)
implicit none
include '../include/parameter'
!
complex*16, intent(in) :: ePt
integer,    intent(in) :: norbs
integer,    intent(in) :: nst
real*8,     intent(in) :: miu
real*8,     intent(in) :: elcoup
real*8,     intent(in) :: eImag
real*8,     intent(in) :: eig(nst)
real*8,    intent(out) :: lambda(norbs,norbs)
!
integer    :: i,j,k,maxdim
real*8     :: dtmp,weight
complex*16 :: ztmp1,ztmp

!
real*8,     allocatable :: dMat(:,:)
complex*16, allocatable :: cMat(:,:,:)
complex*16, allocatable :: green(:)
!

maxdim=max(norbs,nst)
allocate(dMat(maxdim,maxdim),cMat(maxdim,maxdim,3))
allocate(green(nst))
!
lambda=0.d0

ztmp1=czero
green=czero
weight=1.d0/dble(nst)
do i=1,nst
  ztmp=ePt - (eig(i)+miu) + eye*eImag  ! the chemical potential shifts the energy
  green(i) = cunity/ztmp
  !
  ztmp1=ztmp1-green(i)
enddo
!write(700,'(f9.3,e15.7)') dble(ePt), dimag(ztmp1)

! system-lead coupling matrix
dMat=0.d0
do k=1,nst
  do i=1,norbs
    if(i<=(norbs/2).and.k>(nst/2)) cycle
    if(i>(norbs/2).and.k<=(nst/2)) cycle
    dMat(i,k)=elcoup
  enddo
enddo

dMat=dMat*dsqrt(weight)
!
cMat=czero
do i=1,nst
  cMat(i,i,1)=green(i)
enddo

cMat(1:norbs,1:nst,2)=dMat(1:norbs,1:nst)*cunity

call zgemm('n','n',norbs,nst,nst,cunity,cMat(1,1,2),maxdim,cMat(1,1,1),maxdim,czero,cMat(1,1,3),maxdim)
call zgemm('n','c',norbs,norbs,nst,cunity,cMat(1,1,3),maxdim,cMat(1,1,2),maxdim,czero,cMat(1,1,1),maxdim)

lambda(1:norbs,1:norbs)=-dimag(cMat(1:norbs,1:norbs,1))

dtmp=0.d0
do i=1,norbs
  dtmp=dtmp+lambda(i,i)
enddo
write(700,'(f9.3,3e15.7)') dble(ePt), dtmp,lambda(1,1),lambda(2,2)

deallocate(dMat,cMat,green)

end subroutine self1
