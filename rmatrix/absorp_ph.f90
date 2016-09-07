subroutine absorp_ph(ePt,nLead,dim1,dimph,pheigv,coup,fock0,Lambda,miu,telec,weight,absorb)
!
!  first-order BA, the self-energy is constructed by bare electron GF
! 
!  output of greater self-energy is added
!
use transmod,    only: dBlkSeq
use transphonon, only: bocc
implicit none
include '../include/parameter'

complex*16,   intent(in)    :: ePt
integer,      intent(in)    :: nLead,dim1,dimph
real*8,       intent(in)    :: pheigv(dimph)
complex*16,   intent(in)    :: coup(dim1,dim1,dimph)
real*8,       intent(in)    :: fock0(dim1,dim1)
type(dBlkSeq),intent(inout) :: Lambda(nLead)
real*8,       intent(in)    :: miu(nLead), telec
real*8,       intent(in)    :: weight
real*8,     intent(inout)   :: absorb(dimph)
!
real*8,     allocatable :: tlamda(:,:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:),greenl(:,:)
complex*16, allocatable :: greeng0(:,:),greenl0(:,:)
!
integer    :: ipiv(dim1), info
integer    :: i,ni, nj, nk, kd, istat
real*8     :: fermi(nLead),tmp1
complex*16 :: ztmp1,ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(hmat(dim1,dim1),tlamda(dim1,dim1))
allocate(ctmp1(dim1,dim1),ctmp2(dim1,dim1),ctmp3(dim1,dim1), STAT=istat)
allocate(greenr(dim1,dim1),greeng(dim1,dim1),greenl(dim1,dim1)  )
allocate(greeng0(dim1,dim1),greenl0(dim1,dim1)  )
!

!-------------------------
! get G^<(E) and G^>(E)  !
!-------------------------

call self_lead(ePt,nLead,dim1,miu,fock0,Lambda,tlamda)
hmat=dcmplx(fock0,-tlamda)

greenr=-hmat
do ni = 1, dim1
   greenr(ni,ni) = ePt + greenr(ni,ni)
enddo

call zgetrf(dim1,dim1, greenr, dim1, ipiv, info)
call zgetri(dim1, greenr, dim1, ipiv, ctmp2,dim1*dim1, info)! G^r(epsilon+omega)
!
ctmp2 = czero
ctmp3 = czero
do i=1,nLead
  fermi(i)=FermiDirac(telec,dble(ePt),miu(i))
  ctmp2 = ctmp2 + 2.d0*eye*fermi(i)*Lambda(i)%ge
  ctmp3 = ctmp3 + 2.d0*eye*(fermi(i)-1.d0)*Lambda(i)%ge
enddo

! G^<(E)
call zgemm('n','n',dim1,dim1,dim1,cunity,greenr,dim1,ctmp2,dim1,czero,ctmp1,dim1)
call zgemm('n','c',dim1,dim1,dim1,cunity,ctmp1,dim1,greenr,dim1,czero,greenl0,dim1)

! G^>(E) 
call zgemm('n','n',dim1,dim1,dim1,cunity,greenr,dim1,ctmp3,dim1,czero,ctmp1,dim1)
call zgemm('n','c',dim1,dim1,dim1,cunity,ctmp1,dim1,greenr,dim1,czero,greeng0,dim1)

!--------------------------------------
do nk = 1, dimph
  ePtq = ePt - dcmplx(pheigv(nk),0.d0)

  do i=1,nLead
     fermi(i)=FermiDirac(telec,dble(ePtq),miu(i))
  enddo

  !------------------------
  !energy-dependent Gamma !
  !------------------------

  call self_lead(ePtq,nLead,dim1,miu,fock0,Lambda,tlamda)
  hmat=dcmplx(fock0,-tlamda)

  greenr=-hmat
  do ni = 1, dim1
     greenr(ni,ni) = ePtq + greenr(ni,ni)
  enddo
   
  call zgetrf(dim1,dim1, greenr, dim1, ipiv, info)
  call zgetri(dim1, greenr, dim1, ipiv, ctmp2,dim1*dim1, info)! G^r(epsilon+omega)
!
  ctmp2 = czero
  ctmp3 = czero
  do i=1,nLead
    ! sigma^< = i f(epsilon-miu) * Lamda  ! lesser self-energy
    ! greater self-energy
    ctmp2 = ctmp2 + 2.d0*eye*fermi(i)*Lambda(i)%ge
    ctmp3 = ctmp3 + 2.d0*eye*(fermi(i)-1.d0)*Lambda(i)%ge
  enddo

  !---------------
  !  G^<(E-w)    !
  !---------------
  call zgemm('n','n',dim1,dim1,dim1,cunity,greenr,dim1,ctmp2,dim1,czero,ctmp1,dim1)
  call zgemm('n','c',dim1,dim1,dim1,cunity,ctmp1,dim1,greenr,dim1,czero,greenl,dim1)
  
  !---------------
  !  G^>(E-w)    !
  !---------------
  call zgemm('n','n',dim1,dim1,dim1,cunity,greenr,dim1,ctmp3,dim1,czero,ctmp1,dim1)
  call zgemm('n','c',dim1,dim1,dim1,cunity,ctmp1,dim1,greenr,dim1,czero,greeng,dim1)

  !--------------------------------------------
  ! M[G^>(E) M G^<(E-w) - G^<(E) M G^>(E-w)]  !
  !--------------------------------------------
  ztmp1 = dcmplx(bocc(nk), 0.d0)
  !call dcmulti(dim1,dim1,coup(1:dim1,1:dim1,nk),greenl,cunity,ctmp1)
  !call dcmulti(dim1,dim1,coup(1:dim1,1:dim1,nk),greeng,cunity,ctmp2)
  call zgemm('n','n',dim1,dim1,dim1,cunity,coup(1:dim1,1:dim1,nk),dim1,greenl,dim1,czero,ctmp1,dim1)
  call zgemm('n','n',dim1,dim1,dim1,cunity,coup(1:dim1,1:dim1,nk),dim1,greeng,dim1,czero,ctmp2,dim1)

  call zgemm('n','n',dim1,dim1,dim1,cunity,greeng0,dim1,ctmp1,dim1,czero,ctmp3,dim1)
  call zgemm('n','n',dim1,dim1,dim1,-1.d0*cunity,greenl0,dim1,ctmp2,dim1,cunity,ctmp3,dim1)

  !call dcmulti(dim1,dim1,coup(1:dim1,1:dim1,nk),ctmp3,ztmp1,ctmp2)
  call zgemm('n','n',dim1,dim1,dim1,ztmp1,coup(1:dim1,1:dim1,nk),dim1,ctmp3,dim1,czero,ctmp2,dim1)

  do i=1,dim1
    absorb(nk)=absorb(nk) + dble(ctmp2(i,i))*weight
  enddo
  
enddo !nk


deallocate(hmat,tlamda)
deallocate(greenr, greeng,greenl,greenl0,greeng0) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
end subroutine absorp_ph
