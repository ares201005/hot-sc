subroutine self_ph(ePt,nLead,dim1,dimph,pheigv,coup,fock0,Lambda, &
                    miu,telec, cselfr, cselfl,greenl,cselfg)
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
complex*16,   intent(out)   :: cselfr(dim1,dim1), cselfl(dim1,dim1)
complex*16,   intent(out)   :: greenl(dim1,dim1), cselfg(dim1,dim1)
!
integer,    allocatable :: ipiv(:)
real*8,     allocatable :: tlamda(:,:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:)
!
integer    :: i,ni, nj, nk, kd, istat,info
real*8     :: fermi(nLead),tmp1
complex*16 :: ztmp1,ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(ipiv(dim1))
allocate(hmat(dim1,dim1),tlamda(dim1,dim1))
allocate(ctmp1(dim1,dim1),ctmp2(dim1,dim1),ctmp3(dim1,dim1), STAT=istat)
allocate(greenr(dim1,dim1),greeng(dim1,dim1)  )
!
cselfr = czero
cselfl = czero
cselfg = czero

!
do nk = 1, dimph
do kd = 1, 2
  if(kd==1) then
    ePtq = ePt + dcmplx(pheigv(nk),0.d0)
  else 
    ePtq = ePt - dcmplx(pheigv(nk),0.d0)
  endif

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
!
  call zgemm('n','n',dim1,dim1,dim1,cunity,greenr,dim1,ctmp2,dim1,czero,ctmp1,dim1)
  call zgemm('n','c',dim1,dim1,dim1,cunity,ctmp1,dim1,greenr,dim1,czero,greenl,dim1)  ! lesser green's function
  
  call zgemm('n','n',dim1,dim1,dim1,cunity,greenr,dim1,ctmp3,dim1,czero,ctmp1,dim1)
  call zgemm('n','c',dim1,dim1,dim1,cunity,ctmp1,dim1,greenr,dim1,czero,greeng,dim1)  ! greater green's function
!
  if(kd==1) then
    ctmp1 = bocc(nk)*greenr - 0.5d0*greenl
  else
    ctmp1 = (1.d0+bocc(nk))*greenr + 05d0*greenl
  endif
!
  !call dcmulti(dim1,dim1,coup(1:dim1,1:dim1,nk),ctmp1, cunity, ctmp2)
  !call cdmulti(dim1,dim1,ctmp2,coup(1:dim1,1:dim1,nk), cunity, ctmp3)
  call zgemm('n','n',dim1,dim1,dim1,cunity,coup(1:dim1,1:dim1,nk),dim1,ctmp1,dim1,czero,ctmp2, dim1)
  call zgemm('n','n',dim1,dim1,dim1,cunity,ctmp2,dim1,coup(1:dim1,1:dim1,nk),dim1,czero,ctmp3, dim1)

!  retarded and lesser self-energy
  if(kd==1) then
    ztmp1 = dcmplx(1.d0 + bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  endif
  !call dcmulti(dim1,dim1,coup(1:dim1,1:dim1,nk),greenl,ztmp1,  ctmp2)
  !call cdmulti(dim1,dim1,ctmp2,coup(1:dim1,1:dim1,nk), cunity, ctmp1)

  call zgemm('n','n',dim1,dim1,dim1,ztmp1,coup(1:dim1,1:dim1,nk),dim1,greenl,dim1,czero,ctmp2, dim1)
  call zgemm('n','n',dim1,dim1,dim1,cunity,ctmp2,dim1,coup(1:dim1,1:dim1,nk),dim1,czero,ctmp1, dim1)

  cselfr = cselfr + ctmp3
  cselfl = cselfl + ctmp1

  ! greater self-energy
  if(kd==1) then
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk)+1.d0, 0.d0)
  endif
  !call dcmulti(dim1,dim1,coup(1:dim1,1:dim1,nk),greeng,ztmp1,ctmp2)
  !call cdmulti(dim1,dim1,ctmp2,coup(1:dim1,1:dim1,nk),cunity,ctmp1)

  call zgemm('n','n',dim1,dim1,dim1,ztmp1,coup(1:dim1,1:dim1,nk),dim1,greeng,dim1,czero,ctmp2, dim1)
  call zgemm('n','n',dim1,dim1,dim1,cunity,ctmp2,dim1,coup(1:dim1,1:dim1,nk),dim1,czero,ctmp1, dim1)
  
  !cselfg(kd,kd) = cselfg(kd,kd) + ctmp1(kd,kd)  ! RWA

  cselfg = cselfg + ctmp1
  
enddo !kd
enddo !nk

cselfr = 0.5d0 * (cselfg-cselfl)    

deallocate(hmat,tlamda,ipiv)
deallocate(greenr, greeng) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
end
