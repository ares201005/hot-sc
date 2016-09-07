subroutine self_ps(ePt,nst,norbs,dimph,pheigv,coup,metal_ps,miu,telec,cselfr,cselfl,cselfg,lwrite)
!
!  first-order BA, the self-energy is constructed by bare electron GF
! 
!  output of greater self-energy is added
!
use transmod, only: dBlkSeq,mps
use transphonon, only: bocc
implicit none
include '../include/parameter'

complex*16,   intent(in) :: ePt
integer,      intent(in) :: norbs,nst,dimph
real*8,       intent(in) :: pheigv(dimph)
complex*16,   intent(in) :: coup(nst,nst,dimph)
type(mps),    intent(in) :: metal_ps
real*8,       intent(in) :: miu, telec
complex*16,   intent(out):: cselfr(norbs,norbs), cselfl(norbs,norbs)
complex*16,   intent(out):: cselfg(norbs,norbs)
logical,      intent(in) :: lwrite
!
real*8,     allocatable :: rwork(:)
complex*16, allocatable :: msMat(:,:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:),greenl(:,:)
complex*16, allocatable :: dgreenr(:,:),dgreeng(:,:),dgreenl(:,:)
!
integer    :: ipiv(nst), info
integer    :: i,ni, nj, nk, kd, istat
real*8     :: weight
real*8     :: fermi !,bocc
complex*16 :: ztmp1,ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(rwork(2*nst*nst))
allocate(ctmp1(nst,nst),ctmp2(nst,nst),ctmp3(nst,nst))
allocate(hmat(nst,nst),greenr(nst,nst),greeng(nst,nst),greenl(nst,nst))
allocate(dgreenr(nst,nst),dgreeng(nst,nst),dgreenl(nst,nst))
!
dgreenr = czero
dgreenl = czero
dgreeng = czero

weight=1.d0/dble(nst)

hmat=czero
do i=1,nst
  if(dabs(dble(ePt)-miu)<=3.d0) then 
    hmat(i,i)=dcmplx(metal_ps%eig(i),-metal_ps%relax)
  else
    hmat(i,i)=dcmplx(metal_ps%eig(i),-metal_ps%relax)
  endif
enddo
!----------------------
!    \Sigma_ps(E)     !
!----------------------
do nk = 1, dimph
do kd = 1, 2
  if(kd==1) then
    ePtq = ePt + dcmplx(pheigv(nk),0.d0)
  else 
    ePtq = ePt - dcmplx(pheigv(nk),0.d0)
  endif

  !fermi=FermiDirac(1.d-3,dble(ePtq),miu)
  fermi=FermiDirac(telec,dble(ePtq),miu)
  !bocc=BoseEinstein(telec,pheigv(nk))

  !greenr=-hmat
  !do ni = 1, nst
  !   greenr(ni,ni) = ePtq + greenr(ni,ni)
  !enddo
  ! 
  !call zgetrf(nst,nst,greenr,nst, ipiv, info)
  !call zgetri(nst,greenr,nst,ipiv,ctmp2,nst*nst,info)! G^r(epsilon+omega)
 
  greenr=czero
  do ni=1,nst
    greenr(ni,ni)=cunity/(ePtq-hmat(ni,ni))
  enddo

  do ni=1,nst
    do nj=1,nst
      greenl(nj,ni)=-2.d0*fermi*eye*dimag(greenr(nj,ni))
      greeng(nj,ni)=-2.d0*(fermi-1.d0)*eye*dimag(greenr(nj,ni))
    enddo
  enddo

  if(kd==1) then
    ctmp1 = bocc(nk)*greenr - 0.5d0*greenl
  else
    ctmp1 = (1.d0+bocc(nk))*greenr + 0.5d0*greenl
  endif
  !
  call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,ctmp1,nst,czero,ctmp2,nst)
  call zgemm('n','n',nst,nst,nst,cunity,ctmp2,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp3,nst)
  !call zlarcm(nst,nst,coup(1:nst,1:nst,nk),nst,ctmp1,nst,ctmp2,nst,rwork)
  !call zlacrm(nst,nst,ctmp2,nst,coup(1:nst,1:nst,nk),nst,ctmp3,nst,rwork)


  !--------------------------------------
  !  retarded and lesser self-energy    !
  !--------------------------------------
  if(kd==1) then
    ztmp1 = dcmplx(0.d0, 0.d0)
    !ztmp1 = dcmplx(1.d0 + bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  endif
 
  call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greenl,nst,czero,ctmp2,nst)
  call zgemm('n','n',nst,nst,nst,cunity,ctmp2,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp1,nst)
  !call zlarcm(nst,nst,coup(1:nst,1:nst,nk),nst,greenl,nst,ctmp2,nst,rwork)
  !call zlacrm(nst,nst,ctmp2,nst,coup(1:nst,1:nst,nk),nst,ctmp1,nst,rwork)
  ctmp1=ctmp1*ztmp1

  dgreenr = dgreenr + ctmp3
  dgreenl = dgreenl + ctmp1

  !-----------------------
  ! greater self-energy  !
  !-----------------------
  if(kd==1) then
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk)+1.d0, 0.d0)
  endif

  call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greeng,nst,czero,ctmp2,nst)
  call zgemm('n','n',nst,nst,nst,cunity,ctmp2,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp1,nst)
  !call zlarcm(nst,nst,coup(1:nst,1:nst,nk),nst,greeng,nst,ctmp2,nst,rwork)
  !call zlacrm(nst,nst,ctmp2,nst,coup(1:nst,1:nst,nk),nst,ctmp1,nst,rwork)
  ctmp1=ctmp1*ztmp1

  dgreeng = dgreeng + ctmp1
enddo !kd
enddo !nk

dgreenr = 0.5d0 * (dgreeng-dgreenl)    

greenr=czero
do ni=1,nst
  greenr(ni,ni)=cunity/(ePt-hmat(ni,ni))
enddo
!greenr=-hmat
!do ni = 1, nst
!  greenr(ni,ni) = ePt + greenr(ni,ni)
!enddo
!
!call zgetrf(nst,nst,greenr,nst,ipiv,info)
!call zgetri(nst,greenr,nst,ipiv,ctmp2,nst*nst, info)

ztmp1=czero
do ni=1,nst
  ztmp1=ztmp1-greenr(ni,ni)
enddo
if(lwrite) write(70,'(f9.3,e15.7)') dble(ePt), dimag(ztmp1)

if(lwrite) write(80,'(f9.4,8e15.7)') dble(ePt),dgreenl(1,1),dgreenl(2,2)
!-----------------------------
! g^r(E)\Sigma_{ps}(E)g^a(E) !
!-----------------------------
call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,dgreenr,nst,czero,ctmp2,nst)
call zgemm('n','c',nst,nst,nst,cunity,ctmp2,nst,greenr,nst,czero,dgreenr,nst)

call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,dgreenl,nst,czero,ctmp2,nst)
call zgemm('n','c',nst,nst,nst,cunity,ctmp2,nst,greenr,nst,czero,dgreenl,nst)

call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,dgreeng,nst,czero,ctmp2,nst)
call zgemm('n','c',nst,nst,nst,cunity,ctmp2,nst,greenr,nst,czero,dgreeng,nst)

if(lwrite) write(90,'(f9.4,8e15.7)') dble(ePt),dgreenl(1,1),dgreenl(2,2)

!---------------------------------------
! V^\dag(E) g^r\Sigma_{ps}(E) g^a(E) V !
!---------------------------------------
allocate(msMat(norbs,nst))
msMat=czero
do nj=1,nst
  do ni=1,norbs

    if(nj<=(nst/2).and.ni>(norbs/2)) cycle
    if(nj>(nst/2).and.ni<=(norbs/2)) cycle
  !  MsMat(ni,nj)=metal_ps%mscoup*cunity

    if(metal_ps%eig(nj)>=metal_ps%barrier_e) then
      msMat(ni,nj)=metal_ps%mscoup*cunity
    else if(metal_ps%eig(nj)>=0.d0.and.metal_ps%eig(nj)<metal_ps%barrier_e) then
      msMat(ni,nj)=metal_ps%mscoup*dexp((metal_ps%eig(nj)-metal_ps%barrier_e)/telec)*cunity
    else if(metal_ps%eig(nj)<0.d0.and.metal_ps%eig(nj)>(-metal_ps%barrier_h)) then
      msMat(ni,nj)=metal_ps%mscoup*dexp((-metal_ps%eig(nj)-metal_ps%barrier_h)/telec)*cunity
    else
      msMat(ni,nj)=metal_ps%mscoup*cunity
    endif
  enddo
enddo

msMat=msMat*weight*dsqrt(weight)

!------------------------
! V^\dag \delta g^< V   !
!------------------------
deallocate(ctmp1)
allocate(ctmp1(norbs,nst))

call zgemm('n','n',norbs,nst,nst,cunity,msMat,norbs,dgreenr,nst,czero,ctmp1,norbs)
call zgemm('n','c',norbs,norbs,nst,cunity,ctmp1,norbs,msMat,norbs,czero,cselfr,norbs)

call zgemm('n','n',norbs,nst,nst,cunity,msMat,norbs,dgreenl,nst,czero,ctmp1,norbs)
call zgemm('n','c',norbs,norbs,nst,cunity,ctmp1,norbs,msMat,norbs,czero,cselfl,norbs)

call zgemm('n','n',norbs,nst,nst,cunity,msMat,norbs,dgreeng,nst,czero,ctmp1,norbs)
call zgemm('n','c',norbs,norbs,nst,cunity,ctmp1,norbs,msMat,norbs,czero,cselfg,norbs)

deallocate(rwork)
deallocate(msMat)
deallocate(greenr, greeng,greenl,hmat) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
deallocate(dgreenr,dgreenl,dgreeng)

end subroutine self_ps
