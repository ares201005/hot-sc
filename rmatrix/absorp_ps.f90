subroutine absorp_ps(ePt,nst,dimph,pheigv,coup,metal_ps,miu,telec,weight,absorb)
!
!  first-order BA, the self-energy is constructed by bare electron GF
! 
!  output of greater self-energy is added
!
use transmod,    only: dBlkSeq,mps
use transphonon, only: bocc
implicit none
include '../include/parameter'

complex*16,   intent(in)  :: ePt
integer,      intent(in)  :: nst,dimph
real*8,       intent(in)  :: pheigv(dimph)
complex*16,   intent(in)  :: coup(nst,nst,dimph)
type(mps),    intent(in)  :: metal_ps
real*8,       intent(in)  :: miu, telec
real*8,       intent(in)  :: weight
real*8,     intent(inout) :: absorb(dimph)
!
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:),greenl(:,:)
complex*16, allocatable :: greeng0(:,:),greenl0(:,:)
!
integer    :: i,ni, nj, nk, kd, istat
real*8     :: fermi,tmp1
complex*16 :: ztmp1,ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(hmat(nst,nst))
allocate(greenr(nst,nst),greeng(nst,nst),greenl(nst,nst)  )
allocate(greeng0(nst,nst),greenl0(nst,nst))
allocate(ctmp1(nst,nst),ctmp2(nst,nst),ctmp3(nst,nst), STAT=istat)

!-------------------------
! get G^<(E) and G^>(E)  !
!-------------------------

hmat=czero
do i=1,nst
  if(dabs(dble(ePt)-miu)<=3.d0) then 
    hmat(i,i)=dcmplx(metal_ps%eig(i),-metal_ps%relax)
  else
    hmat(i,i)=dcmplx(metal_ps%eig(i),-metal_ps%relax)
  endif
enddo

greenr=czero
do ni = 1, nst
   greenr(ni,ni) = cunity/(ePt-hmat(ni,ni))
enddo

fermi=FermiDirac(telec,dble(ePt),miu)

greenl0=czero
greeng0=czero
do ni=1,nst
   greenl0(ni,ni)=-2.d0*fermi*eye*dimag(greenr(ni,ni))
   greeng0(ni,ni)=-2.d0*(fermi-1.d0)*eye*dimag(greenr(ni,ni))
enddo

!--------------------------------------
do nk = 1, dimph
  ePtq = ePt-dcmplx(pheigv(nk),0.d0)

  fermi=FermiDirac(telec,dble(ePtq),miu)

  greenr=czero
  do ni = 1, nst
     greenr(ni,ni) = cunity/(ePtq-hmat(ni,ni))
  enddo
   
  !---------------------------
  !  G^<(E-w)    G^>(E-w)    !
  !---------------------------
  greenl = czero
  greeng = czero
  do ni=1,nst
    greenl(ni,ni)=-2.d0*fermi*eye*dimag(greenr(ni,ni))
    greeng(ni,ni)=-2.d0*(fermi-1.d0)*eye*dimag(greenr(ni,ni))
  enddo

  !--------------------------------------------
  ! M[G^>(E) M G^<(E-w) - G^<(E) M G^>(E-w)]  !
  !--------------------------------------------
  ztmp1 = dcmplx(bocc(nk), 0.d0)

  ctmp3=czero
  do ni=1,nst
    do nj=1,nst
      ctmp3(nj,ni)=greeng0(nj,nj)*coup(nj,ni,nk)*greenl(ni,ni) - &
                   greenl0(nj,nj)*coup(nj,ni,nk)*greeng(ni,ni)
    enddo
  enddo

  call zgemm('n','n',nst,nst,nst,ztmp1,coup(1:nst,1:nst,nk),nst,ctmp3,nst,czero,ctmp2,nst)
  !call dcmulti(nst,nst,coup(1:nst,1:nst,nk),ctmp3,ztmp1,ctmp2)

  do i=1,nst
    absorb(nk)=absorb(nk) + dble(ctmp2(i,i))*weight
  enddo
  
enddo !nk

deallocate(hmat)
deallocate(greenr, greeng,greenl,greenl0,greeng0) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
end subroutine absorp_ps
