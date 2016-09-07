subroutine LandauerPsLop2(nLead,norbs,nmode,fock0,Lambda,metal_ps,miu,telec,pheigv,coup)
use transmod, only: dBlkSeq,mps
implicit none
include '../include/parameter'
include '../include/trans_data'

!------------------------------------------------------------------------------
! Calculate steady state current within WBL and electron-photon coupling      !
! by integration of transmission function. There are two parts of             !
! transmission function, one is the elastic part:                             !
!  T1(E) = 4 Tr( G^r * Lambda_R * G^a * Lambda_L )                            !
!                                                                             !
! another is the inelastic part:                                              !
!  T2(E) = \Sigma_L G^r(E)\Sigma_{ep}(E)G^a(E) -                              !
!                                                                             !
! Current is  int [(f_L-f_R)T1(E) + T2(E)]                                    !
!                                                                             !
! note: there is not fermi function before T2 part,                           !
!       so integration over total energy range is needed                      !
!                                                                             !
!  Lowest order expansion method is used                                      !
! the effect of system on the metal is consered                               !
!------------------------------------------------------------------------------

integer,       intent(in)    :: nLead,norbs,nmode 
type(dBlkSeq), intent(inout) :: Lambda(nLead)
real*8,        intent(in)    :: fock0(norbs,norbs)
type(mps),     intent(in)    :: metal_ps
real*8,        intent(in)    :: miu(nLead),telec
real*8,        intent(in)    :: pheigv(nmode)
complex*16,    intent(in)    :: coup(metal_ps%nst,metal_ps%nst,nmode)
!
integer :: nsplit,npoint
integer :: i,j,k,ni,nj,nk,iLead,jLead,istat,info,ipiv(norbs)
real*8  :: etop, ebot
real*8  :: eImag,dengy, tmp2,tmp1,telec0
real*8  :: fermi(nLead), volt
real*8  :: tcoef(nLead)
!
integer :: ngrid
real*8  :: downrange=0.d0, uprange=0.d0
namelist /land/ ngrid, downrange, uprange,eImag
!
real*8                  :: flux(2),flux0(2)
real*8,     allocatable :: hot(:,:),hot0(:,:)
real*8,     allocatable :: angPt(:), weight(:)
real*8,     allocatable :: tlamda(:,:)
real*8,     allocatable :: absorb(:), absorb0(:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: green(:,:),dgreen(:,:)
complex*16, allocatable :: ctmp1(:, :), ctmp2(:, :)
complex*16, allocatable :: ctmp3(:, :), ctmp4(:, :)
complex*16, allocatable :: cselfr(:,:), cselfl(:,:),cselfg(:,:)
!
character*10   :: pFile
character*15   :: tFile
complex*16     :: ePt
real*8         :: dtmp
real*8         :: curr(nLead), curr2(nLead)
real*8         :: maxcurr
real*8,external:: FermiDirac
!
write(6,*) 
write(6,*) ' =========== entering LandauerPsLop2 ================'
write(6,*) ' calculate steady state current with e-p coupling'

!-----------------------
! read in parameters   !
!-----------------------
ngrid=10
eImag=1.d-5
rewind(5)
read(5,nml=land)
call flush(6)

ebot = minval(miu) - maxval(pheigv) - downrange  !initial point of energy
etop = maxval(miu) + maxval(pheigv) + uprange    !end point of energy

npoint=mgauleg1
nsplit=int((etop-ebot)/0.1d0*ngrid)
allocate(angPt(nsplit*npoint),weight(nsplit*npoint))
  
write(6,'(A,2f9.3)') ' integration range:', ebot, etop
write(6,'(A,I7)')    ' number of grid:   ', nsplit*npoint

k=0
do i=1,nsplit
  tmp1=(etop-ebot)/dble(nsplit)*dble(i-1) + ebot
  tmp2=(etop-ebot)/dble(nsplit)*dble(i)   + ebot
  do j=1,npoint
    k=k+1
    angPt(k)=((tmp2-tmp1)*dpt(j) + (tmp1+tmp2))*0.5d0
    weight(k)=dw(j)*(etop-ebot)*0.5d0/dble(nsplit)
  enddo
enddo
!
telec0 = telec
!
allocate(hot(metal_ps%nst,2),hot0(metal_ps%nst,2))
allocate(absorb(nmode),absorb0(nmode))
allocate(tlamda(norbs,norbs))
allocate(hmat(norbs,norbs), green(norbs,norbs), dgreen(norbs,norbs),STAT=istat)
allocate(ctmp1(norbs, norbs), ctmp2(norbs, norbs), stat=istat)
allocate(ctmp3(norbs, norbs), ctmp4(norbs, norbs), STAT=istat)
allocate(cselfr(norbs,norbs),cselfl(norbs,norbs),cselfg(norbs,norbs) )
!
pFile  = 'DOS-ps.dat'
curr   = 0.d0
curr2  = 0.d0
absorb = 0.d0
flux   = 0.d0
hot    = 0.d0
do iLead=1,nLead
  write(tFile,'(A9,I1,A4)') "Tcoef-ps-",iLead,".dat"
  open(89, file=tFile, form='formatted', status='replace')
  ! 
  if(iLead==1) open(99, file=pFile, form='formatted', status='replace')
  if(iLead==1) open(79, file="absorb.dat", form='formatted', status='replace')

  do nk=1,nsplit*npoint
    ePt = dcmplx(angPt(nk),eImag)
    dengy=weight(nk)
 
    !--------------------------------------------------
    ! calculate plasmon-induced hot-carrier injection !
    ! in terms of self-energies                       !
    !--------------------------------------------------
    if(iLead==1) then 
      call self_ps_new(ePt,metal_ps%nst,norbs,nLead,nmode,pheigv,coup,metal_ps,miu,telec0,&
                       fock0,Lambda,cselfr,cselfl,cselfg,.true.,.true.,.true.,absorb0,flux0,hot0)
      absorb=absorb+absorb0*dengy
      flux=flux+flux0*dengy
      hot = hot+hot0*dengy
      write(79,'(f9.4,e15.7)') dble(ePt), sum(absorb0)
    else
      call self_ps_new(ePt,metal_ps%nst,norbs,nLead,nmode,pheigv,coup,metal_ps,miu,telec0, &
                       fock0,Lambda,cselfr,cselfl,cselfg,.false.,.false.,.false.,absorb0,flux0,hot0)
    endif

    cSelfr=czero

    !------------------------
    !energy-dependent Gamma !
    !------------------------
    call self_lead(ePt,nLead,norbs,miu,fock0,Lambda,tlamda)

    hmat=dcmplx(fock0,-tlamda)

    green = -hmat
    do ni=1,norbs
       green(ni,ni) = ePt + green(ni,ni)
    enddo
    !
    call zgetrf(norbs, norbs, green, norbs, ipiv, info)
    call zgetri(norbs, green, norbs, ipiv, ctmp2, norbs*norbs, info) 

    !--------
    ! DOS   !
    !--------
    dtmp=0.d0
    tmp1=0.d0
    do ni=1,norbs
      dtmp=dtmp-dimag(green(ni,ni))
      tmp1=tmp1+tLamda(ni,ni)
    enddo
    if(iLead==1) write(99,'(f9.4,6e15.7)') dble(ePt), dtmp, tmp1,tLamda(1,1),tLamda(2,2),tLamda(1,2),tLamda(2,1)

    fermi(iLead)=FermiDirac(telec,dble(ePt),miu(iLead))

    !----------------------------------
    ! delta G = G^r_0 \sigma^r G^r_0  !
    !----------------------------------
    call zgemm('n','n',norbs,norbs,norbs,cunity,green,norbs,cselfr,norbs,czero,ctmp2,norbs)
    call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,green, norbs,czero,dgreen,norbs)

    !----------------------------------------
    ! for the T1 part, right to left        !
    ! LambdaR -> ctmp2 (lower half)         !
    !----------------------------------------

    do jLead=1,nLead
      if(iLead==jLead) cycle
      fermi(jLead)=FermiDirac(telec,dble(ePt),miu(jLead))
      ctmp2=dcmplx(Lambda(jLead)%ge,0.d0)
      !---------------------------------------------------- 
      ! (G^r + delta G)* LambdaR (lower, right)  -> ctmp3 !
      !---------------------------------------------------- 
      call zgemm('n','n',norbs,norbs,norbs,cunity,green,norbs,ctmp2,norbs,czero,ctmp3,norbs)
      call zgemm('n','n',norbs,norbs,norbs,2.d0*cunity,dgreen,norbs,ctmp2,norbs,cunity,ctmp3,norbs)
      !---------------------------------------------------- 
      ! ctmp3 * G^a = ctmp3 * (G^r)^H -> ctmp2            !
      !---------------------------------------------------- 
      call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp3,norbs,green,norbs,czero,ctmp2,norbs)
      !
      ctmp1 = dcmplx(Lambda(iLead)%ge,0.d0)

      !---------------------------------------------------- 
      ! ctmp2 * LambdaL (lower, right) -> ctmp3           !
      !---------------------------------------------------- 
      call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,ctmp1,norbs,czero,ctmp3,norbs)

      tmp1 = 0.d0
      do ni=1,norbs
        tmp1 = tmp1 + dble(ctmp3(ni,ni))/(pi * pi)
      enddo
      if(jLead>iLead) then 
        tCoef(jLead-1)=tmp1*2.d0*pi
      else
        tCoef(jLead)=tmp1*2.d0*pi
      endif

      tmp1 = tmp1 * (fermi(iLead) - fermi(jLead))
      dtmp = 2.d0 * pi * tmp1 * dengy
      curr(iLead) = curr(iLead) + dtmp
    enddo

    !------------------------------------------
    ! for T2 part, to left                    !
    !   G^r * Sigma^>_ep * G^a * Sigma^<_L    !
    ! - G^r * Sigma^<_ep * G^a * Sigma^>_L    !
    !------------------------------------------
    ctmp3 = fermi(iLead) * cselfg - (fermi(iLead)-1.d0) * cselfl

    call zgemm('n','n',norbs,norbs,norbs,cunity,green,norbs,ctmp3,norbs,czero,ctmp1,norbs)
    call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,green,norbs,czero,ctmp2,norbs)

    ctmp1=eye*Lambda(iLead)%ge

    call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,ctmp1,norbs,czero,ctmp3,norbs)

    tmp1 = 0.d0 
    do ni = 1, norbs
      tmp1 = tmp1 + dble(ctmp3(ni,ni))/(pi * pi) 
    enddo
    tcoef(nLead)=tmp1*pi

    curr2(iLead) = curr2(iLead) + 1.d0 * pi * tmp1 * dengy

    write(89,'(5e15.7)') dble(ePt), tCoef(1:nLead)
  enddo
  close(89)
  if(iLead==1) close(99)
  if(iLead==1) close(79)
enddo

curr =curr *1.6022d5*hbarinv
curr2=curr2*1.6022d5*hbarinv

do iLead=1,nLead
  write(6,1000) '  current  through lead: ', iLead, " is ", curr(iLead)
  write(6,1000) '  current2 through lead: ', iLead, " is ", curr2(iLead)
  write(6,*)
enddo
!
hot=hot*hbarinv/(2.d0*pi)
open(555,file='hot.dat',status='replace')
do ni=1,metal_ps%nst
  write(555,'(I5,f12.4,2e15.4)') ni,metal_ps%eig(ni), hot(ni,1), hot(ni,2)
enddo
write(555,'(2e15.4)') sum(hot(:,1)), sum(hot(:,2))
close(555)

absorb=absorb*hbarinv/(2.d0*pi)
flux=flux*hbarinv/(2.d0*pi)
maxcurr= sum(absorb)*1.6022d5
write(6,*)
write(6,'(A,2e15.6,A)') '  hot-carrier injection flux is ', flux(1:2), ' /fs'
write(6,'(A,2e15.6,A)') '  hot-carrier injection curr is ', flux(1:2)*1.6022d5, ' nA'
write(6,*)
write(6,'(A,e15.7,A)') '  absorped photon flux is ', sum(absorb), ' /fs'
write(6,'(A,e15.7,A)') '  maximum photocurrent is ', maxcurr, ' nA'
write(6,'(A,f10.3,A)') '  IQE is                  ', abs(curr2(1)-curr2(2))/2.d0/maxcurr*1d2,'%'

write(6,*) ' ========== leaving LandauerPsLop ============='

1000 format(A, I2, A, e15.7, " nA")

!
deallocate(hot)
deallocate(absorb0,absorb)
deallocate(tlamda)
deallocate(ctmp1, ctmp2, ctmp3, ctmp4, STAT=istat)
deallocate(hmat, green, dgreen, STAT=istat)
deallocate(cselfr, cselfl,cselfg )
!
end subroutine LandauerPsLop2


!----------------------
! self-energy of lead !
!----------------------
subroutine self_lead(ePt,nLead,norbs,miu,fock0,Lambda,tlamda)
use transmod, only: dBlkSeq,mps
implicit none
complex*16,    intent(in)    :: ePt
integer,       intent(in)    :: nLead
integer,       intent(in)    :: norbs
real*8,        intent(in)    :: miu(nLead)
real*8,        intent(in)    :: fock0(norbs,norbs)
type(dBlkSeq), intent(inout) :: Lambda(nLead)
real*8,        intent(out)   :: tlamda(norbs,norbs)
!
integer :: i,j
real*8  :: tmp1,tmp2

tlamda=0.d0

do i=1,nLead
  tmp1=Lambda(i)%width
  tmp1=tmp1*tmp1/((dble(ePt)-miu(i))*(dble(ePt)-miu(i))+tmp1*tmp1)

  tmp1=0.d0
  tmp2=0.d0
  if(dble(ePt)>=fock0(norbs/2+1,norbs/2+1)) then
  !if(dble(ePt)>=fock0(norbs/2+1,norbs/2+1).and.i==2) then
    tmp2= sqrt(dble(ePt)-fock0(norbs/2+1,norbs/2+1))
  else if(dble(ePt)<=fock0(norbs/2,norbs/2)) then
  !else if(dble(ePt)<=fock0(norbs/2,norbs/2).and.i==1) then
    tmp1= sqrt(-dble(ePt)+fock0(norbs/2,norbs/2))
  endif

  Lambda(i)%ge = 0.d0
  Lambda(i)%ge(1:norbs/2,1:norbs/2)=Lambda(i)%g0(1:norbs/2,1:norbs/2)*tmp1
  Lambda(i)%ge(norbs/2+1:norbs,norbs/2+1:norbs)=Lambda(i)%g0(norbs/2+1:norbs,norbs/2+1:norbs)*tmp2
  tlamda=tlamda+Lambda(i)%ge
enddo

end subroutine self_lead

!---------------------------------
! self-energy (system to metal)  !
!---------------------------------
subroutine self_s(ePt,nLead,norbs,nst,miu,telec,fock0,Lambda,tlamda,msMat,selfs,selfl,selfg,tless)
use transmod, only: dBlkSeq
implicit none
include '../include/parameter'
!
complex*16,  intent(in) :: ePt
integer,     intent(in) :: nLead
integer,     intent(in) :: norbs
integer,     intent(in) :: nst
real*8,      intent(in) :: miu(nLead)
real*8,      intent(in) :: telec
real*8,      intent(in) :: fock0(norbs,norbs)
type(dBlkSeq),intent(in):: Lambda(nLead)
real*8,      intent(in) :: tlamda(norbs,norbs)
complex*16,  intent(in) :: msMat(norbs,nst)
complex*16, intent(out) :: selfs(nst,nst)
complex*16, intent(out) :: selfl(nst,nst)
complex*16, intent(out) :: selfg(nst,nst)
logical,    intent(in)  :: tless
!
integer :: i,ni,info
integer,    allocatable :: ipiv(:)
complex*16, allocatable :: grs(:,:)
complex*16, allocatable :: gls(:,:)
complex*16, allocatable :: ctmp1(:,:)
!
real*8                  :: fermi
real*8,     external    :: FermiDirac
!
allocate(ipiv(norbs))
allocate(grs(norbs,norbs),gls(norbs,norbs))
allocate(ctmp1(nst,nst))

grs=-dcmplx(fock0,-tlamda)  
do ni=1,norbs
   grs(ni,ni) = ePt + grs(ni,ni)
enddo

call zgetrf(norbs,norbs,grs,norbs, ipiv, info)
call zgetri(norbs,grs,norbs,ipiv,ctmp1,nst*nst,info)! G^r(epsilon+omega)

if(tless) then
  ! g^<_S
  gls = dcmplx(0.d0,0.d0)
  do i=1,nLead
    !fermi=FermiDirac(telec,dble(ePt),(miu(1)+miu(2))/2.d0)
    fermi=FermiDirac(telec,dble(ePt),miu(i))
    gls= gls+2.d0*eye*fermi*Lambda(i)%ge
  enddo
  call zgemm('n','n',norbs,norbs,norbs,cunity,grs,norbs,gls,norbs,czero,ctmp1,nst)
  call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,nst,grs,norbs,czero,gls,norbs)

  call zgemm('n','n',norbs,nst,norbs,cunity,gls,norbs,msMat,norbs,czero,ctmp1,nst)
  call zgemm('c','n',nst,nst,norbs,cunity,msMat,norbs,ctmp1,nst,czero,selfl,nst)

  ! g^>_S
  gls = dcmplx(0.d0,0.d0)
  do i=1,nLead
    fermi=FermiDirac(telec,dble(ePt),miu(i))
    !fermi=FermiDirac(telec,dble(ePt),(miu(1)+miu(2))/2.d0)
    gls= gls+2.d0*eye*(fermi-1.d0)*Lambda(i)%ge
  enddo
  call zgemm('n','n',norbs,norbs,norbs,cunity,grs,norbs,gls,norbs,czero,ctmp1,nst)
  call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,nst,grs,norbs,czero,gls,norbs)

  call zgemm('n','n',norbs,nst,norbs,cunity,gls,norbs,msMat,norbs,czero,ctmp1,nst)
  call zgemm('c','n',nst,nst,norbs,cunity,msMat,norbs,ctmp1,nst,czero,selfg,nst)
endif

! V_{mi}G_i V_{im} 
call zgemm('n','n',norbs,nst,norbs,cunity,grs,norbs,msMat,norbs,czero,ctmp1,nst)
call zgemm('c','n',nst,nst,norbs,cunity,msMat,norbs,ctmp1,nst,czero,selfs,nst)

! get the imaginary part
!selfs =   dimag(selfg-selfl)*eye/2.d0

!selfs =   dimag(selfs)*eye

!selfs=-selfs

do ni=1,nst
  do i=1,nst
    if((-dimag(selfs(i,ni)))<0.d0) then
      write(6,*) 'warning:, negative line-width function in self_s!'
    endif
  enddo
enddo

!do ni=1,nst
!  do i=1,nst
!    if(i==ni) cycle
!    selfs(i,ni)= cmplx(0.d0,0.d0)
!    selfg(i,ni)= cmplx(0.d0,0.d0)
!    selfl(i,ni)= cmplx(0.d0,0.d0)
!  enddo
!enddo

deallocate(grs,gls,ctmp1,ipiv)
!
end subroutine self_s



!---------------------------------------------------
! plasmon induced self-energy, effect of system on !
! the metal is taken into account                  !
!---------------------------------------------------
subroutine self_ps_new(ePt,nst,norbs,nLead,dimph,pheigv,coup,metal_ps,miu,telec,&
                       fock0,Lambda,cselfr,cselfl,cselfg,lwrite,labsorb,lhot,absorb,flux,hot)
!
use transmod, only: dBlkSeq,mps
use transphonon, only: bocc
implicit none
include '../include/parameter'

complex*16,     intent(in) :: ePt
integer,        intent(in) :: norbs,nst,nLead,dimph
real*8,         intent(in) :: pheigv(dimph)
complex*16,     intent(in) :: coup(nst,nst,dimph)
type(mps),      intent(in) :: metal_ps
real*8,         intent(in) :: miu(nLead), telec
real*8,         intent(in) :: fock0(norbs,norbs)
type(dBlkSeq),intent(inout):: Lambda(nLead)
complex*16,     intent(out):: cselfr(norbs,norbs)
complex*16,     intent(out):: cselfl(norbs,norbs)
complex*16,     intent(out):: cselfg(norbs,norbs)
logical,        intent(in) :: lwrite
logical,        intent(in) :: labsorb
logical,        intent(in) :: lhot
real*8,         intent(out):: absorb(dimph)
real*8,         intent(out):: flux(2)    ! 
real*8,         intent(out):: hot(nst,2) 
!
real*8,     allocatable :: rwork(:)
real*8,     allocatable :: tlamda(:,:)! total lambda of the electrodes
complex*16, allocatable :: grs(:,:)   ! retarded GF of the system
complex*16, allocatable :: gl0(:,:)   ! lesser   GF of the metal
complex*16, allocatable :: gg0(:,:)   ! greater  GF of the metal
complex*16, allocatable :: msMat(:,:) ! metal-system coupling
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: selfsl(:,:),selfsg(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:),greenl(:,:)
complex*16, allocatable :: dgreenr(:,:),dgreeng(:,:),dgreenl(:,:)
!
integer    :: ipiv(nst), info
integer    :: i,ni, nj, nk, kd, istat
real*8     :: miuG,diff(2),dtmp
real*8     :: weight
real*8     :: fermi 
complex*16 :: ztmp1,ztmp2
complex*16 :: ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(rwork(2*nst*nst))
allocate(grs(norbs,norbs))
allocate(gl0(nst,nst),gg0(nst,nst))
allocate(tlamda(norbs,norbs))
allocate(hmat(nst,nst))
allocate(selfsl(nst,nst),selfsg(nst,nst))
allocate(ctmp1(nst,nst),ctmp2(nst,nst),ctmp3(nst,nst))
allocate(greenr(nst,nst),greeng(nst,nst),greenl(nst,nst))
allocate(dgreenr(nst,nst),dgreeng(nst,nst),dgreenl(nst,nst))

!-----------------------------------------------------
! this subroutine assumes number of states in metal  !
! is larger than that of system.                     !
!-----------------------------------------------------

if(nst<norbs) then
  write(6,*) 'warning: too less states in metal!'
  stop
endif
!
!miuG=0.5d0
miuG=(miu(1)+miu(2))/2.d0
!
dgreenr = czero
dgreenl = czero
dgreeng = czero

hmat=czero
do i=1,nst
  hmat(i,i)=dcmplx(metal_ps%eig(i),-metal_ps%relax)
enddo

!-------------------------------
! system-metal coupling matrix !
!-------------------------------
allocate(msMat(norbs,nst))
msMat=czero
do nj=1,nst
  do ni=1,norbs

    !if(nj<=(nst/2).and.ni>(norbs/2)) cycle
    !if(nj>(nst/2).and.ni<=(norbs/2)) cycle
    if(metal_ps%eig(nj)<=0.d0.and.ni>(norbs/2)) cycle
    if(metal_ps%eig(nj)>0.d0.and.ni<=(norbs/2)) cycle
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

!----------------------------------------------------------------------
! G^< and G^> of metal                                                !
! G^r=[E - H_m - Sigma_S]^{-1},  where Sigma_S = V_{mi}G_{S,i} V_{im} !
! G_S is the GF of the system                                         !
!----------------------------------------------------------------------
grs=czero
call self_lead(ePt,nLead,norbs,miu,fock0,Lambda,tlamda)
call self_s(ePt,nLead,norbs,nst,miu,telec,fock0,Lambda,tlamda,msMat,ctmp2,selfsl,selfsg,.true.)

greenr= -hmat - ctmp2
do ni = 1, nst
  greenr(ni,ni) = ePt + greenr(ni,ni)
enddo
!
call zgetrf(nst,nst,greenr,nst,ipiv,info)
call zgetri(nst,greenr,nst,ipiv,ctmp3,nst*nst, info)
  
fermi=FermiDirac(telec,dble(ePt),miuG)

!gl0(1:nst,1:nst)=-2.d0*fermi*eye*dimag(greenr(1:nst,1:nst))
!gg0(1:nst,1:nst)=-2.d0*(fermi-1.d0)*eye*dimag(greenr(1:nst,1:nst))

do ni=1,nst
  selfsl(ni,ni)=selfsl(ni,ni) + 2.d0*fermi*eye*metal_ps%relax
  selfsg(ni,ni)=selfsg(ni,ni) + 2.d0*(fermi-1.d0)*eye*metal_ps%relax
enddo
!
call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,selfsl,nst,czero,ctmp1,nst)
call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,greenr,nst,czero,gl0,nst)
!
call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,selfsg,nst,czero,ctmp1,nst)
call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,greenr,nst,czero,gg0,nst)
!

!----------------------
!    \Sigma_ps(E)     !
!----------------------
absorb=0.d0
hot=0.d0
do nk = 1, dimph
do kd = 1, 2
  if(kd==1) then
    ePtq = ePt + dcmplx(pheigv(nk),0.d0)
  else 
    ePtq = ePt - dcmplx(pheigv(nk),0.d0)
  endif

  fermi=FermiDirac(telec,dble(ePtq),miuG)

  !------------------------------------
  ! self-energy induced by the system !
  !------------------------------------ 
  grs=czero
  call self_lead(ePtq,nLead,norbs,miu,fock0,Lambda,tlamda)
  call self_s(ePtq,nLead,norbs,nst,miu,telec,fock0,Lambda,tlamda,msMat,ctmp2,greenl,greeng,.true.)

  !---------------------
  ! G^r(epsilon+omega) !
  !---------------------
  greenr= - hmat - ctmp2
  do ni = 1, nst
     greenr(ni,ni) = ePtq + greenr(ni,ni)
  enddo

  call zgetrf(nst,nst,greenr,nst, ipiv, info)
  call zgetri(nst,greenr,nst,ipiv,ctmp2,nst*nst,info) 

  !----------------------------------------------
  ! calculate lesser GF from lesser self-energy !
  !----------------------------------------------
  do ni=1,nst
    greenl(ni,ni)=greenl(ni,ni) + 2.d0*fermi*eye*metal_ps%relax
    greeng(ni,ni)=greeng(ni,ni) + 2.d0*(fermi-1.d0)*eye*metal_ps%relax
  enddo

  call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,greenl,nst,czero,ctmp1,nst)
  call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,greenr,nst,czero,greenl,nst)

  call zgemm('n','n',nst,nst,nst,cunity,greenr,nst,greeng,nst,czero,ctmp1,nst)
  call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,greenr,nst,czero,greeng,nst)

  !-----------------------
  ! calculate absorption !
  !-----------------------
  if(kd==2.and.labsorb) then
    ! M* [G^>(E)*M*G^<(E-w) - G^<(E)*M*G^>(E-w)]
    call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greenl,nst,czero,ctmp1,nst)
    call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greeng,nst,czero,ctmp2,nst)

    call zgemm('n','n',nst,nst,nst,cunity,gg0,nst,ctmp1,nst,czero,ctmp3,nst)
    call zgemm('n','n',nst,nst,nst,-1.d0*cunity,gl0,nst,ctmp2,nst,cunity,ctmp3,nst)

    call zgemm('c','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,ctmp3,nst,czero,ctmp1,nst)

    do ni=1,nst
      absorb(nk)=absorb(nk) + dble(ctmp1(ni,ni))*bocc(nk)
    enddo

    if(lhot) then
      ! M G^<(E-w) M
      call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greenl,nst,czero,ctmp1,nst)
      call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp2,nst)

      ! M G^>(E-w) M
      call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greeng,nst,czero,ctmp1,nst)
      call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp3,nst)

      do ni=1,nst
        do nj=1,nst
          hot(ni,1) = hot(ni,1) + dble(gg0(ni,nj)*ctmp2(nj,ni) - gl0(ni,nj)*ctmp3(nj,ni))
          dtmp = dble(gg0(ni,nj)*ctmp2(nj,ni) - gl0(ni,nj)*ctmp3(nj,ni))
          !if(dtmp<0.d0) then
          !  write(6,'(A,f9.3,2I5,6e15.7)') 'ePt=', dble(ePt), nj,ni, dtmp, &
          !         FermiDirac(telec,dble(ePtq),miuG)-FermiDirac(telec,dble(ePt),miuG),&
          !         dimag(gg0(nj,ni)),dimag(greenl(nj,ni)),dimag(gl0(nj,ni)),dimag(greeng(nj,ni))
          !endif
        enddo
      enddo

      call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,gg0,nst,czero,ctmp1,nst)
      call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp2,nst)

      call zgemm('n','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,gl0,nst,czero,ctmp1,nst)
      call zgemm('n','c',nst,nst,nst,cunity,ctmp1,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp3,nst)

      do ni=1,nst
      do nj=1,nst
        hot(ni,2) = hot(ni,2) + dble(greenl(ni,nj)*ctmp2(nj,ni) - greeng(ni,nj)*ctmp3(nj,ni))
      enddo
      enddo

      !do ni=1,nst
      !  do nj=1,nst
      !    hot(nj,1) = hot(nj,1) + dble(gg0(nj,nj)*coup(nj,ni,nk)*greenl(ni,ni)*coup(ni,nj,nk) - &
      !                                 gl0(nj,nj)*coup(nj,ni,nk)*greeng(ni,ni)*coup(ni,nj,nk))
      !    hot(ni,2) = hot(ni,2) + dble(gg0(nj,nj)*coup(nj,ni,nk)*greenl(ni,ni)*coup(ni,nj,nk) - &
      !                                 gl0(nj,nj)*coup(nj,ni,nk)*greeng(ni,ni)*coup(ni,nj,nk))
      !  enddo
      !enddo
    endif
  endif

  if(kd==1) then
    ctmp1 = bocc(nk)*greenr - 0.5d0*greenl
  else
    ctmp1 = (1.d0+bocc(nk))*greenr + 0.5d0*greenl
  endif

  call zgemm('c','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,ctmp1,nst,czero,ctmp2,nst)
  call zgemm('n','n',nst,nst,nst,cunity,ctmp2,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp3,nst)

  !call zlarcm(nst,nst,coup(1:nst,1:nst,nk),nst,ctmp1,nst,ctmp2,nst,rwork)
  !call zlacrm(nst,nst,ctmp2,nst,coup(1:nst,1:nst,nk),nst,ctmp3,nst,rwork)

  dgreenr = dgreenr + ctmp3

  !--------------------------------------
  !  retarded and lesser self-energy    !
  !--------------------------------------
  if(kd==1) then
    !ztmp1 = dcmplx(0.d0, 0.d0)
    ztmp1 = dcmplx(1.d0+bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  endif

  call zgemm('c','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greenl,nst,czero,ctmp2,nst)
  call zgemm('n','n',nst,nst,nst,cunity,ctmp2,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp1,nst)
  !call zlarcm(nst,nst,coup(1:nst,1:nst,nk),nst,greenl,nst,ctmp2,nst,rwork)
  !call zlacrm(nst,nst,ctmp2,nst,coup(1:nst,1:nst,nk),nst,ctmp1,nst,rwork)
  ctmp1=ctmp1*ztmp1

  dgreenl = dgreenl + ctmp1

  !-----------------------
  ! greater self-energy  !
  !-----------------------
  if(kd==1) then
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk)+1.d0, 0.d0)
  endif

  call zgemm('c','n',nst,nst,nst,cunity,coup(1:nst,1:nst,nk),nst,greeng,nst,czero,ctmp2,nst)
  call zgemm('n','n',nst,nst,nst,cunity,ctmp2,nst,coup(1:nst,1:nst,nk),nst,czero,ctmp1,nst)
  !call zlarcm(nst,nst,coup(1:nst,1:nst,nk),nst,greeng,nst,ctmp2,nst,rwork)
  !call zlacrm(nst,nst,ctmp2,nst,coup(1:nst,1:nst,nk),nst,ctmp1,nst,rwork)
  ctmp1=ctmp1*ztmp1

  dgreeng = dgreeng + ctmp1
enddo !kd
enddo !nk

dgreenr = 0.5d0 * (dgreeng-dgreenl)    

!------------------------------
! G^r_M(E)\Sigma_{ps} G^a_M   !
!------------------------------

!------------------------------------
! self-energy induced by the system !
!------------------------------------ 
grs=czero
call self_lead(ePt,nLead,norbs,miu,fock0,Lambda,tlamda)
call self_s(ePt,nLead,norbs,nst,miu,telec,fock0,Lambda,tlamda,msMat,ctmp2,selfsl,selfsg,.true.)

ztmp2=czero
greenr= -hmat - ctmp2
do ni = 1, nst
  greenr(ni,ni) = ePt + greenr(ni,ni)
  ztmp2 = ztmp2 + ctmp2(ni,ni)
enddo
!
call zgetrf(nst,nst,greenr,nst,ipiv,info)
call zgetri(nst,greenr,nst,ipiv,ctmp2,nst*nst, info)

ztmp1=czero
do ni=1,nst
  ztmp1=ztmp1-greenr(ni,ni)
enddo

if(lwrite) write(70,'(f9.3,2e15.7)') dble(ePt), dimag(ztmp1),dimag(ztmp2)
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

if(lhot) then
  !hot=0.d0
  !do i=1,nst
  !   hot(i,1) =  -dble(dgreenl(i,i) * eye)
  !   hot(i,2) =   dble(dgreeng(i,i) * eye)
  !enddo
endif


flux=0.d0
! Sigma^<_S G^>_M
call zgemm('n','n',nst,nst,nst,cunity,selfsl,nst,dgreeng,nst,czero,ctmp2,nst)
call zgemm('n','n',nst,nst,nst,cunity,selfsg,nst,dgreenl,nst,czero,ctmp3,nst)

do ni=1,nst
  flux(1) = flux(1) + ctmp2(ni,ni) 
  flux(2) = flux(2) - ctmp3(ni,ni) 
enddo


if(lwrite) write(90,'(f9.4,8e15.7)') dble(ePt),dgreenl(1,1),dgreenl(2,2)

!------------------------
! V^\dag \delta g^< V   !
!------------------------
call zgemm('n','n',norbs,nst,nst,cunity,msMat,norbs,dgreenr,nst,czero,ctmp1,nst)
call zgemm('n','c',norbs,norbs,nst,cunity,ctmp1,nst,msMat,norbs,czero,cselfr,norbs)

call zgemm('n','n',norbs,nst,nst,cunity,msMat,norbs,dgreenl,nst,czero,ctmp1,nst)
call zgemm('n','c',norbs,norbs,nst,cunity,ctmp1,nst,msMat,norbs,czero,cselfl,norbs)

call zgemm('n','n',norbs,nst,nst,cunity,msMat,norbs,dgreeng,nst,czero,ctmp1,nst)
call zgemm('n','c',norbs,norbs,nst,cunity,ctmp1,nst,msMat,norbs,czero,cselfg,norbs)
!

deallocate(rwork)
deallocate(msMat)
deallocate(selfsl,selfsg)
deallocate(greenr,greeng,greenl,hmat) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
deallocate(dgreenr,dgreenl,dgreeng)
deallocate(grs,tlamda,gl0,gg0)

end subroutine self_ps_new
