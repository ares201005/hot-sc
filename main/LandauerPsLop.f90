subroutine LandauerPsLop(nLead,norbs,nmode,fock0,Lambda,metal_ps,miu,telec,pheigv,coup)
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
real*8,     allocatable :: angPt(:), weight(:)
real*8,     allocatable :: tlamda(:,:)
real*8,     allocatable :: absorb(:),absorb0(:)
complex*16, allocatable :: hmat(:,:), green(:,:),dgreen(:,:)
complex*16, allocatable :: ctmp1(:, :), ctmp2(:, :)
complex*16, allocatable :: ctmp3(:, :), ctmp4(:, :)
complex*16, allocatable :: cselfr(:,:), cselfl(:,:),cselfg(:,:)
!
character*10 :: pFile
character*15 :: tFile
complex*16   :: ePt
real*8       :: pstep, dtmp, dtmpL, dtmpR 
real*8       :: curr(nLead), curr2(nLead)
real*8       :: maxcurr
!
real*8, external :: FermiDirac
!

write(6,*) 
write(6,*) ' =========== entering LandauerPsLop ================'
write(6,*) ' calculate steady state current with e-p coupling'

if(nLead/=3) then
  write(6,'(A)') "Warning: incorrect nLead!"
endif

ngrid=10
eImag=1.d-4
rewind(5)
read(5,nml=land)
call flush(6)

ebot = minval(miu) - maxval(pheigv) - downrange  !initial point of energy
etop = maxval(miu) + maxval(pheigv) + uprange    !end point of energy

npoint=mgauleg1
nsplit=int((etop-ebot)/0.1d0*ngrid)
allocate(angPt(nsplit*npoint),weight(nsplit*npoint))
  
write(6,*) ' integration range:', ebot, etop
write(6,*) ' number of grid:   ', nsplit*npoint

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
allocate(absorb(nmode),absorb0(nmode))
allocate(tlamda(norbs,norbs))
allocate(hmat(norbs,norbs), green(norbs,norbs), dgreen(norbs,norbs),STAT=istat)
allocate(ctmp1(norbs, norbs), ctmp2(norbs, norbs), stat=istat)
allocate(ctmp3(norbs, norbs), ctmp4(norbs, norbs), STAT=istat)
allocate(cselfr(norbs,norbs),cselfl(norbs,norbs),cselfg(norbs,norbs) )
!
pFile='DOS-ps.dat'

curr=0.d0
curr2 = 0.d0
absorb=0.d0
do iLead=1,nLead
  write(tFile,'(A9,I1,A4)') "Tcoef-ps-",iLead,".dat"
  open(89, file=tFile, form='formatted', status='replace')
  ! 
  if(iLead==1) open(99, file=pFile, form='formatted', status='replace')
  if(iLead==1) open(79, file="absorb.dat", form='formatted', status='replace')

  do nk=1,nsplit*npoint
    ePt = dcmplx(angPt(nk),eImag)
    dengy=weight(nk)

    !------------------------
    !energy-dependent Gamma !
    !------------------------
    tlamda=0.d0
    do i=1,nLead
      tmp1=Lambda(i)%width
      tmp1=tmp1*tmp1/((dble(ePt)-miu(i))*(dble(ePt)-miu(i))+tmp1*tmp1)

      !tmp1=1.d0
      !if(i==1.and.dble(ePt)>fock0(1,1)) tmp1=1.d-5
      !if(i==2.and.dble(ePt)<fock0(2,2)) tmp1=1.d-5

      if(abs(dble(ePt))<0.5d0*abs(fock0(2,2)-fock0(1,1)) .and. i<3) then
        tmp1=1.d-5
      else
        tmp1=1.d0
      endif
     
      if(i==3) call self1(ePt,norbs,metal_ps%nst,miu(3),metal_ps%mscoup,metal_ps%relax,metal_ps%eig,Lambda(i)%g0)

      Lambda(i)%ge=Lambda(i)%g0*tmp1
      tlamda=tlamda+Lambda(i)%ge
      !if(i<3) tlamda=tlamda+Lambda(i)%ge
    enddo

    hmat=dcmplx(fock0,-tlamda)

    if(iLead==1) then 
      call self_ps(ePt,metal_ps%nst,norbs,nmode,pheigv,coup,metal_ps,miu(3),telec0,cselfr,cselfl,cselfg,.true.)
      absorb0=0.d0
      call absorp_ps(ePt,metal_ps%nst,nmode,pheigv,coup,metal_ps,miu(3),telec0,1.d0,absorb0)
      absorb=absorb+absorb0*dengy
      write(79,'(f9.4,e15.7)') dble(ePt), sum(absorb0)
    else
      call self_ps(ePt,metal_ps%nst,norbs,nmode,pheigv,coup,metal_ps,miu(3),telec0,cselfr,cselfl,cselfg,.false.)
    endif
    cSelfr=czero

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
    ! G^r * Sigma^>_ep * G^a * \Sigma^<_L     !
    !------------------------------------------
    call zgemm('n','n',norbs,norbs,norbs,cunity,green,norbs,cselfg,norbs,czero,ctmp1,norbs)
    call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,green,norbs,czero,ctmp2,norbs)

    ctmp1=eye*fermi(iLead)*Lambda(iLead)%ge

    call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,ctmp1,norbs,czero,ctmp3,norbs)

    !-------------------------------------
    ! G^r * Sigma^<_ep * G^a * Sigma^>_L !
    !-------------------------------------
    call zgemm('n','n',norbs,norbs,norbs,cunity,green,norbs,cselfl,norbs,czero,ctmp1,norbs)
    call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,green,norbs,czero,ctmp2,norbs)

    ctmp1=eye*(fermi(iLead)-1.d0)*Lambda(iLead)%ge
    call zgemm('n','n',norbs,norbs,norbs,-1.d0*cunity,ctmp2,norbs,ctmp1,norbs,cunity,ctmp3,norbs)
  
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

absorb=absorb*hbarinv/(2.d0*pi)
maxcurr= sum(absorb)*1.6022d5
write(6,*)
write(6,'(A,e15.7,A)') '  absorped photon flux is ', sum(absorb), ' /fs'
write(6,'(A,e15.7,A)') '  maximum photocurrent is ', maxcurr, ' nA'
write(6,'(A,f10.3,A)') '  IQE is                  ', abs(curr2(1))/maxcurr*1d2,'%'

write(6,*) ' ========== leaving LandauerPsLop ============='

1000 format(A, I2, A, e15.7, " nA")

!
deallocate(absorb0,absorb)
deallocate(tlamda)
deallocate(ctmp1, ctmp2, ctmp3, ctmp4, STAT=istat)
deallocate(hmat, green, dgreen, STAT=istat)
deallocate(cselfr, cselfl,cselfg )
!
end subroutine LandauerPsLop
