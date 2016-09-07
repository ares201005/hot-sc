subroutine landauerphscba(nLead,norbs,nmode,fock0,Lambda,miu,telec,pheigv,coup)
use transmod, only: dBlkSeq
use transphonon, only: bocc
implicit none
include '../include/parameter'
!
! Calculate steady state current within WBL and electron-photon coupling
! by integration of transmission function. There are two parts of
! transmission function, one is the elastic part:
!  T1(E) = 4 Tr( G^r * Lambda_R * G^a * Lambda_L )
!
! another is the inelastic part:
!  T2(E) = \Sigma_L G^r(E)\Sigma_{ep}(E)G^a(E) - 
!
! Current is  int [(f_L-f_R)T1(E) + T2(E)]
!
! note: there is not fermi function before T2 part,
!       so integration over total energy range is needed
!!!!!! 
!  Lowest order expansion method is used
!****************************************************************
integer,       intent(in)    :: nLead,norbs,nmode
type(dBlkSeq), intent(inout) :: Lambda(nLead)
real*8,        intent(in)    :: fock0(norbs,norbs)
real*8,        intent(in)    :: miu(nLead),telec
real*8,        intent(in)    :: pheigv(nmode)
complex*16,    intent(in)    :: coup(norbs,norbs,nmode)
!
character*15 :: tFile
!
integer:: iLead,jLead
integer:: i,ni,nj,nk,np,step,iwarn,istat,info,ipiv(norbs)
integer:: maxits, iter
real*8 :: etop,ebot
real*8 :: dengy, tmp2,tmp1,telec0, currL2, coeff(2)
real*8 :: fermi(nLead),curr(nLead),curr2(nLead)
real*8 :: pstep, dtmp, dtmpL, dtmpR 
real*8 :: derr1, derr2, dmcon
real*8 :: tCoef(nLead)
!
real*8, external :: FermiDirac
!
integer :: npoint,nsplit
integer :: ngrid
real*8 :: downrange=0.d0, uprange=0.d0,eImag=1.e-6
namelist /land/ ngrid, downrange, uprange,eImag
!*
complex*16 :: ePt
real*8,     allocatable :: tlamda(:,:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: greenrsav(:,:,:),greenr(:,:,:)
complex*16, allocatable :: greenlsav(:,:,:),greenl(:,:,:)
complex*16, allocatable :: ctmp1(:, :), ctmp2(:, :)
complex*16, allocatable :: ctmp3(:, :), ctmp4(:, :)
complex*16, allocatable :: cselfr(:,:), cselfl(:,:),cselfg(:,:)
!*
write(6,*) 
write(6,*) ' =========== entering landauerphscba ================'
write(6,*) ' calculate steady state current with e-p coupling'
!
!telec0 = 1.5d0 * pheigv(1)
telec0 = telec
!telec0 = 0.d0

eImag=1.d-4
rewind(5)
read(5,nml=land)
call flush(6)
!
etop = maxval(miu) + uprange     !end point of energy
ebot = minval(miu) - downrange  !initial point of energy
ngrid=int((etop-ebot)/0.01d0*ngrid)
dengy = (etop - ebot) / dble(ngrid+1)

write(6,*) ' integration range:', ebot, etop
write(6,*) ' number of grid:   ', ngrid
write(6,*) ' integration step: ', dengy
!
allocate(hmat(norbs,norbs),tlamda(norbs,norbs))
allocate(greenr(norbs,norbs,ngrid), greenrsav(norbs,norbs,ngrid),STAT=istat)
allocate(greenl(norbs,norbs,ngrid), greenlsav(norbs,norbs,ngrid),STAT=istat)
allocate(ctmp1(norbs, norbs), ctmp2(norbs, norbs), stat=istat)
allocate(ctmp3(norbs, norbs), ctmp4(norbs, norbs), STAT=istat)
allocate(cselfr(norbs,norbs),cselfl(norbs,norbs),cselfg(norbs,norbs) )
!
!
if (ngrid .le. 1) then
 write(6,*)
 write(6,*)' error! ngrid = ',ngrid
 stop
endif
!
curr = 0.d0
curr2 = 0.d0


!
dmcon = 1.d-10
maxits = 50
iter = 0
greenr = czero
greenl = czero
scba: do
  do nk=1,ngrid
    ePt = dcmplx(ebot+dble(nk-1)*dengy,eImag)

    !energy-dependent Gamma !
    !------------------------
    tlamda=0.d0
    do i=1,nLead
      tmp1=Lambda(i)%width
      tmp1=tmp1*tmp1/((dble(ePt)-miu(i))*(dble(ePt)-miu(i))+tmp1*tmp1)

      if(abs(dble(ePt))<0.5d0*abs(fock0(2,2)-fock0(1,1))) then
        tmp1=1.d-6
      else
        tmp1=1.d0
      endif

      Lambda(i)%ge=Lambda(i)%g0*tmp1
      tlamda=tlamda+Lambda(i)%ge
    enddo

    hmat=dcmplx(fock0,-tlamda)

    !--------------------------------
    ! get self-energies at nk point !
    !--------------------------------
    cselfr = czero
    cselfl = czero

    if(iter>0) then
      do np = 1, nmode
        step = int(pheigv(np)/dengy)

        ctmp1 = czero
        ctmp2 = czero

        if((nk+step)<=ngrid) then
          ctmp1 = ctmp1 + bocc(np) * greenr(:,:,nk+step) - 0.5d0 * greenl(:,:,nk+step)
          ctmp2 = ctmp2 + (1.d0+bocc(np))*greenl(:,:,nk+step)
        endif

        if((nk-step)>=1) then
          ctmp1 = ctmp1 + (1.d0+bocc(np))*greenr(:,:,nk-step) + 0.5d0 * greenl(:,:,nk-step)
          ctmp2 = ctmp2 + bocc(np)*greenl(:,:,nk-step)
        endif

        !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,np),ctmp1, cunity, ctmp3)
        !call cdmulti(norbs,norbs,ctmp3,coup(1:norbs,1:norbs,np), cunity, ctmp4)
        call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,np),norbs,ctmp1,norbs,czero,ctmp3,norbs)
        call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp3,norbs,coup(1:norbs,1:norbs,np),norbs,czero,ctmp4,norbs)
        cselfr = cselfr + ctmp4

        !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,np),ctmp2, cunity, ctmp3)
        !call cdmulti(norbs,norbs,ctmp3,coup(1:norbs,1:norbs,np), cunity, ctmp4)
        call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,np),norbs,ctmp2,norbs,czero,ctmp3,norbs)
        call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp3,norbs,coup(1:norbs,1:norbs,np),norbs,czero,ctmp4,norbs)
        cselfl = cselfl + ctmp4
      enddo
    endif

    !------------------------------
    !  retarded green's function  ! 
    !------------------------------
    greenrsav(:,:,nk)=-hmat(:,:)-cselfr(:,:)
    do ni=1,norbs
        greenrsav(ni,ni,nk) = ePt + greenrsav(ni,ni,nk)
    enddo

    !---------------------------------
    !zero-order green's function G^r !
    !---------------------------------
    call zgetrf(norbs, norbs, greenrsav(1:norbs,1:norbs,nk), norbs, ipiv, info)
    call zgetri(norbs, greenrsav(1:norbs,1:norbs,nk), norbs, ipiv, ctmp2, norbs*norbs, info) 
    !write(6,*) "green function is obtained"

    !--------------------
    ! total self-energy !
    !--------------------
    ctmp2 = cselfl
    do jLead=1,nLead
      fermi(jLead)=FermiDirac(telec,dble(ePt),miu(jLead))
      ctmp2=ctmp2+2.d0*eye*fermi(jLead)*Lambda(jLead)%ge
    enddo

    !---------------------------
    ! lesser green's function  !
    !---------------------------
    call zgemm('n','n', norbs, norbs, norbs, cunity, greenrsav(1:norbs,1:norbs,nk), norbs, ctmp2, norbs,&
             czero, ctmp1, norbs)
    call zgemm('n','c', norbs, norbs, norbs, cunity, ctmp1, norbs, greenrsav(1:norbs,1:norbs,nk), norbs, &
             czero, ctmp3, norbs)  ! lesser green's function
    greenlsav(1:norbs,1:norbs,nk) = ctmp3(1:norbs,1:norbs)
  enddo 

  !----------------------
  !  check convergence  !
  !----------------------
  derr1 = 0.d0
  derr2 = 0.d0
  do nk = 1, ngrid
    do ni = 1, norbs
      do nj = 1, norbs
        derr1 = max(derr1,cdabs(greenrsav(nj,ni,nk)-greenr(nj,ni,nk)))
        derr2 = max(derr2,cdabs(greenlsav(nj,ni,nk)-greenl(nj,ni,nk)))
      enddo
    enddo
  enddo
  greenr = greenrsav
  greenl = greenlsav
  write(6,'(/,X,A20,I5,2e15.7)') "iter/derr1/derr2:",iter, derr1,derr2
  if(derr1<dmcon.and.derr2<dmcon) then
    write(6,*) 'converged!'
    exit
  endif
  iter = iter + 1
  if(iter>maxits) then
    write(6,*) 'not converged!'
    exit
  endif
enddo scba

!---------------------------------------
! integration for current starts here  !
!---------------------------------------
do iLead=1,nLead

  write(tFile,'(A9,I1,A4)') "Tcoef-pt-",iLead,".dat"
  open(89, file=tFile, form='formatted', status='replace')

  do nk = 1, ngrid
    ePt = dcmplx(ebot+dble(nk-1)*dengy,eImag)
    fermi(iLead)=FermiDirac(telec,dble(ePt),miu(iLead))

    !------------------------
    !energy-dependent Gamma !
    !------------------------
    tlamda=0.d0
    do i=1,nLead
      tmp1=Lambda(i)%width
      tmp1=tmp1*tmp1/((dble(ePt)-miu(i))*(dble(ePt)-miu(i))+tmp1*tmp1)

      if(abs(dble(ePt))<0.5d0*abs(fock0(2,2)-fock0(1,1))) then
        tmp1=1.d-6
      else
        tmp1=1.d0
      endif
      !if(i==1.and.iLead==1) write(6,'(2e15.7)') dble(ePt), tmp1

      Lambda(i)%ge=Lambda(i)%g0*tmp1
      tlamda=tlamda+Lambda(i)%ge
    enddo

    !----------------------------------------
    ! for the T1 part, right to left        !
    ! LambdaR -> ctmp2 (lower half)         !
    !----------------------------------------

    do jLead=1,nLead
      if(jLead==iLead) cycle
      fermi(jLead)=FermiDirac(telec,dble(ePt),miu(jLead))
      ctmp2=dcmplx(Lambda(jLead)%ge,0.d0)
      !
      ! G^r* LambdaR (lower, right)  -> ctmp3
      call zgemm('n','n',norbs,norbs,norbs,cunity,greenr(1:norbs,1:norbs,nk),norbs,ctmp2,norbs,czero,ctmp3,norbs)
      ! ctmp3 * G^a = ctmp3 * (G^r)^H -> ctmp2
      call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp3,norbs,greenr(1:norbs,1:norbs,nk),norbs,czero,ctmp2,norbs)

      ctmp1=dcmplx(Lambda(iLead)%ge,0.d0)

      !
      ! ctmp2 * LambdaL (lower, right) -> ctmp3
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

    !-----------------------------------
    ! for T2 part, to left             !
    ! get cselfg and cselfl first      !
    !-----------------------------------

    cselfr = czero
    cselfl = czero
    do np = 1, nmode
      step = int(pheigv(np)/dengy)
      ctmp1 = czero
      ctmp2 = czero
      if((nk+step)<=ngrid) then
        ctmp1 = ctmp1 + bocc(np) * greenr(:,:,nk+step) - 0.5d0 * greenl(:,:,nk+step)
        ctmp2 = ctmp2 + (1.d0+bocc(np))*greenl(:,:,nk+step)
      endif
      if((nk-step)>=1) then
        ctmp1 = ctmp1 + (1.d0+bocc(np))*greenr(:,:,nk-step) + 0.5d0 * greenl(:,:,nk-step)
        ctmp2 = ctmp2 + bocc(np)*greenl(:,:,nk-step)
      endif

      call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,np),norbs,ctmp1,norbs,czero,ctmp3,norbs)
      call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp3,norbs,coup(1:norbs,1:norbs,np),norbs,czero,ctmp4,norbs)
      !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,np),ctmp1, cunity, ctmp3)
      !call cdmulti(norbs,norbs,ctmp3,coup(1:norbs,1:norbs,np), cunity, ctmp4)
      cselfr = cselfr + ctmp4

      call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,np),norbs,ctmp2,norbs,czero,ctmp3,norbs)
      call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp3,norbs,coup(1:norbs,1:norbs,np),norbs,czero,ctmp4,norbs)
      !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,np),ctmp2, cunity, ctmp3)
      !call cdmulti(norbs,norbs,ctmp3,coup(1:norbs,1:norbs,np), cunity, ctmp4)
      cselfl = cselfl + ctmp4
    enddo
 
    do ni = 1, norbs
      do nj = 1, norbs
        cselfg(nj,ni) = cselfr(nj,ni) - dconjg(cselfr(ni,nj)) + cselfl(nj,ni)
      enddo
    enddo

    !* G^r * Sigma^>_ep * G^a * \Sigma^<_L
    call zgemm('n','n',norbs,norbs,norbs,cunity,greenr(1:norbs,1:norbs,nk),norbs,cselfg,norbs,czero,ctmp1,norbs)
    call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr(1:norbs,1:norbs,nk),norbs,czero,ctmp2,norbs)

    ctmp1=eye*fermi(iLead)*Lambda(iLead)%ge
    call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,ctmp1,norbs,czero,ctmp3,norbs)

    ! G^r * Sigma^<_ep * G^a * Sigma^>_L
    call zgemm('n','n',norbs,norbs,norbs,cunity,greenr(1:norbs,1:norbs,nk),norbs,cselfl,norbs,czero,ctmp1,norbs)
    call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr(1:norbs,1:norbs,nk),norbs,czero,ctmp2,norbs)

    ctmp1=eye*(fermi(iLead)-1.d0)*Lambda(iLead)%ge
    call zgemm('n','n',norbs,norbs,norbs,-1.d0*cunity,ctmp2,norbs,ctmp1,norbs,cunity,ctmp3,norbs)
  
    tmp1 = 0.d0 
    do ni = 1, norbs
      tmp1 = tmp1 + dble(ctmp3(ni,ni))/(pi * pi) 
    enddo
    tCoef(nLead)=tmp1*pi

    curr2(iLead) = curr2(iLead) + 1.d0 * pi * tmp1 * dengy

    write(89,'(5e15.7)') dble(ePt), tCoef(1:nLead)
  enddo
  close(89)
enddo

!
curr =  curr * 1.6022d5 * hbarinv
curr2=  curr2* 1.6022d5 * hbarinv

do iLead=1,nLead
  write(6,1000) '  current  through lead: ', iLead, " is ", curr(iLead)
  write(6,1000) '  current2 through lead: ', iLead, " is ", curr2(iLead)
enddo

write(6,*) ' ========== leaving landauerphscba ============='

1000 format(A, I2, A, e15.7, " nA")

deallocate(ctmp1, ctmp2, ctmp3, ctmp4, STAT=istat)
deallocate(hmat,tlamda, greenr, greenl, greenrsav, greenlsav, STAT=istat)
deallocate(cselfr, cselfl,cselfg )
!

end subroutine landauerphscba
