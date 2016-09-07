subroutine landauer(nLead,norbs,miu,beta,fock0,Lambda)
use transmod, only: dBlkSeq
implicit none

include '../include/parameter'
include '../include/trans_data'

integer, intent(in) :: nLead, norbs
real*8,  intent(in) :: miu(nLead)
real*8,  intent(in) :: beta
real*8,  intent(in) :: fock0(norbs,norbs)
type(dBlkSeq),  intent(inout) :: Lambda(nLead)
!
character*10 :: pFile
character*15 :: tFile

integer :: nsplit,npoint
integer :: i, j, k, l, info, ipiv(norbs), ii, jj, kk, ll, ierror
integer :: iLead, jLead
real*8  :: eImag, ebot,etop, dener
real*8  :: tmp, tmp1, tmp2, tmp3
real*8  :: volt, fermi1,fermi2
real*8, allocatable :: curr(:,:)
complex*16 :: ePt
!
real*8, allocatable :: angPt(:), weight(:)
real*8, allocatable :: tlamda(:,:), deltah(:,:), fock(:,:)
complex*16, allocatable  :: greenr(:, :)
complex*16, allocatable  :: ctmp0(:,:), ctmp(:,:), ctmp1(:, :), ctmp2(:,:)
!
integer :: ngrid=10
real*8 :: downrange=0.d0, uprange=0.d0
namelist /land/ ngrid, downrange, uprange,eImag

!
allocate( tlamda(norbs,norbs), deltah(norbs,norbs) )
allocate( greenr(norbs,norbs), fock(norbs,norbs),ctmp0(norbs,norbs) )
allocate( ctmp(norbs,norbs), ctmp1(norbs,norbs), ctmp2(norbs,norbs) )
allocate(curr(nLead,nLead))
!
curr=0.d0
!
ctmp=0.d0
ctmp1=0.d0
ctmp2=0.d0

write(6,*) ' ========= enter  Landauer.f90 ==========='

ngrid=10
eImag=1.d-4
rewind(5)
read(5,nml=land)
call flush(6)


ebot=  minval(miu) - downrange  !initial point of energy
etop = maxval(miu) + uprange     !end point of energy

npoint=mgauleg1
nsplit=int((etop-ebot)/0.1d0*ngrid)
allocate(angPt(nsplit*npoint),weight(nsplit*npoint))

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


pFile='DOS.dat'

do iLead=1,nLead
  do jLead=1,nLead
    if(iLead==jLead) cycle
    write(6,'(A,I2,A,I2,A)') '  from ', jLead, ' lead to ', iLead,' lead'

    write(tFile,'(A5,I1,A1,I1,A4)') "Tcoef",iLead,"-",jLead,".dat"
    open(89, file=tFile, form='formatted', status='replace')
    if(iLead==1.and.jLead==2) then
      open(99, file=pFile, form='formatted', status='replace')
    endif

    volt = miu(jLead)-miu(iLead)

    ! voltage's effect on fock matrix
    deltah = 0.d0
    fock(1:norbs, 1:norbs) = fock0(1:norbs, 1:norbs) + deltah(1:norbs, 1:norbs)

    ! integration over energy space
    integration:do k=1, nsplit*npoint
  
      ePt=dcmplx(angPt(k),eImag)
      dener=weight(k)

      !energy-dependent Gamma
      tlamda=0.d0
      do i=1,nLead
        tmp=Lambda(i)%width
        tmp=tmp*tmp/((dble(ePt)-miu(i))*(dble(ePt)-miu(i))+tmp*tmp)

        if(abs(dble(ePt))<0.5d0*abs(fock0(2,2)-fock0(1,1)) .and. i<3) then
          tmp =1.d-6
        else
          tmp =1.d0
        endif

        Lambda(i)%ge=Lambda(i)%g0*tmp
        tlamda=tlamda+Lambda(i)%ge
      enddo

      !if(iLead==1.and.jLead==2) write(60,'(f9.4,2e15.7)') dble(ePt), Lambda(3)%ge(1,1),Lambda(3)%ge(2,2)
      !if(iLead==1.and.jLead==2) write(70,'(f9.4,2e15.7)') dble(ePt), tlamda(1,1),tlamda(2,2)

      greenr=-dcmplx(fock,-tlamda)      
      do i=1,norbs
        greenr(i,i) = ePt + greenr(i,i)
      enddo
!
      call zgetrf(norbs, norbs, greenr, norbs, ipiv, info)
      call zgetri(norbs, greenr, norbs, ipiv, ctmp, norbs*norbs, info)
!
      if(iLead==1.and.jLead==2) then
        tmp=0.d0
        do i=1,norbs
          tmp=tmp-dimag(greenr(i,i))
        enddo
        write(99,'(f9.4,e15.7)') dble(ePt), tmp
      endif

      tmp1=(dble(ePt)-miu(iLead))*beta
      tmp2=(dble(ePt)-miu(jLead))*beta

      if(tmp1>100.d0) then
        fermi1=0.d0
      else
        fermi1 = 1.d0/(1.d0 + dexp(tmp1))
      endif
      if(tmp2>100.d0) then
        fermi2=0.d0
      else
        fermi2 = 1.d0/(1.d0 + dexp(tmp2))
      endif
!
      ctmp0(1:norbs, 1:norbs) = dcmplx(Lambda(iLead)%ge(1:norbs,1:norbs), 0.d0) 
      call zgemm('c', 'n', norbs, norbs, norbs, cunity, greenr, norbs, ctmp0, norbs, czero, ctmp2, norbs)

      ctmp0(1:norbs, 1:norbs) = dcmplx(Lambda(jLead)%ge(1:norbs,1:norbs), 0.d0) 
      call zgemm('n', 'n', norbs, norbs, norbs, cunity, ctmp0, norbs, ctmp2, norbs, czero, ctmp1, norbs)
      call zgemm('n', 'n', norbs, norbs, norbs, cunity, greenr, norbs, ctmp1, norbs, czero, ctmp, norbs)

      tmp=0.d0
      do i=1,norbs
       tmp = tmp + (dble(ctmp(i,i))/(PI**2.d0))
      enddo 
! 
      tmp2 = 2.d0*PI*tmp*(fermi1-fermi2)*dener

      write(89,'(2e15.7)') dble(ePt), 2.d0*Pi*tmp
      call flush(89)
      curr(iLead,jLead) = curr(iLead,jLead) + tmp2 
    enddo integration

    close(89)
!
    ! 2.d0 comes from spin 
    curr(iLead,jLead) = 2.d0 * curr(iLead,jLead) * 1.6022d5 * hbarinv ! to nA
    if(jLead/=iLead) curr(iLead,iLead) = curr(iLead,iLead)+curr(iLead,jLead)

    write(6,'(A,I2,A,I2,A,e15.7,A)') '  current from', jLead, ' to ', iLead, ' lead is', curr(iLead,jLead), ' nA'

    if(iLead==1.and.jLead==2) close(99)
  enddo
  write(6,'(A,I2,A,e15.7,A)') '  current through', iLead, ' lead is', curr(iLead,iLead), ' nA'
  write(6,*)
enddo

deallocate(angPt,weight)
deallocate( tlamda, greenr, deltah, fock, ctmp0,ctmp, ctmp1, ctmp2) 
!
write(6,*)
write(6,*) ' ========= leave  Landauer.f90 ==========='


end subroutine landauer
