module transmod
!
real*8, save  :: telec,beta
real*8, save  :: ampL, ampR, ampG, miu0, ebot

integer, save                 :: fieldtyp
!                                fieldtyp=0   : step function
!                                fieldtyp=1   : exponential
!                                fieldtyp=2   : half-sin
!                                fieldtyp=3   : sin
!                                fieldtyp=4   : custom

complex*16, allocatable, save :: newhmat(:,:)
real*8, allocatable, save     :: deltah(:,:)
real*8, allocatable, save     :: LamdaL(:,:), LamdaR(:,:)
real*8, allocatable, save     :: deltapt(:, :)

!
!* external field 
integer, save :: idelta(5), ipot, ncyc, nsig, ndirect(5)
real*8, save :: strength(5), freq(5), tp, tlas, taulas, e0
real*8, save :: ex, ey, ez
logical, save :: lefield

!*
real*8, save :: charge

type :: dBlkSeq
  real*8, allocatable :: g0(:,:)
  real*8, allocatable :: ge(:,:)
  real*8 :: width
end type dBlkSeq
!

type :: mps
  real*8 :: barrier_e  ! barrier for electron
  real*8 :: barrier_h  ! barrier for hole
  real*8 :: mscoup     ! metal-system coupling
  real*8 :: relax      ! relaxation time of metal
  integer:: nst        ! number of states in the metallic contatcts
  !
  real*8, allocatable :: eig(:) ! energy of states
end type mps

end module transmod
