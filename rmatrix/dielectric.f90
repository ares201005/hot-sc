!
! Reference: Applied Optics 37, 5271 (1998).
! 
! A. Lorentz-Drude Model
! 
! epsilon(w) = epsilon_f(w) + epsilon_b(w)
!
! This separates explicitly the intraband effects (ususally referred to as free
! electron effects) from interband effects (ususally referred to as bound electron effects). 
!
! where
!                      \Omega^2_p
! epsilon_f(w) = 1 - --------------, where \Omega_p = \sqrt{f_0}*omegap
!                    w ( w - i G_0 )
!
!                               f_i w^2_p
! epsilon_b(w) = \sum_j  ------------------------
!                        (w^2_j - w^2) + i w G_j
!
! B. Brendel-Bormann Model
!
!                     \Omega^2_p
! epsilon_f(w) = 1 - ------------- + \sum_j x_j(w),
!                     w(w - i G_0)
! where
!
! x_j(w) = 
!
subroutine dielectric(model, metal, omega, epsilonw)
implicit none
include '../include/parameter'
!
integer,     intent(in) :: model         ! choice of dielectric models
character*2, intent(in) :: metal         ! choice of metal
real*8,      intent(in) :: omega         ! frequency 
complex*16,  intent(out):: epsilonw  ! dielectric at given frequency
!
integer, parameter :: nLB = 5
!
real*8  :: dtmp
real*8  :: omegap, f0, G0, omegaj(nLB), fj(nLB), Gj(nLB)
real*8  :: omegac, A, delta
integer :: i, j, k

epsilonw = cZero

! parameters
!
if ( metal=='Ag' ) then

 if(model==1) then
   omegap = 9.01
   f0 = 0.845
   G0 = 0.048
   data fj(1:nLB)     /0.065, 0.124, 0.011, 0.840, 5.646/
   data Gj(1:nLB)     /3.886, 0.452, 0.065, 0.916, 2.419/
   data omegaj(1:nLB) /0.816, 4.481, 8.185, 9.083, 20.29/
 else if(model==2) then
   omegap = 9.08d0
   f0 = 4.18d0 ! dielectric of backgroud
   G0 = 0.06d0 ! damping factor
 endif

else if( metal == 'Au' ) then

 if(model==1) then
   omegap = 9.03
   f0 = 0.760
   G0 = 0.053
   data fj(1:nLB)     /0.024, 0.010, 0.071, 0.601, 4.384/
   data Gj(1:nLB)     /0.241, 0.345, 0.870, 2.494, 2.214/
   data omegaj(1:nLB) /0.415, 0.830, 2.969, 4.304, 12.32/
 else if(model==2) then
   ! ref: acs nano paper 2014
   omegap = 9.08d0
   f0 = 4.18d0 ! dielectric of backgroud
   G0 = 0.06d0 ! damping factor
 else if(model==3) then
   ! ref: DOI 10.1007/s11468-015-0128-7
   omegap = 9.10d0
   f0 = 9.84d0 ! dielectric of backgroud
   G0 = 0.072d0 ! damping factor
   A = 5.6d0
   omegac=2.4d0
   delta=0.17d0
 endif
else
  write(6,*) 'warning: metal is not fund!'
  stop
endif
!
if(model==1) then
  !write(6,*) 'test1'
  epsilonw = cunity - f0 * omegap * omegap / omega / (omega - eye * G0) 
  do i = 1, nLB
    epsilonw = epsilonw + fj(i) * omegap * omegap / (omegaj(i)*omegaj(i) - omega * omega + eye * omega * Gj(i) )
  enddo
else if(model==2) then
  epsilonw = cunity * f0 - omegap * omegap / (omega*(omega+eye*G0))
else if (model==3) then
  epsilonw = cunity * f0 - omegap * omegap / (omega*(omega+eye*G0)) + eye * A/(1.d0+exp(-(omega-omegac)/delta))
endif
!
end subroutine dielectric
