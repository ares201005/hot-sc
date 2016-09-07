module transphonon


logical, save :: lphoton
logical, save :: tplasmon
integer, save :: mode

real*8, save :: pt_telec,pt_beta

real*8, save :: epscoup,epcoup, phfreq         

real*8, allocatable, save :: pheigv(:) !phonon frequency
!real*8, allocatable, save :: coup(:,:,:) ! e-p coupling matrix
complex*16,allocatable, save :: coup(:,:,:)

real*8, allocatable, save :: bocc(:)     ! phonon occupation number

! matrixs needed in the e-p calculation
!real*8, allocatable, save :: sigma0(:,:)
complex*16, allocatable, save :: sigma0(:,:)

!

!
end module transphonon
