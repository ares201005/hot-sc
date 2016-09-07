subroutine dmaxmat(dim1, dim2, dmat, ldm, dmax)
implicit none
!
integer, intent(in) :: dim1, dim2, ldm
double precision, intent(in)  :: dmat(ldm,*)
double precision, intent(out) :: dmax
!
integer ni,nj
!
dmax = 0.d0
do nj=1,dim2
 do ni=1,dim1
  dmax = max(dmax, dabs(dmat(ni,nj)))
 enddo
enddo
!
end subroutine dmaxmat
!----------------------------------------------
subroutine dmaxmat2(dim1, dim2, dmat, ldm, dmax, irow, icol)
implicit none
!
integer, intent(in)  :: dim1, dim2, ldm
integer, intent(out) :: irow, icol
double precision, intent(in)  :: dmat(ldm,*)
double precision, intent(out) :: dmax
!
integer ni,nj
!
dmax = 0.d0
irow = 0
icol = 0
do nj=1,dim2
 do ni=1,dim1
  if ( dmax .le. dabs(dmat(ni,nj)) )  then
   irow = ni
   icol = nj
   dmax = dabs(dmat(ni,nj))
  endif
 enddo
enddo
!
end subroutine dmaxmat2
!-----------------------------------------------
subroutine cmaxmat(dim1, dim2, cmat, ldm, dmax)
implicit none
!
integer, intent(in) :: dim1, dim2, ldm
complex*16, intent(in)  :: cmat(ldm,*)
double precision, intent(out) :: dmax
!
integer ni,nj
!
dmax = 0.d0
do nj=1,dim2
 do ni=1,dim1
  dmax = max(dmax, cdabs(cmat(ni,nj)))
 enddo
enddo
!
end subroutine cmaxmat
