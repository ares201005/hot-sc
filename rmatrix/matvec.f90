subroutine matvec(m, n, const, h, x, y)

implicit none

include '../include/parameter'
integer, intent(in) :: n, m
complex*16, intent(in) :: x(n), h(m,m)
complex*16, intent(in) :: const
complex*16, intent(out) :: y(n)

integer :: i, j, k, nk

complex*16 :: ztmp1

complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:)

if(n/=m*m) then
  write(6,*) ' warning: n is not equal to m * m!'
  stop
endif

allocate( ctmp1(m,m), ctmp2(m,m) )

y = czero
!
do i = 1, m
  do j = 1, m
    nk = (i-1) * m + j
    ctmp1(j,i) = x(nk)
  enddo
enddo


call zgemm('n', 'n', m, m, m, -1.d0*cunity, H, m, ctmp1, m, czero, ctmp2, m)
call zgemm('n', 'c', m, m, m, cunity, ctmp1, m, H, m, cunity, ctmp2, m)

do i = 1, m
  do j = 1, m
    nk = (i-1) * m + j
    ctmp2(j,i) = const * ctmp1(j,i) + ctmp2(j,i)
    y(nk) = ctmp2(j,i)
  enddo
enddo


deallocate( ctmp1, ctmp2 )

end subroutine matvec 
