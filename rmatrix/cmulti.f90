      subroutine cmulti(dim1,dim0,cin1,cin2,const,cout)
      implicit none
!
! complex matrix multiplication subroutine
! cout = cin1 * cin2 * const
!
      integer ni,nj,nk,dim1,dim0
      complex*16 const,tmp
      complex*16 cin1(dim0,dim0),cin2(dim0,dim0),cout(dim0,dim0)
!
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in cmulti'
       call flush(8)
       stop
      endif
!
      do nj=1,dim1
       do ni=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
       enddo
       do nk=1,dim1
        tmp = cin2(nk,nj) * const
        do ni=1,dim1
         cout(ni,nj) = cout(ni,nj) + tmp * cin1(ni,nk)
        enddo
       enddo
      enddo
!
      return
      end

!
      subroutine cmultimm(dim1,dim0,cin1,cin2,const,beta,cout)
      implicit none
!
! complex matrix multiplication subroutine
! cout = cin1 * cin2 * const+beta*cout
!
      integer ni,nj,nk,dim1,dim0
      complex*16 const,tmp,beta
      complex*16 cin1(dim0,dim0),cin2(dim0,dim0),cout(dim0,dim0)
      complex*16,allocatable:: ctmp(:,:)
!
      allocate(ctmp(dim0,dim0))
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in cmulti'
       call flush(8)
       stop
      endif
      ctmp=cout
!
      do nj=1,dim1
       do ni=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
       enddo
       do nk=1,dim1
        tmp = cin2(nk,nj) * const
        do ni=1,dim1
         cout(ni,nj) = cout(ni,nj) + tmp * cin1(ni,nk)
        enddo
       enddo
      enddo
      
!
      do nj=1,dim0
       do ni=1,dim0
        cout(ni,nj)= cout(ni,nj)+beta*ctmp(ni,nj)
       enddo
      enddo
!
      deallocate(ctmp)
      return
      end
!      
      subroutine dmulti(dim1,dim0,din1,din2,dconst,dout)  
      implicit none
!
! double precision matrix multiplication subroutine
! dout = din1 * din2 * dconst
!     
      integer ni,nj,nk,dim1,dim0
      double precision dconst,dtmp
      double precision din1(dim0,dim0),din2(dim0,dim0),dout(dim0,dim0)
!
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in dmulti'
       call flush(8)
       stop
      endif
!
      do nj=1,dim1
       do ni=1,dim1
        dout(ni,nj) = 0.d0
       enddo
       do nk=1,dim1
        dtmp = dconst * din2(nk,nj)
        do ni=1,dim1
         dout(ni,nj) = dout(ni,nj) + dtmp * din1(ni,nk)
        enddo
       enddo
      enddo
!
      return
      end
!
      subroutine ccopym(dim1,dim0,csource,ctarget)
      implicit none
!
      integer dim1,dim0,ni,nj
      complex*16 csource(dim0,dim0),ctarget(dim0,dim0)
!
      do ni=1,dim1
       do nj=1,dim1
        ctarget(ni,nj) = csource(ni,nj) 
       enddo
      enddo
!
      return
      end

      subroutine dcmulti(dim1,dim0,din1,cin2,const,cout)
      implicit none
!
! matrix multiplication subroutine
! cout = din1 * cin2 * const
!
      integer ni,nj,nk,dim1,dim0
      double precision din1(dim0,dim0)
      complex*16 const
      complex*16 cin2(dim0,dim0),cout(dim0,dim0)
!
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in dcmulti'
       call flush(8)
       stop
      endif
!
      do ni=1,dim1
       do nj=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
        do nk=1,dim1
         cout(ni,nj) = cout(ni,nj) + din1(ni,nk) * cin2(nk,nj)
        enddo
        cout(ni,nj) = cout(ni,nj) * const
       enddo
      enddo
!
      return
      end
!
      subroutine cdmulti(dim1,dim0,cin1,din2,const,cout)
      implicit none
!
! complex matrix multiplication subroutine
! cout = cin1 * din2 * const
!
      integer ni,nj,nk,dim1,dim0
      double precision din2(dim0,dim0)
      complex*16 const
      complex*16 cin1(dim0,dim0),cout(dim0,dim0)
!
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in cdmulti'
       call flush(8)
       stop
      endif
!
      do ni=1,dim1
       do nj=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
        do nk=1,dim1
         cout(ni,nj) = cout(ni,nj) + cin1(ni,nk) * din2(nk,nj)
        enddo
        cout(ni,nj) = cout(ni,nj) * const
       enddo
      enddo
!
      return
      end
!
      subroutine cdmultip(dim1,dim0,cin1,din2,const1,const2,cout)
      implicit none
!
! complex matrix multiplication subroutine
! cout = cin1 * din2 * const1 + const2*cout
!
      integer ni,nj,nk,dim1,dim0
      real*8, intent(in) :: din2(dim0,dim0)
      complex*16,intent(in) :: const1,const2
      complex*16,intent(in) :: cin1(dim0,dim0)
      complex*16,intent(inout) :: cout(dim0,dim0)
      complex*16 :: ctmp(dim0,dim0)
!
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in cdmultip'
       call flush(8)
       stop
      endif
!
      ctmp(1:dim0,1:dim0) = cout(1:dim0,1:dim0)
      do ni=1,dim1
       do nj=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
        do nk=1,dim1
         cout(ni,nj) = cout(ni,nj) + cin1(ni,nk) * din2(nk,nj)
        enddo
        cout(ni,nj) = cout(ni,nj) * const1 
       enddo
      enddo
!
      do ni = 1, dim1
       do nj = 1, dim1
         cout(ni, nj) = cout(ni, nj) + ctmp(ni,nj) * const2
       enddo
      enddo
!   
      return
      end
!
      subroutine cdmulti1(dim1,dim0,cin1,dim2,din2,const,dim3,cout)
      implicit none
!
! complex matrix multiplication subroutine
! cout = cin1 * din2 * const
!
      integer ni,nj,nk,dim1,dim0, dim2, dim3
      double precision din2(dim2,dim2)
      complex*16 const
      complex*16 cin1(dim0,dim0),cout(dim3,dim3)
!
      if (dim1 .gt. dim0 .or. dim1 .gt. dim2 .or. dim1 .gt. dim3) then
       write(8,*)'ERROR in cdmulti'
       call flush(8)
       stop
      endif
!
      do ni=1,dim1
       do nj=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
        do nk=1,dim1
         cout(ni,nj) = cout(ni,nj) + cin1(ni,nk) * din2(nk,nj)
        enddo
        cout(ni,nj) = cout(ni,nj) * const
       enddo
      enddo
!
      return
      end
!
      subroutine cmulti1(dim1,dim0,cin1,dim2,cin2,const,dim3,cout)
      implicit none
!
! complex matrix multiplication subroutine
! cout = cin1 * cin2 * const
!
      integer ni,nj,nk,dim1,dim0,dim2,dim3
      complex*16 const,tmp
      complex*16 cin1(dim0,dim0),cin2(dim2,dim2),cout(dim3,dim3)
!
      if (dim1 .gt. dim0 .or. dim1 .gt. dim2 .or. dim1 .gt. dim3) then
       write(8,*)'ERROR in cmulti1'
       call flush(8)
       stop
      endif
!
      do nj=1,dim1
       do ni=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
       enddo
       do nk=1,dim1
        tmp = cin2(nk,nj) * const
        do ni=1,dim1
         cout(ni,nj) = cout(ni,nj) + tmp * cin1(ni,nk)
        enddo
       enddo
      enddo
!
      return
      end
!
      subroutine cmulmvm(dim1,dim0,cmat1,cmat2,cvec,const,cout)
      implicit none
!
! const * cmat1 * cvec * cmat2 (matrix should be complex square)
! 
      integer dim1,dim0,ni,nj,nk
      complex*16 cmat1(dim0,dim0),cmat2(dim0,dim0),cvec(dim0)
      complex*16 const,cout(dim0,dim0),tmp1
!
      if (dim1 .gt. dim0) then
       write(8,*)'ERROR in cmulmvm'
       call flush(8)
       stop
      endif
!
      do nj=1,dim1
       do ni=1,dim1
        cout(ni,nj) = dcmplx(0.d0,0.d0)
       enddo
       do nk=1,dim1
        tmp1 = cmat2(nk,nj) * cvec(nk) * const
        do ni=1,dim1
         cout(ni,nj) = cout(ni,nj) + cmat1(ni,nk) * tmp1
        enddo
       enddo
      enddo
!
      return
      end

