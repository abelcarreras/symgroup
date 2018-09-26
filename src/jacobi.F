c-------------------------------------------------------
        subroutine jacobi(a,n,e,v,iout)
c
c-------------------------------------------------------------------------------
c
c  diagonalize a matrix using the jacobi method
c
c  input:
c          a      matrix to diagonalize. dimension n*(n+1)/2
c          n      dimension of the matrix to diagonalize
c       iout      channel to write the results
c
c  output:
c          e      eigenvalues. dimension n
c          v      eigenvectors. dimension n*n
c
c-------------------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      data toler/1.0d-20/
      dimension a(*), v(n,*), e(*)
      index(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
C--MQL --
! intent only perque es crida amb un 6, no amb variable
!      iout =6
C--MQL --/
      do 10 j=1,n
         do 20 i=1,n
            v(i,j)=0.0d0
   20    continue
         v(j,j)=1.0d0
         jj=j*(j+1)/2
         e(j)=a(jj)
   10 continue
      do 50 iter=1,100
         aijmax=0.0d0
         do 30 i=2,n
            j1=i-1
            do 30 j=1,j1
               ij=i*(i-1)/2+j
               aii=e(i)
               ajj=e(j)
               aij=a(ij)
               aij2=aij**2
               if (aij2 .gt. toler) then
                  if(aij2 .gt. aijmax) aijmax=aij2
                  ded=aii-ajj
                  if (ded .lt. 0.0d0) then
                     dos=-2.0d0
                     ded=-ded
                  else
                     dos=2.0d0
                  endif
                  den=ded+dsqrt(ded**2+4.0d0*aij2)
                  tan=dos*aij/den
                  cos=1.0d0/dsqrt(1.0d0+tan**2)
                  sin=cos*tan
                  do 40 k=1,n
                     x1=cos*v(k,j)-sin*v(k,i)
                     v(k,i)=sin*v(k,j)+cos*v(k,i)
                     v(k,j)=x1
                     if (k.ne.j .and. k.ne.i) then
                        jk=index(j,k)
                        ik=index(i,k)
                        x1=cos*a(jk)-sin*a(ik)
                        a(ik)=sin*a(jk)+cos*a(ik)
                        a(jk)=x1
                     endif
   40             continue
                  e(i)=cos*cos*aii+sin*sin*ajj+2.0d0*sin*cos*aij
                  e(j)=cos*cos*ajj+sin*sin*aii-2.0d0*sin*cos*aij
                  a(ij)=0.0d0
               endif
   30    continue
         if (aijmax .le. toler) then
c
c     order the eigenpairs
c
             do 60 i=1,n
             do 60 j=1,n
                if (e(j) .lt. e(i)) then
                   x1=e(i)
                   e(i)=e(j)
                   e(j)=x1
                  do 70 k=1,n
                     x1=v(k,i)
                     v(k,i)=v(k,j)
                     v(k,j)=x1
   70             continue
               endif
   60       continue
            return
         endif
   50 continue
      write (iout,'('' warning: jacobi diag. did not converge'')')
      end
