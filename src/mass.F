c------------------------------------------------
        subroutine cmasses(a,b,m,ctrl)
c------------------------------------------------
c  retorna l'estructura a centrada al centre de masses
c
	integer           m
      double precision  a(m,3),b(m,3)
      double precision  cm(3)
      logical           ctrl
c
      common /center/cm
c
      if (.NOT.ctrl) then
        cm(:)=0.0d0
        do j=1,3
         do i=1,m
          cm(j)=cm(j)+a(i,j)
         end do
        end do
        cm(:)=cm(:)/dfloat(m)
      end if
c
c       Coloquem la molecula ideal amb el cm a l'origen
c
        do j=1,3
         do i=1,m
	    b(i,j)=a(i,j)-cm(j)
         end do
        end do
c
	return
	end
