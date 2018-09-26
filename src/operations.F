c-------------------------------------------------------
c
      subroutine rotacions(y,b,ms,msf,x)

c Calcula el valor de la mesura per un eix de rotacio cn.

	double precision	y(6),b(3),ms,msf,ms1,x(3),x1(3),modul
	double precision  eigvec(3,3),eigval(3),y0(6),vec(3)

      common  /int2/nconv

	iout  = 6
	small = 1.0d-7
	y0(:) = y(:)
	modul = dsqrt(b(1)**2+b(2)**2+b(3)**2)

      IF((abs(y(1)).lt.small).and.(abs(y(2)).lt.small).and.
     & (abs(y(3)).lt.small).and.(abs(y(4)).lt.small).and.
     & (abs(y(5)).lt.small).and.(abs(y(6)).lt.small)) then
         msf=ms+2.0d0*modul
         x(:)=-b(:)/modul
      ELSE
	  CALL jacobi(y,3,eigval,eigvec,iout)
        IF((abs(b(1)).lt.small).and.(abs(b(2)).lt.small).and.
     &     (abs(b(3)).lt.small)) then
	    x(:)=eigvec(:,3)
          msf=ms-eigval(3)
        ELSE
c Inici de les Iteracions
	    msf=-1.d6
          DO i = 1,3
            vec(:) = eigvec(:,i)
            CALL iterNewton(ms,b,eigval(i),vec,y0,ms1,x1)
            IF (ms1.gt.msf) THEN
              msf=ms1
              x(:)=x1(:)
            ENDIF
 	    ENDDO

c	    i=3
c         if(nconv.eq.1) i=1
c         do while (i.le.3)
c           vec(:)=eigvec(:,i)
c           do ll=1,2
c             call iterNewton(ms,b,eigval(i),vec,y0,ms1,x1)
c             if(ms1.gt.msf) then
c               msf=ms1
c               x(:)=x1(:)
c             end if
c             vec(:)=-vec(:)
c           end do
c           i=i+1
c         end do

      ENDIF
      ENDIF
c
	RETURN
	END

c-------------------------------------------------------
c
        subroutine reflexio(q,m,v,p)

c Calcula les coordenades de l'estructura reflexada
c q es l'estructura inicial i p la final. v es el vector normal al pla

	integer		m
	double precision	v(3),q(1,m,3),p(m,3)

	p = 0.0d0

	do i = 1,3
	do j = 1,3
	do k = 1,m
	   p(k,i) = p(k,i)-2.0d0*v(i)*v(j)*q(1,k,j)
	end do
	end do
	end do

	p(:,:) = q(1,:,:) + p(:,:)

	RETURN
	END
c-------------------------------------------------------
c
      subroutine rotacio(q,n,m,v,p,ctheta,stheta,ttheta)

c Calcula les coordenades de l'estructura rotada
c q es l'estructura inicial i p la final. v es el vector de l'eix de rotacio
c theta es l'angle de rotacio

      integer           m,n
      double precision  v(3),q(n-1,m,3),p(m,3)
	double precision	ctheta(n-1),stheta(n-1),ttheta(n-1)

      p = 0.0d0

      DO i = 1,3
      DO j = 1,3
      DO k = 1,m
	DO l = 1,n-1
        p(k,i)=p(k,i)+ttheta(l)*v(i)*v(j)*q(l,k,j)
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	DO k = 1,m
	DO l = 1,n-1
	  p(k,1)=p(k,1)+stheta(l)*(v(2)*q(l,k,3)-v(3)*q(l,k,2))
        p(k,2)=p(k,2)+stheta(l)*(v(3)*q(l,k,1)-v(1)*q(l,k,3))
        p(k,3)=p(k,3)+stheta(l)*(v(1)*q(l,k,2)-v(2)*q(l,k,1))
	ENDDO
	ENDDO

      DO l = 1,n-1
        p(:,:) = ctheta(l)*q(l,:,:) + p(:,:)
	ENDDO

      RETURN
	END

c-------------------------------------------------------
      subroutine iterNewton(ms,b0,l0,vec,y0,ms1,xf)

c Iteracions per obtenir lambda i la orientacio per la rotacio
c l0 i x0: lambda i vector orientacio inicials, y: matriu en forma y(6),
c b: vector independent. lf i xf: lambda i vector orientacio finals
c cas en que b*x es 0.0

      integer                 ipiv(4),info
      double precision        modul,ms,ms1,lf,l0
      double precision        y0(6),b0(3),vec(3)
      double precision        F(4),J(4,4),x(4),xf(3),b(4,1)

c Definicio del vector inicial x(4)
      x(1) = l0
      DO i = 1,3
        x(i+1)=vec(i)
	ENDDO

c ITERACIONS
      ms1=-1.d6
      ii=0
      jj=0
 100  if (jj.eq.2) goto 2000

c Construccio vector inicial F(4)
      F(1) = y0(1)*x(2)+y0(2)*x(3)+y0(4)*x(4)-x(1)*x(2)+b0(1)
      F(2) = y0(2)*x(2)+y0(3)*x(3)+y0(5)*x(4)-x(1)*x(3)+b0(2)
      F(3) = y0(4)*x(2)+y0(5)*x(3)+y0(6)*x(4)-x(1)*x(4)+b0(3)
      F(4) = x(2)*x(2)+x(3)*x(3)+x(4)*x(4)-1.0d0

c Construccio matriu inicial J(4,4)
      J(1,1) = -x(2)
      J(2,1) = -x(3)
      J(3,1) = -x(4)
      J(4,1) = 0.0d0
      J(1,2) = y0(1)-x(1)
      J(2,2) = y0(2)
      J(3,2) = y0(4)
      J(4,2) = 2.0d0*x(2)
      J(1,3) = y0(2)
      J(2,3) = y0(3)-x(1)
      J(3,3) = y0(5)
      J(4,3) = 2.0d0*x(3)
      J(1,4) = y0(4)
      J(2,4) = y0(5)
      J(3,4) = y0(6)-x(1)
      J(4,4) = 2.0d0*x(4)

      DO i = 1,4
        b(i,1) = -F(i)
      ENDDO
c     call dgesv_(4,1,J,4,ipiv,b,4,info)
      CALL DGESV(4,1,J,4,ipiv,b,4,info)
      DO i = 1,4
        x(i) = x(i)+b(i,1)
      ENDDO
      modul = 0.0d0
      DO i = 2,4
        modul=modul+x(i)*x(i)
      ENDDO

      IF ((abs(b(1,1)).gt.1.d-6).or.(abs(b(2,1)).gt.1.d-6)
     &.or.(abs(b(3,1)).gt.1.d-6).or.(abs(b(4,1)).gt.1.d-6)) THEN
        ii=ii+1
        IF (ii.eq.100) THEN
          DO i = 1,4
            x(i) = x(i)-0.5d0*b(i,1)
          ENDDO
          jj=jj+1
          ii=0
        ENDIF
        GOTO 100
      ENDIF

      lf=x(1)
      do i=1,3
        xf(i)=x(i+1)
      end do
      ms1=b0(1)*xf(1)+b0(2)*xf(2)+b0(3)*xf(3)
      ms1=ms-lf-ms1
c     ms1=2.0d0*ms-lf-ms1  ! To be used when only computing half of the operations.

c2000 if (jj.eq.2) write(3,*)'Newton not converged after 100 steps'
 2000	CONTINUE
c
	RETURN
	END
