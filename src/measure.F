      subroutine mesura(oeix,i_oper,noper,q,p,nll,s,vec,
     .                  ctheta,stheta,ttheta)

c Computes Symmetry Group Measure for each permutation

	integer	      nll,i_oper,dimr,neix,iout,noper,oeix
      double precision  p(noper,nll,3),pp(noper,nll,3)
	double precision  q(nll,3),eigvec(3,3),vecnorm(3)
	double precision  y(6),eigval(3),b(3),vmax(3),vec(3),yy(3,3)
      double precision  ms,ms1,s,sum
      double precision  stheta(noper),ctheta(noper),ttheta(noper)

      common  /int1/dimr
      common  /real2/vecnorm,vmax

      iout=6

      IF (i_oper.eq.1) THEN
	  CALL norma(q,p,noper,nll,ms)
        ms = -ms
      ELSEIF(i_oper.eq.2) THEN
        CALL norma(q,p,noper,nll,ms)
        CALL maty(q,p,noper,nll,y)
        CALL jacobi(y,3,eigval,eigvec,iout)
        IF (dimr.eq.3) THEN
          ms = ms-eigval(3)
          vec(:) = eigvec(:,3)
        ENDIF
c 2D analysis
        IF (dimr.eq.2) THEN
          neix=0
          DO ii = 3,1,-1
            sum = 0.0d0
            DO j = 1,3
              sum = sum+vecnorm(j)*eigvec(j,ii)
            ENDDO
            IF (abs(sum).lt.1e-7) THEN
              ms = ms-eigval(ii)
              vec(:) = eigvec(:,ii)
              GOTO 203
            ELSE
              neix=neix+1
            ENDIF
          ENDDO
          IF (neix.eq.3) ms=0.0d0
 203    ENDIF
c End 2D analysis

c Rotations
	ELSEIF ((i_oper.eq.3).or.(i_oper.eq.4)) THEN
	  DO l = 1,noper
          pp(l,:,:) = ctheta(l)*p(l,:,:)
	  ENDDO
c       IF ((oeix/2)*2.eq.oeix) pp(noper,:,:)=0.5d0*pp(noper,:,:)
        CALL norma(q,pp,noper,nll,ms)
        DO l = 1,noper
          pp(l,:,:) = ttheta(l)*p(l,:,:)
        ENDDO
c       IF ((oeix/2)*2.eq.oeix) pp(noper,:,:)=0.5d0*pp(noper,:,:)
        CALL maty(q,pp,noper,nll,y)
        y(:) = -0.5d0*y(:)
        DO l = 1,noper
          pp(l,:,:) = stheta(l)*p(l,:,:)
        ENDDO
c       IF((oeix/2)*2.eq.oeix) pp(noper,:,:)=0.5d0*pp(noper,:,:)
        CALL gradBk(q,pp,noper,nll,b)
        b(:)=-0.5d0*b(:)
        IF (i_oper.eq.3) THEN
          CALL rotacions(y,b,ms,ms1,vec)
        ELSE
          CONTINUE
c         CALL mat_T(q,p,nll,vmax,ms1)
c         ms1 = ctheta*ms+ttheta*ms1-2.0d0*
c    .    (b(1)*vmax(1)+b(2)*vmax(2)+b(3)*vmax(3))
c         vec(:) = vmax(:)
        ENDIF
        ms = ms1
	ENDIF

 	s = ms

      RETURN
      END

c-------------------------------------------------------
      subroutine mat_T(q,p,m,eix,y)

c Calcula la contribucio a la mesura de rotacio corresponent a la matriu T

      integer                 m
      double precision        q(m,3),p(m,3),eix(3),y,yp,yq

	y = 0.0d0
	do k = 1,m
	  yp=0.0d0
	  yq=0.0d0
	  do j = 1,3
	    yp = yp+p(k,j)*eix(j)
          yq = yq+q(k,j)*eix(j)
	  end do
	  y = y+yp*yq
	end do

	RETURN
	END

c-------------------------------------------------------
	subroutine gradBk(a,b,n,m,c)

c Obtenim el terme independent c de l'equacio de la rotacio
c m es el numero d'atoms i n el numero d'operacions

	integer		m
	double precision 	a(m,3),b(n,m,3),c(3)

	c(:)=0.0d0

	DO l = 1,n
	DO k = 1,m
	  c(1) = c(1) + a(k,3)*b(l,k,2)-a(k,2)*b(l,k,3)
        c(2) = c(2) + a(k,1)*b(l,k,3)-a(k,3)*b(l,k,1)
        c(3) = c(3) - a(k,1)*b(l,k,2)+a(k,2)*b(l,k,1)
	ENDDO
	ENDDO

	RETURN
	END

c-------------------------------------------------------
	subroutine mat_A(q,p,m,v,y)

c Calcula la contribucio a la mesura de rotacio corresponent a la matriu A

	integer		m
 	double precision	q(m,3),p(m,3),v(3),A(3,3),y

	do i=1,3
	 A(i,i)=0.0d0
	end do
	A(1,2) = -v(3)
	A(1,3) = v(2)
	A(2,3) = -v(1)
      A(2,1) = v(3)
      A(3,1) = -v(2)
      A(3,2) = v(1)

	y = 0.0d0
	DO k = 1,m
	DO i = 1,3
	DO j = 1,3
	  y = y + q(k,i)*A(i,j)*p(k,j)
	ENDDO
	ENDDO
	ENDDO

	RETURN
	END
