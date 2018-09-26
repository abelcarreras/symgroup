c-------------------------------------------------------
	subroutine vecnormal(q,m,v,pla_ctrl)
c
c	Calcula el vector normal a una estructura plana
c       pla_ctrl controla que q sigui plana
c
	integer		m,n
	double precision	q(m,3),v(3),q1(3),q2(3)
      double precision  sum
	logical		pla_ctrl
c
	pla_ctrl=.false.
c
	do k=1,m
	 if((q(k,1).gt.1e-7).or.(q(k,2).gt.1e-7).or.(q(k,3).gt.1e-7))then
	  q1(:)=q(k,:)
	  kk=k
	  goto 100
	 end if
    	end do
 100	continue
	do k=k+1,m
         if((q(k,1).gt.1e-7).or.(q(k,2).gt.1e-7).or.
     .      (q(k,3).gt.1e-7))then
          q2(:)=q(k,:)
          goto 200
         end if
	end do
 200	continue
c
	v(1)=q1(2)*q2(3)-q1(3)*q2(2)
        v(2)=q1(3)*q2(1)-q1(1)*q2(3)
        v(3)=q1(1)*q2(2)-q1(2)*q2(1)
c
	n=0
	do k=1,m
	 sum=0.0d0
	 do j=1,3
	  sum=sum+v(j)*q(k,j)
	 end do
	 if(abs(sum).lt.1e-7) n=n+1
	end do
	if (n.eq.m) pla_ctrl=.true.
c
        return
        end

c--------------------------------------------
	subroutine maty(a,b,n,m,y)
c
c	Fa el producte de dues matrius: a(m,3) i b(n,m,3)
c  	El resultat es c(3,3). c=ab+ba
c	Expresem la matriu simetrica com a vector y(6) (triangle superior de c)
c
	integer			m
	double precision 	a(m,3),b(n,m,3),y(6)
c
	y(:)=0.0d0
	ii=0
	do j=1,3
	do i=1,j
	  ii=ii+1
	  do k=1,m
	  do l=1,n
           y(ii)=y(ii)+a(k,i)*b(l,k,j)+a(k,j)*b(l,k,i)
	  end do
	  end do
	end do
	end do
c
	RETURN
	END

c-------------------------------------------------------
        subroutine matyy(a,b,n,m,y)

c Fa el producte de dues matrius: a(m,3) i b(n,m,3)
c El resultat es c(3,3). c=ab+ba
c Expresem la matriu simetrica com a matriu y(3,3)

      integer                 m,n
      double precision        a(m,3),b(n,m,3),y(3,3)

      y(:,:) = 0.0d0

      DO l = 1,n
      DO j = 1,3
      DO i = 1,3
      DO k = 1,m
        y(i,j) = y(i,j) + a(k,i)*b(l,k,j)+a(k,j)*b(l,k,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

	RETURN
      END

	subroutine norma(a,b,n,m,c)

c ---------------------------------------------------------
c Computes the scalar products sum between a and b vectors

	integer	       m,n
	double precision   a(m,3),b(n,m,3),c
c
	c = 0.0d0
	DO l = 1,n
	DO i = 1,m
	DO j = 1,3
	  c = c + a(i,j)*b(l,i,j)
	ENDDO
	ENDDO
	ENDDO

	RETURN
	END
