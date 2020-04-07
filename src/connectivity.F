	subroutine connectivity(m,eti,q,vcon)
c
c Calcula la matriu de connectivitat a partir dels radis atomics (radius.f)
c
	integer		    znum(m)
	double precision    r(85),q(m,3),d,vcon(m*(m+1)/2),thresh
	character	    eti(m)*8,zchar(m)*2,lab*2
	logical              :: ctrl
        character, parameter :: sim(85)*2 = (/
     &        "H ",                                     "HE",
     &        "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE",
     &        "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR",
     &        "K ", "CA", "SC", "TI", "V ", "CR", "MN", "FE",
     &        "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE",
     &        "BR", "KR", "RB", "SR", "Y ", "ZR", "NB", "MO",
     &        "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
     &        "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE",
     &        "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
     &        "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W ",
     &        "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB",
     &        "BI", "PO", "AT" /)
        character, parameter :: fit(85)*2=(/
     &        "1 ","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ","10",
     &        "11","12","13","14","15","16","17","18","19","20",
     &        "21","22","23","24","25","26","27","28","29","30",
     &        "31","32","33","34","35","36","37","38","39","40",
     &        "41","42","43","44","45","46","47","48","49","50",
     &        "51","52","53","54","55","56","57","58","59","60",
     &        "61","62","63","64","65","66","67","68","69","60",
     &        "71","72","73","74","75","76","77","78","79","80",
     &        "81","82","83","84","85" /)
c
      integer nuat(30)
      common /real/r
      common /sumrad/thresh
C--MQL ///
c
C--MQL \\\
       write(6,"(/'Connectivity radius factor: ',f7.4)") thresh
C--MQL ///
c
c Pas de l'etiqueta atomica a Z atomic
c
	do k=1,m
	 read(eti(k)(1:2),"(a2)")lab(1:2)
	 do i=1,2
          l=ichar(lab(i:i))
          if (l.ge.97.and.l.le.122) then
           lab(i:i)=char(l-32)
          endif
	 end do
	 if((ichar(lab(2:2)).lt.65).or.
     &  (ichar(lab(2:2)).gt.90)) lab(2:2)=' '
	 ctrl=.false.
	 ATOM: do i=1,85
	  if(trim(lab).eq.trim(sim(i)))then
	   zchar(k)=fit(i)
	   ctrl=.true.
	   exit ATOM
	  end if
	 end do ATOM
	 if(.not.ctrl) stop 'ERROR. Invalid atomic label'
	 read(zchar(k),'(i2)')znum(k)
	end do
	
C--MQL \\\
        write(3,"('Specie     Radius')")
        ka=0
	do i=1,m
          ctrl=.true.
          CHK: do k=1,ka
            if (nuat(k).eq.znum(i)) then
              ctrl=.false.
              exit CHK
            endif
          enddo CHK
          if (ctrl) then
            ka=ka+1
            nuat(ka)=znum(i)
            write(3,"(i4,6x,f7.4)") nuat(ka), r(nuat(ka))
          endif
        enddo
        write(3,"('Species         R+R       Cutoff')")
        do i=1,ka
        do j=i,ka
          d=r(nuat(i))+r(nuat(j))
          write(3,"(2i4,2(5x,f7.4))") nuat(i), nuat(j), d, d*thresh
        enddo
        enddo
C--MQL ///
c
c Calcul de la matriu de connectivitat
c
	vcon(:)=0.0d0
	kk=0
	do i=1,m
	 do j=1,i
	  kk=kk+1
	  if(i.ne.j) then
	   d=(q(i,1)-q(j,1))**2+(q(i,2)-q(j,2))**2+
     .       (q(i,3)-q(j,3))**2
	   d=dsqrt(d)
	   if(d.lt.thresh*(r(znum(i))+r(znum(j)))) vcon(kk)=1.0d0
	  end if
	 end do
	end do
C--MQL \\\
        write(3,"(/'Atom connectivity')")
        do im=1,m
          write(3,"(i3,'   ')",advance='no') im
          kk=0
          do i=1,m
          do j=1,i
            kk=kk+1
            if(vcon(kk).gt.0.d0.and.i.eq.im)
     &        write(3,"(i4)",advance='no') j
          enddo
          enddo
          write(3,*)
        enddo
        write(3,"(/'Connectivity list')")
        do im=1,m
          kk=0
          do i=1,m
          do j=1,i
            kk=kk+1
            if(vcon(kk).gt.0.d0) then
              if (j.eq.im) write(3,"(2i4)") im, i
            end if
          enddo
          enddo
        enddo
C--MQL ///
c
	return
	end 
