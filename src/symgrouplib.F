
       subroutine symgroup(coordinates, Nmulti_inp,
     &symbols, pm, oper, oeix, m, s_max, NearStruc,
     &v_max, optper, refAxis, sm, vm)
c
c
c   Aquest programa mesura la distancia entre un estructura i l'estructura
c   despres d'aplicar-l'hi una operacio de simetria
c
c-----------------------------------------------------------
c     David Casanova    29/NOV/06     v0c_6 mes flexibilitat en l'input
c	  David Casanova    14/DES/06     opcions input (lectura d'inputs Miquel)
c     David Casanova    17/DES/06     v1.0a friendly version
c     David Casanova    23/NOV/09     adapt to symgroup program
c     David Casanova    23/DES/09
c     David Casanova    26/MAY/11	  remove "poten" keyword
c     Abel Carreras     19/SEP/18	  created library
c-----------------------------------------------------------
c VARIABLES
c
c  iunit unitat d'entrada o sortida de dades
c  m     numero d'atoms
c  nll   numero de lligands
c  pm    posicio del metall
c  coordinates  coordenades dels atoms (inclos el metall)
c  symbols  symbols atomics (o etiquetes dels atoms)
c  oeix  ordre de l'eix de rotacio
c  oper  operacio de simetria: E i r c
c  noper numero d'operacions del grup
c  s     resultat del calcul de cssm
c  dimr  dimensio de la reflexio (2D o 3D)
c  ivecs  posicions dels atoms a partir dels quals es defineix l'orientacio de l'eix (Cn o Sn)
c  nvecs  nombre d'atoms a partir dels quals es defineix l'orientacio de l'eix (Cn o Sn)
c  dvmax  modul al quadrat del l'orientacio (fixada) de l'eix (Cn o Sn)
c  Nmulti numero de mesures per cada estructura (maxim valor 10)
c  ptotal conjunt d'estructures que genera el grup
c-----------------------------------------------------------
	integer	     :: m,pm,nll,es,dimr,i_oper,nvecs,noper,
     .                    oeix,Nmulti,iunit
	double precision :: pi,s,s_max,A
	double precision :: vecnorm(3),vmax(3),v(3),dvmax
      double precision :: sm(20),vm(20,3)
      character        :: oper*1,name*40
      character        :: input*80,title*80,line*80
	logical	     :: keyctrl,pla_ctrl,ctrl_nvecs,group
	logical	     :: label,multi,hconv,connect,permu,
     .            atfix,rfix,dim2D,csd,nolabel,centerfix
c
        character, parameter :: fit(90)*2=(/
     &        "1 ","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ","10",
     &        "11","12","13","14","15","16","17","18","19","20",
     &        "21","22","23","24","25","26","27","28","29","30",
     &        "31","32","33","34","35","36","37","38","39","40",
     &        "41","42","43","44","45","46","47","48","49","50",
     &        "51","52","53","54","55","56","57","58","59","60",
     &        "61","62","63","64","65","66","67","68","69","60",
     &        "71","72","73","74","75","76","77","78","79","80",
     &        "81","82","83","84","85","86","87","88","89","90" /)
c
	logical, allocatable ::	     ctrl(:)
	integer, allocatable ::      ip(:),ipf(:),ivecs(:),ivecs0(:)
      double precision,allocatable :: q(:,:),q0(:,:),qf(:,:)
      double precision,allocatable :: p(:,:,:),ptotal(:,:,:)  ! Group
	double precision,allocatable :: convep(:,:),convap(:),vcon(:)
 	double precision,allocatable :: theta(:),ctheta(:),stheta(:),ttheta(:)
      character*8,allocatable ::       etiq0(:),etiq(:)
C--MQL
        logical :: mute
        integer, allocatable :: kcon(:), jcon(:,:)
c        character :: yn*1    ! DCC: this version always runs
        integer :: npmt(4), khg
        real*8  :: apmt(4)
      double precision  cm(3)

      real*8, dimension(m,3) :: coordinates, NearStruc
      integer :: Nmulti_inp
      character, dimension(m) :: symbols*2
      real*8, dimension(3) :: v_max
      integer, dimension(m) :: optper
      real*8, dimension(4,3) :: refAxis
C--MQL ///
c
	common  /int1/dimr
	common  /int2/nconv
	common  /real2/vecnorm,vmax
	common  /real3/A
      common /intlog/Nmulti
      common /center/cm

C--MQL --
        common /logical/label,multi,hconv,permu,atfix,
     & rfix,dim2D,csd,connect,nolabel,mute,centerfix
C--MQL --/

	title=' '
	pi = 4.0d0*atan(1.0d0)
      dimr  = 3
      nvecs = 0
      k=1
	noper  = 1
	Nmulti = 1
	iunit  = 1
c
	keyctrl = .true.
	label   = .false.
	multi   = .false.
	hconv   = .false.
	connect = .false.
	permu   = .false.
	atfix   = .false.
	rfix    = .false.
      centerfix = .false.
	dim2D   = .false.
	csd     = .false.
	nolabel=.false.
C--MQL
        mute=.true.
C--MQL ///
c
c      Updated things
       if ((pm).eq.0) then
         nll = m
       else
         nll = m - 1
       endif

      title=' '
      name=' '

      Nmulti = nmulti_inp
      if (Nmulti.gt.0) multi = .true.
      label = .true.
      input = 'test'

c Read input file name (first argument)
c      CALL getarg(1,input)
c      IF (trim(input).eq.'') THEN
c        write(0,*) 'No input file provided'
c        write(0,*) 'usage: ./symgroup input_file'
c        stop
c      ENDIF

	IF (trim(input).eq.'+') THEN
c       open(unit=33,file="operations.dat",status='unknown')
	  write(6,*)' '
        write(6,*)'* Identity Symmetry'
        write(6,*)'  E'
        write(6,*)' '
	  write(6,*)'* Inversion Symmetry Group'
	  write(6,*)'  Ci'
        write(6,*)' '
        write(6,*)'* Reflection Symmetry Group'
        write(6,*)'  Cs'
        write(6,*)' '
	  write(6,*)'* Rotational Symmetry Group'
	  write(6,*)'  C n       n: rotation order'
        write(6,*)' '
        write(6,*)'* Rotation-Reflection Symmetry Group'
        write(6,*)'  S n       n: rotation-reflection order'
        write(6,*)' '
        write(6,*)' '
c       close(33)
	  goto 9999
	ENDIF
c
       IF (trim(input).eq.'+keywords') THEN
c      open(unit=33,file="keywords.dat",status='unknown')
        write(6,*)' '
        write(6,"('%LABEL     Comparing vertices with the same label')")
        write(6,"('%MULTI m   m best measures (file .zout2)')")
        write(6,"(
     &'%PERMU     Measure with the permutation in .zper file')")
c       write(6,"('%ATFIX     Fix the symmetry element using m atoms')")
c       write(6,"('%RFIX      Fix the symmetry element orientation')")
        write(6,"(
     &'%CENTERFIX Fix the center of the symmetry element')")
        write(6,"('%2D        2D measure')")
        write(6,"('%CSD       Structures from file .cor (CSD)')")
        write(6,"('%CONNECT   Comparing vertices by the connectivity')")
        write(6,"('%NOLABEL   Structures without labels')")
c       write(6,"('%HCONV     High convergence')")
        write(6,*)' '
c       close(33)
        GOTO 9999
	ENDIF
c
c Open I/O files
      k=max(len_trim(input),5)
      if (input(k-4:k).eq.'.zdat') input(k-4:k)='    '
c      open(unit=1,file=trim(input)//".zdat",status='old',iostat=k)
c      open(unit=2,file=trim(input)//".ztab",status='unknown')
c      open(unit=3,file=trim(input)//".zout",status='unknown')

c Write Header to .zout file
c      write(3,*)' '
c      write(3,*)'  *********************************************'
c      write(3,*)'  *    Welcome to the SYMGROUP Program v0     *'
c      write(3,*)'  *                                           *'
c      write(3,*)'  *   Continuous Symmetry Group calculation   *'
c      write(3,*)'  *********************************************'
c      write(3,*)' '
c      write(3,*)' by David Casanova'
c      write(3,*)' (c) 2009 Electronic Structure Group'
c      write(3,*)'            Universitat de Barcelona'
c      write(3,*)' '

c Write Header to .ztab file
c      write(2,*)' '
c      write(2,*)'  *********************************************'
c      write(2,*)'  *    Welcome to the SYMGROUP Program v0     *'
c      write(2,*)'  *                                           *'
c      write(2,*)'  *   Continuous Symmetry Group calculation   *'
c      write(2,*)'  *********************************************'
c      write(2,*)' '
c      write(2,*)' by David Casanova'
c      write(2,*)' (c) 2009 Electronic Structure Group'
c      write(2,*)'            Universitat de Barcelona'
c      write(2,*)' '

c      if (k.eq.0) then
c       write(0,*)' '
c       write(0,"('Reading data from file ',a)") trim(input)//'.zdat...'
c       write(0,*)' '
c      else
c       write(0,"('File not found')")
c      end if

c	call nextline(iunit)
c	read(iunit,"(a)") line
c	if((line(1:1).eq.'$').and.(len_trim(title).eq.0)) then
c	 title(1:79)=line(2:80)
c	 call nextline(iunit)
c	 read(iunit,"(a)") line
c	end if
c	do while (keyctrl)
c	 if((line(1:1).eq.'%').and.(len_trim(line).gt.1)) then
c	  call options(line)
c	  call nextline(iunit)
c	  read(iunit,"(a)") line
c	 else
c	  keyctrl=.false.
c	 end if
c	end do

c At the moment, fixing orientation (atfix or rfix) not ready.
	IF (atfix .OR. rfix) THEN
	  WRITE(0,*)'ERROR. ATFIX and RFIX not available in this version'
c	  WRITE(3,*)'ERROR. ATFIX and RFIX not available in this version'
	  STOP
	ENDIF
c At the moment, HCONV not ready
	IF (hconv) THEN
        WRITE(0,*)'ERROR. HCONV not available in this version'
c        WRITE(3,*)'ERROR. HCONV not available in this version'
        STOP
      ENDIF

c	backspace(iunit)
c	read(iunit,*)nll,pm
c	read(iunit,*)
c	call m_control(nll,pm)
c
	if(connect) allocate(convep(m,m),convap(m),vcon(m*(m+1)/2))
C--MQL
        allocate (kcon(m),jcon(m,m))
C--MQL ///
c       allocate(q(m,3),q0(m,3),p(m,3),qf(m,3)) Group
        allocate(q(m,3),q0(m,3),qf(m,3))
        allocate(etiq0(m),etiq(m))
        allocate(ip(m),ipf(m),ivecs(m),ivecs0(m))
	allocate(ctrl(m))

c Initialize to 0
        do i = 1,m
           ipf(i) = 0
           ip(i) = 0
           ivecs(i) = 0
           ivecs0(i) = 0
        enddo

c Lectura del fitxer de la permutacio
c       if(permu) then
c         open(unit=55,file=trim(input)//".zper",status='unknown')
c         write(0,"('Reading permutation from file ',a)") trim(input)//'.zper'
c         write(0,"('Only operation measure when permutation is fixed')") ! Group
c         read(55,*)(ipf(j),j=1,5)
c         close(55)
c        end if


c Coordinates assignation
      kka = 1
      do kk = 1,nll
        if (kka.eq.pm) kka = kka + 1
	    etiq0(kk) = symbols(kka)
	    q0(kk,:) = coordinates(kka,:)
	    kka = kka + 1
      enddo
      if (pm.ne.0) then
        etiq0(nll+1) = symbols(pm)
	    q0(nll+1,:) = coordinates(pm,:)
      endif

      atfix = .FALSE.
      rfix = .FALSE.

c------------------------
c Lectura de l'operacio
c------------------------
c	call nextline(iunit)
c	read(1,"(a)") line
c      oper = oper_2

        IF((oper.ne.'e').and.(oper.ne.'i').and.(oper.ne.'r')
     &.and.(oper.ne.'c').and.(oper.ne.'s')) THEN
       write(0,"(/,'ERROR. This is not a defined symmetry operation')")
	    stop
	  ENDIF
c	line(1:1)=' '
c	line=adjustl(line)
c	k=0
c      if ((atfix).or.(rfix)) then
c	 do i=1,len_trim(line)
c	  if (line(i:i).eq.'(') k=i
c	  if (line(i:i).eq.')') kk=i
c	  if (line(i:i).eq.',') line(i:i)=' '
c	 end do
c	 IF(k.eq.0) THEN
c        write(0,"(/,'ERROR. Character ( not found')")
c	   stop
c	 ENDIF
	 if(atfix) then
c	  ctrl_nvecs=.true.
c	  do i=k+1,kk-1
c	   if((ctrl_nvecs).and.(line(i:i).ne.' ')) nvecs=nvecs+1
c	   IF (nvecs.gt.m) THEN
c           write(0,"(/,'ERROR. Too many atoms to fix the axis')")
c	     stop
c	   ENDIF
c	   ctrl_nvecs=(line(i:i).eq.' ')
c	  end do
	  IF (nvecs.eq.0) THEN
          write(0,"(/,'ERROR. There is no atoms to fix the axis')")
	    stop
	  ENDIF
	  IF((oper.eq.'r').and.(nvecs.gt.2)) THEN
           write(0,"(/,'ERROR: Bad fixed plane definition')")
	     stop
	  ENDIF
	 end if
c 	 line(k:k)=' '
c 	 line(kk:)=' '
c	 if(oper.eq.'r') then
c	  if(atfix) read(line,*)(ivecs0(j),j=1,nvecs)
c        if(rfix) read(line,*)(vmax(j),j=1,3)
c	 elseif ((oper.eq.'c').or.(oper.eq.'s')) then
c	  if(atfix) read(line,*)oeix,(ivecs0(j),j=1,nvecs)
c	  if(rfix) read(line,*)oeix,(vmax(j),j=1,3)
c	 end if
c	elseif((oper.eq.'c').or.(oper.eq.'s')) then
c	 read(line,'(i2)')oeix
c	 read(line,'(i2)')
c	end if
c      oeix = oeix_2
c      write(6,*) 'oeix',oeix
c Control atfix definition
	IF (atfix) THEN
        DO i = 1,nvecs
	    IF (ivecs0(i).GT.m) THEN
	      WRITE(0,"(/,'ERROR: Bad atom definition in ATFIX')")
	      STOP
	    ENDIF
	  ENDDO
	ENDIF

	if ((dim2D).and.(oper.eq.'r')) then
	  dimr=2
	endif

C--MQL --
      IF ((oper.eq.'c'.or.oper.eq.'s').and.oeix.gt.90) THEN
        write(0,"(/,'ERROR. Invalid axis order value')")
        stop
      END IF
C--MQL --/
c
c	open(unit=2,file=trim(input)//".ztab",status='unknown')
c	open(unit=3,file=trim(input)//".zout",status='unknown')
c	write(3,"('Writing output files ',a,a)") trim(input)/
c    .  /'.ztab and ',trim(input)//'.zout'
c	if (multi) then
c	 open(unit=4,file=trim(input)//".zout2",status='unknown')
c	 write(3,"('Writing multiple values to file ',a)") trim(input)
c     &   //'.zout2'
C--MQL --
c	 write(4,'(a)') trim(title)
C--MQL --/
c       write(4,*)' '
c	 write(4,*)'**************************************'
c	end if
c      write(2,'(a)') trim(title)
c      write(3,'(a)') trim(title)

c Allocate some space
      IF ((oper.eq.'c').or.(oper.eq.'s')) noper = oeix-1
      allocate(p(noper,m,3))
      allocate(ptotal(noper,m,3))
      allocate(theta(noper),ctheta(noper),stheta(noper),ttheta(noper))

c Write some file header
	line='                                                     '
	IF (oper.eq.'e') THEN
	  i_oper = 0
c	  WRITE(2,800)
c	  WRITE(3,800)
	ELSEIF (oper.eq.'i') THEN
	  i_oper = 1
c        WRITE(2,810)
c        WRITE(3,810)
	ELSEIF (oper.eq.'r') THEN
	  i_oper = 2
c        WRITE(2,820)
c        WRITE(3,820)
	  IF ((atfix).or.(rfix)) THEN
	    theta(1)  =  2.0d0*pi
          ctheta(1) =  1.0d0
          stheta(1) =  0.0d0
          ttheta(1) = -2.0d0
	  ENDIF
	ELSEIF ((oper.eq.'c').or.(oper.eq.'s')) THEN
	  i_oper = 3
	  DO l = 1,noper
	    theta(l)  = dfloat(l)*2.0d0*pi/dfloat(oeix)
          ctheta(l) = cos(theta(l))
          stheta(l) = sin(theta(l))
	  ENDDO
	  IF (oper.eq.'c') THEN
c	    IF (dim2D) THEN
c	      write(2,830) trim(fit(oeix))
c	      write(3,830) trim(fit(oeix))
c	    ELSE
c	      write(2,840) trim(fit(oeix))
c	      write(3,840) trim(fit(oeix))
c	    ENDIF
	    DO l = 1,noper
            ttheta(l) = 1.0d0-ctheta(l)
	    ENDDO
        ELSE
c          if (dim2D) then
c            write(2,850) trim(fit(oeix))
c            write(3,850) trim(fit(oeix))
c          else
c            write(2,860) trim(fit(oeix))
c            write(3,860) trim(fit(oeix))
c          end if
          DO l = 1,noper
            ttheta(l) = (-1.0d0)**l-ctheta(l)
          ENDDO
	  ENDIF
	ENDIF
	if (dim2D) line(1:1)='2'
c
	if((atfix).or.(rfix)) then
	 i_oper=4
c	 write(2,*)'(Fixed Orientation)'
c	 write(3,*)'(Fixed Orientation)'
	end if
c
c	write(2,*)' '
c	write(2,870)
c	write(3,"(/,'**************************************')")

 800	format('Identity Measures')
 810  format('Inversion Measures')
 820  format('Reflection Measures')
 830	format('2D Rotation Measures C',1x,a,1x,'Group')
 840	format('3D Rotation Measures C',1x,a,1x,'Group')
 850	format('2D Rotation-Reflection Measures S',1x,a,1x,'Group')
 860	format('3D Rotation-Reflection Measures S',1x,a,1x,'Group')
 870  format('STRUCTURE',7x,'MEASURE')

c-------------------------------------------------------
c Farem el calcul per les diferents estructures distorsionades
c
c	IF (csd) THEN
c  	  iunit=11
c	  open(unit=iunit,status='old',file=trim(input)//".cor")
c	  write(3,"('Reading data from file ',a)") trim(input)//'.cor'
c	ENDIF
c-------------------------------------------------------

c	STRUCTURES: do es=1,9999999
	STRUCTURES: do es=1,1

c	 call nextline(iunit)
c	 read(iunit,"(a)",end=11, IOSTAT=iostat) name
c	 if (iostat.gt.0) exit
c 	 write(3,'(a10)') name
c	 write(3,*)' '
c	 ivecs(:)=ivecs0(:)
c	 k=0
c	if(atfix) then
c         do j=1,nvecs
c          if(ivecs(j).eq.pm) ivecs(j)=m
c         end do
c	end if
c	 do i=1,m
c	  if (i.eq.pm) then
c	   write(6,*) 'pm', pm
c	   call nextline(iunit)
c	   if(nolabel)then
c	    read(iunit,*)(q0(m,j),j=1,3)
c	   else
c	    read(iunit,"(a)")line
c            line=adjustl(line)
c            l=index(line,' ')
c	    etiq0(m)=line(1:l)
c	    line(1:l)=' '
c	    read(line,*)(q0(m,j),j=1,3)
c	   end if
c	  else
c	   k=k+1
c           call nextline(iunit)
c	   if(nolabel)then
c            read(iunit,*)(q0(k,j),j=1,3)
c           else
c	    read(iunit,"(a)")line
c            line=adjustl(line)
c            l=index(line,' ')
c            etiq0(k)=line(1:l)
c            line(1:l)=' '
c            read(line,*)(q0(k,j),j=1,3)
c	    if (atfix) then
c	     do j=1,nvecs
c	      if(ivecs(j).eq.i) ivecs(j)=k
c	     end do
c	    end if
c	   end if
c	  end if
c	 end do

c      write(6,*)  es
c
c
c   Coloquem la molecula al centre de masses
c
	 call cmasses(q0,q,m,centerfix)
c	 write(3,*)'Centered structure'
c	 if(nolabel) then
c	  do i=1,m
c	   write(3,"(11x,f9.5,2x,f9.5,2x,f9.5)")(q(i,j),j=1,3)
c	  end do
c	 else
c	  do i=1,m
c           write(3,"(1x,a8,2x,f9.5,2x,f9.5,2x,f9.5)")etiq0(i),
c     .    (q(i,j),j=1,3)
c	  end do
c	 end if
c
c CONECTIVITAT
c
	if(connect) then
	 call connectivity(m,etiq0,q,vcon)
C--MQL --
c Escriptura
c	 write(3,"(/2x,'Connectivity matrix')")
c         do ic1=1,m,8
c           ic2=min(m,ic1+7)
c	   write(3,"(/17x,8(i4,4x))") (j,j=ic1,ic2)
c           write(3,"(17x,8(a8))") etiq0(ic1:ic2)
c	   do i=ic1,m
c             kk=((i+1)*i)/2-1
c             kc1=kk+ic1
c             kc2=kk+min(ic2,i)
c	     write(3,"(1x,i7,1x,a8,8(2x,i1,5x))")i,etiq0(i),int(vcon(kc1:kc2))
c	   end do
c         enddo
C--MQL --/
	 call jacobi(vcon,m,convap,convep,6)
	end if

c Group
c	 call norma(q,q,m,A)
c 	 p(:,:)=q(:,:)
       call norma(q,q,1,m,A)
       do l=1,noper
         p(l,:,:)=q(:,:)
       end do
	 etiq(:)=etiq0(:)

c
c Calcul de vmax fixat per alguns atoms
c
	if(atfix)then
	 if((oper.eq.'r').and.(nvecs.eq.2))then
	  dvmax=dot_product(q(ivecs(1),:),q(ivecs(2),:))
	  dvmax=dvmax/dsqrt(q(ivecs(1),1)**2+
     .    q(ivecs(1),2)**2+q(ivecs(1),3)**2)
	  dvmax=dvmax/dsqrt(q(ivecs(2),1)**2+
     .    q(ivecs(2),2)**2+q(ivecs(2),3)**2)
	  if(abs(1.0d0-abs(dvmax)).lt.1.0d-6) then
           write(0,"(1x,a8,
     .   2x,'ERROR. Aligned vertices do not define a plane')") name
c           write(0,*)' '
c           write(0,*)'ERROR. Non planar structure'
c           write(0,*)'**************************************'
	   goto 11
 	  else
	   vmax(1)=q(ivecs(1),2)*q(ivecs(2),3)-
     .     q(ivecs(1),3)*q(ivecs(2),2)
           vmax(2)=q(ivecs(1),3)*q(ivecs(2),1)-
     .     q(ivecs(1),1)*q(ivecs(2),3)
           vmax(3)=q(ivecs(1),1)*q(ivecs(2),2)-
     .     q(ivecs(1),2)*q(ivecs(2),1)
	  end if
	 else
	  vmax(:)=0.0d0
	  do j=1,3
	   do k=1,nvecs
	    vmax(j)=vmax(j)+q(ivecs(k),j)
	   end do
	  end do
	 end if
	end if
c----------------------------------
        if (dim2D) then
         call vecnormal(q,m,vecnorm,pla_ctrl)
         if (.not.pla_ctrl) then
	  write(0,"(1x,a8,2x,'ERROR. Non planar structure')") name
c	  write(3,*)' '
c          write(3,*)'ERROR. Non planar structure'
c	  write(3,*)'**************************************'
          goto 11
         end if
	 if((oper.eq.'c').or.(oper.eq.'s')) vmax(:)=vecnorm(:)
        end if
c----------------------------------
        if (i_oper.eq.4) then
         dvmax=vmax(1)**2+vmax(2)**2+vmax(3)**2
         if(dvmax.lt.1.d-6) then
          write(0,*)'ERROR. Bad fixed axis definition'
          write(0,*)'Defined vector length equal zero'
c         write(0,"(1x,a8,5x,'ERROR. Bad fixed axis definition')") name
          cycle STRUCTURES
         end if
         vmax(:)=vmax(:)/dsqrt(dvmax)
        end if
c----------------------------------
        s_max=-1.0d0
c----------------------------------
      IF (i_oper .EQ. 0) THEN
        s_max=A
        p(1,nll,:)=q(ipf(nll),:)  ! Group
        qf(:,:) = q(:,:)
      ENDIF

c----------------------------------
c Permutacio unica
	IF (permu) THEN
        DO k=1,m
          p(1,k,:)=q(ipf(k),:) ! Group
        ENDDO
c	  call mesura(i_oper,q,p,m,s_max,vmax)
	  call mesura(oeix,i_oper,1,q,p,m,s_max,vmax,
     .                  ctheta,stheta,ttheta) ! Group
	  goto 101
	ENDIF
c----------------------------------
	ctrl(:nll)=.false.
      ip(:nll)=0
      idc=1
	ii=1
	sm(1) = -1.d6
c
C--MQL
        anpmtt=1.d0
        kcon(:m)=0
        jcon(:m,:m)=0
        ncon=0
        if (connect) then
          do ih=1,nll
            kdone=0
            do icon=1,ncon
      if (abs(convep(ih,1)-convep(jcon(icon,1),1)).le.1.d-6) then
                kcon(icon)=kcon(icon)+1
                anpmtt=anpmtt*kcon(icon)
                jcon(icon,kcon(icon))=ih
                kdone=1
                exit
              endif
            enddo
            if (kdone.eq.0) then
              ncon=ncon+1
              kcon(ncon)=1
              jcon(ncon,1)=ih
            endif
          enddo
        elseif (label) then
          do ih=1,nll
            kdone=0
            do icon=1,ncon
              if (trim(etiq(ih)).eq.trim(etiq(jcon(icon,1)))) then
                kcon(icon)=kcon(icon)+1
                anpmtt=anpmtt*kcon(icon)
                jcon(icon,kcon(icon))=ih
                kdone=1
                exit
              endif
            enddo
            if (kdone.eq.0) then
              ncon=ncon+1
              kcon(ncon)=1
              jcon(ncon,1)=ih
            endif
          enddo
        else
          ncon=1
          kcon(1)=nll
          do ih=1,nll
            jcon(1,ih)=ih
            anpmtt=anpmtt*ih
          enddo
        endif
c        if (.not.mute) write(3,"(/'Sets of equivalent points: ',i3)") ncon
c        write(3,"(/'Sets of equivalent points')")
        do ikc=maxval(kcon(:ncon)),1,-1
          nkc=count(kcon(:ncon).eq.ikc)
          if (nkc.ne.0) then
c            write(3,"('Sets of ',i4,' points :',i3)") ikc, nkc
            iapl=0
            do icon=1,ncon
              if (kcon(icon).eq.ikc) then
c                write(3,"(' |',100i4)",advance='no') jcon(icon,1:ikc)
                iapl=iapl+ikc
                if (iapl+ikc.gt.13) then
c                  write(3,"(' |')")
                  iapl=0
                endif
              endif
            enddo
c            if (iapl.ne.0) write(3,"(' |')")
          endif
        enddo
c        if (.not.mute)
c     $ write(3,"('Number of predicted permutations',g10.2)") anpmtt
c    $ write(3,"(/'Number of predicted permutations',g10.2)") anpmtt
        ktim=0
        khg=huge(khg)
        apmt(1)=1.d0
        do im=2,4
          apmt(im)=apmt(im-1)*khg
        enddo
        npmt(:)=0
        kpct0=0
        call system_clock (icl_ini,icl_rat,icl_top)
C--MQL ///

      PERMUT: do while (idc.ne.-1)
c
      ip(idc)=ip(idc)+1
      if (ip(idc).gt.nll) then
        ip(idc)=0
        idc=idc-1
        if (idc.eq.0) exit PERMUT
        ctrl(ip(idc))=.false.
        cycle PERMUT
      endif
      if (ctrl(ip(idc))) cycle PERMUT

	if (label) then
 	  if(trim(etiq(ip(idc))).ne.
     .   trim(etiq(idc))) cycle PERMUT
	end if
	if (connect) then
	 if(abs(convep(ip(idc),1)-convep(idc,1))
     &.gt.1.d-6) cycle PERMUT
	end if

c---------------------------c
c Permutations of the Group c
c---------------------------c
c Inversion
      if (i_oper.eq.1) then
       if(nll/2*2.eq.nll) then
        if (ip(idc).eq.idc) cycle PERMUT
       else
        if (ip(idc).eq.idc) then
         noperm=0
         do k=1,idc
          if (ip(k).eq.k) noperm=noperm+1
         end do
         if (noperm.ne.1) cycle PERMUT
        end if
       end if
      end if

      if (idc.ne.nll) then
        ctrl(ip(idc))=.true.
        idc=idc+1
        cycle PERMUT
      endif

c Inversion and Reflection
      if (i_oper.le.2) then
       do k=1,nll
        if (ip(ip(k)).ne.k) cycle PERMUT
       end do
      end if

c Rotation
      if((i_oper.eq.3).or.(i_oper.eq.4))then
       if (oeix.eq.3) then
        do k=1,nll
        if (ip(ip(ip(k))).ne.k) then
         if((oper.eq.'s ').or.(oper.eq.'sf')) then
          if(ip(ip(ip(ip(ip(ip(k)))))).ne.k) cycle PERMUT
         else
          cycle PERMUT
         end if
        end if
        end do
       end if
       if (oeix.eq.4) then
        do k=1,nll
         if(ip(ip(ip(ip(k)))).ne.k) then
          if((ip(ip(k)).eq.k).and.(ip(k).ne.k)) cycle PERMUT
          if((oper.eq.'s ').or.(oper.eq.'sf')) then
           if(ip(ip(ip(ip(ip(ip(ip(ip(k)))))))).ne.k) cycle PERMUT
          else
           cycle PERMUT
          end if
         end if
        end do
       end if
      end if

c Compute measure for each permutation
c Group
	DO k = 1,nll
	  p(1,k,:) = q(ip(k),:)
	ENDDO
      IF ((i_oper.ge.3).and.(oeix.ge.3)) THEN
        DO l=2,noper
        DO k=1,nll
          p(l,k,:)=p(l-1,ip(k),:)
        ENDDO
        ENDDO
      ENDIF

c	call mesura(i_oper,q,p,m,s,v)
 	CALL mesura(oeix,i_oper,noper,q,p,m,s,v,
     & ctheta,stheta,ttheta)

C--MQL
       if (ktim.eq.0.and.npmt(1).eq.100) then
          call system_clock (icl_fin,icl_rat,icl_top)
          icl_stp=icl_fin-icl_ini
          if (icl_stp.lt.0) icl_stp=icl_stp+icl_top
          cl_seg=(anpmtt*icl_stp/icl_rat)/100.d0
          cl_hor=cl_seg/3600.d0
          cl_dia=cl_hor/24.d0
c          if (.not.mute) then
c           if (cl_dia.gt.1.d0) then
c       write(3,"(/'Estimated calculation time ',g9.2,' days')") cl_dia
c       write(3,"('** WARNING ** Time larger than 1 day')")
c           elseif (cl_hor.gt.1.d0) then
c       write(3,"(/'Estimated calculation time ',f8.3,' hours')") cl_hor
c       write(3,"('** WARNING ** Time larger than 1 hour')")
c           else
c            write(3,"(/'Estimated calculation time ',
c     &            f8.3,' seconds')") cl_seg
c           endif
c DCC: this version allways runs
c           if (cl_seg.gt.300) then
c             write(3,"('Continue (y/n) ? ')",advance='no')
c             read(5,*) yn
c             if (yn.ne.'y') stop "Aborted calculation"
c           endif
c            if (anpmtt.gt.apmt(4)*khg) write(3,
c     &"(/'** WARNING ** Dot counter will stop before the end')")
c          endif
c          ktim=1
c          call system_clock (icl_ini,icl_rat,icl_top)
       endif
       if (ktim.ne.3) then
          npmt(1)=npmt(1)+1
          if (npmt(1).eq.khg) then
            npmt(1)=0
            npmt(2)=npmt(2)+1
            if (npmt(2).eq.khg) then
            npmt(2)=0
            npmt(3)=npmt(3)+1
            if (npmt(3).eq.khg) then
            npmt(3)=0
            npmt(4)=npmt(4)+1
            if (npmt(4).eq.khg) then
            ktim=3
            if (.not.mute) then
             write(0,"(/'** WARNING ** Permutation number too large.')")
             write(0,"( '           ** Counter stopped.')")
             write(0,"( '           ** Program is running.')")
            endif
            endif
            endif
            endif
          endif
          anpmt=sum(npmt(:)*apmt(:))
          kpct=int(50.0d0*anpmt/anpmtt)
          if (ktim.eq.1) then
           if (.not.mute) write(0,"(1a)") repeat('_',50)
           ktim=2
          endif
          if (.not.mute.and.ktim.eq.2.and.kpct.ne.kpct0)
     $     write(0,"(1a)",advance='no') '.'
          kpct0=kpct
        endif
C--MQL ///
c
c	Multiples valors optims
c
	if(s.ge.sm(ii))then
	 if (ii.gt.15) ii=15
	 jj=1
	 do i=1,ii
	  if (s.lt.sm(i)) jj=i+1
	 end do
	 ii=ii+1
	 if (ii.gt.jj) then
	  do ij=ii,jj+1,-1
	   sm(ij)=sm(ij-1)
	   vm(ij,:)=vm(ij-1,:)
	  end do
	 end if
	 sm(jj)=s
	 vm(jj,:)=v(:)
	end if

c Renew optimum value
       IF (s.gt.s_max) THEN
       s_max=s
        vmax=v
       ipf(:)=ip(:)
       ENDIF
       enddo PERMUT

C--MQL
c      if (ktim.ne.3) then
c        if (npmt(2).eq.0) then
c         if (.not.mute) write(6,"(//'Number of permutations done',i20)") npmt(1)
c         write(3,"('Number of permutations done',i20)") npmt(1)
c        write(3,"(/'Number of permutations done',i20)") npmt(1)
c        else
c         if (.not.mute) write(6,"(//'Number of permutations done',g8.2)") anpmt
c         write(3,"('Number of permutations done',g8.2)") anpmt
c        endif
c      endif
      if (ktim.gt.0) then
        sec_max=1.d0*icl_top/icl_rat
        if (cl_seg.gt.sec_max*0.9d0) then
c         if (.not.mute) write(6,"(/'Real calculation time > 1 day')")
          write(0,"(/'Real calculation time > 1 day')")
        else
          call system_clock (icl_end,icl_rat,icl_top)
          icl_tot=icl_end-icl_ini
          if (icl_tot.lt.0) icl_stp=icl_tot+icl_top
          cl_seg=1.0d0*icl_tot/icl_rat
          cl_hor=cl_seg/3600.d0
          cl_dia=cl_hor/24.d0
          if (cl_seg.lt.1000.d0) then
c           if (.not.mute) write(6,"(/'Real calculation time ',f10.4,' seconds')") cl_seg
            write(0,"(/'Real calculation time ',f10.4,' seconds')") cl_seg
          else
c           if (.not.mute) write(6,"(/'Real calculation time ',f10.4,' hours')") cl_hours
            write(0,"(/'Real calculation time ',f10.4,' hours')") cl_hours
          endif
        endif
c        call system('date')
      endif
C--MQL ///
 101	continue
c----------------------------------
c Escriptura estructura final
c
C--MQL
c      write(3,"('Optimal permutation')")
c      do im=1,nll
c        write(3,"(2i5)") im, ipf(im)
c      enddo
c     write(3,"(/)")
C--MQL ///

      DO k=1,nll
c       p(k,:)=q(ipf(k),:)
        ptotal(1,k,:)=q(ipf(k),:) ! Group
	  etiq(k)=etiq0(ipf(k))
      ENDDO
c Group
      IF ((oper.eq.'c').or.(oper.eq.'s')) THEN
        DO l=2,oeix-1
        DO k=1,nll
          ptotal(l,k,:)=ptotal(l-1,ipf(k),:)
        ENDDO
        ENDDO
	ENDIF

	IF (oper.eq.'i') THEN
c       write(3,*)'Inverted structure'
c       qf(:,:)=-p(1,:,:)
c Group
c        write(3,*)'Nearest Ci Structure'
        qf(:,:)=0.5d0*(q(:,:)-ptotal(1,:,:))
	ELSEIF (oper.eq.'r') THEN
c       write(3,*)'Reflected structure'
c       call reflexio(p,m,vmax,qf)
c Group
c        write(3,*)'Nearest Cs Structure'
        call reflexio(ptotal,m,vmax,qf)
        qf(:,:)=0.5d0*(q(:,:)+qf(:,:))
	ELSEIF ((oper.eq.'c').or.(oper.eq.'s')) THEN
c Group
c	  call rotacio(p,m,vmax,qf)
        call rotacio(ptotal,oeix,m,vmax,qf,ctheta,stheta,ttheta)
        qf(:,:)=(qf(:,:)+q(:,:))/dfloat(oeix)
c	  IF (oper.eq.'c') THEN
c Group
c         write(3,*)'Rotated structure'
c          write(3,"('Nearest C',i2,1x,'structure')") oeix
c	  ELSE
c Group
c         write(3,*)'Rotated and reflected structure'
c          write(3,"('Nearest S',i2,1x,'structure')") oeix
c	  ENDIF
	ENDIF
c
c	if (oper.ne.'e') then
c	 do k=1,m
c	  if(nolabel) then
c	   write(3,"(11x,f9.5,2x,f9.5,2x,f9.5)")(qf(k,j),j=1,3)
c	  else
c	   write(3,"(1x,a8,2x,f9.5,2x,f9.5,2x,f9.5)")etiq(k),
c     .    (qf(k,j),j=1,3)
c	  end if
c	 end do
c         write(3,*)' '
c	end if
c
c	if (i_oper.eq.2) then
c         write(3,*)'Normal vector to the optimum reflection plane'
c        end if
c        if (i_oper.eq.3) then
c          write(3,*)'Optimum axis orientation'
c        end if
c        if (i_oper.eq.4) then
c          write(3,*)'Fixed axis orientation'
c        end if

	IF ((i_oper.eq.2).or.(i_oper.eq.3).or.
     $   (i_oper.eq.4))then
c	 write(3,"(3(2x,f9.5))") (vmax(j),j=1,3)
C--MQL
c        write(3,"(/1x,'Reference axis')")
        hz=-(vmax(1)+vmax(2))/vmax(3)
        hx=1.d0/sqrt(2.d0+hz*hz)
        hy=hx
        hz=hz*hx
        h2x=vmax(2)*hz-vmax(3)*hy
        h2y=vmax(3)*hx-vmax(1)*hz
        h2z=vmax(1)*hy-vmax(2)*hx
c        write(3,"(1x,'X',7x,3(2x,f9.5))") 0.d0, 0.d0, 0.d0
c        write(3,"(1x,'X',7x,3(2x,f9.5))") hx, hy, hz
c        write(3,"(1x,'X',7x,3(2x,f9.5))") h2x, h2y, h2z
c        write(3,"(1x,'X',7x,3(2x,f9.5))") (vmax(j),j=1,3)
C--MQL ///
	ENDIF

c Group
c	s_max=50.0d0*(1.0d0-s_max/A)
	IF (i_oper.le.2) THEN
        s_max = 50.0d0 * (1.0d0-s_max/A)
	ELSE
        s_max = (1.0d0+s_max/A)/dfloat(oeix)
	  s_max = 100.0d0 * (1.0d0-s_max)
	ENDIF

c 	write(2,"(1x,a10,',',2x,f9.5)") name,s_max
c	write(3,*)' '
c	write(3,"(1x,'Symmetry measure',2x,f14.8)") s_max
c	write(3,*)'**************************************'
C--MQL
        if (multi) then
C--MQL ///

	DO i = 1,Nmulti
       IF (i_oper.le.2) THEN
         sm(i) = 50.0d0 * (1.0d0-sm(i)/A)
       ELSE
         sm(i) = (1.0d0+sm(i)/A)/dfloat(oeix)
         sm(i) = 100.0d0 * (1.0d0-sm(i))
       ENDIF
	ENDDO
c	write(4,'(1x,a10)') name
c	write(4,*)' '
c	if (Nmulti.le.5) then
c	 write(4,"(5(2x,f9.5))")(sm(j),j=1,Nmulti)
c	 write(4,*)' '
c         do k=1,3
c          write(4,"(5(2x,f9.5))")(vm(j,k),j=1,Nmulti)
c         end do
c	else
c	 write(4,"(5(2x,f9.5))")(sm(j),j=1,5)
c         write(4,*)' '
c         do k=1,3
c          write(4,"(5(2x,f9.5))")(vm(j,k),j=1,5)
c         end do
c         write(4,*)' '
c         write(4,*)' '
c	 if (Nmulti.le.10) then
c          write(4,"(5(2x,f9.5))")(sm(j),j=6,Nmulti)
c          write(4,*)' '
c          do k=1,3
c           write(4,"(5(2x,f9.5))")(vm(j,k),j=6,Nmulti)
c          end do
c	 else
c          write(4,"(5(2x,f9.5))")(sm(j),j=6,10)
c          write(4,*)' '
c          do k=1,3
c           write(4,"(5(2x,f9.5))")(vm(j,k),j=6,10)
c          end do
c          write(4,*)' '
c          write(4,*)' '
c          write(4,"(5(2x,f9.5))")(sm(j),j=11,Nmulti)
c          write(4,*)' '
c          do k=1,3
c           write(4,"(5(2x,f9.5))")(vm(j,k),j=11,Nmulti)
c          end do
c	 end if
c	end if
c        write(4,*)'**************************************'
C--MQL
        endif
C--MQL ///
c
c------------------------------------------------
c Fi del calcul per les diferents estructures
c
C--MQL --
        enddo STRUCTURES

   11   if (.not.mute) write(3,*)' '
C--MQL --

c 10   write(3,"(/,' *** Normal Termination ***')")
c      write(3,*)' '
c      write(3,*)' Thank you for using SYMGROUP '
c      write(3,*)' '

c      close(2)
c      close(3)
c   	if (multi) close(4)
c	if (csd) close(11)
c Recolecta outputs

      iOff = 0
      do i =1, m
        if (i.eq.pm) then
          NearStruc(i,:) = qf(m, :)
          iOff = -1
        else
          NearStruc(i,:) = qf(i+iOff, :)
        endif
      enddo
      v_max = vmax

      do i=1,nll
        optper(i) = ipf(i)
      enddo

      refAxis(1,:) = (/0.d0, 0.d0, 0.d0/)
      refAxis(2,:) = (/hx, hy, hz/)
      refAxis(3,:) = (/h2x, h2y, h2z/)
      refAxis(4,:) = vmax(:)

c      write(3,"(1x,'X',7x,3(2x,f9.5))") 0.d0, 0.d0, 0.d0
c      write(3,"(1x,'X',7x,3(2x,f9.5))") hx, hy, hz
c      write(3,"(1x,'X',7x,3(2x,f9.5))") h2x, h2y, h2z
c      write(3,"(1x,'X',7x,3(2x,f9.5))") (vmax(j),j=1,3)
c      sm, vm

c 9999 close(1)
 9999 END





c********************c
c*                  *c
c*   SUBROUTINES    *c
c*                  *c
c********************c

	subroutine norma(a,b,n,m,c)

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
c-------------------------------------------------------
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
c------------------------------------------------------------------------
c	Subroutines Lectura d'Inputs
c------------------------------------------------------------------------
	subroutine nextline(iunit)
c
c Aquesta subroutine detecta les linies d'input, evitant linies en blanc i comentaris
c
	integer    :: iunit
	character  :: line*80
	logical    :: ctrl
c
	ctrl=.false.
	do while (.not.ctrl)
	 read(iunit,"(a)",end=999) line
	 line=adjustl(line)
	 if((line(1:1).ne.'!').and.(len_trim(line).ne.0)) then
	  backspace(iunit)
	  ctrl=.true.
	 end if
	end do
c
 999 	return
	end subroutine
c------------------------------------------------------------------------
	subroutine options(opt)
c
c Activa les diferents opcions segons les keywords emprades
c
      character  :: opt0*80, opt*80, line*80
      logical    :: label,multi,hconv,permu,atfix,
     & rfix,dim2D,csd,connect,nolabel,centerfix
c
      common /intlog/Nmulti
C--MQL
      logical :: mute
      double precision :: thresh, cm(3)
      common /sumrad/thresh
      common /center/cm
C--MQL ///
C--MQL --
      common /logical/label,multi,hconv,permu,atfix,
     & rfix,dim2D,csd,connect,nolabel,mute,centerfix
C--MQL --/

      opt0=opt
      k=index(opt0,' ')
      opt0=opt0(1:k)//repeat(' ',80-k)
c Passem les minuscules a majuscules
      do i=2,len_trim(opt0)
        k=ichar(opt0(i:i))
        if (k.ge.97.and.k.le.122) then
          opt0(i:i)=char(k-32)
        endif
      enddo

      select case (opt0)

	case('%LABEL')
	 label=.true.
        case('%HCONV')
         hconv=.true.
        case('%CONNECT')
         connect=.true.
C--MQL
         if (len_trim(opt).eq.len_trim(opt0)) then
           thresh=1.1d0
         else
	   backspace(1)
	   read(1,"(a8)",advance='no')line(1:8)
	   read(1,*) thresh
           thresh=thresh/100.d0
         endif
        case('%MUTE')
         mute=.true.
C--MQL ///
        case('%PERMU')
         permu=.true.
	case('%RFIX')
	 rfix = .true.
      case('%CENTERFIX')
        k = len_trim(opt0)
        opt0 = repeat(' ',k)//opt(k+1:80)
        read(opt0,*) cm(1),cm(2),cm(3)
        centerfix = .true.
      case('%ATFIX')
       atfix=.true.
	case('%2D')
	 dim2D=.true.
	case('%CSD')
	 csd=.true.
	case('%NOLABEL')
	 nolabel=.true.
      case('%MULTI')
	  multi=.true.
	  backspace(1)
	  read(1,"(a6)",advance='no')line(1:6)
	  read(1,*)Nmulti

      case default
          write(3,"(/,'* ERROR * Unknown keyword')")
	    stop
      endselect

  999 return

      end subroutine
c------------------------------------------------------------------------
        subroutine m_control(nll,pm)
c
c Controla certs valors d'entrada incorrectes
c
      integer    :: nll,pm
      logical    :: label,multi,hconv,permu,atfix,
     . rfix,dim2D,connect
c
C--MQL --
      logical :: mute
      common /logical/label,multi,hconv,permu,atfix,
     & rfix,dim2D,csd,connect,nolabel,mute, centerfix
C--MQL --/
c
	m=nll
        if(pm.ne.0) m=nll+1
        if (m.lt.3) then
         write(0,999)
         stop
        end if
	if(pm.gt.nll+1) then
	 write(0,*)'ERROR: Not allowed value of the central atom position'
	 stop
	end if
c
	if((permu).and.(pm.ne.0)) then
	 write(0,*)'ERROR. Central atom is not allowed with permu'
	 stop
	end if
c
 999  format('ERROR: Structures with less than 3 points have no sense')
	return
	end subroutine
c------------------------------------------------------------------------
