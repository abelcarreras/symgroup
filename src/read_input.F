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
         write(3,999)
         stop
        end if
	if(pm.gt.nll+1) then
	 write(3,*)'ERROR: Not allowed value of the central atom position'
	 stop
	end if
c
	if((permu).and.(pm.ne.0)) then
	 write(3,*)'ERROR. Central atom is not allowed with permu'
	 stop
	end if
c
 999  format('ERROR: Structures with less than 3 points have no sense')
	return
	end subroutine
c------------------------------------------------------------------------
