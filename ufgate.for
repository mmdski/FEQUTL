C
C
C
      REAL FUNCTION   RSOMY3
     I                      (Y3U)
 
C     + + + PURPOSE + + +
C     Find the unknown depth at section 3 for submerged orifice
C     flow from the momentum balance.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y3U
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y3U    - value of unknown being sought-depth at section 3
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ufgate.cmn'
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTJ
C***********************************************************************
      YT = Y3U
      CALL LKTJ
     I         (DEPTAB,
     M          YT,
     O          J4ATY3)
 
      QSQR = AT**2*(TWOG*(Z1B + Y1 - Y3U - Z3B))/
     A              (1.0 - ALPHA1*(AT/A1)**2)
 
      RSOMY3 = (G*(J4ATY3 - J4) + QSQR*(1.0/(CC*AG) -
     A                      BETA4/A4))/(G*J4)
 
      RETURN
      END
C
C
C
      REAL FUNCTION   RSOMY4
     I                      (Y4U)
 
C     + + + PURPOSE + + +
C     Find the unknown depth at section 4 for submerged orifice
C     flow from the momentum balance.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y4U
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y4U    - value of unknown being sought-depth at section 3
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ufgate.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL DT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKT22
C***********************************************************************
      YT = Y4U
      CALL XLKT22
     I           (DEPTAB,
     M            YT,
     O            A4, TT, DT, J4, KT, DKT, BETA4, DBETA, ALPHA4, DALPHA,
     O            QCT)
 
      RSOMY4 = (G*(J4ATY3 - J4) + QSQR*(1.0/(CC*AG) -
     A                      BETA4/A4))/(G*J4)
 
      RETURN
      END
C
C
C
      REAL FUNCTION   RSWMY3
     I                      (Y3U)
 
C     + + + PURPOSE + + +
C     Find the unknown depth at section 3 for submerged weir flow
C     from the momentum balance.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y3U
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y3U    - value of unknown being sought-depth at section 3
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ufgate.cmn'
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTJ
C***********************************************************************
      YT = Y3U
      CALL LKTJ
     I         (DEPTAB,
     M          YT,
     O          J4ATY3)
 
      QSQR = TWOG*(Z1B + Y1 - Z3B - Y3U)/(1.0/(CD*BG*(Y3U + DZ))**2 -
     A                 ALPHA1/A1**2)
      RSWMY3 = (G*(J4ATY3 - J4) + QSQR*(1.0/((Y3U + DZ)*BG) -
     A                      BETA4/A4))/(G*J4)
 
      RETURN
      END
C
C
C
      REAL FUNCTION   RSWMY4
     I                      (Y4U)
 
C     + + + PURPOSE + + +
C     Find the value of depth at section 4 for submerged weir
C     flow from the momentum balance.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y4U
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y4U    - value of unknown being sought-depth at section 3
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ufgate.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL DT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKT22
C***********************************************************************
C     On entry all constant values must be known and placed in
C     the common block in ufgate.cmn. The cross section elements at
C     section 4 are recomputed here.
      YT = Y4U
      CALL XLKT22
     I           (DEPTAB,
     M            YT,
     O            A4, TT, DT, J4, KT, DKT, BETA4, DBETA, ALPHA4, DALPHA,
     O            QCT)
 
      RSWMY4 = (G*(J4ATY3 - J4) + QSQR*(1.0/((Y3 + DZ)*BG) -
     A                      BETA4/A4))/(G*J4)
      RETURN
      END
C
C
C
      REAL FUNCTION   FINDCC
     I                      (ARG, CONCC, CCTAB)
 
C     + + + PURPOSE + + +
C     Find a contraction coefficient for an underflow gate.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CCTAB
      REAL ARG, CONCC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARG    - argument for lookup of contraction coefficient for an
C              underflow gate.
C     CONCC  - contraction coefficient
C     CCTAB  - address of table giving contraction coefficient
 
C     + + + LOCAL VARIABLES + + +
      INTEGER NTAB
      REAL CC, DF
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTAB
C***********************************************************************
      IF(CCTAB.EQ.0) THEN
C       Take as a constant value.
        FINDCC = CONCC
      ELSE
C       Lookup values in table.
        CALL LKTAB
     I            (CCTAB, ARG, 1,
     O             CC, NTAB, DF)
        FINDCC = CC
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   FOTOSO
     I                   (H1, Q, HDATUM,
     O                    H4, IFLAG)
 
C     + + + PURPOSE + + +
C     For a given upstream head(at section 1) find the
C     tailwater head(at section 4) that defines the
C     boundary between free orifice(FO) flow and submerged orifice(SO)
C     flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER IFLAG
      REAL H1, H4, HDATUM, Q
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     H1     - head at section 1
C     Q      - Flowrate
C     HDATUM - Datum for measuring head
C     H4     - head at section 4
C     IFLAG  - error flag
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ufgate.cmn'
      INCLUDE 'epscom.cmn'

C     + + + LOCAL VARIABLES + + +
      REAL HVC, YLEFT, YRIGHT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL RSOMY4
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTJ, REGFAL, RSOMY4
C***********************************************************************
      HVC = CC*HG + Z2B - HDATUM
      Y3 = HVC + HDATUM - Z3B
      CALL LKTJ
     I         (DEPTAB,
     M          Y3,
     O          J4ATY3)
      QSQR = Q**2
 
C     The water surface elevation at section 4 will be
C     below the elevation at section 1 but above the
C     jet surface elevation at section 2.
 
      YLEFT = HVC + HDATUM - Z4B
      YRIGHT = H1 + HDATUM - Z4B

      CALL REGFAL
     I           (EPSARG, EPSF, RSOMY4,
     M            YLEFT, YRIGHT,
     O            Y4, IFLAG)
 
      H4 = Y4 + Z4B - HDATUM
      RETURN
      END
C
C
C
      SUBROUTINE   FNDFOQ
     I                   (H1, HDATUM,
     O                    Q)
 
C     + + + PURPOSE + + +
C     Find the free orifice flow.  Gate opening and contraction
C     coefficient are in common block
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL H1, HDATUM, Q
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     H1     - head at section 1
C     HDATUM - Datum for measuring head
C     Q      - Flowrate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ufgate.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL DT, HVC
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKT22
C***********************************************************************
      Y1 = H1 + HDATUM - Z1B
      CALL XLKT22
     I           (APPTAB,
     M            Y1,
     O            A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1, DALPHA,
     O            QCT)
 
      AT = CD*CC*AG
      HVC = CC*HG + Z2B - HDATUM
      Q = AT*SQRT(TWOG*(H1 - HVC)/
     A              (1.0 - ALPHA1*(AT/A1)**2))
      RETURN
      END
C
C
C
      SUBROUTINE   LSTOPF 
     I                   (STDOUT, TABLT, TABGT, POWER, OFFSET, A, B,
     I                    PREC, NMAX,
     O                    N, XBRK, EFLAG)
 
C     + + + PURPOSE + + +
C     Define the sequence of breakpoints, XBRK, between A and B
C     so that a linear spline on that breakpoint sequence will
C     interpolate a power function with power=POWER and offset=
C     OFFSET with a relative precision given by PREC.  There
C     can be no more than NMAX points.  TABLT and TABGT are
C     data table addresses that give the argument ratios for a range
C     of powers and precisions.  TABLT is for powers less than 1
C     and TABGT is for powers greater than 1.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, N, NMAX, STDOUT, TABGT, TABLT
      REAL A, B, OFFSET, POWER, PREC, XBRK(NMAX)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     TABLT  - address of function table for powers less than 1.0
C     TABGT  - function table address for powers greater than 1.0
C     POWER  - power in the power function
C     OFFSET - offset to use in the power function
C     A      - lower limit of range of approximation
C     B      - upper limit of range of approximation
C     PREC   - precision of the approximation
C     NMAX   - Maximum number of items allowed in table
C     N      - number of break points
C     XBRK   - breakpoint locations for the linear spline
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IPOW
      REAL ARGRAT, DFCOL, DFROW, POW, RBTOA
 
C     + + + INTRINSICS + + +
      INTRINSIC EXP, FLOAT, LOG
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL TDLK10

C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* A=',F10.4,' => B=',F10.4,' in LSTOPF.')
C***********************************************************************
C     Find the argument ratio for the current power and precision.
C     Note that TABLT and TABGT are addresses of the table and not
C     the table number.

      IF(A.GE.B) THEN
        WRITE(STDOUT,50) A, B
        STOP
      ENDIF 
      IF(POWER.LT.1.0) THEN
        CALL TDLK10
     I             (STDOUT, TABLT, 10, POWER, PREC,
     O              ARGRAT, DFROW, DFCOL)
      ELSE
        CALL TDLK10
     I             (STDOUT, TABGT, 10, POWER, PREC,
     O              ARGRAT, DFROW, DFCOL)
      ENDIF
 
C      WRITE(STDOUT,*) ' LSTOPF: ARGRAT=',ARGRAT
 
C      WRITE(STDOUT,*) ' LSTOPF: A=',A,' B=',B
C      WRITE(STDOUT,*) ' OFFSET=',OFFSET
 
C     Find the ratio of the limits
      RBTOA = (B + OFFSET)/(A + OFFSET)
 
C      WRITE(STDOUT,*) ' LSTOPF: RBTOA=',RBTOA
C     Find the non-integer power of ARGRAT that yields RBTOA
 
      POW = LOG(RBTOA)/LOG(ARGRAT)
      IPOW = POW + 1.0
      N = IPOW + 1
C      WRITE(STDOUT,*) ' LSTOPF: N=',N
C     Recompute ARGRAT so that an integer power of it yields RBTOA
 
      ARGRAT = EXP(LOG(RBTOA)/FLOAT(IPOW))
 
C      WRITE(STDOUT,*) ' LSTOPF: NEW ARGRAT=',ARGRAT
C     Compute the sequence of values in XBRK
 
      IF(N.GT.NMAX) THEN
        EFLAG = 1
        N = -1
      ELSE
        XBRK(1) = A + OFFSET
        XBRK(N) = B + OFFSET
        DO 100 I=2,N-1
          XBRK(I) = ARGRAT*XBRK(I-1)
 100    CONTINUE
        IF(OFFSET.NE.0.0) THEN
          DO 200 I=1,N
            XBRK(I) = XBRK(I) - OFFSET
 200      CONTINUE
        ENDIF
      ENDIF
 
C      DO 300 I=1,N
C        WRITE(STDOUT,*) ' I=',I,' XBRK(I)=',XBRK(I)
C300   CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE READ_UFGATE_ITEMS(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     M                            EFLAG, IDLEN,  
     O                            HG, TAB2D, CCVAL, ANGLE)

C     Get the items of data from input line in UFGATE

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A         EFLAG, IDLEN, TAB2D
      REAL HG, CCVAL, ANGLE
      CHARACTER LINE*(*)

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:759* Only ',I3,' items given in ',
     A   'UFGATE description line.  Need four items.')

C***********************************************************************
      IF(NITEM.LT.4) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

      N = 1
C     Process the gate opening
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') HG

C     Process the table id
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      LKEY = LENSTR(KEY)
      IDLEN = MAX(IDLEN, LKEY)
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 ITAB)
        TAB2D = ITAB
      ELSE
        TAB2D = 0
      ENDIF

C     Process the the contraction coefficient
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') CCVAL

C     Process the angle
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') ANGLE

      RETURN
      END
C
C
C
      SUBROUTINE   UFGATE
     I                   (GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR)
 
C     + + + PURPOSE + + +
C     Compute the function tables needed to define the flow
C     for an underflow gate.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL GRAV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TABDIR - Table directory to remember table numbers
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'ufgate.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER MAXN
      PARAMETER (MAXN=4)
      INTEGER CCTAB, EF, ELFLAG, I, IFLAG, IHG, IHGERR, IHU, IPFD, J,
     A        JBASE, N, NFRAC, NHG, NHU, NRMS, TAB, TABGT, TABLT, TPTAB,
     B        COLWIDTH, NITEM, IDLEN
      INTEGER TAB2D(PMXNHG), ITEM_START(MAXN), ITEM_END(MAXN)
      REAL ANGLE(PMXNHG), BIGERR, BRKPFD, CCFOLL, CCVAL(PMXNHG), CONCC,
     A     DBETA4, DE1TO4, DIV, DROP, DT, FAC, FDROP, FDVEC(PMXNHU),
     B     FINPOW, FWFOTR, H1, H1FOLL, H1FWUL, H1SWSO, H3, H4, H4F,
     C     H4FOLL, H4FWUL, H4SWSO, HDATUM, HERR, HGVEC(PMXNHG),
     D     HSTUFF(PMXNHG,3), HUMAX, HUMIN, HUVEC(PMXNHU), LIMPFD,
     E     LIPREC, LIPVEC(PMXNHG), MAXHU, MINHU, MINPFD, MNHVEC(PMXNHG),
     F     OFFSET, OLDHG, OLDY3, P, PFDTMP(PMXFRC), PFDVEC(PMXFRC), POW,
     G     POWER, PREC, Q, QFO, QFREE, QHAT, QMAT(PMXNHU,PMXFRC), RERR,
     H     RG, RHS, RMS, T4, WORK(PMXFRC), XBRK(PMXFRC), Y, Y2FO, Y2FW,
     I     Y2OLD, Y4F, Y4SW, Y4SWSO, YLEFT, YRIGHT, Z1FWUL, ZW4, ZW4F,
     J     TP, RGFWUL, zrhufd

      real*8 easting, northing
      CHARACTER CHAR6*6, FTYPE*2, HGLAB*10, LAB2*50, LABEL*50, LINE*80,
     A          CQ*10, KEY*16, TABID*16, APPTABID*16, DEPTABID*16,
     B          CCTABID*16, JUST*5, IDOUT*32, ID*16,
     c   zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, FLOAT, LOG, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LENSTR
      REAL FINDCC, RSOMY3, RSOMY4, RSWMY3, RSWMY4
      CHARACTER PUT10*10, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKTAB, FDROOT, FINDCC, FNDELV, FNDFOQ, FOTOSO, inline,
     A         INVTSE, LKTJ, LSTOPF, REGFAL, RSOMY3, RSOMY4, RSWMY3, 
     B         RSWMY4, TABCHK, TWDOUT, XLKT22, PUT10,
     C         GET_INTERNAL_TAB_NUMBER, READ_TABID, LENSTR,
     D         GET_TABID
 
C     + + + INPUT FORMATS + + +
C 1    FORMAT(7X,I5)
 2    FORMAT(6X,A)
 4    FORMAT(7X,I5)
 6    FORMAT(7X,I5)
 8    FORMAT(9X,F10.0)
 10   FORMAT(10X,F10.0)
 12   FORMAT(3X,F10.0)
 14   FORMAT(6X,I5)
 16   FORMAT(7X,F10.0)
 18   FORMAT(7X,F10.0)
 19   FORMAT(9X,F10.0)
 20   FORMAT(8X,F10.0)
 21   FORMAT(8X,F10.0)
 22   FORMAT(10X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' TabId= ',A,' for type 15 table for underflow',
     A     ' gate.')
 52   FORMAT(/,' Label=',A50)
 54   FORMAT(/,' Approach section table id= ',A)
 56   FORMAT(/,' Departure section table id= ',A)
 58   FORMAT(/,' Elevation of gate sill=',F10.3)
 60   FORMAT(/,' Total gate opening width=',F10.3)
 62   FORMAT(/,' Discharge coefficient for approach to gate=',F10.3)
 64   FORMAT(/,' Contraction coefficient table id= ',A)
 65   FORMAT(/,' *ERR:774* Must have at least two gate openings.')
 66   FORMAT(/,' Minimum partial free drop=',F10.5)
 67   FORMAT(/,' Partial free drop at power breakpoint=',F10.5)
 68   FORMAT(/,' Limiting partial free drop=',F10.3)
 69   FORMAT(/,' Final local power for partial free drops=',F10.3)
 70   FORMAT(/,'Two-D table computations for gate opening=',F8.3,
     A    '  Two-D table id= ',A)
71    FORMAT(5X,'Head at section 1 for FW/FO boundary=', F9.4,/,
     A       5X,'Head-relative gate opening at FW/FO boundary=',F9.4)
 72   FORMAT(/,'Upstream head=',F9.4,
     A  ' Elevation=',F10.4,/,2X,'Depth at section 1=',F8.4,
     B  ' Gate opening=',F8.4)
 74   FORMAT(/,
     A '  Partial  Drop    Head    Head   Flow Cont.  Discharge Local',
     B '  Energy',/,
     C '   free    sect.   sect.   sect.  type coef.            power',
     D '   loss',/,
     E '   drop    1->4     3       4           Cc  ',17X,'   1->4',/,
     F ' --------  ------  ------  ------  ---  ---- --------- -----',
     G '  ------')
 75   FORMAT(1X,F8.4,F8.3,F8.3,F8.3,3X,A2,A6,1X,A10,F6.2,F8.3)
 76   FORMAT(1X,F8.4,F8.3,F8.3,F8.3,3X,A2,A6,1X,A10,6X,F8.3)
 77   FORMAT(/,' *WRN:592* Please review results.  One or more',
     A        ' cases with energy gain found.')
 79   FORMAT(/,' Free weir to free orifice transition fraction=',F8.2)
 80   FORMAT(/,' Maximum upstream head=',F8.2)
 81   FORMAT(/,' Minimum non-zero upstream head=',F8.2)
 82   FORMAT(/,' Linear interpolation precision=',F8.3)
 86   FORMAT(/,' Maximum relative error=',F6.3,' Gate opening=',
     A     F8.4,/,'   Upstream head=',F9.4,' Partial free drop=',
     B     F8.5)
 87   FORMAT(/,' *ERR:607* TabId missing or negative.')
 88   FORMAT(/' Root-mean-squared error=',F6.3,' N in sample=',I5)
 89   FORMAT('  Processing UFGATE TabId= ',A)
 90   FORMAT(/,' *ERR:712* Gate sill elevation=',F10.3,' <  elevation',
     A       ' of floor of',/,11X,' departure reach=',F10.3)
 91   FORMAT(/,' *ERR:713* Gate opening width=',F10.3,' <= 0.0')
 92   FORMAT(/,' *ERR:714* Approach loss Cd <=0.0 or > 1.0')
 93   FORMAT(/,' *ERR:715* Contraction coef.=',F10.3,' <= 0.0 or > 1.0')
 94   FORMAT(/,' *ERR:716* In UFGATE: FLAG=',I3,' No solution',A)
 95   FORMAT('TABID= ',A,/,'TYPE=  -15',/,'REFL=',7X,'0.0 LABEL=',A50,
     B /,'   OPENING            TABID   H1FWULR   H4FWULR H4SWSOMDR')
 96   FORMAT(F10.3,1X,A16,3F10.6)
 97   FORMAT(/,' *BUG:XXX* In UFGATE no root in 16 tries for:',A)
 98   FORMAT(/,' *ERR:717* Gate sill elevation=',F10.3,' <  elevation',
     A       ' of floor of',/,11X,' approach reach=',F10.3)
 99   FORMAT(/,' *ERR:* Interpolation precision tables missing.',
     A   ' Check for version',/,11X,'of file: TYPE5.TAB.')
C***********************************************************************
C     Define the linear interpolation precision tables.
      KEY = '10001'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, KEY,
     M                             EFLAG,
     O                             TABLT)
      TABLT = FTPNT(TABLT)
      KEY = '10002'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, KEY,
     M                             EFLAG,
     O                             TABGT)
      TABGT = FTPNT(TABGT)
      IF(TABLT.LT.1.OR.TABGT.LT.1) THEN
        WRITE(STDOUT,99) 
        EFLAG = 1
        RETURN
      ENDIF
      G = GRAV
      TWOG = 2.*GRAV
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'TAB',
     O                EFLAG, TABID, TAB)
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID))
      WRITE(*,89) TABID(1:LENSTR(TABID))
 
      CALL TABCHK
     I           (STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)
 
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) LAB2
      WRITE(STDOUT,52) LAB2
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'APP',
     O                EFLAG, APPTABID, APPTAB)
      WRITE(STDOUT,54) APPTABID(1:LENSTR(APPTABID))
 
      IF(APPTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        EF = 0
        TPTAB = APPTAB
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              APPTAB,
     O              EF)
        IF(EF.EQ.0) THEN
C         Find the invert elevation from the table
          CALL FNDELV
     I               (TPTAB, STDOUT,
     O                EFLAG, Z1B)
        ELSE
          EFLAG = 1
        ENDIF
      ENDIF
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'DEP',
     O                EFLAG, DEPTABID, DEPTAB)
      WRITE(STDOUT,56) DEPTABID(1:LENSTR(DEPTABID))
 
      IF(DEPTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        EF = 0
        TPTAB = DEPTAB
        CALL CHKTAB
     I             (1, STDOUT, FTPNT, MFTNUM,
     M              DEPTAB,
     O              EF)
        IF(EF.EQ.0) THEN
C         Find the invert elevation from the table
          CALL FNDELV
     I               (TPTAB, STDOUT,
     O                EFLAG, Z4B)
        ELSE
          EFLAG = 1
        ENDIF
      ENDIF
 
C     For the moment make Z3B the same as Z4B.  Note that
C     the table called the departure table will give the
C     cross section that defines the upstream hydrostatic pressure
C     force on the control volume downstream of the underflow
C     gate.  An optional table, to be added later will give
C     the table for computing the momentum flux exiting from the
C     control volume.
 
      Z3B = Z4B
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,8,ERR=991) Z2B
 
      WRITE(STDOUT,58)  Z2B
 
      IF(Z2B.LT.Z4B) THEN
        WRITE(STDOUT,90) Z2B, Z4B
        EFLAG = 1
      ENDIF
      IF(Z2B.LT.Z1B) THEN
        WRITE(STDOUT,98) Z2B, Z1B
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,10,ERR=991)  BG
      WRITE(STDOUT,60) BG
 
      IF(BG.LE.0.0)  THEN
        WRITE(STDOUT,91) BG
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,12,ERR=991) CD
      WRITE(STDOUT,62)  CD
      IF(CD.LE.0.0.OR.CD.GT.1.0) THEN
        WRITE(STDOUT,92) CD
        EFLAG = 1
      ENDIF
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'CCTAB',
     O                EFLAG, CCTABID, CCTAB)
      WRITE(STDOUT,64) CCTABID(1:LENSTR(CCTABID))
      IF(CCTAB.LT.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSEIF(CCTAB.GT.0) THEN
        CALL CHKTAB
     I             (2, STDOUT, FTPNT, MFTNUM,
     M              CCTAB,
     O              EFLAG)
      ENDIF
 
C     Get the size of the transition between FW and FO flow.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,19,ERR=991)  FWFOTR
      IF(FWFOTR.EQ.0.0) FWFOTR = 0.1
      WRITE(STDOUT,79) FWFOTR
 
C     Get the maximum upstream head
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,20,ERR=991)  MAXHU
      WRITE(STDOUT,80) MAXHU
 
C     Get the smallest non-zero upstream head
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,21,ERR=991)  MINHU
      IF(MINHU.EQ.0.0) MINHU = 0.1
      WRITE(STDOUT,81) MINHU
 
C     Get the global linear interpolation precisions
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,22,ERR=991)  LIPREC
      IF(LIPREC.EQ.0.0) LIPREC = 0.02
      WRITE(STDOUT,82) LIPREC
 
C     Process the gate opening table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,'(A)') LINE
      JUST = 'RIGHT'
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)
C     Set the  column width
      COLWIDTH =  ITEM_END(2) - ITEM_START(2) + 1
 
 
      OLDHG = 0.0
      NHG = 1
      IDLEN = 0
 100  CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        CALL READ_UFGATE_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, IDLEN,  
     O            HGVEC(NHG), TAB2D(NHG), CCVAL(NHG), ANGLE(NHG))
C        READ(LINE,'(F10.0, I10, 4F10.0)') HGVEC(NHG), TAB2D(NHG),
C     A                      CCVAL(NHG), ANGLE(NHG)
C     B                      MNHVEC(NHG), LIPVEC(NHG)
        MNHVEC(NHG) = 0.0
        LIPVEC(NHG) = 0.0
        IF(HGVEC(NHG).LE.OLDHG) THEN
C         Input complete.
          NHG = NHG - 1
        ELSE
          OLDHG = HGVEC(NHG)
          IF(MNHVEC(NHG).EQ.0.0) THEN
            MNHVEC(NHG) = MINHU
          ENDIF
          IF(LIPVEC(NHG).EQ.0.0) THEN
            LIPVEC(NHG) = LIPREC
          ENDIF
          ID = GET_TABID(TAB2D(NHG))
          IDOUT = ' '
          IDOUT(COLWIDTH-IDLEN+1:COLWIDTH) = ID(1:IDLEN)
          WRITE(STDOUT,
     A        '(F10.3,A,F10.3,F10.1,F10.2,F10.3)') HGVEC(NHG),
     B           IDOUT(1:COLWIDTH), CCVAL(NHG), ANGLE(NHG)
C     C           MNHVEC(NHG), LIPVEC(NHG)
          CALL TABCHK
     I               (STDOUT, PMXTAB,
     M                TAB2D(NHG), TABDIR, EFLAG)
          IF(CCTAB.EQ.0.0) THEN
C           If no table is given for contraction coefficient, then
C           a value must be given by the user.
            IF(CCVAL(NHG).LE.0.0.OR.CCVAL(NHG).GT.1.0) THEN
              WRITE(STDOUT,93) CCVAL(NHG)
               EFLAG = 1
            ENDIF
          ENDIF
          NHG = NHG + 1
          GOTO 100
        ENDIF
 
      IF(NHG.LT.2) THEN
        WRITE(STDOUT,65)
        EFLAG = 1
      ENDIF
 
C     INPUT THE FACTORS CONTROLLING THE DISTRIBUTION OF DROPS
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,'(A)') LINE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) MINPFD
      WRITE(STDOUT,66) MINPFD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) BRKPFD
      WRITE(STDOUT,67) BRKPFD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,18,ERR=991) LIMPFD
      WRITE(STDOUT,68) LIMPFD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) FINPOW
      WRITE(STDOUT,69) FINPOW
 
 
C     COMPUTE THE PROPORTIONS OF FREE DROP
 
      POW = 0.5
      OFFSET = 0.0
      CALL LSTOPF 
     I           (STDOUT, TABLT, TABGT, POW, OFFSET, MINPFD, BRKPFD,
     I            LIPREC, PMXFRC,
     O            N, XBRK, EFLAG)
 
      WORK(1) = 0.0
      DO 195 I=1,N
        WORK(I+1) = XBRK(I)
 195  CONTINUE
      NFRAC = N + 1
C      OFFSET = 0.9*BRKPFD
      OFFSET = 0.0
      CALL LSTOPF 
     I           (STDOUT, TABLT, TABGT, FINPOW, OFFSET, BRKPFD, LIMPFD,
     I            LIPREC, PMXFRC,
     O            N, XBRK, EFLAG)
      DO 196 I=2,N
        WORK(NFRAC+I-1) = XBRK(I)
 196  CONTINUE
      NFRAC = NFRAC + N - 1
 
C      WORK(N+2) = 0.95
C      WORK(N+3) = 0.975
C      WORK(N+4) = 0.98
C      WORK(N+5) = 0.99
      WORK(NFRAC+1) = 1.0
      NFRAC = NFRAC + 1
 
      PFDVEC(1) = WORK(1)
      PFDVEC(2) = WORK(2)
      J = 2
      DO 201 I=3,NFRAC-1
        J = J + 1
        PFDVEC(J) = 0.5*(WORK(I) + WORK(I-1))
        J = J + 1
        PFDVEC(J) = WORK(I)
 201  CONTINUE
      PFDVEC(J+1) = WORK(NFRAC)
      NFRAC = J + 1
C      DO 202 I=1,NFRAC
C        WRITE(STDOUT,*) ' I=',I,' PFDVEC(I)=',PFDVEC(I)
C202   CONTINUE
 
      IF(EFLAG.NE.0) RETURN
 
C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.
 
      HDATUM = Z2B
 
      DZ = Z3B - Z2B
 
      BIGERR = 0.0
      RMS = 0.0
      NRMS = 0
      OLDHG = 0.0
      DO 1000 IHG=1,NHG
        ELFLAG = 0
        HG = HGVEC(IHG)
        WRITE(HGLAB,'(''Hg='',F6.3,'':'')') HG
        LABEL = HGLAB // LAB2
        MINHU = MNHVEC(IHG)
        PREC = LIPVEC(IHG)
        CONCC = CCVAL(IHG)
        ID = GET_TABID(TAB2D(IHG))
        WRITE(STDOUT,70) HG, ID(1:LENSTR(ID))
 
        AG = BG*HG
 
C       Find the free weir flow(FW flow) that exists when critical depth
C       in the gate opening equals the gate opening.
        Y2 = HG
        Q = BG*Y2*SQRT(G*Y2)
        RHS = Z2B + Y2 + (Q/(CD*AG))**2/TWOG - Z1B
C       Make first rough  estimate of the depth at section
C       1.  The final water surface elevation at section 1
C       should be higher than the water surface elevation at section 2
C       at critical flow.
        Y = Y2 + Z2B - Z1B
 
        CALL INVTSE
     I             (G, STDOUT, APPTAB, 22, Q, RHS, Y,
     O              EFLAG)
C        Y1FWUL = Y
        Z1FWUL = Y + Z1B
        H1FWUL = Z1FWUL - HDATUM
        HSTUFF(IHG,1) = H1FWUL/HG
C       Compute the gate-opening ratio for the free-weir flow limit.
        RGFWUL = HG/H1FWUL
 
        WRITE(STDOUT,71) H1FWUL, RGFWUL

C        WRITE(STDOUT,*) ' FW/FO boundary at section 1:'
C        WRITE(STDOUT,*) ' Y1FWUL=',Y1FWUL, ' H1FWUL=', H1FWUL,
C     A                    ' Gate ratio=',RGFWUL
C        WRITE(STDOUT,*) ' H1FWUL/HG=',H1FWUL/HG
 
 
C       Find the non-standard FO flow at the upper limit of
C       FW flow.  This flow should be the same as the flow computed
C       assuming critical flow.  That is why the FO flow is called
C       non-standard.
 
        CC = 1.0
        CALL FNDFOQ
     I             (H1FWUL, HDATUM,
     O              QFO)
 
        WRITE(STDOUT,*) ' FW flow at its upper limit=',Q
        WRITE(STDOUT,*) ' Non-standard FO flow=',QFO
 
C       Compute the water level at section 4 that just initiates
C       submergence of the non-standard FO flow.
        CALL FOTOSO
     I             (H1FWUL, QFO, HDATUM,
     O              H4FWUL, IFLAG)
C        WRITE(STDOUT,*) ' Tailwater head for submergence of ',
C     A               'non-standard FO flow: H4FWUL=',H4FWUL
C        WRITE(STDOUT,*) ' H4FWUL/HG=',H4FWUL/HG
        HSTUFF(IHG,2) = H4FWUL/HG
        IF(IFLAG.NE.0) THEN
          WRITE(STDOUT,94) IFLAG,
     A               ' for nonstandard FO flow tailwater limit.'
          EFLAG = 1
        ENDIF
C       Compute the assumed lower limit of  standard
C       FO flow, that is, FO flow with a contraction coefficient
C       as defined for FO flow only.
        H1FOLL = H1FWUL + FWFOTR*HG
 
C       Find the contraction coefficient at lower limit for FO flow
        IF(ANGLE(IHG).EQ.0.0) THEN
          CCFOLL = FINDCC(HG/H1FOLL, CONCC, CCTAB)
        ELSE
          CCFOLL = FINDCC(ANGLE(IHG), CONCC, CCTAB)
        ENDIF
 
C       Find the tailwater level that just initiates submergence of
C       the FO flow at its lower limit
 
        CC = CCFOLL
        CALL FNDFOQ
     I             (H1FOLL, HDATUM,
     O              QFREE)
 
C       Now find the level at section 4 that defines
C       the limit of free orifice flow.
        CALL FOTOSO
     I             (H1FOLL, QFREE, HDATUM,
     O              H4FOLL, IFLAG)
        IF(IFLAG.NE.0) THEN
          WRITE(STDOUT,94) IFLAG,
     A      ' for lower-limit FO flow tailwater limit.'
          EFLAG = 1
        ENDIF
 
C        WRITE(STDOUT,*) ' Lower limit of standard FO flows:'
C        WRITE(STDOUT,*) ' H1FOLL=',H1FOLL,' H4FOLL=',H4FOLL
C        WRITE(STDOUT,*) ' Standard FO flow=',QFREE
 
C       Find the level at section 4  that is at the boundary
C       between SW and SO when section 1 is midway between
C       head = HG and head = H1FWUL.  Use this level to provide for
C       non-linearity of the boundary between SW and SO.
 
        H1SWSO = 0.5*(HG + H1FWUL)
        TP = H1SWSO+HDATUM-Z1B
        CALL XLKT22
     I             (APPTAB,
     M              TP,
     O              A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1,
     O              DALPHA, QCT)
        QSQR = TWOG*(AG*CD)**2*(H1SWSO + HDATUM - HG - Z2B)/
     A                 (1.0 - ALPHA1*(CD*AG/A1)**2)
 
C       Establish the constant values in RSWMY4.
        Y3 = HG + Z2B - Z3B
        CALL LKTJ
     I           (DEPTAB,
     M            Y3,
     O            J4ATY3)
C       Estimate initial values to bracket the root.
        YLEFT = Z2B + HG - Z4B
        YRIGHT = H1SWSO + HDATUM  - Z4B
        CALL REGFAL
     I             (EPSARG, EPSF, RSWMY4,
     M              YLEFT, YRIGHT,
     O              Y4SWSO, IFLAG)
        IF(IFLAG.NE.0) THEN
          WRITE(STDOUT,94) IFLAG,' for SW to SO boundary at midpoint.'
          EFLAG = 1
        ENDIF
        H4SWSO = Y4SWSO + Z4B - HDATUM
        HSTUFF(IHG,3) = H4SWSO/HG
 
C       Assign the upstream heads for this gate opening.  The
C       ranges are: MINHU to HG, HG to H1FWUL, H1FWUL to H1FOLL,
C       and H1FOLL to MAXHU
        POWER = 1.5
        OFFSET = 0.0
        IF(MINHU.GE.HG) THEN
          MINHU = 0.5*HG
        ENDIF
        HUMIN = MINHU
        CALL LSTOPF 
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, HUMIN, HG,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        DO 210 I=1,N
          HUVEC(I) = XBRK(I)
 210    CONTINUE
        NHU = N
 
        CALL LSTOPF 
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, HG, H1FWUL,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        IF(N.EQ.2) THEN
          XBRK(3) = XBRK(2)
          XBRK(2) = 0.5*(XBRK(1) + XBRK(3))
          N = 3
        ENDIF
        DO 212 I=2,N
          NHU = NHU + 1
          HUVEC(NHU) = XBRK(I)
 212    CONTINUE
        OFFSET = -CC*HG
C       Flow varies as the .5 power of the head - offset.  However,
C       the drop to free flow varies more nearly like the 1.5 power.
C       The spacing for the 1.5 power is smaller than for .5 power;
C       so we use 1.5 power.
        POWER = 1.5
        CALL LSTOPF 
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, H1FWUL, H1FOLL,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        IF(N.EQ.2) THEN
          XBRK(3) = XBRK(2)
          XBRK(2) = 0.5*(XBRK(1) + XBRK(3))
          N = 3
        ENDIF
        DO 213 I=2,N
          NHU = NHU + 1
          HUVEC(NHU) = XBRK(I)
 213    CONTINUE
        IF(OLDHG.GT.0.0) THEN
          HUMAX = MAXHU*HG/OLDHG
        ELSE
          HUMAX = MAXHU
        ENDIF
        IF(HUMAX.LE.H1FOLL) THEN
C         The maximum requested head is less than the head at section
C         1 at the lower limit of FO flow.
          HUMAX = H1FOLL*1.10
        ENDIF
        OLDHG = HG
        CALL LSTOPF
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, H1FOLL, HUMAX,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        DO 214 I=2,N
          NHU = NHU + 1
          HUVEC(NHU) = XBRK(I)
 214    CONTINUE
 
C        WRITE(STDOUT,*) ' NHU=',NHU
 
C       One-shot modification for computing the Elmhurst Q. gate
C       values.

C        IF(IHG.EQ.1) THEN
CC         One head value.
C          NHU = 3
C          HUVEC(1) = 9.52
C          HUVEC(2) = 6.17
C          HUVEC(3) = 8.47
C        ELSEIF(IHG.EQ.2) THEN
CC         One head value
C          NHU = 4
C          HUVEC(1) = 9.40
C          HUVEC(2) = 6.08
C          HUVEC(3) = 8.33
C          HUVEC(4) = 7.5
C        ENDIF
        
C       Now compute the free flow for each of the upstream heads
        DO 900 IHU=1,NHU
          H1 = HUVEC(IHU) + HDATUM - Z2B
          Y1 = H1 + Z2B - Z1B
          WRITE(STDOUT,72) HUVEC(IHU), H1 + Z2B, Y1, HG
          WRITE(STDOUT,74)
C         Set the flow for zero partial free drop to 0.0
          QMAT(IHU,1) = 0.0
C         Find conditions at section 1
          CALL XLKT22
     I               (APPTAB,
     M                Y1,
     O                A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1,
     O                DALPHA, QCT)
          IF(H1.LE.H1FWUL) THEN
C           Flow is free weir flow.   Try linear iteration.
            FTYPE = 'FW'
            FAC = 0.0
            Y2OLD = 0.0
            DIV = 1.0 + 0.5/CD**2
 320        CONTINUE
              Y2 = H1/(DIV - FAC)
              IF(ABS(Y2 - Y2OLD)/Y2.GT.EPSARG) THEN
C               No convergence.  Compute next value of FAC and
C               try again.
                Y2OLD = Y2
                FAC = 0.5*ALPHA1*(Y2*BG/A1)**2
                GOTO 320
              ENDIF
            QFREE = Y2*BG*SQRT(G*Y2)
            Y2FW = Y2
C            WRITE(STDOUT,*) ' H1=',H1,' FW Q=', QFREE,' Y2FW=',Y2FW
 
C           Submerged weir flow can transition to submerged orifice
C           flow if the water level at section 1 is above the
C           elevation of the gate lip.
            IF(Y1 + Z1B.GT.HG+Z2B) THEN
C             Find the tailwater level that causes the water to
C             contact the lip.  That is, when Y3 + Z3B = HG + Z2B, thus
C             Y3 = HG + Z2B - Z3B.  Find square of the flow at this
C             condition.
 
              QSQR = TWOG*(AG*CD)**2*(Z1B + Y1 - HG - Z2B)/
     A                       (1.0 - ALPHA1*(CD*AG/A1)**2)
 
C             Establish the constant values in RSWMY4.
              Y3 = HG + Z2B - Z3B
              CALL LKTJ
     I                 (DEPTAB,
     M                  Y3,
     O                  J4ATY3)
 
C             The water-surface elevation at section 4 will be
C             between the water surface elevation at section 1
C             and the elevation of the gate lip.
              YLEFT = Z2B + HG - Z4B
C              FLEFT = RSWMY4(YLEFT)
              YRIGHT = Z1B + Y1 - Z4B
C              FRIGHT = RSWMY4(YRIGHT)
C              WRITE(STDOUT,*) ' YLEFT=',YLEFT,' FLEFT=',FLEFT
C              WRITE(STDOUT,*) ' YRIGHT=',YRIGHT,' FRIGHT=',FRIGHT
              CALL REGFAL
     I                   (EPSARG, EPSF, RSWMY4,
     M                    YLEFT, YRIGHT,
     O                    Y4SW, IFLAG)
              IF(IFLAG.NE.0) THEN
                WRITE(STDOUT,94) IFLAG,' for SW to SO boundary.'
                EFLAG = 1
              ENDIF
 
C              RGSW = HG/H1
 
C              WRITE(STDOUT,*) ' SW/SO boundary at section 4:'
C              WRITE(STDOUT,*) ' Y4SW=',Y4SW,' H4SW=',Y4SW + Z4B - Z2B,
C     A               ' Gate opening ratio=',RGSW
            ELSE
              Y4SW = Y1 + Z1B - Z4B
            ENDIF
 
C           Now find the level at section 4 that defines
C           the limit of free weir flow, Y4F. Define the constant
C           values.
            Y3 = Y2FW + Z2B - Z3B
            OLDY3 = Y3
            CALL LKTJ
     I               (DEPTAB,
     M                Y3,
     O                J4ATY3)
            QSQR = QFREE**2
 
C           The water surface elevation at section 4 will be
C           below the elevation at section 1 but above the
C           critical elevation at section 2.
            YLEFT = Y2FW + Z2B - Z4B
            YRIGHT = Y1 + Z1B - Z4B
            CALL REGFAL
     I                 (EPSARG, EPSF, RSWMY4,
     M                  YLEFT, YRIGHT,
     O                  Y4F, IFLAG)
            IF(IFLAG.GT.0) THEN
              WRITE(STDOUT,94) IFLAG,' for FW flow tail water limit.'
              EFLAG = 1
            ENDIF
            ZW4F = Y4F + Z4B
            H4F = ZW4F - HDATUM
            FDROP = Y1 + Z1B - ZW4F
 
C            WRITE(STDOUT,*) ' Tailwater at free weir flow limit:'
C            WRITE(STDOUT,*) ' Y4F=',Y4F,' H4F=',H4F
C            WRITE(STDOUT,*) ' Submergence ratio on total head:',
C     A             H4F/(H1 + ALPHA1*(QFREE/A1)**2/TWOG)
 
C           Compute a check on the energy balance from section 1 to 4.
            DE1TO4 = H1 + ALPHA1*(QFREE/A1)**2/TWOG -
     A           (H4F + ALPHA4*(QFREE/A4)**2/TWOG)
            IF(DE1TO4.LT.0.0) ELFLAG = 1
            QMAT(IHU,NFRAC) = QFREE
            FDVEC(IHU) = FDROP
            CHAR6 = '  --- '
            CQ = PUT10(QFREE)
            WRITE(STDOUT,76) 1.00, FDROP,
     A          Y3 + Z3B - HDATUM, H4F, FTYPE, CHAR6, CQ, DE1TO4
C            OLDH3 = Y3 + Z3B - HDATUM
C           Compute the submerged flows.  They may be of two types:
C           SO or SW.
            DO 400 J=NFRAC-1,2,-1
              DROP = FDROP*PFDVEC(J)
              ZW4 = Z1B + Y1 - DROP
              Y4 = ZW4 - Z4B
              H4 = ZW4 - HDATUM
              YLEFT = OLDY3
              YRIGHT = Y4
              CALL XLKT22
     I                   (DEPTAB,
     M                    Y4,
     O                    A4, TT, DT, J4, KT, DKT, BETA4, DBETA, ALPHA4,
     O                    DALPHA, QCT)
              IF(Y1 + Z1B.LE.HG + Z2B) THEN
C               Flow can only be submerged weir flow.
                FTYPE = 'SW'
                CHAR6 = '  --- '
                CALL FDROOT
     I                     (YRIGHT, RSWMY3, EPSF,
     M                      YLEFT,
     O                      IFLAG)
                IF(IFLAG.EQ.1) THEN
                  WRITE(STDOUT, 97)  ' SW flow'
                  EFLAG = 1
                ENDIF
                CALL REGFAL
     I                     (EPSARG, EPSF, RSWMY3,
     M                      YLEFT, YRIGHT,
     O                      Y3, IFLAG)
                IF(IFLAG.NE.0) THEN
                  WRITE(STDOUT,94) IFLAG,' for SW flow.'
                  EFLAG = 1
                ENDIF
                Q = SQRT(QSQR)
C                WRITE(STDOUT,*) ' SW Q=',Q,' Y3=',Y3,' IFLAG=',IFLAG
 
              ELSE
                IF(Y4.GT.Y4SW) THEN
C                 Flow is submerged orifice.
                  RG = HG/H1
                  CC = 1.0
 
                  AT = CD*CC*AG
                  FTYPE = 'SO'
                  WRITE(CHAR6,'(F6.3)') CC
C                  FLEFT = RSOMY3(YLEFT)
C                  WRITE(STDOUT,*) ' YLEFT=',YLEFT,' FLEFT=',FLEFT
C                  FRIGHT = RSOMY3(YRIGHT)
C                  WRITE(STDOUT,*) ' YRIGHT=',YRIGHT,' FRIGHT=',FRIGHT
                  CALL FDROOT
     I                       (YRIGHT, RSOMY3, EPSF,
     M                        YLEFT,
     O                        IFLAG)
                  IF(IFLAG.EQ.1) THEN
                    WRITE(STDOUT, 97)  ' SO flow after SW flow'
                    EFLAG = 1
                  ENDIF
                  CALL REGFAL
     I                       (EPSARG, EPSF, RSOMY3,
     M                        YLEFT, YRIGHT,
     O                        Y3, IFLAG)
                  IF(IFLAG.NE.0) THEN
                   WRITE(STDOUT,94) IFLAG,' for SO flow after SW flow.'
                    EFLAG = 1
                  ENDIF
                  Q = SQRT(QSQR)
                ELSE
C                 Flow is submerged weir.
                  FTYPE = 'SW'
                  CHAR6 = '  --- '
                  CALL FDROOT
     I                       (YRIGHT, RSWMY3, EPSF,
     M                        YLEFT,
     O                        IFLAG)
                  IF(IFLAG.EQ.1) THEN
                    WRITE(STDOUT, 97)  ' SW flow'
                    EFLAG = 1
                  ENDIF
                  CALL REGFAL
     I                       (EPSARG, EPSF, RSWMY3,
     M                        YLEFT, YRIGHT,
     O                        Y3, IFLAG)
                  IF(IFLAG.NE.0) THEN
                    WRITE(STDOUT,94) IFLAG,' for SW flow.'
                    EFLAG = 1
                  ENDIF
                  Q = SQRT(QSQR)
C                  WRITE(STDOUT,*) ' SW Q=',Q,' Y3=',Y3,' IFLAG=',IFLAG
                ENDIF
              ENDIF
              OLDY3 = Y3
              QMAT(IHU,J) = Q
 
              DE1TO4 = H1 + ALPHA1*(Q/A1)**2/TWOG -
     A           (H4 + ALPHA4*(Q/A4)**2/TWOG)
              IF(DE1TO4.LT.0.0) ELFLAG = 1
 
              H3 = Y3 + Z3B - HDATUM
C              P2 = Q/(QFREE*SQRT(PFDVEC(J)))
C              OLDH3 = H3
              P = LOG(Q/QMAT(IHU,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
              CQ = PUT10(Q)
              WRITE(STDOUT,75) PFDVEC(J), DROP,
     A          H3, H4, FTYPE, CHAR6, CQ, P, DE1TO4
 400        CONTINUE
          ELSE
C           Flow is free orifice flow.  Find the contraction
C           coefficient to use.  Find the gate opening ratio for
C           this upstream head
 
            FTYPE = 'FO'
            RG = HG/H1
 
            IF(H1.LT.H1FOLL) THEN
              CC = 1.0 + (H1 - H1FWUL)*(CCFOLL - 1.0)/(H1FOLL - H1FWUL)
            ELSE
              IF(ANGLE(IHG).LE.0.0) THEN
                CC = FINDCC(RG, CONCC, CCTAB)
              ELSE
                CC = FINDCC(ANGLE(IHG), CONCC, CCTAB)
              ENDIF
            ENDIF
C            WRITE(STDOUT,*) ' CC FOR FREE ORIFICE=',CC
            WRITE(CHAR6,'(F6.3)') CC
            AT = CD*CC*AG
            Y2 = CC*HG
            QFREE = AT*SQRT(TWOG*(Z1B + Y1 - Z2B - Y2)/
     A                    (1.0 - ALPHA1*(AT/A1)**2))
            Y2FO = Y2
C           Now find the level at section 4 that defines
C           the limit of free orifice flow, Y4F. Define the constant
C           values.
            Y3 = Y2FO + Z2B - Z3B
            OLDY3 = Y3
            CALL LKTJ
     I               (DEPTAB,
     M                Y3,
     O                J4ATY3)
            QSQR = QFREE**2
 
C           The water surface elevation at section 4 will be
C           below the elevation at section 1 but above the
C           jet surface elevation at section 2.
            YLEFT = Y2FO + Z2B - Z4B
            YRIGHT = Y1 + Z1B - Z4B
            CALL REGFAL
     I                 (EPSARG, EPSF, RSOMY4,
     M                  YLEFT, YRIGHT,
     O                  Y4F, IFLAG)
            IF(IFLAG.NE.0) THEN
              WRITE(STDOUT,94) IFLAG,' for FO tailwater limit.'
              EFLAG = 1
            ENDIF
            ZW4F = Y4F + Z4B
            H4F = ZW4F - HDATUM
            DE1TO4 = H1 + ALPHA1*(QFREE/A1)**2/TWOG -
     A           (H4F + ALPHA4*(QFREE/A4)**2/TWOG)
            IF(DE1TO4.LT.0.0) ELFLAG = 1
C            WRITE(STDOUT,*) ' Tailwater at free orifice flow limit:'
C            WRITE(STDOUT,*) ' Y4F=',Y4F,' H4F=',H4F
            QMAT(IHU,NFRAC) = QFREE
            FDROP = Y1 + Z1B - ZW4F
            FDVEC(IHU) = FDROP
            CQ = PUT10(QFREE)
            WRITE(STDOUT,76) 1.00, FDROP,
     A          Y3 + Z3B - HDATUM, H4F, FTYPE, CHAR6, CQ, DE1TO4
C            OLDH3 = Y3 + Z3B - HDATUM
            FTYPE = 'SO'
            AT = CD*CC*AG
            DO 500 J=NFRAC-1,2,-1
              DROP = FDROP*PFDVEC(J)
              ZW4 = Z1B + Y1 - DROP
              Y4 = ZW4 - Z4B
              H4 = ZW4 - HDATUM
              YLEFT =  CC*HG + Z2B - Z4B
              YRIGHT = Y4
              CALL XLKT22
     I                   (DEPTAB,
     M                    Y4,
     O                    A4, T4, DT, J4, KT, DKT, BETA4, DBETA4,
     O                    ALPHA4, DALPHA, QCT)
              CALL FDROOT
     I                   (YRIGHT, RSOMY3, EPSF,
     M                    YLEFT,
     O                    IFLAG)
              IF(IFLAG.EQ.1) THEN
                WRITE(STDOUT, 97)  ' SO flow'
                EFLAG = 1
              ENDIF
C             Flow is submerged orifice.
              CALL REGFAL
     I                   (EPSARG, EPSF, RSOMY3,
     M                    YLEFT, YRIGHT,
     O                    Y3, IFLAG)
              IF(IFLAG.NE.0) THEN
                WRITE(STDOUT,94) IFLAG,' for SO flow.'
                EFLAG = 1
              ENDIF
              Q = SQRT(QSQR)
              OLDY3 = Y3
              QMAT(IHU,J) = Q
              H3 = Y3 + Z3B - HDATUM
C              P2 = Q/(QFREE*SQRT(PFDVEC(J)))
C              OLDH3 = H3
              P = LOG(Q/QMAT(IHU,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
              DE1TO4 = H1 + ALPHA1*(Q/A1)**2/TWOG -
     A           (H4 + ALPHA4*(Q/A4)**2/TWOG)
              IF(DE1TO4.LT.0.0) ELFLAG = 1
 
C             Compute the rate of change of flow with respect
C             to change in the tailwater level.
C              CALL LKTA(DEPTAB, Y3, A4ATY3)
C              DQDY4 = (G*A4 + (DBETA4 - BETA4*T4/A4)*Q**2/A4)/
C     A                 (Q*(A4ATY3*(ALPHA1/A1**2 - 1.0/(CC*CD*AG)**2)
C     B                  + 2.*(1.0/(CC*AG) - BETA4/A4)))
 
C              DQ = (DQDY4 + 0.5*Q/DROP)/(QFREE*SQRT(PFDVEC(J)))
              CQ = PUT10(Q)
              WRITE(STDOUT,75) PFDVEC(J), DROP,
     A          H3, H4, FTYPE, CHAR6, CQ, P, DE1TO4
 500        CONTINUE
 
          ENDIF
 
C         Compute approximate maximum error and report
          DO 600 J=3,NFRAC-2,2
            QHAT = 0.5*(QMAT(IHU,J-1) + QMAT(IHU,J+1))
            RERR = ABS(QHAT - QMAT(IHU,J))/QMAT(IHU,J)
            RMS = RMS + RERR*RERR
            NRMS = NRMS + 1
            IF(RERR.GT.BIGERR) THEN
              BIGERR = RERR
              HERR = H1
              IPFD = J
              IHGERR = IHG
            ENDIF
 600      CONTINUE
 
C         Eliminate the checking values from QMAT
          JBASE = 3
          DO 700 J=4,NFRAC-1,2
            QMAT(IHU,JBASE) = QMAT(IHU,J)
            JBASE = JBASE + 1
 700      CONTINUE
          QMAT(IHU,JBASE) = QMAT(IHU,NFRAC)
 
 900    CONTINUE
        PFDTMP(1) = PFDVEC(1)
        PFDTMP(2) = PFDVEC(2)
        JBASE = 3
        DO 910 J=4,NFRAC-1,2
          PFDTMP(JBASE) = PFDVEC(J)
          JBASE = JBASE + 1
 910    CONTINUE
        PFDTMP(JBASE) = PFDVEC(NFRAC)
          
        zrhufd = 0.0
        CALL TWDOUT
     I             (STDOUT, STDTAB, TAB2D(IHG), LABEL, NHU, JBASE,
     I              HUVEC, FDVEC, PFDTMP, QMAT, HDATUM,
     I              13, '  UFGATE', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)
 
        IF(ELFLAG.GT.0) THEN
          WRITE(STDOUT,77)
        ENDIF
 1000 CONTINUE
      WRITE(STDOUT,86) BIGERR, HGVEC(IHGERR), HERR,
     A                 PFDVEC(IPFD)
 
      RMS = SQRT(RMS/FLOAT(NRMS))
      WRITE(STDOUT,88) RMS, NRMS
 
C     Output the table of type 15 to the table file
      WRITE(STDTAB,95) TABID(1:LENSTR(TABID)), LAB2
      DO 1005 I=1,NHG
        ID = GET_TABID(TAB2D(I))
        J = LENSTR(ID)
        TABID = ' '
        TABID(16-J+1:16)= ID(1:J)
        WRITE(STDTAB,96) HGVEC(I), TABID, (HSTUFF(I,J), J=1,3)
 1005 CONTINUE
      WRITE(STDTAB,96) -1.0
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
