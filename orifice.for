C
C
C
      SUBROUTINE READ_ORIFICE_ITEMS1(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     O                            NUMBER, SHAPE, EDGE, D, W, 
     O                            INVERT, CD, CW)

C     Get the items of data from the first line of ORIFICE data

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A        NUMBER
      REAL D, W, INVERT, CD, CW
      CHARACTER LINE*(*), SHAPE*6, EDGE*6

C     Local

      INTEGER IE, IS, N
      CHARACTER TPC*20

C     Called program units
      EXTERNAL STRIP_L_BLANKS
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:752* Only ',I3,' items given in ',
     A   'first orifice description line.  Need at least eight items.')
52    FORMAT(/,' *ERR:753* Conversion error in field ',I1,' in:',/,
     A        A)
C***********************************************************************
      IF(NITEM.LT.8) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.'
      ENDIF

      N = 1
C     Process NUMBER
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        NUMBER = 0
      ELSE
        READ(TPC,*,ERR=999) NUMBER
      ENDIF

C     Process the SHAPE description
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      SHAPE = TPC

C     Process the EDGE descriptions
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      EDGE = TPC
      
C     Process the vertical diameter
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        D = 0.0
      ELSE
        READ(TPC,*,ERR=999) D
      ENDIF
      

C     Process horizontal diameter
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        W = 0.0
      ELSE
        READ(TPC,*,ERR=999) W
      ENDIF
     

C     Process the invert elevation
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        INVERT = 0.0
      ELSE
        READ(TPC,*,ERR=999) INVERT
      ENDIF

C     Process discharge coefficient
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        CD = 0.0
      ELSE
        READ(TPC,*,ERR=999) CD
      ENDIF

C     Process the weir coefficient
      N = 8
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        CW = 0.0
      ELSE
        READ(TPC,*,ERR=999) CW
      ENDIF

      RETURN
999   CONTINUE
      WRITE(STDOUT,52) N, LINE
      STOP 'Abnormal stop.  Errors found.' 
      END

C
C
C
      SUBROUTINE READ_ORIFICE_ITEMS2(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     M                            EFLAG,
     O                            IDLENA, IDLENB, 
     O                            APPTAB, ORITAB, MAXZUP, MINHUP,
     O                            LIMPFD, MINPFD, LIPREC)

C     Get the items of data  from second line of ORIFICE data

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A        APPTAB, ORITAB, EFLAG, IDLENA, IDLENB
      REAL MAXZUP, MINHUP, LIMPFD, MINPFD, LIPREC
      CHARACTER LINE*(*)

C     Local

      INTEGER IE, IS, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER LENSTR
      EXTERNAL STRIP_L_BLANKS, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:754* Only ',I3,' items given in ',
     A   'second orifice input line.  Need at least seven items.')
52    FORMAT(/,' *ERR:753* Conversion error in field ',I1,' in:',/,
     A        A)
C***********************************************************************
      IF(NITEM.LT.7) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

C     Process the approach table id
      N = 1
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENA = LENSTR(KEY)
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
C       We have an id given.
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, KEY,
     M                               EFLAG,
     O                               APPTAB)
      ELSE
        APPTAB = 0   
      ENDIF
      

C     Process the orifice table id
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENB = LENSTR(KEY)
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
C       We have an id given.
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, KEY,
     M                               EFLAG,
     O                               ORITAB)
      ELSE
        ORITAB = 0   
      ENDIF

C     Process maximum elevation
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        MAXZUP = 0.0
      ELSE
        READ(TPC,*,ERR=999) MAXZUP
      ENDIF

      
C     Process the minimum upstream head
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        MINHUP = 0.0
      ELSE
        READ(TPC,*,ERR=999) MINHUP
      ENDIF
     

C     Process the limiting partial free drop
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC.EQ.' ') THEN
        LIMPFD = 0.0
      ELSE
        READ(TPC,*,ERR=999) LIMPFD
      ENDIF


C     Process the minimum partial free drop
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        MINPFD = 0.0
      ELSE
        READ(TPC,*,ERR=999) MINPFD
      ENDIF

C     Process linear interpolation precision
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        LIPREC = 0.0
      ELSE
        READ(TPC,*,ERR=999) LIPREC
      ENDIF


      RETURN
999   CONTINUE
      WRITE(STDOUT,52) N, LINE
      STOP 'Abnormal stop.  Errors found.'
      
      END
C     ***********
C     *         *
C     * FUN_CIRC
C     *         *
C     ***********

      DOUBLE PRECISION FUNCTION FUN_CIRC(Y)

C     Return the vertical diameter normalized width of the 
C     opening at the vertical diameter normalized distance 
C     from the invert. 

      DOUBLE PRECISION Y
C*********************************************************************
      FUN_CIRC = 2.D0*DSQRT(Y*(1.D0 - Y))
      RETURN
      END
C     ***********
C     *         *
C     * FUN_RECT
C     *         *
C     ***********

      DOUBLE PRECISION FUNCTION FUN_RECT(Y)

C     Return the vertical diameter normalized width of the 
C     opening at the vertical diameter normalized distance 
C     from the invert. 

      DOUBLE PRECISION Y

      INCLUDE 'orfshape.cmn'
C*********************************************************************
      FUN_RECT = W_OVER_D
      RETURN
      END
C     ***********
C     *         *
C     * FUN_TRI
C     *         *
C     ***********

      DOUBLE PRECISION FUNCTION FUN_TRI(Y)

C     Return the vertical diameter normalized width of the 
C     opening at the vertical diameter normalized distance 
C     from the invert. 

      DOUBLE PRECISION Y

      INCLUDE 'orfshape.cmn'
C*********************************************************************
      FUN_TRI = Y*W_OVER_D
      RETURN
      END
C     ***********
C     *         *
C     * FUN_ODD
C     *         *
C     ***********

      DOUBLE PRECISION FUNCTION FUN_ODD(Y)

C     Return the vertical diameter normalized width of the 
C     opening at the vertical diameter normalized distance 
C     from the invert. 

      DOUBLE PRECISION Y

      INCLUDE 'orfshape.cmn'

C     Local
      INTEGER NTAB

      REAL ARG, PDV, FUNC
C*********************************************************************
      ARG = Y
      CALL LKTAB(ORIFICE_SHAPE_TAB, ARG, 0,
     O                   FUNC, NTAB, PDV)
      FUN_ODD = FUNC
      RETURN
      END
C     ***********
C     *         *
C     * ORF
C     *         *
C     ***********

      DOUBLE PRECISION FUNCTION ORF(Y, A, B, H, FUN)

C     Compute the orifice/weir flow integrand for the orifice 
C     width given by FUN.  We use the W. Kahan rescaling 
C     to improve accuracy of integration. 
      
      DOUBLE PRECISION Y, A, B, H, FUN

      EXTERNAL FUN

C     Local
      DOUBLE PRECISION Z
C**********************************************************************
      Z = 0.25D0*(B - A)*Y*(3.D0 - Y*Y) + .5D0*(B + A)
      ORF = FUN(Z)*DSQRT(H - Z)*(1 - Y*Y)
      RETURN
      END
C     ***********
C     *         *
C     * WIDTH_ORF
C     *         *
C     ***********

      DOUBLE PRECISION FUNCTION WIDTH_ORF(Y, A, B, FUN)

C     Compute the integrand for the normalized area of an
C     orifice opening. 
      
      DOUBLE PRECISION Y, A, B, FUN

      EXTERNAL FUN

C     Local
      DOUBLE PRECISION Z
C**********************************************************************
      Z = 0.25D0*(B - A)*Y*(3.D0 - Y*Y) + .5D0*(B + A)
      WIDTH_ORF = FUN(Z)*(1 - Y*Y)
      RETURN
      END

C     ***********
C     *         *
C     * FIND_FREE_ORIFICE
C     *         *
C     ***********

      SUBROUTINE FIND_FREE_ORIFICE(STDOUT, HW, TWF, FUN, QFREE,
     A                              NATURE)

C     Find the relative free flow for an orifice with the given relative
C     headwater piezometric level.  

      IMPLICIT NONE

      CHARACTER*2 NATURE

      INTEGER STDOUT

      REAL HW, TWF, QFREE

      DOUBLE PRECISION FUN

      EXTERNAL FUN

      INCLUDE 'orfshape.cmn'

C     Local

      INTEGER I

      DOUBLE PRECISION SUM, ORF, A, B, HUB

      EXTERNAL ORF

C**********************************************************************
C     We assume that the free flow is reduced as soon as the 
C     tailwater is at or above the orifice invert.  
      TWF = 0.0

C     If the relative headwater level is <= 1.0, then we have
C     weir flow, else we have orifice flow. 
      HUB = HW
      IF(HUB.LE.1.D0) THEN
C       This is weir flow.  
        A = 0.D0
        B = HUB
        NATURE = 'FW'
      ELSE
C       This is orifice flow
        A = 0.D0
        B = 1.D0
        NATURE = 'FO'
      ENDIF

      SUM = 0.D0
      DO 100 I=1,ORF_NUM
        SUM = SUM + ORF_W(I)*ORF(ORF_X(I), A, B, HUB, FUN)
100   CONTINUE

      SUM = 0.75D0*(B - A)*SUM
      QFREE = SUM
      RETURN
      END
C     ***********
C     *         *
C     * FIND_SUB_ORIFICE
C     *         *
C     ***********

      SUBROUTINE FIND_SUB_ORIFICE(STDOUT, HW, TW, FUN, QSUB,
     A                           NATURE)

C     Find the relative submerged flow for an orifice with the given 
C     relative headwater and tailwater piezometric levels.

      IMPLICIT NONE

      CHARACTER*2 NATURE

      INTEGER STDOUT

      REAL HW, TW, QSUB

      DOUBLE PRECISION FUN

      EXTERNAL FUN

      INCLUDE 'orfshape.cmn'

C     Local

      INTEGER I

      DOUBLE PRECISION SUM, ORF, A, B, HUB, HDB, QFREE, WIDTH_ORF

      EXTERNAL ORF, WIDTH_ORF

C**********************************************************************
C     If the relative headwater level is <= 1.0, then we have
C     weir flow, else we have orifice flow. 
      HUB = HW
      HDB = TW
      IF(HDB.GE.1.0) THEN
C       The opening is completely submerged. 
C       Compute flow using head difference and the area of the
C       opening. 
        NATURE = 'SO'

C       Now compute the full area.   
        A = 0.D0
        B = 1.0
        SUM = 0.D0
        DO 300 I=1,ORF_NUM
          SUM = SUM + ORF_W(I)*WIDTH_ORF(ORF_X(I), A, B, FUN)
300     CONTINUE
        SUM = 0.75D0*(B - A)*SUM*DSQRT(HUB - HDB)
        QSUB =  SUM
      ELSE
        IF(HUB.LE.1.D0) THEN
C         This is weir flow.  Define the free part
          A = HDB
          B = HUB
          NATURE = 'SW'
        ELSE
C         This is orifice flow
C         The opening is only partially submerged.   Define the
C         free part.
          NATURE = 'SO'
          A = HDB
          B = 1.0
        ENDIF
        SUM = 0.D0
        DO 100 I=1,ORF_NUM
          SUM = SUM + ORF_W(I)*ORF(ORF_X(I), A, B, HUB, FUN)
100     CONTINUE
        SUM = 0.75D0*(B - A)*SUM
        QFREE = SUM

C       Now compute the partial area.   
        A = 0.D0
        B = HDB
        SUM = 0.D0
        DO 200 I=1,ORF_NUM
          SUM = SUM + ORF_W(I)*WIDTH_ORF(ORF_X(I), A, B, FUN)
200     CONTINUE
        SUM = 0.75D0*(B - A)*SUM*DSQRT(HUB - HDB)
        
        QSUB = QFREE + SUM
      ENDIF

      QFREE = SUM
      RETURN
      END
C     ***********
C     *         *
C     * FIND_SUB_ORIFICE_FLOW
C     *         *
C     ***********

      SUBROUTINE FIND_SUB_ORIFICE_FLOW(STDOUT, HW_REL, TW_REL, SHAPE,
     A                               QSUB, NATURE)
C     Find submerged flow cases. 

      IMPLICIT NONE

      CHARACTER NATURE*2, SHAPE*6

      INTEGER STDOUT

      REAL HW_REL, TW_REL, QSUB

C     Local
      DOUBLE PRECISION FUN_CIRC, FUN_RECT, FUN_TRI, FUN_ODD
      EXTERNAL FUN_CIRC, FUN_RECT, FUN_TRI, FUN_ODD
C**********************************************************************
      IF(SHAPE.EQ.'CIRC') THEN
        CALL FIND_SUB_ORIFICE(STDOUT, HW_REL, TW_REL, FUN_CIRC, 
     A                        QSUB, NATURE)
      ELSEIF(SHAPE.EQ.'RECT') THEN
        CALL FIND_SUB_ORIFICE(STDOUT, HW_REL, TW_REL, FUN_RECT, 
     A                        QSUB, NATURE)
      ELSEIF(SHAPE.EQ.'TRI') THEN
        CALL FIND_SUB_ORIFICE(STDOUT, HW_REL, TW_REL, FUN_TRI, 
     A                        QSUB, NATURE)
      ELSEIF(SHAPE.EQ.'ODD') THEN
        CALL FIND_SUB_ORIFICE(STDOUT, HW_REL, TW_REL, FUN_ODD, 
     A                        QSUB, NATURE)
      ENDIF
      RETURN
      END
C     ***********
C     *         *
C     * FIND_FREE_ORIFICE_FLOW
C     *         *
C     ***********

      SUBROUTINE FIND_FREE_ORIFICE_FLOW(STDOUT, HW_REL, TWF_REL,
     A                                  SHAPE, QFREE, NATURE)

C     Find free orifice flow.  Could be weir flow if the orifice is
C     not submerged enough. 

      IMPLICIT NONE

      CHARACTER NATURE*2, SHAPE*6

      INTEGER STDOUT

      REAL HW_REL, TWF_REL, QFREE

C     Local
      DOUBLE PRECISION FUN_CIRC, FUN_RECT, FUN_TRI, FUN_ODD
      EXTERNAL FUN_CIRC, FUN_RECT, FUN_TRI, FUN_ODD
C**********************************************************************      

      IF(SHAPE.EQ.'CIRC') THEN
        CALL FIND_FREE_ORIFICE(STDOUT, HW_REL, TWF_REL, FUN_CIRC, 
     A                         QFREE, NATURE)
      ELSEIF(SHAPE.EQ.'RECT') THEN
        CALL FIND_FREE_ORIFICE(STDOUT, HW_REL, TWF_REL, FUN_RECT, 
     A                         QFREE, NATURE)
      ELSEIF(SHAPE.EQ.'TRI') THEN
        CALL FIND_FREE_ORIFICE(STDOUT, HW_REL, TWF_REL, FUN_TRI, 
     A                         QFREE, NATURE)
      ELSEIF(SHAPE.EQ.'ODD') THEN
        CALL FIND_FREE_ORIFICE(STDOUT, HW_REL, TWF_REL, FUN_ODD, 
     A                         QFREE, NATURE)
      ENDIF

      RETURN
      END
C     ***********
C     *         *
C     * ORIFICE
C     *         *
C     ***********

      SUBROUTINE   ORIFICE(GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR, FTP)
 
C     Compute a 2-D table of type 13 for a vertical orifice.
 
      IMPLICIT NONE
C     Dummy arguments
      INTEGER EFLAG, STDIN, STDOUT, STDTAB, FTP
      INTEGER TABDIR(*)
      REAL GRAV
 
C     DEFINITIONS 
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TABDIR - Table directory to remember table numbers
 
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'orfshape.cmn' 
      INCLUDE 'tabupgrade.cmn' 

C     Local
      INTEGER MAXN
      PARAMETER (MAXN=8)
      INTEGER  I, IHUP, IPFD, J,JBASE, NFRAC, NHUP, N_GLOBAL, TAB,
     A         TABLT, TABGT, NN, APPTAB, ORITAB, IT, NITEM,
     B         ITEM_START(MAXN), ITEM_END(MAXN), L, W1, W2, 
     C         NUMBER, IDLENA, IDLENB, W1A, W2A, ftpup, verbose

      REAL  BIGERR, DROP,  FDROP, FDVEC(PMXNHU), HEAD_DATUM,
     A     HUPVEC(PMXNHU), LIMPFD, MAXZUP, LIPREC,  MAXHUP, MINHUP,
     B     MINPFD, PFDTMP(PMXFRC), PFDVEC(PMXFRC), POW,  QFREE, QHAT, 
     C     QMAT(PMXNHU,PMXFRC), RERR, RMS_GLOBAL, WORK(PMXFRC), 
     D     XBRK(PMXFRC), TWE, TWF, FNUMBER, HTEMP, INVERT,
     E     HW_REL, TW_REL, HW, TW, CE, CW, CD, ORIFICE_LIMIT,
     F     SQRT_OF_2G, QSUB, APPZB, YA, APPA, HW_RELT, W, D,
     G     OFFSET, HWE, TWF_REL, TWEF, TWF_RELT, PFD, HERR, zrhufd

      real*8 FUN_CIRC, FUN_RECT, FUN_TRI, FUN_ODD,
     a         easting, northing
      CHARACTER LABEL*50, LINE*80, SHAPE*6, EDGE*6, NATURE*2,
     A          CQ*8, JUST*5, TABID*16, APPTABID*16, ORITABID*16,
     B          TEMPA*25, TEMPB*25, KEY*16,
     c   zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     Intrinsics
      INTRINSIC ABS, FLOAT, LOG, SQRT
 
 
C     External names
      CHARACTER*16 GET_TABID
      INTEGER LENSTR
      EXTERNAL CHKTAB, inline, LSTOPF, TABCHK, TWDOUT, 
     A         FUN_CIRC, FUN_RECT, FUN_TRI, FUN_ODD, VAR_DECIMAL,
     B         GET_INTERNAL_TAB_NUMBER, READ_TABID, LENSTR
C     ************************************FORMATS**********************
 1    FORMAT(7X,I5)
 2    FORMAT(6X,A)
4     FORMAT(I6,2A6,5F6.0)
6     FORMAT(2I6,5F6.0)
 
 50   FORMAT(/,' Table Id= ',A,' for type 13 table for orifice.')
 52   FORMAT(/,' Label=',A50)
54    FORMAT(' ','Number  Shape   Edge VertD HoriD Invrt OrifC WeirC')
55    FORMAT(A,'AppTbId',A,'OriTbId MxZup MnHup MxPFD MnPFD LIPrc')
56    FORMAT(1X,I6,1X,A6,1X,A6,2F6.2,F6.2,2F6.3)
58    FORMAT(A,A,2F6.1,2F6.3,F6.2)
 60   FORMAT(/,' Datum for defining heads=',F10.2)
62    FORMAT(/,' *ERR:580* Shape is ODD but orifice table number,',
     A    ' ORITAB, is blank.')
 72   FORMAT(/,'Upstream head=',F9.4,' Elevation=',F10.4,
     A        ' Free flow=',A8)
 74   FORMAT(/,
     A '  Partial  Drop    Elev.   Head    Flow Discharge Local',/,
     C '   free    sect.   sect.   sect.   Code           power',/,
     E '   drop    1->4     4       4',/,
     F ' --------  ------  ------  ------   --- --------- ------')
 75   FORMAT(1X,F8.4,F8.3,F8.3,F8.3,3X,A3,2X,A8,F7.2)
78    FORMAT(/,'  *ERR:634* Maximum ups head=',F8.2,' <= 0.',
     A ' Max ups elev=',F8.2,' and head datum=',F8.2)
 80   FORMAT(/,' *WRN:595* Minimum non-zero upstream head=',F8.2,
     A     ' <= 0.',' Setting to 0.15')
 86   FORMAT(/,' Maximum relative error=',F6.3,
     A   '  Upstream head=',F9.4,/,'  and partial free drop=', F8.5)
 88   FORMAT(/' Root-mean-squared error=',F6.3,' N in sample=',I5)
 89   FORMAT('  Processing ORIFICE TabId= ',A)
 99   FORMAT(/,' *ERR:635* Interpolation precision tables missing.',
     A   ' Check for version',/,11X,'of file: TYPE5.TAB.')
C**********************************************************************
      SQRT_OF_2G = SQRT(2.*GRAV)
      JUST = 'RIGHT'
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
C     Define the Guass rule values for numerical integration.
      ORF_NUM = 10
      CALL GRULE(ORF_NUM,
     O                   ORF_X, ORF_W)


      CALL inline(STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'TAB',
     O                EFLAG, TABID, TAB)
C      READ(LINE,1,ERR=991) TAB
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID))
      WRITE(*,89) TABID(1:LENSTR(TABID))
 
      CALL TABCHK(STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)
 
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

      CALL inline(STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) LABEL
      WRITE(STDOUT,52) LABEL
 
C     Get the heading line.
      CALL inline(STDIN, STDOUT,
     O           LINE)
C      WRITE(STDOUT,54) LINE
      
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)

      CALL inline(STDIN, STDOUT,
     O           LINE)

C      READ(LINE,4,ERR=991) NUMBER, SHAPE, EDGE, D, W, INVERT, CD, CW
      CALL READ_ORIFICE_ITEMS1(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     O                         NUMBER, SHAPE, EDGE, D, W, 
     O                         INVERT, CD, CW)
      IF(EDGE.EQ.'     ') EDGE = 'SHARP'
      IF(NUMBER.EQ.0) NUMBER = 1
      WRITE(STDOUT,54)
      L = LENSTR(SHAPE)
      TEMPA = ' '
      TEMPB = ' '
      TEMPA(1+6-L:) = SHAPE
      L = LENSTR(EDGE)
      TEMPB(1+6-L:) = EDGE
      WRITE(STDOUT, 56) NUMBER, TEMPA(1:6), TEMPB(1:6), D, W, 
     A                  INVERT, CD, CW
C      CALL STRIP_L_BLANKS(SHAPE)
C      CALL STRIP_L_BLANKS(EDGE)
      IF(SHAPE.EQ.'ROUND'.OR.SHAPE.EQ.'CIRCLE'.OR.SHAPE.EQ.'CIRC') THEN
        SHAPE = 'CIRC'
        W = D
      ENDIF

      IF(SHAPE.EQ.'OTHER') THEN
        SHAPE = 'ODD'
      ENDIF
      IF(SHAPE.EQ.'INVTRI') THEN
        SHAPE = 'TRI'
      ENDIF
      FNUMBER = NUMBER
C     Get the heading line.
      CALL inline(STDIN, STDOUT,
     O           LINE)
C      WRITE(STDOUT,54) LINE
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)

      CALL inline(STDIN, STDOUT,
     O           LINE)
C      READ(LINE,6,ERR=991) APPTAB, ORITAB, MAXZUP, MINHUP,
C     A       LIMPFD, MINPFD, LIPREC
      CALL READ_ORIFICE_ITEMS2(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         IDLENA, IDLENB, 
     O                         APPTAB, ORITAB, MAXZUP, MINHUP,
     O                         LIMPFD, MINPFD, LIPREC)
      APPTABID = GET_TABID(APPTAB)
      ORITABID = GET_TABID(ORITAB)
      W1 = ITEM_END(1) - ITEM_START(1) + 1
      W2 = ITEM_END(2) - ITEM_START(2) + 1
      W1A = W1
      W2A = W2
      IF(W1A.LT.8) W1A = 8
      IF(W2A.LT.8) W2A = 8
      TEMPA = ' '
      TEMPB = ' '
      WRITE(STDOUT,55) TEMPA(1:W1A-8+1), TEMPB(1:W2A-8+1)
      IF(W1.LT.8) W1 = 8
      IF(W2.LT.8) W2 = 8
      TEMPA(1+W1-IDLENA:) = APPTABID
      TEMPB(1+W2-IDLENB:) = ORITABID
      WRITE(STDOUT,58) TEMPA(1:W1), TEMPB(1:W2), MAXZUP, 
     A                 MINHUP, LIMPFD, MINPFD, LIPREC

      IF(APPTAB.GT.0) THEN
C       Check on the existence of the table and its type.
        IT = APPTAB
        CALL CHKTAB(20, STDOUT, FTPNT, PMXTAB,
     M            APPTAB,
     O            EFLAG)
        IF(EFLAG.EQ.0) THEN
          CALL FNDELV(IT, STDOUT,
     O                      EFLAG, APPZB)
        ENDIF
      ENDIF
      IF(SHAPE.EQ.'ODD') THEN
        IF(ORITAB.GT.0) THEN
C         Check on existence of the table for the orifice shape.
          CALL CHKTAB(2, STDOUT, FTPNT, PMXTAB,
     M            ORITAB,
     O            EFLAG)
          ORIFICE_SHAPE_TAB = ORITAB
        ELSE
          WRITE(STDOUT,62) 
          EFLAG = 1
        ENDIF
      ENDIF      
      
      IF(EFLAG.NE.0) RETURN
 
C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.

C     Define the width to vertical diameter ratio
      IF(D.GT.0.0) THEN
        W_OVER_D = DBLE(W)/DBLE(D)
      ENDIF

C     Compute the upstream head sequence.  First, find the head
C     range requested. 
      HEAD_DATUM = INVERT
      WRITE(STDOUT,60) HEAD_DATUM
      MAXHUP = MAXZUP - HEAD_DATUM
      IF(MAXHUP.LE.0.0) THEN
        WRITE(STDOUT,78) MAXHUP, MAXZUP, HEAD_DATUM
        MAXHUP = 2.0
        EFLAG = 1
      ENDIF
      IF(MINHUP.LE.0.0) THEN
        WRITE(STDOUT,80) MINHUP
        MINHUP = 0.125*D
      ENDIF      
      IF(MINHUP.GT.D) THEN
        MINHUP = 0.125*D
      ENDIF
C     Compute the spacing for weir flow.
      HTEMP = D
      POW = 1.5
      OFFSET = 0.0
      CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, MINHUP, HTEMP, 
     I            LIPREC, PMXNHU,
     O            NHUP, XBRK, EFLAG)

      DO 190 I=1,NHUP
        HUPVEC(I) = XBRK(I)
190   CONTINUE
 
C     Compute the orifice flow region.  
      POW = 0.5
      CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, HTEMP+MINHUP, 
     I            MAXHUP, LIPREC, PMXNHU,
     O            NN, XBRK, EFLAG)
      
      DO 191 I=1,NN
        NHUP = NHUP + 1
        HUPVEC(NHUP) = XBRK(I)
191   CONTINUE

C     Compute the proportions of free drop. 
      POW = 0.5
      OFFSET = 0.0
      CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, MINPFD, LIMPFD,
     I            LIPREC, PMXFRC,
     O            NN, XBRK, EFLAG)
 
      WORK(1) = 0.0
      DO 195 I=1,NN
        WORK(I+1) = XBRK(I)
 195  CONTINUE
      NFRAC = NN + 1

      WORK(NFRAC+1) = 1.0
      NFRAC = NFRAC + 1
 
C     Transfer to PFDVEC and insert the intermediate points for computing
C     an estimated interpolation error. 
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

      RMS_GLOBAL = 0.0
      N_GLOBAL = 0.0
      BIGERR = 0.0

      DO 1000 IHUP= 1, NHUP
        HW = HUPVEC(IHUP)
        HW_REL = HW/D
        HWE = HW + HEAD_DATUM

        IF(APPTAB.GT.0) THEN
C         Define the approach area.
          YA = HWE - APPZB
          CALL LKTA(APPTAB,
     M                  YA,
     O                  APPA)
        ELSE
          APPA = -1.0
        ENDIF
                
C       Find the free flow values.
        CALL FIND_FREE_ORIFICE_FLOW(STDOUT, HW_REL, TWF_REL, SHAPE,
     A                           QFREE, NATURE)
        TWF = D*TWF_REL
        TWEF = TWF + HEAD_DATUM
        FDROP = HW - TWF
C       Interpolate for the coefficient to use. 
        ORIFICE_LIMIT = 1.5*D        
        IF(HW.LE.D) THEN
          CE = CW
        ELSEIF(HW.GE.ORIFICE_LIMIT) THEN
          CE = CD
        ELSE
          CE = CW + (HW - D)*(CD - CW)/(ORIFICE_LIMIT - D)
        ENDIF

        QFREE = QFREE*FNUMBER*CE*SQRT(D)*D**2*SQRT_OF_2G

        IF(APPA.GT.0.0) THEN
C         Make an approx. adjustment for approach velocity head
C         as induced by the flow through the orifice itself.  
          HW_RELT = (HW + (QFREE/APPA)**2/(2.*GRAV))/D
          CALL FIND_FREE_ORIFICE_FLOW(STDOUT, HW_RELT, TWF_RELT, SHAPE,
     A                           QFREE, NATURE)
          QFREE = QFREE*FNUMBER*CE*SQRT(D)*D**2*SQRT_OF_2G
        ENDIF
        QMAT(IHUP,NFRAC) = QFREE
C       Set the flow at zero partial free drop to 0.0
        QMAT(IHUP,1) = 0.0
        FDVEC(IHUP) = FDROP

        CALL VAR_DECIMAL(QFREE,
     O                   CQ)
        WRITE(STDOUT,72) HW, HWE, CQ
        WRITE(STDOUT,74)        
        CALL VAR_DECIMAL(QFREE,
     O                   CQ)
        WRITE(STDOUT,75) PFDVEC(NFRAC), FDROP,
     A                    TWEF, TWF, NATURE, CQ
        
        DO 400 J= NFRAC-1, 2, -1
          PFD = PFDVEC(J)
          DROP = FDROP*PFD
          TW = HW - DROP
          TWE = TW + HEAD_DATUM
          TW_REL = TW/D
C         Find the submerged flow values.
          CALL FIND_SUB_ORIFICE_FLOW(STDOUT, HW_REL, TW_REL, SHAPE,
     A                               QSUB, NATURE)
          QSUB = QSUB*FNUMBER*CE*SQRT(D)*D**2*SQRT_OF_2G
          IF(APPA.GT.0.0) THEN
C           Make an approximate adjustment for approach velocity
C           head. 
            HW_REL = (HW + (QSUB/APPA)**2/(2.*GRAV))/D
            CALL FIND_SUB_ORIFICE_FLOW(STDOUT, HW_REL, TW_REL, SHAPE,
     A                               QSUB, NATURE)
            QSUB = QSUB*FNUMBER*CE*SQRT(D)*D**2*SQRT_OF_2G
          ENDIF

          QMAT(IHUP,J) = QSUB
          POW = LOG(QSUB/QMAT(IHUP,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
          CALL VAR_DECIMAL(QSUB,
     O                     CQ)
          WRITE(STDOUT,75) PFDVEC(J), DROP,
     A                    TWE, TW, NATURE, CQ, POW
400     CONTINUE

C       Compute approximate maximum error and report
        DO 600 J=3,NFRAC-2,2
          QHAT = 0.5*(QMAT(IHUP,J-1) + QMAT(IHUP,J+1))
          RERR = ABS(QHAT - QMAT(IHUP,J))/QMAT(IHUP,J)
          RMS_GLOBAL = RMS_GLOBAL + RERR*RERR
          N_GLOBAL = N_GLOBAL + 1
          IF(RERR.GT.BIGERR) THEN
            BIGERR = RERR
            HERR = HW
            IPFD = J
          ENDIF
 600    CONTINUE
 
C       Eliminate the checking values from QMAT
        JBASE = 3
        DO 700 J=4,NFRAC-1,2
          QMAT(IHUP,JBASE) = QMAT(IHUP,J)
          JBASE = JBASE + 1
 700    CONTINUE
        QMAT(IHUP,JBASE) = QMAT(IHUP,NFRAC)
 
1000  CONTINUE

 900    CONTINUE
        PFDTMP(1) = PFDVEC(1)
        PFDTMP(2) = PFDVEC(2)
C       Eliminate checking values of PFD
        JBASE = 3
        DO 910 J=4,NFRAC-1,2
          PFDTMP(JBASE) = PFDVEC(J)
          JBASE = JBASE + 1
 910    CONTINUE
        PFDTMP(JBASE) = PFDVEC(NFRAC)
 
        zrhufd = 0.0
        CALL TWDOUT(STDOUT, STDTAB, TAB, LABEL, NHUP, JBASE,
     I              HUPVEC, FDVEC, PFDTMP, QMAT, HEAD_DATUM,
     I              13, ' ORIFICE', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)
 
      WRITE(STDOUT,86) BIGERR, HERR,
     A                 PFDVEC(IPFD)
 
      RMS_GLOBAL = SQRT(RMS_GLOBAL/FLOAT(N_GLOBAL))
      WRITE(STDOUT,88) RMS_GLOBAL, N_GLOBAL
 
      if(twod_cubic_out.eq.'YES') then
        verbose = 1
        CALL twodfit
     I           (STDOUT, TAB, NHUP, JBASE, 
     I            HUPVEC, FDVEC,
     I            PFDTMP, QMAT, HEAD_DATUM,
     I            13, ' ORIFICE', zrhufd, verbose,
     M            FTP,
     O            EFLAG, ftpup)
      endif

      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END

