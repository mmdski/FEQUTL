C
C
C
      SUBROUTINE READ_FLOODWAY_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         ITAB, OPT, ELEV,  BOT, LEFT, 
     O                         RIGHT, LOSS, flow, IDLENA)

C     Get the items of data a FLOODWAY-table entry line.

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A       ITAB, IDLENA, EFLAG
    
      REAL ELEV, flow

      CHARACTER BOT*8, LEFT*8, LOSS*8, OPT*4, RIGHT*8
      CHARACTER LINE*(*)

C     Local

      INTEGER IE, IS, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:775* Only ',I3,' items given in ',
     A   'floodway description line.  Need st least seven items.')
52    FORMAT(/,' *ERR:753* Conversion error in field# ',I2,' in:',/,
     A        A)
C***********************************************************************
      IF(NITEM.LT.7) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

C     Process the cross-section table id
      N = 1
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENA = LEN_TRIM(KEY)
      IF(KEY(1:1).EQ.'-') THEN
        ITAB = -1
      ELSE
C       Convert from the table id to an internal number.
        IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C         We have an id given.
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 ITAB)
        ELSE
          ITAB = 0   
        ENDIF
      ENDIF

C     Process the adjustment option
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      OPT = TPC(1:4)

C     Process the standard-flood elevation
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        ELEV = 0.0
      ELSE
        READ(TPC,*,ERR=999) ELEV
      ENDIF
      
C     Process the invert-elevation field
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      BOT = TPC(1:8)

C     Process the left-hand encroachment limit
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      LEFT = TPC(1:8)

C     Process the right-hand encroachment limit
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      RIGHT = TPC(1:8)

C     Process the cross-section-loss value
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      LOSS = TPC(1:8)


c     Process the optional item for flow.  Used to compute velocity in the 
c     floodway.
      flow = 0.0 
      if(nitem.eq.8) then
        n = 8
        IS = ITEM_START(N)
        IE = ITEM_END(N)
        TPC = LINE(IS:IE)
        CALL STRIP_L_BLANKS(
     M                    TPC)
        IF(TPC(1:1).EQ.' ') THEN
          flow = 0.0
        ELSE
          READ(TPC,*,ERR=999) flow
        ENDIF
      endif

      RETURN

999   CONTINUE
      WRITE(STDOUT,52) N, LINE
      STOP 'Abnormal stop.  Errors found.' 
      END
C
C
C
      SUBROUTINE   FLDIN
     I                  (STDIN, STDOUT, STDFLD,
     O                   EFLAG, FHEAD)
 
C     + + + PURPOSE + + +
C     Input the floodway specifications.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDFLD, STDIN, STDOUT
      CHARACTER FHEAD*91
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDFLD - Fortran unit number for Floodway file
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FHEAD   - Descriptive heading
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'fldway.cmn'
 
C     + + + SAVED VALUES + + +
      CHARACTER BLANK*8
      SAVE BLANK
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IDLENA, ITABA, MAXN, NITEM

      PARAMETER (MAXN=8)

      INTEGER ITEM_START(MAXN), ITEM_END(MAXN)

      REAL ELEV, flow
      CHARACTER BOT*8, FILNAM*64, LEFT*8, LINE*80, LOSS*8, OPT*4,
     A          RIGHT*8, JUST*5, TABID*16, HEAD*80
      LOGICAL THERE
 
C     + + + EXTERNAL NAMES + + +
      CHARACTER GET_TABID*16
      EXTERNAL INLINE, GET_TABID, os_file_style
 
C     + + + DATA INITIALIZATIONS + + +
      DATA BLANK/'    '/
 
C     + + + OUTPUT FORMATS + + +
 10   FORMAT(/,
     A' Flood way option set.  All tables listed above will be ',
     B   /,' processed with the flood way option.  Tables encountered',
     C   /,' in input but not in the table will be processed normally.')
C 50   FORMAT(/,' Table Identifier  Option  BF Elevation  FEQ Invert',
C     A' Left-----  Right---      Loss')

50    FORMAT(/,A)
 52   FORMAT(1X,A16,4X,A4,F14.2,4X,A8,2X,A8,2X,A8,F10.3,f10.1)

C***********************************************************************
C     Set heading-dependent processing to right justified
      JUST = 'RIGHT'

C     Create output heading
      FHEAD =' Table Identifier  Option  BF Elevation  FEQ Invert Left--
     A---  Right---      Loss      Flow'

      CALL INLINE
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(5X,A)',ERR=991) FILNAM
      IF(FILNAM.NE.' ') THEN
        CALL MAYBE_ADD_HOME(
     M                      filnam)

        call os_file_style(
     m                              filnam)
        WRITE(STDOUT,'('' FLOODWAY TABLE FILE IS:'',A)') FILNAM
        INQUIRE(FILE=FILNAM, EXIST=THERE)
        IF(THERE) THEN
          OPEN(UNIT=STDFLD, FILE=FILNAM, STATUS='OLD')
        ELSE
          WRITE(STDOUT,*) ' FILE NAMED:', FILNAM,' NOT FOUND.'
          WRITE(STDOUT,*) ' CHECK SPELLING OF STANDARD FLOOD FILE.'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
      ENDIF
 
C     Set the flood way flag
 
      FLOOD = 1
 
      CALL INLINE
     I          (STDFLD, STDOUT,
     O           LINE)
      READ(LINE,'(A)') HEAD
      WRITE(STDOUT,'(1X,A80)') HEAD
 
 
C     Input the global value of the loss in conveyance-the fraction of
C     total conveyance to be subtracted from each overbank
 
      CALL INLINE
     I          (STDFLD, STDOUT,
     O           LINE)
      READ(LINE, '(16X, F10.0)') GLBCON
      WRITE(STDOUT,'('' CONVEYANCE LOSS='',F7.2)') GLBCON
 
C     Input the global value of the elevation loss
 
      CALL INLINE
     I          (STDFLD, STDOUT,
     O           LINE)
      READ(LINE,'(15X,F10.0)') GLBELV
      WRITE(STDOUT,'('' ELEVATION LOSS='',F7.2)') GLBELV
 
      CALL INLINE
     I          (STDFLD, STDOUT,
     O           LINE)
C      WRITE(STDOUT,'(1X,A80)') LINE
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)
 

C     Write standard heading to make it simple.  User's heading
C     used only for defining columns.
      WRITE(STDOUT,50) FHEAD

C     Clear the option values
 
      DO 100 I=1,PMXTAB
        FLDOPT(I) ='    '
 100  CONTINUE
 
C     Input the specification until a table number is < 0
 
 200  CONTINUE
        CALL INLINE
     I            (STDFLD, STDOUT,
     O             LINE)


        CALL READ_FLOODWAY_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         ITABA, OPT, ELEV,  BOT, LEFT, 
     O                         RIGHT, LOSS, flow, IDLENA)

 
        IF(ITABA.LE.0) THEN
          WRITE(STDOUT,10)
          CLOSE(STDFLD)
          RETURN
        ENDIF
 
        TABID = GET_TABID(ITABA)

        IF(OPT.NE.'ELEV'.AND.OPT.NE.'CONV'.AND.OPT.NE.'USET'.AND.
     A     OPT.NE.'EQK') THEN
C         INVALID FLOODWAY OPTION
          EFLAG = 1
          WRITE(STDOUT,'(A)') LINE
          WRITE(STDOUT,'(''*ERR:533* INVALID FLOODWAY OPTION IN ABOVE'',
     A              '' LINE'')')
          GOTO 200
        ENDIF
 
        FLDOPT(ITABA) = OPT
        FLDELV(ITABA) = ELEV
        if(nitem.gt.7) then
          fldflow(itaba) = flow
        else
c         A flood flow of 0.0 signals that no flow has been given.  A floodway will never
c         be defined by zero flow so that this should be ok.
          fldflow(itaba) = 0.0
        endif
 
C       Note: elevation must always be given.  other items are optional
        IF(OPT.EQ.'USET') THEN
          IF(LEFT.EQ.BLANK.OR.RIGHT.EQ.BLANK) THEN
            EFLAG = 1
            WRITE(STDOUT,'(A)') LINE
        WRITE(STDOUT,*) ' *ERR:534* IN LINE ABOVE, USET SELECTED BUT',
     A        ' ONE OR BOTH LIMITS MISSING.'
            EFLAG = 1
          ENDIF
        ENDIF
 
        IF(BOT.EQ.BLANK) THEN
C         Bottom profile in feq is same as minimum point in the
C         cross section
 
          FEQBOT(ITABA) = -1.E30
        ELSE
          READ(BOT,'(F8.0)') FEQBOT(ITABA)
        ENDIF
        IF(LEFT.EQ.BLANK) THEN
C         The left side of the cross section has no limit on the
C         flood way
          FLDLT(ITABA) = 1.E30
        ELSE
          READ(LEFT,'(F8.0)') FLDLT(ITABA)
        ENDIF
        IF(RIGHT.EQ.BLANK) THEN
C         The right side of the cross section has no limit on the
C         floodway
          FLDRT(ITABA) = -1.E30
        ELSE
          READ(RIGHT,'(F8.0)') FLDRT(ITABA)
        ENDIF
        IF(LOSS.EQ.BLANK) THEN
C         Take the global loss value
          IF(OPT.EQ.'ELEV') THEN
            FLDLOS(ITABA) = GLBELV
          ELSEIF(OPT.EQ.'EQK') THEN
            FLDLOS(ITABA) = GLBELV
          ELSE
            FLDLOS(ITABA) = GLBCON
          ENDIF
        ELSE
          READ(LOSS,'(F8.0)') FLDLOS(ITABA)
        ENDIF
 
        IF(FLDRT(ITABA).EQ.-1.E30) THEN
          RIGHT = '  -inf  '
        ELSE
          WRITE(RIGHT,'(F8.1)') FLDRT(ITABA)
        ENDIF
        IF(FLDLT(ITABA).EQ.1.E30) THEN
          LEFT = '  +inf  '
        ELSE
          WRITE(LEFT,'(F8.1)') FLDLT(ITABA)
        ENDIF
        IF(FEQBOT(ITABA).EQ.-1.E30) THEN
          BOT = '  same  '
        ELSE
          WRITE(BOT,'(F8.2)') FEQBOT(ITABA)
        ENDIF
        WRITE(STDOUT, 52)
     A    TABID, FLDOPT(ITABA), FLDELV(ITABA), BOT, LEFT,
     B    RIGHT, FLDLOS(ITABA), flow
 
        GOTO 200
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   SETLME
     I                   (ITABA, STDOUT,
     O                    EFLAG, LEFT, RIGHT)
 
C     + + + PURPOSE + + +
C     Set the left and right hand values for the floodway based
C     on a given decrement in elevation from the 100-year flood
C     elevation.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, ITABA, STDOUT
      REAL LEFT, RIGHT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ITABA  - Index for the cross section floodway descriptors
C     STDOUT - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'fldway.cmn'
      INCLUDE 'xscomu.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IL, IR
      REAL DZ, ELEV, FEQBAS
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:529* Too few intersections.  flood-way elevation ',
     A  'likely below section.',/,'   IL=',I5,' IR=',I5)
52    FORMAT(/,' Final encroachment limits: Left=',
     A F10.2,' Right=',F10.2)
C***********************************************************************
C     FIND THE ELEVATION TO USE IN SEARCHING THE CROSS SECTION.
C     MUST ADJUST FOR POTENTIAL DIFFERENCE BETWEEN THE BOTTOM PROFILE
C     AS USED IN FEQ AND IN FEQUTL
 
      IF(FEQBOT(ITABA).LT.-1.E29) THEN
        FEQBAS = ZMINU
      ELSE
        FEQBAS = FEQBOT(ITABA)
      ENDIF
 
      ELEV = FLDELV(ITABA) - FLDLOS(ITABA) - FEQBAS + ZMINU
 
      WRITE(STDOUT,'('' DEFINING ELEVATION='',F10.2)') ELEV
      IF(ELEV.LE.ZMINU) THEN
        WRITE(STDOUT,*) ' *ERR:530* DEFINING ELEVATION BELOW CHANNEL',
     A     ' BOTTOM.'
        EFLAG = 1
        RETURN
      ENDIF
 
C     SEARCH FOR ELEVATION INTERVALS WHICH CONTAIN ELEV
 
      IF(ELEV.GT.ZU(1)) THEN
        WRITE(STDOUT,'('' *ERR:531* FLDWAY ELEV='',F10.2,
     A   '' HIGHER THAN LEFT END ELEV='',F10.2)') ELEV, ZU(1)
        EFLAG = 1
        RETURN
      ENDIF
      IF(ELEV.GT.ZU(NPNTU)) THEN
        WRITE(STDOUT,'('' *ERR:532* FLDWAY ELEV='',F10.2,
     A   '' HIGHER THAN RIGHT END ELEV='',F10.2)') ELEV, ZU(NPNTU)
        EFLAG = 1
        RETURN
      ENDIF
 
C     SEARCH FOR ALL INTERVALS BUT RETAIN ONLY THE FIRST AND LAST BUT
C     ISSUE A WARNING IF THERE ARE ONE OR MORE 'ISLANDS'
 
      IR = 0
      IL = 0
      DO 100 I=2,NPNTU
        IF(ZU(I).LT.ZU(I-1).AND.ELEV.GT.ZU(I).AND.ELEV.LE.ZU(I-1)) THEN
C         FOUND LEFT INTERSECTION
          IF(IL.EQ.0) THEN
            IL = I
C            WRITE(STDOUT,*) ' LEFT INTERSECTION AT:',I
          ELSE
C            WRITE(STDOUT,*) ' LEFT INTERSECTION AT:',I
 
            WRITE(STDOUT,*) ' *WRN:516* MORE THAN ONE LEFT BOUNDARY.',
     A        ' ONLY FIRST ONE RETAINED.'
          ENDIF
        ELSE
          IF(ZU(I).GT.ZU(I-1).AND
     A       .ELEV.LE.ZU(I).AND.ELEV.GT.ZU(I-1)) THEN
C           FOUND A RIGHT INTERSECTION
            IF(IR.EQ.0) THEN
              IR = I
C              WRITE(STDOUT,*) ' RIGHT INTERSECTION AT:',I
            ELSE
C              WRITE(STDOUT,*) ' RIGHT INTERSECTION AT:',I
 
              IR = I
              WRITE(STDOUT,*) ' *WRN:517* MORE THAN ONE RIGHT',
     A            ' BOUNDARY. ONLY LAST ONE RETAINED.'
            ENDIF
          ENDIF
        ENDIF
 100  CONTINUE
 
      IF(IL.EQ.0.OR.IR.EQ.0) THEN
        WRITE(STDOUT,50) IL, IR
        EFLAG = 1
        RETURN
      ENDIF
      DZ = (ZU(IL) - ZU(IL-1))
      IF(DZ.EQ.0.0) THEN
        WRITE(STDOUT,*) ' *BUG:502* ZERO DIFFERENCE IN ELEVATIONS'
        WRITE(STDOUT,*) 'I=',I
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      LEFT = XU(IL-1) +(ELEV - ZU(IL-1))*(XU(IL) - XU(IL-1))/DZ
 
      DZ = (ZU(IR) - ZU(IR-1))
      IF(DZ.EQ.0.0) THEN
        WRITE(STDOUT,*) ' *BUG:502* ZERO DIFFERENCE IN ELEVATIONS'
        WRITE(STDOUT,*) 'I=',I
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      RIGHT = XU(IR-1) + (ELEV - ZU(IR-1))*(XU(IR) - XU(IR-1))/DZ
 
      IF(LEFT.GT.FLDLT(ITABA)) LEFT = FLDLT(ITABA)
      IF(RIGHT.LT.FLDRT(ITABA)) RIGHT = FLDRT(ITABA)
      WRITE(STDOUT,52) LEFT, RIGHT
 
 
      RETURN
      END
C
C
C
      SUBROUTINE   SETLMK
     I                   (TAB, STDOUT, NFAC,
     O                    EFLAG, LEFT, RIGHT)
 
C     + + + PURPOSE + + +
C     Find limits for a conveyance reduction floodway by direct
C     computation of the conveyance of the reduced cross section.
C     There is no unique way of computing the loss of
C     conveyance which is consistent with the intrinsic
C     meaning of conveyance
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDOUT, TAB
      REAL LEFT, NFAC, RIGHT
                   
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     TAB    - Internal table number
C     STDOUT - Fortran unit number for user output and messages
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'fldway.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL ELEV, FEQBAS, FL, FM, FR, KL, KM, KR, KTAR, KTOTAL, LOSS, XL,
     A     XM, XR, XSV(PMXELM), KB, KF, ELEV2, ZSEEK
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSELM

C     ***************************FORMATS********************************
50    FORMAT(/,' *ERR:776* Floodway elevation=',F10.2,
     A ' has conveyance=',1PE12.5,/, 
     B ' < at BFE=',0PF10.2,' with conveyance=',1PE12.5)
52    FORMAT(/,' Desired conveyance-loss fraction on each side=',F10.3)
54    FORMAT(/,' Base-Flood Elev. for Conveyane=',F10.2)
58    FORMAT(/,' Conveyance at base-flood elev.=',1PE12.4)
60    FORMAT(/,' Final encroachment limits: Left=',
     A F10.2,' Right=',F10.2)
C***********************************************************************
C     DETERMINE THE ELEVATION TO USE IN DEFINING THE CONVEYANCE
 
      IF(FEQBOT(TAB).LT.-1.E29) THEN
        FEQBAS = ZMINU
      ELSE
        FEQBAS = FEQBOT(TAB)
      ENDIF
 
      ELEV = FLDELV(TAB) - FEQBAS + ZMINU
 
      WRITE(STDOUT,54) ELEV
 
      IF(ELEV.LE.ZMINU) THEN
        WRITE(STDOUT,*) ' *ERR:537* Base-Flood elevation below',
     A      ' channel invert!'
        EFLAG = 1
        RETURN
      ENDIF
 

C     COMPUTE THE TOTAL CONVEYANCE FOR THE CROSS SECTION BELOW THE
C     base-flood elevataion.
 
      XL = 1.E20
      XR = -1E20
      CALL CXSELM
     I           (ELEV, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I            ZU, SBU, NFAC, XL, XR, LSNU, SNU,
     M            EFLAG,
     O            NU, XSV)
 
      KTOTAL = XSV(5)**2
      ZSEEK = ELEV
      
      
      IF(FLDOPT(TAB).EQ.'EQK ') THEN
C       Use KB for EQK option
        KB = KTOTAL
        WRITE(STDOUT,58) KB

        ELEV2 = ELEV + FLDLOS(TAB)

C       Find conveyance at the new elevation

        XL = 1.E20
        XR = -1E20
        CALL CXSELM
     I           (ELEV2, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I            ZU, SBU, NFAC, XL, XR, LSNU, SNU,
     M            EFLAG,
     O            NU, XSV)
 
        KF = XSV(5)**2
        IF(KF.LE.KB) THEN
          WRITE(STDOUT,50) ELEV2, KF, ELEV, KB
          EFLAG = 1
          RETURN
        ENDIF
        
C       Reset KTOTAL and ZSEEK for this case
        KTOTAL = KF
        ZSEEK = ELEV2
      ENDIF


C     Compute the loss of conveyance that is to be used on 
C     each side
      IF(FLDOPT(TAB).EQ.'EQK ') THEN
        LOSS = 0.5*(KF - KB)/KF
       
      ELSE
        LOSS = FLDLOS(TAB)
      ENDIF
 
      WRITE(STDOUT,52) LOSS
C
C     COMPUTE THE TARGET AMOUNT FOR SETTING THE RIGHT HAND LIMIT
 
 
      KTAR = KTOTAL*(1.0 - LOSS)
 
C     SEARCH FOR AN INTERVAL ON THE CROSS SECTION PERIPHIERY
C     WHICH CONTAINS KTAR . FIND RIGHT HAND LIMIT FIRST.
 
      KR = KTOTAL
      LEFT = XU(1)
      XR = XU(NPNTU)
C      WRITE(STDOUT,*) ' SEARCHING FOR RIGHT LIMIT KTAR=', KTAR
C      WRITE(STDOUT,*) '      XOFF          K'
C      WRITE(STDOUT,'(F10.1,1PE12.5)') XR, KR
      DO 200 I=NPNTU-1, 1, -1
 
        XL = XU(I)
        CALL CXSELM
     I          (ZSEEK, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I           ZU, SBU, NFAC, LEFT, XL, LSNU, SNU,
     M           EFLAG,
     O           NU, XSV)
        KL = XSV(5)**2
 
C        WRITE(STDOUT,'(F10.1,1PE12.5)') XL, KL
 
        IF(KTAR.LE.KR.AND.KTAR.GE.KL) THEN
C         FOUND INTERVAL
          GOTO 210
        ENDIF
 
        XR = XL
        KR = KL
 200  CONTINUE
 
      WRITE(STDOUT,*) ' *ERR:538* NO INTERVAL FOUND FOR KTAR ON RIGHT'
      STOP 'Abnormal stop. Errors found.'
 210  CONTINUE
 
C     FIND THE MATCH VALUE USING BISECTION
 
      FL = KL - KTAR
      FR = KR - KTAR
 
C      WRITE(STDOUT,*) ' BISECTION'
 300  CONTINUE
 
        XM = XL - FL*(XR - XL)/(FR - FL)
 
        CALL CXSELM
     I          (ZSEEK, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I           ZU, SBU, NFAC, LEFT, XM, LSNU, SNU,
     M           EFLAG,
     O              NU, XSV)
 
        KM = XSV(5)**2
 
        FM = KM - KTAR

        IF(ABS(FM/KTAR).LT.1.D-4) GOTO 310
        IF(ABS(XR - XL).LE.EPSDIF) THEN
          WRITE(STDOUT,*) ' Bisection trapped at possible discontinuity'
          WRITE(STDOUT,*) ' at boundary between subsections.'
          WRITE(STDOUT,*) ' Requested reduction cannot be matched'
          GOTO 310
        ENDIF
 
        IF(FM*FL.LT.0.0) THEN
          FR = FM
          XR = XM
        ELSE
          FL = FM
          XL = XM
        ENDIF
        GOTO 300
 
 310  CONTINUE
 
C     FOUND THE BOUNDARY ON THE RIGHT
 
      RIGHT = XM
 
C     NOW FIND THE  TARGET VALUE FOR THE LEFT LIMIT
 
      KTAR = KTOTAL*(1. - 2.*LOSS)
 
C     SEARCH FOR AN INTERVAL ON THE CROSS SECTION PERIPHIERY
C     WHICH CONTAINS KTAR
 
      KL = KM
      XL = XU(1)
C      WRITE(STDOUT,*) ' SEARCHING FOR LEFT LIMIT KTAR=', KTAR
C      WRITE(STDOUT,*) '      XOFF          K'
C      WRITE(STDOUT,'(F10.1,1PE12.5)') XL, KL
      DO 400 I=2, NPNTU
 
        XR = XU(I)
        IF(XR.GE.RIGHT) THEN
          XR = RIGHT - 1.0
        ENDIF
        CALL CXSELM
     I         (ZSEEK, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I          ZU, SBU, NFAC, XR, RIGHT, LSNU, SNU,
     M          EFLAG,
     O          NU, XSV)
        KR = XSV(5)**2
 
C        WRITE(STDOUT,'(F10.1,1PE12.5)') XR, KR
 
        IF(KTAR.LE.KL.AND.KTAR.GE.KR) THEN
C         FOUND INTERVAL
          GOTO 410
        ENDIF
 
        XL = XR
        KL = KR
 400  CONTINUE
 
      WRITE(STDOUT,*) ' *ERR:539* NO INTERVAL FOUND FOR KTAR ON LEFT'
      STOP 'Abnormal stop. Errors found.'
 410  CONTINUE
 
C     FIND THE MATCH VALUE USING BISECTION
 
      FL = KL - KTAR
      FR = KR - KTAR
 
C      WRITE(STDOUT,*) ' BISECTION'
 500  CONTINUE
 
        XM = XL - FL*(XR - XL)/(FR - FL)
 
        CALL CXSELM
     I          (ZSEEK, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I           ZU, SBU, NFAC, XM, RIGHT, LSNU, SNU,
     M           EFLAG,
     O           NU, XSV)
 
        KM = XSV(5)**2
 
        FM = KM - KTAR
        IF(ABS(FM/KTAR).LT.1.D-4) GOTO 510
        IF(ABS(XR - XL).LE.EPSDIF) THEN
          WRITE(STDOUT,*) ' Bisection trapped at possible discontinuity'
          WRITE(STDOUT,*) ' at boundary between subsections.'
          WRITE(STDOUT,*) ' Requested reduction cannot be matched'
          GOTO 510
        ENDIF
 
        IF(FM*FL.LT.0.0) THEN
          FR = FM
          XR = XM
        ELSE
          FL = FM
          XL = XM
        ENDIF
        GOTO 500
 
 510  CONTINUE
 
C     FOUND THE BOUNDARY ON THE LEFT
 
      LEFT = XM
 
C     CHECK FOR THE PRESET LIMITS
 
      IF(LEFT.GT.FLDLT(TAB)) LEFT = FLDLT(TAB)
      IF(RIGHT.LT.FLDRT(TAB)) RIGHT = FLDRT(TAB)
      WRITE(STDOUT,60) LEFT, RIGHT
 
      RETURN
      END
C
C
C
      SUBROUTINE   FNDWAY
     I                   (TAB, STDOUT, EFLAG, NFAC,
     O                    LEFT, RIGHT, area)
 
C     + + + PURPOSE + + +
C     Find the left and right limits for a floodway in the
C     channel.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDOUT, TAB
      REAL LEFT, NFAC, RIGHT, area
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     TAB    - Table number
C     STDOUT - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NFAC   - Factor in Manning's formula
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
c     area   - area of the floodway 
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'fldway.cmn'
      INCLUDE 'xscomu.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL FEQBAS, KPART, KTOTAL, XL, XR, XSV(PMXELM), ZGIVE
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSELM, SETLME, SETLMK, XCHK
C     ***************************FORMAT*********************************
50    FORMAT(/,' Realized conveyence-loss fraction on each side =',
     A           F10.3)
C***********************************************************************
C     CHECK FOR INVALID HORIZONTAL LINE SEGMENTS
 
      CALL XCHK
     I         (STDOUT, NPNTU, ZMINU, XU, NFAC,
     M          ZU)
 
C     SELECT FLOODWAY OPTION
 
C      WRITE(STDOUT,*) ' In FNDWAY TAB=',TAB
C      WRITE(STDOUT,*) ' In FNDWAY FLDOPT(TAB)=', FLDOPT(TAB)
      IF(FLDOPT(TAB).EQ.'ELEV') THEN
        CALL SETLME
     I             (TAB, STDOUT,
     O              EFLAG, LEFT, RIGHT)
      ELSE IF(FLDOPT(TAB).EQ.'CONV'.OR.FLDOPT(TAB).EQ.'EQK') THEN
        CALL SETLMK
     I             (TAB, STDOUT, NFAC,
     O              EFLAG, LEFT, RIGHT)
      ELSE IF(FLDOPT(TAB).EQ.'USET') THEN
        LEFT = FLDLT(TAB)
        RIGHT = FLDRT(TAB)
      ELSE
        RETURN
      ENDIF
 
      IF(EFLAG.GT.0) THEN
        WRITE(STDOUT, *) ' FLOODWAY COMPUTATIONS SUPPRESSED.',
     A    ' ERRORS ENCOUNTERED.'
        RETURN
      ENDIF
 
C     RESET THE LEFT AND RIGHT LIMITS IN THE FLOODWAY TABLE
 
      FLDLT(TAB) = LEFT
      FLDRT(TAB) = RIGHT
 
C     COMPUTE THE REDUCTION IN CONVEYENCE RESULTING FROM THE IMPOSITION
C     OF THE FLOODWAY LIMITS.
 
 
C     DETERMINE THE ELEVATION TO USE IN DEFINING THE CONVEYANCE
 
      IF(FEQBOT(TAB).LT.-1.E29) THEN
        FEQBAS = ZMINU
      ELSE
        FEQBAS = FEQBOT(TAB)
      ENDIF
 
      ZGIVE = FLDELV(TAB) - FEQBAS + ZMINU
      IF(FLDOPT(TAB).EQ.'EQK') THEN
        ZGIVE = ZGIVE + FLDLOS(TAB)
      ENDIF
 
C      WRITE(STDOUT,*) ' ZGIVE=',ZGIVE
C      WRITE(STDOUT,*) ' ZMINU=',ZMINU
C      WRITE(STDOUT,*) ' ZMAXU=',ZMAXU
C      WRITE(STDOUT,*) ' NPNTU=',NPNTU
C      WRITE(STDOUT,*) ' NAVMU=',NAVMU
C      WRITE(STDOUT,*) ' NSUBU=',NSUBU
C      WRITE(STDOUT,*) ' NFAC=',NFAC
C      WRITE(STDOUT,*) 'NU:',(NU(I),I=1,NSUBU)
C
C      DO 9123 I=1,NPNTU
C        WRITE(STDOUT,'(F10.1,F10.2,I5)') XU(I), ZU(I), SBU(I)
C9123  CONTINUE
 
      XL =  999999.
      XR = -999999.
      CALL CXSELM
     I           (ZGIVE, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I            ZU, SBU, NFAC, XL, XR, LSNU, SNU,
     M            EFLAG,
     O            NU, XSV)
 
      KTOTAL = XSV(5)**2
 
      CALL CXSELM
     I           (ZGIVE, STDOUT, NPNTU, NSUBU, NAVMU, ZMINU, ZMAXU, XU,
     I            ZU, SBU, NFAC, LEFT, RIGHT, LSNU, SNU,
     M            EFLAG,
     O            NU, XSV)
      KPART = XSV(5)**2
      area = xsv(3)
 
C      WRITE(STDOUT,*) ' KTOTAL=',KTOTAL
C      WRITE(STDOUT,*) ' KPART=',KPART
 
      FLDLOS(TAB) = (1.0 - KPART/KTOTAL)/2.0
 
      WRITE(STDOUT,50) FLDLOS(TAB)

 
      RETURN
      END
