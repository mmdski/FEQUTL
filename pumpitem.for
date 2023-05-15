C
C
C
      SUBROUTINE READ_AXIALPUMP_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         IDLENA, IDLENB, IDLENC, 
     O                         NET_PUMP_TAB, NUMBER_OF_PUMPS, 
     O                         INLET_LENGTH, INLET_TAB, 
     O                         OUTLET_LENGTH, OUTLET_TAB, 
     O                         INLET_LOSS_FACTOR, H_DES, Q_DES, LABEL)

C     Get the items of data for computing an axial pump table

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A       NET_PUMP_TAB, NUMBER_OF_PUMPS, INLET_TAB, OUTLET_TAB,
     B       IDLENA, IDLENB, IDLENC, EFLAG
      REAL INLET_LENGTH, OUTLET_LENGTH, INLET_LOSS_FACTOR, 
     A      H_DES, Q_DES
      CHARACTER LINE*(*), LABEL*50 

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:755* Only ',I3,' items given in ',
     A   'axial-pump description line.  Need  ten items.')
52    FORMAT(/,' *ERR:753* Conversion error in field ',I1,' in:',/,
     A        A)
C***********************************************************************
      IF(NITEM.LT.10) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

C     Process the net pump table
      N = 1
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENA = LENSTR(KEY)
      IF(KEY(1:1).EQ.'-') THEN
        NET_PUMP_TAB = -1
      ELSE
C       Convert from the table id to an internal number.
        IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C         We have an id given.
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 NET_PUMP_TAB)
        ELSE
          NET_PUMP_TAB = 0   
        ENDIF
      ENDIF

C     Process the number of pumps
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        NUMBER_OF_PUMPS = 0
      ELSE
        READ(TPC,*,ERR=999) NUMBER_OF_PUMPS
      ENDIF

C     Process the inlet conduit length
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        INLET_LENGTH = 0
      ELSE
        READ(TPC,*,ERR=999) INLET_LENGTH
      ENDIF
      
C     Process the inlet conduit table
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENB = LENSTR(KEY)
C     Convert from the table id to an internal number.
      IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C       We have an id given.
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, KEY,
     M                               EFLAG,
     O                               INLET_TAB)
      ELSE
        INLET_TAB = 0   
      ENDIF

C     Process the outlet conduit length
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        OUTLET_LENGTH = 0
      ELSE
        READ(TPC,*,ERR=999) OUTLET_LENGTH
      ENDIF

C     Process the outlet conduit table
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENC = LENSTR(KEY)
C     Convert from the table id to an internal number.
      IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C       We have an id given.
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, KEY,
     M                               EFLAG,
     O                               OUTLET_TAB)
      ELSE
        OUTLET_TAB = 0   
      ENDIF

C     Process the inlet-loss factor
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        INLET_LOSS_FACTOR = 0.0
      ELSE
        READ(TPC,*,ERR=999) INLET_LOSS_FACTOR
      ENDIF

C     Process the design head
      N = 8
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        H_DES = 0.0
      ELSE
        READ(TPC,*,ERR=999) H_DES
      ENDIF

C     Process the design flow
      N = 9
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        Q_DES = 0.0
      ELSE
        READ(TPC,*,ERR=999) Q_DES
      ENDIF

C     Process the table label
      N = 10
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      LABEL = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    LABEL)
      RETURN
999   CONTINUE
      WRITE(STDOUT,52) N, LINE
      STOP 'Abnormal stop.  Errors found.' 
      END
C
C
C
      SUBROUTINE WRITE_AXIALPUMP_ITEMS(
     I                         STDOUT, NITEM, ITEM_START,
     I                         ITEM_END,
     I                         IDLENA, IDLENB, IDLENC, 
     I                         NET_PUMP_TAB, NUMBER_OF_PUMPS, 
     I                         INLET_LENGTH, INLET_TAB, 
     I                         OUTLET_LENGTH, OUTLET_TAB, 
     I                         INLET_LOSS_FACTOR, H_DES, Q_DES, LABEL,
     O                         LINE)

C     Construct a line of output of axial-pump items.

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A       NET_PUMP_TAB, NUMBER_OF_PUMPS, INLET_TAB, OUTLET_TAB,
     B       IDLENA, IDLENB, IDLENC
      REAL INLET_LENGTH, OUTLET_LENGTH, INLET_LOSS_FACTOR, 
     A      H_DES, Q_DES
      CHARACTER LINE*(*), LABEL*50 

C     Local

      INTEGER IE, IS, N, W
      CHARACTER TPC*20, KEY*16

C     Called program units
      CHARACTER GET_TABID*16
      INTEGER LENSTR
      EXTERNAL STRIP_L_BLANKS, LENSTR, GET_TABID
C     ***********************FORMATS************************************
C***********************************************************************
      LINE = ' '
C     Process the net pump table
      N = 1
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      KEY = GET_TABID(NET_PUMP_TAB)
      W = IE - IS + 1
      LINE(IS+W-IDLENA:) = KEY(1:IDLENA)

C     Process the number of pumps
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:5),'(I5)') NUMBER_OF_PUMPS
      W = IE - IS + 1
      LINE(IS+W-5:) = TPC(1:5)
      
C     Process the inlet conduit length
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') INLET_LENGTH
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)
      
C     Process the inlet conduit table
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      KEY = GET_TABID(INLET_TAB)
      W = IE - IS + 1
      LINE(IS+W-IDLENB:) = KEY(1:IDLENB)

C     Process the outlet conduit length
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') OUTLET_LENGTH
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the outlet conduit table
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      KEY = GET_TABID(OUTLET_TAB)
      W = IE - IS + 1
      LINE(IS+W-IDLENC:) = KEY(1:IDLENC)

C     Process the inlet-loss factor
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') INLET_LOSS_FACTOR
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the design head
      N = 8
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') H_DES
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the design flow
      N = 9
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.1)') Q_DES
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the table label
      N = 10
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      LINE(IS+1:) = LABEL
      RETURN
      END
C
C
C
      SUBROUTINE READ_PUMPLOSS_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         IDLENA, IDLENB, IDLENC, 
     O                         LOSS_TAB, NUMBER_OF_PUMPS, 
     O                         INLET_LENGTH, INLET_TAB, 
     O                         OUTLET_LENGTH, OUTLET_TAB, 
     O                         INLET_LOSS_FACTOR, Q_MAX, LABEL)


C     Get the items of data for computing pumploss table

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A       LOSS_TAB, NUMBER_OF_PUMPS, INLET_TAB, OUTLET_TAB,
     B       IDLENA, IDLENB, IDLENC, EFLAG
      REAL INLET_LENGTH, OUTLET_LENGTH, INLET_LOSS_FACTOR, 
     A      Q_MAX
      CHARACTER LINE*(*), LABEL*50 

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:756* Only ',I3,' items given in ',
     A   'pump-loss description line.  Need  nine items.')
52    FORMAT(/,' *ERR:753* Conversion error in field ',I1,' in:',/,
     A        A)
C***********************************************************************
      IF(NITEM.LT.9) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

C     Process the net pump table
      N = 1
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENA = LENSTR(KEY)
      IF(KEY(1:1).EQ.'-') THEN
        LOSS_TAB = -1
      ELSE
C       Convert from the table id to an internal number.
        IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C         We have an id given.
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 LOSS_TAB)
        ELSE
          LOSS_TAB = 0   
        ENDIF
      ENDIF

C     Process the number of pumps
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        NUMBER_OF_PUMPS = 0
      ELSE
        READ(TPC,*,ERR=999) NUMBER_OF_PUMPS
      ENDIF

C     Process the inlet conduit length
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        INLET_LENGTH = 0
      ELSE
        READ(TPC,*,ERR=999) INLET_LENGTH
      ENDIF
      
C     Process the inlet conduit table
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENB = LENSTR(KEY)
C     Convert from the table id to an internal number.
      IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C       We have an id given.
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, KEY,
     M                               EFLAG,
     O                               INLET_TAB)
      ELSE
        INLET_TAB = 0   
      ENDIF

C     Process the outlet conduit length
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        OUTLET_LENGTH = 0
      ELSE
        READ(TPC,*,ERR=999) OUTLET_LENGTH
      ENDIF

C     Process the outlet conduit table
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      IDLENC = LENSTR(KEY)
C     Convert from the table id to an internal number.
      IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C       We have an id given.
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, KEY,
     M                               EFLAG,
     O                               OUTLET_TAB)
      ELSE
        OUTLET_TAB = 0   
      ENDIF

C     Process the inlet-loss factor
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        INLET_LOSS_FACTOR = 0.0
      ELSE
        READ(TPC,*,ERR=999) INLET_LOSS_FACTOR
      ENDIF

C     Process the maximum flow
      N = 8
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        Q_MAX = 0.0
      ELSE
        READ(TPC,*,ERR=999) Q_MAX
      ENDIF


C     Process the table label
      N = 9
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      LABEL = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    LABEL)
      RETURN
999   CONTINUE
      WRITE(STDOUT,52) N, LINE
      STOP 'Abnormal stop.  Errors found.' 
      END
C
C
C
      SUBROUTINE WRITE_PUMPLOSS_ITEMS(
     I                         STDOUT, NITEM, ITEM_START,
     I                         ITEM_END,
     I                         IDLENA, IDLENB, IDLENC, 
     I                         LOSS_TAB, NUMBER_OF_PUMPS, 
     I                         INLET_LENGTH, INLET_TAB, 
     I                         OUTLET_LENGTH, OUTLET_TAB, 
     I                         INLET_LOSS_FACTOR, Q_MAX, LABEL,
     O                         LINE)

C     Construct a line of output of PUMPLOSS items.

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A       LOSS_TAB, NUMBER_OF_PUMPS, INLET_TAB, OUTLET_TAB,
     B       IDLENA, IDLENB, IDLENC
      REAL INLET_LENGTH, OUTLET_LENGTH, INLET_LOSS_FACTOR, 
     A      Q_MAX
      CHARACTER LINE*(*), LABEL*50 

C     Local

      INTEGER IE, IS, N, W
      CHARACTER TPC*20, KEY*16

C     Called program units
      CHARACTER GET_TABID*16
      INTEGER LENSTR
      EXTERNAL STRIP_L_BLANKS, LENSTR, GET_TABID
C     ***********************FORMATS************************************
C***********************************************************************
      LINE = ' '
C     Process the pump loss table id
      N = 1
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      KEY = GET_TABID(LOSS_TAB)
      W = IE - IS + 1
      LINE(IS+W-IDLENA:) = KEY(1:IDLENA)

C     Process the number of pumps
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:5),'(I5)') NUMBER_OF_PUMPS
      W = IE - IS + 1
      LINE(IS+W-5:) = TPC(1:5)
      
C     Process the inlet conduit length
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') INLET_LENGTH
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)
      
C     Process the inlet conduit table
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      KEY = GET_TABID(INLET_TAB)
      W = IE - IS + 1
      LINE(IS+W-IDLENB:) = KEY(1:IDLENB)

C     Process the outlet conduit length
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') OUTLET_LENGTH
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the outlet conduit table
      N = 6
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      KEY = GET_TABID(OUTLET_TAB)
      W = IE - IS + 1
      LINE(IS+W-IDLENC:) = KEY(1:IDLENC)

C     Process the inlet-loss factor
      N = 7
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.2)') INLET_LOSS_FACTOR
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the maximum flow
      N = 8
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      WRITE(TPC(1:8),'(F8.1)') Q_MAX
      W = IE - IS + 1
      LINE(IS+W-8:) = TPC(1:8)

C     Process the table label
      N = 9
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      LINE(IS+1:) = LABEL
      RETURN
      END
C
C
C
      CHARACTER*10 FUNCTION   PUT10
     I                           (X)
 
C     + + + PURPOSE + + +
C     Function to convert a real number into a special compact
C     form of output to retain 6 significant figures for
C     numbers in the range  x < 1e10. 
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL X
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     X      - value to format
 
C     + + + LOCAL VARIABLES + + +
      CHARACTER RESULT*10, WORK*12
C***********************************************************************
      IF(ABS(X).LE.1.01E-9) X = 0.0
      WRITE(WORK,'(1PE12.5)') X
      IF(WORK(10:10).EQ.'+'.AND.WORK(11:11).NE.'0') THEN
        RESULT = ' *********'
      ELSE
        RESULT(1:8) = WORK(1:8)
        RESULT(9:9) = WORK(10:10)
        RESULT(10:10) = WORK(12:12)
      ENDIF
      PUT10 = RESULT
      RETURN
      END
C
C
C
      CHARACTER*10 FUNCTION   PUT10D
     I                           (X)
 
C     + + + PURPOSE + + +
C     Function to convert a real*8 number into a special compact
C     form of output to retain 7 significant figures for
C     numbers in the range  x < 1e10. 
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL*8 X
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     X      - value to format
 
C     + + + LOCAL VARIABLES + + +
      CHARACTER RESULT*10, WORK*12
C***********************************************************************
      IF(ABS(X).LE.1.01E-9) X = 0.D0
      WRITE(WORK,'(1PE12.6)') X
      IF(WORK(10:10).EQ.'+'.AND.WORK(11:11).NE.'0') THEN
C       OVERFLOW
        RESULT = ' *********'
      ELSE
        RESULT(1:8) = WORK(1:8)
        RESULT(9:9) = WORK(10:10)
        RESULT(10:10) = WORK(12:12)
      ENDIF
      PUT10D = RESULT
      RETURN
      END

C
C
C
      SUBROUTINE  OUTPUT_TYPE_234
     I                   (STDOUT, STDTAB, N, TAB, TYPE, FACTOR, 
     I                    HEAD, LABEL, X, F1, F2,
     i                    zone, hgrid, vdatum, unitsys, basis,
     i                    easting, northing)
 
C     + + + PURPOSE + + +
C     Output a function table of type 2, 3, or 4.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, STDOUT, STDTAB, TAB, TYPE
      CHARACTER LABEL*50, HEAD*30, zone*8, hgrid*8, vdatum*8, 
     a        unitsys*8, basis*8
      REAL FACTOR, F1(N), F2(N), X(N)
      real*8 easting, northing
 

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     N      - Number of tabulated levels in the table
C     TAB    - Internal table number
C     TYPE   - type of the table, 2, 3, or 4
C     FACTOR - Factor to put in the FAC position
C     HEAD   - heading to put on the body of the table
C     LABEL  - Label for identification
C     X      - Argument values for table
C     F1      - Function values to output
C     F2     - derivative values to output for types 3 and 4
 
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J, HEAD_LENGTH, LABEL_LENGTH, L
      real term
      CHARACTER ARG_OUT*10, F1_OUT*10, F2_OUT*10, TABID*16
 
C     External functions
      CHARACTER PUT10*10, GET_TABID*16

C     External names
      EXTERNAL VAR_DECIMAL, PUT10, GET_TABID, chk_vdatum_unitsys
            
C     + + + OUTPUT FORMATS + + +
 1    FORMAT(' TABID=',A)
 2    FORMAT('TABID=',A)
 3    FORMAT(' TYPE=',I5)
 4    FORMAT('TYPE=',I5)
 5    FORMAT(' REFL=',F10.3,' FAC=',A10)
 6    FORMAT('REFL=',F10.3,' FAC=',A10)
 7    FORMAT(' ',A,1X,A)
 8    FORMAT(A, 1X,A)
 9    FORMAT(' ',A10,A10,A10)
 10   FORMAT(A10,A10,A10)
 60   format('ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a ' BASIS=',a8,/, 'EASTING=',0PF15.3,' NORTHING=',F15.3)
C***********************************************************************
      call chk_vdatum_unitsys(stdout, vdatum, unitsys, 
     a          ' during output of type 2, 3, or 4 table')
c 
      TABID = GET_TABID(TAB)
      L = len_trim(TABID)      
      WRITE(STDOUT,1) TABID(1:L)
      WRITE(STDTAB,2) TABID(1:L)
 
      WRITE(STDOUT,3) TYPE
      WRITE(STDTAB,4) -TYPE
      if(zone /= 'NONE')  then
c       Output the location information.
        write(stdtab, 60) zone, hgrid, vdatum, unitsys, basis, easting, 
     a        northing
      endif

      ARG_OUT = PUT10(FACTOR)
      WRITE(STDOUT,5) 0.0, ARG_OUT
      WRITE(STDTAB,6) 0.0, ARG_OUT

      IF(TYPE.EQ.2) THEN
        HEAD_LENGTH = 20
      ELSE
        HEAD_LENGTH = 30
      ENDIF 

      LABEL_LENGTH = len_trim(LABEL) 
      IF(LABEL_LENGTH.EQ.0) LABEL_LENGTH = 1
      WRITE(STDOUT,7) HEAD(1:HEAD_LENGTH), LABEL(1:LABEL_LENGTH)
      WRITE(STDTAB,8) HEAD(1:HEAD_LENGTH), LABEL(1:LABEL_LENGTH)
 
      DO 100 J=1,N

        CALL VAR_DECIMAL(X(J),
     O                       ARG_OUT)
        CALL VAR_DECIMAL(F1(J),
     O                       F1_OUT)
        IF(TYPE.GT.2) THEN
          CALL VAR_DECIMAL(F2(J),
     O                       F2_OUT)
          WRITE(STDOUT,9) ARG_OUT, F1_OUT, F2_OUT
          WRITE(STDTAB,10) ARG_OUT, F1_OUT, F2_OUT
        ELSE
          WRITE(STDOUT,9) ARG_OUT, F1_OUT
          WRITE(STDTAB,10) ARG_OUT, F1_OUT
        ENDIF
        
 100  CONTINUE
      if(x(n) > 0.0) then
        term = -1.0
      else
        term = x(n)  - 1.0
      endif
      WRITE(STDTAB,'(F10.1)') term 
      RETURN
      END



C
C
C
      SUBROUTINE   PUMPITEMS
     I                  (OPTION, GRAV, STDIN, STDOUT, STDTAB,
     M                   EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute various functions needed for modeling variable-
C     head and variable-speed pumps. 
C     Options: 
C         SFWMD- Compute the net head for an idealized axial flow
C                pump moving water through a given length of
C                inlet and outlet conduit plus losses at the
C                entrance of the inlet pipe.

C         PUMPLOSS - Compute the loss through the inlet and outlet
C                    conduit plus losses at the entrance of 
C                    the inlet pipe.

C         PUMPADJ- given a pump curve for FEQ and the inlet and
C                  outlet conduits and an entrance loss coeff.,
C                  compute a net pump curve.  
C 
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      CHARACTER OPTION*8
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      REAL GRAV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     OPTION - gives the option
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER NHEAD, MAXN
      PARAMETER (NHEAD=4, MAXN=10)
      INTEGER LOSS_TAB, INLET_TAB, OUTLET_TAB, HEAD_LEN(NHEAD),
     A        EF, INLET_ADRS, OUTLET_ADRS, GROSS_PUMP_TAB,
     B        NET_PUMP_TAB, N, M, TYPE, IT, I,
     C        NUMBER_OF_PUMPS, ITEM_START(MAXN), ITEM_END(MAXN),
     D        NITEM, IDLENA, IDLENB, IDLENC

      REAL INLET_LENGTH, OUTLET_LENGTH, INLET_LOSS_FACTOR,
     A     D_INLET, D_OUTLET, Q_MAX, K_INLET, A_INLET,
     B     K_OUTLET, A, T, DK, DT, B, DB, ARG(PMXOFF), F1(PMXOFF),
     C     F2(PMXOFF), H_DES, Q_DES, Q_FACTOR, Q_UNIT,
     D     FACTOR, AXIAL_H_REL(17),
     E     AXIAL_Q_REL(17)

      real*8 easting,northing

      CHARACTER LABEL*50, LINE*120, HEAD(NHEAD)*120, BODY_HEAD*30,
     A          FLOW_UNITS*5, EXIT_AREA*7, JUST*5, TABID*16,
     b          zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8

C     External Functions
      INTEGER LENSTR
      REAL GETD 
      CHARACTER PUT7*7, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL inline, OUTPUT_TYPE_234, PUT7, LENSTR, GET_TABID,
     A         READ_AXIALPUMP_ITEMS, WRITE_AXIALPUMP_ITEMS

C     Data initializations
      DATA AXIAL_H_REL/1000.,2.96,2.62,2.34,2.10,1.93,1.73,1.59,1.49,
     A                       1.41,1.24,1.00,0.75,0.53,0.18,0.00,-10.0/,
     B     AXIAL_Q_REL/ 0.00,0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,
     C                       0.80,0.90,1.00,1.10,1.20,1.30,1.37,1.37/
C*********************************Formats*******************************
1     FORMAT(6X,F15.0)
2     FORMAT(2I8,F8.0,I8,F8.0,I8,2F8.0,1X,A)
3     FORMAT(2I8,F8.0,I8,F8.0,I8,3F8.0,1X,A)
4     FORMAT(2I8,F8.0,I8,F8.0,I8,F8.0,1X,A)


50    FORMAT(1X,2I8,F8.2,I8,F8.2,I8,F8.2,F8.0,1X,A)
51    FORMAT(/,' Unit for flows in the input is ',F12.8,' ',A5)
52    FORMAT(1X,A)
54    FORMAT(1X,/,' *ERR:726* TabId= ',A,' for inlet conduit cross ',
     A            'section not found.')
56    FORMAT(1X,/,' *ERR:727* TabId=',A,' for outlet conduit cross ',
     A            'section not found.')
58    FORMAT(1X,/,' *BUG: Invalid option in sub. PUMPITEMS.')
60    FORMAT(1X,2I8,F8.2,I8,F8.2,I8,F8.2,F8.1,F8.0, 1X,A)                
62    FORMAT(1X,2I8,F8.2,I8,F8.2,I8,F8.2,F8.0,F8.1,F8.0, 1X,A)                
61    FORMAT(1X,A)
64    FORMAT(' A=',A7)
66    FORMAT(' The values of flow in the pump capacity table are in',
     A       ' the same',/,' units as used for the design flows .',
     B       '  The body of the table gives the',/,' capacity for one',
     C       ' pump and FAC includes the adjustment for flow',/,
     D       ' units as well as number of pumps.')
C***********************************************************************
      IF(GRAV.LT.12.) THEN
C       FEQ wants cubic meters per second for flows.
        FLOW_UNITS = 'm^3/s'
      ELSE
C       FEQ wants cubic feet per second for flows.
        FLOW_UNITS ='f^3/s'
      ENDIF

c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  Global values for zone, hgrid, vdatum, and
c     unitsys are handled in subroutine TABOUT
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

      JUST = 'RIGHT'

C     Get the size of the user's flow unit in the input for flows.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
            
      READ(LINE,1) Q_UNIT

      WRITE(STDOUT,51) Q_UNIT, FLOW_UNITS
      WRITE(STDOUT,*) ' '
      IF(OPTION.EQ.'SFWMD') THEN
        WRITE(STDOUT,66)
        WRITE(STDOUT,*) ' '
      ENDIF

C     Read the heading lines. 
      DO 100 I=1,NHEAD
        CALL inline
     I            (STDIN, STDOUT,
     O             HEAD(I))
        HEAD_LEN(I) = LENSTR(HEAD(I))
100   CONTINUE

      CALL GET_ITEM_LIMITS(
     I                     STDOUT, HEAD(NHEAD), MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)


200   CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)

        IF(OPTION.EQ.'PUMPLOSS') THEN
C          READ(LINE,2) LOSS_TAB, NUMBER_OF_PUMPS, INLET_LENGTH, 
C     A                 INLET_TAB, OUTLET_LENGTH, OUTLET_TAB, 
C     B                 INLET_LOSS_FACTOR,  Q_MAX, LABEL
          CALL READ_PUMPLOSS_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         IDLENA, IDLENB, IDLENC, 
     O                         LOSS_TAB, NUMBER_OF_PUMPS, 
     O                         INLET_LENGTH, INLET_TAB, 
     O                         OUTLET_LENGTH, OUTLET_TAB, 
     O                         INLET_LOSS_FACTOR, Q_MAX, LABEL)

          N = LOSS_TAB
          M = N
        ELSEIF(OPTION.EQ.'SFWMD') THEN
C          READ(LINE,3) NET_PUMP_TAB, NUMBER_OF_PUMPS, 
C     A                INLET_LENGTH, INLET_TAB, 
C     B      OUTLET_LENGTH, OUTLET_TAB, INLET_LOSS_FACTOR, 
C     C      H_DES, Q_DES, LABEL
          CALL READ_AXIALPUMP_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         IDLENA, IDLENB, IDLENC, 
     O                         NET_PUMP_TAB, NUMBER_OF_PUMPS, 
     O                         INLET_LENGTH, INLET_TAB, 
     O                         OUTLET_LENGTH, OUTLET_TAB, 
     O                         INLET_LOSS_FACTOR, H_DES, Q_DES, LABEL)
 
         IF(NUMBER_OF_PUMPS.LE.0) NUMBER_OF_PUMPS = 1
          N = NET_PUMP_TAB
          M = N
        ELSEIF(OPTION.EQ.'PUMPADJ') THEN
          READ(LINE,4) GROSS_PUMP_TAB, NET_PUMP_TAB, INLET_LENGTH,
     A     INLET_TAB, OUTLET_LENGTH, OUTLET_TAB, 
     B     INLET_LOSS_FACTOR, LABEL
          N = GROSS_PUMP_TAB
          M = NET_PUMP_TAB
        ELSE
          WRITE(STDOUT,58) 
          STOP ' Abnormal stop. Bug found.' 
        ENDIF

        IF(N.LE.0) GOTO 300

        WRITE(STDOUT,'(/)')
        DO 210 I=1,NHEAD
          WRITE(STDOUT,52) HEAD(I)(1:HEAD_LEN(I))
210     CONTINUE

        IF(OPTION.EQ.'PUMPLOSS') THEN
C          WRITE(STDOUT,50) LOSS_TAB, NUMBER_OF_PUMPS,
C     A                     INLET_LENGTH, INLET_TAB, 
C     B                     OUTLET_LENGTH, OUTLET_TAB, 
C     C                     INLET_LOSS_FACTOR, Q_MAX, LABEL
          CALL WRITE_PUMPLOSS_ITEMS(
     I                         STDOUT, NITEM, ITEM_START,
     I                         ITEM_END,
     I                         IDLENA, IDLENB, IDLENC, 
     I                         LOSS_TAB, NUMBER_OF_PUMPS, 
     I                         INLET_LENGTH, INLET_TAB, 
     I                         OUTLET_LENGTH, OUTLET_TAB, 
     I                         INLET_LOSS_FACTOR, Q_MAX, LABEL,
     O                         LINE)
          WRITE(STDOUT,61) LINE
        ELSEIF(OPTION.EQ.'SFWMD') THEN
C          WRITE(STDOUT,60) NET_PUMP_TAB, NUMBER_OF_PUMPS,
C     A                    INLET_LENGTH, INLET_TAB, 
C     B           OUTLET_LENGTH, OUTLET_TAB, INLET_LOSS_FACTOR,
C     C           H_DES, Q_DES, LABEL
          CALL WRITE_AXIALPUMP_ITEMS(
     I                         STDOUT, NITEM, ITEM_START,
     I                         ITEM_END,
     I                         IDLENA, IDLENB, IDLENC, 
     I                         NET_PUMP_TAB, NUMBER_OF_PUMPS, 
     I                         INLET_LENGTH, INLET_TAB, 
     I                         OUTLET_LENGTH, OUTLET_TAB, 
     I                         INLET_LOSS_FACTOR, H_DES, Q_DES, LABEL,
     O                         LINE)
          WRITE(STDOUT,61) LINE
        ELSEIF(OPTION.EQ.'PUMPADJ') THEN
          WRITE(STDOUT,62) GROSS_PUMP_TAB, NET_PUMP_TAB, INLET_LENGTH, 
     A      INLET_TAB, OUTLET_LENGTH, OUTLET_TAB, 
     B      INLET_LOSS_FACTOR, LABEL
        ENDIF       

C       Get the values needed by all options.         
        EF = 0      
        IF(INLET_TAB.GT.0) THEN
          INLET_ADRS = INLET_TAB
          CALL CHKTAB
     I               (20, STDOUT, FTPNT, PMXTAB,
     M                INLET_ADRS,
     O                EF)
          IF(EF.NE.0) THEN
            TABID = GET_TABID(INLET_TAB)
            WRITE(STDOUT,54) TABID(1:LENSTR(TABID))
            EFLAG = 1
          ELSE
            D_INLET = GETD(INLET_ADRS, STDOUT)
          ENDIF
        ELSE
          INLET_ADRS = 0
        ENDIF

        EF = 0      
        IF(OUTLET_TAB.GT.0) THEN
          OUTLET_ADRS = OUTLET_TAB
          CALL CHKTAB
     I               (20, STDOUT, FTPNT, PMXTAB,
     M                OUTLET_ADRS,
     O                EF)
          IF(EF.NE.0) THEN
            TABID = GET_TABID(OUTLET_TAB)
            WRITE(STDOUT,56) TABID(1:LENSTR(TABID))
            EFLAG = 1
          ELSE
            D_OUTLET = GETD(OUTLET_ADRS, STDOUT)
            
          ENDIF
        ELSE
          OUTLET_ADRS = 0
        ENDIF

C       Does the table to be created already exist?
        IF(FTPNT(M).NE.0) CALL KIL
     I                             (25,
     M                              M, EFLAG)


        IF(INLET_LOSS_FACTOR.LT.0.0) INLET_LOSS_FACTOR = 0.0

C       Establish the physical parameters. 
        IF(INLET_ADRS.GT.0) THEN
          CALL XLKT20(INLET_ADRS,
     M                D_INLET,
     O                A_INLET, T, DT, K_INLET, DK, B, DB)
        ELSE
          A_INLET = 1.0
          K_INLET = 1.0
          INLET_LENGTH = 0.0
          INLET_LOSS_FACTOR = 0.0
        ENDIF
        IF(OUTLET_ADRS.GT.0) THEN
          CALL XLKT20(OUTLET_ADRS,
     M                D_OUTLET,
     O                A, T, DT, K_OUTLET, DK, B, DB)
          A = FLOAT(NUMBER_OF_PUMPS)*A
          EXIT_AREA = PUT7(A)
          IT = LENSTR(LABEL)
          IF(IT.LT.40) THEN
            WRITE(LABEL(IT+1:IT+11),64) EXIT_AREA
          ENDIF
        ELSE
          K_OUTLET = 1.0
          OUTLET_LENGTH = 0.0
        ENDIF

C       The loss caused by the inlet and outlet conduit and the inlet
C       entrance losses are all proportional to the square of the flow.
C       Compute the factor  to use.

        Q_FACTOR = INLET_LOSS_FACTOR/(2.*GRAV*A_INLET**2) +
     A             INLET_LENGTH/K_INLET**2  +
     B             OUTLET_LENGTH/K_OUTLET**2



        IF(EFLAG.NE.0) GOTO 300


        IF(OPTION.EQ.'PUMPLOSS') THEN
C         Compute the pump loss relationship.  Need only two values for
C         a table of type 4.  Output middle value for checking.
          ARG(1) = 0.0
          F1(1) = 0.0
          F2(1) = 0.0

          ARG(2) = 0.5*Q_MAX*Q_UNIT*FLOAT(NUMBER_OF_PUMPS)
          F1(2) = Q_FACTOR*(0.5*Q_MAX*Q_UNIT)**2
          F2(2) = Q_FACTOR*Q_MAX*Q_UNIT/
     A                            FLOAT(NUMBER_OF_PUMPS)


          ARG(3) = Q_MAX*Q_UNIT*FLOAT(NUMBER_OF_PUMPS)
          F1(3) = Q_FACTOR*(Q_MAX*Q_UNIT)**2
          F2(3) = 2*Q_FACTOR*Q_MAX*Q_UNIT/
     A                            FLOAT(NUMBER_OF_PUMPS)
          TYPE = 4
          BODY_HEAD = '      Flow Head Loss Derivativ'
          FACTOR = 1.0
          N = 3
        ELSEIF(OPTION.EQ.'SFWMD') THEN
C         The pump in question is assumed to be an axial flow pump 
C         with a standard relative performance.  The flow and head are
C         taken relative to those at the point of maximum efficiency and
C         the head and flow at that point are called the design values.

C         The pump curve we need in FEQ is not convenient for deducting
C         the losses.  Therefore we will deduct the losses using the
C         inverse of the FEQ pump curve, the one normally given in 
C         pump handbooks.    However, the values will still be stored
C         in the vectors consistent with the form that FEQ wants.
C         FEQ wants the pump curve to contain the flow in the pump
C         as a function of head across the pump.  However, when
C         we compute the effect of the losses in the conduits and
C         at the inlet conduit entrance, we will treat the 
C         relationship as the head across the pump being a functioin
C         of flow.  

C         Fill the vectors in the order required by FEQ for
C         the head across the pump.
          DO 220 I=1,17
            F1(I) = Q_DES*AXIAL_Q_REL(18-I)
            ARG(I) = H_DES*AXIAL_H_REL(18-I) - 
     A                    Q_FACTOR*(F1(I)*Q_UNIT)**2
220       CONTINUE
          N = 17
          BODY_HEAD = '      Head      Flow'
          TYPE = 2
C         Provide for conversion of GPM to CFS and adjust for the
C         number of pumps.  The body of the table will contain the
C         the flows for ONE pump so that the user can check it
C         more easily.  
          FACTOR = Q_UNIT*FLOAT(NUMBER_OF_PUMPS)
        ENDIF

C       Output the finished table.
                           
        WRITE(STDOUT,*) ' '
        CALL OUTPUT_TYPE_234
     I                   (STDOUT, STDTAB, N, M, TYPE, FACTOR, 
     I                    BODY_HEAD, LABEL, ARG, F1, F2,
     i                    zone, hgrid, vdatum, unitsys, basis,
     i                    easting, northing)


        GOTO 200

300   CONTINUE


      RETURN
      END

