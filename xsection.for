C     Program units for computing cross section tables, checking
C     cross section tables for this or that, interpolating, etc.

c
c
c

      subroutine adj_invert(stdout, npnt, stat, ws_tab, dinvert,
     m                    z,
     o                    eflag)

c     Make adjustments to the "below-water" portion of the cross section. 
c     Only makes sense for certain cross sections.  This only applies
c     in detail to the current main-channel cross sections on the 
c     Nooksack River in Whatcom County, Washington.  However, it can be 
c     used on any river if the pattern is followed.  3 June 2003.

      implicit none

      integer stdout, npnt, eflag, ws_tab
      real stat, dinvert, z(npnt)

c     Local

      integer i, il, ir, ntab
      real ws_elev, zmin, hmax, h, df
c     ***************************formats********************************
50    format(' Adjusting section below approximate time-of-flight',
     a' water surface')
52    format(' Relative change in area below water surface is about: ',
     a   f6.2)
54    format(' *ERR:XXX* DINVERT=',F8.2,' > ',' max. local depth=',f8.2)
c***********************************************************************
c     Lookup the water-surface elevation to use.
      
      CALL LKTAB
     I          (ws_tab, stat, 1,
     O           ws_elev, ntab, df)

      write(stdout,*) ' '
      write(stdout,50) 
      write(stdout,*) ' Time-of-flight surface elevation=',ws_elev

c     Now search z(*) from left to right seeking the first interval 
c     of elevation there that contains ws_elev.  The first point with
c     an elev below ws_elev is remembered as il. 

      il = -1
      do i=2,npnt
        if(ws_elev.le.z(i-1).and.ws_elev.ge.z(i)) then
          il = i
          exit
        endif
      enddo
      if(il.eq.-1) then
        write(stdout,*) ' il undefined in adj_invert'
        stop 'Abnormal stop. Error(s) found.'
      endif

c     Now search z(*) from right to left finding the first point below
c     the water surface again. 

      ir = -1
      do i=npnt,2,-1
        if(ws_elev.le.z(i).and.ws_elev.ge.z(i-1)) then
          ir = i - 1
          exit
        endif
      end do
      if(ir.eq.-1) then
        write(stdout,*) ' ir undefined in adj_invert'
        stop 'Abnormal stop. Error(s) found.'
      endif

      if(ir.le.il) then
        write(stdout,*) ' ir <= il in adj_invert'
        stop 'Abnormal stop. Error(s) found.'
      endif


c     Find the min elevation. 
      zmin = 1.e30
      do i=il,ir
        zmin = min(zmin,z(i))
      end do

      write(stdout,*) ' il=',il, ' z(il)=',z(il)
      write(stdout,*) ' ir=',ir, ' z(ir)=',z(ir)

c     Define the max local depth
      hmax = ws_elev - zmin
      if(hmax.le.0.0) then
        write(stdout,*) ' hmax <= 0.0 in adj_invert'
        stop 'Abnormal stop. Error(s) found.'
      endif

      write(stdout,*) ' Maximum local depth=',hmax
      write(stdout,52) -dinvert/hmax

      if(dinvert.ge.hmax) then
        write(stdout,54) dinvert, hmax
        stop 'Abnormal stop. Error(s) found.'
      endif        
c     Make the adjustments.

      do i=il,ir
        h = ws_elev - z(i)
        if(h.gt.0.0) then
c         Only adjust the part of the cross section that is below
c         the estimated time of flight water level.  We allow whatever
c         changes take place in the line segements that penetrate the water
c         surface.  The typical vertical increment between boundary points
c         from the DTM is so small as to negate any benefit of inserting
c         additional points to make the water surface come out exactly. 
c         Remember that the water surface is only a rough estimate, the 
c         below-water cross section is even more rough so that any noise
c         we introduce by not interpolating to the exact water surface 
c         is minor compared to the large uncertainty we have in the 
c         points on the cross section itself. 

c         The adjustment applies the relative shift at the maximum local 
c         depth to all other positive local depths.  We then get a smooth
c         transition near the water surface at the time of flight.  The 
c         below-water shape, whatever the DTM TIM makes it, is retained.
c         It is approximately a symmetrical triangle but the vagaries of the 
c         TIN consruction and the pattern of points on breaklines near the 
c         Nksk  cause there to be enough deviation from that form to make it
c         difficult to detect the water level at time of flight from the 
c         pattern of points in the cross section. 

          z(i) = z(i) + dinvert*h/hmax
        endif
      end do

      return
      end
      
C
C
C
      SUBROUTINE  PROCESS_FLNAME(STDOUT, FLNAME, FLFILE)

C     Find and input information on a flow line to be used in XSINTERP
C     to define the Easting and Northing fields in an interpolated 
C     cross section. 

      IMPLICIT NONE
      INTEGER STDOUT
      CHARACTER FLNAME*6, FLFILE*64

      INCLUDE 'flowline.cmn'

C     Called program units
      INTEGER GET_UNIT
      EXTERNAL GET_UNIT, FREE_UNIT
C     Local

      INTEGER I, FLUNIT
      LOGICAL THERE
      REAL*8 SUM
      CHARACTER LINE*80
C     ******************************FORMATS*****************************
50    FORMAT(/,'*ERR:761* Flow-line Name= ',A,
     A' not found in file named:',A)
52    FORMAT(/,' RM_ORIGIN=',F10.4,' S_ATRMORG=',F10.1,
     A         ' RM_FACTOR=',F10.1,' S_BEGIN=',F10.1)
C***********************************************************************
      FLUNIT = GET_UNIT(0)
      INQUIRE(FILE=FLFILE, EXIST=THERE)
      IF(THERE) THEN
        OPEN(UNIT=FLUNIT, FILE=FLFILE, STATUS='OLD')
      ELSE
        WRITE(STDOUT,*) ' FILE named:',FLFILE,' not found.'
        WRITE(STDOUT,*) ' check flow-line file (FLFILE). '
        STOP 'Abnormal stop. Errors found.'
      ENDIF

      FL_EASTING = 0.D0
      FL_NORTHING = 0.D0

C     Scan the file seeking the FLNAME starting in column 1. 
100   CONTINUE
        CALL inline
     I            (FLUNIT, STDOUT,
     O             LINE)
        IF(LINE(1:6).NE.'FINISH') THEN
          IF(LINE(1:6).EQ.FLNAME) THEN
C           Found the flow-line name.  Now look for 
C           RM_ORIGIN starting in column 2
200         CONTINUE
              CALL inline
     I                  (FLUNIT, STDOUT,
     O                   LINE)
              IF(LINE(2:10).EQ.'RM_ORIGIN') THEN
C               Now we can begin to read the information. 
                CALL inline
     I                    (FLUNIT, STDOUT,
     O                     LINE)
                READ(LINE,'(F10.0,F10.0,F10.0, F10.0)') 
     A          RM_ORIGIN, S_ATRMORG, RM_FACTOR, S_BEGIN
                WRITE(STDOUT,52) RM_ORIGIN, S_ATRMORG, RM_FACTOR, 
     A                         S_BEGIN


                FL_PNT_KNT = 1
                SUM = S_BEGIN
                CALL inline
     I                    (FLUNIT, STDOUT,
     O                     LINE)
                READ(LINE,'(12X,2F13.0)') FL_EASTING(FL_PNT_KNT), 
     A                                    FL_NORTHING(FL_PNT_KNT)
                FL_DISTANCE(FL_PNT_KNT) = RM_ORIGIN + 
     A                  (SUM - S_ATRMORG)/RM_FACTOR
300             CONTINUE
                  CALL inline
     I                      (FLUNIT, STDOUT,
     O                       LINE)
                  IF(LINE(1:3).NE.'END') THEN
                    FL_PNT_KNT = FL_PNT_KNT + 1
              
                    READ(LINE,'(12X,2F13.0)') FL_EASTING(FL_PNT_KNT), 
     A                                        FL_NORTHING(FL_PNT_KNT)

                    SUM = SUM + SQRT(
     A        (FL_EASTING(FL_PNT_KNT) - FL_EASTING(FL_PNT_KNT-1))**2
     B             +
     C        (FL_NORTHING(FL_PNT_KNT) - FL_NORTHING(FL_PNT_KNT-1))**2)
                    FL_DISTANCE(FL_PNT_KNT) = RM_ORIGIN + 
     A                                  (SUM - S_ATRMORG)/RM_FACTOR
                    GOTO 300
                  ELSE
                    GOTO 400
                  ENDIF
              ELSE
                GOTO 200
              ENDIF
          ELSE
            GOTO 100
          ENDIF
        ELSE   
          WRITE(STDOUT,50) FLNAME, FLFILE
          STOP 'Abnormal stop,  Error found.'
        ENDIF

400   CONTINUE


      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,*) 'Debug dump for x,y location'
      DO 500 I=1,FL_PNT_KNT
        WRITE(STDOUT,60) FL_DISTANCE(I), FL_EASTING(I), FL_NORTHING(I)
60    FORMAT(F10.4,F13.2,F13.2)
500   CONTINUE
      CALL FREE_UNIT(STDOUT, FLUNIT)
      WRITE(STDOUT,*) ' '

      FL_PRESENT = 1

      RETURN
      END



C     ***********
C     *         *
C     * CHK_FOR_GIVEN_TYPE
C     *         *
C     ***********

      SUBROUTINE CHK_FOR_GIVEN_TYPE(IS, IE, TYPE_VALUES, TYPE,
     M                              EFLAG)

C     Check that all types at and between indices IS and IE in the
C     TYPE_VAULES(*) match TYPE.

      IMPLICIT NONE
      INTEGER IS, IE, TYPE_VALUES(IE), TYPE, EFLAG

C     Local

      INTEGER I
C***********************************************************************
      DO 100 I=IS,IE
        IF(TYPE_VALUES(I).NE.TYPE) THEN
          EFLAG = EFLAG + 1
        ENDIF
100   CONTINUE
      RETURN
      END
      
     
C     ***********
C     *         *
C     * GET_XSEC_HEADER
C     *         *
C     ***********

      SUBROUTINE GET_XSEC_HEADER(STDIN, STDOUT,
     M   EFLAG,
     O   TAB, SAVOPT, OUTOPT, MONTON, BETOPT, EXTEND, 
     O   STAT, LEFT, RIGHT, VARN, NAVM, SCALE, SHIFT,
     O   VSCALE, HSHIFT, NSUB, N, zone, hgrid, vdatum,
     O   unitsys, basis)

      IMPLICIT NONE
C     Get the header information for the FEQX and FEQXEXT
C     commands for processing cross sections.  

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG, TAB, EXTEND, NSUB, NAVM

      REAL STAT, LEFT, RIGHT, SCALE, SHIFT, VSCALE, HSHIFT,
     A     N(PMXSUB)

      CHARACTER SAVOPT*8, OUTOPT*8, MONTON*8, BETOPT*8,
     A          VARN*4, zone*8, hgrid*8, vdatum*8, unitsys*8,
     b          basis*8


      INCLUDE 'xtadd.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  NVAL, INTVAL, REAVAL, CONTINUATION_VALUE,
     A         CHRVAL, DPRVAL, EXACT_TYPE, LOWER_TYPE,
     B         N_SYMBOL
      PARAMETER(N_SYMBOL=62, NVAL=PMXSUB, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, CONTINUATION_VALUE=5,
     B          EXACT_TYPE=0,LOWER_TYPE=1)

      INTEGER EFLAG2, I, IS, IE, IT, J, KNT, OPT, IP, SELECTION,
     A                ITEM_KNT
      INTEGER CLEN(NVAL), IVAL(NVAL), TERML(NVAL), TERMCLS(NVAL),
     A        ITEM_TYPE(NVAL)
      REAL RVAL(NVAL)
      REAL*8 DPVAL(NVAL), dnull
      CHARACTER CVAL(NVAL)*256, TERM(NVAL)*1, LINE*120,
     A          KEY*16

      INTEGER LENSTR

      EXTERNAL LENSTR, GET_MULTIPLE_REAL_VALUES, 
     A   GET_INTERNAL_TAB_NUMBER

C     + + + SAVED VALUES + + +
      INTEGER SYMBOL_VALUE(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE SYMBOL_VALUE, SYMBOL_TABLE


      DATA  SYMBOL_TABLE
     1      /'EXTEND  ','MONOTONE','NEWBETA ','NEWBETAE','NEWBETAM',
     A       'NOEXTEND','NOOUT   ','NOSAVE  ','OLDBETA ','OUT1    ',
     B       'OUT12   ','OUT20   ','OUT21   ','OUT22   ','OUT23   ',
     C       'OUT24   ','OUT25   ','SAVE    ','SAVE1   ','SAVE12  ',
     D       'SAVE20  ','SAVE21  ','SAVE22  ','SAVE23  ','SAVE24  ',
     E       'SAVE25  ','TABLE   ','STATION ','LEFT    ','RIGHT   ',
     F       'VARN    ','NAVM    ','SCALE   ','SHIFT   ','VSCALE  ',
     G       'HSHIFT  ','NSUB    ','GISID   ','GIS     ','EASTING',
     H       'NORTHING','TABID','NEWBETAX',   'OUT30   ','OUT31  ',
     I       'OUT32   ','OUT33   ','OUT34   ','OUT35   ','SAVE30',
     J       'SAVE31  ','SAVE32  ','SAVE33  ','SAVE34  ','SAVE35',
     k       'DINVERT ','WS_TABID','ZONE    ','HGRID   ','VDATUM',
     l       'UNITSYS ','BASIS'/

      DATA SYMBOL_VALUE
     1      /         6,         5,         4,         4,         4,
     A                7,         3,         1,         4,         3,
     B                3,         3,         3,         3,         3,
     C                3,         3,         2,         1,         1,
     D                1,         1,         1,         1,         1,
     E                1,         8,         9,        10,        11,
     F               12,        13,        14,        15,        16,
     G               17,        18,        19,        19,        20,
     H               21,         8,         4,         3,         3,
     I                3,         3,         3,         3,         1,
     J                1,         1,         1,         1,         1,
     k               22,        23,        24,        25,        26,
     l               27,        28/

      DATA RESPONSE_TYPE
     A     /26*0, CHRVAL, REAVAL, REAVAL, REAVAL, CHRVAL, INTVAL,
     B            REAVAL, REAVAL, REAVAL, REAVAL, -1, CHRVAL, CHRVAL,
     C            DPRVAL, DPRVAL, CHRVAL, 0, 12*0, REAVAL, CHRVAL,
     e            CHRVAL, CHRVAL, CHRVAL, CHRVAL, CHRVAL/
      DATA CONVERT_RULE
     A     /26*0, LOWER_TYPE, 3*LOWER_TYPE,EXACT_TYPE, EXACT_TYPE,
     B          4*LOWER_TYPE, 3*EXACT_TYPE,2*LOWER_TYPE, LOWER_TYPE, 0,
     C          12*0, LOWER_TYPE, 6*LOWER_TYPE/
      data dnull/-33d6/
C     *****************************FORMATS******************************
 52   FORMAT(/,' *BUG:XXX* Invalid index=',I5,' for name=',A,' in',
     A       ' subroutine GET_XSEC_HEADER.')
54    FORMAT(/,' *ERR:734* Name=',A8,' is unknown in a cross section',
     A   ' description.')
56    FORMAT(/,' Unable to continue due to previous errors.')
 58   FORMAT(' *ERR:504* Number of subsections=',I5,' > ',I5)
60    FORMAT(/,' Processing:',A)
62    FORMAT(/,' Seeking additional roughness values.')
C***********************************************************************
C     Clear the local error flag for subroutine GETVAL and GET_MULTIPLE
C     _REAL_VALUES
      EFLAG2 = 0
C     Read lines of input and process each one until the expected number
C     of lines or an end of header signal is found.  


C     SET DEFAULTS
 
      SAVOPT = 'NOSAVE'
      OUTOPT = 'OUT21'
      MONTON = ' '
      BETOPT = 'OLDBETA'
      LEFT = 1.E30
      RIGHT = -1.E30
      SCALE = 1.0
      VSCALE = 1.0
      SHIFT = 0.0
      HSHIFT = 0.0
      NAVM = 0
      VARN = 'NCON'
      STAT = 0.0
      GISID = ' '
      TABID = ' '
      TAB = 0
      EASTING = dnull
      NORTHING = dnull

      dinvert = 0.0
      WS_TABID= ' '
      zone = 'NONE'
      hgrid = 'NONE'
      vdatum = 'NONE'
      unitsys = 'NONE'
      basis = 'NONE'

C     Select option for GETVAL to return BOTH the string value
C     and the converted value for numeric responses.  
C     In some cases the value will be recomputed. 
      OPT = 2
C     Start a loop over input lines
100   CONTINUE


        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)

        IT = LENSTR(LINE)
        WRITE(STDOUT,60) LINE(1:IT)

        LINE(IT+1:IT+1) = ''''

        CALL GETVAL
     I             (STDOUT, LINE, NVAL, OPT,
     O              ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O              EFLAG2, TERM, TERML, TERMCLS,
     O              ITEM_KNT)
C      WRITE(STDOUT,*) 
C     A  ' Return from GETVAL in GET_XSEC_HEADER: ITEM_KNT=',ITEM_KNT
C      WRITE(STDOUT,97)
C97    FORMAT(1X,12X,'ITEM','   LEN  TYPE T  TCLS')
C      DO 9213 I=1,ITEM_KNT
C        WRITE(STDOUT,99) CVAL(I), CLEN(I), ITEM_TYPE(I), TERM(I),
C     A                  TERMCLS(I)
C99    FORMAT(' ',A16,' ',I5,' ',I5,' ',A1,' ',I5)
C9213  CONTINUE

        IF(EFLAG2.NE.0) THEN
C         Error in parsing the line of input.
          
          WRITE(STDOUT,56) 
          STOP 'Abnormal stop.  Errors found.' 
        ELSE
C         No errors reported.  Process the items found on the current
C         line. 

          CALL CHK_AND_CONVERT_RESPONSE(STDOUT,
     I           N_SYMBOL, SYMBOL_TABLE, RESPONSE_TYPE,
     I           CONVERT_RULE, ITEM_KNT, ITEM_TYPE, CVAL, CLEN,
     M           EFLAG,
     O           IVAL, RVAL, DPVAL)

          I = 1
110       CONTINUE          
C           Find the value for the next item from the 
C           symbol table. 
            KEY = CVAL(I)(1:CLEN(I))
            CALL LSTAB
     I                (KEY, SYMBOL_TABLE, N_SYMBOL,
     O                 IP)
            IF(IP.EQ.0) THEN
C             error-symbol not found
              WRITE(STDOUT,54) KEY
              EFLAG = 1

            ELSE
              SELECTION = SYMBOL_VALUE(IP)            
            ENDIF

            GOTO(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
     A           14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
     b           25, 26, 27, 28),SELECTION

              WRITE(STDOUT,52) I, KEY
              STOP 'Abnormal stop. Errors found.'
 
 1          CONTINUE
C             MOST OF THE SAVE OPTIONS.
              SAVOPT = KEY
              I = I + 1
              GOTO 149
 2          CONTINUE
C             SAVE BY ITSELF WITHOUT A SUFFIX NUMBER
              SAVOPT = 'SAVE21'
              I = I + 1
              GOTO 149
 3          CONTINUE
C             OUTPUT OPTIONS
              OUTOPT = KEY
              I = I + 1
              GOTO 149
 4          CONTINUE
C             VELOCITY COEFFICIENT FLAGS
              BETOPT = KEY
              I = I + 1
              GOTO 149
 5          CONTINUE
C             MONOTONICITY FLAG
              MONTON = KEY
              I = I + 1
              GOTO 149
 6          CONTINUE
C             SET EXTEND OPTION
              EXTEND = 1
              I = I + 1
              GOTO 149
 7          CONTINUE
C             CLEAR THE EXTEND OPTION
              EXTEND = 0
              I = I + 1
              GOTO 149
 8          CONTINUE
C             Set the value for the TABLE id.  
              TABID = CVAL(I+1)(1:CLEN(I+1))
              CALL GET_INTERNAL_TAB_NUMBER
     I                                    (STDOUT, TABID,
     M                                     EFLAG,
     O                                     TAB)

              I = I + 2
              GOTO 149
 9          CONTINUE
              STAT = RVAL(I+1)
              I = I + 2
              GOTO 149
10          CONTINUE
              LEFT = RVAL(I+1)
              I = I + 2
              GOTO 149
11          CONTINUE
              RIGHT = RVAL(I+1)
              I = I + 2
              GOTO 149
12          CONTINUE
              VARN = CVAL(I+1)(1:CLEN(I))
              I = I + 2
              GOTO 149
13          CONTINUE
              NAVM = IVAL(I+1)
              I = I + 2
              GOTO 149
14          CONTINUE
              SCALE = RVAL(I+1)
              I = I + 2
              GOTO 149
15          CONTINUE
              SHIFT = RVAL(I+1)
              I = I + 2
              GOTO 149
16          CONTINUE
              VSCALE = RVAL(I+1)
              I = I + 2
              GOTO 149
17          CONTINUE
              HSHIFT = RVAL(I+1)
              I = I + 2
              GOTO 149
18          CONTINUE
C             Process the number of subsections and the values of Manning's 
C             n for each.  NSUB must be the last entry in the header. 
C             There are two options: the number of subsections is given
C             or the number of subsections is omitted.  When the number of 
C             subsections is omitted, the continuation signal must be given
C             for all lines required for the roughness values but for the 
C             last. 

              CALL  GET_MULTIPLE_REAL_VALUES(
     I          STDIN, STDOUT, INTVAL, REAVAL, CONTINUATION_VALUE, 
     I          OPT, NVAL, 
     M          LINE, ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN,
     M          EFLAG2, TERM, TERML, TERMCLS, ITEM_KNT,
     O          NSUB, N)

              GOTO 200
19          CONTINUE
              GISID = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149
20          CONTINUE
              EASTING = DPVAL(I+1)
              I = I + 2
              GOTO 149
21          CONTINUE
              NORTHING = DPVAL(I+1)
              I = I + 2
              GOTO 149
22          CONTINUE
              dinvert = rVAL(I+1)
              I = I + 2
              GOTO 149
23          CONTINUE
              WS_TABID = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149
24          CONTINUE
              zone = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149
25          CONTINUE
              hgrid = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149
26          CONTINUE
              vdatum = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149
27          CONTINUE
              unitsys = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149
28          CONTINUE
              basis = CVAL(I+1)(1:CLEN(I+1))
              I = I + 2
              GOTO 149

 149        CONTINUE
              IF(I.GT.ITEM_KNT) THEN
C               Get the next line from the input
                GOTO 100
              ELSE
C               Get the next item from the current line
                  GOTO 110
              ENDIF            
        ENDIF

200   CONTINUE
      RETURN
      END

C
C
C
      SUBROUTINE   EXTRAP
     I                   (TYPE, XOFF, MLFT,
     M                    EFLAG, EXT, FTP)
 
C     + + + PURPOSE + + +
C     Extend a cross-section table by extrapolation
 
C     PARAMTERS
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, MLFT, TYPE, XOFF
      REAL EXT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     TYPE   - Table type
C     XOFF   - Offset between successive depth values for cross section
C               function table
C     MLFT   - maximum length of FTAB/ITAB
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     EXT    - extent and nature of extrapolation of the table
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'stdun.cmn'
      INCLUDE 'grvcom.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER JID(25), QCID(25)
      SAVE JID, QCID
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IOFF, it
      REAL A2, A3, ALP2, ALP3, B2, B3, DA, DAEXT, DALP, DALPDA, DB,
     A     DBDA, DKH, DKHDA, DMA, DMADA, DMQ, DMQDA, DQC, DQCDA, DT,
     B     DTDY, DY, KH2, KH3, MA2, MA3, MQ2, MQ3, QC2, QC3, T2, T3, Y2,
     C     Y3, YB2, YB3, MXSLOT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL KIL
 
C     + + + DATA INITIALIZATIONS + + +
      DATA  JID/1,10*0, 1,8*0, 1 ,1, 0, 1 ,1/
      DATA QCID/11*0, 1,9*0, 1, 0, 0, 1/
 
C     + + + OUTPUT FORMATS + + +
 1    FORMAT(/,' *WRN:01* EXTRAPOLATION NOT DONE. TOP WIDTH',
     A      ' < 0. EXT = ',F10.2)
 50   FORMAT(/,'  *WRN:571 Extrapolation not done.',
     A           '  Section is slotted.')
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        MXSLOT = 0.07
      ELSE
        MXSLOT = 0.02134
      ENDIF
C     ON ENTRY FTP POINTS TO THE NEXT OPEN LOCATION IN
C     (FTAB,ITAB).  ON EXIT FTP SHOULD POINT TO THE NEW
C     NEXT OPEN LOCATION IN (FTAB,ITAB).
 
C     DO THE VALUES THAT ARE COMMON TO ALL OF THE 14 TABLE TYPES.
C     THEN BRANCH TO SPECIAL SECTIONS TO FINISH TABLES THAT
C     HAVE MORE TABULATED VALUES IN THEM THAN THE BARE MINIMUM
 
C     GET THE LAST LEVEL IN THE TABLE JUST ADDED TO FTAB/ITAB
 
      Y2 = FTAB(FTP-XOFF)
      T2 = FTAB(FTP-XOFF+1)
      A2 = FTAB(FTP-XOFF+2)
      KH2 = FTAB(FTP-XOFF+3)
      B2 = FTAB(FTP-XOFF+4)
 
C     Do not do extrapolation on slotted cross section tables.
      IF(T2.LT.MXSLOT) THEN
        WRITE(STD6,50)
        RETURN
      ENDIF
C     COMPUTE THE RATE OF CHANGE OF T, KH AND B
 
      IOFF = 2*XOFF
      DY = Y2 - FTAB(FTP-IOFF)
      DT = T2 - FTAB(FTP-IOFF+1)
      DA = A2 - FTAB(FTP-IOFF+2)
      DKH = KH2 - FTAB(FTP-IOFF+3)
      DB = B2 - FTAB(FTP-IOFF+4)
      DTDY = DT/DY
      DKHDA = DKH/DA
      DBDA = DB/DA
 
C     CHECK FOR MODE OF TOP WIDTH EXTRAPOLATION.
 
      IF(EXT.LT.0) DTDY = 0.0
      EXT = ABS(EXT)
 
C     NOW COMPUTE THE NEW LEVEL
 
      Y3 = Y2 + EXT
      T3 = T2 + DTDY*EXT
      IF(T3.LE.0.0) THEN
C       EXTRAPOLATION OF TOP WIDTH LEADS TO NON-POSITIVE VALUE.
        WRITE(STD6,1) EXT
        RETURN
      ENDIF
 
      DAEXT = 0.5*EXT*(T2 + T3)
      A3 = A2 + DAEXT
      KH3 = KH2 + DKHDA*DAEXT
      B3 = B2 + DBDA*DAEXT
      IF(B3.LT.1.0) B3 = 1.0
 
C     NOW STORE AWAY THE VALUES
 
      IF((FTP+XOFF).GE.MLFT) CALL KIL
     I                                (10,
     M                                 FTP, EFLAG)
      FTAB(FTP) = Y3
      FTAB(FTP+1) = T3
      FTAB(FTP+2) = A3
      FTAB(FTP+3) = KH3
      FTAB(FTP+4) = B3
 
 
C     THE CORE SET HAS BEEN DONE.  NOW DO THE TABLES THAT CONTAIN
C     ADDITIONAL ELEMENTS
 
      if(type.gt.25) then
        it = type - 10
      else
        it = type
      endif
      IF(JID(it).EQ.1) THEN
C       FIRST MOMENT EXISTS IN THESE TABLES.
        YB2 = FTAB(FTP-XOFF+5)
        YB3 = YB2 + 0.5*EXT*(A2 + A3) - EXT*EXT*(T3 - T2)/12.
        FTAB(FTP+5) = YB3
      ENDIF
 
      IF(QCID(it).EQ.1) THEN
C       ALP AND QC EXISTS IN THESE TABLES.
        ALP2 = FTAB(FTP-XOFF+6)
        QC2 = FTAB(FTP-XOFF+7)
        DALP =  ALP2 - FTAB(FTP-IOFF+6)
        DQC = QC2 - FTAB(FTP-IOFF+7)
        DALPDA = DALP/DA
        DQCDA = DQC/DA
        ALP3 = ALP2 + DALPDA*DAEXT
        QC3 = QC2 + DQCDA*DAEXT
        IF(QC3.LT.QC2) QC3 = QC2
        IF(ALP3.LT.1.0) ALP3 = 1.0
        FTAB(FTP+6) = ALP3
        FTAB(FTP+7) = QC3
      ENDIF
 
      IF(it.GE.23) THEN
C       WE HAVE THE WEIGHT COEFFICIENTS TO EXTRAPOLATE.  THESE
C       VARY IN POSITION AND MUST BE DONE INDIVIDUALLY.  COMPUTE
C       THE RATES OF CHANGE FOR EACH TYPE.  THEN EXTRAPOLATE
C       AND THEN STORE FOR EACH TYPE
 
        IF(it.EQ.23) THEN
C         MA AND MQ ARE AT OFFSETS OF 5 AND 6
          MA2 = FTAB(FTP-XOFF+5)
          DMA = MA2 - FTAB(FTP-IOFF+5)
          DMADA = DMA/DA
          MA3 = MA2 + DMADA*DAEXT
 
          MQ2 = FTAB(FTP-XOFF+6)
          DMQ = MQ2 - FTAB(FTP-IOFF+6)
          DMQDA = DMQ/DA
          MQ3 = MQ2 + DMQDA*DAEXT
          FTAB(FTP+5) = MA3
          FTAB(FTP+6) = MQ3
        ELSEIF(it.EQ.24) THEN
C         MA AND MQ ARE AT OFFSETS OF 6 AND 7
          MA2 = FTAB(FTP-XOFF+6)
          DMA = MA2 - FTAB(FTP-IOFF+6)
          DMADA = DMA/DA
          MA3 = MA2 + DMADA*DAEXT
 
          MQ2 = FTAB(FTP-XOFF+7)
          DMQ = MQ2 - FTAB(FTP-IOFF+7)
          DMQDA = DMQ/DA
          MQ3 = MQ2 + DMQDA*DAEXT
          FTAB(FTP+6) = MA3
          FTAB(FTP+7) = MQ3
        ELSEIF(it.EQ.25) THEN
C         MA AND MQ ARE AT OFFSETS OF 8 AND 9
          MA2 = FTAB(FTP-XOFF+8)
          DMA = MA2 - FTAB(FTP-IOFF+8)
          DMADA = DMA/DA
          MA3 = MA2 + DMADA*DAEXT
 
          MQ2 = FTAB(FTP-XOFF+9)
          DMQ = MQ2 - FTAB(FTP-IOFF+9)
          DMQDA = DMQ/DA
          MQ3 = MQ2 + DMQDA*DAEXT
          FTAB(FTP+8) = MA3
          FTAB(FTP+9) = MQ3
        ENDIF
      ENDIF
 
c     Make adjustments to the derivatives wrt depth if 30 <= type <=35
      if(type.eq.30) then
        ftab(ftp+5) = (kh3 - kh2)/(y3 - y2)
        ftab(ftp+6) = (b3 - b2)/(y3 - y2)
      elseif(type.eq.31) then
        ftab(ftp+6) = (kh3 - kh2)/(y3 - y2)
        ftab(ftp+7) = (b3 - b2)/(y3 - y2)
      elseif(type.eq.32) then
        ftab(ftp+8) = (kh3 - kh2)/(y3 - y2)
        ftab(ftp+9) = (b3 - b2)/(y3 - y2)
        ftab(ftp+10) = (alp3 - alp2)/(y3 - y2)
      elseif(type.eq.33) then
        ftab(ftp+7) = (kh3 - kh2)/(y3 - y2)
        ftab(ftp+8) = (b3 - b2)/(y3 - y2)
        ftab(ftp+9) = (ma3 - ma2)/(y3 - y2)
        ftab(ftp+10) = (mq3 - mq2)/(y3 - y2)
      elseif(type.eq.34) then
        ftab(ftp+8) = (kh3 - kh2)/(y3 - y2)
        ftab(ftp+9) = (b3 - b2)/(y3 - y2)
        ftab(ftp+10) = (ma3 - ma2)/(y3 - y2)
        ftab(ftp+11) = (mq3 - mq2)/(y3 - y2)
      elseif(type.eq.35) then
        ftab(ftp+10) = (kh3 - kh2)/(y3 - y2)
        ftab(ftp+11) = (b3 - b2)/(y3 - y2)  
        ftab(ftp+12) = (ma3 - ma2)/(y3 - y2)
        ftab(ftp+13) = (mq3 - mq2)/(y3 - y2)
        ftab(ftp+14) = (alp3 - alp2)/(y3 - y2)
      endif

                         



      FTP = FTP + XOFF
 
      RETURN
 
      END
C
C
C
      REAL FUNCTION   FNDYDN
     I                      (DEPTH, M, YVEC, NVEC)
 
C     + + + PURPOSE + + +
C     Find depth-dependent value of Manning's n.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER M
      REAL DEPTH, NVEC(M), YVEC(M)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     DEPTH  - Depth argument for variation of Manning's n with
C               depth(maximum depth or hydraulic depth)
C     M      - Number of defined depth-dependent Manning's n values
C     YVEC   - Depth values for depth dependent Manning's n values
C     NVEC   - Values of Manning's n at each depth
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
C***********************************************************************
C     FIND THE EXTREME VALUES. TABLE UNDERFLOW AND TABLE OVERFLOW ARE
C     NOT ERRORS.  WE ASSIGN THE VALUE AT DEPTH 0.0 TO TABLE UNDERFLOW
C     AND VALUE AT YVEC(M) TO TABLE OVERFLOW.
 
      IF(DEPTH.LE.0.) THEN
        FNDYDN = NVEC(1)
      ELSEIF(DEPTH.GE.YVEC(M)) THEN
        FNDYDN = NVEC(M)
      ELSE
C       DO LINEAR SEARCH FOR INTERVAL AND LINEAR INTERPOLATE
        I = 1
 100    CONTINUE
          IF(DEPTH.LE.YVEC(I+1)) THEN
C           FOUND INTERVAL.
 
            FNDYDN = NVEC(I) + (DEPTH - YVEC(I))
     A               *(NVEC(I+1) - NVEC(I))/(YVEC(I+1) - YVEC(I))
          ELSE
            I = I + 1
            GOTO 100
          ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   CUTTAB
     M                   (NDEP, XST)
 
C     + + + PURPOSE + + +
C     Cut off the part of the table that is in the slot and replace it
C     with a standard slot.
 
      IMPLICIT NONE
C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NDEP
      REAL XST(PMXPNT,PMXELM)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     NDEP   - Number of depth values
C     XST    - Storage table for various elements of cross section
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
      REAL DY, NEWT, OLDT, MXSLOT, MXARG
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        MXSLOT = 0.07
        MXARG = 500.0
      ELSE
        MXSLOT = 0.02134
        MXARG = 150.0
      ENDIF
      OLDT = XST(1,2)
      DO 100 J=2,NDEP
 
        NEWT = XST(J,2)
C       Check for old bridge opening tables that do not have a slot but
C       that close at the top.  Skip these.
        IF(NEWT.EQ.0.0) GOTO 100
C        IF(NEWT.EQ.OLDT.AND.NEWT.LE.MXSLOT) THEN
        IF(ABS(NEWT-OLDT)/NEWT.LE.1.E-3.AND.NEWT.LE.MXSLOT) THEN
C         ASSUME WE ARE IN THE SLOT. KEEP CONVEYANCE, ALPHA,
C         AND BETA CONSTANT. COMPUTE AREA AND FIRST MOMENT
C         OF AREA TO NEW DEPTH OF MXARG.
 
          XST(J,1) = MXARG
          XST(J,5) = XST(J-1,5)
          XST(J,6) = XST(J-1,6)
          XST(J,7) = XST(J-1,7)
          XST(J,13) = XST(J-1,13)
          DY = XST(J,1) - XST(J-1,1)
          XST(J,3) = XST(J-1,3) + 0.5*DY*(NEWT + OLDT)
          XST(J,4) = XST(J-1,4) + 0.5*DY*(XST(J,3) + XST(J-1,3))
          NDEP = J
          RETURN
        ENDIF
        OLDT = NEWT
 100  CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   CHKTAJ
     I                   (STDOUT, NDEP, XST,
     O                    FLAG)
 
C     + + + PURPOSE + + +
C     Check consistency of top width, area, and first moment of area
C     for a cross section.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG, NDEP, STDOUT
      REAL XST(PMXPNT,PMXELM)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NDEP   - Number of depth values
C     XST    - Storage table for various elements of cross section
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      DOUBLE PRECISION AL, AR, DA, DIFF, DJ, DY, TL, TR, YL, YR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* Change in area between depths ',F10.4,' and ',
     A                   F10.4,' is inconsistent.  Rerr=',1PE10.4)
 52   FORMAT(' *BUG:XXX* Change in first moment between depths ',
     A        F10.4,' and ',F10.4,' is inconsistent.  Rerr=',1PE10.4)
C***********************************************************************
      FLAG = 0
      YL = 0.0
      TL = XST(1,2)
      AL = 0.D0
      DO 100 I=2,NDEP
        YR = XST(I,1)
        TR = XST(I,2)
        DY = YR - YL
        DA = 0.5*DY*(TL + TR)
        IF(DA.EQ.0.0) THEN
          WRITE(STDOUT,*) ' PROBLEM IN CHKTAJ: DA = 0.0'
          WRITE(STDOUT,*) ' YL=',YL,' YR=',YR,' DY=',DY
          WRITE(STDOUT,*) ' TL=',TL,' TR=',TR
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        AR = AL + DA
        DJ = 0.5*(DY*(AR + AL - DY*(TR - TL)/6.0))
C       TEST FOR CLOSENESS
 
        DIFF = ABS(DA - (XST(I,3) - XST(I-1,3)))
        IF(DIFF/DA.GT.5.E-2.AND.DA.GT.5.E-1.AND.TL.GT.2.0) THEN
C         FLAG DISCREPANCY
          FLAG =1
          WRITE(STDOUT,50) YL, YR, DIFF/DA
        ENDIF
 
        DIFF = ABS(DJ - (XST(I,4) - XST(I-1,4)))
        IF(DIFF/DJ.GT.8.E-2.AND.DJ.GT.5.E-1.AND.TL.GT.2.0) THEN
C         FLAG DISCREPANCY
          FLAG = 1
          WRITE(STDOUT,52) YL, YR, DIFF/DJ
        ENDIF
 
        YL = YR
        TL = TR
        AL = AR
 100  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   SUBSET
     I                   (STDOUT, NPNT, X, Z, SB, ZMAX, SN, LSN, 
     M                    XL, XR, EFLAG, NPNTS,
     O                    XS, ZS, SBS, SNS, LSNS)
 
C     + + + PURPOSE + + +
C     Extract the part of a cross section contained between
C     offsets XL and XR.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NPNT, NPNTS, STDOUT
      INTEGER SB(NPNT), SBS(NPNTS)
      REAL LSN(NPNT), LSNS(NPNTS), SN(NPNT), SNS(NPNTS), X(NPNT), XL,
     A    XR, XS(NPNTS), Z(NPNT), ZMAX, ZS(NPNTS)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NPNT   - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     ZMAX   - Maximum elevation
C     SN     - Sinuousity at a point on a cross section boundary
C     LSN    - Line segment Manning's n value
C     XL     - Offset at left hand end of segment
C     XR     - Offset at right hand end of segment
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NPNTS  - Number of points on boundary of a cross section
C     XS     - Offsets for subset of a cross section
C     ZS     - Elevation of boundary points in subset
C     SBS    - Subsection numbers for the line segments in the subset
C     SNS    - Sinuousity at a point on a subset of a cross section
C     LSNS   - Line segment Manning's n value for the subset
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IL, ILP1, IR, IRM1
      REAL SNL, SNR, ZL, ZR
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT('0*ERR:518* Subset and section are disjoint')
 52   FORMAT('0*WRN:524* Left-hand subset request of',F10.2,' is left',
     A   /,10X,'of cross section boundary of',F10.2)
 53   FORMAT(10X,'Subset request set to the cross section boundary')
 54   FORMAT('0*WRN:525* Right hand subset request of',F10.2,' is',
     B   ' right',/,10X,'of cross section boundary of',F10.2)
C***********************************************************************
C      WRITE(STDOUT,*) ' '
C      WRITE(STDOUT,*) ' Entering SUBSET'
C      WRITE(STDOUT,*) ' XL=',XL,' XR=',XR

C     CHECK FOR XL AND XR BEING OUTSIDE THE BOUNDS OF THE
C     CROSS SECTION
 
      IF(XL.LT.X(1)) THEN
C       ISSUE WARNING AND RESET TO THE LIMIT OF THE CROSS SECTION
 
        WRITE(STDOUT,52) XL, X(1)
        WRITE(STDOUT,53)
        XL = X(1)
      ENDIF
 
      IF(XR.GT.X(NPNT)) THEN
        WRITE(STDOUT,54) XR, X(NPNT)
        WRITE(STDOUT,53)
        XR = X(NPNT)
      ENDIF
      IF(XL.LT.X(NPNT).AND.XR.GT.X(1)) GOTO 100
        WRITE(STDOUT,51)
        EFLAG = EFLAG + 1
        RETURN
 100  CONTINUE
 

C     IL is the index of the point that is at or to the left 
C     of XL.  In other words the line segment on the cross section
C     boundary between IL  and IL + 1 will contain XL.

C     IR is the index of the point that is at or to the right
C     of XR.  In other words the line segment on the cross section
C     boundary between IR - 1 and IR will contain XR.  
      IL=0
      IR=0
      DO 200 I=2,NPNT
        IF(XL.GE.X(I-1).AND.XL.LT.X(I)) IL = I-1
        IF(XR.LE.X(I).AND.XR.GT.X(I-1)) IR = I
        IF(IR.GT.0.AND.IL.GT.0) GOTO 250
 200    CONTINUE
 250  CONTINUE
 
      ZL = Z(IL) + (XL - X(IL))*(Z(IL+1) - Z(IL))/(X(IL+1) - X(IL))
      SNL = SN(IL) + (XL - X(IL))*(SN(IL+1) - SN(IL))/(X(IL+1) - X(IL))
      ZR = Z(IR-1) + (XR - X(IR-1))*(Z(IR) - Z(IR-1))/(X(IR) - X(IR-1))
      SNR = SN(IR-1) + 
     A          (XR - X(IR-1))*(SN(IR) - SN(IR-1))/(X(IR) - X(IR-1))
 
C     Transfer the subset.  We add vertical walls at each limit. 
C     If these walls are to have resistance, they must be part of
C     a subsection that includes non-sero area.  Thus if we 
C     wish to make the walls have zero resistance, we add two
C     additional subsections and make sure that each vertical
C     end line segment has a unique subsection.  This subsection
C     will then have zero area, the signal to subsequent processing
C     that the line segement contributes no resistance to the flow. 

C     January 29, 2002:  Added walls are frictionless.  Code will need
C     significant revision to allow this to be user selected!

C     Add the two points required to define the vertical wall at 
C     the left-hand extremity of the subset. 
      XS(1) = XL
      ZS(1) = ZMAX
      SNS(1) = SNL
      SBS(1) = SB(IL)
      LSNS(1) = LSN(IL)

      XS(2) = XL
      ZS(2) = ZL
      SNS(2) = SNL
      SBS(2) = SB(IL)
      LSNS(2) = LSN(IL)
 
      NPNTS = 2

C     IL is the index of the point on the boundary that is at 
C     or to the left of XL.  We have already added the point at 
C     XL to the subset.  Thus we want the point just beyond IL 
C     to be the next point on the subset boundary.
      ILP1 = IL + 1

C     IR is the index of the point on the boundary that is at or
C     to the right of XR.  Thus in order to add XR and possibly change
C     its subsection number, we need to stop the transfer of points 
C     just before IR.
      IRM1 = IR - 1
      DO 310 I=ILP1,IRM1
        NPNTS = NPNTS + 1
        XS(NPNTS) = X(I)
        ZS(NPNTS) = Z(I)
        SNS(NPNTS) = SN(I)
        SBS(NPNTS) = SB(I)
        LSNS(NPNTS) = LSN(I)
 310  CONTINUE
 
C     Now add the vertical wall at XR. 
      NPNTS = NPNTS + 1
      XS(NPNTS) = XR
      ZS(NPNTS) = ZR
      SNS(NPNTS) = SNR
      SBS(NPNTS) = SB(IRM1)
      LSNS(NPNTS) = LSNS(NPNTS-1)

      NPNTS = NPNTS + 1
      XS(NPNTS) = XR
      ZS(NPNTS) = ZMAX
      SNS(NPNTS) = SNR
      SBS(NPNTS) = 0
      LSNS(NPNTS) = 0.0


C      WRITE(STDOUT,93) XL,XR
C93    FORMAT('SUBSET results:',/,' XL=',F10.2,' XR=',F10.2)
C      DO 500 I=1,NPNTS
C        WRITE(STDOUT,92) XS(I),ZS(I),SBS(I)
C92      FORMAT(2F10.2,I5)
C500   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   CHKARG
     I                   (STDOUT, MXNDEP, NRZERO, DZLIM, added,
     M                    NDEP, Z)
 
C     + + + PURPOSE + + +
C     Checks and modifies argument sequence in Z(*) to ensure that
C     a near zero point is present and that the argument spacing is
C     not larger than DZLIM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER added, MXNDEP, NDEP, STDOUT
      REAL DZLIM, NRZERO, Z(MXNDEP)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NDEP   - Number of depth values
C     MXNDEP - Maximum number of depth values allowed
C     NRZERO - Value of depth near zero
C     DZLIM  - Maximum value of elevation difference permitted between
C              adjacent entries in a cross section function table
C     Z      - Distinct values of elevation at points on a cross
C               section boundary
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, M, N
      REAL DELZ, DZ, TP
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL RDUP, SORT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT('0*ERR:525*SPACE FOR ELEVATION ARGUMENTS EXHAUSTED')
C***********************************************************************
      N = NDEP

C     Eliminate a too small initial depth.  Can be from roundoff or 
c     may cause computational problems. 
      if(z(2) - z(1).le.nrzero/16.0) then
        do i=3,n
          z(i-1) = z(i)
        end do
        n = n -1
      endif


C     CHECK FOR NEAR ZERO POINT
 
      if(added.eq.0) then
c       only do for sections without an added bottom slot
        IF(Z(2).GT.Z(1) + 1.1*NRZERO) THEN
C         ADD A NEAR ZERO POINT
 
          N = N + 1
          IF(N.LE.MXNDEP) THEN
            Z(N) = Z(1) + NRZERO
          ELSE
            WRITE(STDOUT,50)
            N = N - 1
          ENDIF
        ENDIF
      endif
 
C     CHECK REMAINING ARGUMENT SPACING
 
      DO 500 I=2,NDEP
        DZ = Z(I) - Z(I-1)
        IF(DZ.GT.DZLIM) THEN
C         ADD SOME INTERMEDIATE POINTS
          M = DZ/DZLIM + 1.0
          DELZ = DZ/M
          TP = Z(I-1)
          DO 100 J=1,M-1
            N = N + 1
            TP = TP + DELZ
            IF(N.LE.MXNDEP) THEN
              Z(N) = TP
            ELSE
              WRITE(STDOUT,50)
              N = N - 1
            ENDIF
 100      CONTINUE
        ENDIF
 500  CONTINUE
 
C     N GIVES THE CURRENT NUMBER OF POINTS IN THE LIST AND NO ADDED
C     POINTS DUPLICATE VALUES IN THE LIST.  SORT THE LIST AGAIN TO
C     RETURN IT TO ASCENDING ORDER
 
      CALL SORT
     I         (N,
     M          Z)
 
C     ELIMINATE DUPLICATES.  THE NEAR ZERO POINT SOMETIMES ADDS
C     AN INADVERTANT DUPLICATE.
 
      CALL RDUP
     I         (N,
     M          Z,
     O          NDEP)
 
      RETURN
      END
C
C
C
      SUBROUTINE   XCHK
     I                 (STDOUT, NPNT, ZMIN, X, NFAC,
     M                  Z)
 
C     + + + PURPOSE + + +
C     Check for horizontal segments which are not at ZMIN.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNT, STDOUT
      REAL NFAC, X(NPNT), Z(NPNT), ZMIN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NPNT   - Number of points on boundary of a cross section
C     ZMIN   - Minimum elevation
C     X      - Offsets of points on cross section boundary
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     Z      - Elevation at points on cross section boundary
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL DZ
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + OUTPUT FORMATS + + +
 2    FORMAT(/,' *WRN:504*Line segment ending at (',F10.2,F10.2,
     1 ') is horizontal and not',/,5X,' at minimum elevation.',
     2 '  Right hand end incremented by',F7.3)
C***********************************************************************
C     ESTABLISH STANDARD INCREMENT
 
      IF(NFAC.GT.1.0) THEN
C       English system
        DZ = 0.053
      ELSE
C       Metric sytem
        DZ = 0.0161544
      ENDIF
 
C     SCAN ALL LINE SEGMENTS
 
      DO 100 I=2,NPNT
      IF(ABS(Z(I)-Z(I-1)).GE.0.001) GOTO 100
        IF(ABS(ZMIN-Z(I)).LT.0.001) GOTO 100
          WRITE(STDOUT,2) X(I),Z(I),DZ
          Z(I)=Z(I)+ DZ
 100  CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   REASUB
     I                   (STDOUT, NPNT,
     M                    NSUB, N, SB,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Reassign subsections so that no subsection number is
C     repeated after it has been used one or more times in
C     consecutive sequence.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NPNT, NSUB, STDOUT
      INTEGER SB(*)
      REAL N(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     N      - Manning's n values
C     SB     - Subsection numbers for the line segments
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CURSUB, J, NEWSUB, NXTSUB, OLDSUB
      INTEGER RUN(PMXSUB), SEEN(PMXSUB)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:614* MAXIMUM NUMBER OF SUBSECTIONS=',I5,
     A       ' EXCEEDED IN REASUB.')
C***********************************************************************
      DO 100 J=1,NSUB
        RUN(J) = 0
        SEEN(J) = 0
 100  CONTINUE
 
      OLDSUB = -1
      NXTSUB = NSUB + 1
      IF(NXTSUB.GT.PMXSUB) THEN
        WRITE(STDOUT,50) PMXSUB
        EFLAG = 1
        RETURN
      ENDIF
 
      DO 200 J=1,NPNT-1
        CURSUB = SB(J)
        IF(RUN(CURSUB).EQ.0) THEN
          IF(OLDSUB.NE.-1) THEN
C           RUN OF OLDSUB HAS ENDED
            SEEN(OLDSUB) = SEEN(OLDSUB) + 1
            RUN(OLDSUB) = 0
          ENDIF
 
C         INITIATE NEXT RUN OF CURSUB
          RUN(CURSUB) = 1
          OLDSUB = CURSUB
          IF(SEEN(CURSUB).EQ.0) THEN
            SEEN(CURSUB) = 1
          ELSE
C           ASSIGN NEXT SUBSECTION- START NEW RUN OF CURSUB
            SB(J) = NXTSUB
            N(NXTSUB) = N(CURSUB)
            NEWSUB = NXTSUB
            NXTSUB = NXTSUB + 1
            IF(NXTSUB.GT.PMXSUB) THEN
              WRITE(STDOUT,50) PMXSUB
              EFLAG = 1
              RETURN
            ENDIF
          ENDIF
        ELSE
C         RUN OF CURSUB IS IN PROGRESS
          IF(SEEN(CURSUB).GT.1) THEN
            SB(J) = NEWSUB
          ENDIF
        ENDIF
 200  CONTINUE
 
      NSUB = NXTSUB - 1
      RETURN
      END
C
C
C
      SUBROUTINE   INSPT
     I                  (STDOUT, NPI, XARG, VARTYP,
     M                   NPNT, X, Z, SB, LSN,
     O                   EFLAG)
 
C     + + + PURPOSE + + +
C     Add one or more offsets to the cross section description,
C     interpolating for the elevation, adjust the line segment n,
C     and adjust the subsection number.  The sinuousity is not
C     yet defined.  The line segment n may not be defined in
C     some cases but is adjusted anyway.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NPI, NPNT, STDOUT, VARTYP
      INTEGER SB(PMXPNT)
      REAL LSN(PMXPNT), X(PMXPNT), XARG(NPI), Z(PMXPNT)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NPI    - Number of points to insert
C     XARG   - Vector of offsets of points to be inserted
C     VARTYP - If VARTYP=1 then piecewise linear variation of sinuousity
C               with offset in a cross section is assumed; else
C               if VARTYP=2 then piecewise constant variation of
C               sinuousity with offset in a cross section is assumed.
C     NPNT   - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     LSN    - Line segment Manning's n value
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IS, J, JS, K, SBVAL
      REAL LSNVAL, P, ZVAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:652* Cross section space full in INSPT. Need',
     A       I5,' points but',/,10X,'only',I5,' available.')
 52   FORMAT(/,' *ERR:653* PWC sinuosity offset=',F10.2,' not in the',
     A      ' offsets',/,10X,' for this cross section.')
C***********************************************************************
      IF(NPNT+NPI-2.GT.PMXPNT) THEN
        WRITE(STDOUT,50) PMXPNT, PMXPNT+NPI
        EFLAG = 1
      ELSE
C       WE HAVE ROOM FOR THE ADDED OFFSETS.  THE OFFSETS ARE GIVEN
C       IN ASCENDING ORDER.  THEREFORE SHORTEN THE SEARCH AS WE
C       FIND THE OFFSETS.
 
        IS = 1
 
C       THE FIRST AND LAST OFFSET HAVE BEEN ADDED IN ADJSIN.  DO
C       NOT TRY TO ADD THESE
        DO 500 K=2,NPI-1
C          WRITE(STDOUT,*) ' K=',K,' XARG(K)=',XARG(K)
 
          DO 200 I=IS+1,NPNT
            IF(XARG(K).GE.X(I-1).AND.XARG(K).LE.X(I)) THEN
C             FOUND INTERVAL CONTAINING THE OFFSET. REMEMBER WHERE
C             WE ARE.  START NEXT SEARCH WHERE WE ARE NOW.
 
              IS = I
              JS = I
 
C              WRITE(STDOUT,*) ' FOUND CONTAINING INTERVAL. I=',I
 
C             DO WE NEED TO ADD THE POINT.
              IF(XARG(K).NE.X(I-1).AND.XARG(K).NE.X(I)) THEN
C               ADD THE POINT.  FIND THE NEW VALUES.
                IF(VARTYP.EQ.2) THEN
C                 PIECEWISE CONSTANT SINUOUSITY SHOULD NOT BE HERE.
C                 FOR PWC THE OFFSETS IN THE SINUOUSITY TABLE SHOULD
C                 BE ON A SUBSECTION BOUNDARY.  ALL SUBSECTION
C                 BOUNDARY POINTS SHOULD BE IN X(*).
                  WRITE(STDOUT,52) XARG(K)
 
                  EFLAG = 1
                  RETURN
                ENDIF
 
                P = (XARG(K) - X(JS-1))/(X(JS) - X(JS-1))
                ZVAL = Z(JS-1) + P*(Z(JS) - Z(JS-1))
C                WRITE(STDOUT,*) ' ADDING POINT AT P=',P,' ZVAL=',ZVAL
 
C               THE SUBSECTION IS THE SAME AS THAT FOR THE
C               LEFT-MOST POINT ON THE LINE.  THIS IS TRUE FOR
C               THE LINE SEGMENT N ALSO.
 
                SBVAL = SB(JS-1)
                LSNVAL = LSN(JS-1)
 
C               MOVE VALUES OUT OF THE WAY TO MAKE ROOM FOR THE INSERTION
                DO 100 J=NPNT,JS,-1
                  X(J+1) = X(J)
                  Z(J+1) = Z(J)
                  SB(J+1) = SB(J)
                  LSN(J+1) = LSN(J)
 100            CONTINUE
                NPNT = NPNT + 1
 
C               NOW STORE THE VALUES
 
                X(JS) = XARG(K)
                Z(JS) = ZVAL
                SB(JS) = SBVAL
                LSN(JS) = LSNVAL
 
C               JUMP OUT OF THE LOOP AND GO TO THE NEXT OFFSET
C               TO ADD
                GOTO 500
              ELSE
C               JUMP OUT OF THE LOOP AND GO TO THE NEXT OFFSET
C               TO ADD
                GOTO 500
              ENDIF
            ENDIF
 200      CONTINUE
 
 500    CONTINUE

c        do i=1,npnt
c          write(stdout,54) i, x(i), z(i), sb(i), lsn(i)
c54    format(i5,f12.3,f12.3,i5,f10.3)
c        enddo

      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   TRAN()
 
C     + + + PURPOSE + + +
C     Transfer contents of XSCOMU to XSCOMD
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomd.cmn'
      INCLUDE 'xscomu.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
      INTEGER DVEC(XSCOML), UVEC(XSCOML)
 
C     + + + EQUIVALENCES + + +
      EQUIVALENCE (NPNTU,UVEC(1)), (NPNTD,DVEC(1))
C***********************************************************************
      DO 100 J=1,XSCOML
        DVEC(J)=UVEC(J)
 100    CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   SETOPT
     I                   (STDOUT, CIN,
     O                    SAVOPT, OUTOPT, MONTON, BETOPT)
 
C     + + + PURPOSE + + +
C     Set the options for the output and saving of the cross section
C     tables.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      CHARACTER BETOPT*8, CIN*(*), MONTON*8, OUTOPT*8, SAVOPT*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CIN    - Character string containing options for cross section
C              function table processing
C     SAVOPT - Function table saving option
C     OUTOPT - Output option for the table file for cross section function
C               tables
C     MONTON - Monotonicity flag value
C     BETOPT - Option for computing flux coefficients and critical flow
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'nrdzcm.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER NOPT
      PARAMETER(NOPT=26)
 
C     + + + SAVED VALUES + + +
      INTEGER TABVAL(NOPT)
      CHARACTER TAB(NOPT)*8
      SAVE TAB, TABVAL
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IEND, IS, IT
      CHARACTER NAME*8
 
C     + + + INTRINSICS + + +
      INTRINSIC LEN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BINSER, GETNXT
 
C     + + + DATA INITIALIZATIONS + + +
      DATA
     1TAB   /'EXTEND  ','MONOTONE','NEWBETA ','NEWBETAE','NEWBETAM',
     *       'NOEXTEND','NOOUT   ','NOSAVE  ','OLDBETA ','OUT1    ',
     B       'OUT12   ','OUT20   ','OUT21   ','OUT22   ','OUT23   ',
     C       'OUT24   ','OUT25   ','SAVE    ','SAVE1   ','SAVE12  ',
     D       'SAVE20  ','SAVE21  ','SAVE22  ','SAVE23  ','SAVE24  ',
     E       'SAVE25  '/
      DATA
     1TABVAL/         6,         5,         4,         4,         4,
     A                7,         3,         1,         4,         3,
     B                3,         3,         3,         3,         3,
     C                3,         3,         2,         1,         1,
     D                1,         1,         1,         1,         1,
     E                1/
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:591* ',A8,' is a bad cross section option.')
 52   FORMAT(/,' *BUG:XXX* Invalid index=',I5,' for name=',A8,' in',
     A       ' subroutine SETOPT.')
C***********************************************************************
      IEND = LEN(CIN)
      IS = 1
 
C     SET DEFAULTS
 
      SAVOPT = 'NOSAVE'
      OUTOPT = 'OUT1'
      MONTON = ' '
      BETOPT = 'OLDBETA'
 
 100  CONTINUE
        CALL GETNXT
     I             (CIN, IS,
     O              IT, NAME)
        IS =  IT
        IF(NAME.NE.' ') THEN
C         LOOKUP THE OPTION NAME AND USE THE TABULATED VALUE
C         TO TAKE AN ACTION.
 
          CALL BINSER
     I               (NAME, NOPT, TAB,
     O                I)
          IF(I.EQ.0) THEN
            WRITE(STDOUT,50) NAME
          ELSE
            I = TABVAL(I)
            GOTO(1, 2, 3, 4, 5, 6, 7),I
              WRITE(STDOUT,52) I, NAME
              STOP 'Abnormal stop. Errors found.'
 
 1          CONTINUE
C             MOST OF THE SAVE OPTIONS.
              SAVOPT = NAME
              GOTO 49
 
 2          CONTINUE
C             SAVE BY ITSELF WITHOUT A SUFFIX NUMBER
              SAVOPT = 'SAVE1'
              GOTO 49
 
 3          CONTINUE
C             OUTPUT OPTIONS
              OUTOPT = NAME
              GOTO 49
 
 4          CONTINUE
C             VELOCITY COEFFICIENT FLAGS
              BETOPT = NAME
              GOTO 49
 
 5          CONTINUE
C             MONOTONICITY FLAG
              MONTON = NAME
              GOTO 49
 
 6          CONTINUE
C             SET EXTEND OPTION
              EXTEND = 1
              GOTO 49
 
 7          CONTINUE
C             CLEAR THE EXTEND OPTION
              EXTEND = 0
              GOTO 49
 
 49         CONTINUE
          ENDIF
        ENDIF
        IF(IS.GE.IEND) GOTO 200
        GOTO 100
 200  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   INFQXE
     I                   (STDIN, STDOUT,
     M                    TABDIR, EFLAG,
     O                    TAB, STAT, NPNT, NSUB, X, Z, SB, N, LEFT,
     O                    RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX, LSN,
     O                    NVAR, NATY, YATN, NNY,
     O                    zone, hgrid, vdatum, unitsys, basis)
 
C     + + + PURPOSE + + +
C     Input cross section in FEQ extended format
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NPNT, NSUB, STDIN, STDOUT, TAB
      INTEGER NNY(PMXSUB), NVAR(PMXSUB), SB(PMXPNT), TABDIR(*)
      REAL LEFT, LSN(PMXPNT), N(PMXSUB), NATY(9,PMXSUB), RIGHT, STAT,
     A     X(PMXPNT), YATN(9,PMXSUB), Z(PMXPNT), ZMAX
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8, zone*8, hgrid*8,
     b          vdatum*8, unitsys*8, basis*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TAB    - Table number
C     STAT   - Station value
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     N      - Manning's n values
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
C     SAVOPT - Function table saving option
C     OUTOPT - Output option for the table file for cross section function
C               tables
C     BETOPT - Option for computing flux coefficients and critical flow
C     ZMAX   - Maximum elevation
C     LSN    - Line segment Manning's n value
C     NVAR   - Flag for variation of Manning's n in each subsection
C     NATY   - Mannings's n value at depth in YATN
C     YATN   - Depth values for the Manning's n values in NATY
C     NNY    - Number of values for Manning's n variation with depth
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xtadd.cmn'
      INCLUDE 'ftable.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER NVAL
      PARAMETER(NVAL=20)
 
C     + + + SAVED VALUES + + +
      INTEGER BASNUM, VTYPE(NVAL)
      SAVE BASNUM, VTYPE
 
C     + + + LOCAL VARIABLES + + +
      INTEGER EFLAG2, I, ISUB, IVARN, J, MVAL, NN, NTEMP,
     A        OLDSUB, OPT, NAVM, L
      INTEGER CLEN(NVAL), IVAL(NVAL), TERML(NVAL), TERMCLS(NVAL)
      REAL RVAL(NVAL), HSHIFT, SCALE, SHIFT, XOLD, VSCALE
      DOUBLE PRECISION DPVAL(NVAL)
      CHARACTER  CVAL(NVAL)*256, HEAD*80, LINE*80,
     A          MONTON*8, TERM(NVAL)*1, VARN*4, EXTEND_STRING*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN, MOD
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      CHARACTER GET_TABID*16
      EXTERNAL GETVAL, inline, REASUB, SETOPT, TABCHK, GET_TABID, LENSTR
 
C     + + + DATA INITIALIZATIONS + + +
      DATA VTYPE/2*2,1,17*2/, BASNUM/2/
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,A)
 2    FORMAT(8X,F10.0,6X,F10.0,7X,F10.0)
 3    FORMAT(A4,I5,6F10.0)
 4    FORMAT(A80)
C 6    FORMAT(5X,A4,7X,F10.0,7X,F10.0)
 6    FORMAT(5X,A4,7X,F10.0,7X,F10.0,8X,F10.0,8X,F10.0)
 31   FORMAT(9X,6F10.0)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(/,' TABID=',A,5(1X,A8))
 52   FORMAT(' VARN=',A4,'  Hor. scale factor=',F10.3,
     A ' Vert. shift=',F10.3,/,'  Vert. scale factor=',F10.3,
     B ' Hor. shift=',F10.3)
 53   FORMAT(' STATION=',F10.3,' LEFT=',F10.1,' RIGHT=',F10.1)
 54   FORMAT(' STATION=',F10.3)
 55   FORMAT(1X,A4,I5,6F6.3,1X,/,(10X,6F6.3))
 56   FORMAT(1X,A80)
 57   FORMAT(1X,F10.3,F10.2,I5,F6.3)
 58   FORMAT(' *ERR:504* NUMBER OF SUBSECTIONS=',I5,' > ',I5)
 59   FORMAT(' *ERR:505* SUBSECTION NUMBER TOO LARGE AT OFFSET=',
     1       F10.1)
 60   FORMAT(' *WRN:501* SUBSECTION VALUE FOR FIRST POINT IS',
     1 ' MISSING.  SUBSECTION=1 ASSUMED.')
 61   FORMAT(' *ERR:506* ONLY ONE POINT GIVEN ON BOUNDARY OF THE',
     1  ' CROSS SECTION.')
 62   FORMAT(' *ERR:507* NUMBER OF POINTS IN CROSS SECTION >',I5)
 63   format(' Invert shift=',f10.3,' Water-surface tabid=',a) 
 65   format(' ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a        ' BASIS=',a8)
 66   FORMAT(' Selection of beta option "NEWBETA" implies checking for',
     A      ' monotonicity.')
 67   FORMAT(/,' *WRN:554* Extending left end of cross section',
     A   ' by ',F8.3)
 68   FORMAT(/,' *WRN:555* Extending right end of cross section',
     A   ' by ',F8.3)
 69   FORMAT(/,' *WRN:556* Some point in cross section higher than',
     A  ' either end.',/,10X,'All area above minimum end elevation',
     b  ' is ignored.')
 70   FORMAT(/,' *ERR:647* Manning''s N = 0 in subsection',I5,
     A /,11X,'Value set to 1.0.  Please correct and recompute.')
 72   FORMAT(/,' *ERR:556* Found only',I3,' values on line but ',
     A       'expected',I3)
 74   FORMAT(/,' *ERR:557* Vertical variation of rougness coef.',
     A       ' already defined',/,10X,' in subsection number',I4)
 76   FORMAT(/,' *ERR:564* Missing n-value at depth number ',I3)
 78   FORMAT(/,' *ERR:567*', A4,' is invalid option for vertical ',
     A       'variation of roughness.')
 80   FORMAT(1X,F10.3,F10.2,I5,F6.3,5(F6.2,F6.3))
 82   FORMAT(1X,F10.3,F10.2,I5)
 84   FORMAT(1X,'    OFFSET ELEVATION SUBS    N0    Y1    N1    Y2',
     A '    N2    Y3    N3    Y4    N4    Y5    N5')
86    FORMAT(' GISID=',A16,' EASTING=',F15.3,' NORTHING=',F15.3)
87    FORMAT('  Processing FEQXEXT TabId= ',A)
C***********************************************************************
C     CLEAR THE LOCAL ERROR FLAG FOR THE GETVAL SUBROUTINE
      EFLAG2 = 0
 
      CALL GET_XSEC_HEADER(STDIN, STDOUT,
     M   EFLAG,
     O   TAB, SAVOPT, OUTOPT, MONTON, BETOPT, EXTEND, 
     O   STAT, LEFT, RIGHT, VARN, NAVM, SCALE, SHIFT,
     O   VSCALE, HSHIFT, NSUB, N, zone, hgrid, vdatum, unitsys, basis)
      IF(EXTEND.EQ.0) THEN
        EXTEND_STRING = 'NOEXTEND'
      ELSE
        EXTEND_STRING = 'EXTEND'
      ENDIF
      TABID = GET_TABID(ABS(TAB))
      L = LENSTR(TABID)
      WRITE(*,87) TABID(1:L)
      WRITE(STDOUT,51) TABID(1:L), SAVOPT, OUTOPT, BETOPT, 
     A                 EXTEND_STRING, MONTON
      if(zone /= 'NONE') then
        write(stdout,65) zone, hgrid, vdatum, unitsys, basis
      endif
 
      IF(TAB.LT.0) THEN
        SLOT = 0.03
        TAB = -TAB
        NOCM = 1
      ELSE
         NOCM = 0
         SLOT = 0.0
      ENDIF
      IF(TAB.GE.0) CALL TABCHK
     I                        (STDOUT, PMXTAB,
     M                         TAB, TABDIR, EFLAG)
 
      LEFT = SCALE*LEFT
      RIGHT = SCALE*RIGHT
      IF(VARN.EQ.' ') THEN
        VARN = 'HYDY'
      ENDIF
      IF(VARN.EQ.'HYDY') THEN
        IVARN = 1
      ELSEIF(VARN.EQ.'MAXY') THEN
        IVARN = 2
      ELSEIF(VARN.EQ.'NCON') THEN
        IVARN = 0
      ELSE
        WRITE(STDOUT,78) VARN
        EFLAG = 1
        IVARN = 1
      ENDIF
      IF(GISID.NE.' ') THEN
        WRITE(STDOUT,86) GISID, EASTING, NORTHING
      ENDIF
      IF(LEFT.GE.RIGHT) THEN
        WRITE(STDOUT,54) STAT
      ELSE
        WRITE(STDOUT,53) STAT, LEFT, RIGHT
      ENDIF
      WRITE(STDOUT,52) VARN, SCALE, SHIFT, VSCALE, HSHIFT
      if(dinvert.ne.0.0) then
        write(stdout,63) dinvert, ws_tabid
      endif
      WRITE(STDOUT,55) 'NSUB', NSUB, (N(J),J=1,NSUB)
 
      DO 131 J=1,NSUB
        IF(N(J).LE.0.0) THEN
          WRITE(STDOUT,70) J
          N(J) = 1.0
          EFLAG = 1
        ENDIF
C       SET THE VARIATION OF MANNING'S N FLAG TO INITIAL VALUE
C       Initial value assumes that n is constant with changes in
C       depth.
        NVAR(J) = 0
 131  CONTINUE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) HEAD
C      WRITE(STDOUT,56) HEAD
      WRITE(STDOUT,84)
 
 
C     GET THE CO-ORDINATES ON THE BOUNDARY OF THE CROSS SECTION
 
C     SET THE VALUE FOR CHECKING FOR MONOTONE VARIATION OF OFFSET
      XOLD = -1.E30
 
      NPNT=0
      ZMAX = -1.E30
 
C     SET THE SUB SECTION NUMBER DEFAULT
      DO 100 I=1,PMXPNT
        SB(I) = 0
 100  CONTINUE
 
C     THE OLD SUBSECTION FLAG TO USE DURING OUTPUT
      OLDSUB = 0
 
C     START LOOP OVER INPUT LINES DEFINING THE BOUNDARY OF THE
C     CROSS SECTION
 
 135  CONTINUE
        NPNT=NPNT+1
        IF(NPNT.GT.PMXPNT) THEN
          WRITE(STDOUT,62) PMXPNT
          EFLAG = 1
          NPNT = PMXPNT
          RETURN
        ENDIF
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
 
C       SET THE DEFAULT VALUES FOR THE CURRENT LINE
        DO 137 I=1,NVAL
          RVAL(I) = -1.E30
 137    CONTINUE
        IVAL(3) = -2147483647
 
        OPT = 0
        CALL GETVAL
     I             (STDOUT, LINE, NVAL, OPT,
     O              VTYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O              EFLAG2, TERM, TERML, TERMCLS,
     O              MVAL)
 
C       MVAL GIVES THE NUMBER OF VALUES GIVEN BY THE USER.  THIS
C       COUNT INCLUDES THE ASTERISK FOR A DEFAULT VALUE IN THE FIELD
C       AS WELL AS COMMAS USED TO DEFAULT A FIELD.  THE COUNT GIVES
C       THE NUMBER OF ITEMS, ACTUAL OR DEFAULT UP UNTIL THE END OF LINE
C       CHARACTER OR UNTIL THE END OF LINE.  MVAL GIVES MUCH OF THE
C       INFORMATION ABOUT WHAT THE USER HAS STATED AND INTENDED.
 
C       THERE MUST ALWAYS BE AT LEAST BASNUM VALUES INPUT.  THE EXCESS
C       OVER BASNUM IS THEN USED TO DECIDE WHAT TO DO WITH THE NUMBERS.
C       GETVAL WILL COMPLAIN ABOUT CONVERSION ERRORS.  EFLAG WILL BE
C       SET IF ANY ERROR CONDITION HAS BEEN FOUND.
 
C       FIRST, CHECK FOR ERRORS IN GETVAL. IF FOUND GO TO THE NEXT LINE
C       OF USER INPUT.
        IF(EFLAG2.NE.0) THEN
C         CLEAR THE LOCAL FLAG
          EFLAG2 = 0
          EFLAG = 1
          NPNT = NPNT - 1
          GOTO 135
        ENDIF
 
 
C       NOW CHECK TO MAKE SURE THAT THERE IS THE PROPER MINIMUM NUMBER OF
C       VALUES.
 
        IF(MVAL.LT.BASNUM) THEN
          EFLAG = 1
          WRITE(STDOUT, 72) MVAL, BASNUM
          EFLAG2 = 0
          NPNT = NPNT - 1
          RETURN
        ENDIF
 
C       NO ERRORS DETECTED AT THIS POINT. NOW PROCESS THE VALUES
C       THAT MUST ALWAYS APPEAR
 
        IF(BASNUM.EQ.2) THEN
C         WE HAVE ONLY OFFSET AND ELEVATION
 
          X(NPNT) = SCALE*RVAL(1) + HSHIFT
          Z(NPNT) = VSCALE*RVAL(2) + SHIFT
          
        ELSE
C         WE HAVE PLANE COORDINATES PLUS ELEVATION
          WRITE(STDOUT,*) ' PLANE CO-ORDINATE OPTION MISSING'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
 
 
C       FIND MAXIMUM ELEVATION IN CROSS SECTION FOR LATER CHECKING
        ZMAX = MAX(Z(NPNT), ZMAX)
 
C       CHECK FOR MONOTONICITY OF TOP WIDTH.  THIS REQUIRES THAT
C       THE OFFSET NEVER DECREASE.
 
        IF(MONTON.EQ.'MONOTONE') THEN
          IF(X(NPNT).LT.XOLD) THEN
            WRITE(STDOUT,*)
     A  ' *ERR:508* SECTION VIOLATES MONOTONICITY AT OFFSET=', X(NPNT)
            EFLAG = 1
          ENDIF
        ENDIF
        XOLD = X(NPNT)
 
C       NOW, CHECK THE VARIOUS VALUES AND TAKE THE ACTIONS AS REQUIRED.
C       THE KEY VALUES ARE: SUBSECTION NUMBER, N0 VALUE, AND THE FIRST
C       DEPTH(THAT IS THE FIRST NON-ZERO DEPTH) FOR VERTICAL VARIATION
C       OF MANNING'S N IN A SUBSECTION.
 
        IF(IVAL(BASNUM+1).LT.-1.OR.IVAL(BASNUM+1).EQ.0) THEN
C         USER OMITTED THE SUBSECTION NUMBER.  GET THE PREVIOUS
C         NUMBER AND USE IT IF IT EXISTS.
 
          IF(NPNT.EQ.1) THEN
C           DOES NOT EXIST.  WRITE WARNING AND FORCE SUBSECTION 1.
            WRITE(STDOUT,60)
            SB(1) = 1
          ELSE
            SB(NPNT) = SB(NPNT-1)
          ENDIF
        ELSEIF(IVAL(BASNUM+1).GT.0) THEN
C         STORE THE VALUE.
 
          IF(IVAL(BASNUM+1).GT.NSUB) THEN
            WRITE(STDOUT,59) X(NPNT)
            EFLAG = 1
            IVAL(BASNUM+1) = NSUB
          ENDIF
          SB(NPNT) = IVAL(BASNUM+1)
 
        ELSEIF(IVAL(BASNUM+1).EQ.-1) THEN
C         BOUNDARY DEFINITION TERMINATED.  NO OTHER INFORMATION
C         MAY APPEAR ON THE LINE OF THE FINAL POINT.
 
          SB(NPNT) = -1
          GOTO 145
        ENDIF
 
        ISUB = SB(NPNT)
 
C       NOW CHECK VALUES FOR THE VALUE OF N0.  IF OMITTED BY THE
C       USER IT TAKES ON THE VALUE OF THE MANNING'S N GIVEN FOR THE
C       SUBSECTION.  ANY NON-NEGATIVE VALUE IS VALID.
 
        IF(RVAL(BASNUM+2).LT.0.0) THEN
C         TAKE USER NEGATIVE VALUES TO MEAN THE SAME THING AS
C         OMITTING THE VALUE.  ASSIGN THE LINE SEGMENT N VALUE
C         FROM THE SUBSECTION.
 
          LSN(NPNT) = N(ISUB)
        ELSE
          LSN(NPNT) = RVAL(BASNUM+2)
        ENDIF
 
C       CHECK THE FIRST DEPTH FOR VERTICAL VARIATION OF MANNING'S N.
C       THESE VALUES CAN APPEAR ONLY ONCE PER SUBSECTION.
 
        IF(RVAL(BASNUM+3).GT.0.0) THEN
C         A NON-NEGATIVE VALUE HAS BEEN GIVEN.  USER IS GIVING DATA
C         ON VERTICAL VARIATION OF ROUGHNESS IN A SUBSECTION.
C         HAS THIS BEEN DONE ALREADY FOR THIS SUBSECTION?
          IF(NVAR(ISUB).GT.0) THEN
C           YES. OUTPUT ERROR AND IGNORE THE NEW DATA.
            WRITE(STDOUT,74) ISUB
            EFLAG = 1
          ELSE
C           GET THE NUMBER OF VALUES GIVEN FOR VARIATION OF ROUGHNESS IN
C           THE VERTICAL.  COUNTING FROM THE FIRST DEPTH THE NUMBER OF
C           VALUES MUST BE EVEN AND NOT GREATER THAN 8.
 
 
            NN = MVAL - (BASNUM + 2)
            IF(MOD(NN,2).NE.0) THEN
C             VALUE OF ROUGHNESS IS MISSING.
              WRITE(STDOUT,76) NN/2
              EFLAG = 1
              NN = NN - 1
            ENDIF
 
C           STORE THE VALUES IN THE SUBSECTION NUMBER INDEXED ARRAYS
            IF(NN.GT.0) THEN
              NNY(ISUB) = NN/2 + 1
              YATN(1,ISUB) = 0.0
              NATY(1,ISUB) = LSN(NPNT)
              DO 138 I=1,NN/2
                YATN(I+1,ISUB) = RVAL(BASNUM+3+2*(I-1))
                NATY(I+1,ISUB) = RVAL(BASNUM+4+2*(I-1))
 138          CONTINUE
 
C             SET THE NATURE OF THE VARIATION IN THE VERTICAL
 
              NVAR(ISUB) = IVARN
            ENDIF
          ENDIF
        ENDIF
 
 
 
C       OUTPUT THE VALUES FOR THE CURRENT LINE
 
        IF(NVAR(ISUB).EQ.0) THEN
          WRITE(STDOUT,57) X(NPNT),Z(NPNT), ISUB, LSN(NPNT)
        ELSE
C         VERTICAL VARIATION OF N.  OUTPUT DETAILS ONLY FOR FIRST
C         POINT IN EACH SUBSECTION HAVING SUCH VARIATION
          IF(ISUB.NE.OLDSUB) THEN
C           FIRST LINE OF A NEW SUB-SECTION
            OLDSUB = ISUB
            WRITE(STDOUT,80) X(NPNT), Z(NPNT), ISUB, LSN(NPNT),
     A        (YATN(J,ISUB), NATY(J,ISUB),J=2,NNY(ISUB))
          ELSE
C           NOT FIRST LINE OF A NEW SUB-SECTION
            WRITE(STDOUT,82) X(NPNT), Z(NPNT), ISUB
          ENDIF
        ENDIF
 
 
        GOTO 135
 
 145  CONTINUE
C     OUTPUT THE FINAL LINE
      WRITE(STDOUT,82) X(NPNT), Z(NPNT), SB(NPNT)
 
 
      IF(NPNT.GT.1) GOTO 170
        WRITE(STDOUT,61)
        EFLAG = 1
 170  CONTINUE
 
C     CHECK THE CROSS SECTION FOR NONSENSE BEHAVIOR AT THE END
 
      IF(Z(1).LE.Z(2).AND.X(2).GT.X(1)) THEN
C       THE LEFT MOST LINE SEGMENT HAS UPWARD SLOPE, THEREFORE HIGH
C       POINT IS NOT AT THE LEFT LIMIT.
 
        WRITE(STDOUT,*) ' *WRN:502* UNEXPECTED SLOPE',
     A   ' AT LEFT END.'
        WRITE(STDOUT,*) '  SLOPE EXPECTED TO BE < 0  AT LEFT BOUNDARY'
      ENDIF
      IF(Z(NPNT).LE.Z(NPNT-1).AND.X(NPNT).GT.X(NPNT-1)) THEN
C       THE RIGHT MOST LINE SEGMENT HAS DOWNWARD SLOPE
 
        WRITE(STDOUT,*) ' *WRN:503* UNEXPECTED SLOPE',
     A   ' AT RIGHT END.'
        WRITE(STDOUT,*) '  SLOPE EXPECTED TO BE > 0 AT RIGHT BOUNDARY'
      ENDIF
 
      IF(EXTEND.EQ.1) THEN
C       CHECK FOR ONE END BEING HIGHER THAN THE OTHER
 
        IF(ABS(ZMAX - Z(1)).GT.EPSDIF) THEN
          WRITE(STDOUT,67) ZMAX - Z(1)
        ENDIF
        IF(ABS(ZMAX - Z(NPNT)).GT.EPSDIF) THEN
          WRITE(STDOUT,68) ZMAX - Z(NPNT)
        ENDIF
      ELSE
C       CHECK FOR AN INTERMEDIATE POINT BEING HIGHER THAN EITHER
C       END
        IF(ABS(ZMAX - Z(1)).GT.EPSDIF.AND.
     A     ABS(ZMAX - Z(NPNT)).GT.EPSDIF) THEN
          WRITE(STDOUT,69)
        ENDIF
 
        ZMAX = MIN(Z(1), Z(NPNT))
      ENDIF
 
C     CHECK TO SEE IF SUBSECTION NUMBERS HAVE BEEN REPEATED IN
C     SEPARATE RUNS.
 
      NTEMP = NSUB
      CALL REASUB
     I           (STDOUT, NPNT,
     M            NSUB, N, SB,
     O            EFLAG)
      IF(NSUB.NE.NTEMP) THEN
C       OUTPUT THE NEW ASSIGNMENTS
        WRITE(STDOUT,*) ' '
        WRITE(STDOUT,*) ' *WRN:548* Unwise use of subsection numbers.'
        WRITE(STDOUT,*) '     Subsections have been added to avoid',
     A                  ' repeated usage.'
        WRITE(STDOUT,*) '     Old NSUB=',NTEMP,' New NSUB=',NSUB,
     A                  ' Please check for validity.'
        WRITE(STDOUT,55) 'NSUB', NSUB, (N(J),J=1,NSUB)
        WRITE(STDOUT,56) HEAD
        DO 200 J=1,NPNT
          WRITE(STDOUT,57) X(J),Z(J),SB(J)
 200    CONTINUE
      ENDIF
 

c     Optionally adjust part of the boundary in elevation.  
      if(dinvert.ne.0.0) then
c       Make an adjustment.  Requires a table of type 2 that defines
c       the elevation below which we make an adjustment. 

        CALL GET_INTERNAL_TAB_NUMBER
     I                             (STDOUT, ws_tabid,
     M                              EFLAG,
     O                              ws_tab)
        CALL CHKTAB
     I             (2, STDOUT, FTPNT, MFTNUM,
     M              ws_tab,
     O              EFLAG)
        if(eflag.eq.0) then
          call adj_invert(stdout, npnt, stat, ws_tab, dinvert,
     m                    z,
     o                    eflag)
        endif
      endif                                

      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   FBASEL
     I                   (SNFLG, ZI, NPNT, NSUB, X, Z, SB, LSN, NVAR,
     I                    SN, NATY, YATN, NNY, NFAC, BETOPT,
     O                    N, TS, PS, AS, YBS, NS, SUMQ, SUMFM, SUMFE,
     O                    SUMDQ, SUMDFM, SUMDFE, SUMMA, SUMMQ, YSMX,
     O                    SBSN, QS, KS)
 
C     + + + PURPOSE + + +
C     Find basic elements for cross section at a given elevation.
C     the only horizontal lines segments which should be present
C     are at the minimum elevation  of the cross section.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNT, NSUB, SNFLG
      INTEGER NNY(NSUB), NVAR(NSUB), SB(NPNT)
      REAL KS(NSUB), LSN(NPNT), N(NSUB), NATY(9,NSUB), NFAC,
     A     NS(NSUB), PS(NSUB), QS(NSUB), SBSN(NSUB), SN(NPNT), TS(NSUB),
     B     X(NPNT), YATN(9,NSUB), YSMX(NSUB), Z(NPNT), ZI
      DOUBLE PRECISION AS(NSUB), YBS(NSUB)
      DOUBLE PRECISION SUMDFE, SUMDFM, SUMDQ, SUMFE, SUMFM, SUMMA,
     A                 SUMMQ, SUMQ
      CHARACTER BETOPT*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     SNFLG  - Flag for sinuousity computations
C     ZI     - Current water surface elevation for computing elements
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     LSN    - Line segment Manning's n value
C     NVAR   - Flag for variation of Manning's n in each subsection
C     SN     - Sinuousity at a point on a cross section boundary
C     NATY   - Mannings's n value at depth in YATN
C     YATN   - Depth values for the Manning's n values in NATY
C     NNY    - Number of values for Manning's n variation with depth
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     BETOPT - Option for computing flux coefficients and critical flow
C     N      - Manning's n values
C     TS     - Top width in each subsection
C     PS     - Wetted perimeter for each subsection
C     AS     - Area for each subsection in the cross section
C     YBS    - First moment of area about surface for each subsection
C     NS     - Manning's n for a subsection
C     SUMQ   - Flow(water flux)
C     SUMFM  - Momentum flux
C     SUMFE  - Energy flux
C     SUMDQ  - Derivative w. r. t. depth of the flow(water flux)
C     SUMDFM - Derivative w. r. t. depth of the momentum flux
C     SUMDFE - Derivative w. r. t. depth of the energy flux
C     SUMMA  - Area weighted sinuousity integral
C     SUMMQ  - Flow weighted sinuousity integral
C     YSMX   - Maximum local depth in each subsection
C     SBSN   - Subsection sinuousity
C     QS     - Flow distribution value in a subsection to account for
C               flow path length variations
C     KS     - Conveyance for each subsection
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'gnicom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + SAVED VALUES + + +
      REAL T10D3, T1D3, T2D3, T4D3, T5D3, T7D3, T8D3
      SAVE T10D3, T1D3, T2D3, T4D3, T5D3, T7D3, T8D3
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IGS, IS, J, JMAX, JMIN, NBFLAG, NBX
      REAL  DEPTH, DX, DZ, HALFDX, P, PVEC(PMXPNT), SNI, SNL, SNR, V,
     A     XI, XL, XMID, XR, YL, YR, ZL, ZR
      DOUBLE PRECISION DA
      DOUBLE PRECISION ALOC, BLOC, C, C1, C2, C3, C4, DALOC, DBLOC, DFE,
     A                 DFM, FE, FM, H, K, M, MU, Q, S, SNLOC, SQRTSN, W
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, DBLE, MAX, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FNDYDN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDYDN
 
C     + + + DATA INITIALIZATIONS + + +
      DATA T1D3/0.3333333/, T2D3/0.6666667/, T4D3/1.333333/,
     A  T5D3/1.666667/, T7D3/2.333333/, T8D3/2.666667/, T10D3/3.333333/
C***********************************************************************
C     To accomodate n varying with depth of water, we must make two passes
C     over the cross section boundary.  The first pass computes the
C     pure geometric elements: top width, area, first moment of area, and
C     wetted perimeter.  Once these values are known, the value of n in
C     those subsections having a variable n can be computed.  The second
C     pass computes the values that depend in some way on roughness.  These
C     include the sinuosity elements and the NEWBETA computation of the
C     velocity distribution coefficients.
 
 
c      write(STD6,*) 'In FBASEL ZI=',ZI

C     SET A LOCAL FLAG FOR THE NEW BETA OPTION TO MAKE DECISIONS IN
C     A LOOP FASTER
 
      IF(BETOPT(1:7).EQ.'NEWBETA') THEN
        NBFLAG = 1
        IF(BETOPT(8:8).EQ.'X') THEN
          NBX = 1
        ELSE
          NBX = 0
        ENDIF
      ELSE
        NBFLAG = 0
        NBX = 0
      ENDIF
 
C     CLEAR THE SUMMING VARIABLES
 
      SUMQ = 0.0
      SUMFM = 0.0
      SUMFE = 0.0
      SUMDQ = 0.0
      SUMDFM = 0.0
      SUMDFE = 0.0
      SUMMA = 0.0
      SUMMQ = 0.0
 
      DO 100 J=1,NSUB
        TS(J) = 0.
        AS(J) = 0.D0
        PS(J) = 0.
        YBS(J) = 0.D0
        NS(J) = 0.
        YSMX(J) = 0.
        SBSN(J) = 0.
        QS(J) = 0.
        KS(J) = 0.
 100  CONTINUE
 
C     CLEAR THE INDICES FOR REMEMBERING THE RANGE OF BOUNDARY POINTS FOR
C     THE CURRENT ELEVATION, ZI
 
      JMIN = 0
      JMAX = 0
 
 
      DO 500 J=2,NPNT
        ZL = Z(J-1)
        ZR = Z(J)
 
c        write(std6,*) 'ZL=',ZL,' ZR=',ZR
C       IS THE LINE SEGMENT COMPLETELY ABOVE THE CURRENT WATER LEVEL, ZI?
        IF(ZI.LT.ZL.AND.ZI.LT.ZR) GOTO 500
C         NO, AT LEAST PART OF THE LINE SEGMENT IS WET.
          XL = X(J-1)
          XR = X(J)
          IS = SB(J-1)
c          write(std6,*) ' At least part under water:'
c          write(std6,*) ' XL=',XL,' XR=',XR
C         REMEMBER THE MINIMUM AND THE MAXIMUM VALUE OF J FOR SEGMENTS THAT
C         ARE UNDER WATER IN WHOLE OR IN PART
          IF(JMIN.EQ.0) JMIN = J
          JMAX = J
          IF(ZI.GE.ZL.AND.ZI.GE.ZR) GOTO 300
 
C           INTERSECTION BETWEEN WATER SURFACE AND CURRENT
C           BOUNDARY LINE SEGMENT.  FIND OFFSET OF POINT OF
C           INTERSECTION.
 
            XI = XL + (ZI - ZL)*(XR - XL)/(ZR - ZL)
            IF(ZL.LE.ZR) THEN
              DX = XI - XL
              XR = XI
              ZR = ZI
              DZ = ZI - ZL
              YR = 0.0
              YL = DZ
              YSMX(IS) = MAX(YSMX(IS), YL)
            ELSE
              DX = XR - XI
              XL = XI
              ZL = ZI
              DZ = ZI - ZR
              YL = 0.0
              YR = DZ
              YSMX(IS) = MAX(YSMX(IS), YR)
            ENDIF
 
c            write(std6,*) ' Intersection with water surface:'
c            write(std6,*) ' DX=',DX,' XI=',XI,' DZ=',DZ
c            write(std6,*) ' YL=',YL,' YR=',YR
            TS(IS) = TS(IS) + DX
            DA = DBLE(DX)*DZ/2.D0
            AS(IS) = AS(IS) + DA
c            write(std6,*) ' DA=',DA
c            write(std6,*) ' '
            YBS(IS) = YBS(IS) + DA*DZ/3.D0
            GOTO 490
 300      CONTINUE
 
C           LINE SEGMENT BELOW WATER SURFACE
 
            DX = XR - XL
            DZ = ZR - ZL
            YL = ZI - ZL
            YR = ZI - ZR
c            write(std6,*) 'Segment below water:'
c            write(std6,*) ' DX=',DX,' DZ=',DZ,' YL=',YL,' YR=',YR
            YSMX(IS) = MAX(YSMX(IS), YL, YR)
            TS(IS) = TS(IS) + DX
            AS(IS) = AS(IS) + DBLE(DX)*(YL + YR)/2.
c            write(std6,*) ' DA=',DBLE(DX)*(YL + YR)/2.
c            write(std6,*) ' '

            YBS(IS) = YBS(IS) + DBLE(DX)*(YL*YR + (YR - YL)**2./3.)/2.
 490      CONTINUE
 
          P = SQRT(DX*DX + DZ*DZ)
          IF(LSN(J-1).GT.0.0) THEN
C           COUNT ONLY PERIMETER THAT HAS POSITIVE N VALUES
            PS(IS) = PS(IS) + P
          ENDIF
          PVEC(J-1) = P
 
C         NOW COMPUTE THE VALUE USED TO ESTIMATE THE COMPOSITE N WITHIN
C         A SUB-SECTION WHEN THE N VALUES VARY WITH LINE SEGMENTS BUT
C         NOT WITH DEPTH.  THEY NEED NOT VARY AND OFTEN WILL NOT BUT
C         PROCESS WITH SAME CODE.
 
          IF(NVAR(IS).EQ.0) THEN
            NS(IS) = NS(IS) + P*LSN(J-1)
          ENDIF
 
 500    CONTINUE
 
C       FIND THE VALUES OF ROUGHNESS FOR ANY SUBSECTIONS HAVING DEPTH-
C       DEPENDENT N.  INCLUDES POSSIBILITY THAT N IS CONSTANT.
 
        DO 600 IS=1,NSUB
          IF(NVAR(IS).EQ.1) THEN
C           N DEPENDS ON HYDRAULIC DEPTH IN THE SUBSECTION
            IF(TS(IS).GT.0.0) THEN
              DEPTH = AS(IS)/TS(IS)
            ELSE
              IF(AS(IS).GT.0.0) THEN
                WRITE(STD6,*) ' *BUG:XXX* AREA OR TOP WIDTH INVALID',
     A                     ' IN SUB. FBASEL'
                WRITE(STD6,*) ' IS=',IS, ' T=',TS(IS),' A=',AS(IS)
                STOP 'Abnormal stop. Errors found.'
              ELSE
                DEPTH = 0.
              ENDIF
            ENDIF
          N(IS) = FNDYDN(DEPTH, NNY(IS), YATN(1,IS), NATY(1,IS))
        ELSEIF(NVAR(IS).EQ.2) THEN
C         N DEPENDS ON MAXIMUM DEPTH IN THE SUBSECTION
          DEPTH = YSMX(IS)
          N(IS) = FNDYDN(DEPTH, NNY(IS), YATN(1,IS), NATY(1,IS))
        ELSE
C         COMPUTE COMPOSITE N AND USE IT IN THE SUBSECTION.
          IF(PS(IS).GT.0.0) THEN
            N(IS) = NS(IS)/PS(IS)
            
          ENDIF
        ENDIF
 600  CONTINUE
 
 
      IF(JMIN.LT.2) JMIN = 2
C     START THE SECOND PASS
      DO 900 J=JMIN,JMAX
        ZL = Z(J-1)
        ZR = Z(J)
        IF(ZI.LT.ZL.AND.ZI.LT.ZR) GOTO 900
          XL = X(J-1)
          XR = X(J)
          SNL = SN(J-1)
          SNR = SN(J)
          IS = SB(J-1)
 
          IF(ZI.GE.ZL.AND.ZI.GE.ZR) GOTO 700
 
C           INTERSECTION BETWEEN WATER SURFACE AND CURRENT
C           BOUNDARY LINE SEGMENT.  FIND OFFSET OF POINT OF
C           INTERSECTION.
 
            XI = XL + (ZI - ZL)*(XR - XL)/(ZR - ZL)
            IF(SNFLG.EQ.1) THEN
              IF(XR.GT.XL) THEN
                SNI = SNL + (XI - XL)*(SNR - SNL)/(XR - XL)
              ELSE
C               OVERHANGING BANK?
                SNI = 1.0
              ENDIF
            ELSE
C             IN THIS CASE THE SINUOSITY APPLIES TO THE WHOLE
C             LINE SEGMENT.  ALSO THE SINUOSITY OF THE LINE
C             SEGMENT IS GIVEN BY THE SINUOSITY OF ITS LEFT POINT.
 
              SNI = SNL
            ENDIF
            IF(ZL.LE.ZR) THEN
              DX = XI - XL
              XR = XI
              SNR = SNI
              ZR = ZI
              DZ = ZI - ZL
              YR = 0.0
              YL = DZ
            ELSE
              DX = XR - XI
              XL = XI
              SNL = SNI
              ZL = ZI
              DZ = ZI - ZR
              YL = 0.0
              YR = DZ
            ENDIF
            P =  PVEC(J-1)
            GOTO 890
 700      CONTINUE
 
C           LINE SEGMENT BELOW WATER SURFACE
 
            DX = XR - XL
            DZ = ZR - ZL
            YL = ZI - ZL
            YR = ZI - ZR
            P = PVEC(J-1)
 890      CONTINUE
 
C         COMPUTE THE SINUOSITY VALUES IF NEEDED.
 
C         The following options exist: 0- no sinuousity elements computed,
C           1- a sinuosity value is given at each point on the boundary
C              and linear variation is assumed between adjacent points.
C           2- a sinuosity value is given at each point on the boundary
C              and the value is a constant for each line segment following
C              the point.  Piecewise constant variation.  This option
C              is computed in COMPEL.
 
C
          IF(SNFLG.EQ.1) THEN
            IF(DX.GT.0.0.AND.MAX(YL,YR).GT.0.0) THEN
              M = (ZR - ZL)/DX
              IF(NVAR(IS).GT.0) THEN
C               N DEFINED ON SUBSECTION BASIS
                C = NFAC/(N(IS)*(1. + M**2)**T1D3)
                IF(NBX.EQ.1) THEN
C                 COMPUTE THE CONSTANTS FOR LOCAL FLUX COEF.
                  C1 = NFAC/(N(IS)*(1. + M**2)**0.083333333333333D0)
                  C2 = 75.*GRAV/(4.*C1**2)
                  C3 = -125.*GRAV**1.5/(4.*C1**3)
                  C4 = 6.251*GRAV/C1**2
                ENDIF
              ELSE
                C = NFAC/(LSN(J-1)*(1. + M**2)**T1D3)
                IF(NBX.EQ.1) THEN
C                 COMPUTE THE CONSTANTS FOR LOCAL FLUX COEF.
                  C1 = NFAC/(LSN(J-1)*(1. + M**2)**0.083333333333333D0)
                  C2 = 75.*GRAV/(4.*C1**2)
                  C3 = -125.*GRAV**1.5/(4.*C1**3)
                  C4 = 6.251*GRAV/C1**2
                ENDIF
 
              ENDIF
 
C             COMPUTE THE SUM FOR MA
 
              SUMMA = SUMMA + (XR - XL)*(SNL*YL +
     A                  0.5*(SNL*YR + SNR*YL) + SNR*YR)/3.0
 
              MU = (SNR - SNL)/DX
 
C             SIMPLE CLOSED FORM EXPRESSIONS FOR THE OTHER INTEGRALS
C             ARE NOT AVAILABLE WHEN THE SINUOUSITY VARIATION IS
C             PIECEWISE LINEAR.  THEREFORE USE NUMERICAL INTEGRATION.
C             WE USE A GAUSSIAN RULE.
 
              XMID =.5D0*(DBLE(XL) + XR)
              HALFDX = .5D0*DBLE(DX)
 
C              WRITE(STD6,*) ' M=',M,' MU=',MU
 
              DO 200 IGS=1,NGS
                W = HALFDX*WGS(IGS)
 
C               COMPUTE THE BASIC VALUES FOR THE CURRENT LINE SEGMENT.
 
                S = HALFDX*XGS(IGS) + XMID
                H = YL - M*(S - XL)
                IF(H.LT.0.D0) THEN
                  IF(ABS(DX).LT.0.10) THEN
                    H = 0.D0
                  ELSE
                    WRITE(STD6,*) ' H=',H,' <0'
                    WRITE(STD6,*) ' NGS=',NGS,' IGS=',IGS,
     A                                ' XGS(IGS)=',XGS(IGS)
                    WRITE(STD6,*) ' XR=',XR,' XL=',XL
                    WRITE(STD6,*) ' ZR=',ZR,' ZL=',ZL
                    WRITE(STD6,*) ' S=',S,' YL=',YL
                    WRITE(STD6,*) ' M=',M, ' XMID=',XMID
                    WRITE(STD6,*) ' DX=',DX,' HALFDX=',HALFDX
                    STOP 'Abnormal stop.  Error? found'
                  ENDIF
                ENDIF
                SNLOC = SNL + MU*(S - XL)
                SQRTSN = SQRT(SNLOC)
                K = C*H**1.666666666666666D0
                Q = K/SQRTSN
 
C                WRITE(STD6,*) ' S=',S,' H=',H,' SN=',SNLOC
C                WRITE(STD6,*) ' K=',K,' Q=',Q
 
C               DO SUMS FOR COMPUTATION OF MQ
 
                SUMMQ = SUMMQ + W*SNLOC*Q
                SUMQ = SUMQ + W*Q
 
C               DO SUMS FOR COMPUTATION OF SUBSECTION CONVEYANCE
C               ADJUSTMENT FOR SINUOUSITY
 
                QS(IS) = QS(IS) + W*Q
                KS(IS) = KS(IS) + W*K
 
                IF(NBFLAG.EQ.1) THEN
C                 COMPUTE THE ADDITIONAL SUMS REQUIRED FOR THE
C                 NEWBETA OPTION.  NOTE THAT CONSTANT FACTORS ON
C                 THE DERIVATIVES ARE APPLIED TO THE COMPLETED
C                 SUM OUTSIDE THE LOOP
                  IF(H.GT.0.D0) THEN
                    V = Q/H
                  ELSE
                    V = 0.0
                  ENDIF
C                  WRITE(STD6,*) ' V=',V
                  IF(NBX.EQ.0.OR.H.EQ.0.D0) THEN
C                   ASSUME LOCAL FLUX COEFFICIENTS = 1.0
                    FM = Q*V
                    FE = FM*V
                    DFM = V*V
                    DFE = DFM*V
                  ELSE
C                   ESTIMATE LOCAL FLUX COEFFICIENTS
 
                    ALOC = 1. + C2/H**0.33333333333333D0
     A                        + C3/SQRT(H)
                    DALOC = -C2/(3.*H**1.3333333333333D0)
     A                      -C2/(2.*H**1.50)
 
                    BLOC = 1. + C4/(H**0.33333333333333D0)
                    DBLOC = -C4/(3.*H**1.3333333333333D0)
 
                    FM = Q*V
                    FE = FM*V
                    DFM = V*V
                    DFE = DFM*V
 
                    DFM = DBLOC*FM + DFM*BLOC
                    DFE = DALOC*FE + DFE*ALOC
                    FM = FM*BLOC
                    FE = FE*ALOC
                  ENDIF
 
                  SUMDQ = SUMDQ + W*V
                  SUMFM = SUMFM + W*FM
                  SUMFE = SUMFE + W*FE
                  SUMDFM = SUMDFM + W*DFM
                  SUMDFE = SUMDFE + W*DFE
                ENDIF
 
 200          CONTINUE
 
            ENDIF
          ELSEIF(SNFLG.EQ.2) THEN
C           PICK VALUES OF SINUOSITY AND ASSIGN TO SUBSECTIONS.
 
            SBSN(IS) = SNL
          ENDIF
 
 
          IF(BETOPT(1:7).EQ.'NEWBETA'.AND.SNFLG.NE.1) THEN
C           COMPUTE THE VALUES FOR THE NEW METHOD FOR ESTIMATING
C           BETA AND ALPHA FOR AN OPEN CHANNEL CROSS SECTION.
 
            DX = XR - XL
            IF(DX.GT.0.0) THEN
              M = (ZR - ZL)/DX
              IF(NVAR(IS).GT.0) THEN
C               N DEFINED ON SUBSECTION BASIS
                C = NFAC/(N(IS)*(1. + M**2)**T1D3)
              ELSE
                C = NFAC/(LSN(J-1)*(1. + M**2)**T1D3)
              ENDIF
              IF(SNFLG.EQ.2) THEN
C               WE HAVE POSSIBLY NON-UNITARY VALUES OF SINUOUSITY.
C               THE VALUE FOR THE LINE SEGMENT IS AT THE LEFT END.
C               ADJUST C SO THAT THE LOCAL VELOCITY WILL TAKE INTO
C               ACCOUNT THE SINUOSITY
                C = C/SQRT(SNL)
              ENDIF
 
              IF(ABS(M).GT.1.E-6) THEN
                SUMQ = SUMQ
     A                  - 3.*C*(DBLE(YR)**T8D3 - DBLE(YL)**T8D3)/(8.*M)
 
                SUMFM = SUMFM
     A               - C**2*(0.3*(DBLE(YR)**T10D3 - DBLE(YL)**T10D3))/M
 
                SUMFE = SUMFE
     A                      - C**3*(0.25*(DBLE(YR)**4 - DBLE(YL)**4))/M
 
                SUMDQ = SUMDQ - C*(DBLE(YR)**T5D3 - DBLE(YL)**T5D3)/M
                SUMDFM = SUMDFM - C**2*(DBLE(YR)**T7D3
     A                                 - DBLE(YL)**T7D3)/M
                SUMDFE = SUMDFE - C**3*(DBLE(YR)**3 - DBLE(YL)**3)/M
              ELSE
                SUMQ = SUMQ + DX*C*YL**T5D3
                SUMFM = SUMFM + DX*C**2*YL**T7D3
                SUMFE = SUMFE + DX*C**3*YL**3
                SUMDQ = SUMDQ + DX*C*T5D3*YL**T2D3
                SUMDFM = SUMDFM + DX*C**2*T7D3*YL**T4D3
                SUMDFE = SUMDFE + DX*C**3*3*YL**2
              ENDIF
            ENDIF
          ENDIF
 900  CONTINUE
 
 
C     ADJUST SUMS WHEN SNFLG = 1 AND NBFLAG = 1
 
      IF(SNFLG.EQ.1.AND.NBFLAG.EQ.1) THEN
 
        SUMDQ = T5D3*SUMDQ
        SUMDFM = T7D3*SUMDFM
        SUMDFE = 3*SUMDFE
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   COMPEL
     I                   (ZI, NPNT, NSUB, NAVM, X, Z, SB, NFAC, BETOPT,
     I                    SNFLG, LSN, NVAR, NATY, YATN, NNY, SN, WRN557,
     M                    KOLD, TSOLD,
     O                    N, XSV)
 
C     + + + PURPOSE + + +
C     Compute cross sectional elements at a given elevation.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER NAVM, NPNT, NSUB, SNFLG, WRN557
      INTEGER NNY(NSUB), NVAR(NSUB), SB(NPNT)
      REAL KOLD(NSUB), LSN(NPNT), N(NSUB), NATY(9,NSUB), NFAC, SN(NPNT),
     A     TSOLD(NSUB), X(NPNT), XSV(PMXELM), YATN(9,NSUB), Z(NPNT), ZI
      CHARACTER BETOPT*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ZI     - Current water surface elevation for computing elements
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     NAVM   - Flag for averaging roughness
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     BETOPT - Option for computing flux coefficients and critical flow
C     SNFLG  - Flag for sinuousity computations
C     LSN    - Line segment Manning's n value
C     NVAR   - Flag for variation of Manning's n in each subsection
C     NATY   - Mannings's n value at depth in YATN
C     YATN   - Depth values for the Manning's n values in NATY
C     NNY    - Number of values for Manning's n variation with depth
C     SN     - Sinuousity at a point on a cross section boundary
C     WRN557 - flag to suppress WRN:557 after first time
C     KOLD   - Previous value of conveyance
C     TSOLD  - Old top width in each subsection
C     N      - Manning's n values
C     XSV    - Vector of various elements of cross section
 
C     Added May 22, 1998: definition of XSV(*) contents
c     Extended 13 March 2003

C     Offset    Value
C     1         Maximum depth-y
C     2         Top width
C     3         Area
C     4         First moment of area about water surface
C     5         Square root of conveyance- kh
C     6         Beta
C     7         Alpha
C     8         dBeta/dy
C     9         dAlpha/dy
C     10        Critical flow from momentum
C     11        Critical flow from energy
C     12        Critcal flow assuming Alpha=beta=1
C     13        Critical flow that is selected by user: 10, 11, or 12
C     14        MA- correction of volumes for sinuousity
C     15        MQ- correction of momentum for sinuosity
C     16        Total wetted perimeter- added May 22, 1998
C     17        Average Manning's n value for the cross section
c     18        dkh/dy computed by cubic-spline fit
c     19        dalpha/dy computed by cubic-spline fit
c     20        dbeta/dy computed by cubic-spline fit
c     21        dma/dy computed by cubic-spline fit
c     22        dmq/dy computed by cubic-spline fit
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xscom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL ALPHA, ALPSUM,  ASUM, BETA, BETSUM, DEN, DROP,
     A     KS(PMXSUB), KSUM, KT, KVEC(PMXSUB), NS(PMXSUB), PS(PMXSUB),
     B     PSUM, QS(PMXSUB), R, SBSN(PMXSUB), TEMP, TS(PMXSUB), TSLOT,
     C     TSUM, WN, YBASUM,  YSMX(PMXSUB), MXSLOT
      DOUBLE PRECISION AS(PMXSUB), YBS(PMXSUB)
      DOUBLE PRECISION SUMDFE, SUMDFM, SUMDQ, SUMFE, SUMFM, SUMMA,
     A                 SUMMQ, SUMQ
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FBASEL
 
C     + + + OUTPUT FORMATS + + +
 10   FORMAT(/,' *WRN:505* Decrease in conveyance in subsection',I4,
     A   ' at elevation=',F10.3,/,11X,' decrease=',F6.1,' per cent')
 11   FORMAT(' *ERR:632* At elevation=',F10.3,' subsection area <= 0',
     A   /,11X,' in subsection number=',I5,
     B   /,11X,' Check for input error or switch to NAVM= 1.')
 12   FORMAT(/,' *WRN:557* At elevation=',F10.3,' subsection=',I3,
     A ' is a single vertical',/,11X,' line segment.  Acts as a',
     B ' frictionless wall.')
50    FORMAT(10X,'Above warning given only once but may apply',
     A       ' many times.')
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        MXSLOT = 0.07
      ELSE
        MXSLOT = 0.02134
      ENDIF
      TSLOT = 1.0001*SLOT
C     FIND BASIC ELEMENTS AT THE GIVEN ELEVATION
 
      CALL FBASEL
     I           (SNFLG, ZI, NPNT, NSUB, X, Z, SB, LSN, NVAR, SN, NATY,
     I            YATN, NNY, NFAC, BETOPT,
     O            N, TS, PS, AS, YBS, NS, SUMQ, SUMFM, SUMFE, SUMDQ,
     O            SUMDFM, SUMDFE, SUMMA, SUMMQ, YSMX, SBSN, QS, KS)
 
C      WRITE(STD6,*) ' '
C      WRITE(STD6,*) ' COMPEL: BETOPT=',BETOPT
C      WRITE(STD6,*) ' SUMQ=',SUMQ,' SUMFM=',SUMFM,' SUMFE=',SUMFE
C      WRITE(STD6,'('' COMPEL: ZI='',F10.4)') ZI
C      WRITE(STD6,'(''    OFFSET ELEVATION  SUB'')')
C      DO 40 I=1,NPNT
C        WRITE(STD6,'(2F10.4,I5)') X(I), Z(I), SB(I)
C40    CONTINUE
C      DO  50 I=1,NSUB
C        WRITE(STD6,'(4F10.4)') TS(I), PS(I), AS(I), N(I)
C50    CONTINUE
C      IF(SNFLG.EQ.2) THEN
C        WRITE(STD6,*) ' '
C        WRITE(STD6,*) ' PWC SINU:',(SBSN(J),J=1,NSUB)
C        WRITE(STD6,*) ' SUB AREA:',(AS(J),J=1,NSUB)
C        WRITE(STD6,*) ' N        :',(N(J),J=1,NSUB)
C        WRITE(STD6,*) ' NVAR    :',(NVAR(J),J=1,NSUB)
C        WRITE(STD6,*) ' PS      :',(PS(J),J=1,NSUB)
C      ENDIF
      ASUM = 0.
      TSUM = 0.
      PSUM = 0.
      YBASUM = 0.
      KSUM = 0.
      ALPSUM = 0.
      BETSUM = 0.
      WN = 0.
      DO 700 J=1,NSUB
        KVEC(J) = 0.0
        TSUM = TSUM + TS(J)
        PSUM = PSUM + PS(J)
 
        IF(NAVM.EQ.1) WN = WN + PS(J)*N(J)
 
C       ZERO SUB-SECTION AREA CAN INDICATE NO WATER IN SUBSECTION
C       AT ELEVATION=ZI OR IT CAN MEAN THAT AN OVERHANG EXISTS WHICH
C       EXACTLY BALANCES POSTIVE AND NEGATIVE INCREMENTS OF AREA SO THAT
C       AREA IS ZERO AND PERIMETER IS NOT ZERO.  A PERIMETER BEING ZERO
C       MEANS THAT THE SUBSECTION IS ABOVE THE CURRENT WATER ELEVATION,
C       ZI.
 
        IF(PS(J).NE.0.0) THEN
          ASUM = ASUM + AS(J)
          YBASUM = YBASUM + YBS(J)
          IF(NAVM.EQ.0) THEN
            IF(AS(J).LT.0.0) THEN
C             A SUBSECTION AREA IS <=0.  MUST USE NAVM=1. EITHER
C             INPUT ERROR OR LARGE OVERHANG IN CROSS SECTION.
              WRITE(LOUT,11)  ZI,J
              EFLAG = 1
              AS(J) = ABS(AS(J))
            ELSEIF(AS(J).EQ.0.0.AND.TS(J).EQ.0.0) THEN
C             SUBSECTION CONSISTS OF ONLY ONE VERTICAL LINE SEGMENT.
C             WETTED PERIMETER DOES NOT AFFECT FLOW.
              IF(WRN557.EQ.1) THEN
                WRITE(LOUT,12) ZI,J
                WRN557 = 0
                WRITE(LOUT,50)
              ENDIF
              GOTO 700
            ELSEIF(AS(J).EQ.0.0) THEN
              GOTO 700
            ENDIF
            R = AS(J)/PS(J)
            KT= (NFAC*AS(J)*R**0.666667)/N(J)
 
            IF(SNFLG.EQ.2) THEN
C             COMPUTE THE PIECEWISE CONSTANT VALUES FOR THE
C             SINUOSITY DEPENDENT ELEMENTS.
 
C             ADJUST THE SUBSECTION CONVEYANCE TO REFLECT THE
C             EFFECT OF ASSUMING A CONSTANT DECLINE IN ELEVATION
C             OF THE TOTAL ENERGY LINE INSTEAD OF CONSTANT FRICTION
C             SLOPE FOR SUBSECTION.  AFFECTS ALL CONVEYANCE RELATED
C             VALUES.
 
              KT = KT/SQRT(SBSN(J))
              SUMMA = SUMMA + AS(J)*SBSN(J)
              SUMMQ = SUMMQ + KT*SBSN(J)
            ELSEIF(SNFLG.EQ.1) THEN
C             ADJUST THE SUBSECTION CONVEYANCE FOR THE WEIGHTED
C             VALUE OF 1/SQRT(SINUOUSITY).
 
              IF(KS(J).NE.0.0) THEN
                KT = KT*QS(J)/KS(J)
              ENDIF
        
 
            ENDIF
 
C           ACCOUNT FOR SLOT WHICH MAY BE PRESENT IN THE TOP OF
C           A CLOSED CONDUIT
            IF(KT.LT.KOLD(J)) THEN
              IF(TS(J).LE.TSLOT.AND.TSOLD(J).LE.TSLOT) THEN
                KT = KOLD(J)
              ENDIF
            ENDIF
            KVEC(J) = KT
 
 
            KSUM = KSUM + KT
 
            IF(IUSGS.EQ.1) THEN
C             USE USGS RELATIONSHIP BETWEEN N AND ALPHA FOR EACH
C             SUBSECTION.
              ALPHA = 14.8*N(J) + 0.884
              IF(ALPHA.GT.2.0) ALPHA = 2.0
              IF(ALPHA.LT.1.0) ALPHA = 1.0
              BETA = 1.0 + 0.3467*(ALPHA - 1.0)
            ELSE
              ALPHA = 1.0
              BETA = 1.0
            ENDIF
 
            TEMP = KT*KT/AS(J)
            BETSUM = BETSUM + BETA*TEMP
            ALPSUM = ALPSUM + ALPHA*KT*TEMP/AS(J)
          ENDIF
          TSOLD(J) = TS(J)
        ENDIF
 700  CONTINUE
 
      IF(NAVM.EQ.1) THEN
        IF(PSUM.GT.0.0) THEN
          KSUM = PSUM*(NFAC*ASUM*(ASUM/PSUM)**0.666667)/WN
        ELSE
          KSUM = 0.0
        ENDIF
      ENDIF
      IF(KSUM.LE.0.0) GOTO 710
        TEMP = ASUM/(KSUM*KSUM)
        BETSUM = BETSUM*TEMP
        ALPSUM = ALPSUM*TEMP*ASUM/KSUM
 710  CONTINUE
      IF(BETSUM.LT.1.0) BETSUM = 1.0
      IF(ALPSUM.LT.1.0) ALPSUM = 1.0
 
C     CHECK THE SUBSECTION VARIATION OF CONVEYANCE
 
      DO 800 I=1,NSUB
        IF(KVEC(I).LT.KOLD(I)) THEN
          IF(NOCM.EQ.0)  THEN
            DROP = ((KVEC(I) - KOLD(I))/KOLD(I))*100.0
            IF(TS(I).GT.MXSLOT.AND.ABS(DROP).GT.1.0) THEN
              WRITE(STD6,10) I, ZI, ABS(DROP)
            ENDIF
          ENDIF
        ENDIF
        KOLD(I) = KVEC(I)
 800  CONTINUE
 
C     STORE FINAL VALUES FOR RETURN
 
      XSV(2) = TSUM
      XSV(3) = ASUM
      XSV(4) = YBASUM
      XSV(5) = SQRT(KSUM)
C     Added May 22, 1998
      XSV(16) = PSUM
      IF(ASUM.GT.0.0) THEN
        XSV(17) = NFAC*ASUM*(ASUM/PSUM)**0.666667/KSUM
      ELSE
        XSV(17) = 0.0
      ENDIF
      
      IF(BETOPT.EQ.'OLDBETA') THEN
        XSV(6) = BETSUM
        XSV(7) = ALPSUM
        XSV(8) = 0.0
        XSV(9) = 0.0
        XSV(10) = 0.0
        XSV(11) = 0.0
        XSV(12) = 0.0
        IF(TSUM.GT.0.0) THEN
          XSV(13) = ASUM*SQRT(GRAV*ASUM/TSUM)
        ELSE
          XSV(13) = 0.0
        ENDIF
      ELSE
C       COMPUTE NEW VALUES
        IF(SUMQ.GT.0.0) THEN
          XSV(6) = ASUM*SUMFM/SUMQ**2
          XSV(7) = ASUM**2*SUMFE/SUMQ**3
          XSV(8) = TSUM*SUMFM/SUMQ**2 - 2.*ASUM*SUMFM*SUMDQ/SUMQ**3
     A                                + ASUM*SUMDFM/SUMQ**2
          XSV(9) = 2.*ASUM*TSUM*SUMFE/SUMQ**3
     A                 - 3.*ASUM**2*SUMFE*SUMDQ/SUMQ**4
     B                    + ASUM**2*SUMDFE/SUMQ**3
        ELSE
          XSV(6) = 1.0
          XSV(7) = 1.0
          XSV(8) = 0.0
          XSV(9) = 0.0
          XSV(13) = 0.0
        ENDIF
 
C       COMPUTE CRITICAL FLOW
        IF(ASUM.GT.0.0) THEN
          DEN = XSV(6)*TSUM - XSV(8)*ASUM
          IF(DEN.LE.0.0) THEN
            XSV(10) = -9876541.0
          ELSE
            XSV(10) = ASUM*SQRT(GRAV*ASUM/DEN)
          ENDIF
 
          DEN = XSV(7)*TSUM - 0.5*XSV(9)*ASUM
          IF(DEN.LE.0.0) THEN
            XSV(11) = -9876543.0
          ELSE
            XSV(11) = ASUM*SQRT(GRAV*ASUM/DEN)
          ENDIF
 
          IF(TSUM.GT.0.0) THEN
            XSV(12) = ASUM*SQRT(GRAV*ASUM/TSUM)
          ELSE
            XSV(12) = 0.0
          ENDIF
          IF(XSV(10).GE.0.0.AND.XSV(11).GE.0.0) THEN
            IF(BETOPT.EQ.'NEWBETAM'.OR.BETOPT.EQ.'NEWBETAX') THEN
              XSV(13) = XSV(10)
            ELSEIF(BETOPT.EQ.'NEWBETAE') THEN
              XSV(13) = XSV(11)
            ELSE
              XSV(13) = SQRT(XSV(10)*XSV(11))
            ENDIF
          ELSE
            XSV(13) = XSV(12)
          ENDIF
        ELSE
          XSV(10) = 0.0
          XSV(11) = 0.0
          XSV(12) = 0.0
          XSV(13) = 0.0
        ENDIF
 

        IF(XSV(6).LT.1.0.OR.XSV(7).LT.1.0) THEN
          WRITE(STD6,*) ' *ERR:615* NEWBETA FAILURE. VALUE < 1'
        ENDIF
      ENDIF
 
      IF(SNFLG.EQ.1) THEN
C       COMPUTE THE ELEMENTS BASED ON SINUOSITY THAT VARIES
C       PIECEWISE LINEARLY BETWEEN POINTS ON THE CROSS SECTION
C       BOUNDARY.
 
        IF(ASUM.EQ.0.0) THEN
C         ASSUME ZERO DEPTH
          XSV(14) = 1.0
          XSV(15) = 1.0
        ELSE
          XSV(14) = SUMMA/ASUM
          XSV(15) = SUMMQ/SUMQ
        ENDIF
      ELSEIF(SNFLG.EQ.2) THEN
C       COMPUTE THE ELEMENTS BASED ON SINUOSITY THAT IS
C       PIECEWISE CONSTANT IN EACH SUBSECTION OF THE
C       CROSS SECTION.
 
        IF(ASUM.EQ.0.0) THEN
          XSV(14) = 1.0
          XSV(15) = 1.0
        ELSE
          XSV(14) = SUMMA/ASUM
          XSV(15) = SUMMQ/KSUM
        ENDIF
      ELSE
        XSV(14) = 1.0
        XSV(15) = 1.0
      ENDIF
 
 
c      WRITE(STD6,'(1P7E13.5)') (XSV(J), J=2,7)
 
      RETURN
      END

C
C
C
      SUBROUTINE FIND_GEN_ST_LINE
     I                           (X1, Y1, X2, Y2,
     O                            A, B, C)

C     Find the general equation of a straight line.

      REAL X1, Y1, X2, Y2, A, B, C

C***********************************************************************
      IF(X2.NE.X1) THEN
        IF(Y2.NE.Y1) THEN
C         General case here.  Line is neither horizontal nor vertical.

          A = Y2 - Y1
          B = X1 - X2
          C = X2*Y1 - X1*Y2
        ELSE
C         Line is horizontal.  
          A = 0.0
          B = 1.0
          C = -Y1
        ENDIF
      ELSE
        IF(Y2.NE.Y1) THEN
C         The line is vertical
          A = 1.0
          B = 0.0
          C = -X1
        ELSE
C         The line does not exist!
          A = 0.0
          B = 0.0
          C = 0.0
        ENDIF
      ENDIF
      RETURN
      END        

C
C
C
      SUBROUTINE   LINE_INTERSECT
     I                           (A1, B1, C1, A2, B2, C2,
     O                            EXIST, X, Y)

C     Find the intersection point of two straigt lines given
C     in general linear equation form. 

      INTEGER EXIST

      REAL A1, B1, C1, A2, B2, C2, X, Y,
     A     TP
C***********************************************************************
      TP = A1*B2 - A2*B1
      IF(TP.NE.0.0) THEN

        EXIST = 1
        X = (C2*B1 - C1*B2)/TP
        Y = (C1*A2 - C2*A1)/TP
      ELSE
        EXIST = 0
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   FIND_CROSS
     I                 (STDOUT, NPNT, MIN_AT, DIR, NLIM, X, Z, 
     I                  SLOT_XM, SLOT_ZM, SLOT_X, SLOT_Z,        
     O                 ITYPE, ILCROSS, XCROSS, ZCROSS)      


C     Seek an intersection between the cross section boundary and
C     the given slot boundary. 

      INTEGER ITYPE, MIN_AT, NLIM, NPNT, DIR, ILCROSS, STDOUT

      REAL X(NPNT), Z(NPNT), SLOT_XM, SLOT_ZM, SLOT_X, SLOT_Z,
     A     XCROSS, ZCROSS


C     Local

      INTEGER I, J, EXIST

      REAL X1, Z1, X2, Z2, ASLOT, BSLOT, CSLOT,
     A        A, B, C, XC, ZC, EPS

      DATA EPS/0.0005/

C**********************************************************************
C     Find general equation for the bounding line of the slot.
      CALL FIND_GEN_ST_LINE
     I                     (SLOT_XM, SLOT_ZM, SLOT_X, SLOT_Z,
     O                      ASLOT, BSLOT, CSLOT)

c      WRITE(STDOUT,*) ' SLOT_XM=',SLOT_XM,' SLOT_ZM=',SLOT_ZM
c      WRITE(STDOUT,*) ' SLOT_X=',SLOT_X,' SLOT_Z=',SLOT_Z
c      WRITE(STDOUT,*) ' ASLOT=',ASLOT,' BSLOT=',BSLOT,' CSLOT=',CSLOT
C     Set the intersection type to: no intersection found.
C      WRITE(STDOUT,*) ' FIND_CROSS: MIN_AT=',MIN_AT,' NPNT=',NPNT
      ITYPE = 0
      I = MIN_AT
      X1 = X(I)
      Z1 = Z(I)
      J = I
100   CONTINUE
                
        J = J + DIR
        IF(J.LT.1.OR.J.GT.NPNT) GOTO 500
        X2 = X(J)
        Z2 = Z(J)
            
C       Find general equation for the line segment on the boundary.

c        WRITE(STDOUT,*) ' J=',J
c        WRITE(STDOUT,*) 'X1=',X1,' Z1=',Z1
c        WRITE(STDOUT,*) 'X2=',X2,' Z2=',Z2

        CALL FIND_GEN_ST_LINE
     I                       (X1, Z1, X2, Z2,
     O                         A, B, C)
      
c        WRITE(STDOUT,*) ' A=',A,' B=',B,' C=',C

        CALL LINE_INTERSECT
     I                     (ASLOT, BSLOT, CSLOT, A, B, C,
     O                      EXIST, XC, ZC)

c        WRITE(STDOUT,*) ' EXIST=',EXIST,' XC=',XC,' ZC=',ZC
        IF(EXIST.EQ.0) THEN
C         No useful intersection.

          ITYPE = 0
        ELSE
C         Is the intersection on the boundary line segment?
          IF( MIN(X1,X2)-EPS.LE.XC.AND. XC.LE.MAX(X1,X2)+EPS) THEN
            if(min(z1,z2)-eps.le.zc.and.zc.le.max(z1,z2)+eps) then
              ITYPE = 1
            endif
          ENDIF
        ENDIF

        IF(ITYPE.EQ.1) THEN
          XCROSS = XC
          ZCROSS = ZC
C         Set the index to the left end of the boundary line
C         segment. 
          IF(DIR.GT.0) THEN
            ILCROSS = J - 1
          ELSE
            ILCROSS = J
          ENDIF
          GOTO 500
        ENDIF
        X1 = X2
        Z1 = Z2
        GOTO 100

500   CONTINUE
      RETURN
      END        
          
C
C
C
      SUBROUTINE   SHIFT_VALUES
     I                         (STDOUT, IOLD, INEW, 
     M                          NPNT, X, Z, SB, LSN, SN)

C     Shift cross section boundary items to make room for inserting
C     the slot description.

C     iold  -  index into x(*) for the starting point of the sequence
c              of points that is to be shifted.
c     inew =  index of the location where the contents at iold are to 
c             be after the shift. 
c     npnt -  index to the last point in x(*) on entry.  On exit
c             index to the new location of the last point.  
c     x, z, sb, lsn, and sn -- contain the section boundary descriptive
C            data. 

      INCLUDE 'arsize.prm'

      INTEGER STDOUT, IOLD, INEW, NPNT, SB(PMXPNT)

      REAL  X(PMXPNT), Z(PMXPNT), LSN(PMXPNT), SN(PMXPNT)


C     Local

      INTEGER I, IDIFF, INC
C     ****************************FORMATS*******************************
50    FORMAT(/,' *ERR:723* Maximum number of points on cross section',
     A      ' boundary=',I5, ' exceeded in adding invert slot.')
C***********************************************************************

      IDIFF = INEW - IOLD
      IF(IDIFF.GT.0) THEN
        IS = NPNT
        IE = IOLD
        INC = -1
      ELSEIF(IDIFF.LT.0) THEN
        IS = IOLD
        IE = NPNT
        INC = 1
      ELSE
C       Force loop to be skipped.
        IS = NPNT+1
        IE = IOLD
        INC = 1
      ENDIF
      DO 100 I=IS, IE, INC
        X(I+IDIFF) = X(I)
        Z(I+IDIFF) = Z(I)
        SB(I+IDIFF) = SB(I)
        LSN(I+IDIFF) = LSN(I)
        SN(I+IDIFF) = SN(I)
100   CONTINUE

      NPNT = NPNT + IDIFF
      IF(NPNT.GT.PMXPNT) THEN
        WRITE(STDOUT,50) PMXPNT
        STOP 'Abnormal stop' 
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   ADD_SLOT
     I                      (STDOUT, in_ss_nvar,
     M                       NPNT, NSUB, N, ZMIN, X, Z, SB, LSN, SN,
     o                       added)

C     Potentially add a bottom slot to a cross section description. 
C     The slot size is in common block ABSLOT.

C     See CXSTAB for definition of arguments.

      INCLUDE 'arsize.prm'

      INTEGER added,STDOUT, NSUB, SB(PMXPNT), in_ss_nvar

      REAL X(PMXPNT), Z(PMXPNT), 
     A     LSN(PMXPNT), SN(PMXPNT), N(PMXSUB)
       
      INCLUDE 'abslot.cmn'
      INCLUDE 'xtadd.cmn'

C     Local

      CHARACTER*1  SMATCH, EMATCH

      INTEGER I, J, MIN_KNT, MIN_PNT(PMXPNT), MIN_AT,
     A        ILCROSS, SAVED_SUB, IS, IE, EFLAG

      REAL  EPS, WIDTH(PMXPNT), MINOFF, MAXOFF, BOTTOM, XL, ZL,
     A      XR, ZR, SLOT_XL, SLOT_ZL, SLOT_XR, XM, ZM,
     B      SLOT_XM, SLOT_ZM

      DATA EPS/0.0005/
C     ********************************FORMATS**************************
 50   FORMAT(' *ERR:724* Maximum number of subsections=',I5,
     A       ' exceeded in adding an invert slot.')
 54   format(/,'*PROBLEM* Subset operations on cross sections with',
     a  ' roughness varying with depth not yet supported.')
C**********************************************************************
        
      IF(SLOT_PRESENT.EQ. 0) THEN
        SLOT_DEPTH = 0.0
        added = 0
        RETURN
      elseif(slot_present.ne.1) then
c       Exponential slot being processed.
        return
      ENDIF

      if(in_ss_nvar /= 0)  then
        write(stdout,54)
        stop 'Abnormal stop. Feature not yet supported.'
      endif

      added = 1
C     Compute the slot vertical extent. 
      IF(YSLOT.GT.0.0) THEN
C       Add a fixed size slot to the bottom of each cross section.
C       Redefine ESLOT to be the new invert.
        ESLOT = ZMIN - YSLOT
        SLOT_DEPTH = YSLOT
      ELSE
C       Add a slot that has an invert at a fixed elevation to 
C       each cross section.        
        SLOT_DEPTH = ZMIN - ESLOT
        IF(SLOT_DEPTH.LT.0.0) THEN
C         No slot needed
          SLOT_DEPTH = 0.0
          RETURN
        ENDIF
      ENDIF
C     Change sign of slot so that we can add to get the result.
C     SLOT is now the distance to the bottom of the slot from the min
C     elevation in the cross section. 
      SLOT_DEPTH = -SLOT_DEPTH

C     We know the minimum elevation on entry.  Do another search to 
C     describe the various extremes.

C      WRITE(STDOUT,*) ' In ADD_SLOT: ZMIN=',ZMIN
      MIN_KNT = 0
      MINOFF = 1.E20
      MAXOFF = -1.E30
      DO 100 I=1,NPNT
        IF(ABS(Z(I) - ZMIN).LE.EPS) THEN
C         We have a match. 
          MIN_KNT = MIN_KNT + 1
          MIN_PNT(MIN_KNT) = I
C         Does the point to the right, if it exists, have essentially 
C         the same elevation?
          IF(I.LT.NPNT) THEN
            IF(ABS(Z(I+1) - ZMIN).LE.EPS) THEN
C             Yes.  We have an essentially horizontal segment at
C             the minimum elevation.  
              WIDTH(MIN_KNT) = X(I+1) - X(I)
            ELSE
C             Bottom width to right is zero.
              WIDTH(MIN_KNT) = 0.0
            ENDIF
          ELSE
C           Point to the right does not exist.  Bottom width to
C           right of point is zero. 

            WIDTH(MIN_KNT) = 0.0
          ENDIF
        ENDIF
C       Find the extreme offsets. 
        MAXOFF = MAX(MAXOFF,X(I))
        MINOFF = MIN(MINOFF,X(I))
        
100   CONTINUE
      IF(MIN_KNT.EQ.0) THEN
        WRITE(STDOUT,*) ' BUG in ADD_SLOT.  No minimum match found'
        STOP 'Abnormal stop' 
      ENDIF

      WRITE(STDOUT,89)
89    FORMAT('  Check of min search:',/,'    Pointer     Width')
      DO 1234 I=1,MIN_KNT
        WRITE(STDOUT,90) MIN_PNT(I), WIDTH(I)
1234  CONTINUE
90    FORMAT(' ',I10,F10.3)

      IF(WSLOT.GE.MAXOFF-MINOFF) THEN
        WRITE(STDOUT,*) ' *ERR:725*  Bottom slot wider than section!'
        STOP 'Abnormal stop' 
      ENDIF

C     Now select the point at which to place the slot.  We have no 
C     choice if MIN_KNT = 1!

      IF(MIN_KNT.EQ.1) THEN
        MIN_AT = MIN_PNT(1)
        BOTTOM = WIDTH(1)
      ELSE
        BOTTOM = -1.E30
        DO  110 J=1,MIN_KNT
          IF(WIDTH(J).GT.BOTTOM) THEN
            MIN_AT = MIN_PNT(J)
            BOTTOM = WIDTH(J)
          ENDIF      
110     CONTINUE
      ENDIF

C      WRITE(STDOUT,*) ' MIN_AT=',MIN_AT,' BOTTOM=',BOTTOM

      XM = X(MIN_AT)
      ZM = Z(MIN_AT)

C     Now add the new subsection and adjust the vectors.
      NSUB = NSUB + 1
      IF(NSUB.GT.PMXSUB) THEN
        WRITE(STDOUT,50) PMXSUB
        STOP 'Abnormal stop'
      ENDIF

C     Set the new minimum to the slot invert elevation.

      ZMIN = ESLOT

      IF(BOTTOM.GT.WSLOT) THEN
C       If BOTTOM is > WSLOT put the slot to the right of the 
C       selected minimum point.  It will always fit.  Also the
C       bottom of the cross section is known to be essentially
C       horizontal. 

C       Save the subsection value at the min point because
C       it will be needed after the slot is inserted.
        SAVED_SUB = SB(MIN_AT)
        
        
        CALL SHIFT_VALUES
     I                    (STDOUT, MIN_AT+1, MIN_AT+3, 
     M                     NPNT, X, Z, SB, LSN, SN)

        
        SB(MIN_AT) = NSUB
        LSN(MIN_AT) = N(NSUB)

        X(MIN_AT+1) = XM + 0.5*WSLOT
        Z(MIN_AT+1) = ESLOT
        SB(MIN_AT+1) = NSUB
        SN(MIN_AT+1) = 1.
        LSN(MIN_AT+1) = N(NSUB)

        X(MIN_AT+2) = XM + WSLOT
        Z(MIN_AT+2) = ZM
        SB(MIN_AT+2) = SAVED_SUB
        SN(MIN_AT+2) = 1.
        LSN(MIN_AT+2) = N(SAVED_SUB)
        

      ELSE
C       The bottom is either not horizontal or horizontal and
C       too short to fit the slot requested.  We use a bit of
C       trial and error to fit the slot to the cross section.
C       Start by placing the slot centerline at the minimum
C       point.  Get the  min point local values.


        SLOT_ZM = ESLOT
        SLOT_XM = XM

        SLOT_XL = XM - 0.5*WSLOT
        SLOT_ZL = ZM

C       Seek an intersection of the slot boundary with the cross
C       section boundary on the left. 


C        WRITE(STDOUT,*) ' CALLING FIND_CROSS 1'

        CALL FIND_CROSS
     I                 (STDOUT, NPNT, MIN_AT, -1, 1, X, Z, 
     I                  SLOT_XM, SLOT_ZM, SLOT_XL, SLOT_ZL,
     O                  ITYPE, ILCROSS, XL, ZL)      
        
C        WRITE(STDOUT,*) ' Left cross=',ILCROSS,' ITYPE=',ITYPE

        IF(ITYPE.EQ.0) THEN
C         No intersection.  Move slot to the right so that its
C         left point matches the minimum point and then seek 
C         an intersection on the right. 

          XL = XM
          ZL = ZM
          ILCROSS = MIN_AT

          SLOT_XM = XM + 0.5*WSLOT
          SLOT_XR = XM + WSLOT

          CALL FIND_CROSS
     I                   (STDOUT, NPNT, MIN_AT, 1, NPNT, X, Z, 
     I                    SLOT_XM, SLOT_ZM, SLOT_XR, SLOT_ZL,
     O                    ITYPE, IRCROSS, XR, ZR)      

          IF(ITYPE.EQ.0) THEN
            WRITE(STDOUT,*) ' BUG: No intersection found on right',
     A                     ' when one must exist.'
            STOP ' Abnormal stop' 
          ENDIF

        ELSE
C         Intersection found on the left.  Seek one on the right 
C         with the slot centerline on the minimum point. 
          SLOT_XR = XM + 0.5*WSLOT
          CALL FIND_CROSS
     I                   (STDOUT, NPNT, MIN_AT, 1, NPNT, X, Z, 
     I                    SLOT_XM, SLOT_ZM, SLOT_XR, SLOT_ZL,
     O                    ITYPE, IRCROSS, XR, ZR)      
        
          IF(ITYPE.EQ.0) THEN
C           No intersection on the right when there was one on
C           the left.  Move the the slot to the left so that
C           its right point matches the minimum and then
C           seek a new intersection on the left. 

            XR = XM
            ZR = ZM
            IRCROSS = MIN_AT 

            SLOT_XM = XM - 0.5*WSLOT
            SLOT_XL = XM - WSLOT
            CALL FIND_CROSS
     I                     (STDOUT, NPNT, MIN_AT, -1, 1, X, Z, 
     I                     SLOT_XM, SLOT_ZM, SLOT_XL, SLOT_ZL,
     O                     ITYPE, ILCROSS, XL, ZL)      

            IF(ITYPE.EQ.0) THEN
              WRITE(STDOUT,*) ' BUG: No intersection found on left',
     A                    ' when one must exist.'
              STOP ' Abnormal stop' 
            ENDIF

          ENDIF
        ENDIF

C       Set the starting and ending values with flags to denote
C       essentially exact matches with a boundary point. 

        IF(ABS(X(ILCROSS)-XL).LE.EPS) THEN
          IS = ILCROSS
          SMATCH = 'Y'
        ELSEIF(ABS(X(ILCROSS+1)-XL).LE.EPS) THEN
          IS = ILCROSS + 1
          SMATCH = 'Y'
        ELSE
          IS = ILCROSS
          SMATCH = 'N'
        ENDIF
        IF(ABS(X(IRCROSS)-XR).LE.EPS) THEN
          IE = IRCROSS
          EMATCH = 'Y'
        ELSEIF(ABS(X(IRCROSS+1)-XR).LE.EPS) THEN
          IE = IRCROSS + 1
          EMATCH = 'Y'
        ELSE
          IE = IRCROSS
          EMATCH ='N'
        ENDIF


C        WRITE(STDOUT,*) ' ADD_SLOT: IS=',IS,' IE=',IE
C        WRITE(STDOUT,*) ' SMATCH=',SMATCH,' EMATCH=',EMATCH

        SAVED_SUB = SB(IE)
C        WRITE(STDOUT,*) ' SAVED_SUB=',SAVED_SUB, 'N(SAVED_SUB)=',
C     A                  N(SAVED_SUB)

c     Select the manning's n for the slot:
c     nslot > 0.0 -- use nslot for manning's n in the slot
c     nslot = 0  -- use the mean value of the existing n at the 
c                   start and end points.  
c     nslot < 0.0 -- use the mean value of the existing n at the 
c                    start and end point multiplied by abs(nslot).

      if(nslot.gt.0.0) then
        n(nsub) = nslot
      else
        if(ematch.eq.'Y') then
         i = ie - 1
         if(i.lt.is) i = is
        else
          i = ie
        endif
        n(nsub) = 0.5*(lsn(is) + lsn(i))
        if(nslot.lt.0.0) then
          n(nsub) = n(nsub)*abs(nslot)
        endif
      endif

        IF(SMATCH.EQ.'Y') THEN
          SB(IS) = NSUB
          LSN(IS) = N(NSUB)
        ENDIF
        IF(EMATCH.EQ.'Y') THEN
          IF(SMATCH.EQ.'Y') THEN
            IOLD = IE
            INEW = IS + 2
            CALL SHIFT_VALUES
     I                         (STDOUT, IOLD, INEW, 
     M                          NPNT, X, Z, SB, LSN, SN)
            X(IS+1) = SLOT_XM
            Z(IS+1) = SLOT_ZM
            SB(IS+1) = NSUB
            SN(IS+1) = 1.0
            LSN(IS+1) = N(NSUB)
          ELSE
            IOLD = IE
            INEW = IS + 3
            CALL SHIFT_VALUES
     I                         (STDOUT, IOLD, INEW, 
     M                          NPNT, X, Z, SB, LSN, SN)
            X(IS+1) = XL
            Z(IS+1) = ZL
            SB(IS+1) = NSUB
            SN(IS+1) = 1.
            LSN(IS+1) = N(NSUB)

            X(IS+2) = SLOT_XM
            Z(IS+2) = SLOT_ZM
            SB(IS+2) = NSUB
            SN(IS+2) = 1.
            LSN(IS+2) = N(NSUB)
          ENDIF
        ELSE
          IF(SMATCH.EQ.'Y') THEN
            IOLD = IE + 1
            INEW = IS + 3
            CALL SHIFT_VALUES
     I                         (STDOUT, IOLD, INEW, 
     M                          NPNT, X, Z, SB, LSN, SN)
            X(IS+1) = SLOT_XM
            Z(IS+1) = SLOT_ZM
            SB(IS+1) = NSUB
            SN(IS+1) = 1.
            LSN(IS+1) = N(NSUB)

            X(IS+2) = XR
            Z(IS+2) = ZR
            SB(IS+2) = SAVED_SUB
            SN(IS+2) = 1.
            LSN(IS+2) = N(SAVED_SUB)
          ELSE
            IOLD = IE + 1
            INEW = IS + 4
            CALL SHIFT_VALUES
     I                         (STDOUT, IOLD, INEW, 
     M                          NPNT, X, Z, SB, LSN, SN)
            X(IS+1) = XL
            Z(IS+1) = ZL
            SB(IS+1) = NSUB
            SN(IS+1) = 1.
            LSN(IS+1) = N(NSUB)

            X(IS+2) = SLOT_XM
            Z(IS+2) = SLOT_ZM
            SB(IS+2) = NSUB                            
            SN(IS+2) = 1.
            LSN(IS+2) = N(NSUB)

            X(IS+3) = XR
            Z(IS+3) = ZR
            SB(IS+3) = SAVED_SUB
            SN(IS+3) = 1.
            LSN(IS+3) = N(SAVED_SUB)
          ENDIF
        ENDIF
      ENDIF

C     Adjust the subsection assignments
      CALL  REASUB
     I            (STDOUT, NPNT,
     M             NSUB, N, SB,
     O              EFLAG)

      RETURN
      END

C    
C
C     
      SUBROUTINE SET_SLOT(STDIN, STDOUT, STDTAB)

C     Set the cross section slot values for the triangular slot

      INTEGER STDIN, STDOUT, STDTAB

      INCLUDE 'abslot.cmn'

C     Local

      REAL SLOT
      CHARACTER LINE*80, NAME*5

C     *******************************FORMATS****************************
2     FORMAT(6X,F10.0)
4     FORMAT(6X,F10.0)
6     FORMAT(A5,1X,F10.0)

52    FORMAT(/,' Slot width for cross section invert=',F10.3)
54    FORMAT(/,' Slot n for cross section invert=',F10.3)
56    FORMAT(/,' Elevation of cross section slot invert=',F10.3)
58    FORMAT(/,' Depth of cross section slot invert=',F10.3)
60    FORMAT('; ',A5,F10.2)
C***********************************************************************
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) WSLOT
      WRITE(STDOUT,52) WSLOT
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) NSLOT
      WRITE(STDOUT,54) NSLOT
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) NAME, SLOT
      IF(NAME.EQ.'ESLOT') THEN
C       Set the fixed invert elevation.
        ESLOT = SLOT
C       Set fixed slot depth to invalid negative
C       value
        YSLOT = -1.0
      ELSE
C       Set the fixed slot depth
        YSLOT = SLOT
        ESLOT = 0.0
      ENDIF
      IF(YSLOT.LT.0.0) THEN
        WRITE(STDOUT,56) ESLOT
      ELSE
        WRITE(STDOUT,58) YSLOT
      ENDIF
      
      SLOT_PRESENT = 1

C     Output signal as a comment in the table file.
      IF(YSLOT.LT.0.0) THEN
        WRITE(STDTAB,60) NAME, ESLOT
      ELSE
        WRITE(STDTAB,60) NAME, YSLOT
      ENDIF

      RETURN

 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop' 

      END

C    
C
C     
      SUBROUTINE SET_SLOTE(STDIN, STDOUT, STDTAB)

C     Set the cross section slot values for the exponential slot 

      INTEGER STDIN, STDOUT, STDTAB

      INCLUDE 'abslot.cmn'

C     Local

      REAL SLOT
      CHARACTER LINE*80, NAME*5

C     *******************************FORMATS****************************
2     FORMAT(6X,F10.0)
4     FORMAT(6X,F10.0)
6     FORMAT(A5,1X,F10.0)

52    FORMAT(/,' Slot width for cross section invert=',F10.3)
54    FORMAT(/,' Slot n for cross section invert=',F10.3)
56    FORMAT(/,' Elevation of cross section slot invert=',F10.3)
58    FORMAT(/,' Depth of cross section slot invert=',F10.3)
60    FORMAT('; ',A5,F10.2,' Exponential slot.')
62    format(/,' Relative depth parameter=',f10.7)
64    format(/,' Top width of base section=',f10.6)
66    format(/,' Factor on exponent=',f10.6)
C***********************************************************************
c     Set default values
      rd = 0.37937619
      expfac = 1.0
c     Set  base top width to signal use of default
      tzero = -1.0

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) WSLOT
      WRITE(STDOUT,52) WSLOT


      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) NSLOT
      WRITE(STDOUT,54) NSLOT
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) NAME, SLOT
      IF(NAME.EQ.'ESLOT') THEN
C       Set the fixed invert elevation.
        ESLOT = SLOT
C       Set fixed slot depth to invalid negative
C       value
        YSLOT = -1.0
      ELSE
C       Set the fixed slot depth
        YSLOT = SLOT
        ESLOT = 0.0
      ENDIF
      IF(YSLOT.LT.0.0) THEN
        WRITE(STDOUT,56) ESLOT
      ELSE
        WRITE(STDOUT,58) YSLOT
      ENDIF
      
c     Get the optional parameter input values.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)

      if(line(1:6).eq.'RDEPTH') then
        read(line(8:),'(f10.0)') rd
        write(stdout,62) rd
      else
        backspace (stdin)
      endif
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)

      if(line(1:5).eq.'TZERO') then
        read(line(7:),'(f10.0)') tzero
        write(stdout,64) tzero
      else
        backspace (stdin)
      endif
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)

      if(line(1:6).eq.'EXPFAC') then
        read(line(8:),'(f10.0)') expfac
        write(stdout,66) expfac
      else
        backspace (stdin)
      endif

      SLOT_PRESENT = 2

C     Output signal as a comment in the table file.
      IF(YSLOT.LT.0.0) THEN
        WRITE(STDTAB,60) NAME, ESLOT
      ELSE
        WRITE(STDTAB,60) NAME, YSLOT
      ENDIF

      RETURN

 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop' 

      END

C
C
C
      SUBROUTINE CLEAR_SLOT(STDTAB)

C     Clear the slot values for cross section inverts.

      IMPLICIT NONE
      INTEGER STDTAB
      INCLUDE 'abslot.cmn'
C     *****************************FORMATS******************************
50    FORMAT('; CLRSLOT')
C***********************************************************************
      SLOT_PRESENT = 0
      WSLOT = 0.0
      NSLOT = 0.0
      ESLOT = -999999.
      YSLOT = -1.0
      WRITE(STDTAB,50)
      RETURN
      END

C    
C
C     
      SUBROUTINE   CXSTAB
     I                  (STDOUT, NSUBA, NAVM, NFAC, MAXPNT, LEFT, RIGHT,
     I                   BETOPT, SNFLG, NVAR, NATY, YATN, NNY,
     M                   NPNT, ZMIN, ZMAX, X, Z, SB, EFLAG, LSN, SN,
     O                   N, NDEP, XST)
 
C     + + + PURPOSE + + +
C     Compute cross section table from basic information.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, MAXPNT, NAVM, NDEP, NPNT, NSUBA, SNFLG, STDOUT
      INTEGER NNY(PMXSUB), NVAR(PMXSUB), SB(PMXPNT)
      REAL LEFT, LSN(PMXPNT), N(PMXSUB), NATY(9,PMXSUB), NFAC, RIGHT,
     A     SN(PMXPNT), X(PMXPNT), XST(MAXPNT,PMXELM), YATN(9,PMXSUB),
     B     Z(PMXPNT), ZMAX, ZMIN

      CHARACTER BETOPT*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NSUBA  - Number of subsections
C     NAVM   - Flag for averaging roughness
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     MAXPNT - Maximum number of tabulated values in a cross section
C               function table
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
C     BETOPT - Option for computing flux coefficients and critical flow
C     SNFLG  - Flag for sinuousity computations
C     NVAR   - Flag for variation of Manning's n in each subsection
C     NATY   - Mannings's n value at depth in YATN
C     YATN   - Depth values for the Manning's n values in NATY
C     NNY    - Number of values for Manning's n variation with depth
C     NPNT   - Number of points on boundary of a cross section
C     ZMIN   - Minimum elevation
C     ZMAX   - Maximum elevation
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     LSN    - Line segment Manning's n value
C     SN     - Sinuousity at a point on a cross section boundary
C     N      - Manning's n values
C     NDEP   - Number of depth values
C     XST    - Storage table for various elements of cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'nrdzcm.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER added, I, J, JJ, LIM, NPNTS, NSUB, WRN557
      INTEGER SBSUB(PMXPNT+2)
      REAL KOLD(PMXSUB), LSNS(PMXPNT+2), SNS(PMXPNT+2), SUBMIN,
     A     TSOLD(PMXSUB), XSUB(PMXPNT+2), XSV(PMXELM), ZSUB(PMXPNT+2),
     B     ZT(PMXPNT)
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKARG, COMPEL, RDUP, SORT, SUBSET, XCHK
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' Subset request being processed')
 52   FORMAT(/,
     A' *ERR:633* Only one point remains on the cross section.',
     B  /,10X,' Automatic extension required to compute.')
 54   format(/,'*PROBLEM* Subset operations on cross sections with',
     a  ' roughness varying with depth not yet supported.')
C***********************************************************************
      NSUB = NSUBA
C     Check for a subset operation before doing anything else: Jan.22, 2001
      IF(LEFT.LT.RIGHT) THEN
C       WE HAVE A SUBSET REQUEST
 
        if (nvar(1) == 0) then
          WRITE(STDOUT,50)
          NPNTS = NPNT + 2
          CALL SUBSET
     I             (STDOUT, NPNT, X, Z, SB, ZMAX, SN, LSN,
     M              LEFT, RIGHT, EFLAG, NPNTS, 
     O              XSUB, ZSUB, SBSUB, SNS, LSNS)

C         Change the number of subsections to force the walls added in
C         SUBSET to be frictionless. 
          SBSUB(1) = NSUB + 1
          SBSUB(NPNTS-1) = NSUB + 2
          NSUB = NSUB + 2

          DO J=1,NPNTS
            X(J) = XSUB(J)
            Z(J) = ZSUB(J)
            SB(J) = SBSUB(J)
            SN(J) = SNS(J)
            LSN(J) = LSNS(J)
          end do
          NPNT = NPNTS
        else
          write(stdout,54)
          stop 'Abnormal stop. Feature not yet supported.'
        endif
      ENDIF


c     Each subroutine checks data to determine which should work
      CALL ADD_SLOT
     I             (STDOUT, nvar(1),
     M              NPNT, NSUB, N, ZMIN, X, Z, SB, LSN, SN,
     o              added)

      CALL ADD_SLOTE
     I             (STDOUT, NFAC, nvar(1),
     M              NPNT, NSUB, N, ZMIN, X, Z, SB, LSN, SN,
     o              added)

c     Set all values of nvar to 0, to indicate constant Manning's n in each 
c     subsection.  This is the only option supported for slot insertion 
c     as of 31 August 2007.  If the computations reach here, then the 
c     variation of Manning's n IS constant, so we make sure that all 
c     values in nvar say that!
c      BUG FIX, commented out below line on 03.23.2009 per email from Delbert
c      Franz, documented in fixlist for version 5.93
c      nvar = 0
c      WRITE(STDOUT,*) ' '
c      WRITE(STDOUT,*) ' Check of adjusted cross section boundary.'
c      WRITE(STDOUT,*) ' NPNT=',NPNT
c      write(stdout,*) ' nsub=',nsub
c      WRITE(STDOUT,89)
c89    FORMAT(' Index     Offset  Elevation Subs       Lsn        Sn')
c      DO 1234 I=1,NPNT
c        WRITE(STDOUT,90) I, X(I), Z(I), SB(I), LSN(I), SN(I)
c90    FORMAT(' ',I5,2F11.4,I5,2F10.3)
c1234  CONTINUE

c      WRITE(STDOUT,*) ' '

C     CHECK FOR INVALID HORIZONTAL LINE SEGMENTS
 
      CALL XCHK
     I         (STDOUT, NPNT, ZMIN, X, NFAC,
     M          Z)
 

C     TRANSFER ELEVATIONS TO WORK SPACE, SORT,AND REMOVE DUPLICATES
      DO 100 J=1,NPNT
        ZT(J)=Z(J)
 100  CONTINUE
 
      CALL SORT
     I         (NPNT,
     M          ZT)
      CALL RDUP
     I         (NPNT,
     M          ZT,
     O          NDEP)
 
      IF(EXTEND.EQ.0) THEN
C       SUPPRESS AUTOMATIC EXTENSION OF LOWER END OF CROSS SECTION.
C       ZMAX, ON ENTRY, GIVES THE DESIRED MAXIMUM.
        
        DO 110 J=NDEP,1,-1
          IF(ZMAX.GT.ZT(J)) GOTO 111
 110    CONTINUE
        WRITE(STDOUT,52)
        EFLAG = 1
        RETURN
 
 111    CONTINUE
        NDEP = J + 1
      ENDIF

 
C     CHECK THE ARGUMENT SPACING FOR CONVEYANCE INTERPOLATION ACCURACY
 
      CALL CHKARG
     I           (STDOUT, PMXPNT, NRZERO, DZLIM, added,
     M            NDEP, ZT)
      ZMIN=ZT(1)
      ZMAX=ZT(NDEP)
 
 
      SUBMIN = ZMIN
C      WRITE(STDOUT,*) ' BEFORE CALL TO COMPEL NDEP=',NDEP
 
      JJ = 0
      DO 195 J=1,NSUB
        KOLD(J) = 0.0
        TSOLD(J) = 1.E30
 195  CONTINUE
      WRN557 = 1
      DO 200 J=1,NDEP
        IF(ZT(J).GE.SUBMIN) THEN
          JJ = JJ + 1

          CALL COMPEL
     I               (ZT(J), NPNT, NSUB, NAVM, X, Z, SB, NFAC, BETOPT,
     I                SNFLG, LSN, NVAR, NATY, YATN, NNY, SN, WRN557,
     M                KOLD, TSOLD,
     O                N, XSV)
 
          XSV(1)=ZT(J) - SUBMIN
          LIM = 17
C          IF(BETOPT(1:7).EQ.'NEWBETA') THEN
C            LIM = 13
C          ELSE
C            LIM = 7
C            XST(JJ,13) = XSV(13)
C          ENDIF
          DO 150 I=1,LIM
            XST(JJ,I)=XSV(I)
 150      CONTINUE
C          XST(J,14) = XSV(14)
C          XST(J,15) = XSV(15)
        ENDIF
 200  CONTINUE
      NDEP = JJ
      RETURN
      END
C
C
C
      SUBROUTINE   INFEQX
     I                   (STDIN, STDOUT, MAXPNT, MODE,
     M                    TABDIR, EFLAG,
     O                    TAB, STAT, NPNT, NSUB, NAVM, X, Z, SB, N,
     O                    LEFT, RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX,
     O                    zone, hgrid, vdatum, unitsys, basis)
 
C     + + + PURPOSE + + +
C     Input cross section in FEQX or FEQXLST format.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, MAXPNT, MODE, NAVM, NPNT, NSUB, STDIN, STDOUT, TAB
      INTEGER SB(*), TABDIR(*)
      REAL LEFT, N(*), RIGHT, STAT, X(*), Z(*), ZMAX
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8, zone*8, hgrid*8, vdatum*8,
     a          unitsys*8, basis*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     MAXPNT - Maximum number of tabulated values in a cross section
C               function table
C     MODE   - Mode of processing cross section boundary specification:
C               MODE=1: fixed format.  MODE=2: list format.
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TAB    - Table number
C     STAT   - Station value
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     NAVM   - Flag for averaging roughness
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     N      - Manning's n values
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
C     SAVOPT - Function table saving option
C     OUTOPT - Output option for the table file for cross section function
C               tables
C     BETOPT - Option for computing flux coefficients and critical flow
C     ZMAX   - Maximum elevation

C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xtadd.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IE, IS, J, NTEMP
      REAL HSHIFT, SCALE, SHIFT, XOLD, VSCALE
      CHARACTER CHR10*10, HEAD*80, LINE*80, MONTON*8,
     A          EXTEND_STRING*8, VARN*4
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      EXTERNAL GETNXT, inline, REASUB, SETOPT, TABCHK, LENSTR
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,A)
 2    FORMAT(8X,F10.0,6X,F10.0,7X,F10.0)
 3    FORMAT(A4,I5,6F10.0)
 4    FORMAT(A80)
 5    FORMAT(2F10.0,I5)
 6    FORMAT(5X,I5,7X,F10.0,7X,F10.0,8X,F10.0,8X,F10.0)
 31   FORMAT(9X,6F10.0)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(/,' TabId= ',A,5(1X,A8))
 52   FORMAT(' NAVM=',I5,'  Hor. scale factor=',F10.3,
     A ' Vert. shift=',F10.3,/,'  Vert. scale factor=',F10.3,
     B ' Hor. shift=',F10.3)
 53   FORMAT(' STATION=',F10.3,' LEFT=',F10.1,' RIGHT=',F10.1)
 54   FORMAT(' STATION=',F10.3)
 55   FORMAT(1X,A4,I5,6F6.3,1X,/,(10X,6F6.3))
 56   FORMAT(1X,A80)
 57   FORMAT(1X,F10.1,F10.2,I5)
 58   FORMAT(' *ERR:504* NUMBER OF SUBSECTIONS=',I5,' > ',I5)
 59   FORMAT(' *ERR:505* SUBSECTION NUMBER TOO LARGE AT OFFSET=',
     1       F10.1)
 60   FORMAT(' *WRN:501* SUBSECTION VALUE FOR FIRST POINT IS',
     1 ' MISSING.  SUBSECTION=1 ASSUMED.')
 61   FORMAT(' *ERR:506* ONLY ONE POINT GIVEN ON BOUNDARY OF THE',
     1  ' CROSS SECTION.')
 62   FORMAT(' *ERR:507* NUMBER OF POINTS IN CROSS SECTION >',I5)
 65   format(' ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a        ' BASIS=',a8)
 66   FORMAT(' Selection of beta option "NEWBETA" implies checking for',
     A      ' monotonicity.')
 67   FORMAT(/,' *WRN:554* Extending left end of cross section',
     A   ' by ',F8.3)
 68   FORMAT(/,' *WRN:555* Extending right end of cross section',
     A   ' by ',F8.3)
 69   FORMAT(/,' *WRN:556* Some point in cross section higher than',
     A  ' either end.',/,10X,'All area above minimum end elevation',
     b  ' is ignored.')
 70   FORMAT(/,' *ERR:647* Manning''s N = 0 in subsection',I5,
     A /,11X,'Value set to 1.0.  Please correct and recompute.')
72    FORMAT(' GISID=',A16,' EASTING=',F15.3,' NORTHING=',F15.3)
73    FORMAT('  Processing FEQX TabId= ',A) 
C***********************************************************************
 
      CALL GET_XSEC_HEADER(STDIN, STDOUT,
     M   EFLAG,
     O   TAB, SAVOPT, OUTOPT, MONTON, BETOPT, EXTEND, 
     O   STAT, LEFT, RIGHT, VARN, NAVM, SCALE, SHIFT,
     O   VSCALE, HSHIFT, NSUB, N, zone, hgrid, vdatum, unitsys, basis)
      IF(EXTEND.EQ.0) THEN
        EXTEND_STRING = 'NOEXTEND'
      ELSE
        EXTEND_STRING = 'EXTEND'
      ENDIF
      WRITE(*,73) TABID(1:LENSTR(TABID))
      WRITE(STDOUT,51) TABID(1:LENSTR(TABID)), SAVOPT, OUTOPT, BETOPT,
     A                 EXTEND_STRING, MONTON
      IF(TAB.LT.0) THEN
        SLOT = 0.03
        TAB = -TAB
        NOCM = 1
      ELSE
         NOCM = 0
         SLOT = 0.0
      ENDIF

      if(zone /= 'NONE') then
        write(stdout,65) zone, hgrid, vdatum, unitsys, basis
      endif
      IF(TAB.GE.0) CALL TABCHK
     I                        (STDOUT, PMXTAB,
     M                         TAB, TABDIR, EFLAG)

      LEFT = SCALE*LEFT
      RIGHT = SCALE*RIGHT
      IF(GISID.NE.' ') THEN
        WRITE(STDOUT,72) GISID, EASTING, NORTHING
      ENDIF
      IF(LEFT.GE.RIGHT) THEN
        WRITE(STDOUT,54) STAT
      ELSE
        WRITE(STDOUT,53) STAT, LEFT, RIGHT
      ENDIF
      WRITE(STDOUT,52) NAVM, SCALE, SHIFT, VSCALE, HSHIFT
      WRITE(STDOUT,55) 'NSUB', NSUB, (N(J),J=1,NSUB)
 
      DO 131 J=1,NSUB
        IF(N(J).LE.0.0) THEN
          WRITE(STDOUT,70) J
          N(J) = 1.0
          EFLAG = 1
        ENDIF
 131  CONTINUE
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) HEAD
      WRITE(STDOUT,56) HEAD
 
C     GET THE CO-ORDINATES ON THE BOUNDARY OF THE CROSS SECTION
 
      XOLD = -1.E30
      NPNT=0
      ZMAX = -1.E30
 135  CONTINUE
        NPNT=NPNT+1
        IF(NPNT.LE.MAXPNT) GOTO 137
          WRITE(STDOUT,62) MAXPNT
          EFLAG = EFLAG + 1
          NPNT=MAXPNT
 137    CONTINUE
 
        IF(MODE.EQ.2) THEN
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
C         EMULATE LIST PROCESSING OF ITEMS IN LINE.
 
          IS = 1
          CALL GETNXT
     I               (LINE, IS,
     O                IE, CHR10)
          READ(CHR10,'(F10.0)',ERR=991) X(NPNT)
 
          IS = IE
          CALL GETNXT
     I               (LINE, IS,
     O                IE, CHR10)
          READ(CHR10,'(F10.0)',ERR=991) Z(NPNT)
 
          IS = IE
          CALL GETNXT
     I               (LINE, IS,
     O                IE, CHR10)
          READ(CHR10,'(I10)', ERR=991) SB(NPNT)
 
C          READ(LINE,*,ERR=991) X(NPNT),Z(NPNT),SB(NPNT)
C          WRITE(STDOUT,*) X(NPNT), Z(NPNT), SB(NPNT)
 
        ELSE IF(MODE.EQ.1) THEN
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
          READ(LINE,5,ERR=991) X(NPNT), Z(NPNT),SB(NPNT)
        ENDIF
        X(NPNT) = X(NPNT)*SCALE + HSHIFT
        Z(NPNT) = Z(NPNT)*VSCALE + SHIFT
 
C       FIND MAXIMUM ELEVATION IN CROSS SECTION FOR LATER CHECKING
        ZMAX = MAX(Z(NPNT), ZMAX)
 
C       CHECK FOR MONOTONICITY OF TOP WIDTH.  THIS REQUIRES THAT
C       THE OFFSET NEVER DECREASE.
 
        IF(MONTON.EQ.'MONOTONE') THEN
          IF(X(NPNT).LT.XOLD) THEN
            WRITE(STDOUT,*)
     A  ' *ERR:508* SECTION VIOLATES MONOTONICITY AT OFFSET=', X(NPNT)
            EFLAG = EFLAG + 1
          ENDIF
        ENDIF
 
        XOLD = X(NPNT)
 
        IF(NPNT.GT.1) GO TO 160
          IF(SB(1).GT.0) GOTO 150
            WRITE(STDOUT,60)
            SB(1)=1
            EFLAG = 1
 150      CONTINUE
 160    CONTINUE
C       PROPAGATE SUBSECTION ASSIGNMENTS INTO ZERO VALUES
 
        IF(SB(NPNT).EQ.0) SB(NPNT)=SB(NPNT-1)
 
        WRITE(STDOUT,57) X(NPNT),Z(NPNT),SB(NPNT)
 
        IF(SB(NPNT).LE.NSUB) GOTO 140
          WRITE(STDOUT,59) X(NPNT)
          SB(NPNT) = NSUB
          EFLAG = EFLAG + 1
 140    CONTINUE
 
 
 
        IF(SB(NPNT).GT.0) GOTO 135
 
C     NOTE THAT A NEGATIVE VALUE OF SB(NPNT) TERMINATES INPUT
 
      IF(NPNT.GT.1) GOTO 170
        WRITE(STDOUT,61)
        EFLAG = 1
 170  CONTINUE
 
C     CHECK THE CROSS SECTION FOR NONSENSE BEHAVIOR AT THE END
 
      IF(Z(1).LE.Z(2).AND.X(2).GT.X(1)) THEN
C       THE LEFT MOST LINE SEGMENT HAS UPWARD SLOPE, THEREFORE HIGH
C       POINT IS NOT AT THE LEFT LIMIT.
 
        WRITE(STDOUT,'(/,A,A)') ' *WRN:502* UNEXPECTED SLOPE',
     A   ' AT LEFT END.'
        WRITE(STDOUT,*) '  SLOPE EXPECTED TO BE < 0  AT LEFT BOUNDARY'
      ENDIF
      IF(Z(NPNT).LE.Z(NPNT-1).AND.X(NPNT).GT.X(NPNT-1)) THEN
C       THE RIGHT MOST LINE SEGMENT HAS DOWNWARD SLOPE
 
        WRITE(STDOUT,'(/,A,A)') ' *WRN:503* UNEXPECTED SLOPE',
     A   ' AT RIGHT END.'
        WRITE(STDOUT,*) '  SLOPE EXPECTED TO BE > 0 AT RIGHT BOUNDARY'
      ENDIF
 
      IF(EXTEND.EQ.1) THEN
C       CHECK FOR ONE END BEING HIGHER THAN THE OTHER
 
        IF(ABS(ZMAX - Z(1)).GT.EPSDIF) THEN
          WRITE(STDOUT,67) ZMAX - Z(1)
        ENDIF
        IF(ABS(ZMAX - Z(NPNT)).GT.EPSDIF) THEN
          WRITE(STDOUT,68) ZMAX - Z(NPNT)
        ENDIF
      ELSE
C       CHECK FOR AN INTERMEDIATE POINT BEING HIGHER THAN EITHER
C       END
        IF(ABS(ZMAX - Z(1)).GT.EPSDIF.AND.
     A     ABS(ZMAX - Z(NPNT)).GT.EPSDIF) THEN
          WRITE(STDOUT,69)
        ENDIF
 
        ZMAX = MIN(Z(1), Z(NPNT))
      ENDIF
 
C     CHECK TO SEE IF SUBSECTION NUMBERS HAVE BEEN REPEATED IN
C     SEPERATE RUNS.
 
      NTEMP = NSUB
      CALL REASUB
     I           (STDOUT, NPNT,
     M            NSUB, N, SB,
     O            EFLAG)
      IF(NSUB.NE.NTEMP) THEN
C       OUTPUT THE NEW ASSIGNMENTS
        WRITE(STDOUT,*) ' '
        WRITE(STDOUT,*) ' *WRN:548* Unwise use of subsection numbers'
        WRITE(STDOUT,*) '     Subsections have been added to avoid',
     A                  ' repeated usages.'
        WRITE(STDOUT,*) '     Old NSUB=',NTEMP,' New NSUB=',NSUB,
     A                  ' Please check for validity'
        WRITE(STDOUT,55) 'NSUB', NSUB, (N(J),J=1,NSUB)
        WRITE(STDOUT,56) HEAD
        DO 200 J=1,NPNT
          WRITE(STDOUT,57) X(J),Z(J),SB(J)
 200    CONTINUE
      ENDIF

 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   CHANEL
     I                   (STDIN, STDOUT, STDTAB, NFAC,
     M                    TABDIR, FTP,
     O                    EFLAG, MODE)
 
C     + + + PURPOSE + + +
C     Compute cross sections including the sinuosity elements.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, MODE, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(PMXTAB)
      REAL NFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     TABDIR - Table directory to remember table numbers
C     FTP    - next open location in the function table storage
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     MODE   - Mode of processing cross section boundary specification:
C               MODE=1: fixed format.  MODE=2: list format.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xscomu.cmn'
      INCLUDE 'sincom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER EFLG
      CHARACTER LINE*80, NXTCMD*8
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FEQX, FQXE, inline, STBIN
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:679* Invalid cross section command:',A8,' in ',
     A       'the CHANNEL command.')
 52   FORMAT(1X,A)
 54   FORMAT(/,' *ERR:680* CHANNEL command ended improperly.  FINISH',
     A      ' encountered before ENDCHAN.')
 56   FORMAT(/,' Errors prevent CHANNEL command completion.')
 58   FORMAT(/,' End of file found.  Missing ENDCHAN or FINISH.')
C***********************************************************************
C     DEFINE THE SINUOSITY VALUES FROM THE USER SPECIFICATION.
 
      CALL STBIN
     I          (STDIN, STDOUT,
     O           STL, OFFSET, NUMOFF, SINU, VARTYP, NUMSEC, DIR, EPS,
     O           JAXIS, EFLAG)
 
      IF(EFLAG.NE.0) THEN
C       FLUSH THE REST OF THE CHANNEL INPUT TO AVOID FURTHER ERRORS.
        WRITE(STDOUT,56)
 100    CONTINUE
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
          WRITE(STDOUT,52) LINE
          READ(LINE,'(A8)') NXTCMD
          IF(NXTCMD.EQ.'FINISH') THEN
            WRITE(STDOUT,54)
            RETURN
          ELSEIF(NXTCMD.EQ.'ENDFILE') THEN
            WRITE(STDOUT,58)
            RETURN
          ENDIF
          IF(NXTCMD.NE.'ENDCHAN') GOTO 100
          RETURN
      ENDIF
 
C     READ THE CROSS SECTION INPUT AND PROCESS.  ANY ONE OF THREE
C     COMMANDS IS POSSIBLE: FEQX, FEQXLST, OR FEQXEXT.
 
 200  CONTINUE
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        WRITE(STDOUT,52) LINE
        READ(LINE,'(A8)') NXTCMD
 
        IF(NXTCMD.EQ.'FEQX    ') THEN
C         PROCESS FEQX FORMAT CROSS SECTION WITH SINUOSITY.
 
          MODE= 2
          SNFLGU = VARTYP
          CALL FEQX
     I             (STDIN, STDOUT, STDTAB, NFAC, MODE,
     M              TABDIR, FTP,
     O              EFLG)
 
          EFLAG = EFLAG + EFLG

        ELSEIF(NXTCMD.EQ.'FEQXLST ') THEN
C         PROCESS FEQX FORMAT IN LIST FORM
 
          MODE = 1
          SNFLGU = VARTYP
          CALL FEQX
     I             (STDIN, STDOUT, STDTAB, NFAC, MODE,
     M              TABDIR, FTP,
     O              EFLG)
          EFLAG = EFLAG + EFLG
        ELSEIF(NXTCMD.EQ.'FEQXEXT') THEN
C         PROCESS THE EXTENDED FEQX FORMAT
 
          SNFLGU = VARTYP
          CALL FQXE
     I             (STDIN, STDOUT, STDTAB, NFAC,
     M              TABDIR, FTP,
     O              EFLG)
          EFLAG = EFLAG + EFLG
 
        ELSEIF(NXTCMD.EQ.'ENDCHAN') THEN
C         END THE CHANNEL COMMAND
          RETURN
        ELSEIF(NXTCMD.EQ.'        ') THEN
C         OK. SKIP BLANK LINES
 
        ELSEIF(NXTCMD.EQ.'FINISH  ') THEN
          WRITE(STDOUT,54)
          STOP 'Abnormal stop. Errors found.'
        ELSEIF(NXTCMD.EQ.'SETSLOT') THEN
            CALL SET_SLOT(STDIN, STDOUT, STDTAB)
        ELSEIF(NXTCMD.EQ.'SETSLOTE') THEN
            CALL SET_SLOTE(STDIN, STDOUT, STDTAB)
        ELSEIF(NXTCMD.EQ.'CLRSLOT') THEN
            CALL CLEAR_SLOT(STDTAB)
        ELSE
          WRITE(STDOUT,50) NXTCMD
          EFLAG = 1
        ENDIF
 
        GOTO 200
      END
C
C
C
      SUBROUTINE   FEQX
     I                 (STDIN, STDOUT, STDTAB, NFAC, MODE,
     M                  TABDIR, FTP,
     O                  EFLAG)
 
C     + + + PURPOSE + + +
C     Process a cross section in FEQ format.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, MODE, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL NFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     MODE   - Mode of processing cross section boundary specification:
C               MODE=1: fixed format.  MODE=2: list format.
C     TABDIR - Table directory to remember table numbers
C     FTP    - next open location in the function table storage
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'fldway.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'sincom.cmn'

C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, LOC, MESG, NPI
      REAL LEFT, RIGHT, XARG(PMXNFL), ZMAX, area
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c  zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8

 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LOCSTA
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSTAB, FNDWAY, INFEQX, INSPT, LKUPSN, LOCSTA, TABOUT,
     A         TRAN
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* NAVM > 0 incompatible with curvilinear',
     A       ' element computation.')
C***********************************************************************
C     ENABLE WARNINGS ABOUT SUBSECTION CONVEYANCE VARIATION
      NOCM = 0
      SLOT = 1.E30
      EFLAG = 0
      IF(NPNTU.EQ.0) GOTO 100
        CALL TRAN
 100  CONTINUE
 
C     INPUT NEW CROSS SECTION
 
      IF(MODE.NE.3) THEN
        CALL INFEQX
     I             (STDIN, STDOUT, MXPNTU, MODE,
     M              TABDIR, EFLAG,
     O              TABU, STATU, NPNTU, NSUBU, NAVMU, XU, ZU, SBU, NU,
     O              LEFT, RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX,
     O              zone, hgrid, vdatum, unitsys, basis)
      ELSE
        WRITE(STDOUT,*) ' *ERR:558* DIGITIZED INPUT NOT SUPPORTED'
        STOP 'Abnormal stop. Errors found.'
 
      ENDIF

C     CLEAR AND SET ITEMS USED IN THE COMPUTATION BUT NOT DEFINED BY
C     FEQX INPUT
 
      DO 95 I=1,NSUBU
         NVARU(I) = 0
         NNYU(I) = 0
 95   CONTINUE
 
C     SET DEFAULT FOR  SINUOSITY
      DO 96 I=1,NPNTU
        SNU(I) = 1.0
 96   CONTINUE

      IF(SNFLGU.GT.0) THEN
        IF(NAVMU.NE.0) THEN
          WRITE(STDOUT,50)
          EFLAG = 1
        ENDIF
 
C       LOCATE THE CURRENT STATION IN THE SINUOUSITY TABLE.
        LOC = LOCSTA(STDOUT, STATU, DIR, JAXIS, NUMSEC,
     B                  STL, EPS)
 
        IF(LOC.LE.0) THEN
          EFLAG = 1
        ENDIF
        IF(EFLAG.NE.0) RETURN
 
C       INSERT THE OFFSETS FOR SINUOUSITY DEFINITION INTO THOSE
C       ALREADY DEFINING THE BOUNDARY.
        NPI = NUMOFF(LOC)
 
        DO 98 I=1,NPI
          XARG(I) = OFFSET(LOC,I)
 98     CONTINUE
        CALL INSPT
     I            (STDOUT, NPI, XARG, VARTYP,
     M             NPNTU, XU, ZU, SBU, LSNU,
     O             EFLAG)
 
        IF(EFLAG.NE.0) RETURN
 
C       FIND THE SINUOSITY FOR EACH POINT ON THE CROSS SECTION
C       BOUNDARY
 
        CALL LKUPSN
     I             (LOC, NPNTU, XU, VARTYP, SINU, OFFSET, NPI,
     O              SNU)
 
C        WRITE(STDOUT,*) ' '
C        WRITE(STDOUT,*) ' CHECK OF SINUOSITY ASSIGNMENTS'
C        WRITE(STDOUT,1234)
C1234  FORMAT(1X,4X,'OFFSET ELEVATION SUBS SINUOSITY')
C        DO 1235 I=1,NPNTU
C          WRITE(STDOUT,1236) XU(I), ZU(I), SBU(I), SNU(I)
C1236  FORMAT(1X,F10.2,F10.2,I5,F10.6)
C1235    CONTINUE
      ENDIF
 
      DO 97 I=1,NPNTU-1
        LSNU(I) = NU(SBU(I))
C        WRITE(STDOUT,*) 'I=',I,' LSNU=',LSNU(I)
 97   CONTINUE

C     COMPUTE ELEMENTS FOR CURRENT CROSS SECTION
 
      IF(NPNTU.GT.1) THEN
C       FIND THE MAXIMUM AND MINIMUM ARGUMENT VALUES
        ZMINU= 9999999.
        ZMAXU = -9999999.
        DO 150 J=1,NPNTU
          IF(ZU(J).GT.ZMAXU) ZMAXU = ZU(J)
          IF(ZU(J).LT.ZMINU) ZMINU = ZU(J)
 150    CONTINUE
 
        IF(EXTEND.EQ.0) ZMAXU = ZMAX
 
C       FIND THE VALUE OF RIGHT AND LEFT IF FLOOD WAY OPTION
C       IS ACTIVE
 
        IF(FLOOD.EQ.1) THEN
          CALL FNDWAY
     I               (TABU, STDOUT, EFLAG, NFAC,
     O                LEFT, RIGHT, area)
        ENDIF

        IF(EFLAG.EQ.0) THEN
          CALL CXSTAB
     I               (STDOUT, NSUBU, NAVMU, NFAC, MXPNTU, LEFT, RIGHT,
     I                BETOPT, SNFLGU, NVARU, NATYU, YATNU, NNYU,
     M                NPNTU, ZMINU, ZMAXU, XU, ZU, SBU, EFLAG, LSNU,
     M                SNU,
     O                NU, NDEPU, XSTU)
 
          IF(FLOOD.EQ.1) THEN
C           STORE LEFT AND RIGHT AGAIN TO
C           REFLECT CHANGES MADE IN CSXTAB
            FLDLT(TABU) = LEFT
            FLDRT(TABU) = RIGHT
            fldarea(tabu) = area
          ENDIF
        ENDIF
      ENDIF
 
c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

      call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

C     OUTPUT THE TABLE IF NO ERRORS AND IF OUTPUT IS REQUESTED
 
      IF(EFLAG.NE.0.OR.TABU.EQ.0) GOTO 200
        IF(NOCM.EQ.1) THEN
          MESG = 0
        ELSE
          MESG = 1
        ENDIF
        CALL TABOUT
     I             (STDOUT, STDTAB, TABU, STATU, ZMINU, MESG, SAVOPT,
     I              OUTOPT, BETOPT, zone, hgrid, vdatum, unitsys, basis,
     i                    khflag, alphaflag, betaflag, maflag, mqflag,
     M              NDEPU, XSTU, FTP)
 200  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   XSTMAK
     I                   (STDIN, STDOUT, STDTAB,
     M                    EFLAG, FTP, FTKNT)
 
C     + + + PURPOSE + + +
C     Input a branch like structure and output every interpolated
C     cross section table requested by the user.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTKNT, FTP, STDIN, STDOUT, STDTAB
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
C     FTKNT  - function table counter
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'culcom.cmn'
      INCLUDE 'ftable.cmn'
      INCLUDE 'flowline.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER NBRA
      SAVE NBRA
 
C     + + + LOCAL VARIABLES + + +
      INTEGER BEGTAB, BSHAPE, CHKBAR, ENDTAB, I, IE, IS, J, NBN, NG,
     A        TYPFLG
      INTEGER BRPT(8,1), NEGTAB(MNBN)
      REAL SFAC
      CHARACTER BNODID(MNBN)*8, LINE*80, NODEID*4, FLNAME*6, FLFILE*64
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CLVIN, GETTBN, inline, IXTOUT, os_file_style
 
C     + + + DATA INITIALIZATIONS + + +
      DATA NBRA/1/
 
C     + + + INPUT FORMATS + + +
 2    FORMAT(5X,F10.0)
 4    FORMAT(7X,A4)
 6    FORMAT(7X,A)
 8    FORMAT(7X,A)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' Interpolation of cross section tables not done.',
     a       '  Errors found.')
 52   FORMAT(/,' STATION FACTOR=',F10.2)
 54   FORMAT(/,' NODE ID PRESENT: ',A4)
 56   FORMAT(/,' Flow-line name= ',A,' will be sought in file name: ',
     A  /,5X,A)
C***********************************************************************
C     Clear the flag for flow-line data
      FL_PRESENT = 0
C     Reset the initial index for searching
      L = 1
C     Get the station factor
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) SFAC
      WRITE(STDOUT,52) SFAC
 
C     Get the nodeid option string
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) NODEID
      WRITE(STDOUT,54) NODEID
 
C     Check for the optional input of a river-mile line that 
C     will define the location in plan of the cross-section. 
C     Used to establish the intersection between the flow line
C     and the interpolated cross section.  
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IF(LINE(1:6).NE.'FLNAME')  THEN
C       Put the line back!  Optional input is missing!
        BACKSPACE(STDIN)
      ELSE
C       Optional input is present 
        READ(LINE,6) FLNAME
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,8) FLFILE
        CALL MAYBE_ADD_HOME(
     M                       flfile)
        call os_file_style( 
     m                              flfile)

        WRITE(STDOUT,56) FLNAME, FLFILE
        CALL PROCESS_FLNAME(STDOUT, FLNAME, FLFILE)
      ENDIF

        
      WRITE(STDOUT,*) ' '
C     Set the type flag to a non-culvert case.  TYPFLG is 1 
C     if a culvert barrel is involved.  Note: CHKBAR is set 
C     here even though it is an output value because if TYPFLG is
C     0, CLVIN does not set CHKBAR.  CHKBAR is only used in the
C     CULVERT code.   NEGTAB is only used here and not in the 
C     CULVERT code. We use HLTAB here because it is never used
C     in XSINTERP.  It supports a special means of giving the location
C     of an interpolated cross section that is needed to support the 
C     processing of input for FEQUTL created by the stand alone
C     utility OVERFLOW.  HLTAB is used as work space in this case
C     and does not return useful information here.

      TYPFLG = 0
      CHKBAR = 0
      CALL CLVIN
     I          (SFAC, STDIN, STDOUT, MNBN, NBRA, NODEID, MFTNUM,
     I           MRFTAB, FTPNT, TYPFLG,
     M           EFLAG, FTKNT, FTP,
     O           NBN, BRPT, NSEC, XVEC, ZBVEC, KA, KD, HLTAB, BNODID,
     O           NEGTAB, CHKBAR, IAT3D, IAT6D, BSHAPE, sbkind)
 
      IF(EFLAG.EQ.0) THEN
C       ON RETURN THE INTERPOLATION REQUESTS HAVE BEEN DONE AND
C       THE TABLES ARE STORED IN FTAB() AND NODES AT WHICH
C       INTERPOLATION TABLES BEEN CREATED HAVE A NON-ZERO ENTRY IN
C       NEGTAB.  HOWEVER, NEGTAB DOES NOT  GIVE THE CURRENT
C       TABLE NUMBER.  NSEC() CONTAINS THE TABLE ADDRESS FOR THE
C       TABLE.
 
C       SCAN NEGTAB() FINDING THE REAL TABLES THAT BRACKET ANY
C       INTERPOLATED TABLES, AND OUTPUT THE INTERPOLATED TABLES
C       TO THE STANDARD TABLE FILE.
 
C       CLEAR THE NEGATIVE TABLE FLAG
        NG = 0
 
        I = 1
 100    CONTINUE
          IF(NG.EQ.0) THEN
C           WE HAVE NOT SEEN A NEGATIVE TABLE NUMBER YET
            IF(NEGTAB(I).EQ.0) THEN
C             REMEMBER THIS LOCATION IN CASE IT IS THE INITIAL
C             POINT FOR INTERPOLATION
              IS = I
            ELSE
C             WE HAVE A NEGATIVE TABLE NUMBER. SET THE FLAG
              NG = 1
            ENDIF
          ELSE
C           WE HAVE SEEN AT LEAST ONE NEGATIVE TABLE NUMBER
            IF(NEGTAB(I).EQ.0) THEN
C             FOUND THE FINAL POINT FOR THE INTERPOLATION
              IE = I
 
C             OUTPUT THE TABLES FROM IS+1 TO IE-1
 
              BEGTAB = GETTBN(NSEC(IS))
              ENDTAB = GETTBN(NSEC(IE))
 
C             USE THE BEGINNING AND ENDING TABLE NUMBERS AS A
C             LABEL ON THE TABLE
 
              DO 200 J=IS+1,IE-1
                CALL IXTOUT
     I                     (STDTAB, STDOUT, NSEC(J), BEGTAB, ENDTAB)
 200          CONTINUE
C             CLEAR THE NEGATIVE TABLE FLAG
              NG = 0
              IS = IE
            ENDIF
          ENDIF
 
          I = I + 1
          IF(I.LE.NBN) GOTO 100
 
C       ALL TABLES HAVE BEEN OUTPUT
      ELSE
        WRITE(STDOUT,50)
      ENDIF
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
 
      END
C
C
C
      SUBROUTINE   FQXE
     I                 (STDIN, STDOUT, STDTAB, NFAC,
     M                  TABDIR, FTP,
     O                  EFLAG)
 
C     + + + PURPOSE + + +
C     Process a cross section in FEQ extended format
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL NFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     TABDIR - Table directory to remember table numbers
C     FTP    - next open location in the function table storage
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'fldway.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'sincom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, LOC, MESG, NPI
      REAL LEFT, RIGHT, XARG(PMXNFL), ZMAX, area
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8

 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LOCSTA
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSTAB, FNDWAY, INFQXE, INSPT, LKUPSN, LOCSTA, TABOUT
C***********************************************************************
C     ENABLE WARNINGS ABOUT SUBSECTION CONVEYANCE VARIATION
      NOCM = 0
      SLOT = 1.E30
      EFLAG = 0
 
C     INPUT NEW CROSS SECTION
 
      CALL INFQXE
     I           (STDIN, STDOUT,
     M            TABDIR, EFLAG,
     O            TABU, STATU, NPNTU, NSUBU, XU, ZU, SBU, NU, LEFT,
     O            RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX, LSNU, NVARU,
     O            NATYU, YATNU, NNYU, zone, hgrid, vdatum, unitsys, 
     O            basis)
 
C     SET DEFAULT FOR SINUOUSITY VALUES
      DO 96 I=1,NPNTU
        SNU(I) = 1.0
 96   CONTINUE
      IF(SNFLGU.GT.0) THEN
C       LOCATE THE CURRENT STATION IN THE SINUOUSITY TABLE.
        LOC =   LOCSTA(STDOUT, STATU, DIR, JAXIS, NUMSEC,
     B                  STL, EPS)
        IF(LOC.LE.0) THEN
          EFLAG = 1
        ENDIF
 
        IF(EFLAG.NE.0) RETURN
 
C       INSERT THE OFFSETS FOR SINUOUSITY DEFINITION INTO THOSE
C       ALREADY DEFINING THE BOUNDARY.
        NPI = NUMOFF(LOC)
        DO 98 I=1,NPI
          XARG(I) = OFFSET(LOC,I)
 98     CONTINUE
        CALL INSPT
     I            (STDOUT, NPI, XARG, VARTYP,
     M             NPNTU, XU, ZU, SBU, LSNU,
     O             EFLAG)
 
 
        IF(EFLAG.NE.0) RETURN
 
C       FIND THE SINUOSITY FOR EACH POINT ON THE CROSS SECTION
C       BOUNDARY
C       ELEMENTS
          CALL LKUPSN
     I               (LOC, NPNTU, XU, VARTYP, SINU, OFFSET, NPI,
     O                SNU)
 
C        WRITE(STDOUT,*) ' '
C        WRITE(STDOUT,*) ' CHECK OF SINUOSITY ASSIGNMENTS'
C        WRITE(STDOUT,1234)
C1234  FORMAT(1X,4X,'OFFSET ELEVATION SUBS SINUOSITY')
C        DO 1235 I=1,NPNTU
C          WRITE(STDOUT,1236) XU(I), ZU(I), SBU(I), SNU(I)
C1236  FORMAT(1X,F10.2,F10.2,I5,F10.6)
C1235    CONTINUE
      ENDIF
 
 
C     COMPUTE ELEMENTS FOR CURRENT CROSS SECTION
 
      IF(NPNTU.GT.1) THEN
C       FIND THE MAXIMUM AND MINIMUM ARGUMENT VALUES
        ZMINU= 9999999.
        ZMAXU = -9999999.
        DO 150 J=1,NPNTU
          IF(ZU(J).GT.ZMAXU) ZMAXU = ZU(J)
          IF(ZU(J).LT.ZMINU) ZMINU = ZU(J)
 150    CONTINUE
 
        IF(EXTEND.EQ.0) ZMAXU = ZMAX
 
C       FIND THE VALUE OF RIGHT AND LEFT IF FLOOD WAY OPTION
C       IS ACTIVE
 
        IF(FLOOD.EQ.1) THEN
          CALL FNDWAY
     I               (TABU, STDOUT, EFLAG, NFAC,
     O                LEFT, RIGHT, area)
        ENDIF
        IF(EFLAG.EQ.0) THEN
          navmu = 0    ! This value overridden by contents of navru in this case
          CALL CXSTAB
     I               (STDOUT, NSUBU, NAVMU, NFAC, MXPNTU, LEFT, RIGHT,
     I                BETOPT, SNFLGU, NVARU, NATYU, YATNU, NNYU,
     M                NPNTU, ZMINU, ZMAXU, XU, ZU, SBU, EFLAG, LSNU,
     M                SNU,
     O                NU, NDEPU, XSTU)
 
 
          IF(FLOOD.EQ.1) THEN
C           STORE LEFT AND RIGHT AGAIN TO
C           REFLECT CHANGES MADE IN CSXTAB
            FLDLT(TABU) = LEFT
            FLDRT(TABU) = RIGHT
            fldarea(tabu) = area
          ENDIF
        ENDIF
      ENDIF
 

c     Check if an exponential slot is present and  optionally replace alpha, 
c     beta, ma, and mq from a function table supplied by the user.  The
c     id for this table is the current table id with a lower case a appended
c     to it.  If such a table is not present, then the element values remain
c     as computed. 
      call replace_eslot_elements(stdout, tabu, ndepu,
     m                              xstu)
      
c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

      call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

C     OUTPUT THE TABLE IF NO ERRORS AND IF OUTPUT IS REQUESTED
 
      IF(EFLAG.NE.0.OR.TABU.EQ.0) GOTO 200
        IF(NOCM.EQ.1) THEN
          MESG = 0
        ELSE
          MESG = 1
        ENDIF
        CALL TABOUT
     I             (STDOUT, STDTAB, TABU, STATU, ZMINU, MESG, SAVOPT,
     I              OUTOPT, BETOPT, zone, hgrid, vdatum, unitsys, basis,
     i                    khflag, alphaflag, betaflag, maflag, mqflag,
     M              NDEPU, XSTU, FTP)
 200  CONTINUE
 

      RETURN
      END

C
C
C
      SUBROUTINE   ADD_SLOTE
     I                      (STDOUT, nfac, in_ss_nvar,
     M                       NPNT, NSUB, N, ZMIN, X, Z, SB, LSN, SN,
     O                       ADDED)

C     Potentially add a bottom slot to a cross section description. 
C     This is a slot with sides expanding as an exponential above 
c     a certain depth such that the hydraulic depth is a constant. 

C     See CXSTAB for definition of arguments.

      implicit none
      INCLUDE 'arsize.prm'

      INTEGER ADDED, STDOUT, NSUB, SB(PMXPNT), npnt, in_ss_nvar

      real nfac, zmin

      REAL X(PMXPNT), Z(PMXPNT), 
     A     LSN(PMXPNT), SN(PMXPNT), N(PMXSUB)
       
      INCLUDE 'abslot.cmn'
      INCLUDE 'xtadd.cmn'
                  
      integer np
      parameter (np=23)
      CHARACTER*1 
     a s_slot_match, s_cls_match, s_cls_short,
     b e_slot_match, e_cls_match, e_cls_short,
     c change_sb, retain_sb, remember_sb

      INTEGER I, J, MIN_KNT, MIN_PNT(PMXPNT), MIN_AT,
     A        ILCROSS, remember_SUB, IS, IE, EFLAG, imid, nn,
     b   i_s_cls_match, i_s_slot_match, i_e_cls_match,
     c   i_e_slot_match, i_xs_left, i_xs_right, i_start_slot,
     d   i_end_slot, n_cross, retain_sub, itype, ircross, inew, iold

      REAL  EPS_base, WIDTH(PMXPNT), MINOFF, MAXOFF, BOTTOM, XL, ZL,
     A      XR, ZR, SLOT_XL, SLOT_ZL, SLOT_XR, slot_zr, XM, ZM,
     B      SLOT_XM, SLOT_ZM, y0, ym, t0, top, y,
     c      xslot, pi, min_ls, ls, sum_ls, sum_nls, eps, shift_base,
     d      shift, zmin_new

      real slot_y(np), slot_offset(np)

      DATA EPS_base/0.0005/, pi/3.14159265/, shift_base/0.053/
C     ********************************FORMATS**************************
 50   FORMAT(' *ERR:724* Maximum number of subsections=',I5,
     A       ' exceeded in adding an invert slot.')
 52   format(/,'*ERR:* Minimum line-segment length=',1pe12.4,
     a     ' in slot is smaller',/,5x,
     b     ' than the point tolerance in add_slote=',1pe12.4)
 54   format(/,'*PROBLEM* Subset operations on cross sections with',
     a  ' roughness varying with depth not yet supported.')
C**********************************************************************
c     Set tolerance based on nfac
      if(nfac.lt.1.2) then
c       metric (si more or less)
        eps = 0.3048*eps_base
        shift = 0.3048*shift_base
      else
        eps = eps_base
        shift = shift_base
      endif
      
      IF(SLOT_PRESENT.EQ. 0) THEN
        SLOT_DEPTH = 0.0
        added = 0
        RETURN
      elseif(slot_present.ne.2) then
c       Triangular slot is present.
        return      
      ENDIF

      if(in_ss_nvar /= 0)  then
        write(stdout,54)
        stop 'Abnormal stop. Feature not yet supported.'
      endif

C     We know the minimum elevation on entry.  Do another search to 
C     locate and describe the various extremes, e. g. the width of the 
c     section at the one or more minimum points.

c     Processing a horizontal bottom becomes messy when adding a bottom slot.
c     Especially messy are horizontal bottoms represented by more then one 
c     line segment.  The additional point may be there because there is a change
c     in sub-section number or it may just be an extra point not really needed. 
c     We will try to find these horizontal segments at the current zmin and 
c     adjust them slightly to be non-horizontal. 

      i = 1
      zmin_new = -1.e30
 500  continue      
        if(abs(z(i) - zmin).le.eps) then
c         we have a point at min elevation. 
          if(abs(z(i+1) - zmin).le.eps) then
c           we have a horizontal line segment. 
c           decrement the right end.
            z(i+1) = z(i+1) - shift
            zmin_new = z(i+1)
           
            i = i + 2
            if(i.lt.npnt) goto 500
          endif
        endif
        i = i + 1
        if(i.lt.npnt) goto 500

      if(zmin_new.gt.-1.e30) then
        zmin = zmin_new
      endif
          
c      WRITE(STDOUT,*) ' In ADD_SLOT: ZMIN=',ZMIN
      MIN_KNT = 0
      MINOFF = 1.E30
      MAXOFF = -1.E30
      DO 100 I=1,NPNT
        IF(ABS(Z(I) - ZMIN).LE.EPS) THEN
C         We have a match. 
          MIN_KNT = MIN_KNT + 1
          MIN_PNT(MIN_KNT) = I
C         Does the point to the right, if it exists, have essentially 
C         the same elevation?
          IF(I.LT.NPNT) THEN
            IF(ABS(Z(I+1) - ZMIN).LE.EPS) THEN
C             Yes.  We have an essentially horizontal segment at
C             the minimum elevation.  
              WIDTH(MIN_KNT) = X(I+1) - X(I)
            ELSE
C             Bottom width to right is zero.
              WIDTH(MIN_KNT) = 0.0
            ENDIF
          ELSE
C           Point to the right does not exist.  Bottom width to
C           right of point is zero. 

            WIDTH(MIN_KNT) = 0.0
          ENDIF
        ENDIF
C       Find the extreme offsets. 
        MAXOFF = MAX(MAXOFF,X(I))
        MINOFF = MIN(MINOFF,X(I))
        
100   CONTINUE
      IF(MIN_KNT.EQ.0) THEN
        WRITE(STDOUT,*) ' BUG in ADD_SLOTE.  No minimum match found'
        STOP 'Abnormal stop' 
      ENDIF

      WRITE(STDOUT,89)
89    FORMAT('  Check of min search in add_slote:',
     a      /,'    Pointer     Width')
      DO 1234 I=1,MIN_KNT
        WRITE(STDOUT,90) MIN_PNT(I), WIDTH(I)
1234  CONTINUE
90    FORMAT(' ',I10,F10.3)


      added = 1
C     Compute the slot vertical extent. 
      IF(YSLOT.GT.0.0) THEN
C       Add a fixed size slot to the bottom of each cross section.
C       Redefine ESLOT to be the new invert.
        ESLOT = ZMIN - YSLOT
        SLOT_DEPTH = YSLOT
      ELSE
C       Add a slot that has an invert at a fixed elevation to 
C       each cross section.        
        SLOT_DEPTH = ZMIN - ESLOT
        yslot = slot_depth
        IF(SLOT_DEPTH.LT.0.0) THEN
C         No slot needed
          SLOT_DEPTH = 0.0
          RETURN
        ENDIF
      ENDIF

c     Save positive value of slot depth.  
      ym = slot_depth

c     Set the values for this slot.  The default value of rd has been defined
c     so that we will get a width at the top of the slot equal to WSLOT,
c     the same as for the triangular slot supported earlier.  The 
c     area of the exponential slot is then also rd of the area of the 
c     triangular slot.  These ratios hold for all values of ym.   
c     A value of WSLOT of 1/rd = 2.6359, gives us about the same 
c     area as the triangular slot with a unit top width. 

c     compute the depth of the base cross section
      y0 = rd*ym

c     compute the top width of the base cross section
      if(tzero.lt.0.0) then
c       use default value
        t0 = rd*wslot/10.0
      else
        t0 = tzero
      endif

c      write(stdout,*) ' y0=',y0,' t0=',t0

c     compute the slot shape relative to the zero depth point. 
c     find the middle index for the points on the slot boundary.
c     Note the number of points must always be odd so that there is 
c     a unique middle point. 
      
      nn = (np-1)/2
      imid = nn + 1
      slot_y(imid) = 0.0
      slot_offset(imid) = 0.0
      do i=1, nn

c       use the upper half of the expanded Chebyshev points
        j = nn + i - 1
        y = y0 - (ym - y0)*
     a   cos(pi*(real(2*j-1))/(real(2*(np-2))))/cos(pi/(real(2*(np-2))))
c        y = y0 + real(i-1)*(ym - y0)/real(nn-1)
        top = t0*exp(2.*expfac*(y/y0 - 1.0))
        slot_y(imid - i) = y
        slot_offset(imid - i) = -top/2.0
        slot_y(imid + i) = y
        slot_offset(imid + i) = top/2.0
      end do

      write(stdout,*) ' '
      write(stdout,*) ' Dump of slot boundary'
      write(stdout,2530)
2530  format('index sltoffset elevation')
      do i=1,np
        write(stdout,2351) i, slot_offset(i), slot_y(i)
2351  format(i5,f10.6,f10.5)
      end do

c     Find the smallest line segment length on the slot boundary. 
      min_ls = 1.e30
      do i=imid+1,np
        ls = sqrt( (slot_y(i) - slot_y(i-1))**2 + 
     a             (slot_offset(i) - slot_offset(i-1))**2)
        min_ls = min(ls, min_ls)
      enddo
      write(stdout,*) 
     a ' Minimum-length slot line segment has length=',min_ls

c     Use half of smallest line segment as the tolerance below
      min_ls = 0.5*min_ls
      if(min_ls.le.eps) then
        write(stdout,52) min_ls, eps
        stop 'Abnormal stop.  Error(s) found.'
      endif
      
C     Change sign of slot so that we can add to get the result.
C     SLOT is now the distance to the bottom of the slot from the min
C     elevation in the cross section. 
      SLOT_DEPTH = -SLOT_DEPTH


      IF(WSLOT.GE.MAXOFF-MINOFF) THEN
        WRITE(STDOUT,*) ' *ERR:725*  Bottom slot wider than section!'
        STOP 'Abnormal stop' 
      ENDIF

C     Now select the point at which to place the slot.  We have no 
C     choice if MIN_KNT = 1!

      IF(MIN_KNT.EQ.1) THEN
        MIN_AT = MIN_PNT(1)
        BOTTOM = WIDTH(1)
      ELSE
        BOTTOM = -1.E30
        DO  110 J=1,MIN_KNT
          IF(WIDTH(J).GT.BOTTOM) THEN
            MIN_AT = MIN_PNT(J)
            BOTTOM = WIDTH(J)
          ENDIF      
110     CONTINUE
      ENDIF

C      WRITE(STDOUT,*) ' MIN_AT=',MIN_AT,' BOTTOM=',BOTTOM

      XM = X(MIN_AT)
      ZM = Z(MIN_AT)

C     Now add the new subsection and adjust the vectors.
      NSUB = NSUB + 1
      IF(NSUB.GT.PMXSUB) THEN
        WRITE(STDOUT,50) PMXSUB
        STOP 'Abnormal stop'
      ENDIF

C     Set the new minimum to the slot invert elevation.

      ZMIN = ESLOT

c     clear the various flags for the variety of conditions that may 
c     exist at each end of the slot insertion.  cls-- connecting line segment!
c     s - start, e - end 
      s_slot_match = 'N'
      s_cls_match = 'N'
      s_cls_short = 'N'
      e_slot_match = 'N'
      e_cls_match = 'N'
      e_cls_short = 'N'

      IF(BOTTOM.GT.WSLOT) THEN
C       If BOTTOM is > WSLOT put the slot at the midpoint of
c       of the bottom and then go to the next section to fit
c       the slot.
        
        xm = xm + 0.5*bottom

      endif

C     The bottom is either not horizontal or horizontal and
C     too short to fit the slot requested.  We use a bit of
C     trial and error to fit the slot to the cross section.
C     Start by placing the slot centerline at the minimum
C     point.  Find the coordinates at the end points of the 
c     left-hand line segment at the start of the slot. 

      SLOT_XL = XM + slot_offset(1)
      SLOT_ZL = eslot + slot_y(1)

      SLOT_XM = XM  + slot_offset(2)
      SLOT_ZM = ESLOT + slot_y(2)

      xslot = xm
C     Seek an intersection of the slot boundary with the cross
C     section boundary on the left. 
c      WRITE(STDOUT,*) ' CALLING FIND_CROSS 1'

      CALL FIND_CROSS
     I               (STDOUT, NPNT, MIN_AT, -1, 1, X, Z, 
     I                SLOT_XM, SLOT_ZM, SLOT_XL, SLOT_ZL,
     O                ITYPE, ILCROSS, XL, ZL)      
      
      if(itype.gt.0) then
c       WRITE(STDOUT,*) ' Left cross=',ILCROSS,' ITYPE=',ITYPE
c       Establish the various conditions at this intersection. 
c       The line segment between  the leftmost point of the slot,
c       (slot_xl, slot_zl) and the point of intersection on the 
c       cross-section boundary, (xl, zl) may be shorter than 
c       we want to bother with.  Also the point of intersection may be
c       so close to an existing point on the boundary as may be the 
c       leftmost point on the slot boundary.  

c       Establish match status of the connecting line segment 
c       with respect to points on the cross section boundary. 
        IF(ABS(X(ILCROSS) - XL).LE.EPS) THEN
          i_s_cls_match = ILCROSS
          s_cls_match = 'Y'
        ELSEIF(ABS(X(ILCROSS+1) - XL).LE.EPS) THEN
          i_s_cls_match = ILCROSS + 1
          s_cls_match = 'Y'
        else
          s_cls_match = 'N'
          i_s_cls_match = ilcross
        endif

c       find length of the connecting line segment. 
        ls = sqrt( ( xl - slot_xl)**2 + (zl - slot_zl)**2)
        if(ls.le.min_ls) then
          s_cls_short = 'Y'
        endif


C       Intersection found on the left.  Seek one on the right 
C       with the slot centerline on the minimum point. 
        SLOT_XR = xslot + slot_offset(np)
        slot_zr = slot_zl
        slot_xm = xslot + slot_offset(np-1)
c        write(stdout,*) ' slot_xr=',slot_xr,' slot_xm=',slot_xm
        CALL FIND_CROSS
     I                 (STDOUT, NPNT, MIN_AT, 1, NPNT, X, Z, 
     I                  SLOT_XM, SLOT_ZM, SLOT_XR, SLOT_ZR,
     O                  ITYPE, IRCROSS, XR, ZR)      
        if(itype.gt.0) then
c         Establish conditions for the ending connecting line segment
          IF(ABS(X(IRCROSS) - XR).LE.EPS) THEN
            i_e_cls_match = IRCROSS
            e_cls_match = 'Y'
          ELSEIF(ABS(X(IRCROSS+1)-XR).LE.EPS) THEN
            i_e_cls_match = IRCROSS + 1
            e_cls_match = 'Y'
          else
            e_cls_match = 'N'
            i_e_cls_match = ircross + 1
          endif
          ls = sqrt( (xr - slot_xr)**2 + (zr - slot_zr)**2)
          if(ls.le.min_ls) then
            e_cls_short = 'Y'
          endif

c         WRITE(STDOUT,*) 
c     a    ' right cross after left found=',IRCROSS,' ITYPE=',ITYPE
        else
C         No intersection on the right when there was one on
C         the left.  Move the the slot to the left so that
C         its right point matches the minimum and then
C         seek a new intersection on the left. 

          XR = XM
          ZR = ZM
          IRCROSS = MIN_AT 
          e_slot_match = 'Y'
          i_e_slot_match = ircross

          xslot = xm + slot_offset(1)

          SLOT_XL = xslot + slot_offset(1)
          SLOT_XM = xslot + slot_offset(2)

          
          CALL FIND_CROSS
     I                   (STDOUT, NPNT, MIN_AT, -1, 1, X, Z, 
     I                   SLOT_XM, SLOT_ZM, SLOT_XL, SLOT_ZL,
     O                   ITYPE, ILCROSS, XL, ZL)      
          if(itype.gt.0) then
c           Establish conditions for the starting connecting line segment
            IF(ABS(X(ILCROSS) - XL).LE.EPS) THEN
              i_s_cls_match = ILCROSS
              s_cls_match = 'Y'
            ELSEIF(ABS(X(ILCROSS+1) - XL).LE.EPS) THEN
              i_s_cls_match = ILCROSS + 1
              s_cls_match = 'Y'
            else
              s_cls_match = 'N'
              i_s_cls_match = ilcross
            endif
    
c           find length of the connecting line segment. 
            ls = sqrt( ( xl - slot_xl)**2 + (zl - slot_zl)**2)
            if(ls.le.min_ls) then
              s_cls_short = 'Y'
            endif
          else
            WRITE(STDOUT,*) ' BUG: No intersection found on left',
     A                  ' when one must exist.'
            STOP ' Abnormal stop' 
          ENDIF

        ENDIF
      else
C       No intersection on the left.  Move slot to the right so that its
C       left point matches the minimum point and then seek 
C       an intersection on the right. 

        XL = XM
        ZL = ZM
        ILCROSS = MIN_AT
        s_slot_match = 'Y'
        i_s_slot_match = ilcross
        

        xslot = xm + slot_offset(np)

        SLOT_XM = xslot + slot_offset(np-1)
        SLOT_XR = xslot + slot_offset(np)
        slot_zr = slot_zl
        CALL FIND_CROSS
     I                 (STDOUT, NPNT, MIN_AT, 1, NPNT, X, Z, 
     I                  SLOT_XM, SLOT_ZM, SLOT_XR, SLOT_ZR,
     O                  ITYPE, IRCROSS, XR, ZR) 
        if(itype.gt.0) then
c         Establish conditions for the ending connecting line segment
          IF(ABS(X(IRCROSS) - XR).LE.EPS) THEN
            i_e_cls_match = IRCROSS
            e_cls_match = 'Y'
          ELSEIF(ABS(X(IRCROSS+1)-XR).LE.EPS) THEN
            i_e_cls_match = IRCROSS + 1
            e_cls_match = 'Y'
          else
            e_cls_match = 'N'
            i_e_cls_match = ircross + 1
          endif
          ls = sqrt( (xr - slot_xr)**2 + (zr - slot_zr)**2)
          if(ls.le.min_ls) then
            e_cls_short = 'Y'
          endif
        else
          WRITE(STDOUT,*) ' BUG: No intersection found on right',
     A                   ' when one must exist.'
          STOP ' Abnormal stop' 
        ENDIF

      endif


c     Now setup the various values we need to make room for the slot and to 
c     insert the slot properly. 

c     i_xs_left - index to the last retained point on cross section boundary to the left
c                 of the slot
c     i_xs_right - index to the first retained point on the cross-section boundary to the right 
c                  slot. 

c      both of above refer to the indices as they exist at this point.  Adjustments may have been 
c      made to make better sense of a horizontal bottom.

c     i_start_slot - starting index for the points in the slot to transfer
c     i_end_slot - ending index for the points in the slot to transfer

c     n_cross - number of intersection points to be added.  can be 0, 1, or 2

      retain_sb = 'N'
      change_sb = 'N'
      remember_sb = 'N'

      n_cross = 0
      if(s_slot_match.eq.'Y') then
c       The start point of the slot matches a boundary point exactly. 
        i_xs_left = i_s_slot_match

c       sb and lsn must be changed to slot values at i_xs_left
        change_sb = 'Y'

c       intersection point is not used
        xl = -1.e30

c       skip first point on slot because it is already on the boundary
        i_start_slot = 2                

      else
        if(s_cls_match.eq.'Y') then                                                                
c         starting connecting line segment matches a boundary point.                               
          if(s_cls_short.eq.'N') then                                                              
c           starting connecting line segment has ok length                                         
            i_xs_left = i_s_cls_match                                                              
                                                                                                   
c           sb and lsn must be changed to the slot value at i_xs_left                             
            change_sb = 'Y'
                                                                                                   
c           intersection point is not used                                                         
            xl = -1.e30                                                                            
                                                                                                   
            i_start_slot = 1                                                                       
          else                                                                                     
c           starting connecting line segment is too short.  last retained point shifts 1 to left   
            i_xs_left = i_s_cls_match - 1                                                          
                                                                                                   
c           sb and lsn are to be unchanged at i_xs_left                                            
            change_sb = 'N'
                                                                                                   
            xl = -1.e30                                                                           
            i_start_slot = 1                                                                       
          endif                                                                                    
        else                                                                                       
c         Starting connection line segment does not match a boundary point.  
          if(s_cls_short.eq.'N') then
c           Starting connecting line segment has ok length
            i_xs_left = i_s_cls_match

c           sb and lsn are not changed
            change_sb = 'N'

c           intersection point is used.  Gets slot n.
            n_cross = n_cross + 1

            i_start_slot = 1
          else
c           Starting connecting line segment is too short
            i_xs_left = i_s_cls_match

c           sb and lsn are not changed
            change_sb = 'N'

c           intersection point is not used.
            xl = - 1.e30

            i_start_slot = 1
          endif
        endif  
      endif                            

      if(e_slot_match.eq.'Y') then
c       The end point of the slot matches a boundary point exactly. 
        i_xs_right = i_e_slot_match

c       sb and lsn must be retained  at i_xs_right
        retain_sb = 'Y'

c       intersection point is not used
        xr = 1.e30

c       skip last point on slot because it is already on the boundary
        i_end_slot = np - 1                

      else
        if(e_cls_match.eq.'Y') then                                                                
c         ending connecting line segment matches a boundary point.                               
          if(e_cls_short.eq.'N') then                                                              
c           ending connecting line segment has ok length                                         
            i_xs_right = i_e_cls_match                                                              
                                                                                                   
c           sb and lsn must be retained at i_xs_right                             
            retain_sb = 'Y'
                                                                                                   
c           intersection point is not used                                                         
            xr = 1.e30                                                                            
                                                                                                   
            i_end_slot = np                                                                       
          else                                                                                     
c           ending connecting line segment is too short.   retained point shifts 1 to right   
            i_xs_right = i_e_cls_match + 1                                                          
                                                                                                   
c           remember sb and lsn  at i_xs_right - 1                                            
            remember_sb = 'Y'
                                                                                                   
            xr = 1.e30                                                                           
            i_end_slot = np                                                                       
          endif                                                                                    
        else                                                                                       
c         ending connection line segment does not match a boundary point.  
          if(e_cls_short.eq.'N') then
c           ending connecting line segment has ok length
            i_xs_right = i_e_cls_match

c           remember sb and lsn at i_xs_right - 1
            remember_sb = 'Y'

c           intersection point is used.
            n_cross = n_cross + 1

            i_end_slot = np
          else
c           ending connecting line segment is too short
            i_xs_right = i_e_cls_match

c           remember sb and lsn at i_xs_right - 1
            remember_sb = 'Y'

c           intersection point is not used.
            xr =  1.e30

            i_end_slot = np
          endif
        endif  
      endif                            
                     

      

c     Select the manning's n for the slot:
c     nslot > 0.0 -- use nslot for manning's n in the slot
c     nslot = 0  -- use the mean value of the existing n at the 
c                   start and end points.  
c     nslot < 0.0 -- use the mean value of the existing n at the 
c                    start and end point multiplied by abs(nslot).

      if(nslot.gt.0.0) then
        n(nsub) = nslot
      else
c       compute perimeter-weighted n in the vicinity of the slot
        sum_ls = 0.0
        sum_nls = 0.0
        if(i_xs_right-1.lt.i_xs_left) then
          write(stdout,*) ' problem in computing aver n in add_slote'
          write(stdout,*) 
     a      ' i_xs_left=',i_xs_left,' i_xs_right=',i_xs_right
          stop 'Abnormal stop.  Bug!'
        endif
        do i=i_xs_left,i_xs_right - 1
          ls = sqrt( (x(i+1) - x(i))**2 + (z(i+1) - z(i))**2)
          sum_ls = sum_ls + ls
          sum_nls = sum_nls + ls*n(sb(i))
        enddo
        n(nsub) = sum_nls/sum_ls
        if(nslot.lt.0.0) then
          n(nsub) = n(nsub)*abs(nslot)
        endif
      endif


      write(stdout,145) i_xs_left, i_xs_right, i_start_slot, 
     a   i_end_slot, n_cross
145   format(' i_xs_left=',i5,' i_xs_right=',i5,' i_start_slot=',i5,
     a  ' i_end_slot=',i5,' n_cross=',i5)
      write(stdout,146) xl, xr
146   format(' xl=',1pe12.4,' xr=',1pe12.4)

c     Compute the space needed for the points to be added to the boundary. 
c     we will be adding  (i_end_slot - i_start_slot + 1) + n_cross points 
c     to the boundary.  We have to the contents of x(*), etc. with 
c     iold being i_xs_right.  This point will be the first one to the right 
c     slot that is kept.  We make enough room to place the additional points
c     starting at  the location after i_xs_left.  It should always be true that
c     i_xs_right is larger than i_xs_left.

      if(i_xs_right.le.i_xs_left) then
        write(stdout,*) 
     a  ' i_xs_right=',i_xs_right,' <=',' i_xs_left=',i_xs_left  
        stop 'Abnormal stop.  Bug.'
      endif


c     Remember the sb on the right in certain cases.
      if(remember_sb.eq.'Y') then
        remember_sub = sb(i_xs_right - 1)
      endif

c     Retain sub in certain cases
      if(retain_sb.eq.'Y') then
        retain_sub = sb(i_xs_right)
      endif

c     We know that iold must be i_xs_right.  We will add some points and 
c     in some cases points will be deleted.  Thus 

c     inew = i_xs_right +
c             number of points added                    -     number of points deleted
      inew = i_xs_right + (i_end_slot - i_start_slot + 1) + n_cross 
     a            - (i_xs_right - i_xs_left - 1)

      iold = i_xs_right
      CALL SHIFT_VALUES
     I                   (STDOUT, IOLD, INEW, 
     M                    NPNT, X, Z, SB, LSN, SN)

      if(change_sb.eq.'Y') then
        sb(i_xs_left) = nsub
        lsn(i_xs_left) = n(nsub)
      endif
      is = i_xs_left + 1
      if(xl.gt.-1.e29) then
c       We use the left intersection point. It gets the slot n values 
        x(is) = xl
        z(is) = zl
        sn(is) = 1.0
        sb(is) = nsub
        lsn(is) = n(nsub)
        is = is + 1
      endif

c     transfer the points on the slot boundary 
      do i=i_start_slot,i_end_slot
        x(is+i-i_start_slot) = slot_offset(i) + xslot
        z(is+i-i_start_slot) = slot_y(i) +  eslot
        SN(IS+i-i_start_slot) = 1.0
        SB(IS+i-i_start_slot) = NSUB
        LSN(IS+i-i_start_slot) = N(NSUB)
      end do
c     compute index to the last slot point transfered
      is = is + i_end_slot - i_start_slot      

c     If the intersection point is used it gets the remembered subsection.
c     if the intersection point is not used any rememberd subsection is assigned
c     to the last point on the slot.   Retained subsection values are supossed
c     to happen without additional action!!
      if(remember_sb.eq.'Y') then
        if(xr.gt.1.e29) then
c         intersection point is not used
          sb(is) = remember_sub
          lsn(is) = n(sb(is))
        endif
      endif
      
      if(xr.lt.1.e30) then
c       add the intersection point.  Always gets the remembered subsection
        is = is + 1
        x(is) = xr 
        z(is) = zr
        sn(is) = 1.0
        sb(is) = remember_sub
        lsn(is) = n(sb(is))
      endif

      if(is+1.ne.inew) then
        write(stdout,*) 
     a' Wrong count during transfer. is+1=',is+1,' inew=',inew
        stop 'Abnormal stop. Bug.'
      endif



C     Adjust the subsection assignments
      CALL  REASUB
     I            (STDOUT, NPNT,
     M             NSUB, N, SB,
     O              EFLAG)

      RETURN
      END

c
c
c
      subroutine replace_eslot_elements(stdout, tab, ndep,
     m                              xst)

c     If the cross-section table with internal number, tab, has an 
c     exponential slot, then seek a stored cross section table with the 
c     a tabid that matches the one for tab but with a trailing lower case
c     a appended.  If no such table is available, return without doing 
c     anything.  Otherwise, 

c     1. Replace all values of alpha, beta, ma, and mq in the table
c        tab, at an above the true invert by looking up values in 
c        the stored cross-section function table.  

c     2. Change all values of these elements in the slot to match the 
c        value at that top of the slot. 

c        Write a message to the master output file noting that this 
c        operation has been done.  This operation will produce inconsistent
c        values in the display for the critical flow computations.   However, the 
c        critical flow values will be left unchanged. 

      implicit none

      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'abslot.cmn'

      integer stdout, tab, ndep

      real  xst(pmxpnt,pmxelm)


c     Called routines
      character*16 get_tabid
      external get_tabid, find_internal_tab_number, xlkt25, chktab

c     Local

      integer eflag, i, taba, adrsa

      real yzero, area, top, dt, j, k, dk, b, db, alp, dalp, 
     a    qc, ma, dma, mq, dmq, yarg

      character tabid*16, tabida*16
c     *********************************formatss*************************
50    format(/,' Table id: ',a,' with slot depth=',f10.3,' has its',
     a ' values of alpha,',/,'      beta, ma, and mq',
     b ' reset to match those in the unsloted table given in FTABIN.')
c***********************************************************************
c     Check on slot type.
      if(slot_present.ne.2) then
        return
      endif
c     Construct the name of the table that might have been stored.  
c     Get the tabid of the current table, then append a to it, and see
c     if it exists. 
      eflag = 0
      tabid = get_tabid(tab)
      i = len_trim(tabid)
      if(i.eq.0) then
        write(stdout,*)' Bug: tabid not found for table being processd'
        stop 'Abnormal stop. Bug.'
      endif

      tabida = tabid(1:i)//'a'

c     Now seek the internal table number for tabida
      call find_internal_tab_number
     i                             (tabida,
     o                              taba)

      if(taba.eq.0) then
c       User has not input a table of that name using FTABIN. 
        return
      endif

c     Check to make sure the table is of the correct type.  We now only 
c     support type 25-anything else will not do- maybe later?

      adrsa = taba
      call chktab
     i           (25, stdout, ftpnt, mftnum,
     m            adrsa,
     o            eflag)

      if(eflag.ne.0) then
        stop 'Abnormal stop.  Error(s) found.'
      endif


      write(stdout,50) tabid, yslot

c     Find the initial values from tabida

      yzero = 0.0
      call xlkt25
     I           (adrsa,
     M            yzero,
     O    area, top, dt, j, k, dk, b, db, alp, dalp, 
     o    qc, ma, dma, mq, dmq)
     

C     Added May 22, 1998: definition of XST(i,*) contents
c     Extended 13 March 2003

C     Offset    Value
C     1         Maximum depth-y
C     2         Top width
C     3         Area
C     4         First moment of area about water surface
C     5         Square root of conveyance- kh
C     6         Beta
C     7         Alpha
C     8         dBeta/dy
C     9         dAlpha/dy
C     10        Critical flow from momentum
C     11        Critical flow from energy
C     12        Critcal flow assuming Alpha=beta=1
C     13        Critical flow that is selected by user: 10, 11, or 12
C     14        MA- correction of volumes for sinuousity
C     15        MQ- correction of momentum for sinuosity
C     16        Total wetted perimeter- added May 22, 1998
C     17        Average Manning's n value for the cross section
c     18        dkh/dy computed by cubic-spline fit
c     19        dalpha/dy computed by cubic-spline fit
c     20        dbeta/dy computed by cubic-spline fit
c     21        dma/dy computed by cubic-spline fit
c     22        dmq/dy computed by cubic-spline fit

      do i=1,ndep
        yarg = xst(i,1)
        if(yarg.le.yslot) then
c         we are in the slot.  

          xst(i,6) = b
          xst(i,7) = alp
          xst(i,14) = ma
          xst(i,15) = mq
        else
c         we are out of the slot.  Redefine the values. 
          yarg = yarg - yslot
          call xlkt25
     I               (adrsa,
     M                yarg,
     O        area, top, dt, j, k, dk, b, db, alp, dalp, 
     o        qc, ma, dma, mq, dmq)

          xst(i,6) = b
          xst(i,7) = alp
          xst(i,14) = ma
          xst(i,15) = mq

        endif
      enddo

      return
      end


        


