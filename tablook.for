C     Table lookup not in comprog.for

C
C
C
      SUBROUTINE   FNDELV
     I                   (NUM, LOUT,
     O                    EFLAG, ELEV)
 
C     + + + PURPOSE + + +
C     Find the elevation of the cross section from the table
C     given by the table number in NUM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, LOUT, NUM
      REAL ELEV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     NUM    - Cross section index number
C     LOUT   - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     ELEV   - Elevation of minimum point in the cross section function
C               table
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADRS
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:198* Cross sec. tab. num.=',I5,' not found for',
     A        ' elevation.')
C***********************************************************************
      IF(NUM.LE.0) THEN
        WRITE(LOUT,50) NUM
        EFLAG = 1
        ELEV = 0.0
      ELSE
        ADRS = FTPNT(NUM)
        IF(ADRS.GT.0) THEN
          ELEV = FTAB(ADRS + 5)
        ELSE
C         TABLE DOES NOT EXIST
 
          WRITE(LOUT,50) NUM
          EFLAG = 1
          ELEV = 0.0
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   XLOOKW
     I                   (ADRS, XOFF, ADEPTH,
     O                    AREA, TOP, DTOP, C, W)
 
C     + + + PURPOSE + + +
C     Given depth find area, top-width, first moment, celerity, and
C     escoffier stage variable.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, XOFF
      REAL ADEPTH, AREA, C, DTOP, TOP, W
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADRS   - Address of function table
C     XOFF   - Offset between successive depth values for cross section
C               function table
C     ADEPTH - Depth argument
C     AREA   - Flow area
C     TOP    - Top width for the cross section
C     DTOP   - derivative of the top width with respect to depth
C     C      - Celerity
C     W      - Value of the Escoffier stage variable
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'xscom.cmn'
      INCLUDE 'grvcom.cmn'

C     Called subprograms
      CHARACTER*16 GET_TABID
      EXTERNAL GET_TABID
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HA, LA, LSTA
      REAL AZERO, CZERO, DEPTH, DTP, DY, FACT, H, TONE, TZERO, WZERO,
     A     YDIFF, YONE, YZERO
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + OUTPUT FORMATS + + +
 2000 FORMAT('0','*WRN:522* ESCOF. TABLE BELOW RANGE IN XLOOKW',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',E10.3,
     D      /,1X,' DEPTH          = ',F10.2)
 2010 FORMAT('0','*WRN:523* ESCOF. TABLE ABOVE RANGE IN XLOOKW',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',E10.3,
     D      /,1X,' DEPTH          = ',F10.2)
C***********************************************************************
      DEPTH = ADEPTH
      HA = ITAB(ADRS)
      LA = ADRS + XTIOFF
      LSTA = ITAB(ADRS+3)
 
      YDIFF = DEPTH - FTAB(LSTA)
      IF(YDIFF) 100,300,200
 
C     DEPTH PASSED LESS THAN AT PREVIOUS ACCESS TO TABLE
 
 100  IF( LSTA.GE.LA ) GOTO 110
        WRITE(LOUT,2000) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), DEPTH
        LSTA = LA
        DEPTH = FTAB(LA+XOFF)
        GOTO 300
 
 110  IF( DEPTH.GE.FTAB(LSTA) ) GOTO 300
        LSTA = LSTA - XOFF
        GOTO 100
 
C     DEPTH PASSED GREATER THAN AT PREVIOUS ACCESS TO TABLE
 
 200  IF( LSTA.LT.HA ) GOTO 210
        IF(ABS(DEPTH - FTAB(HA)).GT.0.0051) THEN
          WRITE(LOUT,2010) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), DEPTH
        ENDIF
        DEPTH = FTAB(HA)
        LSTA = HA - XOFF
        GOTO 300
 
 
 210  IF( DEPTH.LE.FTAB(LSTA+XOFF) ) GOTO 300
        LSTA = LSTA + XOFF
        GOTO 200
 
 300  CONTINUE
 
C     FETCH VALUES FROM FTAB
 
      YZERO = FTAB(LSTA)
      TZERO = FTAB(LSTA+1)
      AZERO = FTAB(LSTA+2)
      CZERO = FTAB(LSTA+3)
      WZERO = FTAB(LSTA+4)
 
      YONE = FTAB(LSTA+XOFF)
      TONE = FTAB(LSTA+XOFF + 1)
 
C     DIRECT LINEAR INTERPOLATION FOR TOP WIDTH
 
      DY = YONE - YZERO
      H = DEPTH - YZERO
      DTP =  TONE - TZERO
      FACT = H/DY
      TOP = TZERO + FACT*DTP
      DTOP = DTP/DY
      AREA = AZERO + 0.5*H*(TOP+TZERO)
      IF(DEPTH.EQ.0.0) THEN
        C = 0.0
        W = 0.0
      ELSE
        C = SQRT(GRAV*AREA/TOP)
        W = WZERO + 2.*GRAV*H/(C + CZERO)
      ENDIF
 
C     RESET POINTER FOR LAST ADDRESS
 
      ITAB(ADRS+3) = LSTA
 
      RETURN
      END
C
C
C
      SUBROUTINE   LKTK
     I                 (ADRS,
     M                  YA,
     O                  K)
 
C     + + + PURPOSE + + +
C     Given depth lookup conveyance in a cross section table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS
      REAL K, YA
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADRS   - Address of function table
C     YA     - maximum depth in a cross section
C     K      - conveyance
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'xscom.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER VTYPE(35)
      SAVE VTYPE
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HA, L, LA, TYPE, XOFF
      REAL DK, DY, K0, Y, Y0
 
C     + + + EXTERNAL NAMES + + +
      CHARACTER*16 GET_TABID
      EXTERNAL GET_TABID, XSTYPE
 
C     + + + DATA INITIALIZATIONS + + +
      DATA VTYPE/1,10*0,1,7*0,6*1,4*0,6*1/
 
C     + + + OUTPUT FORMATS + + +
 2000 FORMAT('0','*WRN:41* X-SECTION BELOW RANGE IN LKTK',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',F10.3
     C      /,1X,' TIME           = ',F10.0,
     D      /,1X,' DEPTH          = ',F10.2)
 2010 FORMAT('0','*WRN:42* X-SECTION ABOVE RANGE IN LKTK',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',F10.3
     C      /,1X,' TIME           = ',F10.0,
     D      /,1X,' DEPTH          = ',F10.2)
C***********************************************************************
C     HA = HIGH ADDRESS
C     LA = LOW ADDRESS
C     L = ADDRESS FOUND ON THE LAST CALL TO THIS TABLE
 
      Y = YA
      HA = ITAB(ADRS)
      LA = ADRS + XTIOFF
      L = ITAB(ADRS+3)
 
      TYPE = ITAB(ADRS+2)
      XOFF = OFFVEC(TYPE)
      IF(VTYPE(TYPE).EQ.0) THEN
        CALL XSTYPE
     I             (LOUT, VTYPE, ITAB(ADRS+1))
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      IF(Y.GE.FTAB(L)) THEN
C       CHECK FOR ARGUMENT ABOVE MAX ARG IN THE TABLE
        IF(Y.GT.FTAB(HA)) THEN
          WRITE(LOUT,2010) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), 
     A                     TIME, Y
          L = HA - XOFF
          Y = FTAB(HA)
          YA = Y
        ELSE
 100      CONTINUE
            IF(Y.GT.FTAB(L+XOFF)) THEN
              L = L + XOFF
              GOTO 100
            ENDIF
        ENDIF
      ELSE
C       CHECK FOR ARGUMENT BELOW MIN ARG IN THE TABLE
        IF(Y.LT.FTAB(LA)) THEN
          WRITE(LOUT,2000) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), 
     A                     TIME, Y
          L = LA
          Y = FTAB(L+XOFF)
          YA = Y
        ELSE
 110      CONTINUE
            L = L - XOFF
            IF(Y.LT.FTAB(L)) GOTO 110
        ENDIF
      ENDIF
C     AT THIS POINT L DEFINES THE LOW ARGUMENT END OF THE
C     INTERVAL CONTAINING THE ARGUMENT, PERHAPS ADJUSTED
C     FOR ARGUMENT OUT OF RANGE.
 
C     RESET POINTER FOR LAST ADDRESS
 
      ITAB(ADRS+3) = L
 
C     FETCH VALUES FROM FTAB
 
      Y0 = FTAB(L)
      K0 = FTAB(L+3)
 
      DY = FTAB(L+XOFF) - Y0
 
      DK = (FTAB(L+XOFF+3) - K0)/DY
      K = K0 + (Y - Y0)*DK
      K = K*K
 
      RETURN
      END
C
C
C
      INTEGER FUNCTION   LOCSTA
     I                         (STDOUT, STAT, DIR, JAXIS, NUMSEC, STL,
     I                          EPS)
 
C     + + + PURPOSE + + +
C     Find STAT in STL(*,*).  Must agree closely because no
C     interpolation is defined.  Return value of index.  DIR
C     gives the direction of the station values in STL.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER DIR, JAXIS, NUMSEC, STDOUT
      REAL EPS, STAT, STL(PMXSEC,PMXNFL)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     STAT   - Station being sought
C     DIR    - If DIR > 0 then stations are ascending order, else
C               descending order
C     JAXIS  - Column in which the values for the channel axis are
C               stored
C     NUMSEC - Number of cross sections in sequence
C     STL    - Table for flow line stations
C     EPS    - Tolerance for matching stations
 
C     + + + SAVED VALUES + + +
      INTEGER L
      SAVE L
 
C     + + + LOCAL VARIABLES + + +
      INTEGER LOC
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FSTDEC, FSTINC
 
C     + + + DATA INITIALIZATIONS + + +
      DATA L/1/
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:678* Cross section station=',F10.4,' not found.')
C***********************************************************************
      IF(DIR.GT.0) THEN
C       STATIONS IN STL ARE ASCENDING
 
C       MAKE SURE L IS IN THE RIGHT RANGE
        IF(L.GT.NUMSEC) L = 1
 
        IF(STAT.GT.STL(NUMSEC,JAXIS).OR.STAT.LT.STL(1,JAXIS)) THEN
C         REQUESTED STATION FALLS OUTSIDE THE LIMITS OF THE
C         DEFINED SINUOUSITIES
          WRITE(STDOUT,50) STAT
          LOCSTA = 0
        ELSE
          CALL FSTINC
     I               (STAT, JAXIS, STL, EPS,
     M                L,
     O                LOC)
          IF(LOC.EQ.0) THEN
            WRITE(STDOUT,50) STAT
            LOCSTA = 0
          ELSE
 
            LOCSTA = LOC
          ENDIF
        ENDIF
      ELSE
C       STATIONS IN STL ARE DESCENDING
 
C       MAKE SURE L IS IN THE RIGHT RANGE
        IF(L.GT.NUMSEC) L = 1
 
        IF(STAT.LT.STL(NUMSEC,JAXIS).OR.STAT.GT.STL(1,JAXIS)) THEN
C         REQUESTED STATION FALLS OUTSIDE THE LIMITS OF THE
C         DEFINED SINUOUSITIES
          WRITE(STDOUT,50) STAT
          LOCSTA = 0
        ELSE
          CALL FSTDEC
     I               (STAT, JAXIS, STL, EPS,
     M                L,
     O                LOC)
          IF(LOC.EQ.0) THEN
            WRITE(STDOUT,50) STAT
            LOCSTA = 0
          ELSE
            LOCSTA = LOC
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   XLKA
     I                 (Y, ND, YV, TV, AV,
     O                  T, A)
 
C     + + + PURPOSE + + +
C     Find T and A in a table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ND
      REAL A, AV(ND), T, TV(ND), Y, YV(ND)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - ordinate values on boundary of conduit
C     ND     - Number of tabulated values
C     YV     - Vector of depth values
C     TV     - Vector of top width values
C     AV     - Vector of areas
C     T      - top width of the cross section
C     A      - Cross sectional area
 
C     + + + LOCAL VARIABLES + + +
      INTEGER L, LL
      REAL H, P, PH
C***********************************************************************
      L = 1
 100  CONTINUE
 
        IF(Y.GE.YV(L)) THEN
          LL = L + 1
          IF(LL.GT.ND) THEN
            T = TV(ND)
            A = AV(ND)
            L = ND -1
            RETURN
          ENDIF
          IF(Y.LE.YV(LL)) GOTO 200
          L = LL
          GOTO 100
        ELSE
          L = L - 1
          IF(L.LT.1) THEN
            T = TV(1)
            A = AV(1)
            L = 1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
 
 200  CONTINUE
 
 
C     INTERPOLATE
 
      H=YV(L+1)-YV(L)
      PH=Y-YV(L)
      P=PH/H
      T=TV(L)+P*(TV(L+1)-TV(L))
      A=AV(L)+PH*(T+TV(L))/2.
      RETURN
      END
C
C
C
      SUBROUTINE   XLKALL
     I                   (NDEP, MAXPNT, XST,
     M                    XSV)
 
C     + + + PURPOSE + + +
C     Find elements in table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER MAXPNT, NDEP
      REAL XST(MAXPNT,*), XSV(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     NDEP   - Number of depth values
C     MAXPNT - Maximum number of tabulated values in a cross section
C               function table
C     XST    - Storage table for various elements of cross section
C     XSV    - Vector of various elements of cross section
 
C     + + + SAVED VALUES + + +
      INTEGER L
      SAVE L
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
      REAL H, K, P, PH
 
C     + + + DATA INITIALIZATIONS + + +
      DATA L/1/
C***********************************************************************
      IF(L.GT.NDEP) L=1
 100  CONTINUE
        IF(XSV(1).GE.XST(L,1)) GOTO 120
          IF(L.GT.1) GOTO 110
            DO 102 J=2,7
              XSV(J)=XST(1,J)
 102          CONTINUE
            RETURN
 110      CONTINUE
            L=L-1
            GOTO 100
 120    CONTINUE
 
 130  CONTINUE
        IF(XSV(1).LT.XST(L+1,1)) GOTO 150
          IF(L+1.LT.NDEP) GOTO 140
            DO 132 J=2,7
              XSV(J)=XST(NDEP,J)
 132          CONTINUE
            RETURN
 140      CONTINUE
            L=L+1
            GOTO 130
 150    CONTINUE
 
C      INTERPOLATE FOR THE VALUES
 
      H=XST(L+1,1)-XST(L,1)
      PH=XSV(1)-XST(L,1)
      P=PH/H
 
      XSV(2)=XST(L,2)+P*(XST(L+1,2)-XST(L,2))
      XSV(3)=XST(L,3)+PH*(XST(L,2)+XSV(2))/2.
      XSV(4)=XST(L,4)+PH*((XST(L,3)+XSV(3))/2.-PH*(XSV(2)-XST(L,2))/12.)
      K=XST(L,5)+P*(XST(L+1,5)-XST(L,5))
      XSV(5)=K*K
      XSV(6)=XST(L,6)+P*(XST(L+1,6)-XST(L,6))
      XSV(7)=XST(L,7)+P*(XST(L+1,7)-XST(L,7))
      RETURN
      END
C
C
C
      SUBROUTINE   LKTA
     I                 (ADRS,
     M                  YA,
     O                  A)
 
C     + + + PURPOSE + + +
C     Given depth lookup area in a cross section table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS
      REAL A, YA
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADRS   - Address of function table
C     YA     - depth argument
C     A      - Cross sectional area
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'xscom.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER VTYPE(35)
      SAVE VTYPE
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HA, L, LA, TYPE, XOFF
      REAL A0, DT, DY, H, T, T0, Y, Y0
 
C     + + + EXTERNAL NAMES + + +
      CHARACTER*16 GET_TABID
      EXTERNAL GET_TABID, XSTYPE
 
C     + + + DATA INITIALIZATIONS + + +
      DATA VTYPE/1,10*0,1,7*0,6*1,4*0,6*1/
 
C     + + + OUTPUT FORMATS + + +
 2000 FORMAT('0','*WRN:41* X-SECTION BELOW RANGE IN LKTA',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',F10.3
     C      /,1X,' TIME           = ',F10.0,
     D      /,1X,' DEPTH          = ',F10.2)
 2010 FORMAT('0','*WRN:42* X-SECTION ABOVE RANGE IN LKTA',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',F10.3
     C      /,1X,' TIME           = ',F10.0,
     D      /,1X,' DEPTH          = ',F10.2)
C***********************************************************************
C     HA = HIGH ADDRESS
C     LA = LOW ADDRESS
C     L = ADDRESS FOUND ON THE LAST CALL TO THIS TABLE
 
      Y = YA
      HA = ITAB(ADRS)
      LA = ADRS + XTIOFF
      L = ITAB(ADRS+3)
 
      TYPE = ITAB(ADRS+2)
      XOFF = OFFVEC(TYPE)
      IF(VTYPE(TYPE).EQ.0) THEN
        CALL XSTYPE
     I             (LOUT, VTYPE, ITAB(ADRS+1))
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      IF(Y.GE.FTAB(L)) THEN
C       CHECK FOR ARGUMENT ABOVE MAX ARG IN THE TABLE
        IF(Y.GT.FTAB(HA)) THEN
          WRITE(LOUT,2010) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), 
     A                     TIME, Y
          L = HA - XOFF
          Y = FTAB(HA)
C          YA = Y
        ELSE
 100      CONTINUE
            IF(Y.GT.FTAB(L+XOFF)) THEN
              L = L + XOFF
              GOTO 100
            ENDIF
        ENDIF
      ELSE
C       CHECK FOR ARGUMENT BELOW MIN ARG IN THE TABLE
        IF(Y.LT.FTAB(LA)) THEN
          WRITE(LOUT,2000) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), 
     A                     TIME, Y
          L = LA
          Y = FTAB(L+XOFF)
C          YA = Y
        ELSE
 110      CONTINUE
            L = L - XOFF
            IF(Y.LT.FTAB(L)) GOTO 110
        ENDIF
      ENDIF
C     AT THIS POINT L DEFINES THE LOW ARGUMENT END OF THE
C     INTERVAL CONTAINING THE ARGUMENT, PERHAPS ADJUSTED
C     FOR ARGUMENT OUT OF RANGE.
 
C     RESET POINTER FOR LAST ADDRESS
 
      ITAB(ADRS+3) = L
 
C     FETCH VALUES FROM FTAB
 
      Y0 = FTAB(L)
      T0 = FTAB(L+1)
      A0 = FTAB(L+2)
 
 
      DY = FTAB(L+XOFF) - Y0
      H = Y - Y0
      DT =  (FTAB(L+XOFF+1) - T0)/DY
      T = T0 + H*DT
      A = A0 + 0.5*H*(T + T0)
 
      RETURN
      END
C
C
C
      SUBROUTINE   LKTJ
     I                 (ADRS,
     M                  YA,
     O                  J)
 
C     + + + PURPOSE + + +
C     Given depth lookup first moment of area in a cross section table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS
      REAL J, YA
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADRS   - Address of function table
C     YA     - depth argument
C     J      - first moment of area about water surface
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'xscom.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER VTYPE(35)
      SAVE VTYPE
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HA, L, LA, TYPE, XOFF
      REAL A, A0, DT, DY, H, HH, J0, T, T0, Y, Y0
 
C     + + + EXTERNAL NAMES + + +
      CHARACTER*16 GET_TABID
      EXTERNAL GET_TABID, XSTYPE
 
C     + + + DATA INITIALIZATIONS + + +
      DATA VTYPE/1,10*0,1,7*0,0,1,1,0,1,1,4*0,0,1,1,0,1,1/
 
C     + + + OUTPUT FORMATS + + +
 2000 FORMAT('0','*WRN:41* X-SECTION BELOW RANGE IN LKTJ',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',F10.3
     C      /,1X,' TIME           = ',F10.0,
     D      /,1X,' DEPTH          = ',F10.2)
 2010 FORMAT('0','*WRN:42* X-SECTION ABOVE RANGE IN LKTJ',
     A      /,1X,' TABLE ID   = ',A,
     B      /,1X,' STATION NUMBER = ',F10.3
     C      /,1X,' TIME           = ',F10.0,
     D      /,1X,' DEPTH          = ',F10.2)
C***********************************************************************
C     HA = HIGH ADDRESS
C     LA = LOW ADDRESS
C     L = ADDRESS FOUND ON THE LAST CALL TO THIS TABLE
 
      Y = YA
      HA = ITAB(ADRS)
      LA = ADRS + XTIOFF
      L = ITAB(ADRS+3)
 
      TYPE = ITAB(ADRS+2)
      XOFF = OFFVEC(TYPE)
      IF(VTYPE(TYPE).EQ.0) THEN
        CALL XSTYPE
     I             (LOUT, VTYPE, ITAB(ADRS+1))
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      IF(Y.GE.FTAB(L)) THEN
C       CHECK FOR ARGUMENT ABOVE MAX ARG IN THE TABLE
        IF(Y.GT.FTAB(HA)) THEN
          WRITE(LOUT,2010) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), 
     A                     TIME, Y
          L = HA - XOFF
          Y = FTAB(HA)
          YA = Y
        ELSE
 100      CONTINUE
            IF(Y.GT.FTAB(L+XOFF)) THEN
              L = L + XOFF
              GOTO 100
            ENDIF
        ENDIF
      ELSE
C       CHECK FOR ARGUMENT BELOW MIN ARG IN THE TABLE
        IF(Y.LT.FTAB(LA)) THEN
          WRITE(LOUT,2000) GET_TABID(ITAB(ADRS+1)), FTAB(ADRS+4), 
     A                     TIME, Y
          L = LA
          Y = FTAB(L+XOFF)
          YA = Y
        ELSE
 110      CONTINUE
            L = L - XOFF
            IF(Y.LT.FTAB(L)) GOTO 110
        ENDIF
      ENDIF
C     AT THIS POINT L DEFINES THE LOW ARGUMENT END OF THE
C     INTERVAL CONTAINING THE ARGUMENT, PERHAPS ADJUSTED
C     FOR ARGUMENT OUT OF RANGE.
 
C     RESET POINTER FOR LAST ADDRESS
 
      ITAB(ADRS+3) = L
 
C     FETCH VALUES FROM FTAB
 
      Y0 = FTAB(L)
      T0 = FTAB(L+1)
      A0 = FTAB(L+2)
      J0 = FTAB(L+5)
 
 
      DY = FTAB(L+XOFF) - Y0
      H = Y - Y0
      HH = 0.5*H
      DT =  (FTAB(L+XOFF+1) - T0)/DY
      T = T0 + H*DT
      A = A0 + HH*(T + T0)
      J = J0 + HH*(A + A0 -H*(T - T0)/6.)
 
      RETURN
      END
