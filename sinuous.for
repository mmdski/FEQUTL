C
C
C
      SUBROUTINE   ADJSIN
     I                   (MAXNFL, VARTYP, NUMSEC,
     M                    OFFSET, SINU,
     O                    NUMOFF)
 
C     + + + PURPOSE + + +
C     Adjust the offset and sinuosity values to make later
C     operations simpler.  We will extend the definition of
C     sinuosity to practival infinity in both directions.
C     Infinity for this purpose is defined as 1.E7.  This is far larger
C     than any offset in either english or metric units should
C     ever be.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER NUMSEC, VARTYP
      INTEGER MAXNFL(PMXSEC), NUMOFF(PMXSEC)
      REAL OFFSET(PMXSEC,PMXNFL), SINU(PMXSEC, PMXNFL)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     MAXNFL - Maximum number of flow lines at each section
C     VARTYP - If VARTYP=1 then piecewise linear variation of sinuousity
C               with offset in a cross section is assumed; else
C               if VARTYP=2 then piecewise constant variation of
C               sinuousity with offset in a cross section is assumed.
C     NUMSEC - Number of cross sections in sequence
C     OFFSET - Offsets for the sinuousity values
C     SINU   - Table of sinuousities
C     NUMOFF - Number of offsets
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, M, N
      REAL OFF(PMXNFL), SIN(PMXNFL)
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL SORT2R
C***********************************************************************
      IF(VARTYP.EQ.1) THEN
C       PIECEWISE LINEAR VARIATION.  THE SAME NUMBER OF OFFSETS AS
C       FLOW LINES(SINUOUSITIES) EXISTS.  COLLECT THE VALUES FROM
C       THE ARRAYS, SORT THEM, THEN TRANSFER THEM BACK TO
C       THE ARRAYS ADDING THE EXTRA POINTS AND ADJUST THE VALUE OF
C       NUMOFF().
 
        DO 500 I=1,NUMSEC
C         TRANSFER THE OFFSETS AND SINUOSITIES TO WORK SPACE
          N = 0
          DO 100 J=1,PMXNFL
            IF(OFFSET(I,J).GT.-1.E29) THEN
C             A VALUE EXISTS
              N = N + 1
              OFF(N) = OFFSET(I,J)
              SIN(N) = SINU(I,J)
              IF(N.EQ.MAXNFL(I)) GOTO 101
            ENDIF
 100      CONTINUE
 101      CONTINUE
 
C         NOW SORT THE VALUES TO MAKE SURE THEY ARE IN ASCENDING ORDER
C         OF OFFSET.
          CALL SORT2R
     I               (N,
     M                OFF, SIN)
 
C         TRANSFER THEM BACK TO THE ARRAYS.
C         ADD A FIRST POINT, SMALLER THAN ANY THE USER WILL USE.
          M = 1
          OFFSET(I,M) = -1.E7
          SINU(I,M) = SIN(1) + (OFFSET(I,M) - OFF(1))*
     A                         (SIN(2) - SIN(1))/(OFF(2) - OFF(1))
          IF(SINU(I,M).LT.0.0) THEN
            SINU(I,M) = 0.1
          ELSEIF(SINU(I,M).GT.10.0) THEN
            SINU(I,M) = 10.0
          ENDIF
          DO 200 J=1,N
            M = M + 1
            OFFSET(I,M) = OFF(J)
            SINU(I,M) = SIN(J)
 200      CONTINUE
C         ADD A LAST POINT LARGER THAN ANY THE USER WILL USE
          M = M + 1
          OFFSET(I,M) = 1.E7
          SINU(I,M) = SIN(N) + (OFFSET(I,M) - OFF(N))*
     A                      (SIN(N) - SIN(N-1))/(OFF(N) - OFF(N-1))
          IF(SINU(I,M).LT.0.0) THEN
            SINU(I,M) = 0.1
          ELSEIF(SINU(I,M).GT.10.0) THEN
            SINU(I,M) = 10.0
          ENDIF
          NUMOFF(I) = M
 500    CONTINUE
      ELSEIF(VARTYP.EQ.2) THEN
C       PIECEWISE CONSTANT VARIATION OF SINUOSITY.
C       THE NUMBER OF OFFSETS IS ONE LESS THAN THE NUMBER OF
C       SINUOUSITIES.  THE OFFSET FOR THE SINUOSITY NOT HAVING
C       A OFFSET SHOULD BE SET TO 1.E7.  WE THEN ADD A FIRST POINT
C       AS FOR PIECEWISE LINEAR BUT WE DO NOT NEED TO INTERPOLATE
C       FOR SINUOSITY.
 
        DO 1000 I=1,NUMSEC
C         TRANSFER THE VALUES TO WORK SPACE AND SET THE MISSING OFFSET.
          N = 0
          DO 600 J=1,PMXNFL
            IF(SINU(I,J).GT.-1.E27) THEN
C             FOUND A VALUE
              N = N + 1
              SIN(N) = SINU(I,J)
              IF(OFFSET(I,J).LT.-1.E29) THEN
C               VALUE IS MISSING
                OFF(N) = 1.E7
              ELSE
                OFF(N) = OFFSET(I,J)
              ENDIF
            ENDIF
 600      CONTINUE
 
C         SORT ON OFFSET
          CALL SORT2R
     I               (N,
     M                OFF, SIN)
 
C         TRANSFER THEM BACK TO THE ARRAYS.
C         ADD A FIRST POINT, SMALLER THAN ANY THE USER WILL USE.
          M = 1
          OFFSET(I,M) = -1.E7
          SINU(I,M) = SIN(1)
          DO 700 J=1,N
            M = M + 1
            OFFSET(I,M) = OFF(J)
            SINU(I,M) = SIN(J)
 700      CONTINUE
          NUMOFF(I) = M
 
 1000   CONTINUE
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   CPSINU
     I                   (STDOUT, NUMSEC, JAXIS, FLNTAB, STL, NFLNAM,
     M                    SINU,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Compute the sinuosities not yet known in the sinuosity matrix.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, JAXIS, NFLNAM, NUMSEC, STDOUT
      REAL SINU(PMXSEC,PMXNFL), STL(PMXSEC, PMXNFL)
      CHARACTER FLNTAB(NFLNAM)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NUMSEC - Number of cross sections in sequence
C     JAXIS  - Column in which the values for the channel axis are
C               stored
C     FLNTAB - Flow line name table
C     STL    - Table for flow line stations
C     NFLNAM - Number of flow line names
C     SINU   - Table of sinuousities
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IEND, ISTART, JCOL, VALID
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL PRSRUN
C***********************************************************************
C     FOR EACH NON-AXIS COLUMN, SCAN FOR CONSECUTIVE RUNS OF
C     STATION VALUES.  THEN FOR EACH RUN OF STATION VALUES
C     DEFINE THE LOCAL VALUES AND SCAN FOR CUBIC SPLINE REQUESTS.
C     IF CUBIC SPLINE REQUEST IS FOUND, RESCAN THE RUN OF STATION
C     VALUES AND COMPUTE THE SPLINE VALUES.  THERE MAY BE SPECIFIED
C     SINUOUSITIES INTERIOR TO THE RUN OF STATIONS.
 
      DO 1000 JCOL=1,NFLNAM
        IF(JCOL.NE.JAXIS) THEN
C         NON-AXIS COLUMN.
 
C         CLEAR THE FLAG FOR A VALID RUN OF STATIONS
          VALID = 0
 
          DO 900 I=1,NUMSEC
C           SEARCH FOR A CONSECUTIVE RUN OF STATIONS IN THE CURRENT
C           COLUMN
 
            IF(STL(I,JCOL).GT.-1.E29) THEN
C             A VALID STATION VALUE EXISTS
              IF(VALID.EQ.0) THEN
C               NO RUN IN PROGRESS.  START ONE.
                ISTART = I
                VALID = 1
              ELSE
C               RUN IN PROGRESS.  SEE IF WE ARE AT THE END OF THE
C               CURRENT COLUMN.
                IF(I.EQ.NUMSEC) THEN
C                 WE ARE AT THE END OF THE COLUMN. END THE RUN
C                 AND PROCESS IT.
 
                  IEND = I
 
                  CALL PRSRUN
     I                       (STDOUT, JAXIS, JCOL, FLNTAB, ISTART, IEND,
     I                        STL, NFLNAM,
     M                        SINU,
     O                        EFLAG)
                ENDIF
              ENDIF
            ELSE
C             VALID STATION VALUE DOES NOT EXIST
              IF(VALID.EQ.1) THEN
C               RUN IN PROGRESS.  END OF A RUN.  END IT AND PROCESS IT.
                IEND = I - 1
                VALID = 0
 
                CALL PRSRUN
     I                     (STDOUT, JAXIS, JCOL, FLNTAB, ISTART, IEND,
     I                      STL, NFLNAM,
     M                      SINU,
     O                      EFLAG)
 
              ENDIF
            ENDIF
 900      CONTINUE
        ENDIF
 1000 CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   LKUPSN
     I                   (LOC, NPNT, X, VARTYP, SINU, OFFSET, NOFF,
     O                    SNVEC)
 
C     + + + PURPOSE + + +
C     Lookup the values of sinuosity at the given location.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER LOC, NOFF, NPNT, VARTYP
      REAL OFFSET(PMXSEC,PMXNFL), SINU(PMXSEC,PMXNFL), SNVEC(NPNT),
     A     X(NPNT)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     LOC    - Index giving location of the current sinuousity data
C     NPNT   - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     VARTYP - If VARTYP=1 then piecewise linear variation of sinuousity
C               with offset in a cross section is assumed; else
C               if VARTYP=2 then piecewise constant variation of
C               sinuousity with offset in a cross section is assumed.
C     SINU   - Table of sinuousities
C     OFFSET - Offsets for the sinuousity values
C     NOFF   - Number of offsets
C     SNVEC  - Sinuousity at a point on a cross section boundary
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J, K, L
      REAL ARG, OFF(PMXNFL), SIN(PMXNFL)
C***********************************************************************
C     TRANSFER DATA TO LOCAL WORKSPACE
      DO 100 J=1,NOFF
        OFF(J) = OFFSET(LOC,J)
        SIN(J) = SINU(LOC,J)
 100  CONTINUE
 
      L = 1
      IF(VARTYP.EQ.2) THEN
C       PIECEWISE CONSTANT SINUOSITY
        DO 200 K=1,NPNT-1
C         USE THE MID-POINT OF THE LINE SEGMENT TO AVOID PROBLEMS AT
C         THE BOUNDARY BETWEEN SUBAREAS
 
          ARG = 0.5*(X(K) + X(K+1))
          IF(ARG.GE.OFF(L)) THEN
 210      CONTINUE
            IF(ARG.GT.OFF(L+1)) THEN
              L = L + 1
              GOTO 210
            ENDIF
          ELSE
 220        CONTINUE
              L = L - 1
              IF(ARG.LT.OFF(L)) GOTO 220
          ENDIF
 
C         L POINTS TO THE LEFT END OF THE INTERVAL THAT CONTAINS ARG.
C         FOR PIECEWISE CONSTANT VARIATION THE SINUOSITY VALUE IS
C         GIVEN BY INDEX L+1
          SNVEC(K) = SIN(L+1)
 200    CONTINUE
        SNVEC(NPNT) = SNVEC(NPNT-1)
      ELSE
C       PIECEWISE LINEAR SINUOSITY
        DO 300 K=1,NPNT
          ARG = X(K)
          IF(ARG.GE.OFF(L)) THEN
 310      CONTINUE
            IF(ARG.GT.OFF(L+1)) THEN
              L = L + 1
              GOTO 310
            ENDIF
          ELSE
 320        CONTINUE
              L = L - 1
              IF(ARG.LT.OFF(L)) GOTO 320
          ENDIF
 
C         L POINTS TO THE LEFT END OF THE INTERVAL THAT CONTAINS ARG.
C         INTERPOLATE FOR THE SINUOSITY
          SNVEC(K) = SIN(L) + (ARG - OFF(L))*(SIN(L+1) - SIN(L))
     A                    /(OFF(L+1) - OFF(L))
 
 300    CONTINUE
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   FSTDEC
     I                   (STAT, JAXIS, STL, EPS,
     M                    L,
     O                    I)
 
C     + + + PURPOSE + + +
C     Find station when STL has stations in descending order.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER I, JAXIS, L
      REAL EPS, STAT, STL(PMXSEC,PMXNFL)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STAT   - Station being sought
C     JAXIS  - Column in which the values for the channel axis are
C               stored
C     STL    - Table for flow line stations
C     EPS    - Tolerance for matching stations
C     L      - Starting index for search
C     I      - Index value for the station
C***********************************************************************
      IF(STAT.LE.STL(L,JAXIS)) THEN
 100    CONTINUE
          IF(STAT.LT.STL(L+1,JAXIS)) THEN
            L = L + 1
            GOTO 100
          ENDIF
      ELSE
 110    CONTINUE
          L = L - 1
          IF(STAT.GT.STL(L,JAXIS)) GOTO 110
      ENDIF
 
      IF(STAT.LE.STL(L+1,JAXIS)+EPS.AND.STAT.GE.STL(L+1,JAXIS)-EPS) THEN
        I = L + 1
      ELSE
        IF(STAT.LE.STL(L,JAXIS)+EPS.AND.STAT.GE.STL(L,JAXIS)-EPS) THEN
          I = L
        ELSE
          I = 0
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   FSTINC
     I                   (STAT, JAXIS, STL, EPS,
     M                    L,
     O                    I)
 
C     + + + PURPOSE + + +
C     Find station when STL has stations in ascending order.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER I, JAXIS, L
      REAL EPS, STAT, STL(PMXSEC,PMXNFL)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STAT   - Station being sought
C     JAXIS  - Column in which the values for the channel axis are
C               stored
C     STL    - Table for flow line stations
C     EPS    - Tolerance for matching stations
C     L      - Starting index for search
C     I      - Index value for the station
C***********************************************************************
      IF(STAT.GE.STL(L,JAXIS)) THEN
 100    CONTINUE
          IF(STAT.GT.STL(L+1,JAXIS)) THEN
            L = L + 1
            GOTO 100
          ENDIF
      ELSE
 110    CONTINUE
          L = L - 1
          IF(STAT.LT.STL(L,JAXIS)) GOTO 110
      ENDIF
 
      IF(STAT.LE.STL(L,JAXIS)+EPS.AND.STAT.GE.STL(L,JAXIS)-EPS) THEN
        I = L
      ELSE
        L = L + 1
        IF(STAT.LE.STL(L,JAXIS)+EPS.AND.STAT.GE.STL(L,JAXIS)-EPS) THEN
          I = L
        ELSE
          I = 0
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   PRSRUN
     I                   (STDOUT, JAXIS, JCOL, FLNTAB, ISTART, IEND,
     I                    STL, NFLNAM,
     M                    SINU,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Process a station run.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, IEND, ISTART, JAXIS, JCOL, NFLNAM, STDOUT
      REAL SINU(PMXSEC,PMXNFL), STL(PMXSEC, PMXNFL)
      CHARACTER FLNTAB(NFLNAM)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     JAXIS  - Column in which the values for the channel axis are
C               stored
C     JCOL   - Column in which the values for the channel axis are
C               stored
C     FLNTAB - Flow line name table
C     ISTART - Starting index for a station run
C     IEND   - Ending index for a station run
C     STL    - Table for flow line stations
C     NFLNAM - Number of flow line names
C     SINU   - Table of sinuousities
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CUBIC, J, K, KEND, KSTART, LCODE, N, RCODE, SPFLAG
      REAL H1, H2, HSUM
      DOUBLE PRECISION LVAL, RVAL, X(PMXSEC), Y(PMXSEC), YDOT(PMXSEC)
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL SPLINE
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:568* Only one station exists for flow line:',A8,
     A        'at station:',F10.4)
 52   FORMAT(/,' *ERR:606* Upstream-linear sinuosity impossible for ',
     A       ' flow line:',A8,/,10X,' at station:',F10.4,'  No ',
     B       ' upstream station exists.')
 54   FORMAT(/,' *ERR:649* Down stream-linear sinuosity impossible',
     A       ' for flow line:',A8,/,10X,' at station:',F10.4,'  No ',
     B       ' down stream station exists.')
 56   FORMAT(/,' *ERR:650* `Parabolic sinuosity impossible',
     A       ' for flow line:',A8,/,10X,' at station:',F10.4,
     B       ' Three consecutive stations needed.')
 58   FORMAT(/,' *BUG:XXX* Invalid code for flow line:',A8,' at ',
     A       'station:',F10.4,' code=',1PE15.7)
 60   FORMAT(/,' *ERR:697* Negative sinuosity found at station:',
     A       F10.4)
C***********************************************************************
C     SCAN FOR LOCAL SINUOSITY DEFINITION REQUESTS AND DO THEM.
C     SET THE CUBIC FLAG IF ANY CUBIC SPLINE REQUESTS ARE FOUND.
 
      IF(ISTART.EQ.IEND) THEN
C       INVALID RUN.  ONE POINT IS NOT SUFFICIENT TO DEFINE
C       SINUOSITY.
 
        WRITE(STDOUT,50) FLNTAB(JCOL), STL(ISTART,JCOL)
        EFLAG = 1
        RETURN
      ENDIF
 
 
C     CLEAR THE CUBIC SPLINE REQUEST FLAG
      CUBIC = 0
 
      DO 1000 K=ISTART,IEND
 
        IF(SINU(K,JCOL).EQ.-1.E28) THEN
C         REQUEST FOR LINEAR EVALUATION.  TAKE AVERAGE OF LEFT HAND
C         AND RIGHT HAND SLOPE AS DEFINED BY THE STRAIGHT LINE SEGMENTS
C         ON THE LEFT AND RIGHT OF THE CURRENT POINT.  IF ONLY ONE SEGMENT
C         EXISTS DEFINE THE SLOPE AT THE POINT WITH THE SLOPE OF THAT
C         LINE SEGMENT.
 
          IF(K.EQ.ISTART) THEN
C           NO LEFT HAND SEGMENT EXISTS.
            SINU(K,JCOL) = (STL(K+1,JCOL) - STL(K,JCOL))/
     A                     (STL(K+1,JAXIS) - STL(K,JAXIS))
          ELSEIF(K.EQ.IEND) THEN
C           NO RIGHT HAND SEGMENT EXISTS
            SINU(K,JCOL) = (STL(K,JCOL) - STL(K-1,JCOL))/
     A                     (STL(K,JAXIS) - STL(K-1,JAXIS))
          ELSE
            SINU(K,JCOL) = 0.5*((STL(K+1,JCOL) - STL(K,JCOL))/
     A                     (STL(K+1,JAXIS) - STL(K,JAXIS)) +
     B                     (STL(K,JCOL) - STL(K-1,JCOL))/
     C                     (STL(K,JAXIS) - STL(K-1,JAXIS)))
          ENDIF
        ELSEIF(SINU(K,JCOL).EQ.-2.E28) THEN
C         REQUEST FOR USING LINE SEGMENT ON THE RIGHT.
          IF(K.EQ.IEND) THEN
            WRITE(STDOUT,54) FLNTAB(JCOL), STL(K,JCOL)
            EFLAG = 1
          ELSE
            SINU(K,JCOL) = (STL(K+1,JCOL) - STL(K,JCOL))/
     A                     (STL(K+1,JAXIS) - STL(K,JAXIS))
          ENDIF
        ELSEIF(SINU(K,JCOL).EQ.-3.E28) THEN
C         REQUEST FOR USING LINE SEGMENT ON THE LEFT
          IF(K.EQ.ISTART) THEN
            WRITE(STDOUT,52) FLNTAB(JCOL), STL(K,JCOL)
            EFLAG = 1
          ELSE
            SINU(K,JCOL) = (STL(K,JCOL) - STL(K-1,JCOL))/
     A                     (STL(K,JAXIS) - STL(K-1,JAXIS))
          ENDIF
        ELSEIF(SINU(K,JCOL).EQ.-4.E28) THEN
C         REQUEST TO USE A PARABOLA FITTED TO THREE POINTS TO DEFINE
C         THE SINUOSITY.
          IF(IEND.LT.ISTART+2) THEN
C           TOO FEW POINTS FOR A PARABOLA
            WRITE(STDOUT,56) FLNTAB(JCOL), STL(K,JCOL)
            EFLAG = 1
          ELSE
C           We have three consecutive points: x0, x0+h1, x0+h1+h2.
C           The points used depend on the location of the point
C           at which we wish to estimate the sinuousity.
            IF(K.EQ.ISTART) THEN
C             WE ARE AT X0.
              H1 = STL(K+1,JAXIS) - STL(K,JAXIS)
              H2 = STL(K+2,JAXIS) - STL(K+1,JAXIS)
              HSUM = H1 + H2
              SINU(K,JCOL) =-(2.*H1 + H2)*STL(K,JCOL)/(H1*HSUM)
     A                     + HSUM*STL(K+1,JCOL)/(H1*H2)
     B                     - H1*STL(K+2,JCOL)/(HSUM*H2)
            ELSEIF(K.EQ.IEND) THEN
C             WE ARE AT X0+H1+H2
              H1 = STL(K-1,JAXIS) - STL(K-2,JAXIS)
              H2 = STL(K,JAXIS) - STL(K-1,JAXIS)
              HSUM = H1 + H2
              SINU(K,JCOL) = +H2*STL(K-2,JCOL)/(H1*HSUM)
     A                      -HSUM*STL(K-1,JCOL)/(H1*H2)
     B                      +(2.*H2 + H1)*STL(K,JCOL)/(HSUM*H2)
            ELSE
C             WE ARE AT X0+H1
              H1 = STL(K,JAXIS) - STL(K-1,JAXIS)
              H2 = STL(K+1,JAXIS) - STL(K,JAXIS)
              HSUM = H1 + H2
              SINU(K,JCOL) = -H2*STL(K-1,JCOL)/(H1*HSUM)
     A                      +(H2 - H1)*STL(K,JCOL)/(H1*H2)
     B                      + H1*STL(K+1,JCOL)/(HSUM*H2)
            ENDIF
            IF(SINU(K,JCOL).LE.0.0) THEN
              WRITE(STDOUT,60) STL(K,JAXIS)
              EFLAG = 1
            ENDIF
          ENDIF
        ELSEIF(SINU(K,JCOL).EQ.-5.E28) THEN
C         REQUEST FOR CUBIC SPLINE.  SET THE CUBIC SPLINE REQUEST FLAG
          CUBIC = 1
        ELSEIF(SINU(K,JCOL).LT.-1.E28) THEN
          WRITE(STDOUT,58) FLNTAB(JCOL), STL(K,JCOL), SINU(K,JCOL)
          STOP 'Abnormal stop. Errors found.'
        ENDIF
 1000 CONTINUE
 
 
      IF(CUBIC.EQ.1) THEN
C       ONE OR MORE REQUESTS FOR A CUBIC SPLINE SLOPE IN THIS
C       RUN OF STATIONS.  RESCAN THE RUN AND COMPUTE ONE OR MORE
C       CUBIC SPLINES.  THERE MAY BE ONE OR MORE SINUOUSITIES
C       SPECIFIED INTERAL TO THE RUN.  THUS WE MUST FIND THE
C       BOUNDARY POINTS AND CONDITIONS FOR EACH SPLINE.
C         X GIVES THE STATIONS ALONG THE AXIS
C         Y GIVES THE STATIONS ALONG THE FLOW LINE AT JCOL
C         YDOT GIVES THE DERIVATIVE OF Y WITH RESPECT TO X.
C              THIS IS ALSO THE SINUOSITY.
 
C       CLEAR THE SPLINE FLAG
        SPFLAG = 0
        DO 2000 K=ISTART,IEND
          IF(SINU(K,JCOL).EQ.-5.E28) THEN
C           SPLINE REQUEST.  IS ONE IN PROGRESS?
            IF(SPFLAG.EQ.0) THEN
C             START A SPLINE.
              SPFLAG = 1
 
C             CLEAR THE VALUE OF THE COUNTER FOR POINTS ON
C             THE SPLINE
              N = 0
 
C             CHECK AND SET THE END CONDITION ON THE LEFT.
              IF(K.EQ.ISTART) THEN
C               NO SINUOSITY GIVEN.  USE DEFAULT END CONDITION
C               OF ZERO SECOND DERIVATIVE.
                LCODE = 2
                LVAL = 0.D0
 
C               TRANSFER STATIONS TO WORKSPACE
                N = N + 1
                X(N) = STL(K,JAXIS)
                Y(N) = STL(K,JCOL)
 
C               SET THE STARTING INDEX FOR THE VALUES OF SINUOSITY
C               BEING DEFINED.
 
                KSTART = K
              ELSE
C               START OF SPLINE NOT AT START OF STATION RUN.
C               THE END CONDITION EXISTS AT K-1 AND IS A KNOWN
C               VALUE OF SINUOSITY
                LCODE = 1
                LVAL = SINU(K-1,JCOL)
 
 
C               SET THE STARTING INDEX FOR THE VALUES OF SINUOSITY
C               BEING DEFINED
                KSTART = K - 1
 
C               TRANSFER TWO STATION VALUES TO THE WORKSPACE
                N = N + 1
                X(N) = STL(K-1,JAXIS)
                Y(N) = STL(K-1,JCOL)
 
                N = N + 1
                X(N) = STL(K,JAXIS)
                Y(N) = STL(K,JCOL)
              ENDIF
            ELSE
C             A SPLINE IS ALREADY IN PROGRESS AND WE HAVE ENCOUNTERED
C             ANOTHER SPLINE REQUEST.  CHECK IF WE HAVE REACHED
C             THE END OF THE STATION RUN.
              IF(K.EQ.IEND) THEN
C               AT THE END OF THE RUN.  NO END CONDITION GIVEN BY THE
C               USER.  USE THE ZERO SECOND DERIVATIVE END CONDITION.
                RCODE = 2
                RVAL = 0.D0
 
                N = N + 1
                X(N) = STL(K,JAXIS)
                Y(N) = STL(K,JCOL)
 
C               SET THE INDEX FOR THE LAST SINUOSITY BEING DEFINED
 
                KEND = K
 
C               COMPUTE THE SPLINE AND STORE THE RESULTS IN SINU(*,*)
 
                CALL SPLINE
     I                     (STDOUT, X, Y, N, LCODE, LVAL, RCODE, RVAL,
     O                      YDOT)
 
                DO 1100 J=KSTART,KEND
                  IF(YDOT(J-KSTART+1).LE.0.0) THEN
                    WRITE(STDOUT,60) STL(J,JAXIS)
                    EFLAG = 1
                  ENDIF
                  SINU(J,JCOL) = YDOT(J-KSTART+1)
 1100           CONTINUE
 
C               CLEAR THE SPLINE REQUEST FLAG
 
                SPFLAG = 0
              ELSE
C               NOT AT END AND SPLINE IS IN PROGRESS.  ADD THE POINTS
C               TO THE SPLINE DEFINITION.
                N = N + 1
                X(N) = STL(K,JAXIS)
                Y(N) = STL(K,JCOL)
              ENDIF
 
            ENDIF
          ELSE
C           KNOWN VALUE OF SINUOSITY ENCOUNTERED.  MAY GIVE AN END
C           CONDITION FOR A SPLINE.
            IF(SPFLAG.EQ.1) THEN
C             GIVES END CONDITION FOR THE CURRENT SPLINE.
 
              RCODE = 1
              RVAL = SINU(K,JCOL)
 
 
              N = N + 1
              X(N) = STL(K,JAXIS)
              Y(N) = STL(K,JCOL)
 
 
C             SET THE INDEX FOR THE LAST SINUOSITY VALUE BEING DEFINED.
 
              KEND = K
 
C             COMPUTE THE SPLINE AND STORE THE RESULTS IN SINU(*,*)
 
              CALL SPLINE
     I                   (STDOUT, X, Y, N, LCODE, LVAL, RCODE, RVAL,
     O                    YDOT)
 
              DO 1200 J=KSTART,KEND
                IF(YDOT(J-KSTART+1).LE.0.0) THEN
                  WRITE(STDOUT,60) STL(J,JAXIS)
                  EFLAG = 1
                ENDIF
                SINU(J,JCOL) = YDOT(J-KSTART+1)
 1200         CONTINUE
 
C             CLEAR THE SPLINE REQUEST FLAG
 
              SPFLAG = 0
            ENDIF
          ENDIF
 2000   CONTINUE
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   STBIN
     I                  (STDIN, STDOUT,
     O                   STL, OFFSET, NUMOFF, SINU, VARTYP, NUMSEC, DIR,
     O                   EPS, JAXIS, EFLAG)
 
C     + + + PURPOSE + + +
C     Input the sinuousity definition table, echo results, and
C     compute missing sinuousity values.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER DIR, EFLAG, JAXIS, NUMSEC, STDIN, STDOUT, VARTYP
      INTEGER NUMOFF(PMXSEC)
      REAL EPS, OFFSET(PMXSEC,PMXNFL), SINU(PMXSEC, PMXNFL),
     A     STL(PMXSEC, PMXNFL)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STL    - Table for flow line stations
C     OFFSET - Offsets for the sinuousity values
C     NUMOFF - Number of offsets
C     SINU   - Table of sinuousities
C     VARTYP - If VARTYP=1 then piecewise linear variation of sinuousity
C               with offset in a cross section is assumed; else
C               if VARTYP=2 then piecewise constant variation of
C               sinuousity with offset in a cross section is assumed.
C     NUMSEC - Number of cross sections in sequence
C     DIR    - If DIR > 0 then stations are ascending order, else
C               descending order
C     EPS    - Tolerance for matching stations
C     JAXIS  - Column in which the values for the channel axis are
C               stored
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER NOPT
      PARAMETER(NOPT=22)
 
C     + + + SAVED VALUES + + +
      INTEGER OPTVAL(NOPT)
      CHARACTER OPTTAB(NOPT)*8
      SAVE OPTTAB, OPTVAL
 
C     + + + LOCAL VARIABLES + + +
      INTEGER EFLAG2, EFLAG3, HEDFLG, I, IFLNAM, ISDEF, IT, J, K, MVAL,
     A        MXNSIN, NFL, NFLNAM, NHEAD, NOFF, NSIN, OPT
      INTEGER CLEN(PMXNFL), DIFF(PMXNFL), HDTOCL(PMXNFL), IVAL(PMXNFL),
     A        MAXNFL(PMXSEC), TERMCLS(PMXNFL), TERML(PMXNFL), 
     B        VTYPE(PMXNFL)
      REAL DELTA, MINDEL, OLDOFF, RVAL(PMXNFL)
      DOUBLE PRECISION DPVAL(PMXNFL)
      CHARACTER CRDNAM*4, CVAL(PMXNFL)*256, FLNTAB(PMXNFL)*8, KEY*8,
     A          LINE*196,  TERM(PMXNFL)*1, SINDEF*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL ADJSIN, BINSER, CPSINU, GETVAL, inline, LSATAB, LSTAB,
     A         STBOUT
 
C     + + + DATA INITIALIZATIONS + + +
      DATA OPTTAB/'CSPLINE','CUBIC','LIN','LIND','LINEAR','LINEARD',
     A            'LINEARU','LINU','PARAB','PARABOLA','PARABOLI',
     B            'cspline','cubic','lin','lind','linear','lineard',
     C            'linearu','linu','parab','parabola','paraboli'/
      DATA OPTVAL/5,5,1,2,1,2,3,3,4,4,4,5,5,1,2,1,2,3,3,4,4,4/
 
C     + + + INPUT FORMATS + + +
 2    FORMAT(7X,A8)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:654* No values found on input line.')
 52   FORMAT(/,' *ERR:655* No heading values available. Lines out of',
     A         ' order.')
 54   FORMAT(/,' *ERR:656* Number of sections > ',I5,' Increase PMXSEC',
     A       ' and recompile.')
 56   FORMAT(/,' *ERR:657* LENG specified but no initial station',
     A       ' found.',/,10X,'Specify an initial station before the',
     B       ' first LENG line.')
 58   FORMAT(/,' *ERR:658* Initial station missing for LENG field',
     A       ' containing:',F10.4)
 60   FORMAT(/,' *ERR:659* Station missing for offset:',F10.2)
 62   FORMAT(/,' *ERR:660* Number of sinuosities:',I5,' does not match',
     A      ' number of stations:',I5,/,10X,'Input out of order?')
 64   FORMAT(/,' *ERR:661* Blank heading in field no.:',I3,' invalid.')
 66   FORMAT(/,' *ERR:662* Station missing for sinuosity=',F10.4)
 68   FORMAT(/,' *ERR:663* Number of offsets=',I5,
     A   ' incompatible with',/,10X,' number of flow lines=',I5)
 70   FORMAT(/,' *ERR:671* Invalid name for an input line:',A4, /,
     A 10X,'  Valid names are:STAT, LENG, OFFS, SINU, HEAD, and END.')
 72   FORMAT(/,' *ERR:664* Label for 1-D axis not found.',
     A      '  Must be AXIS, axis, or Axis.')
 74   FORMAT(/,' *ERR:665* LENG value number',I5,' is undefined.')
 76   FORMAT(/,' *ERR:666* STAT value number',I5,' is undefined.')
 78   FORMAT(/,' *ERR:667* OFFS value number',I5,' is undefined.')
 80   FORMAT(/,' *ERR:668* Sinuousity definition:',A8,' unknown.',
     A       ' Assuming CUBIC.')
 82   FORMAT(/,' *ERR:669* Sinuosity option:',A8,' unknown.',
     A      '  Assuming LINEAR.')
 84   FORMAT(/,' *BUG:XXX* Invalid option value=',I5,' for:',A8,
     A        ' in STBIN')
 86   FORMAT(/,' *ERR:670* There are:',I5,' values on the line',
     A     ' but only',I5,' heading values')
 88   FORMAT(/,' Sinuousity table as defined by input.')
 90   FORMAT(/,' Final values of sinuosity.')
 92   FORMAT(/,' *ERR:672 Change in sinuosity variation invalid.',
     A     /,10X,'Check number of offsets or stations.')
 94   FORMAT(/,' *ERR:673* Number of flow lines=',I5,' too small. Must',
     A      ' be at least two flow lines.')
 96   FORMAT(/,' *ERR:674* Change in direction of stationing at',
     A       ' station:',F10.4,' invalid.')
 97   FORMAT(/,' *ERR:675* Station match at station:',F10.4,' invalid.')
 98   FORMAT(/,' *ERR:676* Inconsistent stationing.  Flow line',
     A       ' stations must all',10X,'increase or must all decrease.')
 99   FORMAT(/,' *ERR:677* Offset=',F10.2,' <= previous offset=',F10.2)
C***********************************************************************
C     GET THE SINUOSITY VALUE: LINEAR, PARABOLIC, OR CUBIC.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,'(1X,A80)') LINE
 
      READ(LINE,2) SINDEF
 
      CALL BINSER
     I           (SINDEF, NOPT, OPTTAB,
     O            I)
      IF(I.EQ.0) THEN
        EFLAG = 1
        WRITE(STDOUT,80) SINDEF
        ISDEF = 5
      ELSE
        ISDEF = OPTVAL(I)
      ENDIF
 
C     INITIALIZE THE NUMBER OF FLOW LINE NAMES CURRENTLY KNOWN
      NFLNAM = 0
 
C     CLEAR THE HEADING FLAG.
      HEDFLG = 0
 
C     CLEAR THE COUNTER FOR THE NUMBER OF SECTIONS
      NUMSEC = 0
 
C     CLEAR THE COUNTER FOR THE NUMBER OF HEADING VALUES
      NHEAD = 0
 
C     CLEAR THE LOCAL ERROR FLAGS
      EFLAG2 = 0
      EFLAG3 = 0
 
C     SET THE MAXIMUM NUMBER OF USER FLOW LINES
      MXNSIN = PMXNFL - 2
 
C     SET THE SINUOSITY VARIATION TO UNDEFINED
      VARTYP = -1
 
C     SET THE MINIMUM STATION INTERVAL TO INITIAL VALUE
      MINDEL = 1.E30
 
C     INITIALIZE THE ARRAYS TO DEFAULT VALUES
      DO 110 I=1,PMXSEC
        MAXNFL(I) = 0
        NUMOFF(I) = 0
        DO 100 J=1,PMXNFL
          STL(I,J) = -1.E30
          OFFSET(I,J) = -1.E30
          GOTO(101,102,103,104,105), ISDEF
            WRITE(STDOUT,84) I, SINDEF
            STOP 'Abnormal stop. Errors found.'
 101      CONTINUE
 102      CONTINUE
 103      CONTINUE
C           TAKE AS LINEAR
            SINU(I,J) = -1.E28
            GOTO 109
 104      CONTINUE
C           TAKE AS PARABOLIC
            SINU(I,J) = -4.E28
            GOTO 109
 105      CONTINUE
C           TAKE AS CUBIC
            SINU(I,J) = -5.E28
            GOTO 109
 109      CONTINUE
 100    CONTINUE
 110  CONTINUE
 
C     Clear the station difference flag and the flow line name
C     table.
      DO 112 J=1,PMXNFL
        DIFF(J) = 0.0
        FLNTAB(J) = ' '
 112  CONTINUE
 
 9999 CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        WRITE(STDOUT,'(1X,A80)') LINE
 
C       USE THE FIRST FOUR CHARS OF LINE TO DEFINE THE NATURE OF
C       THE REMAINDER OF THE LINE.
 
        CRDNAM = LINE(1:4)
 
        IF(CRDNAM.EQ.'    ') THEN
C         SKIP BLANK LINES
          GOTO 9999
        ENDIF
 
        IF(CRDNAM.NE.'END'.AND.CRDNAM.NE.'end'.AND.
     A     CRDNAM.NE.'HEAD'.AND.CRDNAM.NE.'head') THEN
 
C         ALL OTHER OPTIONS HAVE ALL VALUES REAL.
 
          IF(CRDNAM.EQ.'SINU'.OR.CRDNAM.EQ.'sinu') THEN
            OPT = 1
            DO 201 I=1,PMXNFL
              RVAL(I) = -1.E30
 201        CONTINUE
          ELSE
            OPT = 0
 
            DO 200 I=1,PMXNFL
              VTYPE(I) = 2
              RVAL(I) = -1.E30
 200        CONTINUE
          ENDIF
 
C          CALL GETVAL
C     I               (STDOUT, LINE(5:80), MXNSIN, OPT,
C     O                VTYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, EFLAG2,
C     O                MVAL)
 
           CALL GETVAL
     I                (STDOUT, LINE(5:), MXNSIN, OPT,
     O                 VTYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O                 EFLAG2, TERM, TERML, TERMCLS,
     O                 MVAL)

          IF(EFLAG2.NE.0) THEN
C           SOME ERROR HAS OCCURRED.  SKIP TO THE NEXT LINE.
            EFLAG3 = 1
            EFLAG2 = 0
          ELSEIF(MVAL.EQ.0) THEN
C           NO VALUES FOUND.
            WRITE(STDOUT,50)
            EFLAG3 = 1
          ELSEIF(HEDFLG.EQ.0) THEN
C           UNABLE TO PROCESS FURTHER.  HEADING LABELS ARE UNDEFINED.
            WRITE(STDOUT,52)
            EFLAG3 = 1
          ELSEIF(MVAL.GT.NHEAD) THEN
C           MORE VALUES ON A LINE THAN WE HAVE DEFINED HEADINGS FOR
C           THE LINE
            WRITE(STDOUT,86) MVAL, NHEAD
            EFLAG3 = 1
          ENDIF
 
          IF(EFLAG3.EQ.0) THEN
 
            IF(CRDNAM.EQ.'STAT'.OR.CRDNAM.EQ.'stat') THEN
C             INPUT OF STATION DATA.
              NUMSEC = NUMSEC + 1
              IF(NUMSEC.GT.PMXSEC) THEN
                WRITE(STDOUT,54) PMXSEC
                EFLAG = 1
                NUMSEC = PMXSEC
              ENDIF
              NFL = 0
              DO 210 J=1,MVAL
                IF(RVAL(J).GT.-1.E29) THEN
C                 TAKE AS A VALID VALUE GIVEN BY THE USER.
                  IT = HDTOCL(J)
                  STL(NUMSEC,IT) = RVAL(J)
                  NFL = NFL + 1
                  IF(NUMSEC.GT.1) THEN
C                   CHECK FOR DIRECTION OF STATIONING
                    IF(STL(NUMSEC-1,IT).GT.-1.E29) THEN
C                     PRECEDING STATION IS DEFINED.
                      DELTA = STL(NUMSEC,IT) - STL(NUMSEC-1,IT)
                      MINDEL = MIN(MINDEL,ABS(DELTA))
                      IF(DELTA.GT.0.0) THEN
C                       STATIONS ARE ASCENDING
                        IF(DIFF(IT).EQ.0) THEN
C                         FIRST STATION INCREMENT FOR THIS FLOWLINE
                          DIFF(IT) = 1
                        ELSEIF(DIFF(IT).EQ.-1) THEN
C                         ERROR.  CHANGE IN DIRECTION OF STATIONING
                          WRITE(STDOUT,96) RVAL(J)
                          EFLAG = 1
                        ENDIF
                      ELSEIF(DELTA.LT.0.0) THEN
C                       STATIONS ARE DESCENDING
                        IF(DIFF(IT).EQ.0.0) THEN
C                         FIRST STATION INCREMENT FOR THIS FLOW LINE
                          DIFF(IT) = -1
                        ELSEIF(DIFF(IT).EQ.1) THEN
                          WRITE(STDOUT,96) RVAL(J)
                        ENDIF
                      ELSE
                        WRITE(STDOUT,97) RVAL(J)
                        EFLAG = 1
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  WRITE(STDOUT,76) J
                  EFLAG = 1
                ENDIF
 210          CONTINUE
 
              MAXNFL(NUMSEC) = MVAL
              IF(MVAL.LT.2) THEN
                WRITE(STDOUT,94) MVAL
                EFLAG = 1
              ENDIF
            ELSEIF(CRDNAM.EQ.'LENG'.OR.CRDNAM.EQ.'leng') THEN
C             FOUND A LINE GIVING THE FLOW LINE LENGTHS.
              NUMSEC = NUMSEC + 1
              IF(NUMSEC.GT.PMXSEC) THEN
                WRITE(STDOUT,54) PMXSEC
                EFLAG = 1
                NUMSEC = PMXSEC
              ENDIF
 
              IF(NUMSEC.LT.1.) THEN
C               NO STATION VALUE GIVEN.  ERROR: IMPROPER INITIALIZATION
                WRITE(STDOUT,56)
                EFLAG = 1
              ELSE
 
                NFL = 0
                DO 310 J=1,MVAL
                  IF(RVAL(J).GT.-1.E29) THEN
C                   TAKE AS A VALID VALUE GIVEN BY THE USER.
                    K = HDTOCL(J)
                    IF(STL(NUMSEC-1,K).LT.-1.E29) THEN
                      WRITE(STDOUT,58) RVAL(J)
                      EFLAG = 1
                      STL(NUMSEC-1,K) = 0.0
                    ENDIF
C                   ADD THE CURRENT LENGTH TO THE PREVIOUS STATION AND
C                   STORE AS A STATION VALUE.
                    STL(NUMSEC,K) = STL(NUMSEC-1,K) + RVAL(J)
                    NFL = NFL + 1
 
 
                    IF(NUMSEC.GT.1) THEN
C                     CHECK FOR DIRECTION OF STATIONING
                      IF(STL(NUMSEC-1,IT).GT.-1.E29) THEN
C                       PRECEDING STATION IS DEFINED.
                        DELTA = STL(NUMSEC,IT) - STL(NUMSEC-1,IT)
                        MINDEL = MIN(MINDEL,ABS(DELTA))
                        IF(DELTA.GT.0.0) THEN
C                         STATIONS ARE ASCENDING
                          IF(DIFF(IT).EQ.0) THEN
C                           FIRST STATION INCREMENT FOR THIS FLOWLINE
                            DIFF(IT) = 1
                          ELSEIF(DIFF(IT).EQ.-1) THEN
C                           ERROR.  CHANGE IN DIRECTION OF STATIONING
                            WRITE(STDOUT,96) RVAL(J)
                            EFLAG = 1
                          ENDIF
                        ELSEIF(DELTA.LT.0.0) THEN
C                         STATIONS ARE DESCENDING
                          IF(DIFF(IT).EQ.0.0) THEN
C                           FIRST STATION INCREMENT FOR THIS FLOW LINE
                            DIFF(IT) = -1
                          ELSEIF(DIFF(IT).EQ.1) THEN
                            WRITE(STDOUT,96) RVAL(J)
                          ENDIF
                        ELSE
                          WRITE(STDOUT,97) RVAL(J)
                          EFLAG = 1
                        ENDIF
                      ENDIF
                    ENDIF
 
                  ELSE
                    WRITE(STDOUT,74) J
                    EFLAG = 1
                  ENDIF
 310            CONTINUE
 
                MAXNFL(NUMSEC) = MVAL
                IF(MVAL.LT.2) THEN
                  WRITE(STDOUT,94) MVAL
                  EFLAG = 1
                ENDIF
              ENDIF
            ELSEIF(CRDNAM.EQ.'OFFS'.OR.CRDNAM.EQ.'offs') THEN
C             WE HAVE OFFSET VALUES GIVEN.
              NOFF = 0
              OLDOFF = -1.E30
              DO 410 J=1,MVAL
                IF(RVAL(J).GT.-1.E29) THEN
C                 TAKE AS A VALID VALUE GIVEN BY THE USER.
                  K = HDTOCL(J)
                  IF(STL(NUMSEC,K).LT.-1.E29) THEN
                    WRITE(STDOUT,60) RVAL(J)
                    EFLAG = 1
                  ENDIF
                  OFFSET(NUMSEC,K) = RVAL(J)
                  IF(RVAL(J).LE.OLDOFF) THEN
                    WRITE(STDOUT,99) RVAL(J), OLDOFF
                    EFLAG = 1
                  ENDIF
                  OLDOFF = RVAL(J)
                  NOFF = NOFF + 1
                ELSE
                  WRITE(STDOUT,78) J
                  EFLAG = 1
                ENDIF
 410          CONTINUE
 
              NUMOFF(NUMSEC) = NOFF
              IF(NOFF.GT.NFL.OR.NOFF.LT.NFL-1) THEN
C               ERROR
                WRITE(STDOUT,68) NOFF, NFL
                EFLAG = 1
              ENDIF
 
              IF(NFL.EQ.NOFF) THEN
                IF(VARTYP.EQ.-1) THEN
C                 PIECEWISE LINEAR VARIATION
                  VARTYP = 1
                ELSEIF(VARTYP.EQ.2) THEN
                  WRITE(STDOUT,92)
                  EFLAG = 1
                ENDIF
              ELSE
                IF(VARTYP.EQ.-1) THEN
C                 PIECEWISE CONSTANT VARIATION
                  VARTYP = 2
                ELSEIF(VARTYP.EQ.1) THEN
                  WRITE(STDOUT,92)
                  EFLAG = 1
                ENDIF
              ENDIF
            ELSEIF(CRDNAM.EQ.'SINU'.OR.CRDNAM.EQ.'sinu') THEN
C             A LINE OF SINUOSITY VALUES HAS BEEN FOUND.
 
              NSIN = 0
              DO 610 J=1,MVAL
 
C               DO THE CONVERSION BASED ON THE TYPE FOUND.
                IF(VTYPE(J).EQ.4) THEN
C                 CHARACTER VALUE RETURNED.  SEE IF WE KNOW IT!
                  CALL BINSER
     I                       (CVAL(J)(1:8), NOPT, OPTTAB,
     O                        I)
                  IF(I.EQ.0) THEN
                    WRITE(STDOUT,82) CVAL(J)(1:8)
                    RVAL(J) = -1.E28
                    GOTO 609
                  ENDIF
                  I = OPTVAL(I)
                  GOTO(601,602,603,604,605), I
                    WRITE(STDOUT,84) I, CVAL(J)(1:8)
                    STOP 'Abnormal stop. Errors found.'
 601              CONTINUE
C                   LINEAR OPTION FOR SINUOSITY AT A POINT
                    RVAL(J) = -1.E28
                    GOTO 609
 602              CONTINUE
C                   LINEAR DOWNSTREAM OPTION FOR SINUOSITY AT A POINT
                    RVAL(J) = -2.E28
                    GOTO 609
 603              CONTINUE
C                   LINEAR UPSTREAM OPTION FOR SINUOSITY AT A POINT
                    RVAL(J) = -3.E28
                    GOTO 609
 604              CONTINUE
C                   PARABOLIC OPTION FOR SINUOSITY AT A POINT
                    RVAL(J) = -4.E28
                    GOTO 609
 605              CONTINUE
C                   CUBIC SPLINE OPTION FOR SINUOSITY AT A POINT
                    RVAL(J) = -5.E28
                    GOTO 609
 609              CONTINUE
                ELSEIF(VTYPE(J).EQ.-1) THEN
C                 ASTERISK OR DUAL COMMAS
                  RVAL(J) = -1.E30
                  NSIN = NSIN + 1
                ELSE
C                 SHOULD BE A NUMBER
                  READ(CVAL(J)(1:15),'(F15.0)') RVAL(J)
                ENDIF
                IF(RVAL(J).GT.-1.E29) THEN
C                 TAKE AS A VALID VALUE GIVEN BY THE USER.
                  K = HDTOCL(J)
                  IF(STL(NUMSEC,K).LT.-1.E29) THEN
                    WRITE(STDOUT,66) RVAL(J)
                    EFLAG = 1
                  ENDIF
                  SINU(NUMSEC,K) = RVAL(J)
                  NSIN = NSIN + 1
                ENDIF
 610          CONTINUE
 
 
              IF(NSIN.NE.NFL) THEN
C               ERROR
                WRITE(STDOUT,62) NSIN, NFL
                EFLAG = 1
              ENDIF
            ENDIF
          ELSE
            EFLAG = 1
            EFLAG3 = 0
          ENDIF
        ELSE
 
          IF(CRDNAM.EQ.'END'.OR.CRDNAM.EQ.'end') THEN
C           INPUT OF SINUOSITY DEFINTION COMPLETED.
            GOTO 10000
          ELSEIF(CRDNAM.EQ.'HEAD'.OR.CRDNAM.EQ.'head') THEN
C           A HEADING LINE HAS BEEN FOUND.  SET TO IDENTIFIERS
 
            OPT = 0
            DO 500 I=1,PMXNFL
              VTYPE(I) = 4
              CVAL(I) = ' '
 500        CONTINUE
C            CALL GETVAL
C     I                 (STDOUT, LINE(5:80), MXNSIN, OPT,
C     O                  VTYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, EFLAG2,
C     O                  MVAL)
           CALL GETVAL
     I                (STDOUT, LINE(5:), MXNSIN, OPT,
     O                 VTYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O                 EFLAG2, TERM, TERML, TERMCLS,
     O                 MVAL)
 
            IF(EFLAG2.NE.0) THEN
C             SOME ERROR HAS OCCURRED.  SKIP TO THE NEXT LINE.
              EFLAG = 1
              EFLAG2 = 0
            ELSEIF(MVAL.EQ.0) THEN
C             NO VALUES FOUND.
              WRITE(STDOUT,50)
              EFLAG = 1
            ELSE
 
C             HERE WE HAVE NO DETECTED ERRORS AND AT LEAST 1 AND NO
C             MORE THAN PMXNFL NON-DEFAULT CHARACTER VALUES IN
C             CVAL(*).
 
              DO 510 J=1,MVAL
                IF(CVAL(J).EQ.' ') THEN
C                 ERROR.  HEADING NAMES MUST BE EXPLICIT-NO DEFAULT
C                 NAMES
                  WRITE(STDOUT,64) J
                  EFLAG = 1
                ENDIF
                KEY = CVAL(J)(1:8)
C               FIND KEY IN TABLE. ADD TO TABLE IF NOT FOUND.
 
                CALL LSATAB
     I                     (STDOUT, KEY, PMXNFL,
     M                      FLNTAB, NFLNAM,
     O                      IFLNAM, EFLAG)
                HDTOCL(J) = IFLNAM
 510          CONTINUE
 
C             SET THE HEADING FLAG
              HEDFLG = 1
 
C             SET THE CURRENT NUMBER OF HEADINGS AVAILABLE
              NHEAD = MVAL
            ENDIF
          ELSE
C           INVALID NAME FOR A LINE
            WRITE(STDOUT,70) CRDNAM
            EFLAG = 1
          ENDIF
        ENDIF
      GOTO 9999
 
 
10000 CONTINUE
 
C     FIND THE AXIS LOCATION
 
      CALL LSTAB
     I          ('AXIS    ', FLNTAB, NFLNAM,
     O           JAXIS)
      IF(JAXIS.EQ.0) THEN
C       TRY ANOTHER SPELLING.
        CALL LSTAB
     I            ('axis    ', FLNTAB, NFLNAM,
     O             JAXIS)
        IF(JAXIS.EQ.0) THEN
C         TRY YET ANOTHER SPELLING
          CALL LSTAB
     I              ('Axis    ', FLNTAB, NFLNAM,
     O               JAXIS)
          IF(JAXIS.EQ.0) THEN
C           GIVE UP AND CALL IT A USER ERROR.
            WRITE(STDOUT,72)
            EFLAG = 1
          ENDIF
        ENDIF
      ENDIF
 
      IF(EFLAG.GT.0) RETURN
 
C     FORCE THE SINUOSITY ON THE AXIS TO BE 1.0
      DO 11000 I=1,NUMSEC
        SINU(I,JAXIS) = 1.0
11000 CONTINUE
 
C     OUTPUT THE PART DEFINED BY INPUT
 
      IF(EFLAG.EQ.0) THEN
        WRITE(STDOUT,88)
        CALL STBOUT
     I             (STDOUT, NUMSEC, STL, OFFSET, SINU, NFLNAM, FLNTAB)
 
      ENDIF
 
C     CHECK DIRECTION OF STATIONING AND SET DIRECTION
      DIR = DIFF(1)
      DO 11010 J=2,NFLNAM
        IF(DIFF(J).NE.DIR) THEN
          WRITE(STDOUT,98)
          EFLAG = 1
        ENDIF
11010 CONTINUE
 
 
C     COMPUTE THE SINUOSITIES NOT DEFINED BY THE USER INPUT
 
      IF(EFLAG.EQ.0) THEN
        CALL CPSINU
     I             (STDOUT, NUMSEC, JAXIS, FLNTAB, STL, NFLNAM,
     M              SINU,
     O              EFLAG)
 
C       OUTPUT THE FINAL SINUOSITY VALUES
 
        WRITE(STDOUT,90)
        CALL STBOUT
     I             (STDOUT, NUMSEC, STL, OFFSET, SINU, NFLNAM, FLNTAB)
 
        IF(EFLAG.EQ.0) THEN
 
C         EXTEND THE SINUOSITY DEFINITION IN CONCORDANCE WITH THE
C         TYPE OF VARIATION.  WE ADD EXTRA POINTS
C         TO MAKE LATER LOOKUP EASIER.
 
 
          CALL ADJSIN
     I               (MAXNFL, VARTYP, NUMSEC,
     M                OFFSET, SINU,
     O                NUMOFF)
        ENDIF
       ENDIF
 
C     COMPUTE THE TOLERANCE FOR FINDING STATIONS.
 
      EPS = 0.1*MINDEL
 
      END
