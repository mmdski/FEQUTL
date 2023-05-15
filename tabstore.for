C
C
C
      SUBROUTINE   TWOD13GET
     I                      (STDOUT, 
     O                       TABLE, NHU, NFRAC, HUVEC, FDROP,
     O                       PFDVEC, QMAT, HDATUM, TYPE, 
     O                       EFLAG)
 
C     + + + PURPOSE + + +
C     Get a 2-D table of type 13 from FTAB/ITAB and put it in
C     the computational matrices.

      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG,  NFRAC, NHU, STDOUT, TABLE, TYPE
      REAL FDROP(PMXNHU), HDATUM, HUVEC(PMXNHU), PFDVEC(PMXFRC),
     A     QMAT(PMXNHU,PMXFRC)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     TABLE  - Table number
C     NHU    - Number of positive upstream heads.
C     NFRAC  - Number of fractions for defining partial free drop
C               values including zero free drop.
C     HUVEC  - Vector of upstream heads
C     FDROP  - Free drop values
C     PFDVEC - Partial free drop vector
C     QMAT   - Flowrate with the rows indexed by upstream head
C              (1:NHU) and the columns by partial free drop values
C              (1:NFRAC)
C     HDATUM - Datum for measuring head
C     TYPE   - Table type
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ftable.cmn'

C     
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IP, J, K, PHU, PPFD, QPNT
C     ****************************************************************** 
C***********************************************************************
      IP = TABLE 
      CALL CHKTAB
     I           (6, STDOUT, FTPNT, MFTNUM,
     M            IP,
     O            EFLAG)
      

      HDATUM = FTAB(IP) 
      TABLE = ITAB(IP+1)
      TYPE = ITAB(IP+2)

      WRITE(STDOUT,*) ' TWOD13GET: TABLE=',TABLE,' TYPE=',TYPE
      PHU = ITAB(IP+3) 
 
      PPFD = ITAB(IP+5)

      WRITE(STDOUT,*) ' TWOD13GET: PHU=',PHU,' PPFD=',PPFD
      
C     Compute the number of nonzero upstream heads
      NHU = (PPFD - PHU)/4 - 1
 
C     Compute the number of partial free drops
      NFRAC = ITAB(PHU+3)
 
      WRITE(STDOUT,*) ' TWOD13GET: NHU=',NHU,' NFRAC=',NFRAC
      K = PHU + 4
      DO 730 I=1, NHU
        HUVEC(I)      = FTAB(K)  
        FDROP(I)      = FTAB(K+1)
        QMAT(I,NFRAC) = FTAB(K+2)
        K = K + 4
730   CONTINUE
 
C     Restore the partial free drops. 

      DO 750 I=1, NFRAC
       PFDVEC(I) = FTAB(PPFD+I-1)
750   CONTINUE 
 
C     Restore the flows for the non-zero upstream heads.  
      QPNT = PPFD + NFRAC
      DO 802 I=1, NHU
        IP = QPNT + I*NFRAC
        DO 801 J=1, NFRAC
          QMAT(I,J) = FTAB(IP+J-1)
 801    CONTINUE
 802  CONTINUE


      RETURN
      END
C
C
C
      SUBROUTINE   TWOD13PUT
     I                      (STDOUT, TABLE, NHU, NFRAC, HUVEC, FDROP,
     I                       PFDVEC, QMAT, HDATUM, TYPE, 
     M                       FTP, 
     O                       EFLAG)
 
C     + + + PURPOSE + + +
C     Store a two-D flow table of  type 13 in the FTAB/ITAB
C     structure given the computational matrices.

      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, NFRAC, NHU, STDOUT, TABLE, TYPE
      REAL FDROP(PMXNHU), HDATUM, HUVEC(PMXNHU), PFDVEC(PMXFRC),
     A     QMAT(PMXNHU,PMXFRC)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     TABLE  - Table number
C     NHU    - Number of positive upstream heads.
C     NFRAC  - Number of fractions for defining partial free drop
C               values including zero free drop.
C     HUVEC  - Vector of upstream heads
C     FDROP  - Free drop values
C     PFDVEC - Partial free drop vector
C     QMAT   - Flowrate with the rows indexed by upstream head
C              (1:NHU) and the columns by partial free drop values
C              (1:NFRAC)
C     HDATUM - Datum for measuring head
C     TYPE   - Table type
C     FTP    - Pointer to next available entry in FTAB/ITAB
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ftable.cmn'

C     
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IP, J, K,  N, PHU, PPFD, QPNT
C     ****************************************************************** 
C***********************************************************************
      IF(FTP+10.GE.MRFTAB) CALL KIL
     I                              (10,
     M                               FTP, EFLAG)
 
C     Store the datum for heads for this table. 
C      WRITE(STDOUT,*) ' TABLE=',TABLE
      IF(FTPNT(TABLE).NE.0) THEN
        CALL KIL
     I           (25,
     M            TABLE, EFLAG)
        STOP 'Abnormal stop. Error found.' 
      ENDIF
      FTPNT(TABLE) = FTP
      FTAB(FTP) = HDATUM
      ITAB(FTP+1) = TABLE
      ITAB(FTP+2) = TYPE


C     Store 2-D table of type 13 in FTAB/ITAB.
C     Set pointer values for the upstream head
      PHU = FTP + 10
 
C      WRITE(STDOUT,*) ' PHU=',PHU
C     Store pointer to the start of the downstream heads.
 
      ITAB(FTP+3) = PHU
 
C     Initialize the value of the last access to the free head
C     list to the zero head value. 
      ITAB(FTP+4) = PHU

C     CHECK FOR SPACE
      N = 10 + 4*(NHU + 1) + (NHU + 2)*NFRAC
      IF(FTP+N.GE.MRFTAB) CALL KIL
     I                             (10,
     M                              FTP, EFLAG)
 
C     Find and store the base pointer for partial free drops. 
 
      PPFD = PHU + 4*(NHU + 1)
C      WRITE(STDOUT,*) ' PPFD=',PPFD
      ITAB(FTP+5) = PPFD
C     Initialize the last access pointer into the partial free drop 
C     sequence
      ITAB(FTP+6) = PPFD

C     Offset locations 7 and 8 are currently not used.
 

C     Store the maximum upstream head
      FTAB(FTP+9) = HUVEC(NHU)
 
C     Compute the base pointer for the flow in the table  
      QPNT = PPFD + NFRAC
C      WRITE(STDOUT,*) ' QPNT=',QPNT

C     Store the values for the upstream heads, and the corresponding
C     drop to free flow, the free flow value, and the offset to 
C     the flows for each upstream head.

C     Store the values for zero upstream head
 
      FTAB(PHU) = 0.0
      FTAB(PHU+1) = 0.0
      FTAB(PHU+2) = 0.0
C     Store the offset between the origin of the partial free
C     drop sequence and the origin of the flows at zero upstream
C     head so that a single addition can be used to find the 
C     flows once the location of the partial free drop argument
C     is found in the partial free drop sequence.
 
      ITAB(PHU+3) = QPNT - PPFD
C      WRITE(STDOUT,*) ' OFFSET AT ORIGIN=',ITAB(PHU+3)

      K = PHU + 4
      DO 730 I=1, NHU
        FTAB(K) = HUVEC(I)
        FTAB(K+1) = FDROP(I)
        FTAB(K+2) = QMAT(I,NFRAC)
        ITAB(K+3) = QPNT - PPFD + I*NFRAC
        K = K + 4
730   CONTINUE
 

C     Store the partial free drops. 

      DO 750 I=1, NFRAC
        FTAB(PPFD+I-1) = PFDVEC(I)
750   CONTINUE 
 
C     Store the zero flows for zero upstream head
 
      DO 800 I=QPNT,QPNT + NFRAC - 1
        FTAB(I) = 0.0
 800  CONTINUE
 
C     Store the flows for the non-zero upstream heads.  Q contains
C     the zero flow at zero partial free drop. 
 
      DO 802 I=1, NHU
        IP = QPNT + I*NFRAC
        DO 801 J=1, NFRAC
          FTAB(IP+J-1) = QMAT(I,J)
 801    CONTINUE
 802  CONTINUE

C     Update the pointer into FTAB/ITAB
      IF(FTP + N.NE.IP + NFRAC) THEN
        WRITE(STDOUT,50) FTP + N, IP + NFRAC
50    FORMAT(/,' *BUG:XXX* In TWOD13PUT: FTP+N=',I10,' IP+NFRAC=',I10)
        STOP 'Abnormal stop: Bug found' 
      ENDIF

      FTP= IP + NFRAC

      RETURN
      END
C
C
C
      SUBROUTINE   STBOUT
     I                   (STDOUT, NUMSEC, STL, OFFSET, SINU, NFLNAM,
     I                    FLNTAB)
 
C     + + + PURPOSE + + +
C     Output the sinuousity table.  If all sinuousity values
C     are undefined for a given cross section, do not output any
C     values.  If one or more are defined output all values
C     giving a blank for values that are not defined.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER NFLNAM, NUMSEC, STDOUT
      REAL OFFSET(PMXSEC, PMXNFL), SINU(PMXSEC,PMXNFL),
     A     STL(PMXSEC, PMXNFL)
      CHARACTER FLNTAB(NFLNAM)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NUMSEC - Number of cross sections in sequence
C     STL    - Table for flow line stations
C     OFFSET - Offsets for the sinuousity values
C     SINU   - Table of sinuousities
C     NFLNAM - Number of flow line names
C     FLNTAB - Flow line name table
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, NFL, SNOUT
      CHARACTER CHR5*5, OUT(PMXNFL)*8, STAOUT*7
 
C     + + + INTRINSICS + + +
      INTRINSIC INDEX
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL MKFMT
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT(/,'  Item',9(1X,A8))
 54   FORMAT(1X,' STAT',9(1X,A8))
 56   FORMAT(1X,' OFFS',9(1X,A8))
 58   FORMAT(1X,' SINU',9(1X,A8))
C***********************************************************************
C     SET UP THE STATION OUTPUT FORMAT
      STAOUT(1:1) = '('
      STAOUT(7:7) = ')'
C     OUTPUT THE TABLE LABEL AND THE HEADING FOR THE TABLE
 
      DO 90 J=1,NFLNAM
        I = INDEX(FLNTAB(J),' ') - 1
        IF(I.EQ.-1) I = 8
        OUT(J) = ' '
        OUT(J)(8-I+1:8) = FLNTAB(J)(1:I)
 90   CONTINUE
      WRITE(STDOUT,52) (OUT(J), J=1,NFLNAM)
 
C     FOR EACH SECTION IN THE SINUOSITY TABLE OUTPUT THE STATIONS,
C     THE OFFSETS, AND IF DEFINED, THE SINUOSITY VALUES.  ALL SECTIONS
C     MUST HAVE AS A MINIMUM ONE OR MORE STATIONS AND EITHER EXACTLY
C     AS MANY OFFSETS AS STATIONS OR ONE LESS STATION THAN OFFSETS.
 
      DO 1000 I=1,NUMSEC
 
C       MUST SEARCH ACROSS ALL COLUMNS TO FIND THE FLOW LINES.
 
        NFL = NFLNAM
        DO 100 J=1,NFL
          IF(STL(I,J).LT.-1.E29) THEN
            OUT(J) = ' '
          ELSE
            CALL MKFMT
     I                (STL(I,J), 8,
     O                 CHR5)
            STAOUT(2:6) = CHR5
            WRITE(OUT(J),STAOUT) STL(I,J)
          ENDIF
 100    CONTINUE
 
        IF(I.GT.1) THEN
          WRITE(STDOUT,*) ' '
        ENDIF
 
        WRITE(STDOUT,54) (OUT(J), J=1,NFL)
 
        DO 200 J=1,NFL
          IF(OFFSET(I,J).LT.-1.E29) THEN
            OUT(J) = ' '
          ELSE
            WRITE(OUT(J),'(F8.2)') OFFSET(I,J)
          ENDIF
 200    CONTINUE
        WRITE(STDOUT,56) (OUT(J), J=1,NFL)
 
C       OUTPUT THE SINUOSITY . CLEAR THE SINUOSITY OUTPUT FLAG.
C       ONLY OUTPUT SINUOSITY AT LOCATIONS WHERE THE STATION IS
C       DEFINED.
        SNOUT = 0
        DO 300 J=1,NFL
          OUT(J) = ' '
          IF(STL(I,J).GT.-1.E29) THEN
            IF(SINU(I,J).GT.-1.E29) THEN
C             A DEFINED VALUE. SET THE FLAG AND PLACE IN OUT
              SNOUT = 1
              IF(SINU(I,J).EQ.-1.E28) THEN
                OUT(J) = '     LIN'
              ELSEIF(SINU(I,J).EQ.-2.E28) THEN
                OUT(J) = '    LIND'
              ELSEIF(SINU(I,J).EQ.-3.E28) THEN
                OUT(J) = '    LINU'
              ELSEIF(SINU(I,J).EQ.-4.E28) THEN
                OUT(J) = 'PARABOLA'
              ELSEIF(SINU(I,J).EQ.-5.E28) THEN
                OUT(J) = '   CUBIC'
              ELSE
                WRITE(OUT(J),'(F8.4)') SINU(I,J)
              ENDIF
            ENDIF
          ENDIF
 300    CONTINUE
 
        IF(SNOUT.GT.0) THEN
          WRITE(STDOUT,58) (OUT(J), J=1,NFL)
        ENDIF
 
 1000 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   SPTOUT
     I                   (STDOUT, STDTAB, N, TAB, LABEL, Z, F)
 
C     + + + PURPOSE + + +
C     Output a standard function table of type 2.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, STDOUT, STDTAB, TAB
      INTEGER LABEL(20)
      REAL F(N), Z(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     N      - Number of tabulated levels in the table
C     TAB    - Table number
C     LABEL  - Label for identification
C     Z      - Argument values for table
C     F      - Function values to output
 
C     + + + SAVED VALUES + + +
      INTEGER TYPE
      SAVE TYPE
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J, L
      REAL Y
      CHARACTER TABID*16
 
C     External names
      INTEGER LENSTR
      CHARACTER GET_TABID*16

      EXTERNAL GET_TABID, LENSTR

C     + + + DATA INITIALIZATIONS + + +
      DATA TYPE/2/
 
C     + + + OUTPUT FORMATS + + +
 1    FORMAT(' TABID= ',A)
 2    FORMAT('TABID= ',A)
 3    FORMAT(' TYPE=',I5)
 4    FORMAT('TYPE=',I5)
 5    FORMAT(' REFL=',F10.3)
 6    FORMAT('REFL=',F10.3)
 7    FORMAT(' ',5X,'DEPTH',4X,'HLCOEF',2X,14A4)
 8    FORMAT(5X,'DEPTH',4X,'HLCOEF',2X,14A4)
 9    FORMAT(' ',F10.3,F10.4)
 10   FORMAT(F10.3,F10.4)
C***********************************************************************
      TABID = GET_TABID(TAB)
      L = LENSTR(TABID)
      WRITE(STDOUT,1) TABID(1:L)
      WRITE(STDTAB,2) TABID(1:L)
 
      WRITE(STDOUT,3) TYPE
      WRITE(STDTAB,4) TYPE
 
      WRITE(STDOUT,5) Z(1)
      WRITE(STDTAB,6) Z(1)
 
      WRITE(STDOUT,7) (LABEL(J),J=1,14)
      WRITE(STDTAB,8) (LABEL(J),J=1,14)
 
      DO 100 J=1,N
        Y=Z(J)-Z(1)
        WRITE(STDOUT,9) Y,F(J)
        WRITE(STDTAB,10) Y,F(J)
 100    CONTINUE
 
C     MAKE LAST ARGUMENT LARGE TO PREVENT TABLE OVERSHOOT
 
      Y = Y + 25.0
      WRITE(STDTAB,10) Y,F(N)
      Y=-1.0
      WRITE(STDTAB,10) Y,Y
      RETURN
      END
C
C
C
      SUBROUTINE   TWDOUT
     I                   (STDOUT, STDTAB, TABLE, LABEL, NHU, NFRAC,
     I                    HUVEC, FDROP, PFDVEC, Q, HDATUM,
     I                    TYPE, SOURCE, zrhufd,
     i                    zone, hgrid, vdatum, unitsys, basis,
     i                    easting, northing,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Output a two-D flow table of the given type. Includes 
C     types 6, 13, and 14.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NFRAC, NHU, STDOUT, STDTAB, TABLE, TYPE
      REAL FDROP(PMXNHU), HDATUM, HUVEC(PMXNHU), PFDVEC(PMXFRC),
     A     Q(PMXNHU,PMXFRC), zrhufd
      real*8 easting, northing
      CHARACTER LABEL*50, SOURCE*8,
     a  zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     TABLE  - Table number
C     LABEL  - Label for identification
C     NHU    - Number of upstream heads
C     NFRAC  - Number of fractions for defining partial free drop
C     HUVEC  - Vector of upstream heads
C     FDROP  - Free drop values
C     PFDVEC - Partial free drop vector
C     Q      - Flowrate
c     zrhufd - gives the drop to free flow at zero upstram head. 
c              Currently is zero for CULVERT, ORIFICE, EMBANKQ, 
c              and may be positive for CHANRAT if the flow surface
c              has a downward slope (sustaining slope).  It is zero
c              in this case also if the slope is zero or retarding. 
c              Only used for fitting cubic splines to the table. 
c              We may add non-zero free drops for other cases if it
c              proves useful. The value is written such that if 
c              it is left out it will be taken as 0.0, the current 
c              value for the drop to free flow when the ups head is zero.
C     HDATUM - Datum for measuring head
C     TYPE   - Table type.
C     SOURCE - source for the table.
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + SAVED VALUES + + +
      INTEGER IOFF
      SAVE IOFF
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, K, L, KNT, LIM
      real*8 hdat
      CHARACTER DUMMY*7, LINE(10)*7, TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LENSTR
      CHARACTER PUT7*7, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL PUT7, GET_TABID, LENSTR
 
C     + + + DATA INITIALIZATIONS + + +
      DATA IOFF/10/
 
C     + + + OUTPUT FORMATS + + +
 21   FORMAT('TYPE=',I5,' HDATUM=',F10.3,1X,A8)
 22   FORMAT('TABID= ',A)
 23   FORMAT('TYPE=',I5,' HDATUM=',F10.3,1X,A8,1x,'zrhufd=',f10.4)
 24   FORMAT('    HEADUP  FREEDROP',1X,A50)
 25   FORMAT(F10.3,F10.5)
 26   FORMAT('    -1.0')
 27   FORMAT(' POINTS ON LINES OF CONSTANT UPSTREAM HEAD')
 28   FORMAT(F10.5,F10.3)
 30   FORMAT('LABEL=',A50)
 32   FORMAT('NHUP=',I5)
 33   FORMAT('NHDN=',I5)
 34   FORMAT('NPFD=',I5)
 35   FORMAT('NPFQ=',I5)
 50   FORMAT('HUP    ',10A7)
 51   FORMAT('QFREE  ',10A7)
 52   FORMAT('FDROP  ',10A7)
 53   FORMAT('HDN    ',10A7)
 54   FORMAT('   PFD',10X,  'Flows for HUP and Proportion of FDROP')
 56   FORMAT(A7,10A7)
 57   FORMAT('   PFQ',10X,'Ups heads for HDN and Proportion of QFREE')
 58   FORMAT(' *ERR:592* TYPE=',I5, ' INVALID FOR 2-D FLOW TABLES.')
 60   format('ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a ' BASIS=',a8,/,'EASTING=',0PF15.3,' NORTHING=',F15.3)
C***********************************************************************
      call chk_vdatum_unitsys(stdout, vdatum, unitsys, 
     a          ' during output of type 6, 13, or 14 table')

c     Apply the global vertical datum shift.  Default must be 0.0
      hdat = hdatum
c      hdat = hdat + g_vdatum_shift

      TABID = GET_TABID(TABLE)
      L = LENSTR(TABID)
       WRITE(STDTAB,22)  TABID(1:L)

      IF(TYPE.EQ.6) THEN
C       OUTPUT THE OLD FORMAT FOR TWO-D FLOW TABLES
 
        WRITE(STDTAB,23)  -TYPE, HDAT, SOURCE, zrhufd
        WRITE(STDTAB,24) LABEL
 
C       PUT IN THE ZERO HEAD POINT
        WRITE(STDTAB,25) 0.0, 0.0
        DO 2000 I=1,NHU
          WRITE(STDTAB,25) HUVEC(I), FDROP(I)
 2000   CONTINUE
        WRITE(STDTAB,26)
        WRITE(STDTAB,27)
        WRITE(STDTAB,28) 0.0, 0.0
        WRITE(STDTAB,28) 1.0, 0.0
        DO 2200 I=1,NHU
          DO 2100 J=1,NFRAC
            WRITE(STDTAB,28) PFDVEC(J), Q(I,J)
 2100     CONTINUE
 2200   CONTINUE
      ELSEIF(TYPE.EQ.13.OR.TYPE.EQ.14) THEN
C       OUTPUT TWO-D FLOW TABLE IN THE RECTANGULAR FORMAT
 
 
        if(type.eq.13) then
          WRITE(STDTAB,23)  -TYPE, HDAT, SOURCE, zrhufd
        else
          WRITE(STDTAB,21)  -TYPE, HDAT, SOURCE
        endif
        if(zone /= 'NONE')  then
c         Output the location information.
          write(stdtab, 60) zone, hgrid, vdatum, unitsys, basis, 
     a             easting,  northing
        endif
        WRITE(STDTAB,30) LABEL
        IF(TYPE.EQ.13) THEN
          WRITE(STDTAB,32) NHU
        ELSE
C         NHU IS REALLY THE NUMBER OF DOWNSTREAM HEADS!
          WRITE(STDTAB,33) NHU
        ENDIF
        IF(TYPE.EQ.13) THEN
          WRITE(STDTAB,34) NFRAC - 1
        ELSE
          WRITE(STDTAB,35) NFRAC - 1
        ENDIF
 
C       OUTPUT THE TABLE WITH IOFF COLUMNS OF DATA UNTIL TABLE IS
C       COMPLETE. THE ZERO HEAD, AND ZERO PROPORTION OF FREE DROP
C       VALUES ARE NOT PRINTED.
 
        KNT = 0
        DO 3000 I=1,NHU,IOFF
          LIM = MIN(IOFF, NHU - KNT)
          DO 100 J=I,I+LIM-1
            LINE(J-KNT) = PUT7(HUVEC(J))
 100      CONTINUE
          IF(TYPE.EQ.13) THEN
            WRITE(STDTAB,50) (LINE(J-KNT), J=I,I+LIM-1)
          ELSE
            WRITE(STDTAB,51) (LINE(J-KNT), J=I,I+LIM-1)
          ENDIF
          DO 101 J=I,I+LIM-1
            LINE(J-KNT) = PUT7(FDROP(J))
 101      CONTINUE
          IF(TYPE.EQ.13) THEN
            WRITE(STDTAB,52) (LINE(J-KNT), J=I,I+LIM-1)
          ELSE
            WRITE(STDTAB,53) (LINE(J-KNT), J=I,I+LIM-1)
          ENDIF
          IF(TYPE.EQ.13) THEN
            WRITE(STDTAB,54)
          ELSE
            WRITE(STDTAB,57)
          ENDIF
 
          DO 2900 J=2,NFRAC
 
            DO 103 K=I,I+LIM-1
              LINE(K-KNT) = PUT7(Q(K,J))
 103        CONTINUE
            DUMMY = PUT7(PFDVEC(J))
            WRITE(STDTAB,56) DUMMY,
     A                       (LINE(K-KNT), K=I,I+LIM-1)
 2900     CONTINUE
          KNT = KNT + IOFF
 3000   CONTINUE
 
      ELSE
        WRITE(STDOUT,58) TYPE
        EFLAG = 1
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   TABOUT
     I                   (STDOUT, STDTAB, TAB, STAT, ZMIN, MESG, SAVOPT,
     I                    OUTOPT, BETOPT, zone, hgrid, vdatum, unitsys, 
     I                    basis, khflag, alphaflag, betaflag, maflag, 
     i                    mqflag,
     M                    NDEP, XST, FTP)
 
C     + + + PURPOSE + + +
C     Output a standard table of elements for a cross section.
      IMPLICIT NONE
C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER FTP, MESG, NDEP, STDOUT, STDTAB, TAB
      REAL STAT, XST(PMXPNT,PMXELM), ZMIN
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c   zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     TAB    - Table number
C     STAT   - Station value
C     ZMIN   - Minimum elevation
C     MESG   - Flag for decreasing conveyance message appearance
C     SAVOPT - Function table saving option
C     OUTOPT - Output option for the table file for cross section function
C               tables
C     BETOPT - Option for computing flux coefficients and critical flow
C     NDEP   - Number of depth values
C     XST    - Storage table for various elements of cross section
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ftable.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'xtadd.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER NSAVE
      PARAMETER(NSAVE=14)
 
C     + + + SAVED VALUES + + +
      INTEGER SOVAL(NSAVE), TYPVEC(NSAVE)
      CHARACTER SOTAB(NSAVE)*8, OUTOTAB(NSAVE)*8
      SAVE SOTAB, SOVAL, TYPVEC
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, I, IT, J, L, N, NSTORE, TYPE, EFLAG,
     a   iy
      REAL DROP, MONE, OLDK, MXSLOT,
     a    yarg,   topf, areaf, jf, kf, betaf,  dkf, dbetaf,
     b       alpf, dalpf, maf, dmaf, mqf, dmqf, dtf, qcf,
     c       kl, km, kr, dkl, dkm, dkr, bl, bm, br, 
     d       dbl, dbm, dbr, alpl, alpm, alpr, dalpl, dalpm, dalpr,
     e       mal, mam, mar, dmal, dmam, dmar, mql, mqm, mqr,
     f       dmql, dmqm, dmqr, max_arg, zmin_out


      CHARACTER ACHR*8, JCHR*8, KCHR*8, QCCHR*8, TCHR*8,
     A          TOPOUT*10, DEPTHOUT*10, AOUT*13,
     B          KHOUT*13, JOUT*13, QCOUT*13,
     c          dkhout*10, dbetaout*10, dalphaout*10, 
     d          dmaout*10, dmqout*10
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LENSTR
      CHARACTER PUT8*8, GET_TABID*16, PUT10*10
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BINSER, CHKTAJ, CUTTAB, PUT8, PUT_STRING_IN_FT,
     A         PUT_DP_IN_FT, GET_TABID, LENSTR, var_decimal,
     B         PUT10
 
C     + + + DATA INITIALIZATIONS + + +
      DATA
     1SOTAB /'SAVE1   ','SAVE12  ','SAVE20  ','SAVE21  ','SAVE22  ',
     A       'SAVE23  ','SAVE24  ','SAVE25  ','SAVE30  ','SAVE31  ',
     B       'SAVE32  ','SAVE33  ','SAVE34  ','SAVE35  '/
      DATA
     1OUTOTAB 
     2   /  'OUT1   ','OUT12  ','OUT20  ','OUT21  ','OUT22  ',
     A       'OUT23  ','OUT24  ','OUT25  ','OUT30  ','OUT31  ',
     B       'OUT32  ','OUT33  ','OUT34  ','OUT35  '/
      DATA
     1SOVAL /      6,         8,         5,         6,         8,
     A             7,         8,        10,         7,         8,
     B            11,        11,        12,        15  /
      DATA
     1TYPVEC/     21,        22,        20,        21,         22,
     A            23,        24,        25,        30,         31,
     B            32,        33,        34,        35 /
 

C     + + + OUTPUT FORMATS + + +
 40   FORMAT(A10,A10,2a13,F10.6,2a10)
 41   FORMAT(A10,A10,2a13,F10.6,A13,2a10)
 42   FORMAT(A10,A10,2a13,F10.6,A13,F10.6,A13,3a10)

 43   FORMAT(A10,A10,2a13,F10.6,2f10.7,4a10)
 44   FORMAT(A10,A10,2a13,F10.6,A13,2F10.7,4a10)
 45   FORMAT(A10,A10,2a13,F10.6,A13,F10.6,A13,2F10.7,5a10)
 50   FORMAT(/,' *BUG:XXX* Save option=',A8,' unknown in TABOUT.')
 51   FORMAT(/,' Elements for section ID: ',A)
 52   FORMAT(/,' TABID= ',A,/,1X,'TYPE=',I5,/,1X,'STATION=',1PE12.5,
     1        /,1X,'ELEVATION=',1PE12.5)
c     English unit system without zone, etc and with zone
 530   FORMAT('TABID= ',A,/,'TYPE=',I5,/,
     2         'STATION=',1PE12.5,
     3         ' GISID=',A,' EASTING=',0PF15.2,' NORTHING=',F15.2,
     4        /,'ELEVATION=',1PE12.5,
     5    17X,' EXT=-99.900000',' FAC=1.000 SLOT=',0PF10.3)
531   FORMAT('TABID= ',A,/,'TYPE=',I5,/,
     1 'ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     2 ' BASIS=',a8,/, 'STATION=', 1PE12.5,
     3        ' GISID=',A,' EASTING=',0PF15.3,' NORTHING=',F15.3,
     4        /,'ELEVATION=',1PE12.5,
     5    17X,' EXT=-99.900000',' FAC=1.000 SLOT=',0PF10.4)

c     Metric unit system without zone, etc and with zone
 532   FORMAT('TABID= ',A,/,'TYPE=',I5,/,
     2         'STATION=',1PE12.5,
     3         ' GISID=',A,' EASTING=',0PF15.2,' NORTHING=',F15.2,
     4        /,'ELEVATION=',1PE12.5,
     5    17X,' EXT=-30.449520',' FAC=1.000 SLOT=',0PF10.3)
533   FORMAT('TABID= ',A,/,'TYPE=',I5,/,
     1 'ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     2 ' BASIS=',a8,/, 'STATION=', 1PE12.5,
     3        ' GISID=',A,' EASTING=',0PF15.3,' NORTHING=',F15.3,
     4        /,'ELEVATION=',1PE12.5,
     5    17X,' EXT=-30.449520',' FAC=1.000 SLOT=',0PF10.4)
 54   FORMAT(/,3X,'Depth Top Wid',4X,'Area',' Sqr(Cv)',3X,
     1       'Beta',7X,'J',2X,'Alpha',3X,'CritQ',5X,'Ma',
     2       5X,'Mq','  AVER n',
     3 '   dkh/dy  dbta/dy  dalp/dy   dma/dy   dmq/dy')
 55   FORMAT(/,' *BUG:XXX* Out option=',A8,' unknown in TABOUT.')
 56   FORMAT(' ',F7.3,3A8,F7.3,A8,F7.3,A8,2F7.3,F8.4,
     a         f8.3,a1,4(f8.4,a1))
 57   FORMAT(F10.0)
 58   FORMAT(/,' *WRN:506* Conveyance non-increasing at depth=',
     A  F10.2,/,11X,' decrease=',F6.1,' per cent')
 59   FORMAT(/,' *BUG:XXX* Invalid index=',I5,' for type branch',
     A       ' in SAVE option in TABOUT.')
 60   FORMAT(' *ERR:569* Insufficient space in ITAB/FTAB to save',
     A     ' cross section table.',/,10X,' need',I5,' more',
     B     ' elements.')
 62   FORMAT(' ','     Depth      Beta   dBetadY     Alpha  DAlphadY',
     A'     Qbeta    Qalpha     Qcrit')
 64   FORMAT(1X,F10.5,4F10.6,3F10.1)


 70   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',6x,
     a       'Beta    dkh/dy  dbeta/dy')
 71   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment    dkh/dy  dbeta/dy')
 72   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',5x,'Alpha',8x,
     b          'Critq    dkh/dy  dbeta/dy dalpha/dy')
 73   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta',8x,'Ma',8x,
     a       'Mq    dkh/dy  dbeta/dy    dma/dy    dmq/dy')
 74   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',8x,'Ma',8x,
     A    'Mq    dkh/dy  dbeta/dy    dma/dy    dmq/dy')
 75   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',5x,'Alpha',8x,'Critq',
     2       8x,'Ma',8x,
     a    'Mq    dkh/dy  dbeta/dy    dma/dy    dmq/dy dalpha/dy')

 80   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',6x,'Beta')
 81   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment')
 82   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',5x,'Alpha',8x,'Critq')
 83   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta',8x,'Ma',8x,'Mq')
 84   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',8x,'Ma',8x,'Mq')
 85   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',5x,'Alpha',8x,'Critq',
     2       8x,'Ma',8x,'Mq')
 90   FORMAT(A10,A10,2a13,F10.6)
 91   FORMAT(A10,A10,2a13,F10.6,A13)
 92   FORMAT(A10,A10,2a13,F10.6,A13,F10.6,A13)

 93   FORMAT(A10,A10,2a13,F10.6,2f10.7)
 94   FORMAT(A10,A10,2a13,F10.6,A13,2F10.7)
 95   FORMAT(A10,A10,2a13,F10.6,A13,F10.6,A13,2F10.7)

 98   format(/,' Note: A caret (^) following a derivative indicates',
     a' imposition of variation limitation at that point.  Continuity',
     b' of',/,7x,' second derivative is lost at these points.',
     c'  Otherwise both the first and second derivative is continuous.')
99    FORMAT(/,' Saving TabId= ',A,'  Internal table number=',I5)
C***********************************************************************
      EFLAG = 0

      IF(GRAV.GT.15.0) THEN
        MXSLOT = 0.07
      ELSE
        MXSLOT = 0.02134
      ENDIF
      MONE=-1.
 
      TABID = GET_TABID(TAB)
      L = LENSTR(TABID)
      WRITE(STDOUT,51) TABID(1:L)
 
C     TRIM TOP OF TABLE IF IT IS A CLOSED CONDUIT AND THEN CHECK FOR
C     CONSISTENCY OF THE TOP WIDTH, AREA, AND FIRST MOMENT OF AREA
 
      CALL CUTTAB
     M           (NDEP, XST)
      CALL CHKTAJ
     I           (STDOUT, NDEP, XST,
     O            FLAG)

 
      if(outopt.eq.' ') then
        type = -25
      elseif(outopt.eq.'NOOUT') then
        type = -25
      else
        CALL BINSER
     I             (OUTOPT, NSAVE, OUTOTAB,
     O              IT)
      
        IF(IT.EQ.0) THEN
          WRITE(STDOUT,55) OUTOPT
          STOP 'Abnormal stop. Errors found.'
        else
          type = -typvec(it)
        endif
      endif
 
      IF(GISID.EQ.' ') THEN
        GISID = TABID
      ENDIF
      WRITE(STDOUT,52) TABID(1:L), TYPE, STAT, zmin
      IF(OUTOPT.NE.'NOOUT'.AND.FLAG.EQ.0) THEN
        IF(GRAV.GT.15.0) THEN
          if(zone =='NONE') then
            WRITE(STDTAB,530) TABID(1:L), TYPE,
     A                      STAT, GISID, EASTING, NORTHING, zmin,
     B                      SLOT_DEPTH
          else
            WRITE(STDTAB,531) TABID(1:L), TYPE,
     1                      zone, hgrid, vdatum, unitsys, basis,  
     A                      STAT, GISID, EASTING, NORTHING, zmin,
     B                      SLOT_DEPTH
          endif
        ELSE
          if(zone == 'NONE') then
            WRITE(STDTAB,532) TABID(1:L), TYPE,                          
     A                       STAT, GISID, EASTING, NORTHING, zmin,
     B                       SLOT_DEPTH
          else
            WRITE(STDTAB,533) TABID(1:L), TYPE,                          
     1                       zone, hgrid, vdatum, unitsys, basis,  
     A                       STAT, GISID, EASTING, NORTHING, zmin,
     B                       SLOT_DEPTH
          endif
        ENDIF
      ENDIF
 
      WRITE(STDOUT,54)
      IF(FLAG.EQ.0) THEN
        IF(OUTOPT.EQ.'OUT20') THEN
          WRITE(STDTAB,80)
        ELSEIF(OUTOPT.EQ.'OUT21'.OR.OUTOPT.EQ.'OUT1') THEN
          WRITE(STDTAB,81)
        ELSEIF(OUTOPT.EQ.'OUT22'.OR.OUTOPT.EQ.'OUT12') THEN
          WRITE(STDTAB,82)
        ELSEIF(OUTOPT.EQ.'OUT23') THEN
          WRITE(STDTAB,83)
        ELSEIF(OUTOPT.EQ.'OUT24') THEN
          WRITE(STDTAB,84)
        ELSEIF(OUTOPT.EQ.'OUT25') THEN
          WRITE(STDTAB,85)
        elseIF(OUTOPT.EQ.'OUT30') THEN
          WRITE(STDTAB,70)
        ELSEIF(OUTOPT.EQ.'OUT31') THEN
          WRITE(STDTAB,71)
        ELSEIF(OUTOPT.EQ.'OUT32') THEN
          WRITE(STDTAB,72)
        ELSEIF(OUTOPT.EQ.'OUT33') THEN
          WRITE(STDTAB,73)
        ELSEIF(OUTOPT.EQ.'OUT34') THEN
          WRITE(STDTAB,74)
        ELSEIF(OUTOPT.EQ.'OUT35') THEN
          WRITE(STDTAB,75)
        ENDIF
      ENDIF

      call chk_vdatum_unitsys(stdout, vdatum, unitsys, 
     a          ' during output/store of a cross-section table')


 
      OLDK=-1.0
      DO 100 J=1,NDEP
 
C       OUTPUT ALL VALUES TO STDOUT EVEN THOUGH NOT REQUESTED FOR THE
C       TABLE TO BE OUTPUT TO STDTAB OR TO BE STORED.
 
        TCHR = PUT8(XST(J,2))
        ACHR = PUT8(XST(J,3))
        KCHR = PUT8(XST(J,5))
        JCHR = PUT8(XST(J,4))
        QCCHR = PUT8(XST(J,13))

        WRITE(STDOUT,56) XST(J,1), TCHR, ACHR, KCHR,XST(J,6),
     A               JCHR, XST(J,7), QCCHR, XST(J,14), XST(J,15),
     B               XST(J,17), 
     c               xst(j,18), khflag(j),
     c               xst(j,20), betaflag(j),
     c               xst(j,19), alphaflag(j),
     c               xst(j,21), maflag(j),
     c               xst(j,22), mqflag(j)
     c                    
C     B               POWERK_VS_J, POWERQC_VS_J, POWERK_VS_A,
C     C               POWERQC_VS_A
 
        IF(FLAG.EQ.0) THEN
          call var_decimal(xst(j,1),
     o                        depthout)
          call var_decimal(xst(j,2),
     o                        topout)
          call var_decimal(xst(j,3),
     o                        aout)
          call var_decimal(xst(j,4),
     o                        jout)
          call var_decimal(xst(j,5),
     o                        khout)
          call var_decimal(xst(j,13),
     o                        qcout)
C
C          call var_decimal(xst(j,18),
C     o                        dkhout)
C          call var_decimal(xst(j,19),
C     o                        dalphaout)
C          call var_decimal(xst(j,20),
C     o                        dbetaout)
C          call var_decimal(xst(j,21),
C     o                        dmaout)
C          call var_decimal(xst(j,22),
C     o                        dmqout)
          
           dkhout     =  put10(xst(j,18))
           dalphaout  =  put10(xst(j,19))
           dbetaout   =  put10(xst(j,20))
           dmaout     =  put10(xst(j,21))
           dmqout     =  put10(xst(j,22))
                      

          IF(OUTOPT.EQ.'OUT20') THEN
            WRITE(STDTAB,90) depthout,topout,aout,khout,
     1                     XST(J,6)
          ELSEIF(OUTOPT.EQ.'OUT21'.OR.OUTOPT.EQ.'OUT1') THEN
            WRITE(STDTAB,91) depthout,topout,aout,khout,
     1                     XST(J,6),jout
          ELSEIF(OUTOPT.EQ.'OUT22'.OR.OUTOPT.EQ.'OUT12') THEN
            WRITE(STDTAB,92) depthout,topout,aout,khout,
     1                     XST(J,6),jout, XST(J,7), qcout
          ELSEIF(OUTOPT.EQ.'OUT23') THEN
            WRITE(STDTAB,93) depthout,topout,aout,khout,
     1                     XST(J,6), XST(J,14), XST(J,15)
          ELSEIF(OUTOPT.EQ.'OUT24') THEN
            WRITE(STDTAB,94) depthout,topout,aout,khout,
     1                     XST(J,6),jout, XST(J,14), XST(J,15)
          ELSEIF(OUTOPT.EQ.'OUT25') THEN
            WRITE(STDTAB,95) depthout,topout,aout,khout,
     1                     XST(J,6),jout, XST(J,7), qcout,
     2                     XST(J,14), XST(J,15)
          elseIF(OUTOPT.EQ.'OUT30') THEN
            WRITE(STDTAB,40) depthout,topout,aout,khout,
     1                     XST(J,6), dkhout, dbetaout
          ELSEIF(OUTOPT.EQ.'OUT31') THEN
            WRITE(STDTAB,41) depthout,topout,aout,khout,
     1                     XST(J,6),jout, dkhout, dbetaout
          ELSEIF(OUTOPT.EQ.'OUT32') THEN
            WRITE(STDTAB,42) depthout,topout,aout,khout,
     1                     XST(J,6),jout, XST(J,7), qcout,
     a                     dkhout, dbetaout, dalphaout
          ELSEIF(OUTOPT.EQ.'OUT33') THEN
            WRITE(STDTAB,43) depthout,topout,aout,khout,
     1                     XST(J,6), XST(J,14), XST(J,15),
     a                    dkhout, dbetaout, dmaout, dmqout
          ELSEIF(OUTOPT.EQ.'OUT34') THEN
            WRITE(STDTAB,44) depthout,topout,aout,khout,
     1                     XST(J,6),jout, XST(J,14), XST(J,15),
     a                    dkhout, dbetaout, dmaout, dmqout
          ELSEIF(OUTOPT.EQ.'OUT35') THEN
            WRITE(STDTAB,45) depthout,topout,aout,khout,
     1                     XST(J,6),jout, XST(J,7), qcout,
     2                     XST(J,14), XST(J,15),
     a                    dkhout, dbetaout, dmaout, dmqout, dalphaout
          ENDIF
        ENDIF
 
C       CHECK FOR NON-INCREASING CONVEYANCE
 
        IF(XST(J,5).LE.OLDK.AND.MESG.EQ.1) THEN
          IF(XST(J,2).GT.MXSLOT) THEN
            DROP = ((XST(J,5)**2 - OLDK**2)/OLDK**2)*100.0
            WRITE(STDOUT,58) XST(J,1), ABS(DROP)
          ENDIF
        ENDIF
        OLDK=XST(J,5)
 
 100  CONTINUE
      IF(FLAG.EQ.0.AND.OUTOPT.NE.'NOOUT') THEN
C       WRITE THE TERMINATING LINE FOR EACH TABLE TO STDTAB.
        WRITE(STDTAB,57) MONE
      ENDIF
 
      write(stdout, 98)
 
      IF(SAVOPT.NE.'NOSAVE'.AND.FLAG.EQ.0) THEN
C       WE HAVE A SAVE OPERATION TO COMPLETE.  FIND THE NUMBER OF
C       VALUES TO STORE FOR EACH TABLE TYPE
 
        CALL BINSER
     I             (SAVOPT, NSAVE, SOTAB,
     O              IT)
        IF(IT.EQ.0) THEN
          WRITE(STDOUT,50) SAVOPT
          STOP 'Abnormal stop. Errors found.'
        ELSE
          NSTORE = SOVAL(IT)
          TYPE = TYPVEC(IT)
   
C         IS THERE ENOUGH ROOM IN FTAB/ITAB FOR THE TABLE?
 
          N = XTIOFF + NDEP*NSTORE
          IF(FTP+N.GT.MRFTAB) THEN
            WRITE(STDOUT,60) FTP + N - MRFTAB
            STOP 'Abnormal stop. Errors found.'
          ENDIF
 
          WRITE(STDOUT,99) TABID(1:L), TAB 
C         STORE THE HEADER INFORMATION FOR THE TABLE
          IF(FTPNT(TAB).NE.0) THEN
            WRITE(STDOUT,*) ' FTPNT(TAB)=',FTPNT(TAB),' TAB=',TAB
            CALL KIL
     I               (25,
     M               TAB, EFLAG)
            STOP 'Abnormal stop. Error found.' 
          ENDIF
 
c         Store the table for subsequent use in this run of FEQUTL.
          max_arg = xst(ndep,1)
          call store_xsec_tab(stdout, ndep, type, tabid, tab,
     i               stat, zmin, slot_depth, gisid, easting, northing,
     i               xst, max_arg, zone, hgrid, vdatum, unitsys, basis,
     m               ftp) 


c         Add test of the lookup in types 30-35.
          if(type.eq.9999) then
            write(stdout,*) ' '
            write(stdout,*) ' Testing lookup for types: 30-35'


            call XLKT25
     I                 (FTPNT(TAB),
     M                  xst(1,1),
     O    areaf, topf, dtf,jf, kl, dkl, bl, dbl, alpl, dalpl, 
     o    qcf, mal, dmal, mql, dmql)

            do iy=2,ndep 
              yarg = 0.5*(xst(iy-1,1) + xst(iy,1))

              call XLKT25
     I                   (FTPNT(TAB),
     M                    yarg,
     O    areaf, topf, dtf, jf, km, dkm, bm, dbm, alpm, dalpm, 
     o    qcf, mam, dmam, mqm, dmqm)

              call XLKT25
     I                   (FTPNT(TAB),
     M                    xst(iy,1),
     O    areaf, topf, dtf, jf, kr, dkr, br, dbr, alpr, dalpr, 
     o    qcf, mar, dmar, mqr, dmqr)

              write(stdout,9124) xst(iy,1),
     a((xst(iy,1) - xst(iy-1,1))*(dkl + 4.*dkm + dkr)/6. - (kr - kl))/
     a     max(abs(kr-kl), 0.01),
     b((xst(iy,1) - xst(iy-1,1))*(dbl + 4.*dbm + dbr)/6. - (br - bl))/
     b    max(abs(br-bl),0.002),
     b((xst(iy,1) - xst(iy-1,1))*(dalpl + 4.*dalpm + dalpr)/6. - 
     b (alpr - alpl))/max(abs(alpr-alpl),0.002),
     b((xst(iy,1) - xst(iy-1,1))*(dmal + 4.*dmam + dmar)/6. -
     c   (mar - mal))/max(abs(mar-mal),0.002),
     b((xst(iy,1) - xst(iy-1,1))*(dmql + 4.*dmqm + dmqr)/6. -
     c   (mqr - mql))/max(abs(mqr-mql),0.002)

9124  format(f10.4, 6(1pe12.4))
              

                kl = kr
                dkl = dkr
                bl = br
                dbl = dbr
                alpl = alpr
                dalpl = dalpr
                mal = mar
                dmal = dmar
                mql = mqr
                dmql = dmqr
               write(stdout,9123) xst(iy,1), topf - xst(iy,2),
     a           areaf - xst(iy,3),
     a           jf - xst(iy,4),
     a           kr - xst(iy,5)**2,
     b           dkr - 2.*xst(iy,5)*xst(iy,18),
     b           qcf - xst(iy,13),
     c           br - xst(iy,6), dbr - xst(iy,20),
     c           alpr - xst(iy,7), dalpr - xst(iy,19),
     d           mar - xst(iy,14), dmar - xst(iy,21),
     d           mqr - xst(iy,15), dmqr - xst(iy,22)

9123  format(f10.4,15(1pe12.4))

            end do
          endif



        ENDIF
      ENDIF
 
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


      IF(BETOPT(1:7).EQ.'NEWBETA'.AND.FLAG.EQ.0) THEN
        WRITE(STDOUT,*) ' '
        WRITE(STDOUT,*) ' Estimates of critical flow:'
        WRITE(STDOUT,62)
        DO 900 J=1,NDEP
C                          Y,        Beta,    dBeta/dy   Alpha   
C      dAlpha/dy  Qcm           Qce        Qcsimple
          WRITE(STDOUT,64) XST(J,1), XST(J,6), XST(J,8), XST(J,7), 
     A XST(J,9),  XST(J,10),    XST(J,11), XST(J,12)
 900    CONTINUE
      ENDIF

      IF(FLAG.EQ.1) THEN
        STOP 'Abnormal stop. Potential BUG found.'
      ENDIF
      RETURN
      END

C
C
C
      SUBROUTINE GET_XS_LOCATION(STDOUT, RM,
     O                           EASTING, NORTHING)
C     Get the location of the cross section from the stored flow-line data

      IMPLICIT NONE
      INTEGER STDOUT
      REAL*8 RM, EASTING, NORTHING

      INCLUDE 'flowline.cmn'

C     Local

C     **********************FORMATS*************************************
50    FORMAT(/,' *ERR:757* Flow-line station=',F12.4,
     A   ' overflows flow-line data')
52    FORMAT(/,' *ERR:758* Flow-line station=',F12.4,
     A     ' underflows flow-line data')
C***********************************************************************
C     Find interval in flow-line distance that contains RM.

      IF(RM.GE.FL_DISTANCE(L)) THEN
C       CHECK FOR ARGUMENT ABOVE MAX ARG IN THE TABLE
        IF(RM.GT.FL_DISTANCE(FL_PNT_KNT)) THEN
          WRITE(STDOUT,50) RM
          STOP 'Abnormal stop. Error found.'
        ELSE
 100      CONTINUE
            IF(RM.GT.FL_DISTANCE(L+1)) THEN
              L = L + 1
              GOTO 100
            ENDIF
        ENDIF
      ELSE
C       CHECK FOR ARGUMENT BELOW MIN ARG IN THE TABLE
        IF(RM.LT.FL_DISTANCE(1)) THEN
          WRITE(STDOUT,52) RM
          STOP 'Abnormal stop.  Error found.'
        ELSE
 110      CONTINUE
            L = L - 1
            IF(RM.LT.FL_DISTANCE(L)) GOTO 110
        ENDIF
      ENDIF
C     At this point L defines the low argument end of the
C     interval containing the argument.

C     Interpolate linearly for the values

      EASTING = FL_EASTING(L) + (RM - FL_DISTANCE(L))*
     A (FL_EASTING(L+1) - FL_EASTING(L))/
     B(FL_DISTANCE(L+1) - FL_DISTANCE(L))
      NORTHING = FL_NORTHING(L) + (RM - FL_DISTANCE(L))*
     A(FL_NORTHING(L+1) - FL_NORTHING(L))/
     B(FL_DISTANCE(L+1) - FL_DISTANCE(L))

      RETURN
      END



C
C
C
      SUBROUTINE   IXTOUT
     I                   (STDTAB, STDOUT, ADRS, BEGTAB, ENDTAB)
 
C     + + + PURPOSE + + +
C     Output the cross section table at ADRS. Output a label
C     to show that it is an interpolated table.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, BEGTAB, ENDTAB, STDOUT, STDTAB
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDTAB - Fortran unit number for output of function tables
C     STDOUT - Fortran unit number for standard user output
C     ADRS   - Address of function table
C     BEGTAB - Upstream table number for cross section used to define
C               the interpolated cross section table being output
C     ENDTAB - Downstream table number for cross section used to define
C               the interpolated cross section table being output
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'flowline.cmn'
      include 'grid_datum.cmn'

C     + + + SAVED VALUES + + +
      CHARACTER PRFORM(12)*80
      SAVE PRFORM
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HA, I, J, LA, TABN, TYPE, XOFF,
     A        LTAB, LBEG, LEND
      REAL ELEV, STAT

      REAL*8 EASTING, NORTHING, RM

      CHARACTER*16 TABID, BEGTABID, ENDTABID, GISID
      character*8 zone, hgrid, vdatum, unitsys, basis

C     Called program units
      INTEGER LENSTR
      CHARACTER GET_TABID*16, get_chr8_from_ft*8
      EXTERNAL LENSTR, GET_TABID, GET_XS_LOCATION,
     a    get_chr8_from_ft

C     + + + DATA INITIALIZATIONS + + +
      DATA PRFORM/
     1'(F10.5,F10.3,1P,2E13.5,0PF10.4)',
     2'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5)',
     3'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,0PF10.4,1PE13.5)',
     4'(F10.5,F10.3,1P,2E13.5,0P3F10.4)',
     5'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,0P2F10.4)',
     6'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,0PF10.4,1PE13.5,0P2F10.4)'
     b,
     7'(F10.5,F10.3,1P,2E13.5,0PF10.4,2(1pe10.4e1))',
     8'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,2(1pe10.4e1))',
     9'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,0PF10.4,1PE13.5,3(1pe10.4e
     *1))',            
     a'(F10.5,F10.3,1P,2E13.5,0P3F10.4,4(1pe10.4e1))',
     b'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,0P2F10.4,4(1pe10.4e1))',
     c'(F10.5,F10.3,1P,2E13.5,0PF10.4,1PE13.5,0PF10.4,1PE13.5,0P2F10.4,5
     *(1pe10.4e1))'/
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT('TABID= ',A,/,'TYPE=',I5,/,'STATION=',1PE12.5,
     A         ' GISID=',A,' EASTING=',0PF15.2,' NORTHING=',F15.2,
     B       /,';', 5X,'Interpolated between tables: ',A,' and ',A,
     C        /,'ELEVATION=',1PE12.5,
     D    17X,' EXT=-99.900000',' FAC=1.000 SLOT=',0PF10.3)
521   FORMAT('TABID= ',A,/,'TYPE=',I5,/,'STATION=',1PE12.5,
     A         ' GISID=',A,' EASTING=',0PF15.2,' NORTHING=',F15.2,
     B       /,';', 5X,'Interpolated between tables: ',A,' and ',A,
     C        /,'ELEVATION=',1PE12.5,
     D    17X,' EXT=-30.449520',' FAC=1.000 SLOT=',0PF10.3)

 53   FORMAT('TABID= ',A,/,'TYPE=',I5,/,
     1 'ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     2    ' BASIS=',a8,/, 'STATION=',1PE12.5,
     A         ' GISID=',A,' EASTING=',0PF15.2,' NORTHING=',F15.2,
     B       /,';', 5X,'Interpolated between tables: ',A,' and ',A,
     C        /,'ELEVATION=',1PE12.5,
     D    17X,' EXT=-99.900000',' FAC=1.000 SLOT=',0PF10.3)
531   FORMAT('TABID= ',A,/,'TYPE=',I5,/,
     1 'ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     2 ' BASIS=',a8,/, 'STATION=',1PE12.5,
     A         ' GISID=',A,' EASTING=',0PF15.2,' NORTHING=',F15.2,
     B       /,';', 5X,'Interpolated between tables: ',A,' and ',A,
     C        /,'ELEVATION=',1PE12.5,
     D    17X,' EXT=-30.449520',' FAC=1.000 SLOT=',0PF10.3)

 70   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',6x,
     a       'Beta    dkh/dy  dbeta/dy')
 71   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment    dkh/dy  dbeta/dy')
 72   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',5x,'Alpha',8x,
     b          'Critq    dkh/dy  dbeta/dy dalpha/dy')
 73   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta',8x,'Ma',8x,
     a       'Mq    dkh/dy  dbeta/dy    dma/dy    dmq/dy')
 74   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',8x,'Ma',8x,
     A    'Mq    dkh/dy  dbeta/dy    dma/dy    dmq/dy')
 75   format(5x,'Depth Top_width',9x,'Area   Sqrt(Conv)',
     1       6x,'Beta First_moment',5x,'Alpha',8x,'Critq',
     2       8x,'Ma',8x,
     a    'Mq    dkh/dy  dbeta/dy dalpha/dy    dma/dy    dmq/dy')


 80   FORMAT(5X,'DEPTH TOP WIDTH',9X,'AREA   SQRT(CONV)',6X,'BETA')
 81   FORMAT(5X,'DEPTH TOP WIDTH',9X,'AREA   SQRT(CONV)',
     1       6X,'BETA FIRST MOMENT')
 82   FORMAT(5X,'DEPTH TOP WIDTH',9X,'AREA   SQRT(CONV)',
     1       6X,'BETA FIRST MOMENT',5X,'ALPHA',8X,'CRITQ')
 83   FORMAT(5X,'DEPTH TOP WIDTH',9X,'AREA   SQRT(CONV)',
     1       6X,'BETA',8X,'MA',8X,'MQ')
 84   FORMAT(5X,'DEPTH TOP WIDTH',9X,'AREA   SQRT(CONV)',
     1       6X,'BETA FIRST MOMENT',8X,'MA',8X,'MQ')
 85   FORMAT(5X,'DEPTH TOP WIDTH',9X,'AREA   SQRT(CONV)',
     1       6X,'BETA FIRST MOMENT',5X,'ALPHA',8X,'CRITQ',
     2       8X,'MA',8X,'MQ')
C***********************************************************************
C     GET KEY VALUES
 
      TABN = ITAB(ADRS+1)
      TYPE = ITAB(ADRS+2)
      HA = ITAB(ADRS)
      LA = ADRS + XTIOFF
      STAT = FTAB(ADRS+4)
      ELEV = FTAB(ADRS+5)

c     30 October 2007:  NO, NO, NO do not apply the datum shift
c     when interpolating cross sections.  If the datum shift is applied
c     to all cross sections when they are input to fequtl.  Therefore
c     applying it to interpolated sections is incorrect.  The datum shift
c     requested by the user is already present in the invert elevations
c     of the cross sections used as the basis for interpolation!
c     Apply the datum shift to the datum for depth
c      elev = elev + g_vdatum_shift

      TABID = GET_TABID(TABN)
      GISID = TABID
      LTAB = LENSTR(TABID)
      BEGTABID = GET_TABID(BEGTAB)
      LBEG = LENSTR(BEGTABID)
      ENDTABID = GET_TABID(ENDTAB)
      LEND = LENSTR(ENDTABID)

c     Get the location values (other than easting and northing),
c     from the ADRS. 
      zone = get_chr8_from_ft(adrs+22)
      hgrid = get_chr8_from_ft(adrs+24)
      vdatum = get_chr8_from_ft(adrs+26)
      unitsys = get_chr8_from_ft(adrs+28)
      basis = get_chr8_from_ft(adrs+30)


      IF(FL_PRESENT.EQ.1) THEN
C       We have a defining flow line.  Compute the location for the 
C       interpolated table.
        RM = DBLE(STAT)
        CALL GET_XS_LOCATION(STDOUT, RM,
     O                       EASTING, NORTHING)
      ELSE
        EASTING = 0.D0
        NORTHING = 0.D0
      ENDIF
 
      if(zone == 'NONE') then
        IF(GRAV.GT.15.0) THEN
          WRITE(STDTAB,52) TABID(1:LTAB), -TYPE, STAT,
     A           GISID, EASTING, NORTHING, BEGTABID(1:LBEG), 
     B                         ENDTABID(1:LEND), ELEV,
     C                         FTAB(ADRS+7)
        ELSE
          WRITE(STDTAB,521) TABID(1:LTAB), -TYPE, STAT,
     A           GISID, EASTING, NORTHING, BEGTABID(1:LBEG), 
     B                         ENDTABID(1:LEND), ELEV,
     C                         FTAB(ADRS+7)
        ENDIF
      else
        IF(GRAV.GT.15.0) THEN
          WRITE(STDTAB,53) TABID(1:LTAB), -TYPE,
     a      zone, hgrid, vdatum, unitsys, basis, STAT,
     A           GISID, EASTING, NORTHING, BEGTABID(1:LBEG), 
     B                         ENDTABID(1:LEND), ELEV,
     C                         FTAB(ADRS+7)
        ELSE
          WRITE(STDTAB,531) TABID(1:LTAB), -TYPE, 
     a      zone, hgrid, vdatum, unitsys, basis, STAT,
     A           GISID, EASTING, NORTHING, BEGTABID(1:LBEG), 
     B                         ENDTABID(1:LEND), ELEV,
     C                         FTAB(ADRS+7)
        ENDIF
      endif
      
      XOFF = OFFVEC(TYPE)
 
      if(type.ge.30) then
         type = type - 4
      endif
      GOTO(20, 21, 22, 23, 24, 25, 30, 31, 32, 33, 34, 35), TYPE-19
 
 20   CONTINUE
        WRITE(STDTAB,80)
        GOTO 40
 21   CONTINUE
        WRITE(STDTAB,81)
        GOTO 40
 22   CONTINUE
        WRITE(STDTAB,82)
        GOTO 40
 23   CONTINUE
        WRITE(STDTAB,83)
        GOTO 40
 24   CONTINUE
        WRITE(STDTAB,84)
        GOTO 40
 25   CONTINUE
        WRITE(STDTAB,85)
        GOTO 40

 30   CONTINUE
        WRITE(STDTAB,70)
        GOTO 40
 31   CONTINUE
        WRITE(STDTAB,71)
        GOTO 40
 32   CONTINUE
        WRITE(STDTAB,72)
        GOTO 40
 33   CONTINUE
        WRITE(STDTAB,73)
        GOTO 40
 34   CONTINUE
        WRITE(STDTAB,74)
        GOTO 40
 35   CONTINUE
        WRITE(STDTAB,75)
        GOTO 40


 40   CONTINUE
 
 
      DO 100 I=LA,HA,XOFF
          WRITE(STDTAB,PRFORM(TYPE-19)) (FTAB(I+J),J=0,XOFF-1)
 100  CONTINUE
      WRITE(STDTAB,'(F10.0)') -1.0
 
      RETURN
      END
