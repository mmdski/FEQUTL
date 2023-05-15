C
C
C
      REAL FUNCTION   FISE
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute the residual for inversion of specific energy.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'isecom.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER NOALP(6)
      SAVE NOALP
 
C     + + + LOCAL VARIABLES + + +
      REAL ALP, AT, B, DALP, DB, DK, DT, J, K, QC, T, YLOC
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKT20, XLKT22
 
C     + + + DATA INITIALIZATIONS + + +
      DATA NOALP/1, 1, 0, 1, 1, 0/
C***********************************************************************
      YLOC = Y
      IF(NOALP(TYPE-19).EQ.1) THEN
C       NO ALPHA IN THE TABLE
        CALL XLKT20
     I             (ADRST,
     M              YLOC,
     O              AT, T, DT, K, DK, B, DB)
 
         FISE = (Y + 0.5*QT**2/(GRAVT*AT**2))/ET - 1.0
      ELSE
C       TABLE DOES HAVE ALPHA.  USE IT.
        CALL XLKT22
     I             (ADRST,
     M              YLOC,
     O              AT, T, DT, J, K, DK, B, DB, ALP, DALP, QC)
 
        FISE = (Y + 0.5*ALP*QT**2/(GRAVT*AT**2))/ET - 1.0
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   INVTSE
     I                   (GRAV, STDOUT, ADRS, ATYPE, Q, E, Y,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     For a given flow rate, Q, specific energy, E, and
C     cross section table address, ADRS, compute the depth of
C     water in the cross section using the initial guess, Y,
C     as a starting point.  The goal is to find a subcritical
C     solution.   There should always be at least one subcritical
C     solution if there is any solution at all.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, ATYPE, EFLAG, STDOUT
      REAL E, GRAV, Q, Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     GRAV   - value of acceleration due to gravity
C     STDOUT - Fortran unit number for user output and messages
C     ADRS   - Address of function table
C     ATYPE  - Type of the cross section table involved
C     Q      - Flowrate
C     E      - Specific energy value.
C     Y      - maximum depth in a cross section
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'isecom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL YD, YU
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FISE, REGFAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:503* GUESSED DEPTH IN SUB. CRITQ IS NOT POSITIVE')
 56   FORMAT(' *ERR:562* SOLUTION DOES NOT EXIST AT DEPTH=',F10.3,
     A /,10X,'CONSTRICTED SECTION IS NOT A CONSTRICTION!')
C***********************************************************************
      IF(Y.LE.0.0) THEN
        WRITE(STDOUT,50)
        EFLAG = 1
        RETURN
      ENDIF
 
      QT = Q
      ET = E
      GRAVT = GRAV
      ADRST = ADRS
      TYPE = ATYPE
      OPUNIT = STDOUT
 
      YD = Y
      YU = ET
      CALL REGFAL
     I           (EPSARG, EPSF, FISE,
     M            YD, YU,
     O            Y, FLAG)
 
      IF(FLAG.GE.1) THEN
C       SOLUTION DOES NOT EXIST
        WRITE(STDOUT,56) Y
        EFLAG = 1
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   CRITQ
     I                  (GRAV, STDIN, STDOUT, STDTAB,
     M                   EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute critical flow table for flow through a constriction
C     when the velocity head of the approaching flow is included
C     in the analysis.  A user supplied discharge coefficient is
C     used to account for contraction losses and friction losses.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      REAL GRAV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER NOQC(6)
      SAVE NOQC
 
C     + + + LOCAL VARIABLES + + +
      INTEGER APPTAB, ATYPE, CONTAB, CTYPE, I, IOFF, J, N, NDEPTH,
     A        TABLE, WFLAG, XOFF, type
      REAL AC, ALPHA, CD, EA, QAVEC(PMXPNT), QC, TC, YA, YAVEC(PMXPNT),
     A     YC, ZA, ZC, factor
      real*8 easting, northing
      CHARACTER BETOPT*8, CIN*50, LABEL*50, LINE*80, MONTON*8, NAME*8,
     A          OUTOPT*8, SAVOPT*8, TABID*16, APPTABID*16, CONTABID*16,
     b          body_head*30, zone*8, hgrid*8, vdatum*8, unitsys*8,
     c          basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      EXTERNAL CHKCFC, inline, INVTSE, KIL, SETOPT, LENSTR,
     A         TAB_IN_USE, READ_TABID_PLUS, READ_TABID
 
C     + + + DATA INITIALIZATIONS + + +
      DATA NOQC/1,1,0,1,1,0/
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:559* COEFFICIENT OF DISCHARGE <= 0.0 OR > 1.0')
 52   FORMAT('TABID= ',A)
 54   FORMAT(/,' TABID= ',A)
 56   FORMAT('TYPE=    2')
 58   FORMAT(' TYPE=    2')
 60   FORMAT('REFL=0.0')
 62   FORMAT(' REFL=0.0')
 64   FORMAT('     DEPTH DISCHARGE',2X,A50)
 66   FORMAT('      DEPTH DISCHARGE',2X,A50)
 68   FORMAT(F10.3,F10.1)
 70   FORMAT(1X,F10.3,F10.1)
 72   FORMAT(' *ERR:560* APPROACH SECTION TABLE DOES NOT ',
     A     ' REPRESENT',/,' A CROSS SECTION.')
 74   FORMAT(' *ERR:561* CONTRICTED SECTION TABLE DOES NOT ',
     B     ' REPRESENT',/,10X,' A CROSS SECTION.')
 76   FORMAT(/,' TESTING TabId=',A,' FOR MONOTONE CELERITY AND ',
     A       'CRITICAL FLOW.')
 78   FORMAT(' *ERR:565* BOTTOM ELEV. OF APPROACH SECTION ABOVE ',
     A  'BOTTOM ELEV.',/,10X,' OF CONSTRICTED SECTION.')
 80   FORMAT(' *ERR:566* INSUFFICIENT SPACE IN ITAB/FTAB TO SAVE TABLE',
     A /,10X,'IN CRITQ. NEED',I6,' MORE ELEMENTS.')
 82   FORMAT(' *ERR:540* TabId=',A,' NOT FOUND')
 84   FORMAT(/,' Approach Section  Constricted Section',/,
     A         ' Elevation  Depth  Elevation     Depth      Flow')
 85   FORMAT(F10.2,F7.2,F11.2,F10.2,F10.1)
 86   FORMAT(/,' TabId= ',A,2X,'Internal number=',I6,2X,A)
C***********************************************************************
      body_head = '      DEPTH DISCHARGE'
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID_PLUS
     I                   (STDOUT, LINE,
     O                    EFLAG, TABID, TABLE, CIN)
 
      WRITE(STDOUT,86) TABID(1:LENSTR(TABID)), TABLE, CIN

      CALL SETOPT
     I           (STDOUT, CIN,
     O            SAVOPT, OUTOPT, MONTON, BETOPT)
 
C     MAKE SURE TABLE NUMBER IS NOT ALREADY USED IN THIS INPUT
 
      IF(FTPNT(TABLE).NE.0) CALL  TAB_IN_USE
     I                                      (TABID,
     M                                       EFLAG)
 
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

C     INPUT THE TABLE NUMBER FOR THE APPROACH CROSS SECTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'APPTAB',
     O                EFLAG, APPTABID, APPTAB)
      WRITE(STDOUT,'(1X,A,A)') 'APPTAB=', APPTABID
 
      IF(FTPNT(APPTAB).LE.0) THEN
        WRITE(STDOUT,82)  APPTABID
        EFLAG = 1
      ENDIF
 
C     INPUT THE TABLE NUMBER FOR THE CONSTRICTED CROSS SECTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'CONTAB',
     O                EFLAG, CONTABID, CONTAB)
      WRITE(STDOUT,'(1X,A,A)') 'CONTAB=', CONTABID
 

      IF(FTPNT(CONTAB).EQ.0) THEN
        WRITE(STDOUT,82)  CONTABID
        EFLAG = 1
      ENDIF
 
C     INPUT THE DISCHARGE COEFFICIENT FOR CONTRACTION AND APPROACH SECTION
C     LOSSES
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(A8,F10.0)',ERR=991) NAME, CD
      WRITE(STDOUT,'(1X,A8,F10.3)') NAME, CD
      IF(CD.GT.1.0.OR.CD.LE.0.0) THEN
        WRITE(STDOUT,50)
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(A50)',ERR=991) LABEL
      WRITE(STDOUT,'(1X,A50)') LABEL
 
      IF(EFLAG.GT.0) RETURN
 
C     CONVERT THE TABLE NUMBERS TO TABLE ADDRESSES FOR CONVENIENCE
C     TABLE NUMBERS ARE RETAINED IN THE TABLES FOR ERROR REPORTING
 
      APPTAB = FTPNT(APPTAB)
      CONTAB = FTPNT(CONTAB)
 
C     ENSURE THAT THE TABLES ARE OF THE CORRECT TYPE.
 
      ATYPE = ITAB(APPTAB+2)
      CTYPE = ITAB(CONTAB+2)
      IF(ATYPE.LT.20.OR.ATYPE.GT.25) THEN
        WRITE(STDOUT,72)
        EFLAG = 1
      ENDIF
      IF(CTYPE.LT.20.OR.CTYPE.GT.25) THEN
        WRITE(STDOUT,74)
        EFLAG = 1
      ENDIF
 
C     CHECK FOR NON-MONOTONE CELERITY AND CRITICAL FLOW
 
      IF(EFLAG.GT.0) RETURN
 
      WRITE(STDOUT,76) CONTABID(1:LENSTR(CONTABID))
      CALL CHKCFC
     I           (GRAV, STDOUT, CONTAB,
     O            WFLAG)
 
C     GET THE BOTTOM ELEVATIONS FOR EACH CROSS SECTION
 
      ZA = FTAB(APPTAB + 5)
      ZC = FTAB(CONTAB + 5)
 
      IF(ZA.GT.ZC) THEN
        WRITE(STDOUT,78)
        EFLAG = 1
      ENDIF
 
      IF(EFLAG.GT.0) RETURN
 
C     BASIC DATA HAS BEEN INPUT IN VALID FORM.  NOW FOR EACH POSITIVE
C     DEPTH VALUE IN THE CONSTRICTED SECTION TABLE, COMPUTE THE
C     CRITICAL FLOW AND THE CORRESPONDING DEPTH IN THE APPROACH
C     CROSS SECTION
 
      NDEPTH = 1
c     Set the point at zero depth in approach section. Always present.
      yavec(ndepth) = 0.0
      qavec(ndepth) = 0.0

c     Check if zero flow persists for non-zero depth in approach section.
      IF(ZC.GT.ZA) THEN
        ndepth = ndepth + 1
        yavec(ndepth) = zc -za
        qavec(ndepth) = 0.0
      ENDIF

 
C     IOFF GIVES THE OFFSET FROM THE TABLE ADDRESS TO THE FIRST
C     POSITIVE DEPTH  ENTRY IN THE TABLE.
 
      XOFF = OFFVEC(CTYPE)
      IOFF = XTIOFF + XOFF
 
      WRITE(STDOUT,84)
 
      WRITE(STDOUT,85) ZA + ZC - ZA, ZC - ZA, ZC, 0.0, 0.0
 100  CONTINUE
 
        YC = FTAB(CONTAB + IOFF)
        TC = FTAB(CONTAB + IOFF + 1)
        AC = FTAB(CONTAB + IOFF + 2)
        IF(NOQC(CTYPE-19).EQ.1) THEN
C         TABLE DOES NOT HAVE CRITICAL FLOW.  THEREFORE COMPUTE
C         IT FROM THE SIMPLE FORMULA ASSUMMING UNIFORM VELOCITY
C         DISTRIBUTION
 
          QC = AC*SQRT(GRAV*AC/TC)
          ALPHA = 1.0
        ELSE
C         TABLE DOES HAVE QC.  GET IT
          QC = FTAB(CONTAB + IOFF + 7)
          ALPHA = FTAB(CONTAB + IOFF + 6)
        ENDIF
 
C       NOW COMPUTE THE SPECIFIC ENERGY WHICH MUST EXIST IN THE
C       APPROACH SECTION IN ORDER TO PRODUCE THE CRITICAL FLOW IN
C       THE CONSTRICTED SECTION
 
        EA = ZC + YC + 0.5*ALPHA*QC**2/(GRAV*(AC**2)*(CD**2)) - ZA
 
C       NOW FIND THE DEPTH IN THE APPROACH SECTION WHICH CORRESPONDS
C       TO THE SPECIFIC ENERGY, EA, AND THE FLOWRATE, QC.
C       INVTSE- INVerT Specific Energy
 
        YA = ZC + YC - ZA
 
        CALL INVTSE
     I             (GRAV, STDOUT, APPTAB, ATYPE, QC, EA, YA,
     O              EFLAG)
 
        NDEPTH = NDEPTH + 1
 
        YAVEC(NDEPTH) = YA
        QAVEC(NDEPTH) = QC
 
        WRITE(STDOUT,85) ZA + YA, YA, ZC + YC, YC, QC
C       INCREMENT TO THE NEXT LEVEL IN THE CROSS SECTION TABLE
 
        IOFF = IOFF + XOFF
 
        IF(CONTAB + IOFF.LE.ITAB(CONTAB)) GOTO 100
 
 
C     TABLE IS COMPLETE. OUTPUT THE TABLE TO BOTH STDOUT AND
C     STDTAB
 
      IF(SAVOPT(1:4).EQ.'SAVE') THEN
C       CHECK ON SPACE IN THE TABLE SYSTEM
 
        N = OFF234 + 2*(2 + NDEPTH)
        IF((N + FTP).GT.MRFTAB) THEN
          WRITE(STDOUT,80) N + FTP - MRFTAB
          EFLAG = 1
        ENDIF
 
        
        FTPNT(TABLE) = FTP
        ITAB(FTP+1) = TABLE
        ITAB(FTP+2) = 2
        ITAB(FTP+3) = FTP + OFF234
        FTAB(FTP+4) = 0.0
        FTAB(FTP+5) = 1.0
        itab(ftp+18) = 0     !datum is not used.
        ftab(ftp+19) = 0.0
        J = ITAB(FTP+3)
      ENDIF
c      WRITE(STDTAB,52) TABID(1:LENSTR(TABID))
      WRITE(STDOUT,54) TABID(1:LENSTR(TABID))
 
c      WRITE(STDTAB,56)
      WRITE(STDOUT,58)
 
c      WRITE(STDTAB,60)
      WRITE(STDOUT,62)
 
c      WRITE(STDTAB,64) LABEL
      WRITE(STDOUT,66) LABEL
 
c      WRITE(STDTAB,68) 0.0, 0.0
      WRITE(STDOUT,70) 0.0, 0.0
      IF(ZC.GT.ZA) THEN
c        WRITE(STDTAB,68) ZC - ZA, 0.0
        WRITE(STDOUT,70) ZC - ZA, 0.0
      ENDIF
 
c      IF(SAVOPT(1:4).EQ.'SAVE') THEN
c        FTAB(J) = 0.
c        FTAB(J+1) = 0.0
c        J = J + 2
c        IF(ZC.GT.ZA) THEN
c          FTAB(J) = ZC - ZA
c          FTAB(J+1) = 0.0
c          J = J + 2
c        ENDIF
c      ENDIF
 
      DO 2300 I=1,NDEPTH
c        WRITE(STDTAB,68) YAVEC(I), QAVEC(I)
        WRITE(STDOUT,70) YAVEC(I), QAVEC(I)
        IF(SAVOPT(1:4).EQ.'SAVE') THEN
          FTAB(J) = YAVEC(I)
          FTAB(J+1) = QAVEC(I)
          J = J + 2
        ENDIF
 2300 CONTINUE
 
c      WRITE(STDTAB,68) -1.0, 0.0

      type = 2
      factor = 1.0
      CALL OUTPUT_TYPE_234
     I                 (STDOUT, STDTAB, ndepth, table, TYPE, FACTOR, 
     I                  BODY_HEAD, LABEL, yavec, qavec, qavec,
     i                  zone, hgrid, vdatum, unitsys, basis,
     i                  easting, northing)

 
        IF(SAVOPT(1:4).EQ.'SAVE') THEN
          ITAB(FTP) = J - 2
          FTP = J
        ENDIF
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
