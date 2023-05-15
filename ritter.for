C
C
C
      REAL FUNCTION   FRIT
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute the generalized ritter residual.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ritcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER NTAB
      REAL C, DT, PDV, T
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTAB, XLOOKW
C***********************************************************************
C     FIND THE ESCOFFIER VARIABLE AND THE AREA
 
 
      CALL XLOOKW
     I           (WTAB, OFF, Y,
     O            A, T, DT, C, W)
 
C     FIND THE FLOW THROUGH THE CONSTRICTION
 
      CALL LKTAB
     I          (QTAB, Y, 0,
     O           Q, NTAB, PDV)
 
C     COMPUTE THE RESIDUAL
 
      IF(A.LE.0.0) A = 1.0
      FRIT = ((W1 + V1 - W) -  Q/A)/(W1 + V1)

C      WRITE(RITOUT,*) ' Y=',Y,' FRIT=',FRIT 
      RETURN
      END
C
C
C
      SUBROUTINE   FNDRIT
     I                   (STDOUT, APPTAB, XOFF, CONTAB, Y1, Q1,
     O                    YPEAK, QPEAK)
 
C     + + + PURPOSE + + +
C     Use modified regula falsi to find the generalized Ritter
C     flood peak for the given conditions.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER APPTAB, CONTAB, STDOUT, XOFF
      REAL Q1, QPEAK, Y1, YPEAK
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     APPTAB - Address of table describing reservoir cross section
C     XOFF   - Offset between successive depth values for cross section
C               function table
C     CONTAB - Address of the table giving the flow through the
C               failure opening
C     Y1     - Initial maximum depth in reservoir
C     Q1     - Flowrate at section 1
C     YPEAK  - Depth at peak outflow from reservoir
C     QPEAK  - Peak flow rate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ritcom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL A1, C1, DT1, T1, YD, YU
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FRIT, REGFAL, XLOOKW
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT('0*BUG:504* NO SOLUTION FOR COMMAND GRITTER.',/,
     A  10X,' UPPER DEPTH=',F10.3,' LOWER DEPTH=',F10.3)
C***********************************************************************
      WTAB = APPTAB
      QTAB = CONTAB
      OFF = XOFF
      RITOUT = STDOUT
 
C     FIND THE VALUES AT THE INITIAL DEPTH IN THE RESERVOIR
 
      CALL XLOOKW
     I           (APPTAB, XOFF, Y1,
     O            A1, T1, DT1, C1, W1)
 
      V1 = Q1/A1
 
C     IF THE CRITICAL FLOW RELATIONSHIP IS MONOTONE INCREASING
C     WITH DEPTH, THEN THERE SHOULD BE ONLY ONE POSITIVE ROOT
C     AND IT SHOULD BE BETWEEN 0 AND Y1
 
 
      YD = 0.2*Y1
      YU = Y1
 
      CALL REGFAL
     I           (EPSARG, EPSF, FRIT,
     M            YD, YU,
     O            YPEAK, FLAG)
      IF(FLAG.GE.1) THEN
C       NO SOLUTION
        WRITE(STDOUT,50) YD, YU
        RETURN
      ENDIF
      QPEAK = Q
      RETURN
 
      END
C
C
C
      SUBROUTINE   RITTER
     I                   (GRAV, STDIN, STDOUT,
     M                    EFLAG)
 
C     + + + PURPOSE + + +
C     Compute a generalized ritter solution for the instantaneous
C     dam break peak flow.  The reservoir is still  assumed to be prismatic,
C     horizontal, and frictionless.  however, the cross section
C     is non-rectangular and the failure need not be complete.
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT
      REAL GRAV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER APPTAB, ATYPE, CONTAB, HA, IOFF, J, WFLAG, XOFF
      REAL A, C, COLD, Q1, QPEAK, T, W, WOLD, Y, Y1, YOLD, YPEAK
      CHARACTER LABEL*79, LINE*80, NAME*8, APPTABID*16, CONTABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKCFC, FNDRIT, inline, KIL, READ_TABID
 
C     + + + INPUT FORMATS + + +
 2    FORMAT(2F10.0)
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT(/,'      DEPTH       TOP      AREA  CELERITY ESCOFFIER')
 53   FORMAT(/,' TABLE FOR CELERITY AND ESCOFFIER VARIABLE FOR ',
     A       'RESERVOIR')
 54   FORMAT(1X,F10.3,F10.1,1PE10.3,0PF10.3,1PE10.3)
 56   FORMAT(1X,F10.3,F10.1,F10.3,1PE12.5)
 58   FORMAT(/,'  INITDEPTH  INITFLOW PEAKDEPTH    PEAKFLOW',A)
 72   FORMAT(' *ERR:560* Approach section table id does not ',
     A     ' represent',/,' a cross section.')
 74   FORMAT(' *ERR:563* Constricted flow table id is invalid.',
     A      /,10X,'Type must be 2.')
 76   FORMAT(/,' Testing TABID= ',A,' for monotone celerity and ',
     A       'critical flow.')
 78   FORMAT(' *ERR:540* TABID= ',A,' not found')
 82   FORMAT(/,' *ERR:682* Cross section type=',I3,' not supported',
     A      ' in GRITTER. ',/,10X,' Use types: 20, 21, 23, or 24.')
C***********************************************************************
C     INPUT THE TABLE ID FOR THE APPROACH CROSS SECTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'APPTAB',
     O                EFLAG, APPTABID, APPTAB)
C      READ(LINE,'(A8,A)',ERR=991) NAME, APPTABID
      WRITE(STDOUT,'(1X,A,A)') 'APPTAB= ', APPTABID
 
      IF(FTPNT(APPTAB).LE.0) THEN
        WRITE(STDOUT,78) APPTABID
        EFLAG = 1
      ENDIF
 
C     INPUT THE TABLE ID FOR THE CONSTRICTED FLOW
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'CONTAB',
     O                EFLAG, CONTABID, CONTAB)
      WRITE(STDOUT,'(1X,A,A)') 'CONTAB= ', CONTABID
 
      IF(FTPNT(CONTAB).LE.0) THEN
        WRITE(STDOUT,78) CONTABID
        EFLAG = 1
      ENDIF
 
C     INPUT THE BREACH OFFSET VALUE
 
C      CALL inline(STDIN, STDOUT, LINE)
C      READ(LINE,'(A8,F10.0)',ERR=991) NAME, Z0
C      WRITE(STDOUT,'(1X,A8,F10.3)') NAME, Z0
C      IF(Z0.LT.0.0) THEN
C        WRITE(STDOUT,50)
C        EFLAG = 1
C      ENDIF
 
      IF(EFLAG.GT.0) RETURN
 
C     CONVERT THE TABLE NUMBERS TO TABLE ADDRESSES FOR CONVENIENCE
C     TABLE NUMBERS ARE RETAINED IN THE TABLES FOR ERROR REPORTING
 
      APPTAB = FTPNT(APPTAB)
      CONTAB = FTPNT(CONTAB)
 
C     ENSURE THAT THE TABLES ARE OF THE CORRECT TYPE.
 
      ATYPE = ITAB(APPTAB+2)
      IF(ATYPE.LT.20.OR.ATYPE.GT.25) THEN
        WRITE(STDOUT,72)
        EFLAG = 1
      ENDIF
 
      IF(ATYPE.EQ.22.OR.ATYPE.EQ.25) THEN
        WRITE(STDOUT,82) ATYPE
        EFLAG = 1
      ENDIF
      IF(ITAB(CONTAB+2).NE.2) THEN
        WRITE(STDOUT,74)
        EFLAG = 1
      ENDIF
 
C     CHECK FOR NON-MONOTONE CELERITY AND CRITICAL FLOW
 
      IF(EFLAG.GT.0) RETURN
      WRITE(STDOUT,76) APPTABID
      CALL CHKCFC
     I           (GRAV, STDOUT, APPTAB,
     O            WFLAG)
 
 
      IF(EFLAG.GT.0) RETURN
 
C     SET THE OFFSET FOR THE CROSS SECTION TABLE
      XOFF = OFFVEC(ATYPE)
 
C     SET THE HIGH ADDRESS OF THE APPROACH TABLE
      HA = ITAB(APPTAB)
 
C     BASIC DATA HAS BEEN INPUT IN VALID FORM.
C     NOW MODIFY THE APPROACH CROSS SECTION TABLE SO THAT
C     IT WILL REPRESENT THE CELERITY AND THE ESCOFFIER STAGE VARIABLE
C     IN THE RESERVOIR. WE WILL REPLACE THE SQRT OF CONVEYANCE AND
C     THE MOMENTUM FLUX CORRECTION FACTOR IN THE TABLE.
 
C     WRITE HEADING FOR THE ESCOFFIER TABLE PRINTOUT
 
      WRITE(STDOUT,53)
      WRITE(STDOUT,52)
 
C     SET THE TYPE TO ZERO TO MAKE IT UNUSABLE FOR OTHER OPERATIONS
 
      ITAB(APPTAB+2) = 0
 
C     CLEAR THE ZERO DEPTH ENTRIES
 
      FTAB(APPTAB+XTIOFF+3) = 0.0
      FTAB(APPTAB+XTIOFF+4) = 0.0
 
      WRITE(STDOUT,54) (FTAB(APPTAB+XTIOFF+J),J=0,4)
 
C     IOFF GIVES THE OFFSET FROM THE TABLE ADDRESS TO THE FIRST
C     POSITIVE DEPTH  ENTRY IN THE TABLE.
 
 
      IOFF = XTIOFF + XOFF
 
      COLD = 0.0
      YOLD = 0.0
      WOLD = 0.0
 
 100  CONTINUE
 
        Y = FTAB(APPTAB + IOFF)
        T = FTAB(APPTAB + IOFF + 1)
        A = FTAB(APPTAB + IOFF + 2)
        C = SQRT(GRAV*A/T)
 
        W =  WOLD + 2.*GRAV*(Y - YOLD)/(C + COLD)
        WRITE(STDOUT,54) Y, T, A, C, W
 
        COLD = C
        YOLD = Y
        WOLD = W
 
        FTAB(APPTAB + IOFF + 3) = C
        FTAB(APPTAB + IOFF + 4) = W
 
C       INCREMENT TO THE NEXT LEVEL IN THE CROSS SECTION TABLE
 
        IOFF = IOFF + XOFF
 
        IF(APPTAB + IOFF.LE.HA) GOTO 100
 
C     THE ESCOFFIER TABLE IS COMPLETE.  NOW COMPUTE THE PEAK FLOW
C     FOR EACH CASE GIVEN IN THE INPUT.
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(A79)',ERR=991) LABEL
C      WRITE(STDOUT,'(1X,A79)') LABEL
      WRITE(STDOUT,58) LABEL(21:79)
 
 200  CONTINUE
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
        READ(LINE,2,ERR=991) Y1, Q1
        IF(Y1.LE.0.0) RETURN
 
C        IF(Y1.LE.Z0) THEN
C          WRITE(STDOUT,80)
C          EFLAG =1
C        ENDIF
 
        CALL FNDRIT
     I             (STDOUT, APPTAB, XOFF, CONTAB, Y1, Q1,
     O              YPEAK, QPEAK)
 
 
        WRITE(STDOUT,56) Y1, Q1, YPEAK, QPEAK
 
        GOTO 200
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
