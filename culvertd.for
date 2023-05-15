C
C
C
      REAL FUNCTION   NDRSD
     I                     (Y)
 
C     + + + PURPOSE + + +
C     Find the normal flow residual at depth Y in the cross
C     section  with the flow given in the labeled common block
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'ndrsdc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL K, YLOC
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTK
C***********************************************************************
      YLOC = Y
      CALL LKTK
     I         (ADR,
     M          YLOC,
     O          K)
 
      NDRSD = (RTSBOT*K - FLOW)/FLOW
 
C      WRITE(STD6,51) Y, K, FLOW, NDRSD
 
      RETURN
      END
C
C
C
      SUBROUTINE   FNDND
     I                  (STDOUT, ADRS, Q, SBOT,
     M                   D,
     O                   YE)
 
C     + + + PURPOSE + + +
C     Find normal depth for the given bottom slope
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, STDOUT
      REAL D, Q, SBOT, YE
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     ADRS   - Address of function table
C     Q      - Flowrate
C     SBOT   - Bottom slope for normal depth
C     D      - Maximum vertical extent of a closed conduit
C     YE     - Normal depth
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'ndrsdc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL FL, FR, K, QMAX, Y, YL, YMAX, YR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL FMXARG, NDRSD
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FMXARG, GETTBN, LKTK, NDRSD, REGFLT, LENSTR, GET_TABID
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:625* TABID=',A,' overflow seeking norm. depth for',
     A       ' flow=',F10.2,11X,' Table MaxArg=',F10.2)
 52   FORMAT(' *ERR:626* TABID=',A,' underflow seeking norm. depth for',
     A       ' flow=',F10.3)
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' FNDND.')
 60   FORMAT(' *BUG:XXX* REGFLT ITERATION>100 IN FNDND')
C***********************************************************************
C     SET THE VALUES FOR THE FUNCTION, NDRSD
      ADR = ADRS
      FLOW = Q
      RTSBOT = SQRT(SBOT)
      CALL LKTK
     I         (ADR,
     M          D,
     O          K)
C     COMPUTE THE MAXIMUM FLOW IN THE CONDUIT FOR THE GIVEN SLOPE
      QMAX = K*RTSBOT
 
C     GET THE MAXIMUM DEPTH IN THE TABLE
      YMAX = FMXARG(ADRS)
      IF(Q.GE.QMAX) THEN
        YE = YMAX
C       WRITE(STDOUT,*) ' NORMAL DEPTH=',YMAX,' D=',D,' K=',K
        RETURN
      ENDIF
 
C     SEARCH FOR A POSITIVE RESIDUAL
      YR = YE
      YL = -1.0
 100  CONTINUE
        FR = NDRSD(YR)
        IF(FR.GE.0.0) THEN
          GOTO 110
        ELSE
          YL = YR
          FL = FR
          IF(ABS(FR).LT.EPSF) THEN
            YR = 1.01*YR
            GOTO 100
          ELSE
            YR = 0.7*YR + 0.3*D
            IF(ABS(YR - YMAX).LE.EPSARG) THEN
              TABID = GET_TABID(GETTBN(ADRS))
              WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), Q, YMAX
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
C     POSITIVE RESIDUAL FOUND- SEARCH FOR NEGATIVE RESIDUAL
 
      IF(YL.LT.0.0) THEN
        YL = 0.6*YR
 120    CONTINUE
          FL = NDRSD(YL)
          IF(FL.LE.0.0) THEN
            GOTO 130
          ELSE
            FR = FL
            YR = YL
            YL = 0.6*YL
            IF(ABS(YL).LT.1.E-6) THEN
              TABID = GET_TABID(GETTBN(ADRS))
              WRITE(STDOUT, 52) TABID(1:LENSTR(TABID)), Q
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF
 
C     WE HAVE A SIGN CHANGE OR ONE OR BOTH POINTS HAVE ZERO RESIDUAL
      CALL REGFLT
     I           (EPSARG, EPSF, NDRSD,
     M            YL, YR, FL, FR,
     O            Y, FLAG)
 
C      WRITE(STDOUT, 56) YL, Y, YR
C      WRITE(STDOUT, 62) FL, FR
 
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      YE = Y
 
      RETURN
      END
C
C
C
      REAL FUNCTION   CRQRE
     I                     (Y)
 
C     + + + PURPOSE + + +
C     Find the critical flow residual at depth Y in the cross
C     section and with the flow given in the labeled common block
C     using energy principles
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'crqrec.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL QC, YLOC
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTQC
C***********************************************************************
      YLOC = Y
      CALL LKTQC
     I          (ADR,
     M           YLOC,
     O           QC)
 
      CRQRE = (QC - FLOW)/QC
 
C      WRITE(STD6,51) FLOW, QC, Y, CRQRE
 
      RETURN
      END
C
C
C
      SUBROUTINE   FNDCDE
     I                   (STDOUT, ADRS, Q,
     M                    YE)
 
C     + + + PURPOSE + + +
C     Find critical depth using energy relationship in the cross
C     section at ADRS at the flow Q.  YE contains an estimate
C     of the critical depth.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, STDOUT
      REAL Q, YE
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     ADRS   - Address of function table
C     Q      - Flowrate
C     YE     - Critical depth
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'crqrec.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL EPF, FL, FR, Y, YL, YMAX, YR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL CRQRE, FMXARG
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CRQRE, FMXARG, GETTBN, RGF5, LENSTR, GET_TABID
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:595* TABID=',A,' overflow seeking crit. depth for',
     A       ' flow=',F10.2,11X,' Table MaxArg=',F10.2)
 52   FORMAT(' *ERR:596* TABID=',A,' underflow seeking crit. depth for',
     A       ' flow=',F10.3)
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF5 CLAIMS NONE IN',
     A       ' FNDCDE.')
 60   FORMAT(' *BUG:XXX* RGF5 TAKES MORE THAN 100 ITERATIONS.')
C***********************************************************************
C     SET THE VALUES FOR THE FUNCTION, CRQRE
      ADR = ADRS
      FLOW = Q
      EF = 0
 
      EPF = EPSF
C     GET THE MAXIMUM DEPTH IN THE TABLE
      YMAX = FMXARG(ADRS)
 
C     SEARCH FOR A POSITIVE RESIDUAL
      YR = YE
      YL = -1
 100  CONTINUE
        FR = CRQRE(YR)
        IF(FR.GE.0.0) THEN
          GOTO 110
        ELSE
          YL = YR
          FL = FR
          IF(ABS(FR).LT.EPSF) THEN
            YR = 1.01*YR
            GOTO 100
          ELSE
            YR = 0.9*YR + 0.1*YMAX
            IF(ABS(YR - YMAX)/YMAX.LE.1.E-6) THEN
              TABID = GET_TABID(GETTBN(ADRS))
              WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), Q, YMAX
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
C     POSITIVE RESIDUAL FOUND- SEARCH FOR NEGATIVE RESIDUAL
 
      IF(YL.LT.0.0) THEN
        YL = 0.5*YE
 120    CONTINUE
          FL = CRQRE(YL)
          IF(FL.LE.0.0) THEN
            GOTO 130
          ELSE
            FR = FL
            YR = YL
            YL = 0.5*YL
            IF(ABS(YL).LT.EPSABS) THEN
              TABID = GET_TABID(GETTBN(ADRS))
              WRITE(STDOUT, 52) TABID(1:LENSTR(TABID)), Q
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF
C     WE HAVE A SIGN CHANGE OR ONE OR BOTH POINTS HAVE ZERO RESIDUAL
C     FORCE CONVERGENCE ON RESIDUAL ONLY.
 
C      WRITE(STDOUT,*) ' FNDCDE: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL RGF5
     I         (0.0, EPF, CRQRE,
     M          YL, YR, FL, FR,
     O          Y, FLAG)
 
C      WRITE(STDOUT, 56) YL, Y, YR, Q
 
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      YE = Y
 
      RETURN
      END
C
C
C
      REAL FUNCTION   INVSE
     I                     (Q, E, Y)
 
C     + + + PURPOSE + + +
C     Compute the inverse of the unit width specific energy seeking
C     a subcritical solution.  INVSE assumes that a subcritical solution
C     exists.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL E, Q, Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Q      - Flowrate
C     E      - Specific energy value.
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'grvcom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL DY, EMIN, F, FP, QSQR, VSQR, YC, YT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' More than 100 iterations in INVSE. Q=',F13.3,
     A   ' E=',F10.4,/,5X,'Y=',F10.4,' F=',F10.4,' DY=',F10.3,
     B   ' YC=',F10.4)
C***********************************************************************
C      WRITE(STD6,*) ' INVSE: Q=',Q,' E=',E,' Y=',Y
 
      QSQR = Q*Q
      YC = (QSQR/GRAV)**0.3333333
      EMIN = YC + QSQR/(GRAV2*YC*YC)
      IF(E.LT.EMIN) THEN
C       No solution exists.  Return critical depth.
        INVSE = YC
        RETURN
      ENDIF
      IF(ABS(E-EMIN)/E.LT.1.E-3) THEN
C       Take critical depth as the solution
        INVSE = YC
        RETURN
      ENDIF
      YT = Y
      KNT = 0
 100  CONTINUE
        IF(YT.EQ.0.0) THEN
          WRITE(STD6,*) ' Y=0 IN INVSE. Q=',Q,' E=',E
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        VSQR = QSQR/YT**2
        F = YT + VSQR/GRAV2 - E
        IF(ABS(F/E).LT.5.E-4) THEN
          INVSE = YT
          RETURN
        ELSEIF(ABS(YT - YC)/YC.LT.0.001) THEN
          INVSE = YT
          RETURN
        ENDIF
        FP = 1.0 - VSQR/(YT*GRAV)
        DY = -F/FP
        YT = YT + DY
        IF(YT.LT.YC) THEN
          YT = YT - DY
          YT = 0.5*(YT + YC)
        ELSE
          IF(ABS(DY/YT).LT.1.E-4) THEN
            INVSE = YT
            RETURN
          ENDIF
        ENDIF
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          WRITE(STD6,50) Q, E, YT, F, DY, YC
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        GOTO 100
      END
C
C
C
      REAL FUNCTION   DEGCON
     I                      (C123, A1, A)
 
C     + + + PURPOSE + + +
C     Make final adjustment to the coef for discharge for types 1, 2, and
C     3 for the degree of contraction
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL A, A1, C123
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     C123   - Discharge coefficient for culvert flow types 1, 2, and 3
C     A1     - Area of approach section
C     A      - Area of control section for culvert
 
C     + + + LOCAL VARIABLES + + +
      REAL M
C***********************************************************************
C      WRITE(STD6,*) ' DEGCON: A1=',A1,' A=',A
      M = 1.0 - A/A1
C      WRITE(STD6,*) ' DEGCON: M=',M
      IF(M.LT.0.0) M = 0.0
      IF(M.GT.0.80) THEN
        DEGCON = C123
      ELSE
        DEGCON = 0.98 - (0.98 - C123)*M/0.80
      ENDIF
C      WRITE(STD6,*) ' DEGCON=',DEGCON
      RETURN
      END
C
C
C
      SUBROUTINE   DISLSS
     I                   (IU, ID, IAT3D, IAT6D, DUP, TY6LSS, XVEC,
     O                    SEVEC)
 
C     + + + PURPOSE + + +
C     Distribute the losses for computation of a type 5 profile
C     in a culvert.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER IAT3D, IAT6D, ID, IU
      REAL DUP, SEVEC(ID), TY6LSS, XVEC(ID)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     IAT3D  - index to vena contracta location
C     IAT6D  - index to vena contracta location
C     DUP    - vertical diameter of culvert barrel at upstream end
C     TY6LSS - estimated type 6 loss
C     XVEC   - Stations for nodes on a branch or along a culvert
C     SEVEC  - vector or slopes that represent eddy losses
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL S
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
C***********************************************************************
      DO 122 I=IU, ID
        SEVEC(I) = 0.0
 122  CONTINUE
 
      IF(IAT6D.GT.0) THEN
C       Spread the losses over an approximate 3*DUP distance
        S = TY6LSS/(ABS(XVEC(IAT3D) - XVEC(IAT6D)))
        DO 125 I=IAT3D+1,IAT6D
          SEVEC(I) = S
 125    CONTINUE
      ELSE
C       Spread the losses over 3*DUP distance but allow reduction of
C       loss because the culvert is so short.
        S = TY6LSS/(3.*DUP)
        DO 127 I=IAT3D+1,ID
          SEVEC(I) = S
 127    CONTINUE
      ENDIF
      RETURN
      END
C
C
C
      REAL FUNCTION   TY6RAT
     I                      (Q, CLASS, D, AEXIT)
 
C     + + + PURPOSE + + +
C     Compute the piezometric level ratio for a type 6 outlet.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL AEXIT, D, Q
      CHARACTER CLASS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Q      - Flowrate
C     CLASS  - Class for culvert shape, BOX, PIPE, ..
C     D      - Maximum vertical extent of a closed conduit
C     AEXIT  - exit area for the culvert
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL NE, QE, RAT, FAC
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        FAC = 1.0
      ELSE
        FAC = 1.8113089
      ENDIF      
      IF(CLASS.EQ.'BOX') THEN
        RAT = Q/(AEXIT*SQRT(GRAV*D))
        IF(RAT.LT.1.0) THEN
          TY6RAT = 1.0 - 0.2*RAT
        ELSEIF(RAT.LT.5.0) THEN
          TY6RAT = 0.875 - 0.075*RAT
        ELSE
          TY6RAT = 0.5
        ENDIF
      ELSE
C       ALL OTHER CLASSES TREATED LIKE A PIPE CULVERT
C       ESTIMATE THE EQUIVALENT NUMBER OF PIPES OF DIAMETER D
C       TO ADJUST THE FLOW RATE
 
        NE = 4.*AEXIT/(3.1416*D**2)
        QE = Q/NE
        RAT = FAC*QE/D**2.5
        IF(RAT.LT.1.0) RAT = 1.0
 
        IF(TY6SUP.EQ.1) THEN
C         Discharge is supported.
          IF(RAT.LE.2.0) THEN
            TY6RAT = 1.0 - 0.1*(RAT - 1.0)
          ELSEIF(RAT.LE.3.0) THEN
            TY6RAT = 1.01 - 0.055*RAT
          ELSEIF(RAT.LE.4.0) THEN
            TY6RAT = 1.061 - 0.0720*RAT
          ELSEIF(RAT.LE.5.0) THEN
            TY6RAT = 1.129 - 0.089*RAT
          ELSEIF(RAT.LE.6.0) THEN
            TY6RAT = 0.899 - 0.043*RAT
          ELSEIF(RAT.LE.7.0) THEN
            TY6RAT = 0.797 - 0.026*RAT
          ELSEIF(RAT.LE.8.0) THEN
            TY6RAT = 0.72 - 0.0152364*RAT
          ELSE
            TY6RAT = 0.5 + 0.52788/RAT**0.80925
          ENDIF
        ELSE
C         Discharge is unsupported.
          IF(RAT.LE.2.0) THEN
            TY6RAT = 1.058 - 0.058*RAT
          ELSEIF(RAT.LE.3.0) THEN
            TY6RAT = 1.056 - 0.057*RAT
          ELSEIF(RAT.LE.4.0) THEN
            TY6RAT = 1.188 - 0.101*RAT
          ELSEIF(RAT.LE.5.0) THEN
            TY6RAT = 1.240 - 0.114*RAT
          ELSEIF(RAT.LE.6.0) THEN
            TY6RAT = 0.940 - 0.054*RAT
          ELSEIF(RAT.LE.7.0) THEN
            TY6RAT = 0.808 - 0.032*RAT
          ELSEIF(RAT.LE.8.0) THEN
            TY6RAT = 0.85 - 0.038*RAT
          ELSE
            TY6RAT = 0.5 + 1.560/RAT**1.546
          ENDIF
        ENDIF
      ENDIF
 
C      WRITE(STD6,*) ' TY6RAT: RAT=',RAT,' Q=',Q,' TY6RAT=',TY6RAT
C      WRITE(STD6,*) ' AEXIT=',AEXIT,' D=',D
      RETURN
      END
C
C
C
      REAL FUNCTION   FNDSZL
     I                      (STDOUT, CULCLS, LOVERD, NBAR, RGHFAC,
     I                       RBVAL, TB15AD, TB16AD)
 
C     + + + PURPOSE + + +
C     Compute the value of the So limit for type 6 flow.  Uses
C     Figures 15 and 16 in the USGS TWRI on culvert flows as
C     a basis.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT, TB15AD
      INTEGER TB16AD(4)
      REAL LOVERD, NBAR, RBVAL, RGHFAC
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     LOVERD - ratio of culvert barrel length to vertical diameter
C     NBAR   - average Manning's n
C     RGHFAC - roughness factor for pipe culverts
C     RBVAL  - rounding/beveling value
C     TB15AD - address of function table for figure 15
C     TB16AD - addresses of function tables for figure 16
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IL, IR
      REAL DFDC, DFDR, F, FL, FR, LOVD, RDL, RDR, ROVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL TDLK10
C***********************************************************************
      LOVD = LOVERD
      IF(LOVD.GT.500) THEN
        LOVD = 500.0
      ENDIF
      IF(CULCLS.EQ.'BOX'.OR.NBAR.LT.0.019) THEN
 
C       Take the barrel to be smooth.
        IF(RBVAL.GT.0.06) THEN
          ROVERD = 0.06
        ELSE
          ROVERD = RBVAL
        ENDIF
        CALL TDLK10
     I             (STDOUT, TB15AD, 10, LOVD, ROVERD,
     O              F, DFDR, DFDC)
        FNDSZL = F
      ELSE
C       Take the barrel to be rough.  Select the pair of tables that
C       contain the value of rounding and beveling.
        IF(RBVAL.GT.0.03) THEN
          ROVERD = 0.03
        ELSE
          ROVERD = RBVAL
        ENDIF
        IF(ROVERD.LE.0.01) THEN
          IR = 2
          RDR = 0.01
        ELSEIF(ROVERD.LT.0.02) THEN
          IR = 3
          RDR = 0.02
        ELSE
          IR = 4
          RDR = 0.03
        ENDIF
        IL = IR - 1
        RDL = RDR - 0.01
 
        CALL TDLK10
     I             (STDOUT, TB16AD(IL), 10, LOVD, RGHFAC,
     O              FL, DFDR, DFDC)
        CALL TDLK10
     I             (STDOUT, TB16AD(IR), 10, LOVD, RGHFAC,
     O              FR, DFDR, DFDC)
 
        FNDSZL = FL +   (ROVERD - RDL)*(FR - FL)/(RDR - RDL)
      ENDIF
 
      RETURN
      END
C
C
C
      REAL FUNCTION   GETD
     I                    (ADRS, STDOUT)
 
C     + + + PURPOSE + + +
C     Find the maximim vertical diameter of the cross section
C     given in the table at address adrs.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, STDOUT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADRS   - Address of function table
C     STDOUT - Fortran unit number for user output and messages
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
      INCLUDE 'grvcom.cmn'
 
C     Called program units
      CHARACTER GET_TABID*16
      EXTERNAL GET_TABID

C     + + + LOCAL VARIABLES + + +
      INTEGER HIGHA, I, LOWA, TYPE, XOFF
      REAL TOP, TOPOLD, DIFF1, DIFF2
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT('0*ERR:590* TabId: ',A,' is not of a',
     A        ' closed conduit.')
 52   FORMAT('0*WRN:529* TabId: ',A,' has slot wddth=',F7.3,
     A       '.',/,11X,'may not be of a closed conduit.')
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        DIFF1 = 0.005
        DIFF2 = 0.15
      ELSE
        DIFF1 = 0.00152
        DIFF2 = 0.04572
      ENDIF
C     GET THE BOUNDING ADDRESSES FOR THE DEPTH VALUES
 
C      WRITE(STD6,*) ' SEARCHING TABLE#=',ITAB(ADRS+1)
      TYPE = ITAB(ADRS+2)
      XOFF = OFFVEC(TYPE)
 
      HIGHA = ITAB(ADRS)
      LOWA = ADRS + XTIOFF
 
C     SCAN THE TOP WIDTH VALUES FROM THE MAXIMUM DEPTH IN THE TABLE
C     DOWNWARD UNTIL THE TOP WIDTH INCREASES BY MORE THAN A PRESET
C     TOLERANCE.
 
      TOPOLD = FTAB(HIGHA+1)
      DO 100 I=HIGHA-XOFF,LOWA,-XOFF
        TOP = FTAB(I+1)
C        WRITE(STD6,*) 'I =',I,' Y=',FTAB(I),' TOP=',TOP,' TOPOLD=',TOPOLD
        IF(TOP.GT.TOPOLD + DIFF1) THEN
C         ASSUME THAT TOP OF CONDUIT IS AT TOPOLD
          IF(TOPOLD.GT.DIFF2) THEN
            WRITE(STDOUT,52) ITAB(ADRS+1), TOPOLD
          ENDIF
          GETD = FTAB(I+XOFF)
          RETURN
        ENDIF
        TOPOLD = TOP
 100  CONTINUE
 
      WRITE(STDOUT,50) GET_TABID(ITAB(ADRS+1))
      GETD = -1.0
      RETURN
      END
C
C
C
      REAL FUNCTION   YOVERD
     I                      (ADR, D, AFULL, ARATIO)
 
C     + + + PURPOSE + + +
C     Find the ratio of depth to vertical height for a closed conduit
C     that corresponds to the ratio of area to full area given by
C     ARATIO.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADR
      REAL AFULL, ARATIO, D
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADR    - Address of function table
C     D      - vertical diameter for the culvert barrel
C     AFULL  - full area of conduit
C     ARATIO - ratio of partial area to full area
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
      INCLUDE 'stdun.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HIGH, KNT, LAST, LOW, XOFF
      REAL A, AL, DT, DY, F, FP, P, PT, T, TL, TR, Y, YL, YR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX No convergence in YOVERD. A=',F8.3,' P=',F8.3)
C***********************************************************************
      A = AFULL*ARATIO
 
      XOFF = OFFVEC(ITAB(ADR+2))
C     Get the address of the argument last accessed in the table.
      LAST = ITAB(ADR+3)
 
      LOW = ADR + XTIOFF
C     Get the address of the highest argument of the table.
      HIGH = ITAB(ADR)
 
C     Adjust addresses to point to the area and not the depth
C     in the table.
 
      LAST = LAST + 2
      LOW = LOW + 2
      HIGH = HIGH + 2
      IF(A.GE.FTAB(LAST)) THEN
 100    CONTINUE
          IF(A.GT.FTAB(LAST+XOFF)) THEN
            LAST = LAST + XOFF
            GOTO 100
          ENDIF
      ELSE
 110    CONTINUE
          LAST = LAST - XOFF
          IF(A.LT.FTAB(LAST)) GOTO 110
      ENDIF
 
C     Save the last access location.
      ITAB(ADR+3) = LAST - 2
 
C     Find the value of depth that would give the area A.
      AL = FTAB(LAST)
      IF(A.EQ.AL) THEN
        Y = FTAB(LAST - 2)
      ELSE
        TL = FTAB(LAST - 1)
        YL = FTAB(LAST - 2)
        TR = FTAB(LAST + XOFF - 1)
        YR = FTAB(LAST + XOFF - 2)
        DT = TR - TL
        DY = YR - YL
        IF(DT.EQ.0.0) THEN
C         Rectangular
          Y = YL + (A - AL)/TL
        ELSEIF(TL.EQ.0.0) THEN
C         Triangular. YL  is zero and AL is zero
          Y = SQRT(2.*YR*A/TR)
        ELSE
C         Trapezoidal
          P = (A - AL)/(DY*TL)
          KNT = 0
 120      CONTINUE
            T = TL + P*DT
            F  = AL + 0.5*P*DY*(T + TL) - A
            FP = DY*T
 
            PT = P -  F/FP
            IF(ABS(F)/A.GT.0.1*EPSF.AND.ABS(P - PT).GT.0.1*EPSF) THEN
              KNT = KNT + 1
              IF(KNT.GT.20) THEN
                WRITE(STD6,50) A, P
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              IF(PT.GT.1.0) THEN
                PT = 0.5*(P + 1.)
              ENDIF
              IF(PT.LT.0.0) THEN
                PT = 0.5*P
              ENDIF
              P = PT
              GOTO 120
            ENDIF
          Y = YL + PT*DY
        ENDIF
      ENDIF
 
      YOVERD = Y/D
      RETURN
      END
C
C
C
      REAL FUNCTION   FCD123
     I                      (STDOUT, TYPE, CULCLS, DUP, Z1TRUE)
 
C     + + + PURPOSE + + +
C     Find the coef of discharge for flow types 1, 2, and 3 making
C     all adjustments except for the degree of contraction adjustment.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT, TYPE
      REAL DUP, Z1TRUE
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     TYPE   - Culvert flow type
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     DUP    - vertical diameter of culvert barrel at upstream end
C     Z1TRUE - known elevation at section 1
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADRTAB
      REAL FR, QC, R, Y
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTQC
 
C     + + + OUTPUT FORMATS + + +
C 50   FORMAT(' *WRN:533* Invalid headwater ratio:',F7.3,
C     A  ' in FCD123 at Z1=',F10.2,/,11X,' Reset to maximum of:',F5.2)
 56   FORMAT(' *BUG:XXX* CULCLS=',A8,' IS NOT SUPPORTED IN FCD123.')
C***********************************************************************
      IF(CULCLS.EQ.'PIPE'.OR.CULCLS.EQ.'FLARED') THEN
        R = (Z1TRUE - ZB2)/DUP
C          WRITE(STDOUT,*) ' R=',R,' Z1TRUE=',Z1TRUE
        IF(R.LE.0.0.OR.R.GT.1.6) THEN
          IF(RATFLG.EQ.0) THEN
            RATFLG = 1
C            WRITE(STDOUT,50) R, Z1TRUE, 1.6
          ENDIF
          R = 1.599
        ENDIF
        IF(R.LT.0.4) THEN
          FCD123 = 0.93
        ELSE
          FCD123 = 0.86773 + R*(0.31564 + R*(-0.48463 + R*(
     A               0.22841 - R*0.04183)))
        ENDIF
C       Apply adjustments for pipe culverts.
        IF(CULCLS.EQ.'PIPE')  FCD123 = KRB*KPROJ*FCD123
      ELSEIF(CULCLS.EQ.'BOX') THEN
        IF(TYPE.EQ.3) THEN
C         SELECT WHICH END TO USE
          IF(A2.LT.A3) THEN
            ADRTAB = ADRXS2
            Y = Y2
          ELSE
            ADRTAB = ADRXS3
            Y = Y3
          ENDIF
C         COMPUTE FROUDE NUMBER
 
          CALL LKTQC
     I              (ADRTAB,
     M               Y,
     O               QC)
          FR = Q3/QC
          IF(FR.GT.1.1) THEN
C            WRITE(STDOUT,54) FR, Y3
            FR = 1.0
          ENDIF
C          WRITE(STDOUT,*) ' FCD123: FR=',FR,' A3=',A3
          FCD123 = 0.7269 + FR*(0.3333 - FR*0.1102)
        ELSE
          FCD123 = 0.95
        ENDIF
C       Apply adjustments for box culverts.
        FCD123 = KRB*KWING*KPROJ*FCD123
 
      ELSEIF(CULCLS.EQ.'MITER') THEN
        R = (Z1TRUE - ZB2)/DUP
        IF(R.LE.0.0.OR.R.GT.1.6) THEN
          IF(RATFLG.EQ.0) THEN
            RATFLG = 1
C            WRITE(STDOUT,50) R, Z1TRUE
          ENDIF
          R = 1.599
        ENDIF
        IF(R.LT.0.4) THEN
          FCD123 = 0.88
        ELSE
          FCD123 = 0.74542 + R*(0.50912 + R*(-0.46277 + R*0.077665))
        ENDIF
C       Apply adjustments for mitered pipe.
        FCD123 = KRB*KPROJ*FCD123
 
      ELSEIF(CULCLS.EQ.'RCPTG') THEN
        FCD123 = 0.95
      ELSE
        WRITE(STDOUT,56) CULCLS
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      IF(FCD123.GT.0.98) THEN
        FCD123 = 0.98
      ENDIF
C      WRITE(STDOUT, *) ' FCD123 =', FCD123
 
      RETURN
      END
C
C
C
      REAL FUNCTION   FDFRFC
     I                      (IU, ID)
 
C     + + + PURPOSE + + +
C     Find the factor on square of flow that will yield the
C     friction and eddy losses when the culvert barrel is flowing
C     full over its length.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IU
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'culcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL AL, AR, DX, KL, KR
      DOUBLE PRECISION SUM
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SNGL
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTA, LKTK
C***********************************************************************
      SUM = 0.D0
      CALL LKTA
     I         (NSEC(IU),
     M          DVEC(IU),
     O          AL)
      CALL LKTK
     I         (NSEC(IU),
     M          DVEC(IU),
     O          KL)
      DO 200 I=IU+1, ID
        CALL LKTA
     I           (NSEC(I),
     M            DVEC(I),
     O            AR)
        CALL LKTK
     I           (NSEC(I),
     M            DVEC(I),
     O            KR)
        DX = ABS(XVEC(I) - XVEC(I-1))
        IF(AL.LE.AR) THEN
C         The flow is expanding.
          SUM = SUM + DX/(KL*KR) + KD(I)*(1.0/AL**2 - 1.0/AR**2)/
     A                                                         GRAV2
        ELSE
C         The flow is contracting.
          SUM = SUM + DX/(KL*KR) + KA(I)*(1.0/AR**2 - 1.0/AL**2)/
     A                                                         GRAV2
        ENDIF
        KL = KR
        AL = AR
 200  CONTINUE
      FDFRFC = SNGL(SUM)
      RETURN
      END
C
C
C
      SUBROUTINE   FNDCD5
     I                   (STDOUT, CULCLS, HDRAT, RBV,
     O                    CDIS)
 
C     + + + PURPOSE + + +
C     Find the discharge coefficient for type 5 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      REAL CDIS, HDRAT, RBV
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     HDRAT  - ratio of head to vertical diameter
C     RBV    - rounding/beveling value
C     CDIS   - discharge coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'cdcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER NTAB
      REAL DF, DFDC, DFDR, WWCD
 
C     + + + INTRINSICS + + +
      INTRINSIC MAX
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTAB, TDLK10
C***********************************************************************
C     Get the basic coefficient.
      IF(CFHWA.GT.0.0) THEN
C       User has selected the FHWA definition of the type 5 discharge
C       coefficient.  Note that the FHWA coefficient c, CFHWA, has already
C       been multiplied by the value of 2g
        CDIS = SQRT((1.0 - (YFHWA - AFHWA*SFHWA)/HDRAT)/CFHWA)
      ELSE
        IF(CULCLS.EQ.'FLARED') THEN
C         Special 1-D table lookup
          CALL LKTAB
     I              (TB8ADR, HDRAT, 1,
     O               CDIS, NTAB, DF)
        ELSE
C         Get the basic value for all other types.
          CALL TDLK10
     I               (STDOUT, TB6ADR, 10, HDRAT, RBV,
     O                CDIS, DFDR, DFDC)
 
          IF(CULCLS.EQ.'MITER') THEN
            CDIS = 0.92*CDIS
          ELSEIF(CULCLS.EQ.'BOX') THEN
            IF(WWANGL.GT.0.0) THEN
C             Get the coefficient for a box culvert with wing walls.
 
              CALL TDLK10
     I                   (STDOUT, TB7ADR, 10, HDRAT, WWANGL,
     O                    WWCD, DFDR, DFDC)
 
              CDIS = MAX(CDIS, WWCD)
            ENDIF
          ENDIF
        ENDIF
 
        IF(CULCLS.NE.'BOX') THEN
          CDIS = CDIS*KPROJ
        ENDIF
      ENDIF
      IF(CDIS.GT.1.0) THEN
        CDIS = 0.98
      ENDIF
      CD = CDIS
      RETURN
      END
C
C
C
      SUBROUTINE   FSCMAT
     I                   (STDOUT, ADR, SBOT,
     O                    EFLAG, NSCMAT, YATMAT, SATMAT)
 
C     + + + PURPOSE + + +
C     Find the list of matches bewteen critical slope and bottom
C     slope for a cross section for a conduit. There may be none, one,
C     two, and maybe more matches.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADR, EFLAG, NSCMAT, STDOUT
      REAL SATMAT(7), SBOT, YATMAT(7)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     ADR    - Address of function table
C     SBOT   - bottom slope
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NSCMAT - number of critical slope matches
C     YATMAT - depth at match between critical slope and bottom slope
C     SATMAT - sign of the residual function slope when critical slope
C               and bottom slope match
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HIGH, I, KNT, KNT2, LOW, NSC, TYPE, XOFF
      REAL DIV, DYOVYH, K, NUM, QC, SC, SCHAT, SCL, SCMAX, SCMIN, SCR,
     A     SCVEC(MNDEP+3), TP, Y, YASCMN, YASCMX, YHAT, YL, YR,
     B     YVEC(MNDEP+3)
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, EXP, LOG
 
C     + + + EXTERNAL NAMES + + +
      CHARACTER GET_TABID*16
      EXTERNAL LKTK, LKTQC, GET_TABID
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* More than seven matches in FSCMAT.')
 52   FORMAT(/,' *ERR:688* TabId=',A,' has type=',I3,' but type=',
     A    ' 22 or 25 required.')
 54   FORMAT(/,' Smallest critical slope found=',F10.4,' at depth=',
     A    F10.3,/,' Largest critical slope found=',F10.4,' at ',
     B  'depth=',F10.3)
 56   FORMAT(/,' *WRN:534* No convergence on depth for So=Sc.  Last',
     A  ' relative correction=',1PE12.4,/,11X,'Continuing with latest',
     B  ' estimate of depth.')
 58   FORMAT(/,
     A' Depths for which critical slope, Sc, match, So, invert',
     B     ' slope=',F10.4)
 59   FORMAT('      Depth  (Sc-So)/So')
 60   FORMAT(1X,F10.3,1PE12.4)
C***********************************************************************
      TYPE = ITAB(ADR+2)
      IF(TYPE.NE.22.AND.TYPE.NE.25) THEN
        WRITE(STDOUT,52) GET_TABID(ITAB(ADR+1)), TYPE
        EFLAG = 1
        RETURN
      ENDIF
      XOFF = OFFVEC(TYPE)
 
      LOW = ADR + XTIOFF
C     Get the address of the highest argument of the table.
      HIGH = ITAB(ADR)
 
C     Try to force computation of the critical slope at
C     small depths.  Note that the definition of critical flow
C     in the table may be inaccurate at this point.
C     Compute critical slope for all but the first(zero depth) and
C     last(top of slot) depths for the cross section of a closed conduit.
 
C      YMIN = FTAB(LOW + XOFF)
 
      SCMAX = -1.E30
      SCMIN = 1.E30
C      DO 90 I=1,3
C        Y = YMIN*FLOAT(I)/4.0
C        CALL LKTK(ADR, Y, K)
C        CALL LKTQC(ADR, Y, QC)
C        SC = (QC/K)**2
C        SCVEC(I) = SC
C        YVEC(I) = Y
C        WRITE(STDOUT,*) ' Y=',Y,' SC=',SC
C        SCMAX = MAX(SC, SCMAX)
C        SCMIN = MIN(SC, SCMIN)
C90    CONTINUE
      NSC = 0
C     Compute the remaining points at the tabulated depths.
      DO 100 I=LOW+XOFF,HIGH-XOFF, XOFF
        Y = FTAB(I)
        K = FTAB(I+3)**2
        QC = FTAB(I+7)
        NSC = NSC + 1
        SC = (QC/K)**2
        IF(SC.GT.SCMAX) THEN
          SCMAX = SC
          YASCMX = Y
        ENDIF
        IF(SC.LT.SCMIN) THEN
          SCMIN = SC
          YASCMN = Y
        ENDIF
        SCVEC(NSC) = SC
        YVEC(NSC) = Y
C        WRITE(STDOUT,*) ' Y=',YVEC(NSC),' SC=',SCVEC(NSC)
 100  CONTINUE
 
 
      WRITE(STDOUT,54) SCMIN, YASCMN, SCMAX, YASCMX
 
      IF(SBOT.LT.SCMIN) THEN
C       There is no match.
        NSCMAT = 0
      ELSEIF(SBOT.GT.SCMAX) THEN
C       There is a match but we need not find it.   Note that critical
C       slope approaches infinity  as depth approaches zero and
C       as depth approaches D.   Very small depths are of no interest
C       and the match above the maximum in the table is not needed.
C       Signal that the slope will be steep for any critical flow
C       possible in the culvert.
        NSCMAT = -1
      ELSE
C       Do a search here and find all the matches with their
C       depths.  There should be at least one match.
        WRITE(STDOUT,58) SBOT
        WRITE(STDOUT,59)
        KNT = 0
        SCL = SCVEC(1)
        YL = YVEC(1)
        DO 200 I=2,NSC
          SCR = SCVEC(I)
          YR = YVEC(I)
          IF(SCR.LE.SBOT.AND.SCL.GE.SBOT.OR.
     A       SCL.LE.SBOT.AND.SCR.GE.SBOT) THEN
C           Interpolate for the value of depth using linear interpolation
C           in logarithms.
            KNT = KNT + 1
            IF(KNT.GT.7) THEN
              WRITE(STDOUT,50)
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            DIV = LOG(SCR/SCL)
            NUM = LOG(YR/YL)
            IF(DIV.EQ.0.0) THEN
C             The critical slope is constant.  Depth is undefined
C             in the interval.  Take the right hand value.
              YATMAT(KNT) = YR
              SATMAT(KNT) = 0.0
            ELSE
              YHAT = YL*EXP(LOG(SBOT/SCL)*NUM/DIV)
              IF(SCR.GT.SCL) THEN
                SATMAT(KNT) = 1.0
              ELSE
                SATMAT(KNT) = -1.0
              ENDIF
C             Refine the interpolation.  Errors can be large
C             enough to cause problems later.
              KNT2 = 0
 150          CONTINUE
                CALL LKTQC
     I                    (ADR,
     M                     YHAT,
     O                     QC)
                CALL LKTK
     I                   (ADR,
     M                    YHAT,
     O                    K)
                SCHAT = (QC/K)**2
                DYOVYH = -(1.0 - SBOT/SCHAT)*NUM/DIV
                IF(ABS(DYOVYH).GT.EPSF) THEN
                  TP = YHAT*(1.0 + DYOVYH)
                  IF(TP.LT.YL) THEN
                    YHAT = 0.5*(YHAT + YL)
                  ELSEIF(TP.GT.YR) THEN
                    YHAT = 0.5*(YHAT + YR)
                  ELSE
                    YHAT = TP
                  ENDIF
                  KNT2 = KNT2 + 1
                  IF(KNT2.GT.30) THEN
                    WRITE(STDOUT,56)  DYOVYH
                    GOTO 160
                  ENDIF
                  GOTO 150
                ENDIF
 160          CONTINUE
              WRITE(STDOUT,60) YHAT, (SCHAT - SBOT)/SBOT
              YATMAT(KNT) = YHAT
            ENDIF
          ENDIF
          SCL = SCR
          YL = YR
 200    CONTINUE
        IF(KNT.EQ.0) THEN
          WRITE(STDOUT,*) ' *BUG:XXX No match found in FSCMAT when',
     A           ' at least one should exist.'
          STOP 'Abnormal stop. Errors found.'
        ELSE
          NSCMAT = KNT
        ENDIF
      ENDIF
C      WRITE(STDOUT,*) ' NSCMAT=',NSCMAT
      RETURN
      END
C
C
C
      SUBROUTINE   MKXVEC
     I                   (LOVERD, D,
     M                    NHAT,
     O                    XVEC, IAT3D, IAT6D)
 
C     + + + PURPOSE + + +
C     Construct a pattern of points along the culvert barrel to
C     make sure that the profile computations are sufficiently
C     precise.  Also insert any special fixed locations.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER IAT3D, IAT6D, NHAT
      REAL D, LOVERD, XVEC(MNBN)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     LOVERD - ratio of culvert barrel length to vertical diameter
C     D      - vertical diameter of culvert barrel
C     NHAT   - number of nodes for culvert barrel
C     XVEC   - Stations for nodes on a branch or along a culvert
C     IAT3D  - index to vena contracta location
C     IAT6D  - index to vena contracta location
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IP
      DOUBLE PRECISION W(MNBN), X(MNBN)
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL GRULE
C***********************************************************************
      IF(NHAT.GT.MNBN) THEN
        NHAT = MNBN - 3
      ENDIF
 
C     Get the Gauss points.
 
      CALL GRULE
     I          (NHAT,
     O           X, W)
 
      XVEC(1) = 0.0
      IP = 2
      IAT3D = 0
      IAT6D = 0
      DO 100 I=1,NHAT
        XVEC(IP) = 0.5D0*(X(I) + 1.D0)*LOVERD
        IF(XVEC(IP).GT.3.1.AND.XVEC(IP-1).LT.2.9) THEN
C         Insert 3.0
          XVEC(IP+1) = XVEC(IP)
          XVEC(IP) = 3.0
          IAT3D = IP
          IP = IP + 1
        ELSEIF(XVEC(IP).GE.2.9.AND.XVEC(IP).LE.3.1) THEN
          IAT3D = IP
        ENDIF
        IF(XVEC(IP).GT.6.1.AND.XVEC(IP-1).LT.5.9) THEN
C         Insert 6.0
          XVEC(IP+1) = XVEC(IP)
          XVEC(IP) = 6.0
          IAT6D = IP
          IP = IP + 1
        ELSEIF(XVEC(IP).GE.5.9.AND.XVEC(IP).LE.6.1) THEN
          IAT6D = IP
        ENDIF
        IP = IP + 1
 100  CONTINUE
      NHAT = IP
      XVEC(NHAT) = LOVERD
 
      DO 200 I=2,NHAT
        XVEC(I) = D*XVEC(I)
 200  CONTINUE
 
 
      RETURN
      END
C
C
C
      REAL FUNCTION   R44TO4
     I                      (YA)
 
C     + + + PURPOSE + + +
C     Residual function for computing values at section 4 from
C     those at section 44.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL YA
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     YA     - unknown being sought
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs4com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'depmc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL Y
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
      Y = YA
      CALL XLKTAL
     I           (ADRXS4,
     M            Y,
     O            A4, T4, DT4, J4, K4, DK4, BET4, DBET4, ALP4, DALP4)
 
      R44TO4 = (Y + ALP4*(Q4/A4)**2/GRAV2) - E4 
 
      RETURN
      END
C
C
C
      REAL FUNCTION   R4TO44
     I                      (YA)
 
C     + + + PURPOSE + + +
C     Residual function for computing values at section 44 from
C     those at section 4.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL YA
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     YA     - unknown being sought
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'x44com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'depmc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL Y
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
      Y = YA
      CALL XLKTAL
     I           (ADRS44,
     M            Y,
     O            A44, T44, DT44, J44, K44, DK44, BET44, DBET44, ALP44,
     O            DALP44)
 
      R4TO44 = (Y + ALP44*(Q44/A44)**2/GRAV2) - E44 
 
      RETURN
      END
C
C
C
      REAL FUNCTION   RAPP
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute residual function for approach reach solution.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'rappc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL APPFAC, ARATIO, CD, EXLOSS, VH1, VH2
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
      Y1 = Y
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1,
     O            A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1)
 
      Z1 = ZB1  + Y1
      VH1 = (Q1/A1)**2/GRAV2
      VH2 = (Q2/A2)**2/GRAV2
 
C      IF(VH2.LT.VH1) THEN
      IF(A1.LE.A2) THEN
C       EXPANSION(NEGATIVE ACCELERATION) IN FLOW INSTEAD OF CONTRACTION
C        WRITE(STD6,*) ' APPRO: EXPANDING FLOW'
        ARATIO = A1/A2
        IF(ARATIO.GT.0.95) THEN
C         Interpolate coefficients to make the transition between
C         the two cases smooth.
          CD = CDIN + 20.0*(1.0 - ARATIO)*(0.98 - CDIN)
          APPFAC = 20.0*(1.0 - ARATIO)*APPEXP
        ELSE
          CD = 0.98
          APPFAC = APPEXP
        ENDIF
        CONF = 0
        EXLOSS = APPFAC*(ALP1*VH1 - ALP2*VH2)
        IF(EXLOSS.LT.0.0) THEN
          EXLOSS = 0.0
        ENDIF
        RAPP = (ALP1 - APPLOS)*VH1 + ZB1 + Y1 -
     A    (ALP2*VH2 + ZB2 + Y2  + APPLEN*Q1*Q2/(K1*K2)
     B     + EXLOSS + (1./CD**2 - 1.)*VHLOSS)
      ELSE
C        WRITE(STD6,*) ' APPRO: CONTRACTING FLOW'
C      WRITE(STD6,*) ' CDIN=',CDIN,' VHLOSS=',VHLOSS,' APPLEN=',APPLEN
C      WRITE(STD6,70)
C70    FORMAT(6X,7X,'Y',7X,'A',9X,'K',7X,'Q',6X,'VH',3X,'ALPHA',
C     A        5X,'TEL')
C      WRITE(STD6,50) Y1, A1, K1, Q1, VH1, ALP1,
C     A                        ZB1 + Y1 + ALP1*(Q1/A1)**2/GRAV2
C50    FORMAT(' SEC1:',F8.4,F8.1,F10.1,F8.2,F8.4,F8.4,F8.4)
C      WRITE(STD6,52) Y2, A2, K2, Q2, VH2, ALP2,
C     A                       ZB2 + Y2 + ALP2*(Q2/A2)**2/GRAV2
C52    FORMAT(' SEC2:',F8.4,F8.1,F10.1,F8.2,F8.4,F8.4,F8.4)
        CONF = 1
        CD = CDIN
        RAPP = ALP1*VH1 + ZB1 + Y1 -
     A    (ALP2*VH2 + ZB2 + Y2  + APPLEN*Q1*Q2/(K1*K2) + APPLOS*VH1
     B      + (1./CD**2 - 1.)*VHLOSS)
      ENDIF
 
C      WRITE(STD6,*) ' AT Y1=',Y1,' RAPP=',RAPP,' CONF=',CONF
      RETURN
      END
C
C
C
      REAL FUNCTION   RCON
     I                    (CONCOF)
 
C     + + + PURPOSE + + +
C     Residual function to find the contraction coefficient for
C     type 5 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL CONCOF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     CONCOF - contraction coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'rconc.cmn'
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL YOVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL YOVERD
C***********************************************************************
      RCON = (CONCOF/CDIS)**2*
     A   (1.0 - YOVERD(ADR, D, AFULL, CONCOF)/HOVERD + DZVC/H)  - 1.0
      RETURN
      END
C
C
C
      REAL FUNCTION   RCON2
     I                     (CONCOF)
 
C     + + + PURPOSE + + +
C     Residual function to find the contraction coefficient for
C     type 5 flow given the conditions at the limit of type 1 or
C     type 2 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL CONCOF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     CONCOF - contraction coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'rconc.cmn'
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL YOVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL YOVERD
C***********************************************************************
      RCON2 = (CONCOF**2*(EVC - D*YOVERD(ADR, D, AFULL, CONCOF)) - VHF)/
     A             VHF
      RETURN
      END
C
C
C
      REAL FUNCTION   RDPM26
     I                      (Y)
 
C     + + + PURPOSE + + +
C     Compute residual function for idealized departure reach
C     to define the section 44 values given the section 3 and
C     roadway flow values for flow types 2 and 6.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'x44com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'dpm26c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL M44
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
C     VARIABLES WITH SUFFIX 43 ARE DEFINED ON ENTRY AND GIVE
C     THE CROSS SECTION ELEMENTS FROM THE DEPARTURE SECTION AT
C     THE PIEZOMETRIC LEVEL AT THE EXIT OF THE CULVERT FOR FLOW
C     TYPES 2 AND 6. VALUES IN THE CULVERT AT THE EXIT OF THE
C     CULVERT ARE ALSO DEFINED ON ENTRY.
 
      Y44 = Y
      Z44 = ZB44 + Y44
 
      CALL XLKTAL
     I           (ADRS44,
     M            Y44,
     O            A44, T44, DT44, J44, K44, DK44, BET44, DBET44, ALP44,
     O            DALP44)
 
      M44 = GRAV*J44 + BET44*Q44**2/A44
      RDPM26 =  1.0 - M43/M44
 
C      WRITE(OUTUN,*) ' AT Y=',Y,' RDPM26 =', RDPM26
C      WRITE(OUTUN,*) ' Q44=',Q44,' Q3=',Q3,' J44=',J44,' J43=',J43,
C     A    ' A3=',A3,' A44=',A44,' RMFLUX=',RMFLUX,
C     B    ' BET44=',BET44,' BET3=',BET3,' M44=',M44,' M43=',M43,
C     C    ' ALP3=',ALP3,' ALP44=',ALP44,' J3=',J3,' J3Z43=',
C     D    J3Z43
 
      RETURN
      END
C
C
C
      REAL FUNCTION   RROVD
     I                     (RB)
 
C     + + + PURPOSE + + +
C     Residual function for finding the relative rounding/beveling
C     that causes the slope of the culvert to become the boundary
C     between types 5 and 6.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL RB
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     RB     - rounding/beveling value
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'rrovdc.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL RBV, TP
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FNDSZL
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDSZL
C***********************************************************************
      RBV = RB
      TP = FNDSZL(STD6, CLASS, LOVERD, NBAR,
     A                       RGHFAC, RBV, TB15AD, TB16AD)
 
      RROVD = 100.0*(TP - SZERO)
C      WRITE(STD6,*) ' RROVD:AT RB=',RB,' SZERO=',SZERO,' SLIM=',TP,
C     A              ' RROVD=',RROVD
      RETURN
      END
C
C
C
      SUBROUTINE   FNDCC
     I                  (STDOUT, CD, DUP, ADRS, AIN, HDRAT, HEAD, DZ,
     O                   CC)
 
C     + + + PURPOSE + + +
C     Estimate the contraction coefficient for type 5 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, STDOUT
      REAL AIN, CC, CD, DUP, DZ, HDRAT, HEAD
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CD     - discharge coefficient
C     DUP    - vertical diameter of culvert barrel at upstream end
C     ADRS   - Address of function table
C     AIN    - full area of conduit
C     HDRAT  - ratio of head to vertical diameter
C     HEAD   - head at section 1
C     DZ     - change in elevation from entrance to vena contracta
C     CC     - contraction coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'rconc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL CCHIGH, CCLOW, CCT, F, FH, FL
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL RCON
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL RCON, RGF
C***********************************************************************
C     Set values in common for finding the contraction coefficient.
      CDIS = CD
      HOVERD = HDRAT
      H = HEAD
      DZVC = DZ
      AFULL = AIN
      D = DUP
      ADR = ADRS
 
C     Make the first estimate.
      CC = 0.6
      CCLOW = 0.0
      CCHIGH = 0.0
 100  CONTINUE
        F = RCON(CC)
C        WRITE(STDOUT,*) ' FNDCC: CC=',CC,' F=',F
        IF(ABS(F).LE.EPSF) THEN
C         Close enough
        ELSE
          IF(F.GT.0.0) THEN
            CCHIGH = CC
            FH = F
            IF(CCLOW.EQ.0.0) THEN
              CC = 0.95*CC
              IF(CC.LE.EPSF) THEN
                WRITE(STDOUT,*) ' FNDCC: NO NEGATIVE RESIDUAL'
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 100
            ELSE
              GOTO 110
            ENDIF
          ELSEIF(F.LT.0.0) THEN
            CCLOW = CC
            FL = F
            IF(CCHIGH.EQ.0.0) THEN
              CCT = 1.05*CC
              IF(CCT.GT.1.0) THEN
                CCT = 0.5*(1.0 + CC)
              ENDIF
              CC = CCT
              GOTO 100
            ENDIF
          ENDIF
 110      CONTINUE
 
          CALL RGF
     I            (1.E-6, EPSF, RCON,
     M             CCLOW, CCHIGH, FL, FH,
     O             CC, FLAG)
 
        ENDIF
 
C      WRITE(STDOUT,*) ' CC=',CC,' FL=',FL
 
      RETURN
      END
C
C
C
      SUBROUTINE   FNDCC2
     I                   (STDOUT, DUP, ZBVC,
     O                    CC)
 
C     + + + PURPOSE + + +
C     Estimate the contraction coefficient for type 5 flow at the
C     limit of type 1 or type 2 flow.  Results from the limit
C     computations are assumed to be in common block TYPLIM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      REAL CC, DUP, ZBVC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     DUP    - vertical diameter of culvert barrel at upstream end
C     ZBVC   - bottom elevation at the vena contracta
C     CC     - contraction coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'rconc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, KNT
      REAL CCHIGH, CCLOW, CCT, EPSFL, F, FH, FL
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL RCON2
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL RCON2, RGF
C***********************************************************************
      EPSFL = 0.5*EPSF
C     Set values in common for finding the contraction coefficient.
      AFULL = A2FULL
      D = DUP
      ADR = ADRXS2
C     Compute the specific energy at the vena contracta and the
C     full flow velocity head at section 2.
C      WRITE(STDOUT,*) ' FNDCC2: Z1L=',Z1L,' Q1L=',Q1L,' A1L=',A1L
C      WRITE(STDOUT,*) ' A2FULL=',A2FULL,' Q2L=',Q2L,' ZBVC=',ZBVC
C      WRITE(STDOUT,*) ' K1L=',K1L,' K2FULL=',K2FULL
      VHF = (Q2L/A2FULL)**2/GRAV2
      EVC = (ALP1L - APPLOS)*(Q1L/A1L)**2/GRAV2 + Z1L - ZBVC
     A       - APPLEN*(Q1L*Q2L)/(K1L*K2FULL)
 
C     Make the first estimate.
      CC = 0.6
      CCLOW = 0.0
      CCHIGH = 0.0
      KNT = 0
 100  CONTINUE
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          WRITE(STDOUT,*) ' *BUG:XXX* No convergence in FNDCC2'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        F = RCON2(CC)
C        WRITE(STDOUT,*) ' FNDCC2: CC=',CC,' F=',F
        IF(ABS(F).LE.EPSFL) THEN
C         Close enough
        ELSE
          IF(F.GT.0.0) THEN
            CCHIGH = CC
            FH = F
            IF(CCLOW.EQ.0.0) THEN
              CC = 0.95*CC
              IF(CC.LE.EPSF) THEN
                WRITE(STDOUT,*) ' FNDCC2: NO NEGATIVE RESIDUAL.'
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 100
            ELSE
              GOTO 110
            ENDIF
          ELSEIF(F.LT.0.0) THEN
            CCLOW = CC
            FL = F
            IF(CCHIGH.EQ.0.0) THEN
              CCT = 1.05*CC
              IF(CCT.GT.1.0) THEN
                CCT = 0.5*(1.0 + CC)
                IF(CCT.GE.1.0 - EPSF) THEN
                  WRITE(STDOUT,*) ' FNDCC2: NO POSITIVE RESIDUAL.'
                  STOP 'Abnormal stop. Errors found.'
                ENDIF
              ENDIF
              CC = CCT
              GOTO 100
            ENDIF
          ENDIF
 110      CONTINUE
 
          CALL RGF
     I            (1.E-6, EPSFL, RCON2,
     M             CCLOW, CCHIGH, FL, FH,
     O             CC, FLAG)
 
        ENDIF
 
C      WRITE(STDOUT,*) ' CC=',CC,' FL=',FL
 
      RETURN
      END
C
C
C
      SUBROUTINE   FROVD
     I                  (STDOUT, CULCLS,
     M                   ROVD)
 
C     + + + PURPOSE + + +
C     Find the value of relative rounding/beveling that causes
C     the type 5/6 limiting slope to match the slope of the culvert.
C     Used to assist in estimating the full-flow-inducing condition
C     in the culvert for type 5 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      REAL ROVD
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     ROVD   - relative rounding/beveling found
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rrovdc.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL F, FL, FR, RB, RBL, RBMAX, RBMIN, RBR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL RROVD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL RGF, RROVD
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF CLAIMS NONE IN',
     A       ' FROVD.')
 60   FORMAT(' *BUG:XXX* RGF: MORE THAN 100 ITERATIONS IN FROVD.')
C***********************************************************************
      CLASS = CULCLS
 
C     Since the flow is known to be type 5 it must be true that
C     the value of relative rounding/beveling that causes a
C     match with the culvert slope is greater than the value
C     that currently exists for the culvert.  That is the culvert
C     slope is already known to be larger than the limit slope
C     for type 6 flow.
 
      RBMIN = RBVAL
 
      IF(CULCLS.EQ.'BOX') THEN
        RBMAX = 0.06
      ELSE
        IF(NBAR.LT.0.019) THEN
C         Barrel is smooth.
          RBMAX = 0.06
        ELSE
          RBMAX = 0.03
        ENDIF
      ENDIF
 
c     25 August 2004:  Supply an initial value.

      ROVD= 0.5*(RBMIN + RBMAX)
      RB = ROVD
 
      FL = 0.0
      FR = 0.0
      F = RROVD(RB)
      IF(ABS(F).LE.EPSF) THEN
C       DONE
      ELSE
        IF(F.LT.0.0) THEN
          RBL = RB
          FL = F
          IF(FR.EQ.0.0) THEN
            F = RROVD(RBMAX)
            IF(F.LT.0.0) THEN
C             Maximum value reached and no match.  Take maximum value.
              ROVD = RBMAX
              RETURN
            ELSE
              FR = F
              RBR = RBMAX
            ENDIF
          ENDIF
        ELSEIF(F.GT.0.0) THEN
          RBR = RB
          FR = F
          IF(FL.EQ.0.0) THEN
            F = RROVD(RBMIN)
            IF(F.GT.0.0) THEN
C             Minimum value reached.  Should not happen.
        WRITE(STDOUT,*) ' *BUG:XXX Min. RB value problem in FROVD'
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            FL = F
            RBL = RBMIN
          ENDIF
        ENDIF
C        WRITE(STDOUT,*) ' FROVD: RBL=',RBL,' RBR=',RBR
C        WRITE(STDOUT,*) ' FROVD: FL=',FL,' FR=',FR
        CALL RGF
     I          (EPSARG, EPSF, RROVD,
     M           RBL, RBR, FL, FR,
     O           ROVD, FLAG)
 
        IF(FLAG.EQ.1) THEN
          WRITE(STDOUT, 54)
          STOP 'Abnormal stop. Errors found.'
        ELSEIF(FLAG.EQ.2) THEN
          WRITE(STDOUT,60)
          STOP 'Abnormal stop. Errors found.'
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   APPRO
     I                  (STDOUT, CD, VHL, QROAD,
     O                   CONFLG, NSFLAG)
 
C     + + + PURPOSE + + +
C     Compute the elevation at section 1 given the flow and the elevation
C     at section 2.  We assume that all values at section 2 are known.
C     CD gives the discharge coefficient for entrance losses and
C     VHL is the velocity head to which it is applied.  Thus APPRO
C     can be used for all types of flow.  CD is based on the known
C     target area at section 1.  APPRO will be used to compute the
C     residual at section 1 between the target elevation and the
C     estimated elevation.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER CONFLG, NSFLAG, STDOUT
      REAL CD, QROAD, VHL
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CD     - Discharge coefficient
C     VHL    - Provides the velocity head to use in computed the
C               entrance losses
C     QROAD  - Flow over the roadway
C     CONFLG - CONFLG=0: flow contracts as it enters the culvert and
C              CONFLG=1: flow expands as it enters the culver
C     NSFLAG - No solution flag
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'rappc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, MINFLG
      REAL DY, F, FH, FL, FOLD, Y, YC, YHIGH, YLOW, YMAX, YOLD, YSTART
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL FMXARG, RAPP
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FMXARG, FNDCDE, GETTBN, RAPP, REGFLT
 
C     + + + OUTPUT FORMATS + + +
 53   FORMAT(' *WRN:537* TABID=',A,' overflow seeking XS1 depth for',
     a       ' flow=',F10.1,/,10X,' in subroutine APPRO.')
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' APPRO.')
 60   FORMAT(' *BUG:XXX* REGFLT: MORE THAN 100 ITERATIONS IN APPRO.')
C***********************************************************************
C     SET VALUES FOR LOCAL COMMON BLOCK
      VHLOSS = VHL
      CDIN = CD
 
      Q1 = Q2 + QROAD
 
C     CLEAR NO SOLUTION FLAG INDICATING THAT THERE IS A SOLUTION
      NSFLAG = 0
 
      YMAX = FMXARG(ADRXS1)
      MINFLG = 0
C     CHECK FOR A SOLUTION
C     Try the value at the known depth
      F = RAPP(YUPTRU)
C      WRITE(STDOUT,*) ' APPRO: YUPTRU=',YUPTRU,' F=',F,' Q2=',Q2,
C     A   ' Y2=',Y2,' CONF=',CONF
      IF(ABS(F).LE.EPSABS) THEN
C       Close enough.
C        WRITE(STDOUT,*) ' APPRO: FIRST RETURN. Y1=', Y1
        CONFLG = CONF
        RETURN
      ELSE
        IF(F.LT.0.0) THEN
C         Start search for sign change from YUPTRU
          YSTART = YUPTRU
        ELSE
C         Evaluate at the critical depth if critical depth
C         differs from YUPTRU.
          YC = YUPTRU
          CALL FNDCDE
     I               (STDOUT, ADRXS1, Q1,
     M                YC)
          IF(ABS(YUPTRU-YC).LT.EPSARG) THEN
C           Start search at YUPTRU and stop at minimum.
            MINFLG = 1
            YSTART = YUPTRU
          ELSE
            F = RAPP(YC)
C            WRITE(STDOUT,*) ' APPRO: YC=',YC,' F=',F, ' CONF=',CONF
            IF(F.LT.0.0) THEN
C             Start search at YC
              YSTART = YC
            ELSE
C             Start search at YUPTRU and stop at minimum
              MINFLG = 1
              YSTART = YUPTRU
            ENDIF
          ENDIF
        ENDIF
      ENDIF
 
      IF(MINFLG.EQ.0) THEN
C       At least one subcritical solution exists.  Search upward for it.
        YLOW = 0.0
        YHIGH = 0.0
        Y = YSTART
        DY =  YUPTRU/8.0
C       SEARCH FOR SIGN CHANGE
 100    CONTINUE
          IF(F.LT.0.0) THEN
            FL = F
            YLOW = Y
          ELSE
            FH = F
            YHIGH = Y
          ENDIF
          IF(YLOW.EQ.0.0.OR.YHIGH.EQ.0.0) THEN
            Y = Y + DY
            DY = 1.25*DY
            IF(Y.GT.YMAX) THEN
              TABID = GET_TABID(GETTBN(ADRXS1))
              WRITE(STDOUT,53) TABID(1:LENSTR(TABID)), Q2
              WRITE(STDOUT,*) ' Maximum depth=',YMAX
              NSFLAG = 2
              CONFLG = CONF
              Z1 = ZB1
              RETURN
            ELSE
              F = RAPP(Y)
C              WRITE(STDOUT,*) ' MINFLG=',MINFLG,' Y=',Y,' F=',F
              GOTO 100
            ENDIF
          ENDIF
C         Change of sign found here.
      ELSE
C       There may be no solution.  Search from YSTART and stop
C       when there is an increase in the residual.  Then if the
C       residual is negative, continue search for a sign change.
C       If not negative, there is no subcritical solution.
        FOLD = F
        YOLD = YSTART
C        WRITE(STDOUT,*) ' YOLD=',YOLD,' FOLD=',FOLD
        Y = YSTART
        DY = 0.0125
 200    CONTINUE
          Y = Y + DY
          F = RAPP(Y)
C          WRITE(STDOUT,*) ' MIN. SEARCH. Y=',Y,' F=',F
          IF(F.GT.FOLD) THEN
C           Stop the search.  FOLD is the minimum value found.
            IF(ABS(FOLD).LE.EPSABS) THEN
C             Close enough.  Reset the values.
              F = RAPP(YOLD)
              CONFLG = CONF
C              WRITE(STDOUT,*) ' APPRO: SECOND RETURN. Y1=',Y1
              RETURN
            ELSE
              IF(FOLD.LT.0.0) THEN
                YLOW = 0.0
                YHIGH = 0.0
                F = FOLD
                Y = YOLD
 210            CONTINUE
                  IF(F.LT.0.0) THEN
                    YLOW = Y
                    FL = F
                  ELSE
                    YHIGH = Y
                    FH = F
                  ENDIF
                  IF(YLOW.EQ.0.0.OR.YHIGH.EQ.0.0) THEN
                    Y = Y + DY
                    IF(Y.GT.YMAX) THEN
                      WRITE(STDOUT,53) GETTBN(ADRXS1), Q2
                      NSFLAG = 2
                      CONFLG = CONF
                      Z1 = ZB1
                      RETURN
                    ELSE
                      F = RAPP(Y)
C              WRITE(STDOUT,*) ' MINFLG=',MINFLG,' Y=',Y,' F=',F
                      GOTO 210
                    ENDIF
                  ENDIF
              ELSE
C               No subcritical solution.
                Z1 = ZB1 + 0.5*YSTART
C                WRITE(STDOUT,*) ' APPRO: No solution. Y1=',
C     A                      0.5*YSTART
                NSFLAG = 1
                CONFLG = CONF
                RETURN
              ENDIF
            ENDIF
          ELSE
            FOLD = F
            YOLD = Y
            GOTO 200
          ENDIF
      ENDIF
 
C      WRITE(STDOUT,*) ' SIGN CHANGE HERE.'
C      WRITE(STDOUT,*) ' YLOW=',YLOW,' FL=',FL
C      WRITE(STDOUT,*) ' YHIGH=',YHIGH,' FH=',FH
 
      CALL REGFLT
     I           (EPSARG, EPSABS, RAPP,
     M            YLOW, YHIGH, FL, FH,
     O            Y1, FLAG)
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
C      WRITE(STDOUT,*) ' APPRO: Y1=',Y1
 
C     FINAL RESULTS IN XS1COM
C      WRITE(STDOUT,*) ' APPRO: CONF=', CONF
      CONFLG = CONF
      RETURN
      END
C
C
C
      SUBROUTINE   F44TO4
     I                   (STDOUT,
     O                    NSFLAG)
 
C     + + + PURPOSE + + +
C     Compute the values at section 4 given values at section 44.
C     All values for section 44 in its common block are known.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NSFLAG, STDOUT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NSFLAG - No solution flag
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs4com.cmn'
      INCLUDE 'x44com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'depmc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL FL, FR, VH44, YC, YL, YR
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL R44TO4
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDCDE, R44TO4, REGFLT
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' F44TO4.')
 60   FORMAT(' *BUG:XXX* REGFLT: MORE THAN 100 ITERATIONS IN F44TO4.')
C***********************************************************************
      NSFLAG = 0
C     Flow at the two sections is always the same.  Road flow
C     contribution is assumed to take place between sections 43 and
C     sections 44 where we use a simple momentum balance.
 
      Q4 = Q44
 
      IF(ADRXS4.EQ.ADRS44.AND.ZB4.EQ.ZB44) THEN
C       The departure reach is horizontal and prismatic.  Transfer
C       values from section 44 to section 4.
        Z4 = Z44
        Y4 = Y44
        A4 = A44
        T4 = T44
        J4 = J44
        BET4 = BET44
        ALP4 = ALP44
        K4 = K44
      ELSE
C       Departure reach is either non-prismatic or non-horizontal.
C       Compute values at section 4.   Find critical depth at
C       section 4.  Estimate the depth at section 4.
 
        YC = Y44 + ZB44 - ZB4
        IF(YC.LE.0.0) THEN
          YC = 0.1
        ENDIF
        CALL FNDCDE
     I             (STDOUT, ADRXS4, Q4,
     M              YC)
 
C       Critical flow at section 4 gives a lower bound for the
C       flow depth at section 4.  Estimate the specific energy
C       available at section 4.
 
 
        VH44 = ALP44*(Q44/A44)**2/GRAV2
        ZTEL44 = ZB44 + Y44 + VH44
 
        E4 = ZTEL44 - ZB4
 
C       Upper bound for depth at section 4 is the specific energy.
 
        Y4 = E4
 
        IF(Y4.LE.YC) THEN
C         No  solution exists.
          NSFLAG = 1
          Y4 = 0.0
        ELSE
C         A subcritical solution exists.  Try to find it.
 
          YL = YC
          FL = R44TO4(YL)
          YR = Y4
          FR = R44TO4(YR)
 
          CALL REGFLT
     I               (EPSARG, EPSABS, R44TO4,
     M                YL, YR, FL, FR,
     O                Y4, FLAG)
 
          Z4 = Y4 + ZB4
 
C          WRITE(STDOUT,*) ' F44TO4: ZB4=',ZB4,' Z4=',Z4
          IF(FLAG.EQ.1) THEN
            WRITE(STDOUT, 54)
            STOP 'Abnormal stop. Errors found.'
          ELSEIF(FLAG.EQ.2) THEN
            WRITE(STDOUT,60)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
        ENDIF
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   F4TO44
     I                   (STDOUT,
     O                    NSFLAG)
 
C     + + + PURPOSE + + +
C     Compute the values at section 44 given values at section 4.
C     All values for section 4 in its common block are known.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NSFLAG, STDOUT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NSFLAG - No solution flag
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs4com.cmn'
      INCLUDE 'x44com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'depmc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG
      REAL E44C, FL, FR, VH4, YL, YR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL R4TO44
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDCDE, R4TO44, REGFLT, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' F4TO44.')
 60   FORMAT(' *BUG:XXX* REGFLT: MORE THAN 100 ITERATIONS IN F4TO44.')
C***********************************************************************
      NSFLAG = 0
C     Flow at the two sections is always the same.  Road flow
C     contribution is assumed to take place between sections 43 and
C     sections 44 where we use a simple momentum balance.
 
      Q44 = Q4
 
      IF(ADRXS4.EQ.ADRS44.AND.ZB4.EQ.ZB44) THEN
C       The departure reach is horizontal and prismatic.  Transfer
C       values from section 4 to section 44.
        Z44 = Z4
        Y44 = Y4
        A44 = A4
        T44 = T4
        J44 = J4
        BET44 = BET4
        ALP44 = ALP4
        K44 = K4
      ELSE
C       Departure reach is either non-prismatic or non-horizontal.
C       Compute values at section 44.   Find critical depth at
C       section 44.  Estimate starting depth.
 
        Y44C = Y4 + ZB4 - ZB44
        IF(Y44C.LE.0.0) THEN
          Y44C = 0.1
        ENDIF
        CALL FNDCDE
     I             (STDOUT, ADRS44, Q44,
     M              Y44C)
 
C        WRITE(STDOUT,*) ' F4TO44: YC AT 44=',Y44C
        Y44 = Y44C
        CALL XLKTAL
     I             (ADRS44,
     M              Y44,
     O              A44, T44, DT44, J44, K44, DK44, BET44, DBET44,
     O              ALP44, DALP44)
 
C       Find the minimum specific energy at section 44, E44C
 
        E44C = Y44C + ALP44*(Q44/A44)**2/GRAV2
 
C        WRITE(STDOUT,*) ' E44C=',E44C
 
C       Find the specific energy available at section 44, E44, that
C       is, the energy as transfered from section 4.
 
        VH4 = ALP4*(Q4/A4)**2/GRAV2
        ZTEL4 = ZB4 + Y4 + VH4
 
C        WRITE(STDOUT,*) ' ZTEL4=',ZTEL4,' VH4=',VH4,' ZB4=',ZB4
 
        E44 = ZTEL4 - ZB44
 
C        WRITE(STDOUT,*) ' F4TO44: E44=',E44,' ZB44=',ZB44
 
        IF(E44.LT.E44C) THEN
C         No  solution exists.
          NSFLAG = 1
          Y44 = 0.0
        ELSE
C         A solution exists.
          YL = Y44C
          FL = R4TO44(YL)
          IF(ABS(FL).LE.EPSABS) THEN
C           Take the solution to be at critical depth at section 44.
C           Values already exist in the common block for section 44.
          ELSE
C           A subcritical solution exists.  Try to find it.  Upper
C           bound for depth is the specific energy.
 
            YR = E44
            FR = R4TO44(YR)
 
C            WRITE(STDOUT,*) ' Before call to REGFLT: YL=',YL,
C     A     ' FL=',FL,' YR=',YR,' FR=',FR
 
            CALL REGFLT
     I                 (EPSARG, EPSABS, R4TO44,
     M                  YL, YR, FL, FR,
     O                  Y44, FLAG)
 
            Z44 = Y44 + ZB44
 
            IF(FLAG.EQ.1) THEN
              WRITE(STDOUT, 54)
              STOP 'Abnormal stop. Errors found.'
            ELSEIF(FLAG.EQ.2) THEN
              WRITE(STDOUT,60)
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
        ENDIF
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   DPM26
     I                  (STDOUT, QRF, MRF, Z3T,
     O                   EXPFLG)
 
C     + + + PURPOSE + + +
C     Find the section 4 elevation given values at section 3 for
C     flow types 0, 1, 2, 5, and 6 using idealized momentum balance.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EXPFLG, STDOUT
      REAL MRF, QRF, Z3T
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     QRF    - Flow over the roadway
C     MRF    - Momentum flux from flow over the road
C     Z3T    - Elevation of water surface at section 3
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'x44com.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'dpm26c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, NSFLAG
      REAL FH, FL, TP, YHIGH, YL, YLOW, YMAX, YT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FMXARG, RDPM26
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F44TO4, F4TO44, FMXARG, FNDCDE, LKTJ, RDPM26, RGF,
     A         XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF CLAIMS NONE IN',
     A       ' DPM26.')
 56   FORMAT(' *BUG:XXX* NO POSITIVE RESIDUAL FOR DPM26.')
 60   FORMAT(' *BUG:XXX* RGF: MORE THAN 100 ITERATIONS IN DPM26.')
 64   FORMAT(/,' *WRN:536* No expansion of flow in departure reach.',
     A  '  Free flow becomes type 7.')
 68   FORMAT(' *BUG:XXX* No solution in F44TO4 in DPM26.')
 70   FORMAT(' Negative residual impossible in DPM26. Try increasing',
     A       ' momentum flux over roadway.')
C***********************************************************************
C     SET EXPANSION FLAG TO TRUE
 
      EXPFLG = 1
 
C     SET VALUES IN SPECIAL COMMON
      OUTUN = STDOUT
      RMFLUX = MRF
 
C     FIND VALUES IN DEPARTURE REACH AT UPSTREAM END AT THE CRITICAL
C     ELEVATION SUPPLIED BY Z3T.
 
      Z43 = Z3T
      Y43 = Z3T - ZB43
 
      CALL XLKTAL
     I           (ADRS43,
     M            Y43,
     O            A43, T43, DT43, J43, K43, DK43, BET43, DBET43, ALP43,
     O            DALP43)
 
C     Find the first moment of area in the barrel at the
C     water level at section 43.  This is often the same as
C     the level at section 3 but not for type 5 flow.
 
      TP = Z43-ZB3
      CALL LKTJ
     I         (ADRXS3,
     M          TP,
     O          J3Z43)
C      WRITE(STDOUT,*) ' DPM26: J3Z43=',J3Z43,' J3=',J3,' J43=',J43
C      WRITE(STDOUT,*) ' Y3=',Y3,' A3=',A3,' Q3=',Q3,' RMFLUX=',
C     A                 RMFLUX
C      WRITE(STDOUT,*) ' Z43 - ZB3=',Z43 - ZB3
C     Compute the momentum flux + impulse function for section
C     43.  Flux always comes from conditions in the barrel and from
C     the flow over the roadway(if there is any).  If the piezometric
C     levels at section 43 and section 3 are the same, then J3 and
C     J3Z43 have the same value.
 
C     There are cases in which the exit area has been falsified in
C     order to compute a transition between free flow types.
C     If this has been done, then BETAF in TYPTRN.COM will be
C     > 0.0.
      IF(BETAF.GT.0.0) THEN
        BETA3 = BETAF
        ALPHA3 = ALPHAF
      ELSE
        BETA3 = BET3
        ALPHA3 = ALP3
      ENDIF
      M43 = GRAV*(J43 + J3 - J3Z43) + BETA3*Q3**2/A3 + RMFLUX
 
C     Remember key values in order to solve for type 7 if that becomes
C     necessary.
      Z43OLD = Z43
      Q43OLD = Q3
C     COMPUTE THE FLOW AT SECTION 44
 
      Q44 = Q3 + QRF
      Q4 = Q44
 
C     WILL CRITICAL DEPTH IN SECTION 4 DROWN THE CONTROL AT SECTION 3?
C     FIND CRITICAL DEPTH IN SECTION 4 - GIVE ESTIMATED DEPTH TO
C     START PROCESS
 
      Y4C = Z3T - ZB4
      IF(Y4C.LT.0.0) THEN
        Y4C = 0.1
      ENDIF
 
      CALL FNDCDE
     I           (STDOUT, ADRXS4, Q4,
     M            Y4C)
 
      Z4 = ZB4 + Y4C
 
C      WRITE(STDOUT,*) ' CRIT DEPTH AT 4 IN DPM26=',Y4C,' AT FLOW=',Q4
 
C     Now compute the values at section 44 given that section 4 is
C     critical.  Need to lookup values at section 4 because not all of
C     them are defined by FNDCDE.
 
      Y4 = Y4C
      CALL XLKTAL
     I           (ADRXS4,
     M            Y4C,
     O            A4, T4, DT4, J4, K4, DK4, BET4, DBET4, ALP4, DALP4)
 
C      WRITE(STDOUT,*) ' Y4=',Y4,' ZB4=',ZB4,' ALP4=',ALP4,
C     A                 ' A4=',A4
      CALL F4TO44
     I           (STDOUT,
     O            NSFLAG)
 
      IF(NSFLAG.EQ.1) THEN
C       No solution at section 44.  Therefore, critical depth
C       at section 4 will not drown the flow at the exit of the
C       culvert.  Search for a negative residual for a depth at
C       section 44.  F4TO44 has computed critical depth at section 44.
C       Try for negative residual at that depth.
 
        Y44 = Y44C
        YL = Y44
        FL = RDPM26(YL)
        IF(FL.GT.0.0) THEN
          WRITE(STDOUT,70)
          WRITE(STDOUT,*) ' No solution in F4TO44 in DPM26:'
          WRITE(STDOUT,*)  ' YL=',YL,' FL=',FL
          STOP 'Abnormal stop. Errors found.'
        ENDIF
      ELSE
C       Solution at section 44.  Now see if the flow at culvert
C       exit(section 3) is drowned.
 
        FL = RDPM26(Y44)
        YLOW = Y44
         IF(FL.GT.0.0) THEN
C         CRITICAL DEPTH AT THE END OF THE DEPARTURE REACH WOULD
C         DROWN THE CONTROL.  WE ASSUME THAT THERE IS THEREFORE
C         NO EXPANSION OF THE FLOW FROM THE CULVERT TO THE DEPARTURE
C         REACH.
 
          WRITE(STDOUT,64)
          EXPFLG = 0
C          WRITE(STDOUT,*) ' Solution in F4TO44 in DPM26: Y44=',Y44
C          WRITE(STDOUT,*) ' FL at Y44=',FL
          RETURN
        ENDIF
      ENDIF
 
C     CRITICAL DEPTH WILL NOT DROWN THE CONTROL BUT SOME HIGHER ELEVATION
C     WILL.
 
 
C     SEARCH FOR SIGN CHANGE IN THE RESIDUAL FUNCTION
C     ALREADY IS <= 0.0 AT Y44
 
      YMAX = FMXARG(ADRS44)
 
C     With expansion in the departure reach the water surface elevation
C     at section 44 should be at least as high as at section 43.
 
      YHIGH = 1.1*Y44
      IF(YHIGH + ZB44.LT. Z43) THEN
        YHIGH = Z43 - ZB44
      ENDIF
 110  CONTINUE
        FH = RDPM26(YHIGH)
        IF(FH.LT.0.0) THEN
          FL = FH
          YLOW = YHIGH
          YHIGH = 1.2*YHIGH
          IF(YHIGH.GT.YMAX) THEN
            YHIGH = .5*(YLOW + YMAX)
            IF(ABS(YHIGH - YMAX).LT.EPSABS) THEN
              WRITE(STDOUT,56)
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
          GOTO 110
        ENDIF
 
C     WE HAVE A SIGN CHANGE HERE
 
C      XMIN = YLOW
C      XMAX = YHIGH
C      CALL SPSCNT(YLOW, YHIGH, EPSARG, EPSF, 100, FLAG, XMIN,
C     A                  XMAX, DPM26)
C      YT = YHIGH
      CALL RGF
     I        (EPSARG, EPSF, RDPM26,
     M         YLOW, YHIGH, FL, FH,
     O         YT, FLAG)
 
      Y44 = YT
      Z44 = Y44 + ZB44
C      WRITE(STDOUT,*) ' Solution for Y44 at free flow limit in DPM26:',
C     A               ' Y44=',Y44,' A44=',A44
 
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C     Now find the values in section 4 given the values in section 44
 
      CALL F44TO4
     I           (STDOUT,
     O            NSFLAG)
      IF(NSFLAG.EQ.1) THEN
        WRITE(STDOUT,68)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C      WRITE(STDOUT,*) ' DPM26: Solution at 4: Y4=',Y4,' Z4=',Z4
      RETURN
      END
C
C
C
      SUBROUTINE   LOCJMP
     I                   (STDOUT, IU, ID, YUP, YDN, Q, DUP,
     M                    TY6LSS,
     O                    YHIGH, JMPLOC, YEND, IERR, PROTYP)
 
C     + + + PURPOSE + + +
C     Locate a jump, if one exists, in the culvert barrel for the
C     given flow and initial depths at the starting node and ending
C     node.  Return details about the jump, the exit depth of the
C     culvert, and so forth.  JMPLOC= 0 if jump does not exist
C     and index of the node downstream of its location otherwise.
C     YHIGH gives the estimated depth on the high side of the jump
C     and is 0.0 otherwise.  YEND gives the end depth at the culvert
C     exit in all cases for which an end depth exists. IERR is 0
C     if the returned values have meaning and > 0 otherwise.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IERR, IU, JMPLOC, STDOUT
      REAL DUP, Q, TY6LSS, YDN, YEND, YHIGH, YUP
      CHARACTER PROTYP*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     YUP    - depth at vena contracta
C     YDN    - depth at exit of culvert for defining subcritical profile
C     Q      - Flowrate
C     DUP    - vertical diameter of culvert barrel at upstream end
C     TY6LSS - estimated type 6 loss
C     YHIGH  - depth on high side of hydraulic jump
C     JMPLOC - index to node below the jump
C     YEND   - depth found to exist at the exit of the culvert
C     IERR   - error flag
C     PROTYP - descripter for the nature of the profile
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'culcom.cmn'
      INCLUDE 'grvcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IC, ICSB, ICSP, IJUMP, ISB, ISP, ISUB, ISUP, SBFLAG,
     A        SPFLAG
      REAL ALPT, AT, BETT, DALPT, DBETT, DH, DKT, DTT, JT, KT, MSUBL,
     A     MSUBR, MSUPL, MSUPR, P, QCT, SBZDN, SBZUP, SPZDN, SPZUP, TT,
     B     YSUB, YSUBL, YSUBR, YSUPL, YSUPR, YT
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DISLSS, LKTJ, SFPSBE, SFPSPE, XLKT22
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* Invalid case=',I5,' in subroutine LOCJMP.')
 56   FORMAT(/,' *WRN:593* The super critical profile for Type 5 flow',
     1    ' is too short.',/,' This suggests that the',
     B   ' expansion loss applied after the vena contracta',/,' is',
     C   ' too large.  The loss will be reduced one or more times to',
     D   ' find',/,' a longer super critical profile.  The value of KD',
     E    ', the',/,' expansion loss in the barrel may be too large.',
     F    '  Values > 0.4 should be avoided.')
 57   FORMAT(/,' Initial expansion loss=',F8.3)
 58   FORMAT(/,' *WRN:594* No super critical profile of adequate ',
     A     'length exists.',
     B    /,' Culvert barrel may be non-prismatic or have significant',
     C    /,' changes in slope.')
 59   FORMAT(/,' Final expansion loss=',F8.3)
C***********************************************************************
C     Approach:  Attempt to compute a supercritical profile from the
C     starting node with the given depth.  Then attempt to compute
C     a subcritical profile from the ending node with the given depth.
C     The attempt to compute each profile has three
C     outcomes:  1. The profile could not be started, that is , not
C     even one distance step was possible.  This is failure.
C     2. A partial profile was computed.  3. A complete profile was
C     computed.  Thus there are nine outcomes possible.
 
C      Super     Sub-     #          Significance
C      critical  critical
C      profile   profile
C      -------   -------   ---------------------------------------------
C      fail      fail     1 Should not happen.  May indicate a strange
C                           culvert barrel slope variation.  Reduce the
C                           estimated loss, TY6LSS, and try the
C                           supercritical profile until the loss is
C                           essentially zero.
C      fail      partial  2 Should not happen.  May indicate a strange
C                           culvert barrel slope variation.  Reduce the
C                           estimated loss, TY6LSS, and try the
C                           supercritical profile until the loss is
C                           essentially zero.
C      fail      complete 3 Jump impossible. Profile is subcritical.
C
C      partial   fail     4 Should not happen.  May indicate a strange
C                           culvert barrel slope variation.  Reduce the
C                           estimated loss, TY6LSS, and try the
C                           supercritical profile until the loss is
C                           essentially zero.
C      partial  partial   5 Two outcomes:  no overlap- should not
C                           happen;  overlap-locate jump and estimate
C                           high-side depth.
C      partial  complete  6 Locate jump and estimate high-side depth.
C                           If no jump possible profile is subcritical.
C      complete fail      7 No jump possible and profile is supercritical
C
C     complete partial    8 Locate jump and estimate high-side depth.
C                           If no jump possible then profile is
C                           supercritical.
C     complete complete   9 Locate jump and estimate high-side depth.
C                           If no jump possible then something is
C                           probably wrong!
 
 
C     Attempt to compute super-critical profile to the end of
C     the barrel.
 
      SPZUP = YUP + ZBVEC(IU)
 
C     Distribute losses for computing the supercritical profile.
      CALL DISLSS
     I           (IU, ID, IAT3D, IAT6D, DUP, TY6LSS, XVEC,
     O            SEVEC)
 
      CALL SFPSPE
     I           (STDOUT, IU, ID, Q, SPZUP,
     O            ISP, SPZDN, SPFLAG)
      IF(SPFLAG.EQ.1) THEN
        ICSP = 2
      ELSE
        IF(ISP.GT.IU) THEN
          ICSP = 1
        ELSE
          ICSP = 0
        ENDIF
      ENDIF
C     Attempt to compute a subcritical profile.
 
      DH = 0.0
      SBZDN = YDN + ZBVEC(ID)
      CALL SFPSBE
     I           (STDOUT, IU, ID, DH, Q, SBZDN,
     O            ISB, SBZUP, SBFLAG)
C      WRITE(STDOUT,*) ' LOCJMP: SBFLAG=',SBFLAG
      IF(SBFLAG.EQ.1) THEN
        ICSB = 2
      ELSE
        IF(ISB.LT.ID) THEN
          ICSB = 1
        ELSE
          ICSB = 0
        ENDIF
      ENDIF
 
      IC = ICSB +3*ICSP + 1
C      WRITE(STDOUT,*) ' LOCJMP: ICSB=',ICSB,' ICSP=',ICSP,' IC=',IC
      GOTO(100, 100, 300, 100, 500, 600, 700, 800, 900), IC
 
        WRITE(STDOUT,50) IC
        STOP 'Abnormal stop. Errors found.'
 
 100    CONTINUE
C         Error condition- can indicate a bug in the software or
C         a condition in the culvert that cannot be computed.
          WRITE(STDOUT,56)
C         The indications are that the super critical profile should
C         have completed because the subcritical profile was either
C         a failure or incomplete.  The loss, TY6LSS, can only be
C         a rough approximation to the actual loss.  If the culvert
C         barrel is steep, the usual cause of failure for the
C         subcritical profile, then the supercritical profile should
C         be able to complete.  Thus reduce the loss by stages until
C         the supercritical profile is complete or the loss is
C         reduced to essentially zero.
 
          WRITE(STDOUT,57)  TY6LSS
 110      CONTINUE
            TY6LSS = 0.95*TY6LSS
            CALL DISLSS
     I                 (IU, ID, IAT3D, IAT6D, DUP, TY6LSS, XVEC,
     O                  SEVEC)
            CALL SFPSPE
     I                 (STDOUT, IU, ID, Q, SPZUP,
     O                  ISP, SPZDN, SPFLAG)
C            WRITE(STDOUT,*) ' TY6LSS=',TY6LSS,' ISP=',ISP,' ID=',ID
            IF(SPFLAG.EQ.1) THEN
C             The supercritical profile is complete.
              ICSP = 2
            ELSEIF(ISP.GT.ISB) THEN
C             The profiles overlap.
              ICSP = 1
            ELSE
              IF(TY6LSS/YUP.LT.0.001) THEN
C               Process has failed.  Super critical profile and
C               perhaps subcritical profile have major problems.
                WRITE(STDOUT,58)
                IERR = 1
                JMPLOC = 0
                YHIGH = 0.0
                YEND = 0.0
                PROTYP = ' '
                GOTO 9000
              ELSE
                GOTO 110
              ENDIF
            ENDIF
            IC = ICSB +3*ICSP + 1
C            WRITE(STDOUT,*) ' LOCJMP: ICSB=',ICSB,' ICSP=',ICSP,
C     A                        ' IC=',IC
            WRITE(STDOUT,59) TY6LSS
            GOTO(100, 100, 300, 100, 500, 600, 700, 800, 900), IC
 
            WRITE(STDOUT,50) IC
            STOP 'Abnormal stop. Errors found.'
 
 300    CONTINUE
C         Supercritical profile fails but subcritical is complete.
C          WRITE(STDOUT,*) ' LOCJMP: Found option 3'
          IERR = 0
          JMPLOC = 0
          YHIGH = 0.0
          YEND = YDN
          PROTYP= 'SUB'
          GOTO 9000
 
 500    CONTINUE
C         Possible mixed profile.
C          WRITE(STDOUT,*) ' LOCJMP: Found option 5'
          IF(ISP.LE.ISB) THEN
C           No overlap.  Unable to connect the two profiles.
C           Possible bug or error.
            IERR = 1
            JMPLOC = 0
            YHIGH = 0.0
            YEND = 0.0
            PROTYP = ' '
            GOTO 9000
          ELSE
C           Overlap. Find the jump location and high-side value
            PROTYP = 'MIXED'
            GOTO 2000
          ENDIF
 
 600    CONTINUE
C         Super critical profile partial and sub-critical complete.
C         Set the profile type to correct value if a jump
C         fails to exist.  Then find the  jump location and high-side
C         value.
C          WRITE(STDOUT,*) ' LOCJMP: Found option 6'
          PROTYP = 'SUB'
          YEND = YDN
          JMPLOC = 0
          YHIGH = 0.0
          GOTO 2000
 
 700    CONTINUE
C         Subcritical profile fails but supercritical is complete.
C          WRITE(STDOUT,*) ' LOCJMP: Found option 7'
          IERR = 0
          JMPLOC = 0
          YHIGH = 0.0
          YEND =  SPZDN - ZBVEC(ID)
          PROTYP= 'SUP'
          GOTO 9000
 
 800    CONTINUE
C         Supercritical profile complete but subcritical partial.
C         Set the profile type to correct value if a jump fails to
C         exist.  Then find the jump location and high-side value.
C          WRITE(STDOUT,*) ' LOCJMP: Found option 8'
          PROTYP = 'SUP'
          YEND =  SPZDN - ZBVEC(ID)
          YHIGH = 0
          JMPLOC = 0
          GOTO 2000
 
 900    CONTINUE
C         Both profiles complete.
C          WRITE(STDOUT,*) ' LOCJMP: Found option 9'
          PROTYP = 'BOTH'
          GOTO 2000
 
 2000 CONTINUE
C     Find the jump location and high-side value here.  ISP points
C     to the end of the supercritical profile and ISB points to
C     the end of the subcritical profile and ISP > ISB.   Start the
C     search at ISB.  Adjust for adding critical depth at end of
C     failed profiles.
 
      IF(ISP.LT.ID) THEN
        ISP = ISP + 1
      ENDIF
      IF(ISB.GT.IU) THEN
        ISB = ISB - 1
      ENDIF
 
      ISUP = 0
      ISUB = 0
      YSUPL = YVECSP(ISB)
      YT = MIN(DVEC(ISB),YSUPL)
      CALL XLKT22
     I           (NSEC(ISB),
     M            YT,
     O            AT, TT, DTT, JT, KT, DKT, BETT, DBETT, ALPT, DALPT,
     O            QCT)
      MSUPL = BETT*Q**2/AT + GRAV*JT
 
      YSUBL = YVECSB(ISB)
      YT = MIN(DVEC(ISB), YSUBL)
      CALL XLKT22
     I           (NSEC(ISB),
     M            YT,
     O            AT, TT, DTT, JT, KT, DKT, BETT, DBETT, ALPT, DALPT,
     O            QCT)
      CALL LKTJ
     I         (NSEC(ISB),
     M          YVECSB(ISB),
     O          JT)
      MSUBL = BETT*Q**2/AT + GRAV*JT
 
      IF(MSUPL.GT.MSUBL) THEN
        ISUP = 1
      ELSEIF(MSUPL.LT.MSUBL) THEN
        ISUB = 1
      ENDIF
 
C      WRITE(STDOUT,52)
C      WRITE(STDOUT,54) ISB, YSUPL, MSUPL, YSUBL, MSUBL
      IJUMP = 0
      DO 130 I=ISB+1,ISP
 
        YSUPR = YVECSP(I)
        YT = MIN(DVEC(I), YSUPR)
        CALL XLKT22
     I             (NSEC(I),
     M              YT,
     O              AT, TT, DTT, JT, KT, DKT, BETT, DBETT, ALPT, DALPT,
     O              QCT)
        MSUPR = BETT*Q**2/AT + GRAV*JT
 
        YSUBR = YVECSB(I)
        YT = MIN(DVEC(I), YSUBR)
        CALL XLKT22
     I             (NSEC(I),
     M              YT,
     O              AT, TT, DTT, JT, KT, DKT, BETT, DBETT, ALPT, DALPT,
     O              QCT)
        CALL LKTJ
     I           (NSEC(I),
     M            YSUBR,
     O            JT)
        MSUBR = BETT*Q**2/AT + GRAV*JT
 
C        WRITE(STDOUT,54) I, YSUPR, MSUPR, YSUBR, MSUBR
 
        IF(MSUPL.GE.MSUBL.AND.MSUPR.LE.MSUBR) THEN
C         There is a jump in the interval.
C         Find relative distance from the upstream end of the
C         current distance increment
 
          P = (MSUPL - MSUBL)/(MSUPL - MSUBL - MSUPR + MSUBR)
          YSUB = YSUBL + P*(YSUBR - YSUBL)
          IJUMP = I
C          GOTO 131
        ENDIF
        IF(MSUPL.GT.MSUBL) THEN
          ISUP = 1
        ELSEIF(MSUPL.LT.MSUBL) THEN
          ISUB = 1
        ENDIF
 
        MSUPL = MSUPR
        YSUPL = YSUPR
 
        MSUBL = MSUBR
        YSUBL = YSUBR
 
 130  CONTINUE
 131  CONTINUE
      IF(IJUMP.EQ.0) THEN
C       No jump even with overlap.
        IF(PROTYP.EQ.'BOTH') THEN
          IF(ISUP.EQ.1) THEN
            IF(ISUB.EQ.0) THEN
C             Super critical profile always prevails.
              PROTYP = 'SUP'
              YEND = YVECSP(ID)
            ELSE
C             Pattern does not show a jump but there is a
C             change in relationship.  Take supercritical flow
C             as prevailing.
              PROTYP= 'SUP'
              YEND = YVECSP(ID)
            ENDIF
          ELSE
            IF(ISUB.EQ.0) THEN
             WRITE(STDOUT,*) ' *BUG:XXX* LOCJMP: ISUB=0 & ISUP=0'
              STOP 'Abnormal stop. Errors found.'
            ELSE
C             Subcritical profile always prevails
              PROTYP = 'SUB'
              YEND = YVECSB(ID)
            ENDIF
          ENDIF
        ENDIF
        IERR = 0
        JMPLOC = 0
        YHIGH = 0.0
      ELSE
        JMPLOC = IJUMP
        YHIGH = YSUB
        YEND = YDN
        IERR = 0
        PROTYP = 'MIXED'
      ENDIF
 
 
 9000 CONTINUE
C      WRITE(STDOUT,*) ' LOCJMP at exit: PROTYP=',PROTYP,' YHIGH=',YHIGH,
C     A                ' YEND=',YEND,' JMPLOC=',JMPLOC,' IERR=',IERR
      RETURN
      END
C
C
C
      REAL FUNCTION   RQVSTW
     I                      (Q)
 
C     + + + PURPOSE + + +
C     Residual function for finding the flow in the culvert barrel
C     given the upstream head and a fixed tailwater level, Z43FIX.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Q
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Q      - Flowrate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rqvstw.cmn'
      INCLUDE 'typtrn.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CONFLG, IS, NSFLAG, SFLAG
      REAL DDN, DH, DUP
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL APPRO, DEGCON, FCD123, FNDCDE, SFPSBE, XLKTAL
C***********************************************************************
C      WRITE(OUTUN,*) ' RQVSTW ENTRY: Q=',Q,' Z43FIX=',Z43FIX
      Y3P = Z43FIX - ZB3
      DUP = DVEC(IUP)
      DDN = DVEC(IDN)
 
C     Establish the flows
 
      Q1 = Q + WFRD
      Q2 = Q
      Q3 = Q
 
C     Assume free surface flow at least part way.  Find solution
C     under this assumption.  The calling program unit must then
C     decide if that assumption is valid or not.  We cannot make the
C     decision here because it may be premature.  That is, this function
C     is called as part of an iterative solution.  In that process
C     there may be values of flow that cause the entrance soffit
C     to be submerged but that flow is not the final flow; it is only
C     one of many flows encountered in the process of finding
C     the final flow.  This final flow may be such that the entrance
C     soffit is free of water.
C     Find discharge coefficient without interpolation of any kind.
 
      C123 = FCD123(OUTUN, 3, CLASS, DUP, Z1T)
C      WRITE(OUTUN,*) ' RQVSTW: C123=',C123
C      IF(FQTYPE.EQ.1) THEN
C        AVH = A2
C      ELSE
        AVH = A3
C      ENDIF
      CD = DEGCON(C123, ABASE, AVH)
C      WRITE(OUTUN,*) ' RQVSTW: CD=',CD, ' AVH=',AVH
 
      IF(BETAF.EQ.-1.0) THEN
C       Computing submergence when the coef. of discharge and the
C       area for velocity head vary with the depth in the barrel exit.
C       Assumes that the flow is moving toward full flow at the
C       exit with full flow at the entrance.
        CD = CDF + (Y3P - Y3PF)*(C46 - CDF)/(DDN - Y3PF)
        AVH = AVHF + (Y3P - Y3PF)*(A2FULL - AVHF)/(DDN - Y3PF)
      ENDIF
 
      DH = (1.0/CD**2 - 1.0)*(Q/AVH)**2/GRAV2
C      WRITE(OUTUN,*) ' RQVSTW: DH=',DH,' Q=',Q
      VHL = 0.0
 
      CALL SFPSBE
     I           (OUTUN, IUP, IDN, DH, Q3, Z43FIX,
     O            IS, Z2, SFLAG)
 
C      WRITE(OUTUN,*) ' RQVSTW: SFPSBE:SFLAG=',SFLAG,' Z2=',Z2,' DH=',DH,
C     A             ' IS=',IS,' Z43FIX=',Z43FIX
      Y2 = Z2 - ZB2
      CALL XLKTAL
     I           (ADRXS2,
     M            Y2,
     O            A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2)
 
      IF(SFLAG.EQ.1) THEN
        SBFLAG = 0
      ELSE
C       Try continuing the profile towards critical depth.
        CALL FNDCDE
     I             (OUTUN, ADRXS2, Q3,
     M              Y2)
C        YCRIT = Y2
C        YSTART = YVECSB(IDN)
C        SLOPE = (ZBVEC(IUP) - ZBVEC(IDN))/ABS(XVEC(IDN) - XVEC(IUP))
C        SE = DH/ABS(XVEC(IUP) - XVEC(IDN))
C        SLOPE = SLOPE - SE
C        OFFSET = XVEC(IUP) - XVEC(IDN)
C        OFF = OFFSET
C        FFAC = 1.0 - KD(IUP+1)
C        CALL ENDSFP(OUTUN, ADRXS2, FFAC, YSTART, YCRIT, Q3, OFF,
C     A                  SLOPE, Y2)
C        WRITE(OUTUN,*) ' OFF=',OFF,' OFFSET=',OFFSET,' Y=',Y2

        Z2 = Y2 + ZB2
        CALL XLKTAL
     I             (ADRXS2,
     M              Y2,
     O              A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2)
C        SBFLAG = -1
        SBFLAG = 0
C        WRITE(OUTUN,*) ' RQVSTW: TEL entrance=',
C     A             Z2 + ALP2*(Q2/A2)**2/GRAV2,' IS=',IS,' IUP=',IUP
      ENDIF
C      WRITE(OUTUN,*) ' RQVSTW: Y3P=',Y3P,' Y2=',Y2,' CD=',CD,
C     A    ' AVH=',AVH,' VHL=',VHL
C     Compute the water level at section 1.
 
      IF(SBFLAG.EQ.0) THEN
        CALL APPRO
     I            (OUTUN, CD, VHL, WFRD,
     O             CONFLG, NSFLAG)
        IF(NSFLAG.EQ.1) THEN
C         Take special action in the calling routine.  Flow is
C         probably too large.
          SBFLAG = -2
        ENDIF
      ENDIF
C     Compute the residual
 
      RQVSTW = Z1 - Z1T 
C      WRITE(OUTUN,*) ' AT EXIT RQVSTW=',RQVSTW,' Y2=',Y2,
C     A            ' SBFLAG=',SBFLAG
C      WRITE(OUTUN,*) ' Z1=',Z1,' Z1T=',Z1T
      RETURN
      END
C
C
C
      REAL FUNCTION   RTY0
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute the residual for the approach reach when
C     the conditions at section 1 are given and a sub-critical
C     solution at section 2 is sought.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rty0c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IS
      REAL APPFAC, ARATIO, CDIN, DH, VH1, VH2, YT
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, LKTA, SFPTY1, XLKTAL
C***********************************************************************
      Y2 = Y
      YT = MIN(Y, DUP0)
      CALL XLKTAL
     I           (ADRXS2,
     M            YT,
     O            A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2)
 
      Z2 = ZB2  + Y2
      VH1 = (Q1/A1)**2/GRAV2
      VH2 = (Q2/A2)**2/GRAV2
 
      C123 = FCD123(OUTUN0, 3, CLASS0, DUP0, Z1T0)
 
      CD = DEGCON(C123, A1T0, A3)
 
      CDIN = CD
      IF(A1.LE.A2) THEN
C       EXPANSION(NEGATIVE ACCELERATION) IN FLOW INSTEAD OF CONTRACTION
        ARATIO = A1/A2
        IF(ARATIO.GT.0.95) THEN
C         Interpolate coefficients to make the transition between
C         the two cases smooth.
          CD = CDIN + 20.0*(1.0 - ARATIO)*(0.98 - CDIN)
          APPFAC = 20.0*(1.0 - ARATIO)*APPEXP
        ELSE
          CD = 0.98
          APPFAC = APPEXP
        ENDIF
C        CONF = 0
        RTY0 = (ALP1 - APPLOS)*VH1 + ZB1 + Y1 -
     A    (ALP2*VH2 + ZB2 + Y2  + APPLEN*Q1*Q2/(K1*K2)
     B     + APPFAC*(ALP1*VH1 - ALP2*VH2))
      ELSE
C        CONF = 1
        CD = CDIN
        RTY0 = ALP1*VH1 + ZB1 + Y1 -
     A    (ALP2*VH2 + ZB2 + Y2  + APPLEN*Q1*Q2/(K1*K2) + APPLOS*VH1)
      ENDIF
 
C      WRITE(OUTUN0,*) ' AT Y1=',Y1,' RTY0',RTY0,' CONF=',CONF
 
C     COMPUTE THE PROFILE TO THE END OF THE BARREL
      DH = 0.0
      CALL SFPTY1
     I           (OUTUN0, IUP0, IDN0, DH, Q2, Z2,
     O            IS, Z3, SFLAG0)
      IF(SFLAG0.EQ.0) THEN
        RETURN
      ENDIF
C     Estimate type 3 losses.  Only rough estimate is possible or
C     needed.
      Y3 = Z3 - ZB3
C     Define new value of area at section 3.
      YT = MIN(Y3,DDN0)
      CALL LKTA
     I         (ADRXS3,
     M          YT,
     O          A3)
      DH = (1.0/CD**2 - 1.0)*(Q3/A3)**2/GRAV2
 
 100  CONTINUE
        CALL SFPTY1
     I             (OUTUN0, IUP0, IDN0, DH, Q2, Z2,
     O              IS, Z3, SFLAG0)
        IF(SFLAG0.EQ.0) THEN
          DH = 0.8*DH
          GOTO 100
        ENDIF
      Y3 = Z3 - ZB3
C     Define new value of area at section 3.
      YT = MIN(Y3,DDN0)
      CALL LKTA
     I         (ADRXS3,
     M          YT,
     O          A3)
      RETURN
      END
C
C
C
      REAL FUNCTION   RTY2
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute residual for type 2 flow
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'rty2c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FTYPE, IS, NSFLAG, SFLAG
      REAL DH
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL APPRO, DEGCON, FCD123, SFPSBE, XLKT22, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT(/,' *WRN:572* Unable to compute subcritical profile for ',
     A  'Type 2.', /,11X,'Adjusting critical depth and trying ',
     B  ' again.',/,11X,'Y3=',F10.3,' Q3=',F10.2,' IS=',I5,' SFLAG=',
     C    I5)
 53   FORMAT(/,' *WRN:573*  Unable to compute approach reach for ',
     A  'Type 2.',/,11X,'Increasing critical depth and trying ',
     B  'again.',/,11X,'Y3=',F10.3,' Q3=',F10.2)
C***********************************************************************
C      WRITE(STD6,*) ' RTY2: Y=',Y
      SBFLAG = 0
      FTYPE = 2
      Y3 = Y
C     COMPUTE CRITICAL FLOW AT SECTION 3
C     GET VALUES AT SECTION 3
      CALL XLKT22
     I           (ADRXS3,
     M            Y3,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3,
     O            Q3C)
 
      Q3 = Q3C
      Q2 = Q3
 
C      WRITE(STD6,*) ' RTY2: Y3=',Y,' Q3C=',Q3
C     ADD IN THE FREE FLOW OVER THE ROADWAY
      Q1 = Q3 + WFRDF
      Z3 = ZB3 + Y3
 
C     DEFINE THE COEF OF DISCHARGE
 
      C123 = FCD123(OUTUN, FTYPE, CLASS, DVEC(IUP), Z1TRUE)
 
C     MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C     CHANNEL CONTRACTION
 
      CD = DEGCON(C123, A1TRUE, A3)
      DH = (1.0/CD**2 - 1.0)*(Q3/A3)**2/GRAV2
      VHL = 0.0
 
C      WRITE(OUTUN,*) ' RTY2: CD=',CD,' A for VH=',A3,' DH=',DH
C      WRITE(OUTUN,*) ' RTY2: Q3=',Q3,'  Y3=',Y3
      CALL SFPSBE
     I           (OUTUN, IUP, IDN, DH, Q3, Z3,
     O            IS, Z2, SFLAG)
C      WRITE(OUTUN,*) ' RTY2: AFTER SFPSBE: Z2=',Z2,' Z3=',Z3
      IF(IS.NE.IUP) THEN
C       UNABLE TO COMPUTE SUBCRITICAL PROFILE TO ENTRANCE
        IF(SFLAG.EQ.0) THEN
          RTY2 = 1.
          WRITE(OUTUN,52) Y3, Q3, IS, SFLAG
          SBFLAG = 2
          RETURN
        ENDIF
      ENDIF
 
C     DEFINE VALUES AT SECTION 2 BECAUSE SFPSBE DOES NOT USE
C     XS2COM
 
      Y2 = YVECSB(IUP)
      IF(Y2.GE.MAXARG) THEN
C       Take type 2 flow to be impossible.  Set flag for the
C       RGF routine.
        SBFLAG = 3
        RTY2 = -1.E31
        RETURN
      ENDIF
      CALL XLKTAL
     I           (ADRXS2,
     M            Y2,
     O            A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2)
 
      CALL APPRO
     I          (OUTUN, CD, VHL, WFRDF,
     O           CONF, NSFLAG)
 
C      WRITE(OUTUN,*) ' RTY2: NSFLAG=',NSFLAG
      NS = NSFLAG
      IF(NSFLAG.EQ.1) THEN
C       Approach reach problem.  Assume that the flow has been made
C       too small.
        RTY2 = 1.
        SBFLAG = 2
        WRITE(OUTUN,53) Y3, Q3
        RETURN
      ELSEIF(NSFLAG.EQ.2) THEN
C       Take type 2 flow to be impossible
        SBFLAG = 3
        RETURN
      ENDIF
 
      RTY2 =  Z1 - Z1TRUE 
 
C      WRITE(OUTUN,*) ' AT EXIT RTY2=',RTY2,' Q3=',Q3,' Y3=',Y3,
C     A               ' NSFLAG=',NSFLAG
      RETURN
      END
C
C
C
      SUBROUTINE   SUPSUB
     I                   (STDOUT, IU, ID, YUP, Y3LIM, Q, DUP, YCAT3,
     M                    TY6LSS,
     O                    ZAT3, ZAT43, JMPLOC)
 
C     + + + PURPOSE + + +
C     Find the conditions at the exit of the culvert required to
C     estimate drop to free flow for type 5 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IU, STDOUT
      REAL DUP, Q, TY6LSS, Y3LIM, YCAT3, YUP, ZAT3, ZAT43
      CHARACTER JMPLOC*32
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     YUP    - depth at vena contracta
C     Y3LIM  - full-flow-inducing depth at culvert exit for type 5 flow
C     Q      - Flowrate
C     DUP    - vertical diameter of culvert barrel at upstream end
C     YCAT3  - critical depth at section 3
C     TY6LSS - estimated type 6 loss
C     ZAT3   - elevation of water surface at section 3
C     ZAT43  - elevation of water surface at section 43
C     JMPLOC - character string giving the jump location
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs3com.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IERR, IJUMP
      REAL YEND, YHIGH
      CHARACTER PROTYP*8
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LOCJMP
C***********************************************************************
 
C     Seek jump in the barrel.
 
      CALL LOCJMP
     I           (STDOUT, IU, ID, YUP, YCAT3, Q, DUP,
     M            TY6LSS,
     O            YHIGH, IJUMP, YEND, IERR, PROTYP)
 
      IF(IERR.NE.0) THEN
        WRITE(STDOUT,*) ' *BUG:XXX SUPSUB: Error in locating jump.'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      IF(PROTYP.EQ.'MIXED') THEN
C       Jump is in the barrel.
        JMPLOC = 'Jump in the barrel'
        ZAT3 = YEND + ZB3
        ZAT43 = ZAT3
 
      ELSEIF(PROTYP.EQ.'SUP') THEN
C       Jump is in the exit.
        JMPLOC = 'Jump in the barrel exit'
        ZAT3 = YEND + ZB3
        ZAT43 = ZB3 + Y3LIM
C        WRITE(STDOUT,50)
 
      ELSE
C       Some problem here.  Should not be here.
        WRITE(STDOUT,*) ' *BUG:XXX* SUPSUB: Invalid profile outcome=',
     A        PROTYP
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      END
C
C
C
      SUBROUTINE   GETFRF
     I                   (Z,
     O                    ZSBRDF, FDRDW)
 
C     + + + PURPOSE + + +
C     Compute the free flow over the road and find its
C     submergence limit.   Results are stored in RDFCOM and EMBCOM
C     and in ZSBRDF.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL FDRDW, Z, ZSBRDF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Z      - elevation at the approach to the roadway
C     ZSBRDF - water surface elevation at section 43 that begins
C              submergence of flow over the roadway
C     FDRDW  - free drop for flow over roadway
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'embcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'embwrq.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL HRDFD, HRDFU, ZT
      REAL*8 FD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL EMBSUB, FRFEMB
C***********************************************************************
      IF(Z.GT.MINCRS) THEN
C       THERE IS FLOW OVER THE ROADWAY. FIND ITS FREE FLOW VALUE AND
C       THE DOWNSTREAM ELEVATION AT THE FREE FLOW BOUNDARY
 
C       Remember that FRFEMB was called.  This means that
C       road flow values must be reestablished.
 
C       FIND THE HEAD ON THE ROADWAY
 
        HRDFU = Z - MINCRS
 
        CALL FRFEMB
     I             (HRDFU, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I              OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I              HLMAX,
     M              HLFLAG, HPFLAG,
     O              ZT, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL, QRDFM,
     O              QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR, APPL,
     O              APPM, APPR, WL, WM, WR, AELL, AELM, AELR, WFRDF,
     O              MFRDF, EFRDF)
        CALL EMBSUB
     I             (MINLOC, MINCRS, NOFF, SURF, TOTHL, TOTHR, HRDFU,
     O              HRDFD, FD, ZSBRDF)
        FDRDW = FD 
 
C       ZSBRDF GIVES THE TAILWATER ELEVATION WHICH MUST BE
C       REACHED FOR SUBMERGENCE OF THE FLOW OVER THE ROADWAY TO OCCUR.
      ELSE
        HRDFU = 0.0
        WFRDF = 0.0
        MFRDF = 0.0
        EFRDF = 0.0
        ZSBRDF = MINCRS
        FDRDW = 0.0
      ENDIF
 
      RETURN
      END
