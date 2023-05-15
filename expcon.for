C
C
C
      REAL FUNCTION   FRLRES
     I                      (QT)
 
C     + + + PURPOSE + + +
C     Find froude number residual at the left section
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL QT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     QT     - Flow rate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'feccom.cmn'
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDHPL
C***********************************************************************
      Q = QT
      CALL FNDHPL
C      WRITE(STD6,*) ' FRLRES: Q=',Q,' QCL=',QCL,' HPL=',HPL
      IF(FLG.EQ.0.OR.FLG.EQ.5) THEN
C       COMPUTE THE FROUDE NUMBER RESIDUAL
        FRLRES = 1.0 - Q/QCL
      ELSE
C       SET TO A NEGATIVE VALUE TO FIND A SOLUTION.
 
        FRLRES = -1.0
      ENDIF
C      WRITE(STD6,*) 'FRLRES=',FRLRES
      RETURN
      END
C
C
C
      REAL FUNCTION   ECECHK()
 
C     + + + PURPOSE + + +
C     Compute a check on the energy balance for flow through
C     an expansion-contraction
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'feccom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL ECLOSS, FLOSS, HTR
 
C     + + + INTRINSICS + + +
      INTRINSIC MAX, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FACDC, GMEAN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FACDC, GMEAN
C***********************************************************************
C     COMPUTE THE MECHANICAL ENERGY LOSS DUE TO BOUNDARY FRICTION
 
      FLOSS = DX*Q**2/GMEAN(KL, KR, TGM)
 
C     COMPUTE THE MECHANICAL ENERGY LOSS DUE TO ADDITIONAL TURBULENCE
C     INTRODUCED BY EITHER EXPANSION OR CONTRACTION OF THE FLOW
 
      ECLOSS = Q*FACDC(Q*(SQRT(ALPHAR)/AR - SQRT(ALPHAL)/AL),
     A                    SMOOTH, KA, KD)*
     B          (SQRT(ALPHAR)/AR + SQRT(ALPHAL)/AL)/GRV2
 
      HTR = HPR + ALPHAR*Q**2/(GRV2*AR**2)
      ECECHK = HTL - HTR - ECLOSS - FLOSS
C     WRITE(STD6,*) ' ECECHK=',ECECHK
C     WRITE(STD6,*) ' HTL-HTR=', HTL - HTR
      ECECHK = ECECHK/MAX(HTL-HTR,EPSDIF)
C     WRITE(STD6,*) ' ECECHK=',ECECHK
      RETURN
      END
C
C
C
      REAL FUNCTION   FACDC
     I                     (X, SMOOTH, KA, KD)
 
C     + + + PURPOSE + + +
C     Smoothing function for losses based on velocity head
C     difference.  We smooth near a difference of zero and in
C     terms of inverse area(approximately).
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL KA, KD, SMOOTH, X
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     X      - Argument defining transition losses
C     SMOOTH - Smoothing parameter for expansion and contraction
C              losses
C     KA     - Acceleration loss.
C     KD     - Deceleration loss.
 
C     + + + LOCAL VARIABLES + + +
      REAL H, P
C***********************************************************************
      IF(SMOOTH.EQ.0.0) THEN
C       NO SMOOTHING OPTION
        IF(X.LE.0.0) THEN
          FACDC = -KD*X
        ELSE
          FACDC = KA*X
        ENDIF
      ELSE
C       FIT A CUBIC OVER INTERVAL (-SMOOTH, SMOOTH).  TURNS OUT TO BE
C       A PARABOLA!
        IF(X.LE.-SMOOTH) THEN
          FACDC = -KD*X
        ELSEIF(X.LT.SMOOTH) THEN
C         USE FITTED FUNCTION
          H = SMOOTH + SMOOTH
          P = (X + SMOOTH)/H
          FACDC = SMOOTH*(KD +P*(-2.*KD + P*(KD + KA)))
        ELSE
          FACDC = KA*X
        ENDIF
      ENDIF
C        WRITE(STD6,50) X, FACDC
C50    FORMAT(' FACDC EXIT: X=',1PE10.3,' FACDC=',E10.3)
      RETURN
      END
C
C
C
      REAL FUNCTION   GMEAN
     I                     (XA, YA, TA)
 
C     + + + PURPOSE + + +
C     Compute the  square of the generalized mean of two numbers
C     using the parameter TA.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL TA, XA, YA
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     XA     - One of two arguments for generalized mean
C     YA     - One of two arguments for generalized mean
C     TA     - Parameter defining generalized mean value
 
C     + + + LOCAL VARIABLES + + +
      DOUBLE PRECISION R, T, X, Y
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
C***********************************************************************
      IF(ABS(TA).LT.0.01) THEN
C       COMPUTE SQUARE OF GEOMETRIC MEAN AND RETURN
        GMEAN = XA*YA
      ELSE
        X = XA
        Y = YA
        T = TA
        R = X/Y
        GMEAN = (Y*(.5D0*R**T + .5D0)**(1.D0/T))**2
C        GMEAN = ((0.5D0*(X**T + Y**T))**(1.D0/T))**2
      ENDIF
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION   FHPL
     I                                (H)
 
C     + + + PURPOSE + + +
C     Compute the residual function defining the piezometric head
C     at the left end of the transition when the flow
C     is fixed and the piezometric head on the right is fixed.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      DOUBLE PRECISION H
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     H      - Piezometric head
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'feccom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL B, DALPHA, DB, DKL, DTL, JL, TL, X, Y
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FACDC, GMEAN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FACDC, GMEAN, XLKT22
C***********************************************************************
C      WRITE(STD6,*) ' FHPL: H=',H,' Q=',Q
 
C     THE ELEMENTS AT R ARE CONSTANT AND MUST BE KNOWN ON ENTRY--
C     CONTAINED IN FECCOM.FOR.
 
C     FIND ELEMENTS AT L FOR THE CURRENT PIEZOMETRIC HEAD
 
      HPL = H
      YL = HPL + HDATUM - ZBL
      IF(XTABL.GT.0) THEN
        CALL XLKT22
     I             (XTABL,
     M              YL,
     O              AL, TL, DTL, JL, KL, DKL, B, DB, ALPHAL, DALPHA,
     O              QCL)
C       WRITE(STD6,*) ' FHPL: HPL=',HPL,' YL=',YL,' QCL=',QCL
        HTL = HPL + ALPHAL*(Q/AL)**2/GRV2
      ELSE
C       RESERVOIR ON LEFT
        AL =1.E10
        ALPHAL = 0.0
        KL = 1.E10
        HTL = HPL
        QCL = 1.E20
      ENDIF
 
 
      X = Q*(SQRT(ALPHAR)/AR - SQRT(ALPHAL)/AL)
      Y = Q*(SQRT(ALPHAR)/AR + SQRT(ALPHAL)/AL)
 
      FHPL = (HTL - HPR - Q**2*(ALPHAR/(GRV2*AR**2) + DX/
     A       GMEAN(KL, KR, TGM)) -
     B       FACDC(X, SMOOTH, KA, KD)*Y/GRV2)/HTL
 
C      WRITE(STD6,50) HPL, FHPL
C50    FORMAT(' FHPL EXIT: HPL=',F10.4,' FHPL=',F10.4)
      RETURN
      END
C
C
C
      SUBROUTINE   FNDHPL()
 
C     + + + PURPOSE + + +
C     Find the piezometric head at left given the piezometric head
C     at right and the flow rate.  If no solution set FLG in FECCOM.
C     flow rate is in FECCOM.COM. Solution is returned in FECCOM also.
 
      IMPLICIT NONE

C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'feccom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL ARGL, ARGR, HPLEST, HPLMAX, HPLMIN, HTR, SUM, SUM1, YMAX
      DOUBLE PRECISION ALEFT, ARIGHT, DARG, FLEFT, FRIGHT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SNGL
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FMXARG
      DOUBLE PRECISION FHPL
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FHPL, FMXARG, SECANT
C***********************************************************************
C     USE EXHAUSTIVE SEARCH FOR SUBCRITICAL SOLUTION.
 
      HTR = HPR + ALPHAR*(Q/AR)**2/GRV2
 
      IF(XTABL.EQ.0) THEN
C       RESERVOIR ON THE LEFT.  A SUBCRITICAL SOLUTION ALWAYS EXISTS
C       BECAUSE THE FLOW WILL BE CRITICAL AT THE RIGHT AT MOST.
C       MINIMUM VALUE AT LEFT IS THE TOTAL ENERGY LINE ELEV AT RIGHT
 
 
        HPL = 1.05*HTR
        HPLEST = HTR
        CALL SECANT
     I             (HPLEST, EPSARG, EPSF, EPSABS, 30, HTR, 1.E5, FHPL,
     M              HPL,
     O              FLG)
        IF(FLG.GT.0) THEN
          WRITE(STD6,*) ' FNDHPL: NO SOLUTION WITH RESERVOIR ON LEFT.'
          STOP 'Abnormal stop. Errors found.' 
        ENDIF
      ELSE
C       DO A SEARCH ON VALID UPSTREAM HEADS FOR A SIGN CHANGE IN
C       FHPL
 
        YMAX = FMXARG(XTABL)
        HPLMAX = ZBL + YMAX - HDATUM
        IF(HPLMAX.GT.2.*HTR) THEN
          HPLMAX = 2.*HTR
        ENDIF
        HPLMIN = ZBL - HDATUM
        IF(HPLMIN.LT.0.0) THEN
          HPLMIN = 0.0
        ENDIF
        ARIGHT = HPLMAX
        FRIGHT = FHPL(ARIGHT)
        DARG = (HPLMAX - HPLMIN)/32.
 
C       WRITE(STD6,*) ' FNDHPL: ARIGHT=',ARIGHT,' FRIGHT=',FRIGHT
        SUM = 1.E10
        DO 100 I=1,31
          ALEFT = ARIGHT - DARG
          FLEFT = FHPL(ALEFT)
C       WRITE(STD6,*) ' FNDHPL: ALEFT=',ALEFT,' FLEFT=',FLEFT
 
          IF(FLEFT*FRIGHT.LE.0D0) GOTO 105
          SUM1 = ABS(FLEFT) + ABS(FRIGHT)
          IF(SUM1.LT.SUM) THEN
            SUM = SUM1
            ARGL = ALEFT
            ARGR = ARIGHT
          ENDIF
          ARIGHT = ALEFT
          FRIGHT = FLEFT
 100    CONTINUE
C       NO SIGN CHANGE HERE.
C       TRY FOR A SOLUTION BECAUSE THE SEARCH MAY HAVE MISSED
C       WRITE(STD6,*) ' FNDHPL: NO SIGN CHANGE.'
        HPL = (ARGR)
        CALL SECANT
     I         ((ARGL), EPSARG, EPSF, EPSABS, 30, HPLMIN, HPLMAX, FHPL,
     M              HPL,
     O              FLG)
        IF(FLG.GT.0) THEN
          FLG = 3
C         WRITE(STD6,*) ' FNDHPL: NO SOLUTION WITH NO SIGN CHANGE.'
          RETURN
        ENDIF
 
        GOTO 110
 105    CONTINUE
C       SOLUTION EXISTS.  USE SECANT TO FIND IT IN THE INTERVAL.
        HPL = SNGL(ARIGHT)
        CALL SECANT
     I    (SNGL(ALEFT), EPSARG, EPSF, EPSABS, 30, HPLMIN, HPLMAX, FHPL,
     M              HPL,
     O              FLG)
        IF(FLG.GT.0) FLG = 4
 110    CONTINUE
C       CHECK FOR SUBCRITICAL SOLUTION
        IF(FLG.EQ.0.AND.Q/QCL.GE.1.0)  THEN
          FLG = 5
C          WRITE(STD6,*) ' FNDHPL: FROUDE > 1: Q=',Q,' QCL=',QCL
        ENDIF
      ENDIF
C      WRITE(STD6,50) HTL, HPL, HPR, FLG
C50    FORMAT('0FNDHPL EXIT: HTL=',F10.4,' HPL=',F10.4,' HPR=',F10.4,
C     A         ' FLG=',I5)
      RETURN
      END
C
C
C
      SUBROUTINE   FRFTRN
     I                   (STDOUT,
     O                    EFLAG, QFREE, CNTLOC)
 
C     + + + PURPOSE + + +
C     Find free flow and upstream piezometric head at free flow.
C     also report the control location.  the upstream piezometric
C     head value, HPL, is in FECCOM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDOUT
      REAL QFREE
      CHARACTER CNTLOC*4
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     QFREE  - Free flow
C     CNTLOC - String giving location of critical control
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'feccom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER DNFLAG, FLAG, I
      REAL ALD, ALPHAD, BR, DALPHA, DBR, DKR, DQ, DTR, FL, FR, HPLD,
     A     HPRD, HTLD, JR, QC, QD, QL, QMAX, QR, TR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FRLRES
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDHPL, FRLRES, REGFLT, XLKT22
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT('0*ERR:551 No control found for transition in FRFTRN.')
C***********************************************************************
C      WRITE(STDOUT,*) ' ENTERING FRFTRN: HPR=', HPR
 
C     RESERVOIRS AT R MUST BE TREATED SEPERATELY
 
      IF(XTABR.GT.0) THEN
C       CROSS SECTION PRESENT AT R
 
C       GET ELEMENTS FOR THE SECTION AT THE RIGHT.
 
        YR = HPR + HDATUM - ZBR
        CALL XLKT22
     I             (XTABR,
     M              YR,
     O              AR, TR, DTR, JR, KR, DKR, BR, DBR, ALPHAR, DALPHA,
     O              QCR)
 
C       TRY TO FIND SUBCRITICAL VALUE AT LEFT FOR CRITICAL FLOW AT
C       RIGHT SECTION
 
        DNFLAG = 0
        Q = QCR
        CALL FNDHPL
        IF(FLG.EQ.0) THEN
C         SUBCRITICAL SOLUTION WAS FOUND.
 
C         CONTROL IS ON THE RIGHT AND WE HAVE A SOLUTION.
 
          CNTLOC = '  DN'
          QFREE = Q
C         SAVE VALUES FOR LATER
          HPRD = HPR
          HPLD = HPL
          ALPHAD = ALPHAL
          ALD = AL
          HTLD = HTL
          QD = Q
          DNFLAG = 1
C         WRITE(STDOUT,*) ' CONTROL ON RIGHT: Q=',Q
        ELSE
C         NO SOLUTION ASSUMING CONTROL ON THE RIGHT
C          WRITE(STDOUT,54)
        ENDIF
C         WRITE(STDOUT,*) ' CHECK FOR CONTROL ON LEFT.'
C         SEARCH FOR SOLUTION AT OTHER END.  MAY BE MORE THAN ONE
C         CONTROL POINT.
C         SEARCH ON FLOW UNTIL WE HAVE
C         CLOSELY BRACKETED CRITICAL FLOW.  THEN USE REGULA FALSI
C         TO FIND THE CONDITIONS AT CRITICAL FLOW.
 
          QMAX = QCR
          DQ = QMAX/64.0
C         WRITE(STDOUT,*) ' QMAX=',QMAX,' DQ=',DQ
C         WE WILL START AT Q = 0.0 BECAUSE THE FROUDE NUMBER RESIDUAL
C         IS KNOWN TO BE 1.0.
 
          QL = 0.0
          FL = 1.0
C         SEARCH FOR A SUBCRITICAL SOLUTION
          DO 100 I=1,64
            QR = QL + DQ
            FR = FRLRES(QR)
C           WRITE(STDOUT,*) 'FRFTRN  Q=',Q,' FR=',FR,' HPL=',HPL,' QCL=',QCL
            IF(FR*FL.LE.0.0) GOTO 105
            QL = QR
            FL = FR
 100      CONTINUE
 
C         WRITE(STDOUT,*) ' NO SIGN CHANGE FOUND IN SEARCH'
          IF(DNFLAG.EQ.1) THEN
C           SOLUTION DOWNSTREAM EXISTS
            HPR   = HPRD
            HPL   = HPLD
            ALPHAL= ALPHAD
            AL    = ALD
            HTL   = HTLD
            Q    = QD
            RETURN
          ELSE
            WRITE(STDOUT,50)
            EFLAG = 1
            RETURN
          ENDIF
 105      CONTINUE
C         FOUND SIGN CHANGE.  SOLVE FOR THE CRITICAL CONDITION.
 
C          WRITE(STDOUT,60) QL, FL, QR, FR
C60    FORMAT(' FRFTRN BEFORE REGFLT: QL=',F10.2,' FL=',F10.4,
C     A      ' QR=',F10.2,' FR=',F10.4)
 
          CALL REGFLT
     I               (EPSARG, EPSF, FRLRES,
     M                QL, QR, FL, FR,
     O                QC, FLAG)
 
C         WRITE(STDOUT,62) QC, FL, QL, QR, FR
C62    FORMAT(' FRFTRN AFTER REGFLT: QC=',F10.3,' F(QC)=',1PE10.3,
C     A        ' QL=',0PF10.3,' QR=',F10.3,' FR=',F10.3)
 
C         COMPARE THE RESULTS WITH THOSE OBTAINED, IF ANY, ASSUMING
C         CONTROL DOWNSTREAM.  TAKE THE SET WITH THE SMALLEST FLOW.
 
          IF(DNFLAG.EQ.1) THEN
            IF(QD.LT.Q) THEN
              HPR   = HPRD
              HPL   = HPLD
              ALPHAL= ALPHAD
              AL    = ALD
              HTL   = HTLD
              Q    = QD
              CNTLOC = '  DN'
              RETURN
            ENDIF
          ENDIF
 
          IF(ABS(FL).GT.0.01) THEN
            WRITE(STDOUT,50)
            EFLAG = 1
            RETURN
          ENDIF
          CNTLOC = '  UP'
          QFREE = QC
C        ENDIF
 
      ELSE
C       RESERVOIR ON THE RIGHT.  SEEK VALUES USING SPECIAL METHOD
C       BECAUSE CRITICAL FLOW CAN ONLY OCCUR ON THE LEFT.
 
        WRITE(STDOUT,*) ' CODE FOR RESERVOIR ON THE RIGHT NOT DONE YET.'
      WRITE(STDOUT,*) ' A LARGE CROSS SECTION ON RIGHT YIELDS THE SAME'
        WRITE(STDOUT,*) ' RESULTS AS A RESERVOIR ON THE RIGHT.'
        STOP 'Abnormal stop. Errors found.'
 
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   FNDECT
     I                   (STDOUT, XTABU, ZBU, XTABD, ZBD, DXX, SMTH,
     I                    TGMEAN, GRV, KACC, KDEC, NHD, HDAT, HDVEC,
     I                    NFRAC, PFQVEC,
     O                    EFLAG, HUFVEC, HUMAT, QFVEC)
 
C     + + + PURPOSE + + +
C     Compute a 2-d table for an expansion-contraction.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NFRAC, NHD, STDOUT, XTABD, XTABU
      REAL DXX, GRV, HDAT, HDVEC(PMXNHU), HUFVEC(PMXNHU),
     A     HUMAT(PMXNHU,PMXFRC), KACC, KDEC, PFQVEC(PMXFRC),
     B     QFVEC(PMXNHU), SMTH, TGMEAN, ZBD, ZBU
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     XTABU  - Address of upstream cross section function table
C     ZBU    - Bottom elevation at upstream section
C     XTABD  - Address of downstream cross section function table
C     ZBD    - Bottom elevation at downstream section
C     DXX    - Length of the transition.  Used for boundary friction
C               estimates
C     SMTH   - Smoothing parameter for expansion and contraction
C              losses
C     TGMEAN - Generalized mean parameter
C     GRV    - value of acceleration due to gravity
C     KACC   - Loss coefficient when flow is accelerating
C     KDEC   - Loss coefficient when flow is decelerating
C     NHD    - Number of downstream heads
C     HDAT   - Datum for measuring head
C     HDVEC  - Vector of prescribed downstream heads
C     NFRAC  - Number of fractions for defining partial free flows
C     PFQVEC - Partial free flow vector
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     HUFVEC - Vector of computed upstream heads at free flow
C     HUMAT  - Matrix of computed upstream heads
C     QFVEC  - Free flow vector
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'feccom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL EBAL, MAXEE, MAXPFQ, QFREE
      CHARACTER CNTLOC*4, QTYPE*4,CQ*8
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL ECECHK
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL ECECHK, FNDHPL, FRFTRN
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,'  Dwns piezometric head=',F10.4,' Ups head at',
     A       ' free flow=',F10.4)
 52   FORMAT('  Fraction   Upstream    Downstream      Flow',
     A     '     Location  Sense',/,
     B  '  of Free   Piezometric  Piezometric   through       ',
     C  'of      of',/,
     D  '    Flow       Head         Head      Transition  Control',
     E  '   Flow')
 54   FORMAT(' ',1X,F7.4,2X,2X,F8.4,5X,F8.4,3X,3X,A8,1X,2X,A4,
     A           4X,1X,A4,3X,1X)
 56   FORMAT(' Max. error in energy relative to energy loss=',1PE8.1,
     A    '  at PFQ=',0PF10.4)
C***********************************************************************
C     ESTABLISH THE VALUES IN THE COMMON BLOCK WHICH REMAIN CONSTANT
 
      GRV2 = 2.*GRV
      XTABL = XTABU
      XTABR = XTABD
      ZBL = ZBU
      ZBR = ZBD
      SMOOTH = SMTH
      DX = DXX
      KA = KACC
      KD = KDEC
      HDATUM = HDAT
      TGM = TGMEAN
 
C     SET VALUES FOR RESERVOIRS
 
      IF(XTABL.EQ.0) THEN
        ALPHAL = 0.0
        AL = 1.E10
      ENDIF
      IF(XTABR.EQ.0) THEN
        ALPHAR = 0.0
        AR = 1.E10
      ENDIF
 
      DO 2000 I=1,NHD
 
C       INITIALIZE THE MAXIMUM ENERGY ERROR VALUE
 
        MAXEE = -1.0
 
        HPR = HDVEC(I)
 
C       COMPUTE FREE FLOW AND UPSTREAM HEAD THROUGH THE TRANSITION
 
        CALL FRFTRN
     I             (STDOUT,
     O              EFLAG, QFREE, CNTLOC)
 
 
        IF(EFLAG.NE.0) THEN
          RETURN
        ENDIF
C       FOR EACH OF THE PARTIAL FREE FLOWS(EXCLUDING 0.00 AND 1.0)
C       COMPUTE THE UPSTREAM PIEZOMETRIC HEAD
 
        WRITE(STDOUT,50) HPR, HPL
 
 
        WRITE(STDOUT,52)
 
        IF(ALPHAL/AL**2.GT.ALPHAR/AR**2) THEN
          QTYPE = ' EXP'
        ELSE
          QTYPE = ' CON'
        ENDIF
        EBAL = ECECHK()
        IF(EBAL.GT.MAXEE) THEN
          MAXEE = EBAL
          MAXPFQ = PFQVEC(NFRAC)
        ENDIF
        CALL VAR_DECIMAL(QFREE,
     O                   CQ)
        WRITE(STDOUT,54) PFQVEC(NFRAC), HPL, HPR, CQ, CNTLOC, QTYPE
 
        CNTLOC = 'BOTH'
 
        HUMAT(I,1) = HPR
        HUMAT(I,NFRAC) = HPL
        HUFVEC(I) = HPL
        QFVEC(I) = QFREE
 
        DO 500 J=NFRAC-1,2,-1
          Q = QFREE*PFQVEC(J)
 
          IF(XTABR.GT.0) THEN
C           CROSS SECTION EXISTS AT R
            CALL FNDHPL
 
          ELSE
C           NO CROSS SECTION EXISTS AT R. RESERVOIR ON THE RIGHT.
            WRITE(STDOUT,*) ' RESERVOIR ON RIGHT NOT YET DONE.'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          IF(ALPHAL/AL**2.GT.ALPHAR/AR**2) THEN
            QTYPE = ' EXP'
          ELSE
            QTYPE = ' CON'
          ENDIF
          EBAL = ECECHK()
          IF(EBAL.GT.MAXEE) THEN
            MAXEE = EBAL
            MAXPFQ = PFQVEC(J)
          ENDIF
          CALL VAR_DECIMAL(Q,
     O                     CQ)
          WRITE(STDOUT,54) PFQVEC(J), HPL, HPR, CQ, CNTLOC, QTYPE
          HUMAT(I,J) = HPL
 
 500    CONTINUE
 
        WRITE(STDOUT,56) MAXEE, MAXPFQ
 
 2000 CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE READ_EXPCON_ITEMS1(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     M                            EFLAG, IDLEN,  
     O                            CHAR4, XTAB, X, ZB)

C     Get the items of data from the first set of EXPCON input

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A         EFLAG, IDLEN, XTAB
      REAL X, ZB
      CHARACTER LINE*(*), CHAR4*4

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:741* Only ',I3,' items given in ',
     A   'EXPCON-1 description line.  Need at least four items.')

C***********************************************************************
      IF(NITEM.LT.4) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

      N = 1
C     Process the location
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      CHAR4 = TPC

C     Process the table id
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      LKEY = LENSTR(KEY)
      IDLEN = MAX(IDLEN, LKEY)
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 ITAB)
        XTAB = ITAB
      ELSE
        XTAB = 0
      ENDIF

C     Process the station 
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') X

C     Process the invert elevation
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') ZB

      RETURN
      END
C
C
C
      SUBROUTINE READ_EXPCON_ITEMS2(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     M                            EFLAG, IDLEN,  
     O                            CHAR4, TAB, KACC, KDEC, LAB)

C     Get the items of data from the second set of EXPCON input

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A         EFLAG, IDLEN, TAB
      REAL KACC, KDEC
      CHARACTER LINE*(*), CHAR4*4, LAB*50

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:742* Only ',I3,' items given in ',
     A   'EXPCON-2 description line.  Need at least five items.')

C***********************************************************************
      IF(NITEM.LT.5) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

      N = 1
C     Process the direction
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      CHAR4 = TPC

C     Process the table id
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      LKEY = LENSTR(KEY)
      IDLEN = MAX(IDLEN, LKEY)
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 ITAB)
        TAB = ITAB
      ELSE
        TAB = 0
      ENDIF

C     Process the the acceleration loss
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') KACC

C     Process the deceleration loss
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F10.0)') KDEC

C     Process the label
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      LAB = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    LAB)

      RETURN
      END
C
C
C
      SUBROUTINE   EXPCON
     I                   (STDIN, STDOUT, STDTAB, GRV,
     M                    TABDIR,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Compute a two-d table for flow through an expansion-contraction.
C     use the table with downstream piezometric head and flow as
C     argument.  The table gives the upstream piezometric head.

      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL GRV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     GRV    - value of acceleration due to gravity
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER MAXN
      PARAMETER(MAXN=5)
      INTEGER DUTAB, I, NFRAC, NHD, TABTYP, UDTAB, WFLAG, XTABL, XTABR,
     A        IDLEN, COLWIDTH, NITEM, ITEM_START(MAXN), ITEM_END(MAXN)
      REAL DX, HDATUM, HDOLD, HDVEC(PMXNHU), HUFVEC(PMXNHU),
     A     HUMAT(PMXNHU,PMXFRC), KACCDU, KACCUD, KDECDU, KDECUD,
     B     PFQVEC(PMXFRC), POWER, QFVEC(PMXNHU), QOLD, SMOOTH, TGMEAN,
     C     XL, XR, ZBL, ZBR, zrhufd
      real*8 easting, northing, dnull, eastingl, northingl, eastingr,
     a       northingr
      CHARACTER CHAR4*4, CHAR5*5, HEAD*80, LABDU*50, LABUD*50, LINE*80,
     A          JUST*5, TABID*16, IDOUT*32, TABIDU*16, TABIDD*16,
     b  zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, FLOAT, MAX
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      CHARACTER GETTOK*4, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKCFC, CHKTAB, FNDECT, GETTBN, GETTOK, inline, TABCHK,
     A         TWDOUT, GET_TABID, READ_EXPCON_ITEMS1,
     B         READ_EXPCON_ITEMS2, LENSTR
 
      data dnull/-33d6/

C     + + + INPUT FORMATS + + +
 1    FORMAT(A4,I5,F5.0,F10.0)
 2    FORMAT(A4,I5,2F5.0,1X,A50)
 4    FORMAT(7X,F10.0)
 5    FORMAT(6X,F10.0)
 16   FORMAT(A5,1X,F10.0)
 24   FORMAT(A80)
 26   FORMAT(A5,1X,I5)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' ',A80)
 52   FORMAT(' ',A4,A,F5.0,F10.2)
 54   FORMAT(' ',A4,A,2F5.2,1X,A50)
 56   FORMAT(' SMOOTH=',F10.5)
 57   FORMAT(' GENERALIZED MEAN VALUE PARAMETER=',F10.2,/,
     A  '  MEAN VALUE IS THEREFORE:')
 66   FORMAT(' ',A5,'=',F10.2)
 70   FORMAT(/,' *ERR:618* Upstream section out of order in EXPCON.')
 71   FORMAT(/,' *ERR:619* Downstream section ouf of order in EXPCON.')
 72   FORMAT(/,' *ERR:620* Loss coefficient < 0 or > 1:',F8.3)
 74   FORMAT(/,' *ERR:621* U to D coefficients and table not first.')
 75   FORMAT(/,' *ERR:622* D to U coefficients and table not second.')
 76   FORMAT(/,'  Transition length is zero.  Boundary friction',
     A   ' losses ignored.')
 78   FORMAT(/,' COMPUTING TABLE FOR FLOW UP -> DOWN')
 80   FORMAT(/,' COMPUTING TABLE FOR FLOW DOWN -> UP')
 82   FORMAT(/,' *ERR:623* EXPCON requires at least one cross section.',
     A   '  None were found.')
 88   FORMAT(/,' Checking cross section tables for possible critical',
     A       ' flow',/,5X,' or celerity problems.')
 89   FORMAT(/,' *WRN:550* EXPCON command may not converge because',
     A       ' critical',/,10X,' flow/celerity decreases in TabId= ',
     B         A)
 90   FORMAT('  Processing EXPCON TabIds= ',A,' and ', A)
 92   FORMAT(/,' ',A80)
 94   FORMAT(/,' ',A5,'=',I5)
 95   FORMAT(/,' ','INPUT COMPLETE. BEGIN COMPUTATIONS')
 96   FORMAT(/,' ','DATUM FOR HEADS IS:',F10.2)
 98   FORMAT(/,' *WRN:519* Free flow=',F10.3,' decreases at downstream',
     A        ' head=',F10.4)
 99   FORMAT(/,' *ERR:624* DOWNSTREAM HEADS NON-INCREASING AT:',F10.2)
C***********************************************************************
      TABTYP = 14
      JUST = 'RIGHT'

c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

C     INPUT THE BASIC DATA FROM USER.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24) HEAD
      WRITE(STDOUT,50) HEAD
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24) HEAD
      WRITE(STDOUT,50) HEAD
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, HEAD, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)
 
C     Set the  column width
      COLWIDTH =  ITEM_END(2) - ITEM_START(2) + 1
      
C     INPUT THE UPSTREAM TABLE INFORMATION
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IDLEN = 0
      CALL READ_EXPCON_ITEMS1(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, IDLEN,  
     O                         CHAR4, XTABL, XL, ZBL)

      TABID = GET_TABID(XTABL)
      IDOUT = ' '
      IDOUT(COLWIDTH-IDLEN+1:COLWIDTH) = TABID(1:IDLEN)
      WRITE(STDOUT,52) CHAR4, IDOUT(1:COLWIDTH), XL, ZBL
      CHAR4 = GETTOK(CHAR4)
      IF(CHAR4.NE.'UP') THEN
        WRITE(STDOUT,70)
        EFLAG = 1
      ELSE
C       CHECK THE TABLE FOR EXISTENCE
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              XTABL,
     O              EFLAG)
      ENDIF
 
C     INPUT THE DOWNSTREAM TABLE INFORMATION
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IDLEN = 0
      CALL READ_EXPCON_ITEMS1(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, IDLEN,  
     O                         CHAR4, XTABR, XR, ZBR)
      TABID = GET_TABID(XTABR)
      IDOUT = ' '
      IDOUT(COLWIDTH-IDLEN+1:COLWIDTH) = TABID(1:IDLEN)
      WRITE(STDOUT,52) CHAR4, IDOUT(1:COLWIDTH), XR, ZBR
      CHAR4 = GETTOK(CHAR4)
      IF(CHAR4.EQ.'DN'.OR.CHAR4.EQ.'DOWN') THEN
C       CHECK THE TABLE FOR EXISTENCE
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              XTABR,
     O              EFLAG)
      ELSE
        WRITE(STDOUT,71)
        EFLAG = 1
      ENDIF
 
C     MAKE SURE THAT AT LEAST ONE CROSS SECTION EXISTS
      IF(XTABL.EQ.0.AND.XTABR.EQ.0) THEN
        WRITE(STDOUT,82)
        EFLAG = 1
      ENDIF
 
C     INPUT THE COEFFICIENT AND OUTPUT TABLES
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24) HEAD
      WRITE(STDOUT,50) HEAD
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24) HEAD
      WRITE(STDOUT,50) HEAD
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, HEAD, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)
C     Set the  column width
      COLWIDTH =  ITEM_END(2) - ITEM_START(2) + 1
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IDLEN = 0
      CALL READ_EXPCON_ITEMS2(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, IDLEN,  
     O                         CHAR4, UDTAB, KACCUD, KDECUD, LABUD)
C      READ(LINE,2,ERR=991) CHAR4, UDTAB, KACCUD, KDECUD, LABUD
      TABIDU = GET_TABID(UDTAB)
      IDOUT = ' '
      IDOUT(COLWIDTH-IDLEN+1:COLWIDTH) = TABIDU(1:IDLEN)
      WRITE(STDOUT,54) CHAR4, IDOUT(1:COLWIDTH), KACCUD, KDECUD, LABUD
      CHAR4 = GETTOK(CHAR4)
 
 
C     CHECK TABLE NUMBER FOR VALIDITY
      CALL TABCHK
     I           (STDOUT, PMXTAB,
     M            UDTAB, TABDIR, EFLAG)
C     CHECK LOSS COEFFICIENTS FOR VALIDITY.
      IF(KACCUD.GT.1.0.OR.KACCUD.LT.0.0) THEN
        WRITE(STDOUT,72) KACCUD
        EFLAG = 1
      ENDIF
      IF(KDECUD.GT.1.0.OR.KDECUD.LT.0.0) THEN
        WRITE(STDOUT,72) KDECUD
        EFLAG = 1
      ENDIF
      IF(CHAR4.NE.'UD'.AND.CHAR4.NE.'UTOD') THEN
        WRITE(STDOUT,74)
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IDLEN = 0
      CALL READ_EXPCON_ITEMS2(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, IDLEN,  
     O                         CHAR4, DUTAB, KACCDU, KDECDU, LABDU)
C      READ(LINE,2,ERR=991) CHAR4, DUTAB, KACCDU, KDECDU, LABDU
      TABIDD = GET_TABID(DUTAB)
      IDOUT = ' '
      IDOUT(COLWIDTH-IDLEN+1:COLWIDTH) = TABIDD(1:IDLEN)
      WRITE(STDOUT,54) CHAR4, IDOUT(1:COLWIDTH), KACCDU, KDECDU, LABDU
      CHAR4 = GETTOK(CHAR4)
 
C     CHECK TABLE NUMBER FOR VALIDITY
      IF(DUTAB.GT.0) THEN
        CALL TABCHK
     I             (STDOUT, PMXTAB,
     M              DUTAB, TABDIR, EFLAG)
C       CHECK LOSS COEFFICIENTS FOR VALIDITY.
        IF(KACCDU.GT.1.0.OR.KACCDU.LT.0.0) THEN
          WRITE(STDOUT,72) KACCDU
          EFLAG = 1
        ENDIF
        IF(KDECDU.GT.1.0.OR.KDECDU.LT.0.0) THEN
          WRITE(STDOUT,72) KDECDU
          EFLAG = 1
        ENDIF
        IF(CHAR4.NE.'DU'.AND.CHAR4.NE.'DTOU') THEN
          WRITE(STDOUT,75)
          EFLAG = 1
        ENDIF
      ENDIF
 
      WRITE(*,90)  TABIDU(1:LENSTR(TABIDU)), TABIDD(1:LENSTR(TABIDD))
C     INPUT THE SMOOTHING CONSTANT FOR LOSSES WHEN AREAS ARE NEARLY
C     EQUAL.
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) SMOOTH
      WRITE(STDOUT,56) SMOOTH
 
C     INPUT THE MEAN VALUE PARAMETER
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,5,ERR=991) TGMEAN
      WRITE(STDOUT,57) TGMEAN
      IF(TGMEAN.LT.-1.0) THEN
        WRITE(STDOUT,*) '  MINIMUM < MEAN VALUE < HARMONIC MEAN.'
      ELSEIF(TGMEAN.EQ.-1.0) THEN
        WRITE(STDOUT,*) '   HARMONIC MEAN.'
      ELSEIF(TGMEAN.LT.0.0) THEN
        WRITE(STDOUT,*) '  HARMONIC MEAN < MEAN VALUE < GEOMETRIC MEAN.'
      ELSEIF(TGMEAN.EQ.0.0) THEN
        WRITE(STDOUT,*) '   GEOMETRIC MEAN.'
      ELSEIF(TGMEAN.LT.1.0) THEN
        WRITE(STDOUT,*) '  GEOMETRIC MEAN < MEAN VALUE < ARITHMETIC',
     A                  ' MEAN.'
      ELSEIF(TGMEAN.EQ.1.) THEN
        WRITE(STDOUT,*) '   ARITHMETIC MEAN.'
      ELSEIF(TGMEAN.GT.1.0) THEN
        WRITE(STDOUT,*) '  ARITHMETIC MEAN < MEAN VALUE <  MAXIMUM.'
      ENDIF
      WRITE(STDOUT,*) ' '
 
C     INPUT THE DOWNSTREAM HEADS AND THE FACTORS CONTROLLING THE
C     DISTRIBUTION OF PARTIAL FREE FLOWS.
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,26,ERR=991) CHAR5, NFRAC
      WRITE(STDOUT,94) CHAR5, NFRAC
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) CHAR5, POWER
      WRITE(STDOUT,66) CHAR5,POWER
 
C     COMPUTE THE PARTIAL FREE FLOWS
 
      DO 200 I=1,NFRAC
        PFQVEC(I) = (FLOAT(I-1)/FLOAT(NFRAC-1))**POWER
 200  CONTINUE
 
C     INPUT THE DOWNSTREAM HEAD SEQUENCE
 
      I = 1
      HDOLD = -1.0
 300  CONTINUE
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
        READ(LINE,'(F10.0)',ERR=991) HDVEC(I)
        WRITE(STDOUT,'(1X,F10.2)') HDVEC(I)
        IF(HDVEC(I).LE.0.0) THEN
          NHD = I - 1
          GOTO 310
        ELSE
          IF(HDVEC(I).LE.HDOLD) THEN
            WRITE(STDOUT,99) HDVEC(I)
            EFLAG = 1
          ENDIF
          HDOLD = HDVEC(I)
          I = I + 1
          IF(I.GT.PMXNHU) THEN
            WRITE(STDOUT,
     A '('' *ERR:548* MORE THAN '',I5,'' UPSTREAM HEADS'')')
     B           PMXNHU
            I = PMXNHU
            EFLAG = 1
          ENDIF
 
          GOTO 300
        ENDIF
 310  CONTINUE
 
      IF(EFLAG.NE.0) RETURN
 
C     CHECK FOR THE CRITICAL FLOW AND CELERITY IN THE TABLE.
      WRITE(STDOUT,88)
      IF(XTABL.GT.0) THEN
        CALL CHKCFC
     I             (GRV, STDOUT, XTABL,
     O              WFLAG)
        IF(WFLAG.NE.0) THEN
          WRITE(STDOUT,89) GET_TABID(GETTBN(XTABL))
        ENDIF
      ENDIF
      IF(XTABR.GT.0) THEN
        CALL CHKCFC
     I             (GRV, STDOUT, XTABR,
     O              WFLAG)
        IF(WFLAG.NE.0) THEN
          WRITE(STDOUT,89) GET_TABID(GETTBN(XTABR))
        ENDIF
      ENDIF
 
c     Check for existence of easting-northing in the cross-section tables
      if(easting <= dnull) then
c       User did not input a value of easting.
c       Get values from each of the tables.
        call get_east_north(
     i                      stdout, xtabl,
     o                      eastingl, northingl)
        call get_east_north(
     i                      stdout, xtabr,
     o                      eastingr, northingr)
c       If both are defined, take the average, else, take 
c       the value that is defined.  
        if(eastingl > dnull) then
          if(eastingr > dnull) then
c           Both are defined.
            easting = 0.5*(eastingl + eastingr)
            northing = 0.5*(northingl + northingr)
          else
c           Only left values are defined.
            easting = eastingl
            northing = northingl
          endif
        else
c         left values undefined
          if(eastingr > dnull) then
c           right values are defined
            easting = eastingr
            northing = northingr
          endif
        endif
      endif

C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.
 
      WRITE(STDOUT,95)
 
C     SELECT THE DATUM FOR HEADS.
      HDATUM = MAX(ZBL, ZBR)
      WRITE(STDOUT,96) HDATUM
 
C     COMPUTE THE DISTANCE BETWEEN THE SECTIONS
 
      DX = ABS(XL - XR)
      IF(DX.EQ.0.0) THEN
        WRITE(STDOUT,76)
      ENDIF
 
      WRITE(STDOUT,78)
 
      CALL FNDECT
     I           (STDOUT, XTABL, ZBL, XTABR, ZBR, DX, SMOOTH, TGMEAN,
     I            GRV, KACCUD, KDECUD, NHD, HDATUM, HDVEC, NFRAC,
     I            PFQVEC,
     O            EFLAG, HUFVEC, HUMAT, QFVEC)
 
C     CHECK FOR DECREASING FREE FLOW
      QOLD = 0.0
      DO 320 I=1,NHD
        IF(QFVEC(I).LT.QOLD) THEN
          WRITE(STDOUT,98) QFVEC(I), HDVEC(I)
        ENDIF
        QOLD = QFVEC(I)
 320  CONTINUE
 
C     OUTPUT THE TABLE
      zrhufd = 0.0
      CALL TWDOUT
     I           (STDOUT, STDTAB, UDTAB, LABUD, NHD, NFRAC, QFVEC,
     I            HDVEC, PFQVEC, HUMAT, HDATUM,
     I            TABTYP,'  EXPCON', zrhufd,
     i            zone, hgrid, vdatum, unitsys, basis,
     i            easting, northing,
     O            EFLAG)
 
      IF(DUTAB.GT.0) THEN
        WRITE(STDOUT,80)
 
        CALL FNDECT
     I             (STDOUT, XTABR, ZBR, XTABL, ZBL, DX, SMOOTH, TGMEAN,
     I              GRV, KACCDU, KDECDU, NHD, HDATUM, HDVEC, NFRAC,
     I              PFQVEC,
     O              EFLAG, HUFVEC, HUMAT, QFVEC)
 
 
C       CHECK FOR DECREASING FREE FLOW
        QOLD = 0.0
        DO 330 I=1,NHD
          IF(QFVEC(I).LT.QOLD) THEN
            WRITE(STDOUT,98) QFVEC(I), HDVEC(I)
          ENDIF
          QOLD = QFVEC(I)
 330    CONTINUE
 
C       OUTPUT THE TABLE
        zrhufd = 0.0
        CALL TWDOUT
     I             (STDOUT, STDTAB, DUTAB, LABDU, NHD, NFRAC, QFVEC,
     I              HDVEC, PFQVEC, HUMAT, HDATUM,
     I              TABTYP,'  EXPCON', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)
 
      ENDIF
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
