C     ***********
C     *         *
C     * QCULCIR
C     *         *
C     ***********

      SUBROUTINE QCULCIR (YX, D, A, R, HD)

C     THIS SUBROUTINE COMPUTES AREA, HYDRAULIC RADIUS AND HYDRAULIC DEPTH 
C     FOR THE NON-OVERLAP OF TWO CIRCLES


C     Dummy arguments
      REAL YX, D, A, R, HD

C**********************************************************************
      Z = D/2.
      DEP = YX
      IF(DEP.GE.D) DEP = 0.9999*D
      IF(DEP.LE.0.) DEP = 0.0001*D 
      DAB = DEP - Z 
      Y = DAB/Z 
      Y1 = ABS(Y) 
C *****  ARCSIN APPROXIMATION 
      PHIY = 1.570796 +(-0.214512 + (0.0878763 + (-0.0449589 + 
     A            (0.0193499 -0. 00433777*Y1)*Y1)*Y1)*Y1)*Y1 
      ANGLE = 1.570796 - (1.0 - Y1)**0.5*PHIY 
      IF (Y) 10,20,20 
   10 ANGLE = -ANGLE
   20 DAC = ANGLE + 1.570796
      A = (DAB*(D*DEP - DEP*DEP)**0.5) + (Z*Z*DAC)
      IF(A.LE.0.) A = 0.000001 
      T = 2.*(Z*Z - DAB*DAB)**0.5 
      WP = D*DAC
      HD = A/T
      R = A/WP
      RETURN
      END
C     ***********
C     *         *
C     * QCULREC
C     *         *
C     ***********

      SUBROUTINE QCULREC (Y, D, W, A, R, HD) 
C 
C     This subroutine computes area, depth, and hydraulic radius for
C     rectangular channel .
C 

C     Dummy arguments. 

      REAL Y, D, W, A, R, HD

C     Local values

      REAL DEP
C***********************************************************************
      DEP = Y 
      IF(DEP.GE.D) DEP = 0.9999*D
      IF(DEP.LE.0.) DEP = 0.0001*D 
      A = W*DEP 
      R = W*DEP/(W + 2.*DEP)
      HD = DEP
      RETURN
      END 
C     ***********
C     *         *
C     * QCULWEI
C     *         *
C     ***********

      SUBROUTINE QCULWEI (QW, CODE,
     A           HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL, L, D, 
     B           W, N, K, C, HW, TW, KE, WB, WE, SWB, SWE, CW, A, AW) 

C     This subroutine computes culvert flow when weir control exists
C     (flashboard and riser) 


C     Dummy arguments

      INTEGER BARREL, GTYPE
      REAL INEL, L, N, K, KE, KWE
      CHARACTER CODE*3

C**********************************************************************
      H = HWE - WE
     
      QW = CW*WB*H**1.5 
      IF(TWE.GT.WE) QW = QW*(1. - ((TWE - WE)/(HWE - WE))**1.5)**0.385 
      IF(HWE.GT.SWE) THEN
         QSW = CW*SWB*(HWE - SWE)**1.5
         IF(TWE.GT.SWE) QSW = QSW*
     A            (1. - ((TWE - SWE)/(HWE - SWE))**1.5)**0.385 
         QW = QW + QSW
      ENDIF 
C *****  MODIFY ENTRANCE LOSS COEFFICIENT FOR CULVERT FLOW COMPUTATION
C *****      ASSUME ENTRANCE LOSS COEFFICIENT OVER WEIR TO BE 0.10
      KWE = 0.10
      AW = WB*H 
      IF(HWE.GT.SWE) AW = AW + (HWE - SWE)*SWB 
      IF(BARREL.EQ.0) CALL QCULCIR(HW,D,A,R,HD)
      IF(BARREL.EQ.1) CALL QCULREC(HW,D,W,A,R,HD)
      KE = KWE*(A/AW)**2 + K
      CODE(2:2) = 'U' 
      CODE(3:3) = 'F' 
      IF (TWE.GT.WE) CODE(3:3) = 'S'
      RETURN
      END 

C     ***********
C     *         *
C     * QCULGAT
C     *         *
C     ***********

      SUBROUTINE QCULGAT (GTYPE, BARREL, HW, GGAP, D, W, K, KE, A, AG) 
C 
C     THIS SUBROUTINE COMPUTES GATE CONTROL PARAMETERS FOR GENERAL 
C     CULVERT ROUTINE <QCULV> 
C 
C     Dummy arguments
      
      REAL K,KE 
      INTEGER BARREL,GTYPE

C***********************************************************************
      IF (GTYPE.EQ.0.AND.BARREL.EQ.0) THEN
         Z = D/2. 
         GGAP1 = Z - GGAP/2.
         CALL QCULCIR(GGAP1, D, AG1, R, HD) 
         A1 = 3.1416*Z**2 
         AG = A1 - 2.*AG1 
         CALL QCULCIR(HW, D, A, R, HD)
         KE = K*(A/AG)**2 
      ELSEIF (GTYPE.GE.1.AND.BARREL.EQ.0) THEN 
         CALL QCULCIR(GGAP, D, AG, R, HD) 
         CALL QCULCIR(HW,D,A,R,HD)
         IF(GTYPE.NE.2) KE = K*(A/AG)**2 
      ELSE
         AG = GGAP*W
         CALL QCULREC(HW, D, W, A, R, HD)
         IF(GTYPE.NE.2) KE = K*(A/AG)**2 
      ENDIF
C 
      RETURN
      END 

C     ***********
C     *         *
C     * QCULDIT
C     *         *
C     ***********

      SUBROUTINE QCULDIT (QC, CODE,
     A           HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL, L, D,
     B           W, N, K, C, HW, TW, KE) 
C 
C     This subroutine computes open channel flow through culvert
C 
      INTEGER BARREL, GTYPE
      REAL INEL, L, N, K, KE, FROUDE 
      CHARACTER CODE*3
C***********************************************************************
      A = -999.0
      AG = -99.0
      Q1 = 9999.0 
      Q2 = 9999.0 
      Q3 = 9999.0 
      IF(GTYPE.EQ.2) THEN
         KE = 0.25*KE 
      ELSEIF(HW .LT. 1.3*GGAP) THEN 
         KE = 0.25*K
      ELSE
         CALL QCULGAT(GTYPE, BARREL, HW, GGAP, D, W, K, KE, A, AG) 
         KE = 0.25*KE 
      ENDIF 
C ***** 
C *****   INLET CONTROL 
C ***** 
      IF (BARREL.EQ.1) GO TO 30 
      YC1 = 0.75*HW 
C *****  ITERATION FOR INLET CRITICAL DEPTH YC1 
      DO 10 I = 1,250 
         CALL QCULCIR(YC1, D, A1, R1, HD1)
         HW1 = YC1 + (1. + KE)*HD1/2. 
         DEV = (HW - HW1)
         IF (ABS(DEV).LE.0.005) GO TO 20 
         YC1 = YC1 + DEV*0.1
   10 CONTINUE
C ***** 
   20 Q = SQRT(32.2*HD1)*A1 
      GO TO 40
   30 YC1 = 2.*HW/(3. + KE) 
      CALL QCULREC(YC1, D, W, A1, R1, HD1) 
      Q = SQRT(32.2*HD1)*A1 
   40 Q1 = Q
      V1 = Q1/A1
      FROUDE = V1/SQRT(32.2*HD1)
      IF(TWE .GT. YC1 + INEL) GO TO 120 
      SC = (Q*N/(1.49*A1*R1**.6667))**2 
      SB = (INEL - OUTEL)/L 
      IF(SB.GE.SC) THEN
         QC = Q1
         CODE(1:1) = 'H'
         CODE(2:2) = 'U'
         CODE(3:3) = 'F'
         RETURN 
      ENDIF 
C ***** 
C *****  FREE FALL AT OUTLET
C ***** 
      YC2 = 0.8*YC1 
      Y1 = YC1 + OUTEL - INEL 
      V1 = 1E-20
      V2 = 1E-20
C *****  ITERATION FOR OUTLET CRITICAL DEPTH YC2
      DO 80 I = 1,1000
         IADJ = I/100 + 1 
         IF(Y1.GE.HW) Y1 = 0.999*HW
         IF(YC2.GE.(Y1 + INEL - OUTEL)) YC2 = (Y1 + INEL - OUTEL)*0.999
         IF (YC2.LE.0.0) YC2 = 0.001
         IF (Y1.LE.YC2 + OUTEL - INEL) Y1 = (YC2 + OUTEL - INEL)*1.001
         IF(BARREL.EQ.0) THEN
            CALL QCULCIR(Y1, D, A1, R1, HD1)
            CALL QCULCIR(YC2, D, A2, R2, HD2) 
         ELSE 
            CALL QCULREC(Y1, D, W, A1, R1, HD1)
            CALL QCULREC(YC2, D, W, A2, R2, HD2) 
         ENDIF
         CONVEY = 1.49/N*SQRT(A1*R1**.6667*A2*R2**.6667)
         F = (Y1 + INEL) - (YC2 + OUTEL)
         S = F/L
         IF(S.LE.0.) S = 0.0001
         Q = CONVEY*SQRT(S) 
         V1 = Q/A1
         V2 = Q/A2
         FROUDE = V2/SQRT(32.2*HD2) 
         HD22 = V2*V2/32.2
         DEV2 = HD22 - HD2
         IF (ABS(DEV2).LE.0.005) GO TO 70
         YC2 = YC2 + DEV2*0.1/IADJ
         GO TO 80 
C *****  ITERATION FOR Y1 
   70    Y11 = HW - (1. + KE)*V1**2/64.4 
         DEV1 = Y11 - Y1
         IF (ABS(DEV1).LE.0.005) GO TO 90
         Y1 = Y1 + DEV1*0.1/IADJ
   80 CONTINUE
C ***** 
   90 Q2 = Q
      CODE(2:2) = 'U' 
      IF (AG.LT.A) CODE(2:2) = 'C'
      CODE(3:3) = 'F' 
      IF(I.GE.1000) GO TO 170 
      IF (TW.GT.YC2) GO TO 120
      IF(Q2.LE.Q1) THEN
         QC = Q2
         CODE(1:1) = 'T'
      ELSE
         QC = Q1
         CODE(1:1) = 'H'
      ENDIF 
      RETURN
C ***** 
C *****   TAILWATER EFFECT
C ***** 
  120 Y1 = 1.01*(TWE - INEL)
      CODE(2:2) = 'U' 
      IF (AG.LT.A) CODE(2:2) = 'C'
      CODE(3:3) = 'S' 
      Y2 = TW 
      V1 = 1E-20
      V2 = 1E-20
C *****  ITERATION FOR Y1 
      DO 130 I = 1,250
         IADJ = I/50 + 1
         IF(Y1.GE.HW) Y1 = HW*0.999
         IF(Y1.LE.(TWE - INEL)) Y1 = (TWE - INEL)*1.001
         IF(BARREL.EQ.0) THEN
            CALL QCULCIR(Y1, D, A1, R1, HD1)
            CALL QCULCIR(Y2, D, A2, R2, HD2)
         ELSE 
            CALL QCULREC(Y1, D, W, A1, R1, HD1)
            CALL QCULREC(Y2, D, W, A2, R2, HD2)
         ENDIF
         CONVEY = 1.49/N*SQRT(A1*R1**.6667*A2*R2**.6667)
         F = (Y1 + INEL) -(Y2 + OUTEL) 
         S = F/L
         IF(S.LE.0.) S = 0.0001
         Q = CONVEY*SQRT(S) 
         V1 = Q/A1
         V2 = Q/A2
         Y11 = HW - (1 + KE)*V1**2/64.4 
         DEV = Y11 - Y1 
         IF (ABS(DEV).LE.0.005) GO TO 140
         Y1 = Y1 + DEV*0.1/IADJ 
  130 CONTINUE
C ***** 
  140 Q3 = Q
      IF(I.GE.250) GO TO 170
      IF(Q3.LE.Q1) THEN
         QC = Q3
         CODE(1:1) = 'T'
      ELSE
         QC = Q1
         CODE(1:1) = 'H'
      ENDIF 
      RETURN
C *****  ERROR DETECTION
  170 QC = MIN(Q1,Q2,Q3)
      CODE(1:1) = '?' 
      RETURN
      END
C     ***********
C     *         *
C     * QCULORI
C     *         *
C     ***********
 
      SUBROUTINE QCULORI(QC, CODE,
     A           HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL, L, D,
     B           W, N, K, C, HW, TW, KE) 
C 
C     This subroutine computes culvert flow when only part of pipe is 
C     full.
C 
      INTEGER BARREL,GTYPE
      REAL INEL,L,N,K,KE,KF 
      CHARACTER CODE*3

C**********************************************************************
      CALL QCULGAT(GTYPE, BARREL, HW, GGAP, D, W, K, KE, A, AG)
      H = HW - 0.6*GGAP 
      Q = C*AG*SQRT(64.4*H) 
      Q1 = Q
      IF(HW .LT. 1.3*D) THEN 
         QC = Q1
         CODE(1:1)  =  'O'
         CODE(2:2) = 'U'
         CODE(3:3) = 'F'
         RETURN 
      ENDIF 

C     Note by ddf:  Appears to assume that the piezometric level
C     at the pipe outlet is always 0.7*D.  Therefore, the tailwater
C     level must exceed the elevation at 0.7*D at the outlet
C     before the flow defined by an the full pipe relationship is affected.

      H = MIN(HWE - (OUTEL + 0.7*D), HWE - TWE) 
      IF(BARREL.EQ.0) R = D/4. 
      IF(BARREL.EQ.1) R = (D*W)/(2*D + W)
      KF = 29.1*N**2*L/R**1.3333
      Q = A*SQRT(64.4*H/(1. + KE + KF)) 
      Q2 = Q
      IF(Q2.LE.Q1) THEN
         QC = Q2
         CODE(1:1) = 'P'
         IF(HWE - (OUTEL + 0.7*D).LT. HWE - TWE) THEN
           CODE(3:3) = 'F'
         ELSE
           CODE(3:3) = 'S'
         ENDIF
      ELSE
         QC = Q1
         CODE(1:1) = 'O'
C        Improper class given here in orginal. Change to free flow.
C         CODE(3:3) = 'S'
         CODE(3:3) = 'F'
      ENDIF 
      CODE(2:2) = 'U' 
      IF (AG.LT.A) CODE(2:2) = 'C'
      RETURN
      END
C     ***********
C     *         *
C     * QCULPIP
C     *         *
C     ***********

      SUBROUTINE QCULPIP (QC, CODE,
     A           HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL, L, D, 
     B           W, N, K, C, HW, TW, KE) 
C 
C     THIS SUBROUTINE COMPUTES CULVERT FLOW WHEN PIPE IS FULL 
C 
      INTEGER BARREL, GTYPE
      REAL INEL, L, N, K, KE, KF 
      CHARACTER CODE*3
C***********************************************************************
      CALL QCULGAT(GTYPE, BARREL, HW, GGAP, D, W, K, KE, A, AG)
      H = HWE - TWE 
      IF(HW.LE.(1.3*GGAP)) KE = 0.25*KE
      IF(BARREL.EQ.0) R = D/4. 
      IF(BARREL.EQ.1) R = (D*W)/(2*D + W)
      KF = 29.1*N**2*L/R**1.3333
      Q = A*SQRT(64.4*H/(1. + KE + KF)) 
      QC = Q
      CODE(1:1) = 'F' 
      CODE(2:2) = 'U' 
      IF (AG.LT.A) CODE(2:2) = 'C'
      CODE(3:3) = 'S' 
      RETURN
      END 

C     ***********
C     *         *
C     * SFWMD_QCULV  *
C     *         *
C     ***********
 
      SUBROUTINE SFWMD_QCULV(CODE, GTYPE, BARREL, BOARD, C, CW, D, HWE,
     A                  INEL, K, L, N, OUTEL, SWB, SWE, TWE,
     B                  W, WB, QA)

C     Dummy arguments 

      CHARACTER CODE*3
      INTEGER GTYPE, BARREL
      REAL  BOARD, C, CW, D,  HWE, INEL, K, L,     
     A     N, OUTEL, SWB, SWE, TWE,  W, WB, QA
C************************************************************************
C                      I D E N T I F I C A T I O N
C                      ---------------------------
C SUBROUTINE NAME          - QCULV (ORIGINALLY PROGRAM <CULVERT>)

C SUBROUTINE DESCRIPTION   - COMPUTES THE DISCHARGE THROUGH GATED 
C                            CULVERTS 

C ORIGINAL PROGRAMMER      - ANDREW FAN 
C                            WATER RESOURCES DIVISION 
C                            RESOURCE PLANNING DEPARTMENT 
C                            SOUTH FLORIDA WATER MANAGEMENT DISTRICT

C REFERENCE                - FAN, ANDREW, PROGRAM DOCUMENTATION,
C                            A GENERAL PROGRAM TO COMPUTE FLOW THROUGH
C                            GATED CULVERTS, TECHNICAL MEMORANDUM,
C                            SFWMD, OCTOBER 1985. 

C ORIGINAL PRODUCTION DATE - 10/31/83 (ORIGINAL PROGRAM DONE 10/8/82) 

C************************************************************************

C                 L I S T   O F   IDENTIFIERS
C                 --------------------------------- 

C VARIABLE   TYPE            DESCRIPTION OF USE 
C --------   -------------   -------------------------------------------
C 
C BARREL     INTEGER         CULVERT SHAPE (0 = CIRCLE; 1 = BOX)
C BOARD      REAL            BOARD ELEVATION (GATE-OPENING), IN FEET
C 
C C          REAL            GATE CONSTRICTION FLOW COEFFICIENT 
C                               (RANGE = 0.6 (SQUARE-EDGE) TO 
C                                        0.9 (ROUND-EDGE))
C CODE       CHAR*3          FLOW REGIME INDICATOR, OUTPUT VARIABLE,
C                               WHERE:  
C                               CODE(1:1) = TYPE CULVERT FLOW 
C                                              F = FULL PIPE FLOW 
C                                              H = HEADWATER-CONTROLLED 
C                                                     (OPEN CHANNEL - 
C                                                        INLET CONTROL) 
C                                                        FLOW 
C                                              O = ORIFICE-CONTROLLED 
C                                                     (PARTIAL PIPE - 
C                                                        NON-INLET
C                                                        CONTROL) FLOW
C                                              P = PARTIAL PIPE FLOW -
C                                                     INLET CONTROL 
C                                              T = TAILWATER CONTROLLED 
C                                                     (OPEN CHANNEL - 
C                                                        NON-INLET
C                                                        CONTROL) FLOW
C                                              W = WEIR CONTROL 
C                                              ? = CHECK "ERRORF" FOR 
C                                                     POSSIBLE PROBLEM
C                                                     IN CALCULATING
C                                                     FLOW
C                               CODE(2:2) = CONTROLLED ("C") OR 
C                                              UNCONTROLLED ("U") FLOW
C                                              INDICATOR
C                               CODE(3:3) = SUBMERGED ("S") OR
C                                              FREE ("F") FLOW INDICATOR
C CW         REAL            WEIR COEFFICIENT (APPROX. 3.3) 
C 
C D          REAL            CULVERT DIAMETER (CIRCLE) OR CULVERT 
C                               VERTICAL HEIGHT (BOX) 
C 
C GTYPE      INTEGER         GATE TYPE CODE (0 = CIRCLE; 1 = RECTANGLE; 
C                               2 = WEIR) 
C 
C HWE        REAL            HEADWATER (UPSTREAM STAGE) ELEVATION,
C                               IN FEET M.S.L.

C INEL       REAL            INLET INVERT ELEVATION (FEET M.S.L.) 

C K          REAL            INLET LOSS COEFFICIENT 
C                               (FLUSH HEADWALL </= 0.5;
C                                  ROUNDED INLET = 0.04 TO 0.20;
C                                  PROJECTING INLET = 0.8 TO 0.9) 

C L          REAL            CULVERT LENGTH 

C N          REAL            MANNINGS' "N" COEFFICIENT
C                               (CONCRETE = 0.012 +/- 0.002;
C                                  CMP = 0.021 +/- 0.006) 

C OUTEL      REAL            OUTLET INVERT ELEVATION (FEET M.S.L.)
C 
C QA         REAL            CALCULATED DISCHARGE (OUTPUT VARIABLE) 

C SWB        REAL            LENGTH OF OVERFLOW PORTION OF RISER
C SWE        REAL            ELEVATION OF OVERFLOW PORTION OF RISER 

C TWE        REAL            TAILWATER (DOWNSTREAM STAGE) ELEVATION,
C                               IN FEET M.S.L.

C W          REAL            CULVERT DIAMETER (CIRCLE) OR CULVERT 
C                               HORIZONTAL WIDTH (BOX)
C WB         REAL            WEIR WIDTH 

C***********************************************************************

C* THIS PROGRAM COMPUTES THE DISCHARGE THROUGH GATED CULVERTS 

C* CODE: W=WEIR CONTROL; F=FULL PIPE FLOW; O = ORIFICE CONTROL; 
C*       P = PARTIAL PIPE FLOW; H = HEAD WATER CONTROL; 
C*       T = TAILWATER CONTROL; ? = CHECK "ERRORF"

C*                   A. FAN  10/8/82
C*************************************************************

C     Local variables

      LOGICAL WEIRCTRL
 
      REAL KE, QW, QC, WE, GGAP, HW, TW, ELMAX

C***----- INITIALIZATION FOR EACH NEW COMPUTATION -----
 
      QW = 0.0
      QC = 0.0
      QA = 0.0
 
      CODE(1:3) = '   '
      IF(GTYPE.EQ.2) THEN
         WE = BOARD
         GGAP = D 
      ELSE
         WE = 0.0 
         GGAP = MIN(BOARD,D)
      ENDIF 
 
 
C *****
      HW = HWE - INEL
      TW  =  TWE - OUTEL
      ELMAX = MAX( INEL, OUTEL ) 
      IF(GGAP.LE.0.) GO TO 110
      IF(GTYPE.EQ.2 .AND. HWE.LE.WE .AND. TWE.LE.WE) GO TO 110
      IF (ABS(HWE - TWE).LE.0.0001) GO TO 110 
      IF(HWE.LE.ELMAX .AND. TWE.LE.ELMAX) GO TO 110
 
C --- Set weir control flag. (JMO 1/4/92)
 
      WEIRCTRL = .FALSE.
 
      IF(GTYPE.EQ.2) THEN
         CALL QCULWEI(QW, CODE,
     1           HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL, L,
     +           D, W, N, K, C,
     2           HW, TW, KE, 
     3           WB, WE, SWB, SWE, CW, A, AW)
 
 
C ------ Flow through culvert with flashboard is weir flow if flashboard control
C        is significant; or if flow is free, not submerged.  (JMO 1/4/92)

C         IF (AW/A .LT. 0.2 .OR. TWE .LT. BOARD) THEN
C            WEIRCTRL = .TRUE.
C            GO TO 100
C         ENDIF
 
      ENDIF 
      IF(TW.GE.D) THEN 
         CALL QCULPIP(QC, CODE,
     1           HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL, L,
     +           D, W, N, K, C,
     2           HW, TW, KE) 
      ELSEIF(HW.GE.1.3*D .OR. HW.GE.2.0*GGAP) THEN
              CALL QCULORI(QC, CODE, 
     1                HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL,
     +                L, D,W, N,
     +                K, C,
     2                HW, TW, KE)
           ELSE 
              CALL QCULDIT(QC, CODE, 
     1                HWE, TWE, GGAP, BARREL, GTYPE, INEL, OUTEL,
     +                L, D, W, N,
     +                K, C, 
     2                HW, TW, KE)
      ENDIF 
 
C --- Flow through culvert is weir flow if flow is controlled by flashboard
C     either completely or primarily.  (JMO 1/4/92)
 
C  100 IF (WEIRCTRL .OR. (GTYPE .EQ. 2 .AND. QW. LE. QC) ) THEN
 
      IF(QW.LE.QC) THEN
C       Assume weir control.
        QA = QW
        CODE(1:1) = 'W'
        IF(GTYPE.EQ.2) THEN
          IF(TWE.LE.BOARD) THEN
C           Make sure that the flow control state is free.
            CODE(3:3) = 'F'
          ELSE
            CODE(3:3) = 'S'
          ENDIF
        ENDIF
      ELSE
         QA = QC
      ENDIF 
  110 RETURN
      END
C     ***********
C     *         *
C     * RES_SFWMD_FIND_FREE
C     *         *
C     ***********
      REAL FUNCTION RES_SFWMD_FIND_FREE(Z)

C     Residual function for finding the free flow values for 
C     the SFWMD gated culvert routine.  We will use a bisection
C     method.  This routine then returns two distinct values. 
C     +1 if the flow is called free and -1 if the flow is called
C     submerged. 

      REAL Z

      INCLUDE 'sfwmdqcv.cmn'
50    FORMAT(' CHK: HWE=',F8.2,' TWE=',F8.2,' FLOW=',F8.2,' CODE=',
     A       A3)
C**********************************************************************

      CALL SFWMD_QCULV(CODE, GTYPE, BARREL, BOARD, C, CW, D, HWE,
     A                  INEL, K, L, N, OUTEL, SWB, SWE, Z,
     B                  W, WB, QA)
      IF(CODE(3:3).EQ.'F') THEN
        RES_SFWMD_FIND_FREE = 1
      ELSE
        RES_SFWMD_FIND_FREE = -1
      ENDIF
      RETURN
      END
C     ***********
C     *         *
C     * BISECT
C     *         *
C     ***********

      SUBROUTINE BISECT(EPS, FBASE, ZBASE, 
     M                  ZL, FL, ZR, FR, F, 
     O                  ZROOT, EFLAG)

C     Find a root of a function by bisection.  EFLAG returned as 0
C     if convergence.  Otherwise returned as > 0.

      INTEGER EFLAG

      REAL EPS, FBASE, ZBASE, ZL, FL, ZR, FR, ZROOT

      REAL F
      EXTERNAL F

C     Local

      INTEGER KNT
      REAL ZM, FM
C***********************************************************************
C     On entry FL and FR must be of differing sign and ZL and ZR must
C     be distinct.  The function, F, may never become small.  Therefore,
C     we must use convergence testing that checks for smallness of
C     the interval as well as the function.
      
      EFLAG = 0
      KNT = 0
100   CONTINUE
        ZM = 0.5*(ZL + ZR)        
        FM = F(ZM)
        IF(FM*FL.LT.0.0) THEN
C         Retain the the sign change in the interval.
          FR = FM
          ZR = ZM
        ELSE
          FL = FM
          ZL = ZM
        ENDIF

        IF(ABS(FM/FBASE).LT.EPS.OR.ABS(ZL - ZR)/ZBASE.LT.EPS) THEN
C         Return the L value in order to get free flow case

          ZROOT = ZL
C         Evaluate the residual once more to make sure we
C         have the flow that agrees
          FL = F(ZL)

          RETURN
        ENDIF
        KNT = KNT + 1
        IF(KNT.GT.25) THEN
          EFLAG = 1
          ZROOT = ZM
          RETURN
        ENDIF
        GOTO 100
      END        
C     ***********
C     *         *
C     * SFWMD_FIND_FREE
C     *         *
C     ***********
      
      SUBROUTINE SFWMD_FIND_FREE(STDOUT, ZUP, 
     O                               ZDN, QFREE)

C     Find the free flow values for a SFWMD gated culvert. 

      INTEGER STDOUT

      REAL ZUP, ZDN, QFREE

      INCLUDE 'sfwmdqcv.cmn'

C     Local

      INTEGER EFLAG

      REAL  EPS, FBASE, ZBASE, FL, FR, ZL, ZR, DROP

      REAL RES_SFWMD_FIND_FREE
      EXTERNAL RES_SFWMD_FIND_FREE

C***********************************************************************
      HWE = ZUP

C     Find a sign change in the residual function and then use 
C     bisection to find the boundary between free and submerged flow.
C     The flow will always be free if the drop from the headwater
C     level to the tailwater level is small enough.  Put the tailwater
C     at the exit level from the culvert barrel.

      ZL =  OUTEL
      FL = RES_SFWMD_FIND_FREE(ZL)

C     Culvert will probably be submerged if the tailwater is 
C     within .005 of the drop from the high point of the culvert 
C     invert.

      DROP = 0.005*(ZUP -  MAX(INEL, OUTEL))
      ZR = ZUP - DROP
      FR = RES_SFWMD_FIND_FREE(ZR)

      IF(FR*FL.GE.0.0) THEN
        WRITE(STDOUT,*) ' Problem in sign change in SFWMD_FIND_FREE.'
        STOP 'Abnormal stop. Errors found.'
      ENDIF

      FBASE = 1
      ZBASE =  ZUP - MAX(INEL, OUTEL, BOARD)
      EPS = 0.001
      CALL BISECT(EPS, FBASE, ZBASE, 
     M                  ZL, FL, ZR, FR, RES_SFWMD_FIND_FREE, 
     O                  ZROOT, EFLAG)

      IF(EFLAG.GT.0) THEN
        WRITE(STDOUT,*) ' No convergence in SFWMD_FIND_FREE.'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      ZDN = ZROOT
      QFREE = QA
      RETURN
      END
C     ***********
C     *         *
C     * RISERCLV
C     *         *
C     ***********
      SUBROUTINE   RISERCLV(GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR)
 
C     Compute a 2-D table of type 13 for one or more SFWMD riser
C     culverts. Uses SFWMD methods for culvert flow.
 
C     Dummy arguments
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL GRAV
 
C     DEFINITIONS 
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TABDIR - Table directory to remember table numbers
 
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'sfwmdqcv.cmn'
 
C     Local
      INTEGER  I, IHUP, IPFD, J,JBASE, NFRAC, NHUP, N_GLOBAL, TAB,
     A         TABLT, TABGT, NN

      REAL  BIGERR, DROP,  FDROP, FDVEC(PMXNHU), HEAD_DATUM,
     A     HUPVEC(PMXNHU), LIMPFD, MAXZUP, LIPREC,  MAXHUP, MINHUP,
     B     MINPFD, PFDTMP(PMXFRC), PFDVEC(PMXFRC), POW,  QFREE, QHAT, 
     C     QMAT(PMXNHU,PMXFRC), RERR, RMS_GLOBAL, WORK(PMXFRC), 
     D     XBRK(PMXFRC), TWE, TWEF, FNUMBER, HTEMP, zrhufd

      real*8 easting, northing
      CHARACTER LABEL*50, LINE*80, SHAPE*6, TABID*16, KEY*16,
     a     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     Intrinsics
      INTRINSIC ABS, FLOAT, LOG, SQRT
 
 
C     External names
      INTEGER LENSTR
      EXTERNAL CHKTAB, inline, LSTOPF, TABCHK, TWDOUT,
     A         READ_TABID, LENSTR
C     ************************************FORMATS**********************
 1    FORMAT(7X,I5)
 2    FORMAT(6X,A)
4     FORMAT(I6,A6,7F6.0)
6     FORMAT(10F6.0)
 
 50   FORMAT(/,' Table id= ',A,' for type 13 table for riser',
     A     ' culvert.')
 52   FORMAT(/,' Label=',A50)
54    FORMAT(' ',A80)
56    FORMAT(1X,I6,A6,F6.2,2F6.1,2F6.1,F6.1,F6.3)
58    FORMAT(1X,F6.1,F6.1,F6.1,F6.1,F6.1,F6.1,F6.2,F6.2,F6.3,F6.3)
 60   FORMAT(/,' Datum for defining heads=',F10.2)
 72   FORMAT(/,'Upstream head=',F9.4,' Elevation=',F10.4)
 74   FORMAT(/,
     A '  Partial  Drop    Elev.   Head    Flow Discharge Local',/,
     C '   free    sect.   sect.   sect.   Code           power',/,
     E '   drop    1->4     4       4',/,
     F ' --------  ------  ------  ------   --- --------- ------')
 75   FORMAT(1X,F8.4,F8.3,F8.3,F8.3,3X,A3,F10.1,F7.2)
78    FORMAT(/,'  *ERR:634* Maximum ups head=',F8.2,' <= 0.',
     A ' Max ups elev=',F8.2,' and head datum=',F8.2)
 80   FORMAT(/,' *WRN:595* Minimum non-zero upstream head=',F8.2,
     A     ' <= 0.',' Setting to 0.15')
 86   FORMAT(/,' Maximum relative error=',F6.3,
     A   '  Upstream head=',F9.4,/,'  and partial free drop=', F8.5)
 88   FORMAT(/' Root-mean-squared error=',F6.3,' N in sample=',I5)
 89   FORMAT('  Processing RISERCLV TabId= ',A)
 99   FORMAT(/,' *ERR:635* Interpolation precision tables missing.',
     A   ' Check for version',/,11X,'of file: TYPE5.TAB.')
C**********************************************************************
C     Define the linear interpolation precision tables.
      KEY = '10001'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, KEY,
     M                             EFLAG,
     O                             TABLT)
      TABLT = FTPNT(TABLT)
      KEY = '10002'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, KEY,
     M                             EFLAG,
     O                             TABGT)
      TABGT = FTPNT(TABGT)
      IF(TABLT.LT.1.OR.TABGT.LT.1) THEN
        WRITE(STDOUT,99) 
        EFLAG = 1
        RETURN
      ENDIF
      CALL inline(STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'TAB',
     O                EFLAG, TABID, TAB)
C      READ(LINE,1,ERR=991) TAB
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID))
      WRITE(*,89) TABID(1:LENSTR(TABID))
 
      CALL TABCHK(STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)

c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, easting, basis,
     O             northing)

 
      CALL inline(STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) LABEL
      WRITE(STDOUT,52) LABEL
 
C     Get the heading line.
      CALL inline(STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,54) LINE

      CALL inline(STDIN, STDOUT,
     O           LINE)

      READ(LINE,4,ERR=991) NUMBER, SHAPE, K, D, W, INEL, OUTEL, L, N
      WRITE(STDOUT, 56) NUMBER, SHAPE, K, D, W, INEL, OUTEL, L, N
      CALL STRIP_L_BLANKS(SHAPE)
      IF(SHAPE.EQ.'PIPE') THEN
        W = D
      ENDIF

      FNUMBER = NUMBER
C     Get the heading line.
      CALL inline(STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,54) LINE

      CALL inline(STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CW, BOARD, WB, SWE, SWB, MAXZUP, MINHUP,
     A       LIMPFD, MINPFD, LIPREC
      WRITE(STDOUT,58)  CW, BOARD, WB, SWE, SWB, MAXZUP, MINHUP,
     A       LIMPFD, MINPFD, LIPREC

C     Gate option not used.  Set value of contraction coefficient
C     to avoid side effects. 
      C = 0.6
      IF(EFLAG.NE.0) RETURN
 
C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.

      IF(SHAPE.EQ.'PIPE') THEN
        BARREL = 0
      ELSE
        BARREL = 1
      ENDIF

C     Force to riser culvert only.  Gated culverts done elsewhere.
      GTYPE = 2

C     Compute the upstream head sequence.  First, find the head
C     range requested. 
      HEAD_DATUM = MAX(BOARD, INEL, OUTEL)
      WRITE(STDOUT,60) HEAD_DATUM
      MAXHUP = MAXZUP - HEAD_DATUM
      IF(MAXHUP.LE.0.0) THEN
        WRITE(STDOUT,78) MAXHUP, MAXZUP, HEAD_DATUM
        MAXHUP = 2.0
        EFLAG = 1
      ENDIF
      IF(MINHUP.LE.0.0) THEN
        WRITE(STDOUT,80) MINHUP
        MINHUP = 0.15
      ENDIF      
C     Use power of 1.5 because is requires closer spacing than 
C     power of .5. 
      POW = 1.5
      OFFSET = 0.0
      HTEMP = SWE - HEAD_DATUM
      IF(HTEMP + MINHUP.LT.MAXHUP) THEN
C       Define ups head sequence in two parts: MINHUP to HTEMP,
C       and then from HTEMP+MINHUP to MAXHUP
        CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, MINHUP, HTEMP,
     I              LIPREC, PMXNHU,
     O              NHUP, XBRK, EFLAG)

        DO 190 I=1,NHUP
          HUPVEC(I) = XBRK(I)
190     CONTINUE
 
        CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, HTEMP+MINHUP,
     I              MAXHUP, LIPREC, PMXNHU,
     O              NN, XBRK, EFLAG)
        
        DO 191 I=1,NN
          NHUP = NHUP + 1
          HUPVEC(NHUP) = XBRK(I)
191     CONTINUE
      ELSE
C       Define head sequence for whole range: MINHUP to MAXHUP.
        CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, MINHUP, MAXHUP,
     I              LIPREC, PMXNHU,
     O              NHUP, HUPVEC, EFLAG)
      ENDIF

C     Compute the proportions of free drop. 
      POW = 0.5
      OFFSET = 0.0
      CALL LSTOPF(STDOUT, TABLT, TABGT, POW, OFFSET, MINPFD, LIMPFD,
     I            LIPREC, PMXFRC,
     O            NN, XBRK, EFLAG)
 
      WORK(1) = 0.0
      DO 195 I=1,NN
        WORK(I+1) = XBRK(I)
 195  CONTINUE
      NFRAC = NN + 1

      WORK(NFRAC+1) = 1.0
      NFRAC = NFRAC + 1
 
C     Transfer to PFDVEC and insert the intermediate points for computing
C     an estimated interpolation error. 
      PFDVEC(1) = WORK(1)
      PFDVEC(2) = WORK(2)
      J = 2
      DO 201 I=3,NFRAC-1
        J = J + 1
        PFDVEC(J) = 0.5*(WORK(I) + WORK(I-1))
        J = J + 1
        PFDVEC(J) = WORK(I)
 201  CONTINUE
      PFDVEC(J+1) = WORK(NFRAC)
      NFRAC = J + 1
C      DO 202 I=1,NFRAC
C        WRITE(STDOUT,*) ' I=',I,' PFDVEC(I)=',PFDVEC(I)
C202   CONTINUE

      RMS_GLOBAL = 0.0
      N_GLOBAL = 0.0
      BIGERR = 0.0

      DO 1000 IHUP= 1, NHUP
        HW = HUPVEC(IHUP)

        HWE = HW + HEAD_DATUM

C       Find the free flow values.
        CALL SFWMD_FIND_FREE(STDOUT, HWE,
     O                             TWEF, QFREE)
        FDROP = HWE - TWEF

        QFREE = QFREE*FNUMBER
        QMAT(IHUP,NFRAC) = QFREE
C       Set flow to zero at this ups head for zero partial free drop
        QMAT(IHUP,1) = 0.0

        FDVEC(IHUP) = FDROP

        WRITE(STDOUT,72) HW, HWE
        WRITE(STDOUT,74)        
        WRITE(STDOUT,75) PFDVEC(NFRAC), FDROP,
     A                    TWEF, TWEF - HEAD_DATUM, CODE, QFREE
        
        DO 400 J= NFRAC-1, 2, -1
          PFD = PFDVEC(J)
          DROP = FDROP*PFD
          TWE = HWE - DROP
          CALL SFWMD_QCULV(CODE, GTYPE, BARREL, BOARD, C, CW, D, HWE,
     A                  INEL, K, L, N, OUTEL, SWB, SWE, TWE,
     B                  W, WB, QA)
          QA = QA*FNUMBER
          QMAT(IHUP,J) = QA
          POW = LOG(QA/QMAT(IHUP,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
          WRITE(STDOUT,75) PFDVEC(J), DROP,
     A                    TWE, TWE - HEAD_DATUM, CODE, QA, POW
400     CONTINUE

C       Compute approximate maximum error and report
        DO 600 J=3,NFRAC-2,2
          QHAT = 0.5*(QMAT(IHUP,J-1) + QMAT(IHUP,J+1))
          RERR = ABS(QHAT - QMAT(IHUP,J))/QMAT(IHUP,J)
          RMS_GLOBAL = RMS_GLOBAL + RERR*RERR
          N_GLOBAL = N_GLOBAL + 1
          IF(RERR.GT.BIGERR) THEN
            BIGERR = RERR
            HERR = HW
            IPFD = J
          ENDIF
 600    CONTINUE
 
C       Eliminate the checking values from QMAT
        JBASE = 3
        DO 700 J=4,NFRAC-1,2
          QMAT(IHUP,JBASE) = QMAT(IHUP,J)
          JBASE = JBASE + 1
 700    CONTINUE
        QMAT(IHUP,JBASE) = QMAT(IHUP,NFRAC)
 
1000  CONTINUE

 900    CONTINUE
        PFDTMP(1) = PFDVEC(1)
        PFDTMP(2) = PFDVEC(2)
C       Eliminate checking values of PFD
        JBASE = 3
        DO 910 J=4,NFRAC-1,2
          PFDTMP(JBASE) = PFDVEC(J)
          JBASE = JBASE + 1
 910    CONTINUE
        PFDTMP(JBASE) = PFDVEC(NFRAC)
 
        zrhufd = 0.0
        CALL TWDOUT(STDOUT, STDTAB, TAB, LABEL, NHUP, JBASE,
     I              HUPVEC, FDVEC, PFDTMP, QMAT, HEAD_DATUM,
     I              13, ' RISERCLV', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)
 
      WRITE(STDOUT,86) BIGERR, HERR,
     A                 PFDVEC(IPFD)
 
      RMS_GLOBAL = SQRT(RMS_GLOBAL/FLOAT(N_GLOBAL))
      WRITE(STDOUT,88) RMS_GLOBAL, N_GLOBAL
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
