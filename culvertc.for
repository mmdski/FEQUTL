C
C
C
      REAL FUNCTION   SBER
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute the subcritical residual for the energy based
C     steady flow profile computation
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'grvcom.cmn'
      INCLUDE 'sberc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL A, ALP, BET, DALP, DBET, DK, DT, E, EDDY, J, K, T, YLOCAL,
     A     YT
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
      YLOCAL = Y
      YT = MIN(YLOCAL, D)
      CALL XLKTAL
     I           (ADRS,
     M            YT,
     O            A, T, DT, J, K, DK, BET, DBET, ALP, DALP)
 
C      WRITE(STD6,*) ' SBER: Y=',YLOCAL,' K=',K,' KR=',KR,' DX=',DX
      IF(A.LT.AR) THEN
C       Flow is decelerating. Area on left is less than on the right.
        EDDY = KDEC*(1.0/A**2 - 1.0/AR**2)*QT**2/GRAV2
      ELSE
C       Flow is accelerating.  Area on left is larger than on the right.
        EDDY = KACC*(1.0/AR**2 - 1.0/A**2)*QT**2/GRAV2
      ENDIF
 
      E = Y + ALP*(QT/A)**2/GRAV2
      SBER = (E  - (DX*(QT*QT/(K*KR) + SE) + RHS + EDDY))/E
 
C      WRITE(STD6,*) ' SBER=',SBER,' SFM + SE=',QT*QT/(K*KR) + SE
      RETURN
      END
C
C
C
      SUBROUTINE   FULBAR
     I                   (STDOUT, A1TRUE, ALP1T, K1TRUE, Z1TRUE, A2FULL,
     I                    K2FULL, A3FULL, CDIS, AVH, Z3P,
     O                    Q, ZAT2)
 
C     + + + PURPOSE + + +
C     Compute the flow in the culvert given the elevation at
C     section 1, Z1TRUE, and the piezometric level at section 3,
C     Z3P.  The barrel is flowing full all the way.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      REAL A1TRUE, A2FULL, A3FULL, ALP1T, AVH, CDIS, K1TRUE, K2FULL, Q,
     A     Z1TRUE, Z3P, ZAT2
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     A1TRUE - area at section 1 for the known water surface elevation
C     ALP1T  - value of energy flux correction coefficient at section 1
C     K1TRUE - conveyance at section 1 for the known elevation there
C     Z1TRUE - known elevation at section 1
C     A2FULL - area at section 2 with full flow in barrel
C     K2FULL - full barrel conveyance at section 2
C     A3FULL - area at section 3 with full flow in barrel
C     CDIS   - discharge coefficient
C     AVH    - area for computing velocity head
C     Z3P    - elevation of piezometric surface at section 3
C     Q      - Flowrate
C     ZAT2   - elevation of water surface at section 2
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'appcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL CCON, DIV, NUM, NUM1, QT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* No convergence on Q in subroutine FULBAR.')
C***********************************************************************
C     Make the first estimate as if roadflow were zero.  If it is
C     zero, the result is valid.  Otherwise, iterate
C     to include the effect of velocity head of approach from the
C     flow over the road.
 
 
      NUM = GRAV2*(Z1TRUE - Z3P)
      DIV = (1.0 + CDIS**2*((AVH/A3FULL)**2 - 1.0 - (ALP1T - APPLOS)*
     A        (AVH/A1TRUE)**2 + GRAV2*AVH**2*(APPLEN/(K1TRUE*K2FULL)
     B        + FRCFAC)))
C      IF(FRCFAC.EQ.0.0) THEN
CC       Prismatic barrel.
C        DIV = (1.0 + GRAV2*(CDIS*A2FULL)**2*(APPLEN/(K1TRUE*K2FULL)
C     B       + L23/K2FULL**2 - (ALP1T - APPLOS)/(GRAV2*A1TRUE**2)))
C     A
C      ELSE
C        DIV = (1.0  + CDIS**2*((A2FULL/A3FULL)**2 - 1.0 +
C     B      GRAV2*A2FULL**2*(APPLEN/(K1TRUE*K2FULL)
C     C       + FRCFAC - (ALP1T - APPLOS)/(GRAV2*A1TRUE**2))))
C      ENDIF
C      WRITE(STDOUT,*) ' DIV=',DIV,' NUM=',NUM
      Q = CDIS*AVH*SQRT(NUM/DIV)
 
      IF(WFRD.GT.0.0) THEN
C       Add terms to the numerator to represent values affected
C       by flow over the road.
 
        KNT = 0
 100    CONTINUE
          IF(KNT.GT.100) THEN
            WRITE(STDOUT,50)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          NUM1 = NUM + WFRD*((ALP1T - APPLOS)*
     A                         (WFRD + Q + Q)/(GRAV2*A1TRUE**2)
     B                          - APPLEN*Q/(K1TRUE*K2FULL))
          QT = CDIS*AVH*SQRT(NUM1/DIV)
 
C          WRITE(STDOUT,*) ' KNT=',KNT,' Q=',Q,' QT=',QT
          IF(ABS(QT - Q)/Q.GT.EPSF) THEN
            Q = QT
            GOTO 100
          ENDIF
      ENDIF
 
C     Compute estimate of the piezometric level at section 2.   Estimate
C     the coefficient of contraction from the discharge coefficient.
 
      CCON = 1.0/(SQRT(1.0/CDIS**2 -1.0) + 1.0)
      ZAT2 = ALP1T*((Q + WFRD)/A1TRUE)**2/GRAV2 + Z1TRUE -
     A                (Q/(CCON*A2FULL))**2/GRAV2
C      WRITE(STDOUT,*) ' ZAT2=',ZAT2,' CCON=',CCON
      RETURN
      END
C
C
C
      SUBROUTINE   SFPCC
     I                  (IU, ID, DH, Q, ZDN,
     O                   IS, ZUP, SFLAG)
 
C     + + + PURPOSE + + +
C     Compute steady flow hydraulic gradeline for a closed
C     conduit flowing full. ZDN gives the elevation of the
C     hydraulic grade line at the downstream end of the
C     conduit.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IS, IU, SFLAG
      REAL DH, Q, ZDN, ZUP
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     DH     - Entrance loss for the culvert
C     Q      - Flowrate
C     ZDN    - Water surface elevation at downstream end of barrel
C     IS     - Index of node for the last value computed
C     ZUP    - Water surface elevation at section 2
C     SFLAG  - Flag describing results:
C                0-results incomplete.  Super critical flow encountered
C                  and flow is thought to be caused by a steep slope.
C                1-results complete
C               -1-results incomplete.  Super critical flow encountered
C                  and DX may be too large.
C               -2-results incomplete.  Super critical flow encountered
C                  and DX is too large.
C               -3-no solution found at the last point.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'culcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADRS, I
      REAL AL, ALPL, ALPR, AR, BETL, BETR, DALPL, DALPR, DBETL, DBETR,
     A     DKL, DKR, DTL, DTR, DX, JL, JR, KL, KR, QL, QR, SE, TL, TR,
     B     YA, YL, YR, ZBL, ZBR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
C      WRITE(STDOUT,60)
      QR = Q
      QL = Q
      ZBR = ZBVEC(ID)
      YR = ZDN - ZBR
C     DH GIVES THE LOSS OF ENERGY HEAD CAUSED BY THE CULVERT IN
C     EXCESS OF BARREL FRICTION. DISTRIBUTE ALONG LENGTH OF THE
C     BARREL
 
      SE = DH/ABS(XVEC(IU) - XVEC(ID))
C      WRITE(STDOUT,*) ' SFPCC: DH=',DH,' SE=',SE,' Q=',Q
C      ZR = ZDN
      ADRS = NSEC(ID)
      YA = DVEC(ID)
      CALL XLKTAL
     I           (ADRS,
     M            YA,
     O            AR, TR, DTR, JR, KR, DKR, BETR, DBETR, ALPR, DALPR)
 
C      WRITE(STDOUT,62) ID, XVEC(ID), YR, ZR, DVEC(ID)
 
      YVECSB(ID) = YR
      DO 500 I=ID-1,IU,-1
        ADRS = NSEC(I)
        DX = ABS(XVEC(I+1) - XVEC(I))
        ZBL = ZBVEC(I)
        YA = DVEC(I)
        CALL XLKTAL
     I             (ADRS,
     M              YA,
     O              AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL, DALPL)
 
        YL = YR + ZBR + ALPR*(QR/AR)**2/GRAV2 + DX*(QL*QR/(KL*KR) + SE)
     A      - ZBL - ALPL*(QL/AL)**2/GRAV2
C        ZL = ZBL + YL
        ZBR = ZBL
        YR = YL
        QR = QL
        KR = KL
        AR = AL
        ALPR = ALPL
        YVECSB(I) = YL
C        WRITE(STDOUT,62) I, XVEC(I), YL, ZL, DVEC(I)
 500  CONTINUE
 
      ZUP = YL + ZBL
      IS = IU
      SFLAG = 1
      RETURN
      END
C
C
C
      REAL FUNCTION   SPER
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute the super critical residual for the energy based
C     steady flow profile computation
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'grvcom.cmn'
      INCLUDE 'sperc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL A, ALP, BET, DALP, DBET, DK, DT, E, EDDY, J, K, T, YLOC, YT
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKTAL
C***********************************************************************
      YLOC = Y
      YT = MIN(YLOC, D)
      CALL XLKTAL
     I           (ADR,
     M            YT,
     O            A, T, DT, J, K, DK, BET, DBET, ALP, DALP)
 
      IF(AL.LT.A) THEN
C       Flow is decelerating. Area on left is less than on the right.
        EDDY = KDEC*(1.0/AL**2 - 1.0/A**2)*QR**2/GRAV2
      ELSE
C       Flow is accelerating.  Area on left is larger than on the right.
        EDDY = KACC*(1.0/A**2 - 1.0/AL**2)*QR**2/GRAV2
      ENDIF
      E = Y + ALP*(QR/A)**2/GRAV2
      SPER = (RHS - (E + DX*(QL*QR/(K*KL) + SE) + EDDY))/E
 
      RETURN
      END
C
C
C
      SUBROUTINE   SBFEBC
     I                   (STDOUT, EU, ED, PLCWTB, GLCWTB, PHCWTB,
     I                    GHCWTB, PSUBTB, GSUBTB, NOFF, SURF, HLCRIT,
     I                    XL, XR, HL, HM, HR, QL, QM, QR, TOTHL, TOTHM,
     I                    TOTHR, YFL, YFM, YFR, APPL, APPM, APPR, WL,
     I                    WM, WR, AELL, AELM, AELR, RMFFAC,
     O                    QROAD, MROAD, EROAD)
 
C     + + + PURPOSE + + +
C     SuBmerged Flow EmBankment with Culvert - Find the submerged flow,
C     estimated momentum flux, and estimated energy flux
C     over a submerged embankment.  The
C     free flow computations establish most of the argument values
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER GHCWTB, GLCWTB, GSUBTB, NOFF, PHCWTB, PLCWTB, PSUBTB,
     A        STDOUT
      REAL AELL(*), AELM(*), AELR(*), APPL(*), APPM(*), APPR(*), ED,
     A     EROAD, EU, HL(*), HLCRIT, HM(*), HR(*), MROAD, QL(*), QM(*),
     B     QR(*), QROAD, RMFFAC, TOTHL(*), TOTHM(*), TOTHR(*), WL(*),
     C     WM(*), WR(*), XL(*), XR(*), YFL(*), YFM(*), YFR(*)
      CHARACTER SURF(*)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     EU     - Water surface elevation upstream.
C     ED     - Downstream water surface elevation
C     PLCWTB - Paved low-head weir coefficient table
C     GLCWTB - Gravel surface low head weir coefficient table
C     PHCWTB - Paved high-head weir coefficient table
C     GHCWTB - Gravel surface high head weir coefficient table
C     PSUBTB - Paved submergence table
C     GSUBTB - Gravel surface submergence table
C     NOFF   - Number of offsets
C     SURF   - Nature of the embankment surface
C     HLCRIT - Ratio of piezometric head to crest breadth at boundary
C              between low head and high head flow
C     XL     - Offset at left hand end of segment
C     XR     - Offset at right hand end of segment
C     HL     - Piezometric head on left end of segment
C     HM     - Piezometric head on middle of segment
C     HR     - Piezometric head on right end of segment
C     QL     - Flow on left hand end of segment
C     QM     - Flow at middle of the segment
C     QR     - Flow on right hand end of segment
C     TOTHL  - Total head at left hand end of segment
C     TOTHM  - Total head at middle of segment
C     TOTHR  - Total head at right hand end of segment
C     YFL    - Estimated depth at crest at left hand end of segment
C     YFM    - Estimated depth at crest at middle of segment
C     YFR    - Estimated depth at crest at right hand end of segment
C     APPL   - Elevation of approach at left end of line segment
C     APPM   - Elevation of approach at middle of line segment
C     APPR   - Elevation of approach at right end of line segment
C     WL     - Breadth of the crest at left hand end of segment
C     WM     - Breadth of the crest at middle of segment
C     WR     - Breadth of the crest at right hand end of segment
C     AELL   - Estimated ratio of head loss to head loss at incipient
C              submergence at left end of line segment
C     AELM   - Estimated ratio of head loss to head loss at incipient
C              submergence at middle of line segment
C     AELR   - Estimated ratio of head loss to head loss at incipient
C              submergence at right end of line segment
C     RMFFAC - Adjustment factor for roadway momentum flux
C     QROAD  - Flow over the roadway
C     MROAD  - Momentum flux from flow over the road
C     EROAD  - Energy flux over the roadway
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'grvcom.cmn'
      INCLUDE 'subcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HCWTAB, ISEG, LCWTAB, SUBTAB
      REAL CREST, DEPTH, DX, EL, ELEFT, EM, EMID, ER, ERIGHT, ESEG, HD,
     A     MLEFT, MMID, MRIGHT, MSEG, QLEFT, QMID, QRIGHT, QSEG, SELL,
     B     SELM, SELR, SHTOT, SRAT, YC, YCRIT
      DOUBLE PRECISION ESUB, MSUB, QSUB
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL INVSE
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL INVSE, STOTHQ
C***********************************************************************
C     COMPUTE THE SUBMERGED FLOW GIVEN THAT THE FREE FLOW VALUES
C     HAVE BEEN ESTABLISHED AND GIVEN A DOWNSTREAM WATER SURFACE
C     ELEVATION.
 
      QSUB = 0.D0
      MSUB = 0.D0
      ESUB = 0.D0
 
C      WRITE(STDOUT,*) ' '
C      WRITE(STDOUT,*) ' SBFEBC'
 
      IF(ED.GT.EU) THEN
        WRITE(STDOUT,*) ' *BUG:XXX* DOWNSTREAM ELEV > UPSTREAM ELEV.'
        WRITE(STDOUT,'('' EU='',F10.4,'' ED='',F10.4)') EU, ED
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(ED.EQ.EU) THEN
        QROAD = 0.0
        MROAD = 0.0
        EROAD = 0.0
        RETURN
      ENDIF
 
C      WRITE(STDOUT,'('' EU='',F10.4,'' ED='',F10.4)') EU, ED
 
C     FOR EACH SEGMENT OF THE WEIR FIND THE EFFECT OF THIS
C     VALUE OF DOWNSTREAM PIEZOMETRIC HEAD
 
      DO 4100 ISEG=1,NOFF-1
 
C        WRITE(STDOUT,*) ' '
C        WRITE(STDOUT,*) ' ISEG=',ISEG
 
C       SELECT THE TABLES FOR SUBMERGENCE AND WEIR COEF
        IF(SURF(ISEG).EQ.'PAVED') THEN
          SUBTAB = PSUBTB
          HCWTAB = PHCWTB
          LCWTAB = PLCWTB
          SRAT = PRAT
        ELSE
          SUBTAB = GSUBTB
          HCWTAB = GHCWTB
          LCWTAB = GLCWTB
          SRAT = GRAT
        ENDIF
 
        IF(HL(ISEG).GT.0.0.OR.HR(ISEG).GT.0.0) THEN
C         SEGMENT HAS NON-ZERO FREE FLOW
          QLEFT = QL(ISEG)
          SHTOT = TOTHL(ISEG)
          IF(QLEFT.GT.0.0) THEN
C           NOTE THAT THE UPSTREAM WATER SURFACE ELEVATION
C           LESS THE PIEZOMETRIC HEAD GIVES THE WEIR CREST
C           ELEVATION
 
            CREST = EU - HL(ISEG)
            HD = ED - CREST
 
            IF(HD.GT.SRAT*SHTOT) THEN
C             FLOW SUBMERGED AT THIS POINT. FIND NEW VALUES
              DEPTH = EU - APPL(ISEG)
              CALL STOTHQ
     I                   (HCWTAB, LCWTAB, SUBTAB, HLCRIT, WL(ISEG),
     I                    HL(ISEG), HD, DEPTH,
     M                    SHTOT,
     O                    QLEFT)
 
C             COMPUTE ESTIMATE OF  SUBMERGED ENERGY LOSS
 
              SELL = SHTOT - (QLEFT/(HD + DEPTH - HL(ISEG)))**2/GRAV2
     A                 - HD
C             ESTIMATE CREST DEPTH ASSUMING SOME ENERGY LOSS
C             The over riding rule is that the crest depth used to 
C             estimate energy and momentum fluxes cannot be less than
C             the depth estimated for free flow.  We assume this is
C             true because submergence would affect both the flow and
C             the crest depth.  
              EL = SHTOT - AELL(ISEG)*SELL
              YC = INVSE(QLEFT, EL, HD)
              YC = MAX(HD, YC)
              YC = MAX(YC, YFL(ISEG))
            ELSE
              YC = YFL(ISEG)
            ENDIF
            MLEFT = QLEFT**2/YC
            ELEFT = QLEFT*((QLEFT/YC)**2/GRAV2 + CREST + YC)
          ELSE
            QLEFT = 0.0
            MLEFT = 0.0
            ELEFT = 0.0
          ENDIF
 
C         WRITE(STDOUT,*) ' LEFT: DEPTH=',DEPTH,' HD=',HD,' QLEFT=',
C     A               QLEFT,' SHTOT=',SHTOT, ' HEAD=',HL(ISEG),
C     B               ' SELL=', SELL,' YC=',YC
 
          QMID = QM(ISEG)
          SHTOT = TOTHM(ISEG)
          IF(QMID.GT.0.0) THEN
            CREST = EU - HM(ISEG)
            HD = ED - CREST
            IF(HD.GT.SRAT*SHTOT) THEN
              DEPTH = EU - APPM(ISEG)
              CALL STOTHQ
     I                   (HCWTAB, LCWTAB, SUBTAB, HLCRIT, WM(ISEG),
     I                    HM(ISEG), HD, DEPTH,
     M                    SHTOT,
     O                    QMID)
 
C             COMPUTE ESTIMATE OF  SUBMERGED ENERGY LOSS
 
              SELM = SHTOT - (QMID/(HD + DEPTH - HM(ISEG)))**2/GRAV2
     A                 - HD
C             ESTIMATE CREST DEPTH ASSUMING SOME ENERGY LOSS
              EM = SHTOT - AELM(ISEG)*SELM
              YC = INVSE(QMID, EM, HD)
              YC = MAX(HD, YC)
              YC = MAX(YC, YFM(ISEG))
            ELSE
              YC = YFM(ISEG)
            ENDIF
            MMID = QMID**2/YC
            EMID = QMID*((QMID/YC)**2/GRAV2 + CREST + YC)
          ELSE
            QMID = 0.0
            MMID = 0.0
            EMID = 0.0
          ENDIF
 
C         WRITE(STDOUT,*) ' MID: DEPTH=',DEPTH,' HD=',HD,' QMID=',
C     A               QMID,' SHTOT=',SHTOT, ' HEAD=',HM(ISEG),
C     B               ' SELM=', SELM  ,' YC=',YC
          QRIGHT = QR(ISEG)
          SHTOT = TOTHR(ISEG)
          IF(QRIGHT.GT.0.0) THEN
            CREST = EU - HR(ISEG)
            HD = ED - CREST
            IF(HD.GT.SRAT*SHTOT) THEN
              DEPTH = EU - APPR(ISEG)
              CALL STOTHQ
     I                   (HCWTAB, LCWTAB, SUBTAB, HLCRIT, WR(ISEG),
     I                    HR(ISEG), HD, DEPTH,
     M                    SHTOT,
     O                    QRIGHT)
 
C             COMPUTE ESTIMATE OF  SUBMERGED ENERGY LOSS
 
              SELR = SHTOT - (QRIGHT/(HD + DEPTH - HR(ISEG)))**2/GRAV2
     A                 - HD
C             ESTIMATE CREST DEPTH ASSUMING SOME ENERGY LOSS
              ER = SHTOT - AELR(ISEG)*SELR
              YC = INVSE(QRIGHT, ER, HD)
              YC = MAX(YC, HD)
              YC = MAX(YC, YFR(ISEG))
            ELSE
              YC = YFR(ISEG)
            ENDIF
            MRIGHT = QRIGHT**2/YC
            ERIGHT = QRIGHT*((QRIGHT/YC)**2/GRAV2 + CREST + YC)
          ELSE
            QRIGHT = 0.0
            MRIGHT = 0.0
            ERIGHT = 0.0
          ENDIF
 
C         WRITE(STDOUT,*) ' RIGHT: DEPTH=',DEPTH,' HD=',HD,' QRIGHT=',
C     A               QRIGHT,' SHTOT=',SHTOT, ' HEAD=',HR(ISEG),
C     B               ' SELR=',SELR,' YC=',YC
 
C         NOW COMPUTE THE FLOW AS AFFECTED BY SUBMERGENCE
 
          DX = ABS(XR(ISEG) - XL(ISEG))
 
C         INTEGRATE OVER THE WETTED LENGTH USING SIMPSON'S RULE
 
          QSEG = DX*(QLEFT + 4.*QMID + QRIGHT)/6.
          MSEG = DX*(MLEFT + 4.*MMID + MRIGHT)/6.
          ESEG = DX*(ELEFT + 4.*EMID + ERIGHT)/6.
 
C          WRITE(STDOUT,*) ' '
C          WRITE(STDOUT,*) ' DX =',DX
C          WRITE(STDOUT,'('' SEG FLOWS:'',3F12.4)') QLEFT, QMID,
C     A                                  QRIGHT
 
C          WRITE(STDOUT,'('' QSEG='',F12.4)')
C     A                    QSEG
 
          QSUB = QSUB + QSEG
          MSUB = MSUB + MSEG
          ESUB = ESUB + ESEG
        ENDIF
 4100 CONTINUE
 
C     STORE THE VALUE AT THIS FRACTION OF THE FREE DROP
 
      QROAD = QSUB
      MROAD = RMFFAC*MSUB
      EROAD = ESUB
C       WRITE(STDOUT,*) ' SBFEBC: QROAD=',QROAD,' MROAD=',MROAD,
C     A        ' EROAD=',EROAD
      RETURN
      END
C
C
C
      SUBROUTINE   SFPSBE
     I                   (STDOUT, IU, ID, DH, Q, ZDN,
     O                    IS, ZUP, SFLAG)
 
C     + + + PURPOSE + + +
C     Compute a sub-critial steady flow profile in the channel defined
C     by NSEC, XVEC, ZBVEC, IU, AND ID with the starting water surface
C     elevation at the downstream end(ID) given in ZDN.  The sequence
C     of depths is contained in YVECSB.  The elevation at the upstream end
C     of the channel is in ZUP and IS points to the last node at which
C     a subcritical flow was found.  SFLAG is used to indicate the
C     results:  SFLAG = 1: computations completed and ZUP and YVECSB
C     contain the final results.  SFLAG = 0: computations not completed
C     supercritical flow was encountered.  It is believed that the
C     supercritical flow is physically caused by a steep slope.
C     The last valid depth is contained in YVECSB(IS).  SFLAG = -1:
C     computations not completed, supercritical flow was encountered
C     and it is believed to be caused by DX being too large. The last
C     valid depth is again in YVECSB(IS).  SFLAG = -2: computations
C     not completed supercritical flow was encountered and it
C     is clear that DX was too large.  The last valid depth is
C     in YVECSB(IS). SFLAG = -3: no solution found at the last
C     point.  SFLAG = 2: initial condition is super critical.
C     DH gives estimate of the loss in head over the reach from
C     other sources.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IS, IU, SFLAG, STDOUT
      REAL DH, Q, ZDN, ZUP
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     DH     - Entrance loss for the culvert
C     Q      - Flowrate
C     ZDN    - Water surface elevation at downstream end of barrel
C     IS     - Index of node for the last value computed
C     ZUP    - Water surface elevation at section 2
C     SFLAG  - Flag describing results:
C                0-results incomplete.  Super critical flow encountered
C                  and flow is thought to be caused by a steep slope.
C                1-results complete
C               -1-results incomplete.  Super critical flow encountered
C                  and DX may be too large.
C               -2-results incomplete.  Super critical flow encountered
C                  and DX is too large.
C               -3-no solution found at the last point.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'sberc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADL, ADR, FLAG, I, NDSIDE
      REAL AL, ALPL, BETL, BETR, DALPL, DALPR, DBETL, DBETR, DKL, DKR,
     A     DL, DR, DTL, DTR, EL, ER, FL, FR, FROUDE, JL, JR, KL, QC,
     B     QCL, QCR, QN, SB, SBOLD, TL, TR, XD, XL, XR, YC, YL, YMAX,
     C     YNL, YNR, YR, YSUB, YT, ZBL, ZBR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL FMXARG, SBER
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FMXARG, FNDCDE, FNDND, GETTBN, LKTQC, REGFLT, SBER,
     A         XLKT22, XLKTAL, LENSTR, GET_TABID
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:597* INITIAL DEPTH=',F10.2,' <= 0 IN SFPSBE.')
 56   FORMAT(' *ERR:598* TABID=',A,' overflow seeking subcritical',
     A       ' solution. MaxArg=',F10.2)
 58   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' SFPSBE.')
 70   FORMAT(' *BUG:XXX* REGFLT TAKES MORE THAN 100 ITERATIONS:SFPSBE.')
C***********************************************************************
      IF(DH.GT.0.0) THEN
        SE = DH/ABS(XVEC(IU) - XVEC(ID))
      ELSEIF(DH.EQ.0.0) THEN
        SE = 0.0
      ENDIF
C      WRITE(STDOUT,*) ' SFPSBE: DH=',DH,' SE=',SE
C      WRITE(STDOUT,60)
 
      QT = Q
 
      ZBR = ZBVEC(ID)
      XD = XVEC(ID)
      XR = XVEC(ID)
      XL = XVEC(ID-1)
      YR = ZDN - ZBR
C      WRITE(STDOUT,*) ' SFPSBE: YR=',YR
      YVECSB(ID) = YR
C      ZR = ZDN
      IF(YR.LE.0.0) THEN
        WRITE(STDOUT, 50) YR
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      ADR = NSEC(ID)
      D = DVEC(ID)
      YT = MIN(YR, D)
      CALL XLKT22
     I           (ADR,
     M            YT,
     O            AR, TR, DTR, JR, KR, DKR, BETR, DBETR, ALPR, DALPR,
     O            QCR)
 
      AVEC(ID) = AR
      KVEC(ID) = KR
C     COMPUTE THE FROUDE NUMBER
 
      FR = Q/QCR
C      WRITE(STDOUT,62) ID, XVEC(ID), YR, ZR, DVEC(ID)
 
C      WRITE(STDOUT,*) ' SFPSBE: FR=',FR
      IF(FR.GT.1.025) THEN
        SFLAG = 2
C        WRITE(STDOUT, 54) FR
C        WRITE(STDOUT,*) ' Q=',Q,' AR=',AR,' YR=',YR,' TR=',TR,
C     A                  ' ALP4=',ALPR,' DALPR=',DALPR,' QCR=',QCR
C        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      DX = ABS(XD - XL)
      ZBL = ZBVEC(ID-1)
      SBOLD = (ZBL - ZBR)/DX
C      WRITE(STDOUT,*) ' SBOLD=',SBOLD
      IF(FR.GT.0.975.AND.FR.LT.1.025) THEN
C       COMPUTE NORMAL FLOW AT INITIAL DEPTH
        IF(SBOLD.GT.0.0) THEN
          QN = KR*SQRT(SBOLD)
        ELSE
          QN = 0.0
        ENDIF
C        WRITE(STDOUT,*) ' SFPSBE: AT START QCR=',QCR,' QN=',QN
        IF(QN.GT.QCR) THEN
C         INITIAL SLOPE IS SUPER CRITICAL.  NO SUB-CRITICAL SOLUTION
C         EXISTS WITH INITIAL CRITICAL DEPTH.
          WRITE(STDOUT,*) ' SFPSBE: Dns invert slope is super critical.'
c          SFLAG = 0   Change to 2 on 27 Feb. 2008.  
          SFLAG = 2
          IS = ID
C          WRITE(STDOUT,*) ' ERROR RETURN: IS=',IS,' SFLAG=',SFLAG
          RETURN
        ENDIF
      ENDIF
C     AT THIS POINT THE INITIAL CONDITION IS VALID. START THE LOOP
 
C     MAKE FIRST ESTIMATE OF YC AND YN
      YC = YR
      YNR= YR
 
C     COMPUTE NORMAL DEPTH FOR THE INITIAL POINT.  SBOLD IS ASSUMED
C     TO BE CONSTANT OVER DX.
 
C      WRITE(STDOUT,*) ' SBOLD=',SBOLD
      IF(SBOLD.GT.0.0) THEN
        CALL FNDND
     I            (STDOUT, ADR, Q, SBOLD,
     M             DVEC(ID),
     O             YNR)
      ELSE
        YNR = 0.0
      ENDIF
C      WRITE(STDOUT,*) ' YNR=',YNR
      YNVEC(ID) = YNR
      YNL = YNR
      ADL = ADR
 
C     FIND CRITICAL DEPTH AT START POINT
C      WRITE(STDOUT,*) ' ESTIMATED CRITICAL DEPTH=',YC,' Q=',Q
      CALL FNDCDE
     I           (STDOUT, ADR, Q,
     M            YC)
C      WRITE(STDOUT,*) ' SFPSBE: CRITICAL DEPTH AT START=',YC,' Q=',Q
 
      ER = YR + ALPR*(Q/AR)**2/GRAV2
C      WRITE(STDOUT,*) ' SFPSBE: TEL at exit=', ER +  ZBR
C      EVEC(ID) = ER
C      YCVEC(ID) = YC
C      SFVEC(ID) = (Q/KR)**2
      DR = DVEC(ID)
      DO 500 I=ID-1,IU,-1
        D = DVEC(I)
        DL = D
        IS = I
        XL = XVEC(I)
        DX = ABS(XR - XL)
        ZBL = ZBVEC(I)
        KACC = KA(I+1)
        KDEC = KD(I+1)
        SB = (ZBL - ZBR)/DX
C        WRITE(STDOUT,*) ' SB=',SB
 
        ADL = NSEC(I)
 
C       SET ADDRESS IN COMMON BLOCK FOR RESIDUAL COMPUTATIONS
        ADRS = ADL
 
        IF(SB.GT.0.0.AND.2.*ABS(SB - SBOLD)/(SB + SBOLD).GT.5.E-3) THEN
C         RECOMPUTE THE NORMAL DEPTH AT RIGHT END.
 
          CALL FNDND
     I              (STDOUT, ADR, Q, SB,
     M               DVEC(I+1),
     O               YNR)
C          WRITE(STDOUT,*) ' RECOMPUTED YNR=',YNR
        ELSE
C         USE NORMAL DEPTH FROM UPSTREAM END OF PREVIOUS ELEMENT
 
          YNR = YNL
        ENDIF
 
        IF(ADR.NE.ADL) THEN
          CALL FNDCDE
     I               (STDOUT, ADL, Q,
     M                YC)
C          WRITE(STDOUT,*) ' SFPSBE: CRITICAL DEPTH=', YC
          IF(SB.GT.0.0) THEN
            CALL FNDND
     I                (STDOUT, ADL, Q, SB,
     M                 DVEC(I),
     O                 YNL)
C            WRITE(STDOUT,*) ' RECOMPUTED YNL=',YNL
          ELSE
            YNL = 0.0
          ENDIF
        ELSE
          YNL = YNR
        ENDIF
C        WRITE(STDOUT,*) ' YNL=',YNL
C       SET THE SIDE OF NORMAL DEPTH FOR THE INITIAL DEPTH IN THIS
C       ELEMENT
 
        IF(YR.GE.YNR) THEN
          NDSIDE = 1
        ELSE
          NDSIDE = -1
        ENDIF
 
C        WRITE(STDOUT,*) ' NDSIDE=',NDSIDE
 
C       COMPUTE THE CONSTANT PART OF THE EQUATION
 
        RHS = ER + ZBR - ZBL
C        RHS = ER + ZBR - ZBL + DX*SE
 
C         WRITE(STDOUT,*) ' RHS=',RHS
 
        YMAX = FMXARG(ADL)
 
        YT = YC
        FL = SBER(YT)
C        WRITE(STDOUT,*) ' YT=',YT,' FL=',FL
        IF(FL.LE.0.0) THEN
C         SUBCRITICAL SOLUTION EXISTS.  SEARCH FOR A POSITIVE
C         RESIDUAL
 
          IF(YR.GT.YC) THEN
            YSUB = YR
          ELSE
            YSUB = 1.01*YC
          ENDIF
 
 110      CONTINUE
            FR = SBER(YSUB)
C            WRITE(STDOUT,*) ' YSUB=',YSUB,' FR=',FR
            IF(FR.LT.0.0) THEN
              FL = FR
              YT = YSUB
              YSUB = 1.05*YSUB
              IF(ABS(YMAX - YSUB).LE.EPSARG) THEN
                TABID = GET_TABID(GETTBN(ADL))
                WRITE(STDOUT,56) TABID(1:LENSTR(TABID)), YMAX
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 110
            ENDIF
        ELSE
C         NO SUBCRITICAL SOLUTION
          IF(ABS(FL).LT.EPSF) THEN
C           Take critical depth as the result.
            YL = YC
            CALL XLKT22
     I                 (ADL,
     M                  YL,
     O                  AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL,
     O                  DALPL, QCL)
            FROUDE = 1.0
            GOTO 200
          ENDIF
C         USE NORMAL DEPTH.  IF THE SUPERCRITICAL DEPTH IS SPURIOUS
C         NORMAL DEPTH SHOULD PASS THE TESTS FOR A VALID SOLUTION
C          WRITE(STDOUT,*) ' TAKING NORMAL DEPTH=',YNL
          IF(YNL.EQ.0.0) THEN
            WRITE(STDOUT,*) ' *ERR:698* No subcritical solution',
     A       ' in steady flow profile.'
            WRITE(STDOUT,*) ' Try reducing step length'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          YL = YNL
          YT = MIN(YL, D)
          CALL XLKT22
     I               (ADL,
     M                YT,
     O                AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL,
     O                DALPL, QCL)
C
          FROUDE = Q/QCL
C          WRITE(STDOUT,*) ' NORMAL DEPTH FROUDE=',FROUDE
          GOTO 200
        ENDIF
 
C       WE HAVE A SIGN CHANGE IN THE INTERVAL (YC, YSUB) OR
C       AT LEAST ONE RESIDUAL IS ZERO
 
C        WRITE(STDOUT,*) ' ROOT BRACKETED'
C        WRITE(STDOUT,*) ' YLOW=',YT,' FL=',FL
C        WRITE(STDOUT,*) ' YHIGH=',YSUB,' FH=',FR
 
        CALL REGFLT
     I             (0.0, 5.E-6, SBER,
     M              YT, YSUB, FL, FR,
     O              YL, FLAG)
        IF(FLAG.EQ.1) THEN
          WRITE(STDOUT,58)
          STOP 'Abnormal stop. Errors found.'
        ELSEIF(FLAG.EQ.2) THEN
          WRITE(STDOUT,70)
          STOP 'Abnormal stop. Errors found.'
        ENDIF
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        YT = MIN(YL, D)
        CALL XLKTAL
     I             (ADL,
     M              YT,
     O              AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL, DALPL)
 
C       CHECK IF THE SOLUTION IS VALID.  ALLOW SLIGHTLY SUPERCRITICAL
C       FLOW TO ACCOMMADATE THE DISCREPENCIES IN TABULATED CRITICAL
C       FLOWS.
 
        CALL LKTQC
     I            (ADL,
     M             YT,
     O             QC)
        FROUDE = Q/QC
C        WRITE(STDOUT,*) ' AT YL=',YL,' FROUDE=',FROUDE
 
 200    CONTINUE
        IF(FROUDE.GT.1.05) THEN
C         SOLUTION IS NO LONGER VALID.
C         NO SUBCRITICAL SOLUTION EXISTS.
C          WRITE(STDOUT,*) ' AT YL=',YL,' FROUDE=',FROUDE,' IS=',IS
 
          SB = (ZBL - ZBR)/DX
          IS = I + 1
C         Store current critical depth in current node's slot.
          YVECSB(I) = YC
          IF(SB.LE.0.0) THEN
C           SUPER CRITICAL FLOW IS CLEARLY A COMPUTATIONAL ARTIFACT
            SFLAG = -2
          ELSE
            QN = KL*SQRT(SB)
C            WRITE(STDOUT,*) ' QN AT LEFT=',QN
            IF(QN.GT.QC) THEN
C             SUPER CRITICAL FLOW MAY HAVE A PHYSICAL BASIS
              SFLAG = 0
            ELSE
C             CHANNEL SLOPE IS NOT STEEP. COMPUTATIONAL ARTIFACT.
              SFLAG = -2
            ENDIF
          ENDIF
C          WRITE(STDOUT,*) ' ERROR RETURN: IS=',IS,' SFLAG=',SFLAG
C          GOTO 501
          ZUP = YL + ZBL
          RETURN
        ENDIF
 
        IF(YNL.LT.DL.AND.YNR.LT.DR) THEN
C         CHECK IF ESTIMATED PROFILE HAS CROSSED NORMAL DEPTH
C         when flow is free surface at both ends.
          IF(NDSIDE.EQ.1) THEN
C           PROFILE WAS ABOVE NORMAL DEPTH AT START
            IF(YL.LT.YNL) THEN
C             CROSSING NORMAL DEPTH.  FORCE YL TO BE AT NORMAL DEPTH
              YL = YNL
              YT = MIN(D, YL)
              CALL XLKTAL
     I                   (ADL,
     M                    YT,
     O                    AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL,
     O                    DALPL)
            ENDIF
          ELSE
C           PROFILE AT START WAS BELOW NORMAL DEPTH
            IF(YL.GT.YNL) THEN
C             CROSSING NORMAL DEPTH.  FORCE YL TO BE AT NORMAL DEPTH
              YL = YNL
              YT = MIN(D, YL)
              CALL XLKTAL
     I                   (ADL,
     M                    YT,
     O                    AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL,
     O                    DALPL)
            ENDIF
          ENDIF
        ENDIF
C        ZL = ZBL + YL
        EL = YL + ALPL*(Q/AL)**2/GRAV2
        ER = EL
        ZBR = ZBL
        YR = YL
        XR = XL
 
        AR = AL
        ALPR = ALPL
        KR = KL
        ADR = ADL
        SBOLD = SB
        DR = DL
        AVEC(I) = AL
        KVEC(I) = KL
C        EVEC(I) = EL
        YVECSB(I) = YL
C        YNVEC(I) = YNL
C        YCVEC(I) = YC
C        SFVEC(I) = (Q/KL)**2
 500  CONTINUE
C501   CONTINUE
C      WRITE(STDOUT,*) ' SE=',SE,' Q=',Q
C      WRITE(STDOUT,60)
C      SUM = 0.0
C      DO 700 I=IS,ID
C        WRITE(STDOUT,62) I, XVEC(I), YVECSB(I), YCVEC(I), YNVEC(I),
C     A      ZBVEC(I), YVECSB(I) + ZBVEC(I), DVEC(I), SFVEC(I),
C     B      EVEC(I) + ZBVEC(I)
C        SUM = SUM + SFVEC(I)
C700   CONTINUE
C      SUM = SUM - .5*(SFVEC(IU) + SFVEC(ID)) -
C     A      (SFVEC(ID) - SFVEC(ID-1) - SFVEC(2) + SFVEC(1))/12.0
C      DE = SUM*ABS(XVEC(ID) - XVEC(IU))/FLOAT(ID - IU) + DH
C      WRITE(STDOUT,*) ' DE BY INTEGRATION=',DE, ' DE BY SUBTRACTION=',
C     A       EVEC(IU) + ZBVEC(IU) - EVEC(ID) - ZBVEC(ID)
 
C      WRITE(STDOUT,*) ' SFPSBE: TEL at entrance=', EL +  ZBL
      
      ZUP = YL + ZBL
      IS = IU
      SFLAG = 1
      RETURN
      END
C
C
C
      SUBROUTINE   SFPSPE
     I                   (STDOUT, IU, ID, Q, ZUP,
     O                    IS, ZDN, SFLAG)
 
C     + + + PURPOSE + + +
C     Compute a super critial steady flow profile in the channel defined
C     by NSEC, XVEC, ZBVEC, IU, and ID with the starting water surface
C     elevation at the upstream end(IU) given in ZUP.  The sequence
C     of depths is contained in YVECSP.  the elevation at the downstream end
C     of the channel is in ZDN and IS points to the last node at which
C     a super critical flow was found.  SFLAG is used to indicate the
C     results:  SFLAG = 1: computations completed and ZDN and YVECSP
C     contain the final results.  SFLAG = 0: computations not completed
C     subcritical flow was encountered.  It is believed that the
C     subcritical flow is physically caused by a less than steep slope.
C     The last valid depth is contained in YVECSP(IS).  SFLAG = -1:
C     computations not completed subcritical flow was encountered
C     and it is believed to be caused by DX being too large. The last
C     valid depth is again in YVECSP(IS).  SFLAG = -2: computations
C     not completed subcritical flow was encountered and it
C     is clear that DX was too large.  The last valid depth is
C     in YVECSP(IS).
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IS, IU, SFLAG, STDOUT
      REAL Q, ZDN, ZUP
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     Q      - Flowrate
C     ZUP    - Water surface elevation at section 2
C     IS     - Index of node for the last value computed
C     ZDN    - Water surface elevation at downstream end of barrel
C     SFLAG  - Flag describing results:
C                0-results incomplete.  Super critical flow encountered
C                  and flow is thought to be caused by a steep slope.
C                1-results complete
C               -1-results incomplete.  Super critical flow encountered
C                  and DX may be too large.
C               -2-results incomplete.  Super critical flow encountered
C                  and DX is too large.
C               -3-no solution found at the last point.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'sperc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADROLD, ADRS, FLAG, I
      REAL AC, ALPC, ALPR, AR, BETC, BETL, BETR, DALPC, DALPL, DALPR,
     A     DBETC, DBETL, DBETR, DKC, DKL, DKR, DTC, DTL, DTR, EC, EL,
     B     ER, FL, FR, JC, JL, JR, KC, KR, QCL, QN, SB, SFM, TC, TL, TR,
     C     YC, YL, YR, YSUP, YT, ZBL, ZBR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL SPER
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDCDE, GETTBN, REGFLT, SPER, XLKT22, XLKTAL, LENSTR,
     A         GET_TABID
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:599* INITIAL DEPTH=',F10.2,' <= 0 IN SFPSPE.')
 54   FORMAT(' *BUG:XXX* INITIAL CONDITION IN SFPSPE HAS FR=',F7.2,
     A       ' < 1.')
 56   FORMAT(' *BUG:XXX* TABID=',A,' UNDERFLOW SEEKING SUPER',
     A      ' CRITICAL SOLUTION.')
 58   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' SFPSPE.')
 60   FORMAT(' *BUG:XXX* REGFLT TAKES MORE THAN 100 ITERATIONS:SFPSPE.')
C***********************************************************************
      QR = Q
      QL = Q
      ADROLD = -1
      ZBL = ZBVEC(IU)
      YL = ZUP - ZBL
      IF(YL.LE.0.0) THEN
        WRITE(STDOUT, 50) YL
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      ADRS = NSEC(IU)
      D = DVEC(IU)
      YT = MIN(D, YL)
      CALL XLKT22
     I           (ADRS,
     M            YT,
     O            AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL, DALPL,
     O            QCL)
 
      FR = (QL/QCL)**2
      IF(FR.LT.0.995) THEN
        WRITE(STDOUT, 54) FR
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C     AT THIS POINT THE INITIAL CONDITION IS VALID. START THE LOOP
 
      EL = YL + ALPL*(QL/AL)**2/GRAV2
      YVECSP(IU) = YL
C     MAKE FIRST ESTIMATE OF YC
      YC = YL
      DO 500 I=IU+1,ID
        D = DVEC(I)
        SE = SEVEC(I)
        DX = ABS(XVEC(I-1) - XVEC(I))
        ZBR = ZBVEC(I)
        KACC = KA(I)
        KDEC = KD(I)
C       FIND CRITICAL DEPTH AT THE CURRENT LOCATION. USE LAST
C       AVAILABLE DEPTH AS THE FIRST ESTIMATE.
 
        ADRS = NSEC(I)
        ADR = ADRS
        IF(ADR.NE.ADROLD.OR.QL.NE.QR) THEN
          CALL FNDCDE
     I               (STDOUT, ADRS, QR,
     M                YC)
 
C         COMPUTE SPECIFIC ENERGY AT CRITICAL DEPTH
 
          CALL XLKTAL
     I               (ADRS,
     M                YC,
     O                AC, TC, DTC, JC, KC, DKC, BETC, DBETC, ALPC,
     O                DALPC)
 
          EC = YC + ALPC*(QR/AC)**2/GRAV2
        ENDIF
C       COMPUTE THE AVERAGE FRICTION SLOPE FOR THE ELEMENT
C       ASSUMING CRITICAL FLOW AT THE CURRENT NODE
 
        SFM = QL*QR/(KC*KL)
 
C       COMPUTE THE CONSTANT PART OF THE EQUATION
 
        RHS = EL + ZBL - ZBR
 
C       CHECK FOR EXISTENCE OF A SUPER CRITICAL SOLUTION
 
        IF(EC + DX*(SFM+SE).GT.RHS) THEN
C         NO SUPER CRITICAL SOLUTION EXISTS.
C          WRITE(STDOUT,*)  ' EC + DX*(SFM+SE)=',EC + DX*(SFM+SE)
C          WRITE(STDOUT,*) ' RHS=',RHS
C
C         Store current value of critical depth in the current
C         node's position.  Needed in some applications.
          YVECSP(I) = YC
C         Do diagnostic dump to check for existence.
C          WRITE(STDOUT,*) ' CHECKING FOR A SOLUTION'
C          DO 199 YT = 0.125*YC, 1.25*YC, YC/64.0
C            FT = SPER(YT)
C            WRITE(STDOUT,64) YT, FT
C199       CONTINUE
          SB = (ZBL - ZBR)/DX
          IS = I - 1
          IF(SB.LE.0.0) THEN
C           SUBCRITICAL FLOW IS PHYSICALLY POSSIBLE
            SFLAG = 0
          ELSE
            QN = KC*SQRT(SB)
            IF(QN.LT.QR) THEN
C             SUBCRITICAL FLOW IS PHYSICALLY POSSIBLE
              SFLAG = 0
            ELSE
C             CHANNEL SLOPE IS  STEEP. SUBCRITICAL FLOW IS
C             A COMPUTATIONAL ARTIFACT.
              SFLAG = -2
            ENDIF
          ENDIF
          RETURN
        ELSE
C         SUPER CRITICAL SOLUTION EXISTS-- FIND IT.
C         AT YC THE RESIDUAL FUNCTION IS > 0. SEARCH FOR RESIDUAL
C         FUNCTION < 0 STARTING AT YL IF IT IS LESS THAN YC.
 
          YT = YC
          FR = SPER(YT)
          IF(FR.LT.0.0) THEN
 
C            WRITE(STDOUT,*) ' SFPSPE: I=',I
C            WRITE(STDOUT,*) ' FR =',FR,' < 0.0 AT YC=',YC
C            WRITE(STDOUT,*) ' CHECKING FOR A SOLUTION'
C            DO 200 YT = 0.125*YC, 1.25*YC, YC/64.0
C              FT = SPER(YT)
C              WRITE(STDOUT,64) YT, FT
C200         CONTINUE
 
            YVECSP(I) = YC
C            IF(ABS(FR).LE.0.01) THEN
C              IS = I
C            ELSE
              IS = I - 1
C            ENDIF
            SFLAG = 0
            RETURN
          ENDIF
          IF(YL.LT.YC) THEN
            YSUP = YL
          ELSE
            YSUP = 0.5*YC
          ENDIF
 
 100      CONTINUE
            FL = SPER(YSUP)
            IF(FL.LT.0.0) THEN
              GOTO 110
            ELSE
              FR = FL
              YT = YSUP
              YSUP = 0.5*YSUP
              IF(YSUP.LE.EPSARG) THEN
                TABID = GET_TABID(GETTBN(ADRS))
                WRITE(STDOUT,56) TABID(1:LENSTR(TABID))
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 100
            ENDIF
 110      CONTINUE
 
C         WE HAVE A SIGN CHANGE IN THE INTERVAL (YSUP, YC)
 
          CALL REGFLT
     I               (EPSARG, EPSF, SPER,
     M                YSUP, YT, FL, FR,
     O                YR, FLAG)
          IF(FLAG.EQ.1) THEN
            WRITE(STDOUT,58)
            STOP 'Abnormal stop. Errors found.'
          ELSEIF(FLAG.EQ.2) THEN
            WRITE(STDOUT,60)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          YT = MIN(D, YR)
          CALL XLKTAL
     I               (ADRS,
     M                YT,
     O                AR, TR, DTR, JR, KR, DKR, BETR, DBETR, ALPR,
     O                DALPR)
 
          ER = YR + ALPR*(QR/AR)**2/GRAV2
          EL = ER
          ZBL = ZBR
          YL = YR
          QL = QR
          AL = AR
          ALPL = ALPR
          ADROLD = ADR
          KL = KR
          YVECSP(I) = YR
        ENDIF
 500  CONTINUE
 
      ZDN = YR + ZBR
      IS = ID
      SFLAG = 1
      RETURN
      END
C
C
C
      SUBROUTINE   SFPTY1
     I                   (STDOUT, IU, ID, DH, Q, ZUP,
     O                    IS, ZDN, SFLAG)
 
C     + + + PURPOSE + + +
C     Compute a steady flow profile in a culvert with type 1 flow
C     starting at the critical depth at the entrance and
C     computing subcritical flow to the exit of the culvert.
C     Used to define the limit of type 1 flow.
C     DH gives estimated head loss  that occurs in the
C     submerging flow.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IS, IU, SFLAG, STDOUT
      REAL DH, Q, ZDN, ZUP
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     DH     - Entrance loss for the culvert
C     Q      - Flowrate
C     ZUP    - Water surface elevation at section 2
C     IS     - Index of node for the last value computed
C     ZDN    - Water surface elevation at downstream end of barrel
C     SFLAG  - Flag describing results:
C                0-results incomplete.  Super critical flow encountered
C                  and flow is thought to be caused by a steep slope.
C                1-results complete
C               -1-results incomplete.  Super critical flow encountered
C                  and DX may be too large.
C               -2-results incomplete.  Super critical flow encountered
C                  and DX is too large.
C               -3-no solution found at the last point.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'sperc.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADROLD, ADRS, FLAG, I
      REAL AC, ALPC, ALPR, AR, BETC, BETL, BETR, DALPC, DALPL, DALPR,
     A     DBETC, DBETL, DBETR, DKC, DKL, DKR, DTC, DTL, DTR, EC, EL,
     B     ER, FL, FR, JC, JL, JR, KC, KR, QCL, SFM, TC, TL, TR, YC, YL,
     C     YMAX, YR, YSUB, YT, ZBL, ZBR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, DBLE, MAX, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL FMXARG, SPER
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FMXARG, FNDCDE, GETTBN, REGFLT, SPER, XLKT22, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* INITIAL DEPTH=',F10.2,' <= 0 IN SFPTY1.')
 56   FORMAT(' *BUG:XXX* TABID=',A,' OVERFLOW SEEKING ',
     A      ' SUBCRITICAL SOLUTION.')
 58   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' SFPTY1.')
 60   FORMAT(' *BUG:XXX* REGFLT TAKES MORE THAN 100 ITERATIONS:SFPTY1')
C***********************************************************************
      SFLAG = 1
      SE = DH/ABS(XVEC(IU) - XVEC(ID))
C      WRITE(STDOUT,*) ' SFPTY1: SE=',SE
 
      QR = Q
      QL = Q
      ADROLD = -1
      ZBL = ZBVEC(IU)
      YL = DBLE(ZUP) - DBLE(ZBL)
      IF(YL.LE.0.0) THEN
        WRITE(STDOUT, 50) YL
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      ADRS = NSEC(IU)
      CALL XLKT22
     I           (ADRS,
     M            YL,
     O            AL, TL, DTL, JL, KL, DKL, BETL, DBETL, ALPL, DALPL,
     O            QCL)
 
C     COMPUTE SQUARE OF THE FROUDE NUMBER
 
      FR = (QL/QCL)**2
C      IF(FR.LT.0.990) THEN
C        WRITE(STDOUT, 54)  FR
C        STOP 'Abnormal stop. Errors found.'
C      ENDIF
 
C     AT THIS POINT THE INITIAL CONDITION IS VALID. START THE LOOP
 
      EL = YL + ALPL*(QL/AL)**2/GRAV2
      YVECSP(IU) = YL
C      EVEC(IU) = EL
C      SFVEC(IU) = (QL/KL)**2
C     MAKE FIRST ESTIMATE OF YC
      YC = YL
      DO 500 I=IU+1,ID
        D = DVEC(I)
        IS = I
        DX = ABS(XVEC(I-1) - XVEC(I))
        ZBR = ZBVEC(I)
        KACC = KA(I)
        KDEC = KD(I)
 
C       FIND CRITICAL DEPTH AT THE CURRENT LOCATION. USE LAST
C       AVAILABLE DEPTH AS THE FIRST ESTIMATE.
 
        ADRS = NSEC(I)
        ADR = ADRS
        IF(ADR.NE.ADROLD.OR.QL.NE.QR) THEN
          CALL FNDCDE
     I               (STDOUT, ADRS, QR,
     M                YC)
 
C         COMPUTE SPECIFIC ENERGY AT CRITICAL DEPTH
 
          CALL XLKTAL
     I               (ADRS,
     M                YC,
     O                AC, TC, DTC, JC, KC, DKC, BETC, DBETC, ALPC,
     O                DALPC)
 
          EC = YC + ALPC*(QR/AC)**2/GRAV2
 
        ENDIF
C       COMPUTE THE AVERAGE FRICTION SLOPE FOR THE ELEMENT
C       ASSUMING CRITICAL FLOW AT THE CURRENT NODE
 
        SFM = QL*QR/(KC*KL)
 
C       COMPUTE THE CONSTANT PART OF THE EQUATION
 
        RHS = EL + ZBL - ZBR - DX*SE
 
C       FORCE A UNIQUE SUBCRITICAL SOLUTION.  DX MAY HAVE TO BE
C       REDUCED TO OBTAIN IT.
C        WRITE(STDOUT,*) 'SFPTY1: RHS=',RHS,' EC=',EC,' SFM=',SFM
C        WRITE(STDOUT,*) ' DIFF=',EC + DX*SFM - RHS
C        WRITE(STDOUT,*) ' EL=',EL,' ZBL-ZBR=',ZBL-ZBR
C        WRITE(STDOUT,*) ' YC=',YC,' YL=',YL
        IF(ABS((EC + DX*SFM - RHS)/MAX(1.0,ABS(RHS))).LE.EPSF) THEN
C         Take the solution to be critical depth.
          YR = YC
        ELSEIF(EC + DX*SFM.GT.RHS) THEN
C         NO UNIQUE SUBCRITICAL SOLUTION EXISTS OR NO SUBCRITICAL
C         SOLUTION EXISTS.
 
C          WRITE(STDOUT,*) ' SFPTY1: EC + DX*SFM=',EC + DX*SFM,
C     A                ' RHS=',RHS,' AT I=',I
C          WRITE(STDOUT,*) ' SFPTY1: NO SOLUTION.'
          SFLAG = 0
          IS = I - 1
          RETURN
        ELSE
C         UNIQUE SUBCRITICAL SOLUTION EXISTS-- FIND IT.
C         AT YC THE RESIDUAL FUNCTION IS > 0. SEARCH FOR RESIDUAL
C         FUNCTION < 0 STARTING AT YL IF IT IS GREATER THAN YC.
 
          FR = SPER(YC)
          IF(FR.LT.0.0) THEN
C           Treat as no solution.
            SFLAG = 0
            IS = I - 1
            RETURN
          ENDIF
          IF(YL.GT.YC) THEN
            YSUB = YL
          ELSE
            YSUB = YC
          ENDIF
          YMAX = FMXARG(ADRS)
 100      CONTINUE
            FL = SPER(YSUB)
            IF(FL.GT.0.0) THEN
              YSUB = 0.9*YSUB + 0.1*YMAX
              IF(ABS(YSUB - YMAX).LE.EPSARG) THEN
                TABID = GET_TABID(GETTBN(ADRS))
                WRITE(STDOUT,56) TABID(1:LENSTR(TABID))
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 100
            ENDIF
 
C         WE HAVE A SIGN CHANGE IN THE INTERVAL (YSUB, YC)
 
          YT = YC
          CALL REGFLT
     I               (0.0, 5.E-6, SPER,
     M                YSUB, YT, FL, FR,
     O                YR, FLAG)
          IF(FLAG.EQ.1) THEN
            WRITE(STDOUT,58)
            STOP 'Abnormal stop. Errors found.'
          ELSEIF(FLAG.EQ.2) THEN
            WRITE(STDOUT,60)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
        ENDIF
        YT = MIN(YR, D)
        CALL XLKTAL
     I             (ADRS,
     M              YT,
     O              AR, TR, DTR, JR, KR, DKR, BETR, DBETR, ALPR, DALPR)
 
        ER = YR + ALPR*(QR/AR)**2/GRAV2
        EL = ER
        ZBL = ZBR
        YL = YR
        QL = QR
        AL = AR
        ALPL = ALPR
        ADROLD = ADR
        KL = KR
        YVECSP(I) = YR
C        EVEC(I) = ER
C        YCVEC(I) = YC
C        YNVEC(I) = 0.0
C        SFVEC(I) = (QR/KR)**2
 500  CONTINUE
C      YCVEC(IU) = YCVEC(IU+1)
C      YNVEC(IU) = YNVEC(IU+1)
 
 
C      WRITE(STDOUT,*) ' SE=',SE,' Q=',Q
C      WRITE(STDOUT,61)
C      SUM = 0.0
C      DO 700 I=IU,IS
C        WRITE(STDOUT,62) I, XVEC(I), YVECSP(I), YCVEC(I), YNVEC(I),
C     A      ZBVEC(I), YVECSP(I) + ZBVEC(I), DVEC(I), SFVEC(I),
C     B      EVEC(I) + ZBVEC(I)
C        SUM = SUM + SFVEC(I)
C700   CONTINUE
C      SUM = SUM - .5*(SFVEC(IU) + SFVEC(ID)) -
C     A      (SFVEC(ID) - SFVEC(ID-1) - SFVEC(2) + SFVEC(1))/12.0
C      DE = SUM*ABS(XVEC(ID) - XVEC(IU))/FLOAT(ID - IU) + DH
C      WRITE(STDOUT,*) ' DE BY INTEGRATION=',DE, ' DE BY SUBTRACTION=',
C     A       EVEC(IU) + ZBVEC(IU) - EVEC(ID) - ZBVEC(ID)
 
      ZDN = YR + ZBR
      IS = ID
      SFLAG = 1
      RETURN
      END
