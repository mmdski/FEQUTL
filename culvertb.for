C
C
C
      REAL FUNCTION   RTY61
     I                     (Z)
 
C     + + + PURPOSE + + +
C     Compute residual funtion for finding the boundary of type 2 and
C     type 61 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Z
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Z      - water surface elevation being sought
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'rty6c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ISB, KNT, SFLAG
      REAL ARATIO, DH, DIV, DUP, EXPFAC, FDRDW, NUM, Q, QT, YT, ZSBRDF,
     A     ZT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, FNDCDE, GETFRF, LKTA, SFPSBE, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* No convergence on flow in fun. RTY61')
 52   FORMAT(/,' *ERR:696* Initial submergence of roadway flow starts',
     A    ' at elev.=',F10.3,/11X,'but critical-flow elev.=',F10.3,
     B    ' in culvert exit.',/,11X,'Flow over the road cannot be',
     C    ' combined with culvert flow when',/,11X,'critical flow in',
     D    ' the culvert exit causes submergence of flow over the road.')
C***********************************************************************
C     Compute the flow in the culvert.  Flow over the road may be
C     involved.
      CALL GETFRF
     I           (Z,
     O            ZSBRDF, FDRDW)
      Y1L = Z - ZB1
      Z1L = Z
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1L,
     O            A1L, T1L, DT1L, J1L, K1L, DK1L, BET1L, DBET1L, ALP1L,
     O            DALP1L)
 
      IF(A1L.GE.A2FULL) THEN
        EXPFAC = 0.0
      ELSE
        ARATIO = A1L/A2FULL
        IF(ARATIO.GT.0.95) THEN
          EXPFAC = 20.0*(1.0 - ARATIO)*APPEXP
        ELSE
          EXPFAC = APPEXP
        ENDIF
      ENDIF
      Q = 0
      KNT = 0
      DIV = (1.0  - EXPFAC - (ALP1L*(1.0 - EXPFAC) - APPLOS)*
     A      (A2FULL/A1L)**2 + APPLEN*GRAV2*A2FULL**2/(K1L*K2FULL))
C      WRITE(OUTUN,*) ' RTY61: DIV=',DIV,' Z=',Z,' A1L=',A1L,
C     A               ' EXPFAC=',EXPFAC,' A2FULL=',A2FULL,
C     B               ' WFRDF=',WFRDF
      IF(DIV.LE.0.0) THEN
        SBFLAG = -2
        RETURN
      ENDIF
 100  CONTINUE
        NUM = GRAV2*(Z - YAT2 - ZB2 - WFRDF*(
     A     Q*APPLEN/(K1L*K2FULL) - (2.*Q + WFRDF)/(GRAV2*A1L**2)*
     B     (ALP1L*(1.0 - EXPFAC) - APPLOS)))
C        WRITE(OUTUN,*) ' RTY61: YAT2=',YAT2,' ZB2=',ZB2,' K1L=',K1L,
C     A                 ' K2FULL=',K2FULL,' APPLOS=',APPLOS
C        WRITE(OUTUN,*) ' NUM=',NUM
        IF(NUM.LT.0.0) THEN
          SBFLAG = -1
          RETURN
        ENDIF
C        WRITE(OUTUN,*) 'RTY61: A2FULL=',A2FULL,' NUM=',NUM,' DIV=',DIV
C        WRITE(OUTUN,*) ' Z=',Z,' YAT2+ZB2=',YAT2+ZB2,' WFRDF=',WFRDF,
C     A          ' A1L=',A1L,' K1L=',K1L,' K2FULL=',K2FULL,
C     B          ' ALP1L=',ALP1L,' APPLOS=',APPLOS,' GRAV2=',GRAV2,
C     C          ' APPLEN=',APPLEN,' Q=',Q
        QT = A2FULL*SQRT(NUM/DIV)
        IF(WFRDF.GT.0.0) THEN
          KNT = KNT + 1
          IF(KNT.GT.100) THEN
            WRITE(OUTUN,50)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
C         Iterate to find the effect of roadflow.
          IF(ABS(Q - QT)/QT.GT.EPSF) THEN
            Q = QT
            GOTO 100
          ENDIF
        ENDIF
C      WRITE(OUTUN,*) ' QT=',QT
      Q1L = QT + WFRDF
      Q2L = QT
      Q3L = QT
 
      DUP = DVEC(IUP)
C     COMPUTE CRITICAL DEPTH AT SECTION 3
 
      Y3PART = 0.6*DUP
      CALL FNDCDE
     I           (OUTUN, ADRXS3, Q3L,
     M            Y3PART)
      Z3PART = ZB3 + Y3PART
      IF(Z3PART.GT.ZSBRDF) THEN
C       Critical depth in the culvert exit drowns free flow
C       over the road.  Cannot handle at this time.
        WRITE(OUTUN,52) ZSBRDF + ZDAT, Z3PART + ZDAT
        EF = 1
        RETURN
      ENDIF
 
C      WRITE(OUTUN,*) ' RTY61:  Q2L=',Q2L,' Y3PART=',Y3PART
 
      CALL LKTA
     I         (ADRXS3,
     M          Y3PART,
     O          A3PART)
      IF(FTYPE.EQ.2) THEN
        C123 = FCD123(OUTUN, 2, CLASS, DUP, Z)
        CD = DEGCON(C123, A1L, A3PART)
        AVH = A3PART
      ELSE
        CD = C46
        AVH = A2FULL
      ENDIF
 
      DH = (1.0/CD**2 - 1.0)*(Q3L/AVH)**2/GRAV2
 
C      WRITE(OUTUN,*) ' RTY61: A1=',A1L,' Z3PART=',Z3PART,
C     A              ' CD=',CD,' A3PART=',A3PART,' DH=',DH
 
      CALL SFPSBE
     I           (OUTUN, IUP, IDN, DH, QT, Z3PART,
     O            ISB, ZT, SFLAG)
      IF(SFLAG.EQ.1) THEN
        YT = ZT - ZB2
        SBFLAG = 0
      ELSE
        SBFLAG = -1
C        WRITE(OUTUN,*) ' SFP FAILED: ISB=',ISB
        RETURN
      ENDIF
 
      RTY61 = (YT - YAT2)/YAT2
 
C      WRITE(OUTUN,*) ' RTY61 EXIT: RTY61=',RTY61,' FTYPE=',FTYPE,
C     A               ' Q=',QT,' Z=',Z,' SFLAG=',SFLAG
      RETURN
      END
C
C
C
      REAL FUNCTION   RTY1
     I                    (Y)
 
C     + + + PURPOSE + + +
C     Compute the residual function for type 1 flow
 
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
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'rty1c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL APPFAC, ARATIO, CDIN, Q1C, Q2C, VH1, VH2
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, LKTQC, XLKT22
C***********************************************************************
C     GET VALUES AT SECTION 2
      Y2 = Y
      CALL XLKT22
     I           (ADRXS2,
     M            Y2,
     O            A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2,
     O            Q2C)
 
      Q2 = Q2C
 
C     ADD IN THE FREE FLOW OVER THE ROADWAY
      Q1 = Q2 + WFRDF
 
C     COMPUTE SQUARE OF FROUDE NUMBER IN APPROACH SECTION
      CALL LKTQC
     I          (ADRXS1,
     M           Y1,
     O           Q1C)
 
      FRSQ = (Q1/Q1C)**2
 
      Z2 = ZB2 + Y2
 
      VH1 = (Q1/A1)**2/GRAV2
      VH2 = (Q2/A2)**2/GRAV2
 
C      WRITE(OUTUN1,*) ' VH1=',VH1,' VH2=',VH2
C      WRITE(OUTUN1,*) ' A1=',A1,' A2=',A2,' GRAV2=',GRAV2
C     CHECK FOR EXPANSION INSTEAD OF CONTRACTION.
C      IF(VH2.GT.VH1) THEN
      IF(A1.GT.A2) THEN
C       WE HAVE A CONTRACTION(I.E. ACCELERATION OF FLOW)
C       DEFINE THE COEF OF DISHCHARGE
C        CONF = 1
        C123 = FCD123(OUTUN1, FTYPE, CLASS1, DU, Z1)
C        WRITE(OUTUN1,*) ' RTY1: C123=',C123,' FTYPE=',FTYPE 
C       MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C       CHANNEL CONTRACTION
 
        CD = DEGCON(C123, A1, A2)
C        WRITE(OUTUN1,*) ' RTY1: CD=',CD,' A1=',A1,' A2=',A2 
        RTY1 = (ALP1 - APPLOS)*VH1 - (ALP2 + (1./CD**2 -1.))*VH2
     A         + Z1 - Z2 - APPLEN*Q1*Q2/(K1*K2)
      ELSE
C        CONF = 0
C       WE HAVE AN EXPANSION(I.E. DECELERATION OF FLOW)
        ARATIO = A1/A2
        IF(ARATIO.GT.0.95) THEN
C         Interpolate coefficients to make the transistion between
C         the two cases smooth.  Define the coefficient for standard
C         type 1 case.
 
          C123 = FCD123(OUTUN1, FTYPE, CLASS1, DU, Z1)
 
C         MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C         CHANNEL CONTRACTION
 
          CDIN = DEGCON(C123, A1, A2)
 
          CD = CDIN + 20.0*(1.0 - ARATIO)*(0.98 - CDIN)
          APPFAC = 20.0*(1.0 - ARATIO)*APPEXP
        ELSE
          CD = 0.98
          APPFAC = APPEXP
        ENDIF
 
        RTY1 = (ALP1 - APPLOS)*VH1 - (ALP2 + (1./CD**2 -1.))*VH2
     A         + Z1 - Z2 - APPLEN*Q1*Q2/(K1*K2) -
     B              APPFAC*(ALP1*VH1 - ALP2*VH2)
 
      ENDIF
      DHTY1 = (1./CD**2 -1.)*VH2
C      WRITE(OUTUN1,*) ' AT Y=',Y,' RTY1=',RTY1,' CD=',CD
C      WRITE(OUTUN1,*) ' Q1=',Q1,' Q2=',Q2,' FRSQ=',FRSQ
C      WRITE(OUTUN1,*) ' TYPE1 LOSS =',DHTY1,' CONF=',CONF

      RETURN
      END
C
C
C
      REAL FUNCTION   RY2GY1
     I                      (Y)
 
C     + + + PURPOSE + + +
C     Compute the residual function for type 1 flow
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - value of the unknown being sought
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'rty1c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL APPFAC, ARATIO, CDIN, Q1C, Q2C, VH1, VH2
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, LKTQC, XLKT22
C***********************************************************************
C     GET VALUES AT SECTION 2
      Y2L = Y
      CALL XLKT22
     I           (ADRXS2,
     M            Y2L,
     O            A2L, T2L, DT2L, J2L, K2L, DK2L, BET2L, DBET2L, ALP2L,
     O            DALP2L, Q2C)
 
      Q2L = Q2C
 
C     ADD IN THE FREE FLOW OVER THE ROADWAY
      Q1L = Q2L + WFRDF
 
C     COMPUTE SQUARE OF FROUDE NUMBER IN APPROACH SECTION
      CALL LKTQC
     I          (ADRXS1,
     M           Y1L,
     O           Q1C)
 
      FRSQ = (Q1L/Q1C)**2
 
      Z2L = ZB2 + Y2L
 
      VH1 = (Q1L/A1L)**2/GRAV2
      VH2 = (Q2L/A2L)**2/GRAV2
 
C      WRITE(OUTUN1,*) ' RY2GY1: VH1=',VH1,' VH2=',VH2
C      WRITE(OUTUN1,*) ' RY2GY1; A1L=',A1L,' A2L=',A2L,' GRAV2=',GRAV2
C     CHECK FOR EXPANSION INSTEAD OF CONTRACTION.
      IF(A1L.GT.A2L) THEN
C       WE HAVE A CONTRACTION(I.E. ACCELERATION OF FLOW)
C       DEFINE THE COEF OF DISHCHARGE
        C123 = FCD123(OUTUN1, FTYPE, CLASS1, DU, Z1L)
 
C       MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C       CHANNEL CONTRACTION
 
        CD = DEGCON(C123, A1L, A2L)
 
        RY2GY1 = (ALP1L - APPLOS)*VH1 - (ALP2L + (1./CD**2 -1.))*VH2
     A         + Z1L - Z2L - APPLEN*Q1L*Q2L/(K1L*K2L)
      ELSE
C       WE HAVE AN EXPANSION(I.E. DECELERATION OF FLOW)
        ARATIO = A1L/A2L
        IF(ARATIO.GT.0.95) THEN
C         Interpolate coefficients to make the transistion between
C         the two cases smooth.  Define the coefficient for standard
C         type 1 case.
 
          C123 = FCD123(OUTUN1, FTYPE, CLASS1, DU, Z1L)
 
C         MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C         CHANNEL CONTRACTION
 
          CDIN = DEGCON(C123, A1L, A2L)
 
          CD = CDIN + 20.0*(1.0 - ARATIO)*(0.98 - CDIN)
          APPFAC = 20.0*(1.0 - ARATIO)*APPEXP
        ELSE
          CD = 0.98
          APPFAC = APPEXP
        ENDIF
 
        RY2GY1 = (ALP1L - APPLOS)*VH1 - (ALP2L + (1./CD**2 -1.))*VH2
     A         + Z1L - Z2L - APPLEN*Q1L*Q2L/(K1L*K2L) -
     B              APPFAC*(ALP1L*VH1 - ALP2L*VH2)
 
      ENDIF

C      WRITE(OUTUN1,*) ' AT Y=',Y,' RY2GY1=',RY2GY1,' CD=',CD
C      WRITE(OUTUN1,*) ' Q1L=',Q1L,' Q2L=',Q2L,' FRSQ=',FRSQ
      
      RETURN
      END
C
C
C
      SUBROUTINE   FRFT0
     I                  (STDOUT, HDATUM, HUP, IU, ID, CULCLS,
     O                   EFLAG, TYPE, EXPFLG, ZSBRDF, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute a flow of type 0 for the given upstream head, HUP.
C     Primary purpose of type 0 is to eliminate computational
C     problems with subsequent types.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL FREED, HDATUM, HUP, QFREE, ZSBRDF
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TYPE   - Culvert flow type
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     ZSBRDF - water surface elevation at section 43 that begins
C              submergence of flow over the roadway
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rty1c.cmn'
      INCLUDE 'rty0c.cmn'
      INCLUDE 'rqvstw.cmn'
      INCLUDE 'typtrn.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CONFLG, FLAG, IS, KNT, NSFLAG, SFLAG
      REAL A1T, DDN, DH, DUP, F, FL, FR, QC1, QCMAX, RES, TEL1, TEL2,
     A     TEL3, Y, Y2MIN, YHIGH, YLOW, YT, TP
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123, RTY0, RTY1
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL APPRO, DEGCON, DPM26, FCD123, FNDCDE, LKTJ, LKTQC, RGF3,
     A         RTY0, RTY1, SFPSBE, XLKT22, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF3 CLAIMS NONE IN',
     A       ' FRFT0.')
 60   FORMAT(' *BUG:XXX* RGF3: MORE THAN 100 ITERATIONS IN FRFT0.')
 66   FORMAT(/,' Type 0 flow drowned by control at section 2.')
 68   FORMAT(/,' Type 0 flow drowned by control at section 3.')
 70   FORMAT(/,' Flow is type 0. Flow in barrel=',F10.2,' Road flow=',
     A         F10.2)
 72   FORMAT(/,' *ERR:699* No result at section 2 for flow type 0',
     A         ' after 100 tries.')
 76   FORMAT(/,' *WRN:575* Minimum depth at section 2=',F8.2,' reached',
     A          ' and no',/,11X,'positive residual for flow type 0',
     B          /,11X,'when holding road flow fixed.')
 78   FORMAT(/,' *WRN:576* Type 0 flow has tailwater at section 43',
     A    ' higher',/,11X,'than road crest or submergence limit by',
     B     F8.3)
 80   FORMAT(/,' *WRN:577* Type 0 flow also has tailwater at section',
     A  ' 43 higher',/,11X,'than elevation at section 1.  This is',
     B    ' peculiar',/,11X,' but not impossible.  Please review',
     C    ' input carefully.')
 82   FORMAT(/,' *WRN:578 Unable to complete profile when roadflow',
     A   ' is constant.')
 84   FORMAT(/,'*ERR:700* Barrel flows full with type 0 flow. ',
     A ' Culvert representation MUST be',/,11X,' changed.  Type 0 flow',
     B ' invalid at this flow level.',/,
     C   11X,'  See error message summary.')
 86   FORMAT(/,' Elev. of energy grade line at section 1=',F7.3,
     A ' > ',F7.3,', elev. at section 2.')
 88   FORMAT(' Elev. of energy grade line at section 1=',F7.3,
     A ' > ',F7.3,', elev. at section 3.')
 90   FORMAT(/,' *ERR:701* Flow over the road=',F10.1,' > critical ',
     A   'flow at section 1=',F10.1,/,11X,' for type 0 flow.  Culvert ',
     B   ' representation MUST be changed.',/,11X,
     C   '  See error message summary.')
 92   FORMAT(/,' Testing for flow type 0 at upstream head =',F8.3)
C***********************************************************************
C     DEFINE THE SPECIAL COMMON BLOCK VALUES
      OUTUN1 = STDOUT
      OUTUN = STDOUT
      OUTUN0 = STDOUT
      DU = DVEC(IU)
      DUP = DU
      DUP0 = DU
      CLASS0 = CULCLS
      CLASS1 = CULCLS
      CLASS = CULCLS
      IUP = IU
      IDN = ID
      IUP0 = IU
      IDN0 = ID
      DDN = DVEC(ID)
      DDN0 = DDN
 
C     Clear the flag for special departure reach treatment
      BETAF = 0.0
      BETA3 = 0.0
      ALPHA3 = 0.0
      EXPFLG = 1
      TYPE = 0
C     FIND FLOW RATE AT SECTION 1
 
      Y1 = HDATUM + HUP - ZB1
      Z1T= Y1 + ZB1
      CALL XLKT22
     I           (ADRXS1,
     M            Y1,
     O            A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1,
     O            QC1)
      TEL1 = Z1T + ALP1*(QC1/A1)**2/GRAV2
      A1T = A1
      Z1T0 = Z1T
C     Set approach area in common blocks
      ABASE = A1T
      A1T0 = A1T
      QFIXED = QC1
C      Z43MAX = 0.0
C      WRITE(STDOUT,*) ' FRFT0: Y1=',Y1,' QC1=',QC1,' WFRD=',WFRD
      WRITE(STDOUT,92) HUP
      Q2 = QC1 - WFRD
      IF(Q2.LE.0.0) THEN
        WRITE(STDOUT,90)  WFRD, QC1
        EFLAG = 1
        RETURN
      ELSE
 
        Q1 = QC1
        Q3 = Q2
        Q4 = Q1
 
C       TEST 1
C       Compute critical depth at the culvert entrance.  Test for
C       possibility of flow type 1.
 
        TP = TY1YTD*DUP
        CALL LKTQC
     I            (ADRXS2,
     M             TP,
     O             QCMAX)
        IF(Q2.GE.QCMAX) THEN
          TYPE = 1
        ELSE
C          Y2MAX = DUP
          Y2 = 0.5*DUP
          CALL FNDCDE
     I               (STDOUT, ADRXS2, Q2,
     M                Y2)
          CALL XLKTAL
     I               (ADRXS2,
     M                Y2,
     O                A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2,
     O                DALP2)
          Y2MIN = Y2
          TEL2 = ZB2 + Y2 + ALP2*(Q2/A2)**2/GRAV2
          IF(TEL2.GE.TEL1) THEN
C            Type 0 flow is drowned by control at section 2.
            TYPE = 1
          ELSE
            WRITE(STDOUT,86) TEL1, TEL2
C           Type 0 flow may be drowned by control at section 2.
            Z1 = Z1T
 
            RES = RTY1(Y2)
            IF(RES.LT.0.0) THEN
C             Control at section 2 does drown control at section 1.
              TYPE = 1
            ELSE
              TYPE = 0
            ENDIF
          ENDIF
        ENDIF
        IF(TYPE.EQ.1) THEN
C         Type 0 flow drowned by control at section 2.
          WRITE(STDOUT,66)
          RETURN
        ENDIF
 
C       TYPE 0 PASSED TEST 1.
 
C       TEST 2.  Compute critical depth at section 3 and check for
C       flow type 2.
        Y3 = 0.5*DDN
        CALL FNDCDE
     I             (STDOUT, ADRXS3, Q3,
     M              Y3)
        CALL XLKTAL
     I             (ADRXS3,
     M              Y3,
     O              A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3)
        Z3 = ZB3 + Y3
 
        TEL3 = Z3 + ALP3*(Q3/A3)**2/GRAV2
        IF(TEL3.GE.TEL1) THEN
C         Type 0 flow is drowned by control at section 3.
          TYPE = 1
        ELSE
          WRITE(STDOUT,88) TEL1, TEL3
C         Type 0 flow may be drowned by control at section 3.
C         Find the water level at section 2
 
C         DEFINE THE COEF OF DISCHARGE
 
          C123 = FCD123(STDOUT, 2, CULCLS, DVEC(IUP), Z1T)
 
          CD = DEGCON(C123, A1T, A3)
          DH = (1.0/CD**2 - 1.0)*(Q3/A3)**2/GRAV2
          VHL = 0.0
 
          CALL SFPSBE
     I               (STDOUT, IU, ID, DH, Q3, Z3,
     O                IS, Z2, SFLAG)
 
          IF(SFLAG.NE.1) THEN
C           Subcritical profile could not be completed.  Take this
C           to mean that the barrel slope is steep for the given
C           flow and that control at section 3 is not possible.
            TYPE = 0
          ELSE
C           Check the energy line values.
            Y2 = Z2 - ZB2
            CALL XLKTAL
     I                 (ADRXS2,
     M                  Y2,
     O                  A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2,
     O                  DALP2)
            TEL2 = ZB2 + Y2 + ALP2*(Q2/A2)**2/GRAV2
            IF(TEL2.GT.TEL1) THEN
C             Control at section 3 is possible and drowns type 0
              IF(Y2.GT.DUP) THEN
                TYPE = 6
              ELSE
                TYPE = 1
              ENDIF
            ELSE
C             Control at section 3 is possible and may drown type 0.
C             Return values  are in XS1COM set of variables.
              CALL APPRO
     I                  (STDOUT, CD, VHL, 0.0,
     O                   CONFLG, NSFLAG)
              IF(NSFLAG.EQ.1) THEN
C               No solution in the approach reach if the flow at
C               section 3 is critical and the flow is critical flow at
C               section 1.
                TYPE = 0
              ELSE
C               Solution possible for section 1 elevation.  Does it
C               drown the control at section 1?
                IF(Z1.GT.Z1T) THEN
                  IF(Y2.GT.DUP) THEN
                    TYPE = 6
                  ELSE
                    TYPE = 1
                  ENDIF
                ELSE
                  TYPE = 0
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
 
 
        IF(TYPE.NE.0) THEN
C         Type 0 flow drowned by control at section 3.
          WRITE(STDOUT,68)
          RETURN
        ENDIF
 
C       TYPE 0 PASSED TEST 2.  The only control left is at section 4.
C       However, this is type 7, computed only after a truly free flow
C       type has been found.  Therefore, we find the level at section
C       43, the immediate tailwater location for flow through
C       the culvert, that for the given flow, QC1, will cause
C       the upstream water level to be matched.  This will be
C       taken to drown the control at section 1.  Once this level
C       is determined, the drop to free flow is computed using
C       DPM26.  At this time it is possible that there is no
C       expansion in the departure reach.   This identifies flow
C       type 7 just as for the other flow types.
 
        WRITE(STDOUT,70) Q2, WFRD
        TYPE = 0
 
        F = -1.0
        IF(F.LE.0.0) THEN
C         There is a solution.
 
          YLOW = 0.0
          YHIGH = 0.0
          Y = Y2MIN
          KNT = 0
 100      CONTINUE
            KNT = KNT + 1
            IF(KNT.GT.100) THEN
              WRITE(OUTUN,72)
              GOTO 800
            ENDIF
 
            F = RTY0(Y)
 
C            WRITE(STDOUT,*) ' FRFT0: RTY0=',F,' AT Y=',Y
            IF(ABS(F).LE.EPSABS) THEN
C             Close enough.  Results in XS2COM
 
            ELSE
              IF(F.GT.0.0) THEN
                YHIGH = Y
                FR = F
                IF(YLOW.EQ.0.0) THEN
C                 Seek a negative residual.
                  YT = 1.1*Y
C                  IF(YT.GT.Y2MAX) THEN
C                    YT = 0.5*(Y + Y2MAX)
C                    IF(YT.EQ.Y2MAX) THEN
C                      WRITE(OUTUN,74) Y2MAX
C                      STOP 'Abnormal stop. Errors found.'
C                    ENDIF
C                  ENDIF
                  Y = YT
                  GOTO 100
                ELSE
                  GOTO 110
                ENDIF
              ELSE
                YLOW = Y
                FL = F
                IF(YHIGH.EQ.0.0) THEN
C                 Seek a positive residual.
                  YT = 0.9*Y
                  IF(YT.LT.Y2MIN) THEN
                    YT = 0.5*(Y + Y2MIN)
                    IF(YT.EQ.Y2MIN) THEN
                      WRITE(OUTUN,76) Y2MIN
                      GOTO 800
                    ENDIF
                  ENDIF
                  Y = YT
                  GOTO 100
                ELSE
                  GOTO 110
                ENDIF
              ENDIF
 110          CONTINUE
 
C             A sign change exists here.
              CALL RGF3
     I                 (1.E-6, EPSABS, RTY0,
     M                  YLOW, YHIGH, FL, FR,
     O                  Y2, FLAG)
              IF(FLAG.EQ.1) THEN
                WRITE(STDOUT, 54)
                STOP 'Abnormal stop. Errors found.'
              ELSEIF(FLAG.EQ.2) THEN
                WRITE(STDOUT,60)
                STOP 'Abnormal stop. Errors found.'
              ENDIF
            ENDIF
 
C          WRITE(STDOUT,*) ' FRFT0: Y2=',Y2
          IF(SFLAG0.EQ.0) THEN
C           Profile problems in RTY0
            WRITE(STDOUT,82)
            IF(WFRD.EQ.0.0) THEN
              ZSBRDF = Z43
            ENDIF
            GOTO 800
          ENDIF
          Z43 = Z3
C          WRITE(STDOUT,*) ' Z43=',Z43,' ZSBRDF=',ZSBRDF,' Z1T=',Z1T
          IF(Z43.GT.ZSBRDF) THEN
C           Tailwater is higher than the roadway crest.
            WRITE(STDOUT,78) Z43 - ZSBRDF
            IF(Z43.GT.Z1T) THEN
C             Tailwater is higher than the upstream head
              WRITE(STDOUT,80)
            ENDIF
C            Z43MAX = Z43
            GOTO 800
          ENDIF
          GOTO 900
        ENDIF
      ENDIF
 
      IF(TYPE.EQ.1) RETURN
 
 800  CONTINUE
        IF(WFRD.GT.0.0) THEN
          WRITE(STDOUT,*) '  Type 0 fails with flow over the road.'
        ELSE
          WRITE(STDOUT,*)  '  Type 0 fails with no flow over the road.'
        ENDIF
        RETURN
 
 900  CONTINUE
 
 
      Y43 = Z43 - ZB43
      Z43OLD = Z43
      Z3P = Z43
      Z3 = Z3P
      Y3 = Z3 - ZB3
      IF(Y2.GE.DUP.AND.Y3.GE.DDN) THEN
        WRITE(STDOUT,84)
        EFLAG = 1
        RETURN
      ENDIF
      TP = MIN(Y3,DDN)
      CALL XLKTAL
     I           (ADRXS3,
     M            TP,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3)
      CALL LKTJ
     I         (ADRXS3,
     M          Y3,
     O          J3)
 
C      WRITE(STDOUT,*) ' FRFT0: Z43=',Z43,' Q3=',Q3
      QFREE = Q3
      Q3FREE = Q3
      Y2FREE = Y2
      Y3FREE = Y3
      LFTYPE = 0
      LSTYPE = -1
 
 
C     Now find the drop to the water level at section 4.
 
      CALL DPM26
     I          (STDOUT, WFRD, MFRD, Z3,
     O           EXPFLG)
      IF(EXPFLG.EQ.0) THEN
        FREED = 0.0
        RETURN
      ENDIF
 
C     DEPARTURE SECTION VALUES ARE IN XS4COM
 
      FREED = Z1T - Z4
C      WRITE(OUTUN,*) ' FRFT0: FREED=',FREED
      RETURN
 
      END
C
C
C
      SUBROUTINE   TY1GY2
     I                   (STDOUT, IU, ID, DUP, CULCLS, YAT2, TRUEA1,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Compute results for flow type 1 when the depth at section 2
C     is given.  Compute the conditions to section 3 as well.
C     Used to establish limits for type 1 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, ID, IU, STDOUT
      REAL DUP, TRUEA1, YAT2
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     DUP    - vertical diameter of culvert barrel at upstream end
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     YAT2   - depth at section 2
C     TRUEA1 - area at section 1
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'y1gy2.cmn'
      INCLUDE 'grvcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, IS, SFLAG
      REAL AVH, C123, CD, DH, F, FH, FL, MAXARG, Q1CL, Q2CL, Y, Y3L,
     A     YHIGH, YLOW, YT, ZBEG, ZDN
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123, FMXARG, RY1GY2
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, FMXARG, LKTA, REGFLT, RY1GY2, SFPTY1,
     A         XLKT22
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:570* Zero depth before negative residual while',
     A         ' seeking Type 1 limit.')
 52   FORMAT(/,' *ERR:594* Maximum depth before positive residual',
     A         ' while seeking Type 1 limit.')
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' TY1GY2.')
 60   FORMAT(' *BUG:XXX* REGFLT: MORE THAN 100 ITERATIONS IN TY1GY2.')
 62   FORMAT(/,' *ERR:683* Flow at section 1 has Froude number=',F8.3,
     A    /,11X,'seeking Type 1 flow limit.')
 72   FORMAT('  Initial loss for type 3 submergence=',1PE10.3)
 74   FORMAT('  Final loss for type 3 submergence=',1PE10.3)
C***********************************************************************
C     Place values in common block for residual function.
      CLASS = CULCLS
      OUTUN = STDOUT
      D = DUP
 
C     Define values at section 2.  The flow is critical at this
C     section, so define flow in the culvert as well.
      Y2L = YAT2
      Z2L = ZB2 + Y2L
      CALL XLKT22
     I           (ADRXS2,
     M            Y2L,
     O            A2L, T2L, DT2L, J2L, K2L, DK2L, BET2L, DBET2L, ALP2L,
     O            DALP2L, Q2CL)
 
      Q2L = Q2CL
 
C     Make first estimate of the elevation at section 1
      Z1L = ZB2 + 1.5*Y2L
      Y1L = Z1L - ZB1
      IF(Y1L.LT.0.0) THEN
        Y1L = 0.2*DUP
      ENDIF
      Y = Y1L
      MAXARG = FMXARG(ADRXS1)
C     Search for a sign change in the residual function.
      YLOW = 0.0
      YHIGH = 0.0
 
 100  CONTINUE
        F = RY1GY2(Y)
        IF(ABS(F).GT.EPSABS) THEN
C         Residual is too large.  Search for a sign change.
          IF(F.GT.0.0) THEN
            YHIGH = Y
            FH = F
            IF(YLOW.EQ.0.0) THEN
C             Make depth at section 1 smaller to find a negative
C             residual.
              Y = 0.95*Y
              IF(Y.LE.EPSABS) THEN
                WRITE(STDOUT,50)
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 100
            ENDIF
          ELSE
            YLOW = Y
            FL = F
            IF(YHIGH.EQ.0.0) THEN
C             Make depth at section 1 larger to find a positive
C             residual.
              YT = 1.05*Y
              IF(YT.GT.MAXARG) THEN
                YT = 0.5*(Y + MAXARG)
                IF(MAXARG - YT.LE.EPSABS) THEN
                  WRITE(STDOUT,52)
                  STOP 'Abnormal stop. Errors found.'
                ENDIF
              ENDIF
              Y = YT
              GOTO 100
            ENDIF
          ENDIF
 
C         We have a sign change in the residual at this point.
C         Find a root.
 
          CALL REGFLT
     I               (0.0, EPSABS, RY1GY2,
     M                YLOW, YHIGH, FL, FH,
     O                Y, FLAG)
 
          IF(FLAG.EQ.1) THEN
            WRITE(STDOUT, 54)
            STOP 'Abnormal stop. Errors found.'
          ELSEIF(FLAG.EQ.2) THEN
            WRITE(STDOUT,60)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
        ENDIF
C     Find the final values at section 1.
      Y1L = Y
      Z1L = ZB1 + Y1L
      CALL XLKT22
     I           (ADRXS1,
     M            Y1L,
     O            A1L, T1L, DT1L, J1L, K1L, DK1L, BET1L, DBET1L, ALP1L,
     O            DALP1L, Q1CL)
C     Check for supercritical flow at section 1
      IF(Q1L/Q1CL.GT.1.02) THEN
        EFLAG = 1
        WRITE(STDOUT, 62) Q1L/Q1CL
        RETURN
      ENDIF
 
C     Find the submergence conditions at section 3.  Same method
C     as in FRFT1.
 
C     Find the depth at the first distance step using direct
C     integration to try to avoid continuing problems with
C     finding free drop for type 1 flow.
C      OFFSET = ABS(XVEC(IU+1) - XVEC(IU))
C      SB = (ZBVEC(IU) - ZBVEC(IU+1))/OFFSET
C      FFAC = 1.0 - KD(IU+1)
C      CALL BEGSFP(NSEC(IU), FFAC, Y2L, Q2L, OFFSET, SB, YSTART)
 
C      WRITE(STDOUT,*) ' TY1GY2: OFFSET=',OFFSET,' YSTART=',YSTART,
C     A                ' Y2=',Y2
C      ZBEG = YSTART + ZBVEC(IU+1)
      ZBEG = Z2L
      DH = 0.0
      CALL SFPTY1
     I           (STDOUT, IU, ID, DH, Q2L, ZBEG,
     O            IS, ZDN, FLAG)
      IF(FLAG.EQ.0) THEN
        WRITE(STDOUT,*) ' TY1GY2: Problem in SFPTY1 with DH=0.0'
        WRITE(STDOUT,*) ' IU=',IU,' ID=',ID,' IS=',IS
        WRITE(STDOUT,*) ' Z2L=',Z2L,' Y2L=',Y2L,' Q2L=',Q2L
        EFLAG = 1
      ELSE
C       Now that we have an estimate of the conditions at section
C       3 when the critical control at section 2 is just being
C       submerged, estimate the type 3 losses that should take
C       place.
        Y3L = ZDN - ZB3
        CALL LKTA
     I           (ADRXS3,
     M            Y3L,
     O            AVH)
C       Select the losses as in RQVSTW.
        C123 = FCD123(STDOUT, 3, CULCLS, DUP, Z1L)
        CD = DEGCON(C123, TRUEA1, AVH)
        DH = (1.0/CD**2 - 1.0)*(Q2L/AVH)**2/GRAV2
        WRITE(STDOUT,72) DH
 200    CONTINUE
          CALL SFPTY1
     I               (STDOUT, IU, ID, DH, Q2L, ZBEG,
     O                IS, ZDN, SFLAG)
          IF(SFLAG.EQ.0) THEN
C           If estimated losses cause computational problems, reduce
C           the losses and try again until the computations are
C           successful.
            DH = 0.9*DH
            GOTO 200
          ENDIF
        WRITE(STDOUT,74) DH
        Y3LTY1 =  ZDN - ZB3
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   TY1BDY
     I                   (STDOUT, DUP, DDN, IU, ID, CULCLS, TRUEA1,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Find the type 1 boundary, if it exists, and the values
C     needed for transition from type 1 to type 5 or from type 1 to
C     type 6.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, ID, IU, STDOUT
      REAL DDN, DUP, TRUEA1
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     DUP    - vertical diameter of culvert barrel at upstream end
C     DDN    - vertical diameter of culvert barrel at downstream end
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     TRUEA1 - area at section 1
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'rty1c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, IS, KNT, NSCMAT, SFLAG
      REAL A, DH, F, FDRDW, FH, FL, HTD, Q2LC, QC, SATMAT(8), TP,
     A     TY6LSS, Y, Y2LIM, Y3L, YAT2, YATMAT(8), YCAT3, YHIGH, YLOW,
     B     YT, YVC, ZAT3, ZAT43, ZBEG, ZDN, ZSBRDF
      CHARACTER JMPLOC*32
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123, RY2GY1, YOVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, FNDCC2, FNDCDE, FSCMAT, GETFRF, LKTA,
     A         LKTQC, REGFLT, RY2GY1, SFPTY1, SUPSUB, TY1GY2, XLKT22,
     B         XLKTAL, YOVERD
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:684* Negative depth at section 1 for type 1 flow',
     A     ' when head ratio limit',/,11X,'there=',F8.3,' and the',
     B     ' culvert vertical diameter=',F8.3)
 52   FORMAT(/,' *ERR:685* Culvert soffit reached in TY1BDY seeking',
     A          ' a negative residual.')
 54   FORMAT(/,' *ERR:686* Culvert invert reached in TY1BDY seeking',
     A            ' a positive residual.')
 56   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' TY1BDY.')
 60   FORMAT(' *BUG:XXX* REGFLT: MORE THAN 100 ITERATIONS IN TY1BDY.')
 62   FORMAT(/,' The depth limit at section 2 causes the',
     A  ' limiting head ratio=',F8.3,/,'  to exceed the',
     B   ' limiting ratio for type 1 flow=',F8.3)
 64   FORMAT(/,' The depth limit at section 2 causes the limiting',
     A        ' head ratio to be:',F8.3)
 66   FORMAT(/,' *BUG:XXX* More than 100 iterations and no sign change',
     A         ' in TY1BDY.')
 68   FORMAT(/,' Type 1 depth limit at section 2=',F8.3,', taken at',
     A         ' maximum value allowed.')
 70   FORMAT(/,' Type 1 depth limit at section 2=',F8.3,', defined by',
     A      ' match of bottom',/,'   slope and critical slope.')
 72   FORMAT(/,' Lower limit for type 1 flow at section 1=',F8.3)
 74   FORMAT(/,' Lower limit for type 1 flow does not exist.')
 76   FORMAT(/,' *WRN:534* Unable to force type 6 Cd to match type 1',
     A         ' flow at its limit.',/,11X,' Using type 6 Cd. ',
     B     'Manual adjustment of 2-D table may be needed.')
 78   FORMAT(/,' *ERR:687* Upper limit for type 1 flow does not',
     A  ' exist.  Type 0 flow found.',/,11X,'Culvert representation',
     B  ' MUST be changed.  See error message summary.')
 80   FORMAT(/,' Type 1 limits:',/,
     A       4X,'Depth at section 1 for lower type 1=',F8.3,/,
     B       4X,'Depth at section 1 for upper type 1=',F8.3,/,
     C       4X,'Depth at section 3 for upper type 1=',F8.3)
 82   FORMAT(/,' Parameters for transition from type 1 to type 5:',/,
     A        4X,'Contraction coefficient=',F8.3,/,
     B        4X,'Momentum-flux coefficient=',F8.3,/,
     C        4X,'Kinetic-energy flux coefficient=',F8.3)
 84   FORMAT(/,' Parameters for transition from type 1 to type 6:',/,
     A        4X,'Coefficient of discharge=',F8.3,/,
     B        4X,'Momentum-flux coefficient=',F8.3,/,
     C        4X,'Kinetic-energy flux coefficient=',F8.3)
 86   FORMAT('  Initial loss for type 3 submergence=',1PE10.3)
 88   FORMAT('  Final loss for type 3 submergence=',1PE10.3)
89    FORMAT(' Minimum depth at section 2 for critical flow match >',
     A       ' limit at section 2=',F8.3,'.',/,
     B       '  Type 1 flow not possible.')
90    format(/,' Warning: Sequence of conditions to this point',
     a /,'indicates possible computational failure.  If computations',
     b /,'fail, be sure to read all descriptions of error messages',
     c /,'that appear.  Representation of the culvert may need',
     d /,'adjustment to allow flow computations to continue.')
C***********************************************************************
      EFLAG = 0
C     Find if type 1 flow even exists and if so find values
C     for type 1 at its limit.
 
      CALL FSCMAT
     I           (STDOUT, ADRXS2, SZERO,
     O            EFLAG, NSCMAT, YATMAT, SATMAT)
 
      Y2LIM = DUP*TY1YTD
      IF(NSCMAT.EQ.0) THEN
C       Type 1 flow does not exist for this culvert.  This should not
C       happen often because we are called from FRFT1.
        CD1 = -1.0
        WRITE(STDOUT,*) ' Type 1 flow not possible.'
        RETURN
      ELSEIF(NSCMAT.LT.0) THEN
C       The bottom slope is so large that any flow will be
C       type 1.  Set the depth at section 2 to the limit allowed
C       for type 1 flow.
        YAT2 = Y2LIM
        WRITE(STDOUT,68) YAT2
      ELSE
        IF(NSCMAT.GT.1) THEN
C         More than one depth.  Assume the first one is the lower
C         limit for type 1 flow.  Find the value of elevation at section
C         1 for this limit to provide a basis for controlling computation
C         of type 2 flow at shallow depths.
          YAT2 = YATMAT(1)
          IF(YAT2.GT.Y2LIM) THEN
            CD1 = -1.0
            WRITE(STDOUT,89) Y2LIM
            RETURN
          ELSE
            IF(SATMAT(1).GT.0.0) THEN
C             Something wrong somewhere.
              WRITE(STDOUT,*) ' *BUG:XXX* Inconsistent slope for ',
     A                 'lower critical slope in TY1BDY.'
              STOP 'Abnormal stop. Errors found.'
            ELSE
              CALL TY1GY2
     I                   (STDOUT, IU, ID, DUP, CULCLS, YAT2, TRUEA1,
     O                    EFLAG)
              IF(EFLAG.EQ.0) THEN
                Z1TY1L = Z1L
                WRITE(STDOUT,72) Z1TY1L
              ELSE
                EFLAG = 0
                WRITE(STDOUT,74)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C     Take the maximum depth, check the slope
C     and then check against the maximum allowed.
      YAT2 = YATMAT(NSCMAT)
      IF(SATMAT(NSCMAT).LT.0.0) THEN
C       Something wrong somewhere.
        WRITE(STDOUT,*) ' *BUG:XXX* Inconsistent slope for critical',
     A         ' slope in TY1BDY.'
        STOP 'Abnormal stop. Errors found.'
      ELSE
        IF(YAT2.GT.Y2LIM) THEN
C         The uppermost match with critical slope is above
C         the maximum depth permitted for type 1.
          YAT2 = Y2LIM
          WRITE(STDOUT,68) YAT2
        ELSE
          WRITE(STDOUT,70) YAT2
        ENDIF
      ENDIF
 
C     Find the type 1 conditions that apply when the depth at section
C     2 is given.
 
      CALL TY1GY2
     I           (STDOUT, IU, ID, DUP, CULCLS, YAT2, TRUEA1,
     O            EFLAG)
 
      IF(EFLAG.NE.0) THEN
C       Invalid result.
        WRITE(STDOUT,78)
        RETURN
      ENDIF
 
C     Results are in section 1.  Check to make sure that the
C     limit on depth at section 2 does not lead to an elevation at
C     section 1 that exceeds the limit set there.
      HTD = (Z1L - ZB2)/DUP
      IF(HTD.GT.TY1HTD) THEN
C       The limit imposed at section 2 does exceed the limit imposed at
C       section 1 for type 1 flow.  Set the head at section 1 to its
C       limit and compute type 1 flow for that condition.  We cannot
C       use FRFT1 because we are called from FRFT1.
 
        WRITE(STDOUT,62) HTD, TY1HTD
        Z1L =  DUP*TY1HTD + ZB2
        Y1L = Z1L - ZB1
        IF(Y1.LT.0.0) THEN
          WRITE(STDOUT,50) TY1HTD, DUP
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        CALL XLKTAL
     I             (ADRXS1,
     M              Y1L,
     O              A1L, T1L, DT1L, J1L, K1L, DK1L, BET1L, DBET1L,
     O              ALP1L, DALP1L)
C       Compute the free flow over the roadway.
 
        CALL GETFRF
     I             (Z1L,
     O              ZSBRDF, FDRDW)
C       Set values in the common block for the type 1 residual function.
 
        CLASS1 = CULCLS
        DU = DUP
        FTYPE = 1
        OUTUN1 = STDOUT
 
        YLOW = 0.0
        YHIGH = 0.0
        Y = YAT2
        KNT = 0
 100    CONTINUE
          F = RY2GY1(Y)
C          WRITE(STDOUT,*) ' TY1BDY: Y=',Y,' F=',F
          KNT = KNT + 1
          IF(KNT.GT.100) THEN
            WRITE(STDOUT,66)
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          IF(ABS(F).GT.EPSABS) THEN
C           Residual too large.  Seek a change in sign for the
C           residual function.
            IF(F.GT.0.0) THEN
              YHIGH = Y
              FH = F
              IF(YLOW.EQ.0.0) THEN
C               Make depth at section 2 larger to find a negative
C               residual.
                YT =  1.05*Y
                IF(YT.GE.DUP) THEN
                  YT = 0.5*(Y + DUP)
                  IF(DUP - YT.LE.EPSABS) THEN
                    WRITE(STDOUT,52)
                    STOP 'Abnormal stop. Errors found.'
                  ENDIF
                ENDIF
                Y = YT
                GOTO 100
              ENDIF
            ELSE
              YLOW = Y
              FL = F
              IF(YHIGH.EQ.0.0) THEN
C               Make depth at section 2 smaller to find a positive
C               residual.
                Y = 0.95*Y
                IF(Y.LE.EPSABS) THEN
                  WRITE(STDOUT,54)
                  STOP 'Abnormal stop. Errors found.'
                ENDIF
                GOTO 100
              ENDIF
            ENDIF
C           We have a change in sign of the residual function here.
C           Find a root.
            CALL REGFLT
     I                 (EPSARG, EPSABS, RY2GY1,
     M                  YLOW, YHIGH, FL, FH,
     O                  Y, FLAG)
            IF(FLAG.EQ.1) THEN
              WRITE(STDOUT, 56)
              STOP 'Abnormal stop. Errors found.'
            ELSEIF(FLAG.EQ.2) THEN
              WRITE(STDOUT,60)
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
        Y2L = Y
        Z2L = Y2L + ZB2
        CALL XLKT22
     I             (ADRXS2,
     M              Y2L,
     O              A2L, T2L, DT2L, J2L, K2L, DK2L, BET2L, DBET2L,
     O              ALP2L, DALP2L, Q2LC)
 
C       Find the submergence conditions at section 3.  Same method
C       as in FRFT1.
C       Find the depth at the first distance step using direct
C       integration to try to avoid continuing problems with
C       finding free drop for type 1 flow.
C        OFFSET = ABS(XVEC(IU+1) - XVEC(IU))
C        SB = (ZBVEC(IU) - ZBVEC(IU+1))/OFFSET
C        FFAC = 1.0 - KD(IU+1)
C        CALL BEGSFP(NSEC(IU), FFAC, Y2L, Q2L, OFFSET, SB, YSTART)
C
C        WRITE(STDOUT,*) ' TY1GY2: OFFSET=',OFFSET,' YSTART=',YSTART,
C     A                  ' Y2=',Y2
C        ZBEG = YSTART + ZBVEC(IU+1)
        ZBEG = Z2L
        DH = 0.0
        CALL SFPTY1
     I             (STDOUT, IU, ID, DH, Q2L, ZBEG,
     O              IS, ZDN, FLAG)
 
        IF(FLAG.EQ.0) THEN
          WRITE(STDOUT,*) ' TY1BDY: Problem in SFPTY1 with DH=0.0'
          WRITE(STDOUT,*) ' ID=',ID,' IS=',IS
          STOP 'Abnormal stop. Errors found.'
        ELSE
C         Now that we have an estimate of the conditions at section
C         3 when the critical control at section 2 is just being
C         submerged, estimate the type 3 losses that should take
C         place.
          Y3L = ZDN - ZB3
          CALL LKTA
     I             (ADRXS3,
     M              Y3L,
     O              AVH)
C         Select the losses as in RQVSTW.
          C123 = FCD123(STDOUT, 3, CULCLS, DUP, Z1L)
          CD = DEGCON(C123, TRUEA1, AVH)
          DH = (1.0/CD**2 - 1.0)*(Q2L/AVH)**2/GRAV2
          WRITE(STDOUT,86) DH
 200      CONTINUE
            CALL SFPTY1
     I                 (STDOUT, IU, ID, DH, Q2L, ZBEG,
     O                  IS, ZDN, SFLAG)
            IF(SFLAG.EQ.0) THEN
C             If estimated losses cause computational problems, reduce
C             the losses and try again until the computations are
C             successful.
              DH = 0.9*DH
              GOTO 200
            ENDIF
          WRITE(STDOUT,88) DH
 
        ENDIF
        Y3LTY1 =  ZDN - ZB3
      ELSE
        WRITE(STDOUT,64) HTD
      ENDIF
 
      Z1TY1 = Z1L
      CD1 = CD
      TP = MIN(Y3LTY1,DDN)
      CALL LKTA
     I         (ADRXS3,
     M          TP,
     O          A3PART)
      Z3PART = Y3LTY1 + ZB3
      Q3L = Q2L
      WRITE(STDOUT,80) Z1TY1L - ZB1, Z1TY1 - ZB1, Y3LTY1
 
C     Now we have finally established a limit for type 1 flow and
C     the conditions that exist at that limit.  Now try to
C     establish the values needed for the transition to type 5 or
C     type 6 flow.
 

      IF(HTD.GT.1.0) THEN
C       Compute transition parameters only when head to opening
C       height ratio is greater than 1.0.  When the ratio is
C       small then type 2 is likely the predominant type and
C       its parameters will be used in the transition.
        IF(HHTYPE.EQ.5) THEN
C         High head type is type 5.  Compute parameters for
C         transition from type 1 to type 5.  Force the type 5
C         flow equation to match the limiting type 1 flow by
C         computing a contraction coefficient for that purpose.

          CALL FNDCC2
     I               (STDOUT, DUP, ZBVEC(IAT3D),
     O                CC1T5)
C          WRITE(STDOUT,*) ' TY1BDY: CC1T5=',CC1T5
 
C         Compute the profile to the end of the barrel
 
          YVC = DUP*YOVERD(NSEC(IAT3D), DUP, A2FULL, CC1T5)
 
          TY6LSS = ((1.0/C46**2 - 1.0))*(Q2L/A2FULL)**2/GRAV2
 
          CALL LKTQC
     I              (NSEC(IAT3D),
     M               YVC,
     O               QC)
C          WRITE(STDOUT,*) ' TY1BDY: DEPTH AT VC=',YVC,' FROUDE=',Q2L/QC
C         Estimate the critical depth.
          YCAT3 = YVC
          CALL FNDCDE
     I               (STDOUT, ADRXS3, Q2L,
     M                YCAT3)
C          WRITE(STDOUT,*) ' CRITICAL DEPTH AT SEC. 3=',YCAT3
          CALL SUPSUB
     I               (STDOUT, IAT3D, ID, YVC, Y3LTY1, Q2L, DUP, YCAT3,
     M                TY6LSS,
     O                ZAT3, ZAT43, JMPLOC)
 
C          WRITE(STDOUT,*) ' JMPLOC=',JMPLOC
          Z1TY51 = Z1TY1
C          Y3LTY5 = DDN*TY5SBF
          Y3TY51 = ZAT3 - ZB3
 
C         Compute a value of beta to apply to the type 5 equation
C         results that will closely match the type 1 equation
C         at the boundary between them.
          CALL LKTA
     I             (ADRXS3,
     M              Y3TY51,
     O              A)
          BT1T5 = A/A3PART
          AP1T5 = (A/A3PART)**2
          WRITE(STDOUT,82) CC1T5, BT1T5, AP1T5
 
        ELSE
C         Find parameters for transition from type 1 to type 6.
          TP = ((ALP1L - APPLOS)*(Q1L/A1L)**2/GRAV2 + Z1L -
     A      ((Q3L/A3FULL)**2/GRAV2 + Z3PART) - L23*(Q2L/K2FULL)**2
     B     - APPLEN*(Q1L*Q2L)/(K1L*K2FULL))*(GRAV2*(A2FULL/Q2L)**2)
          IF(TP.LT.0.0) THEN
C           Problem.  Cannot find a Cd that matches type 1 flow.
            WRITE(STDOUT,76)
            CD1T6 = C46
          ELSE
            CD1T6 = SQRT(1.0/(1.0 + TP))
          ENDIF
          BT1T6 = A3FULL/A3PART
          AP1T6 = (BT1T6)**2
          WRITE(STDOUT,84) CD1T6, BT1T6, AP1T6
 
        ENDIF
      else
c       CULVERT has computed the free-flow limits for type 1 flow.  However, the
c       limiting condition for the water-surface elevation at section 1 is such that 
c       the ratio of approach head to vertical diameter is < 1.0!  
        write(stdout, 90) 

      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   FRFT1
     I                  (STDOUT, HDATUM, HUP, DUP, DDN, IU, ID, CULCLS,
     I                   TRUEA1,
     O                   EFLAG, FG, TYPE, CONFLG, EXPFLG, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute flow type 1 for the given upstream head, HUP, and
C     the approach and hydraulic data are given in the
C     labelled common blocks.
C     FG gives the first guess of the depth for type 2 if
C     it follows.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CONFLG, EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL DDN, DUP, FG, FREED, HDATUM, HUP, QFREE, TRUEA1
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     HUP    - Head upstream
C     DUP    - vertical diameter of culvert barrel at upstream end
C     DDN    - vertical diameter of culvert barrel at downstream end
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     TRUEA1 - area at section 1
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FG     - First guess
C     TYPE   - Culvert flow type
C     CONFLG - CONFLG=0: flow contracts as it enters the culvert and
C              CONFLG=1: flow expands as it enters the culver
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'rty1c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, IS, ISB, SFLAG
      REAL DH, F, FAC, FDRDW, FH, FL, FRH, FRL, QN, SB, SBE, TP, 
     A     YARG, YE, YHIGH, YLOW, YL, YR, ZBEG, ZSBRDF, ZUP
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123, RTY1
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, DPM26, FCD123, FNDCDE, GETFRF, LKTJ, REGFLT,
     A         RTY1, SFPSBE, SFPTY1, TY1BDY, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT REGFLT CLAIMS NONE IN',
     A       ' FRFT1.')
 56   FORMAT(' *WRN:539* NO ROOT FOR TYPE 1 FLOW. TRYING TYPE 2.')
 60   FORMAT(' *BUG:XXX* REGFLT: MORE THAN 100 ITERATIONS IN FRFT1.')
 62   FORMAT(/,' Type 1 rejected.',/,
     A '  Normal flow in barrel at critical',
     B ' depth=',F8.2,' <  critical flow=',F8.2)
 64   FORMAT(/,' Type 1 accepted.',/,
     A '  Normal flow in barrel at critical',
     B ' depth=',F8.2,' > critical flow=',F8.2)
 66   FORMAT(/,' *WRN:579* Type 1 flow drowned assuming flow is type',
     A  ' 2.',/,11X,' Slope reduction in culvert barrel.')
 68   FORMAT(/,' *BUG:XXX* FRFT1: Type 2 flow successful when it',
     A  ' should be impossible. ',/,5X,' Barrel slope at entrance=',
     B  F10.4,' Barrel slope at exit=',F10.4)
 70   FORMAT(/,' *WRN:589* Flows are within 1 per cent.  Unavoidable',
     A   ' convergence differences',/,11X,'may cause CULVERT to make',
     B   ' the wrong branch and fail.  If this',/,11X,'occurs change',
     C   ' the upstream head to avoid the near flow match.')
 72   FORMAT('  Initial loss for type 3 submergence=',1PE10.3)
 74   FORMAT('  Final loss for type 3 submergence=',1PE10.3)
C***********************************************************************
      BETAF = 0.0
      BETA3 = 0.0
      ALPHA3 = 0.0
c      write(stdout,*) 'TY1BDY: CD1=',cd1
      IF(CD1.EQ.0.0) THEN
C       Type 1 limit not set yet.  Try to find it.
        CALL TY1BDY
     I             (STDOUT, DUP, DDN, IU, ID, CULCLS, TRUEA1,
     O              EFLAG)
 
C       Recompute flow over the roadway because TY1BDY may have
C       computed it for a different upstream head.
        CALL GETFRF
     I             (HUP+HDATUM,
     O              ZSBRDF, FDRDW)
 
        IF(CD1.LT.0.0) THEN
          TYPE = 2
          FG = 0.0
          RETURN
        ENDIF
      ENDIF
      IF(HUP + HDATUM.GT.Z1TY1) THEN
C        WRITE(STDOUT,*) ' FRFT1: Limit elevation=',Z1TY1,' exceeded.'
C        WRITE(STDOUT,*) ' Current elevation=',HUP + HDATUM
        IF(CD2.GT.0.0) THEN
          IF(HUP + HDATUM.GT.Z1TY2) THEN
            IF(HHTYPE.EQ.5) THEN
              TYPE = 5
            ELSE
              TYPE = 6
            ENDIF
          ELSE
            TYPE = 2
            IF(Z1TY2.GT.ZB1) THEN
C             We have a type 1 limit that  has a type 2 flow above it.
C             This means that the type 1 limit is critical flow over the
C             length of the culvert.   Thus the submergence level for
C             type 1 flow at this limit is critical depth.  In this case
C             the value of Y3LTY1 is critical depth.   Start the
C             process at a slightly higher depth.
              write(stdout,*) ' Y3LTY1=',y3lty1
              FG = Y3LTY1 + 0.1*(HUP + HDATUM - Z1TY2)
            ELSE
              FG = 0.0
            ENDIF
          ENDIF
        ELSE
C         We have not tried type 2 yet.  Do so to check on the
C         type 2 limit.
          WRITE(STDOUT,*) ' FRFT1: Trying to find type 2 limit.'
          TYPE = 2
          FG = Y3LTY1*1.01
        ENDIF
        RETURN
      ENDIF
C      WRITE(STDOUT,*) ' FRFT1: Z1=',Z1
      CLASS1 = CULCLS
      DU = DUP
      FTYPE = 1
      CONFLG = 1
      OUTUN1 = STDOUT
 
      Y3FREE = 0.0
      Y2FREE = 0.0
      Q3FREE = 0.0
C     FIND AN INTERVAL CONTAINING A ROOT.
 
      FAC = 0.01
      YHIGH = 0.0
      YLOW = 0.0

      FRL = -1.0
      FRH = -1.0
 100  CONTINUE
        YARG = FAC*(Z1 - ZB2)
        F = RTY1(YARG)
        IF(F.LE.0.0) THEN
          YLOW = YARG
          FL = F
          FRL = FRSQ
        ELSE
          YHIGH = YARG
          FH = F
          FRH = FRSQ
        ENDIF
        IF(YHIGH.EQ.0.0.OR.YLOW.EQ.0.0) THEN
          FAC = 1.1*FAC
          IF(FAC.GE.1.0) THEN
C           TRY TYPE 2
            WRITE(STDOUT,56)
            FREED = 0.0
            QFREE = Q2
            TYPE = 2
            FG = 0.0
            RETURN
          ENDIF
          IF(FRL.GT.0.0.AND.FRH.GT.0.0) THEN
            IF(MIN(FRL, FRH).GT.2.0) THEN
C             ASSUME NO SUBCRITICAL SOLUTION EXISTS
             WRITE(STDOUT,*) ' Type 1 rejected while seeking root. ',
     A              ' Minimum FRSQ=',MIN(FRL, FRH)
              TYPE = 2
              FREED = 0.0
              QFREE = Q2
              FG = 0.0
              RETURN
            ENDIF
          ENDIF 
          GOTO 100
        ENDIF
 
C     AT THIS POINT WE HAVE A SIGN CHANGE IN THE RESIDUAL.  FIND THE
C     ROOT IN THE INTERVAL. CHECK THE FROUDE NUMBER IN THE APPROACH
C     SECTION.
 
      CALL REGFLT
     I           (EPSARG, EPSABS, RTY1,
     M            YLOW, YHIGH, FL, FH,
     O            Y2, FLAG)
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 50)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      IF(ABS(FL).GT.EPSDIF) THEN
        WRITE(STDOUT,*) ' Flow type 1 residual at convergence=',
     A   ABS(FL), ' > ',EPSDIF,'  Result may be invalid.'
      ENDIF
C     THE FINAL VALUES ARE IN XS2COM VARIABLES
 
      QFREE = Q2
      Z2 = ZB2 + Y2
      FG = Y2
      Y2TY1 = Y2
C     CHECK FOR VALIDITY OF TYPE 1 FLOW AND IF VALID DETERMINE THE
C     DROP TO FREE FLOW.  FLOW MUST PASS 4 TESTS BEFORE IT IS
C     TYPE 1.

C      WRITE(STDOUT,*) ' FRFT1: CD=',CD,' for Type 1 flow' 
C      WRITE(STDOUT,*) ' FRFT1: A2=',A2,' A1=',A1
      IF(FRSQ.GT.1.0) THEN
C       Flow at section 1 is supercritical.  Flow type not 1.
        WRITE(STDOUT,*) ' Type 1 rejected after solution: FRSQ=',FRSQ
        TYPE = 2
        FREED = 0.0
        QFREE = Q2
        FG = 0.0
        RETURN
      ENDIF
 
      SB = (ZBVEC(IU) - ZBVEC(IU+1))/ABS(XVEC(IU) - XVEC(IU+1))
      IF(SB.LE.0.0) THEN
C        WRITE(STDOUT,*) ' TYPE 1 REJECTED. SB=',SB
        TYPE = 2
      ELSE
        QN = K2*SQRT(SB)
        IF(QN.LT.Q2) THEN
C         Slope is mild.  Super critical flow not supported.
          WRITE(STDOUT,62) QN, Q2
          IF(ABS(QN-Q2)/QN.LE.0.01) THEN
            WRITE(STDOUT,70)
          ENDIF
          TYPE = 2
        ELSE
          WRITE(STDOUT,64) QN, Q2
          IF(ABS(QN-Q2)/QN.LE.0.01) THEN
            WRITE(STDOUT,70)
          ENDIF
C         COMPUTE PROFILES TO CHECK FOR SLOPE CHANGES
 
 
C         FIND CRITICAL DEPTH AT EXIT AND COMPUTE SUBCRITICAL
C         PROFILE
 
          YE = Y2
          CALL FNDCDE
     I               (STDOUT, ADRXS3, Q2,
     M                YE)
          Y3 = YE
          Z3 = ZB3 + Y3
          CALL XLKTAL
     I               (ADRXS3,
     M                Y3,
     O                A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3,
     O                DALP3)
 
C         DEFINE THE COEF OF DISCHARGE
 
          C123 = FCD123(STDOUT, 2, CULCLS, DVEC(IU), Z1)
 
C         MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C         CHANNEL CONTRACTION
 
          CD = DEGCON(C123, A1, A3)
 
          Q3 = Q2
          DH = (1.0/CD**2 - 1.0)*(Q3/A3)**2/GRAV2
          CALL SFPSBE
     I               (STDOUT, IU, ID, DH, Q2, Z3,
     O                ISB, ZUP, SFLAG)
 
          IF(ISB.EQ.IU) THEN
C           SUB CRITICAL PROFILE HAS REACHED THE ENTRANCE.  THUS
C           CRITICAL DEPTH AT THE ENTRANCE IS DROWNED
C           BECAUSE THE CULVERT CANNOT SUSTAIN SUPERCRITICAL
C           FLOW THROUGHOUT ITS LENGTH
            SBE = (ZBVEC(ID-1) - ZBVEC(ID))/ABS(XVEC(ID-1) - XVEC(ID))
            IF(SBE.LT.SB) THEN
              WRITE(STDOUT,66)
              TYPE = 2
              FG = YE
            ELSE
              WRITE(STDOUT,68) SB, SBE
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(TYPE.EQ.2) THEN
        RETURN
      ELSE

C       DETERMINE THE DROP TO FREE FLOW.
C       COMPUTE SPECIAL STEADY FLOW PROFILE FROM ENTRANCE TO
C       EXIT AND THEN USE THE DEPARTURE REACH FROM EXIT TO SECTION
C       4 TO DETERMINE THE LEVEL IN SECTION 4 REQUIRED TO
C       DROWN TYPE 1 FLOW

       
C       Find the elevation of water surface at section 3 that
C       just drowns critical flow at section 2.  ZBEG gives the
C       elevation at the entrance for critical flow.

        ZBEG = Z2

        IF(1.EQ.2) THEN


C       Use bisection to find the tailwater that comes close to drowning
C       critical flow at the entrance.  
C       Search for a sign change.  We know that critical flow at the
C       culvert exit does not drown critical flow at the entrance. 
C       Therefore, the subcritical profile from the exit to the entrance
C       did not complete. 

        YL = YE
        Y3 = YE  

        Q3 = Q2
300     CONTINUE
          Y3 = 1.1*Y3
          Z3 = ZB3 + Y3
          CALL XLKTAL
     I               (ADRXS3,
     M                Y3,
     O                A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3,
     O                DALP3)
 
C         DEFINE THE COEF OF DISCHARGE
 
          C123 = FCD123(STDOUT, 3, CULCLS, DVEC(IU), Z1)
          CD = DEGCON(C123, A1, A3)
          DH = (1.0/CD**2 - 1.0)*(Q3/A3)**2/GRAV2
          CALL SFPSBE
     I               (STDOUT, IU, ID, DH, Q2, Z3,
     O                ISB, ZUP, SFLAG)
C          WRITE(STDOUT,*) ' FRFT1: In search of bracket Y3=',Y3
          IF(ISB.NE.IU) THEN
C           Profile did not reach entrance.
            YL = Y3
            GOTO 300
          ENDIF

C       Profile did reach the entrance.  Set YR and start the
C       bisection process.
        YR = Y3
310     CONTINUE
          Y3 = 0.5*(YL + YR)
          Z3 = ZB3 + Y3
          CALL XLKTAL
     I               (ADRXS3,
     M                Y3,
     O                A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3,
     O                DALP3)
 
C         DEFINE THE COEF OF DISCHARGE
 
          C123 = FCD123(STDOUT, 3, CULCLS, DVEC(IU), Z1)
          CD = DEGCON(C123, A1, A3)
          DH = (1.0/CD**2 - 1.0)*(Q3/A3)**2/GRAV2
          CALL SFPSBE
     I               (STDOUT, IU, ID, DH, Q2, Z3,
     O                ISB, ZUP, SFLAG)
          IF(ISB.NE.IU) THEN
C           Profile did not reach entrance.
            YL = Y3
          ELSE
C           Profile did reach entrance.
            YR = Y3
          ENDIF
C          WRITE(STDOUT,*) ' FRFT1: Sqeezing YL=',YL,' YR=',YR
          IF(ABS(YL - YR)/YR.GT.EPSF) GOTO 310

        ELSE

        DH = 0.0
        CALL SFPTY1
     I             (STDOUT, IU, ID, DH, Q2, ZBEG,
     O              IS, Z3, SFLAG)
        IF(SFLAG.EQ.0) THEN
          WRITE(STDOUT,*) ' FRFT1: Problem in SFPTY1 with DH=0.0'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
 
C       Now that we have an estimate of the conditions at section
C       3 when the critical control at section 2 is just being
C       submerged, estimate the type 3 losses that should take
C       place.

        Y3 = Z3 - ZB3
        TP = MIN(Y3,DDN)
        CALL XLKTAL
     I             (ADRXS3,
     M              TP,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3)

C       Select the losses as in RQVSTW.
        C123 = FCD123(STDOUT, 3, CULCLS, DUP, HDATUM + HUP)
C        WRITE(STDOUT,*) ' FRFT1: C123=',C123
        AVH = A3
C        AVH = A2
        CD = DEGCON(C123, TRUEA1, AVH)
C        WRITE(STDOUT,*) ' FRFT1: CD=',CD
        DH = (1.0/CD**2 - 1.0)*(Q2/AVH)**2/GRAV2
        WRITE(STDOUT,72) DH
 200    CONTINUE
          CALL SFPTY1
     I               (STDOUT, IU, ID, DH, Q2, ZBEG,
     O                IS, Z3, SFLAG)
          IF(SFLAG.EQ.0) THEN
C           If estimated losses cause computational problems, reduce
C           the losses and try again until the computations are
C           successful.
            DH = 0.9*DH
            GOTO 200
          ENDIF

        WRITE(STDOUT,74) DH

        ENDIF
        Y3 =  Z3 - ZB3
        TP = MIN(Y3,DDN)
        CALL XLKTAL
     I             (ADRXS3,
     M              TP,
     O              A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3)
        CALL LKTJ
     I           (ADRXS3,
     M            Y3,
     O            J3)
        Q3FREE = Q3
        Y3FREE = Y3
        Y2FREE = Y2
        LFTYPE = 1
        LSTYPE = -1
C       FIND ELEVATION IN THE DEPARTURE REACH
 
        Z3P = Z3
 
        CALL DPM26
     I            (STDOUT, WFRDF, MFRDF, Z3,
     O             EXPFLG)
 
        IF(EXPFLG.EQ.0) THEN
          FREED = 0.0
        ELSE
          FREED = Z1 - Z4
        ENDIF
        RETURN
      ENDIF
      END
C
C
C
      SUBROUTINE   FRFT2
     I                  (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O                   EFLAG, FG, TYPE, EXPFLG, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute flow of type 2 for the given upsteam head, HUP
 
C     FG is a first guess passed from type 1 flow or 0.0
C     if type 1 was not computed
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL FG, FREED, HDATUM, HUP, QFREE, ZDATUM
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     ZDATUM - datum for local elevation
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FG     - First guess
C     TYPE   - Culvert flow type
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'rty2c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, ICASE, IFLAG, KNT
      REAL A, DDN, DUP, FHIGH, FLOW, HW, QC, QN, R, TP, TY6LSS, YAT2,
     A     YCAT3, YHIGH, YLOW, YT, YVC, ZAT3, ZAT43
      CHARACTER JMPLOC*32
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FMXARG, RTY2, YOVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DPM26, F61BDY, FMXARG, FNDCC2, LKTA, LKTQC, RGF, RTY2,
     A         SUPSUB, XLKT22, XLKTAL, YOVERD
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF CLAIMS NONE IN',
     A       ' FRFT2.')
 52   FORMAT(/,' Lower type 1 limit does not exist.  Type 2 limit ',
     A       'does not exist.',/,' Increasing the approach length ',
     B       'may help if the approach',/,' area is about the same ',
     C       'size as the culvert flow area,',/,' or if the kinetic-',
     D       'energy-flux coefficient in the',/,' approach section ',
     E       'is significantly > 1.0.',//,' Trying high-head option.')
 56   FORMAT(' *WRN:540* No positive residual for type 2 flow.')
 58   FORMAT(' *WRN:541* No negative residual for type 2 flow.')
 60   FORMAT(' *BUG:XXX* RGF: More than 100 iterations in FRFT2.')
 66   FORMAT(/,' *WRN:552* Residual at convergence=',F10.3,' > ',F8.4)
 67   FORMAT(/,'  Type 2 may not be possible.  Trying high head.')
 68   FORMAT('  Current residual=',1PE10.3,' trying to continue.')
 70   FORMAT(/,' Type 2 flow no longer possible.  Trying high-head',
     A         ' flow.')
 76   FORMAT(/,' *WRN:580* Unable to force type 6 Cd to match type 2',
     A         ' flow at its limit.',/,11X,' Using type 6 Cd. ',
     B     'Manual adjustment of 2-D table may be needed.')
 78   FORMAT(/,' *WRN:581* Unable to force type 61 Cd to match type 2',
     A         ' flow at its limit.',/,11X,' Using type 6 Cd. ',
     B     'Manual adjustment of 2-D table may be needed.')
 80   FORMAT(/,'*WRN:590* Type 2 failure at ups head=',F8.3,' <= head',
     A   ' at lower limit',/,11X,'of Type 1=',F8.3,'  If heads are',
     B   ' close, increase ups head to',/,11X,
     C   'exceed Type 1 lower limit.')
 82   FORMAT(/,' *WRN:591* Type 2 failure at ups head=',F8.3,' => head',
     A   ' at lower limit',/,11X,'of Type 1=',F8.3,'  If heads are',
     B   ' close, increase ups head to',/,11X,'exceed Type 1 lower',
     C   ' limit by a larger amount.')
 84   FORMAT(/,' Type 2 limits:',/,
     A         4X,'Local elevation at section 1=', F8.3,/,
     B         4X,'Head/vertical diameter ratio at section 2=',F8.3,/,
     C         4X,'Depth at section 3=',F8.3)
 86   FORMAT(/,' Parameters for transition from type 2 to type 6:',/,
     A        4X,'Coefficient of discharge=',F8.3,/,
     B        4X,'Momentum-flux coefficient=',F8.3,/,
     C        4X,'Kinetic-energy flux coefficient=',F8.3)
 88   FORMAT(/,' Type 61 limits:',/,
     A        4X,'Local elevation at section 1=',F8.3,/,
     B        4X,'Head/vertical diameter ratio at sec. 2=', F8.3)
 90   FORMAT(/,' Parameters for transition from type 61 to type 6:',/,
     A        4X,'Coefficient of discharge=',F8.3,/,
     B        4X,'Momentum-flux coefficient=',F8.3,/,
     C        4X,'Kinetic-energy flux coefficient=',F8.3)
 92   FORMAT(/,' Parameters for transition from type 2 to type 5:',/,
     A        4X,'Contraction coefficient=',F8.3,/,
     B        4X,'Momentum-flux coefficient=',F8.3,/,
     C        4X,'Kinetic-energy flux coefficient=',F8.3,/,
     D        4X,'Depth at section 3=',F8.3)
C***********************************************************************
      BETAF = 0.0
      BETA3 = 0.0
      ALPHA3 = 0.0
      DUP = DVEC(IU)
      DDN = DVEC(ID)

C      WRITE(STDOUT,*) ' FRFT2: CD2=',CD2 
      IF(CD2.EQ.0.0) THEN
C       Type 2 limit has not been computed.  Try to find the limit
C       here.
        CALL F61BDY
     I             (STDOUT, 2, HDATUM, ZDATUM, HUP, IU, ID, CULCLS, DUP,
     O              EFLAG)
C        WRITE(STDOUT,*) ' FRF2: AFTER F61BDY: EFLAG=',EFLAG 
        IF(EFLAG.EQ.-1) THEN
C         Type 2 limit does not exist.  Set a value in CD2 to show
C         that the limit does not exist.
          WRITE(STDOUT,*) ' FRFT2: Type 2 limit above upper Type 1',
     A       ' limit does not exist.'
          CD2 = -1.0
          EFLAG = 0
          CD2T6 = -1.0
          IF(Z1TY1L.GT.ZB1) THEN
C           Lower type 1 limit exists.  Therefore, type 2 flow
C           may exist below the lower type 1 limit.
            WRITE(STDOUT,*) ' Lower type 1 limit exists.  Will',
     A         ' attempt Type 2 computations.'
 
          ELSE
C           Type 2 limit above type 1 limit does not exist. Lower
C           type 1 limit does not exist.  Not clear what is
C           happening here.  Try a  high-head option.
            WRITE(STDOUT,52) 
            TYPE = 6
            RETURN
          ENDIF
        ELSE
          CD2 = CD
          AVH2 = AVH
          Y3LTY2 = Y3PART
          WRITE(STDOUT,84) Z1TY2, (Z1TY2 - ZB2)/DUP, Y3LTY2
          IF(HHTYPE.EQ.6) THEN
            IF(ZB2.GE.ZB3) THEN
C             Type 61 limit does not exist but type 2 limit does.
C             Find parameters for transition from type 2 to type 6.
              TP = ((ALP1L - APPLOS)*(Q1L/A1L)**2/GRAV2 + Z1L -
     A          ((Q3L/A3FULL)**2/GRAV2 + Z3PART) - L23*(Q2L/K2FULL)**2
     B         - APPLEN*(Q1L*Q2L)/(K1L*K2FULL))*(GRAV2*(A2FULL/Q2L)**2)
              IF(TP.LT.0.0) THEN
                WRITE(STDOUT,76)
                CD2T6 = C46
              ELSE
                CD2T6 = SQRT(1.0/(1.0 + TP))
              ENDIF
              BT2T6 = A3FULL/A3PART
              AP2T6 = (BT2T6)**2
              WRITE(STDOUT,86) CD2T6,  BT2T6, AP2T6
              CD61 = -1.0
            ELSE
C             Type 2 flow exists and the barrel slope is adverse.
C             Therefore, type 61 exists!
              YAT2 = ZB3 + DDN - ZB2
              CALL F61BDY
     I                   (STDOUT, 61, HDATUM, ZDATUM, HUP, IU, ID,
     I                    CULCLS, YAT2,
     O                    IFLAG)
              IF(IFLAG.EQ.-1) THEN
                WRITE(STDOUT,*) ' Type 61 limit not found when it must',
     A                  ' exist.'
                STOP 'Abnormal stop. Errors found.'
              ELSE
                CD61 = CD
                AVH61 = AVH
C               Type 61 exists.  Therefore transition is between 61
C               and 6.
 
                TP = ((ALP1L - APPLOS)*(Q1L/A1L)**2/GRAV2 + Z1L -
     A           ((Q3L/A3FULL)**2/GRAV2 + Z3PART) - L23*(Q2L/K2FULL)**2
     B          - APPLEN*(Q1L*Q2L)/(K1L*K2FULL))*(GRAV2*(A2FULL/Q2L)**2)
                IF(TP.LT.0.0) THEN
                  WRITE(STDOUT,78)
                  CD61T6 = C46
                ELSE
                  CD61T6 = SQRT(1.0/(1.0 + TP))
                ENDIF
                BT61T6 = A3FULL/A3PART
                AP61T6 = (BT61T6)**2
                WRITE(STDOUT,88) Z1TY61, (Z1TY61 - ZB2)/DUP
                WRITE(STDOUT,90) CD61T6, BT61T6, AP61T6
              ENDIF
            ENDIF
          ELSE
C           High head type is type 5.  Compute parameters for
C           transition from type 2 to type 5.  Force the type 5
C           flow equation to match the limiting type 2 flow by
C           computing a contraction coefficient for that purpose.
 
            CALL FNDCC2
     I                 (STDOUT, DUP, ZBVEC(IAT3D),
     O                  CC2T5)
C            WRITE(STDOUT,*) ' FRFT2: CC2T5=',CC2T5
 
C           Compute the profile to the end of the barrel
 
            YVC = DUP*YOVERD(NSEC(IAT3D), DUP, A2FULL, CC2T5)
 
            TY6LSS = ((1.0/C46**2 - 1.0))*(Q2L/A2FULL)**2/GRAV2
 
            CALL LKTQC
     I                (NSEC(IAT3D),
     M                 YVC,
     O                 QC)
C            WRITE(STDOUT,*) ' FRFT2: DEPTH AT VC=',YVC,' FROUDE=',Q2L/QC
            YCAT3 = Y3LTY2
C            WRITE(STDOUT,*) ' CRITICAL DEPTH AT SEC. 3=',YCAT3
            CALL SUPSUB
     I                 (STDOUT, IAT3D, ID, YVC, Y3LTY2, Q2L, DUP, YCAT3,
     M                  TY6LSS,
     O                  ZAT3, ZAT43, JMPLOC)
 
C            WRITE(STDOUT,*) ' JMPLOC=',JMPLOC
            Z1TY52 = Z1TY2
C            Y3LTY5 = DDN*TY5SBF
            Y3TY52 = ZAT3 - ZB3
 
C           Compute a value of beta to apply to the type 5 equation
C           results that will closely match the type 2 equation
C           at the boundary between them.
            CALL LKTA
     I               (ADRXS3,
     M                Y3TY52,
     O                A)
            BT2T5 = A/A3PART
            AP2T5 = (A/A3PART)**2
 
            WRITE(STDOUT, 92) CC2T5, BT2T5, AP2T5, Y3TY52
          ENDIF
        ENDIF
      ENDIF
 
 
      MAXARG = FMXARG(ADRXS2)
C     DEFINE THE TRUE VALUES AT SECTION 1
 
      Z1TRUE = HDATUM + HUP
 
C     Type 2 flow may occur under a variety of cases.  There
C     are three levels at section 1 that are involved.  Z1TY1L
C     gives the lower limit for type 1 flow.  This is taken
C     to be the upper limit for type 2 flow.  Z1TY1 gives the
C     upper limt for type 1 flow.  Between Z1TY1L and Z1TY1 the
C     flow is type 1.  The last level is Z1TY2 giving the upper
C     limit for type 2 flow.  The upper level for type 1 flow is
C     complex.  If it is defined by critical slope, then there
C     may be type 2 flow at levels above it.  However, the critical
C     slope definition can lead to nonsense for pipe culverts
C     because the depth in the barrel is so close to the soffit
C     and the critical flow is so large that the water surface
C     elevation at section 1 becomes so large that the entrance
C     will be flowing full.  This happens because critical flow
C     loses its physical meaning as the soffit is approached in
C     a closed conduit that has converging walls.  The mathematical
C     meaning still exists and with the slot in the top, FEQUTL
C     should always be able to find a critical slope that will
C     exceed any bottom slope encounted in the field.
C     In principle the lower limit for type 1 flow should always
C     exist also because the critical slope approaches infinity
C     as the depth in the barrel approaches zero.  The minimum
C     depth tabulated for a closed conduit is 0.08 or less.
C     Flow depths this small are of little interest.  Therefore,
C     the lower limit for type 1 flow may not always be computed
C     even though it may exist.  The lower limit may not
C     exist because the flow at that level is type 0 and not
C     type 1.  Thus non-existence of the lower level is not
C     an error!
 
      ICASE = 0
      IF(Z1TY1L.GT.ZB1) THEN
C       Lower limit for type 1 flow exists.
        ICASE = ICASE + 1
      ENDIF
      IF(Z1TY1.GT.ZB1) THEN
C       Upper limit for type 1 flow exists.
        ICASE = ICASE + 2
      ENDIF
      IF(Z1TY2.GT.ZB1) THEN
C       Upper limit for type 2 flow exists.
        ICASE = ICASE + 4
      ENDIF
C      WRITE(STDOUT,*) ' FRFT2: CASE=',ICASE
C     There are 8 outcomes: 0 through 7.  Not all are possible here
C     but check for them all to find bugs.
      IF(ICASE.EQ.0) THEN
C       None of the limits exist.
        TYPE = 6
        RETURN
      ELSEIF(ICASE.EQ.1) THEN
C       Should not happen.  No upper limit for type 1 but there
C       is a lower limit.  Nonsense!
        WRITE(STDOUT,*) ' *BUG:XXX* ICASE=1 IN FRFT2.'
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(ICASE.EQ.2) THEN
C       Should not happen.  Only upper limit for type 1 exists.
C       Type 2 flow should not occur but here we are in FRFT2!
        WRITE(STDOUT,*) ' *BUG:XXX* ICASE=2 IN FRFT2.'
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(ICASE.EQ.3) THEN
C       Both type 1 limits exist. Type 2 upper limit does not.
        IF(Z1TRUE.GT.Z1TY1L) THEN
C         Type 1 elevation in FRFT2!  Happens if we go above
C         the upper limit for type 1 but we cannot find a type 2
C         limit above the upper type 1 and a lower type 1 limit
C         exists.
          WRITE(STDOUT,82) HUP, Z1TY1L - HDATUM
          WRITE(STDOUT,*) ' FRFT2: Type 2 appears impossible.',
     A      '   Trying high-head option'
          TYPE = 6
          RETURN
        ENDIF
      ELSEIF(ICASE.EQ.4) THEN
C       Only type 2 upper limit exists.
        IF(Z1TRUE.GT.Z1TY2) THEN
          TYPE = 6
          RETURN
        ENDIF
      ELSEIF(ICASE.EQ.5) THEN
C       Should not happen.  No upper limit for type 1 but there is
C       a lower limit.
        WRITE(STDOUT,*) ' *BUG:XXX* ICASE=5 IN FRFT2.'
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(ICASE.EQ.6) THEN
C       Both the type 1 upper limit and the type 2 upper limit exist.
C       The lower type 1 limit does not.
        IF(Z1TY2.LT.Z1TY1) THEN
C         Should not happen.
          WRITE(STDOUT,*) ' Z1TY2=',Z1TY2,' Z1TY1=', Z1TY1
          WRITE(STDOUT,*) ' *BUG:XXX* ICASE=6 IN FRFT2. INVALID',
     A           ' LIMIT RELATIONSHIP.'
          STOP 'Abnormal stop. Errors found.'
C        ELSEIF(Z1TRUE.LT.Z1TY1) THEN
C         WRITE(STDOUT,*) ' *BUG:XXX* ICASE=6 IN FRFT2. INVALID',
C     A           ' ELEVATION.'
C         STOP 'Abnormal stop. Errors found.'
        ELSEIF(Z1TRUE.GE.Z1TY2) THEN
          TYPE = 6
          RETURN
        ENDIF
C       If the current Z1 is between type1 and type 2 limits change
C       the first guess.
        IF(Z1TRUE.GE.Z1TY1.AND.Z1TRUE.LE.Z1TY2) THEN
          FG = Y3LTY1 + (Z1TRUE - Z1TY1)*(Y3LTY2 - Y3LTY1)/
     A                                     (Z1TY2 - Z1TY1)
        ENDIF
      ELSEIF(ICASE.EQ.7) THEN
C       All three limits exist.
        IF(Z1TY2.LT.Z1TY1) THEN
C         Should not happen.
          WRITE(STDOUT,*) ' *BUG:XXX* ICASE=7 IN FRFT2. INVALID',
     A           ' LIMIT RELATIONSHIP.'
          STOP 'Abnormal stop. Errors found.'
        ELSEIF(Z1TRUE.GT.Z1TY1L.AND.Z1TRUE.LE.Z1TY1) THEN
C         Should not be here.  Type 1 exists and the elevation
C         is in type 1 range.
          WRITE(STDOUT,*) ' *BUG:XXX* ICASE=7 IN FRFT2. TYPE 1',
     A           ' ELEVATION APPEARS.'
          STOP 'Abnormal stop. Errors found.'
        ELSEIF(Z1TRUE.GT.Z1TY2) THEN
          TYPE = 6
          RETURN
        ENDIF
C       If the current Z1 is between type1 and type 2 limits change
C       the first guess.
        IF(Z1TRUE.GE.Z1TY1.AND.Z1TRUE.LE.Z1TY2) THEN
          FG = Y3LTY1 + (Z1TRUE - Z1TY1)*(Y3LTY2 - Y3LTY1)/
     A                                     (Z1TY2 - Z1TY1)
        ENDIF
      ELSE
        WRITE(STDOUT,*) ' *BUG:XXX* ICASE ERROR. OUTSIDE RANGE 0:7'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      HW = Z1TRUE - ZB2
      R = HW/DUP
      IF(R.GT.1.5) THEN
C       HEAD WATER RATIO TOO LARGE FOR TYPE 2
C        WRITE(STDOUT,*) ' FRFT2: HEAD WATER RATIO=',R,' > 1.5'
        TYPE = 6
        RETURN
      ENDIF
      Y1 = Z1TRUE - ZB1
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1,
     O            A1TRUE, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1TR,
     O            DALP1)
 
C     SET THE SPECIAL COMMON BLOCK VALUES
      OUTUN = STDOUT
      IUP = IU
      IDN = ID
      CLASS = CULCLS
      SBOT = (ZBVEC(ID-1) - ZBVEC(ID))/ABS(XVEC(ID-1) - XVEC(ID))
 
C     ESTIMATE THE CRITICAL DEPTH AT SECTION 3 AND THEN SEARCH FOR
C     A SIGN CHANGE IN THE RESIDUAL. IF Y3 IS TOO LARGE THEN
C     THE RESIDUAL IS POSITIVE.
 
 
      HW = MIN(Z1TRUE - ZB2, Z1TRUE - ZB3)
C      WRITE(STDOUT,*) ' FRFT2: FG=',FG,' Y3FREE=',Y3FREE
C     SEARCH FOR A NEGATIVE RESIDUAL
      IF(FG.GT.0.0) THEN
        YLOW = FG
      ELSE
        IF(Y3FREE.GT.0.0) THEN
          YLOW = Y3FREE
        ELSE
          YLOW = .3*HW
        ENDIF
      ENDIF
      IF(YLOW.GE.DDN) THEN
        YLOW = 0.8*DDN
      ENDIF
C     CHK TO MAKE SURE THAT THIS VALUE OF DEPTH CAN SUSTAIN A
C     CRITICAL FLOW
 
      IF(SBOT.GT.0.0) THEN
        SBOT = SQRT(SBOT)
C       WRITE(STDOUT,*) ' SQRT(SBOT)=',SBOT
 105    CONTINUE
C         FIND CRITICAL FLOW AT SECTION 3
C         GET VALUES AT SECTION 3
          CALL XLKT22
     I               (ADRXS3,
     M                YLOW,
     O                A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3,
     O                DALP3, Q3C)
 
C         COMPUTE THE NORMAL FLOW AT THIS DEPTH.
          QN = K3*SBOT
 
C         WRITE(STDOUT,*) ' YLOW=',YLOW,' Q3C=',Q3C,' QN=',QN
          IF(QN.GT.Q3C) THEN
            YLOW = 0.9*YLOW
            IF(YLOW.LT.EPSARG) THEN
              WRITE(STDOUT,*) ' CRITICAL FLOW IMPOSSIBLE FOR TYPE 2.'
              WRITE(STDOUT,*) ' POSSIBLE ERROR IN TYPE 1 CHK.'
              TYPE = 1
              RETURN
            ENDIF
            GOTO 105
          ENDIF
          FG = YLOW
      ENDIF
 
      YHIGH = 0.0
      FHIGH = 1000.
      KNT = 0
 100  CONTINUE
        FLOW = RTY2(YLOW)
        IF(SBFLAG.EQ.3) THEN
          WRITE(STDOUT,70)
          IF(Z1TRUE.LE.Z1TY1L) THEN
C           We are below the lower limit for the type 1 flow
C           region.  We may be close to that limit and no
C           solution may exist for type 2 flow close to this
C           limit.
            WRITE(STDOUT,80) HUP, Z1TY1L - HDATUM
          ENDIF
          TYPE = 6
          RETURN
        ENDIF
        IF(SBFLAG.EQ.2) THEN
C         Profile failure.  Assume that the flow has been made too
C         small.
          IF(YHIGH.GT.0.0) THEN
C           A residual has already been found and it is positive.
            YT = 0.5*(YHIGH + YLOW)
          ELSE
            YT = 1.05*YLOW
          ENDIF
          IF(YT.GT.DDN) THEN
            YT = .5*(YLOW + DDN)
          ENDIF
          YLOW = YT
          KNT = KNT + 1
          IF(KNT.GT.100) THEN
            WRITE(STDOUT,*) ' Profile failure for Type 2 flow.'
            WRITE(STDOUT,*) ' While increasing critical depth.'
            IF(YHIGH.GT.0.0.AND.ABS(FHIGH).LT.EPSDIF) THEN
              WRITE(OUTUN,68) FHIGH
              GOTO 200
            ENDIF
            IF(Z1TRUE.LE.Z1TY1L) THEN
C             We are below the lower limit for the type 1 flow
C             region.  We may be close to that limit and no
C             solution may exist for type 2 flow close to this
C             limit.
              WRITE(STDOUT,80) HUP, Z1TY1L - HDATUM
            ENDIF
            IF(Y2.GE.DUP) THEN
              WRITE(STDOUT,70)
              TYPE = 6
              RETURN
            ENDIF
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          GOTO 100
        ENDIF
        IF(FLOW.GE.0.0) THEN
          YHIGH = YLOW
          FHIGH = FLOW
          YLOW = 0.95*YLOW
          IF(YLOW.LE.EPSARG) THEN
C           A NEGATIVE RESIDUAL SHOULD BE FOUND.
            WRITE(STDOUT,58)
            FREED = 0.0
            QFREE = Q3
            TYPE = 5
            RETURN
          ENDIF
          GOTO 100
        ENDIF
      IF(YHIGH.EQ.0.0) THEN
C       MUST SEARCH FOR POSTIVE RESIDUAL BECAUSE NONE FOUND WHEN
C       SEARCHING FOR A NEGATIVE RESIDUAL.  YLOW is defined here
C       so start with it.
 
        YHIGH = 1.05*YLOW
        IF(YHIGH.GT.DDN) THEN
          YHIGH = 0.5*(YLOW + YHIGH)
        ENDIF
C        WRITE(STDOUT,*) ' FRFT2: SEARCH FOR POSITIVE RESID. YHIGH=',
C     A                   YHIGH
        KNT = 0
 110    CONTINUE
          FHIGH = RTY2(YHIGH)
          IF(SBFLAG.EQ.3) THEN
            WRITE(STDOUT,70)
            IF(Z1TRUE.LE.Z1TY1L) THEN
C             We are below the lower limit for the type 1 flow
C             region.  We may be close to that limit and no
C             solution may exist for type 2 flow close to this
C             limit.
              WRITE(STDOUT,80) HUP, Z1TY1L - HDATUM
            ENDIF
            TYPE = 6
            RETURN
          ENDIF
          IF(SBFLAG.EQ.2) THEN
C           Profile failure.  Assume that flow is still too small.
            YT = 1.05*YHIGH
            IF(YT.GT.DDN) THEN
              YT = 0.5*(YHIGH + DDN)
            ENDIF
            YHIGH = YT
            KNT = KNT + 1
            IF(YHIGH-DDN.GE.EPSABS.OR.KNT.GT.100) THEN
              WRITE(STDOUT,*) ' Profile failure for Type 2 flow.'
              WRITE(STDOUT,*) ' While increasing critical depth'
              WRITE(STDOUT,*) ' searching for positive residual.'
              WRITE(STDOUT,*) ' Probable program failure.'
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            GOTO 110
          ENDIF
          IF(FHIGH.LE.0.0) THEN
            FLOW = FHIGH
            YLOW = YHIGH
            YHIGH = 0.98*YHIGH + 0.02*DDN
            IF(ABS(YHIGH - DDN).LE.EPSARG) THEN
              WRITE(STDOUT,56)
              FREED = 0.0
              QFREE = Q3
              TYPE = 5
              RETURN
            ENDIF
          GOTO 110
          ENDIF
      ENDIF
 
C     WE NOW HAVE A SIGN CHANGE IN (YLOW, YHIGH)
 
 
C      WRITE(STDOUT,*) ' YLOW=',YLOW,' FLOW=',FLOW
C      WRITE(STDOUT,*) ' YHIGH=',YHIGH,' FHIGH=',FHIGH
      CALL RGF
     I        (EPSARG, EPSABS, RTY2,
     M         YLOW, YHIGH, FLOW, FHIGH,
     O         Y3, FLAG)
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 50)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.3) THEN
C       SUBCRITICAL PROFILE NOT POSSIBLE. THUS TRY TYPE 6
        WRITE(STDOUT,70)
        IF(Z1TRUE.LE.Z1TY1L) THEN
C         We are below the lower limit for the type 1 flow
C         region.  We may be close to that limit and no
C         solution may exist for type 2 flow close to this
C         limit.
          WRITE(STDOUT,80) HUP, Z1TY1L - HDATUM
        ENDIF
        TYPE = 6
        RETURN
      ENDIF
 
 200  CONTINUE
      QFREE = Q3
      Q3FREE = Q3
      Y2FREE = Y2
      Z2 = ZB2 + Y2
      Y3FREE = Y3
      LFTYPE = 2
      LSTYPE = -1
C      WRITE(STDOUT,*) ' FRFT2: Y3=',Y3,' Q3=',Q3
 
C     CHK THAT THE RESIDUAL AT POINT OF CONVERGENCE IS SMALL.  RESIDUAL
C     IS A DIFFERENCE IN ELEVATION AT SECTION 1 AND IS GIVEN BY FLOW.
      IF(ABS(FLOW).GT.EPSDIF) THEN
        WRITE(STDOUT,66) ABS(FLOW), EPSDIF
        IF(Z1TY1.GT.ZB1) THEN
          WRITE(STDOUT,67)
C         Try a high head flow.  We may be trapped near the soffit of
C         the conduit where a narrow band of type 2 flow may exist.
          TYPE = 6
          RETURN
        ENDIF
      ENDIF
 
C     NOW COMPUTE THE ELEVATION IN THE DEPARTURE SECTION REQUIRED
C     TO DROWN CRITICAL DEPTH AT THE CULVERT EXIT.
 
      Z3P = Z3
      Z3PEST = Z3P
      CALL DPM26
     I          (STDOUT, WFRDF, MFRDF, Z3,
     O           EXPFLG)
 
 
      IF(EXPFLG.EQ.0) THEN
        RETURN
      ENDIF
 
C     DEPARTURE SECTION VALUES ARE IN XS4COM
 
      FREED = Z1TRUE - Z4
 
      RETURN
 
      END
C
C
C
      SUBROUTINE   FNDQ5
     I                  (STDOUT, HDATUM, HUP, IU, CULCLS, RBVAL, QRF,
     O                   YVC, Q5, CC)
 
C     + + + PURPOSE + + +
C     Find the type 5 flow rate for the culvert.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER IU, STDOUT
      REAL CC, HDATUM, HUP, Q5, QRF, RBVAL, YVC
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     RBVAL  - rounding/beveling value
C     QRF    - Flow over the roadway
C     YVC    - depth at vena contracta
C     Q5     - free flow of type 5
C     CC     - contraction coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'appcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL AVC, CD, DEN, DUP, DZVC, HDRAT, HEAD, HEADVC, NUM, NUMFAC, Q,
     A     QC1, QC2, QT, Z1TRUE
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL YOVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDCC, FNDCD5, XLKT22, YOVERD
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(/,' *BUG:XXX* No convergence on Q in subroutine FNDQ5.')
C***********************************************************************
      DUP = DVEC(IU)
      CALL XLKT22
     I           (ADRXS2,
     M            DUP,
     O            A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2,
     O            QC2)
      Z1TRUE = HDATUM + HUP
      Y1 = Z1TRUE - ZB1
      CALL XLKT22
     I           (ADRXS1,
     M            Y1,
     O            A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1,
     O            QC1)
      HEAD =  HUP + HDATUM - ZB2
 
      HDRAT = HEAD/DUP
C     Find the discharge coefficient for type 5 flow.
 
      CALL FNDCD5
     I           (STDOUT, CULCLS, HDRAT, RBVAL,
     O            CD)
 
C     Deduce the contraction coefficient from the CD given by
C     the USGS tables.  These appear to ignore approach velocity,
C     approach reach friction, and approach contraction losses.
C     Do the same in estimating the contraction coefficient.  The
C     change in bottom elevation to the approximate vena contracta
C     location is included because it was likely present in the
C     model studies.
 
      DZVC = ZB2 - ZBVEC(IAT3D)
      CALL FNDCC
     I          (STDOUT, CD, DUP, ADRXS2, A2, HDRAT, HEAD, DZVC,
     O           CC)
 
      YVC = DUP*YOVERD(NSEC(IAT3D), DUP, A2, CC)
 
      HEADVC = HUP + HDATUM - ZBVEC(IAT3D)
C      WRITE(STDOUT,*) ' YVC=',YVC,' YRATIO=',YVC/DUP,' QRF=',QRF
C      WRITE(STDOUT,*) ' HEADVC=',HEADVC
      
      AVC = CC*A2
      DEN = 1.0 + (AVC)**2*(APPLEN*GRAV2/(K1*K2) - (ALP1 - APPLOS)/
     A                                                   A1**2)
C      WRITE(STDOUT,*) 'AVC=',AVC,' APPLEN=',APPLEN,' GRAV2=',GRAV2,
C     A             ' K1=',K1,' K2=',K2,' ALP1=',ALP1,' APPLOS=',APPLOS,
C     B             ' A1=',A1
      NUM = GRAV2*(HEADVC - YVC) + (ALP1 - APPLOS)*(QRF/A1)**2
C      WRITE(STDOUT,*) ' NUM=',NUM, ' DEN=',DEN,' CC=',CC
 
      
C     Starting estimate ignores the cross product terms. This estimate
C     is correct if QRF = 0.
      Q = AVC*SQRT(NUM/DEN)
      IF(QRF.GT.0.0) THEN
        KNT = 0
        NUMFAC = QRF*((ALP1 - APPLOS)/A1**2 - GRAV2*APPLEN/(K1*K2))
 100    CONTINUE
          QT = AVC*SQRT((NUM + Q*NUMFAC)/DEN)
          IF(ABS(Q - QT)/QT.GE.EPSF) THEN
            Q = QT
            KNT = KNT + 1
            IF(KNT.GT.100) THEN
              WRITE(STDOUT,54)
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            GOTO 100
          ENDIF
          Q = QT
      ENDIF
 
C      WRITE(STDOUT,*) ' FNDQ5: TYPE 5 FLOW=',Q,' with CC=',CC
      Q5 = Q
      RETURN
      END
C
C
C
      SUBROUTINE   FRFT5
     I                  (STDOUT, HDATUM, HUP, IU, ID, CULCLS,
     O                   TYPE, EXPFLG, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute the free flow and the drop to free flow for
C     type 5 flow.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EXPFLG, ID, IU, STDOUT, TYPE
      REAL FREED, HDATUM, HUP, QFREE
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     TYPE   - Culvert flow type
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'culcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL AVC, BETA, DDN, DEN, DUP, DZVC, HDRAT, HEAD, HEADVC, NUM,
     A     NUMFAC, Q, QC, QC1, QC2, QC3, QRF, QT, TP, TY6LSS, Y3LIM, YC,
     B     YVC, Z1TRUE, Z43T
      CHARACTER JMPLOC*32
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL YOVERD
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DPM26, FNDCC, FNDCD5, FNDCDE, LKTJ, LKTQC, SUPSUB,
     A         XLKT22, YOVERD
 
C     + + + OUTPUT FORMATS + + +
 54   FORMAT(/,' *BUG:XXX* No convergence on Q in subroutine FRFT5.')
 56   FORMAT(/,' *WRN:583* BETAF < 0.0 for type 5 flow.')
 58   FORMAT(/,' Parameters for transition from type=',I3,
     A       ' to type 4:',/,
     B       4X,'Coefficient of discharge=',F8.3,/,
     C       4X,'Piezometric depth at section 3=',F8.3,/,
     D       4X,'Momentum-flux coefficient=',F8.3,/,
     E       4x,'Kinetic-energy flux coefficient=',F8.3)
 60   FORMAT(/,' Jump location:',A)
 62   FORMAT('  Note: With the jump in the barrel exit, the depth at',
     A      ' section 3',/,2X,'given in the summary tables is the ',
     B      ' depth upstream of the jump.',/,2X,'The value for the',
     C      ' head at section 3 & 43 is the piezometric',/,2X,
     D      'level in the exit that causes full flow.  These water',
     E      ' levels differ.')
C***********************************************************************
C      WRITE(STDOUT,*) ' FRFT5: Z1TY51=',Z1TY51,' Z1TY52=',Z1TY52,
C     A                ' Z1TY5=',Z1TY5

      DUP = DVEC(IU)
      DDN = DVEC(ID)
      CALL XLKT22
     I           (ADRXS2,
     M            DUP,
     O            A2, T2, DT2, J2, K2, DK2, BET2, DBET2, ALP2, DALP2,
     O            QC2)
      Z1TRUE = HDATUM + HUP
C      WRITE(STDOUT,*) ' FRFT5: Z1TRUE=',Z1TRUE
      Y1 = Z1TRUE - ZB1
      CALL XLKT22
     I           (ADRXS1,
     M            Y1,
     O            A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1,
     O            QC1)
      HEAD =  HUP + HDATUM - ZB2
 
      IF(Z1TRUE.GE.Z1TY5) THEN
        TYPE = 5
        LFTYPE = 5
        LSTYPE = -1
        BETAF = 0.0
        BETA3 = 0.0
        ALPHA3 = 0.0
        Y3LIM = Y3LTY5
        HDRAT = HEAD/DUP
C        WRITE(STDOUT,*) ' FRFT5: HDRAT=',HDRAT
C       Find the discharge coefficient for type 5 flow.
 
        CALL FNDCD5
     I             (STDOUT, CULCLS, HDRAT, RBVAL,
     O              CD)
 
C        WRITE(STDOUT,*) ' FRFT5: CD=',CD
 
 
C        Deduce the contraction coefficient from the CD given by
C       the USGS tables.  These appear to ignore approach velocity,
C       approach reach friction, and approach contraction losses.
C       Do the same in estimating the contraction coefficient.  The
C       change in bottom elevation to the approximate vena contracta
C       location is included because it was likely present in the
C       model studies.
 
        DZVC = ZB2 - ZBVEC(IAT3D)
        CALL FNDCC
     I            (STDOUT, CD, DUP, ADRXS2, A2, HDRAT, HEAD, DZVC,
     O             CC)
        Y3LIM = Y3LTY5
      ELSE
        IF(Z1TY52.GT.Z1TY51) THEN
          TYPE = 52
          LFTYPE = 52
          LSTYPE = -1
          IF(Z1TRUE.GE.Z1TY52.AND.Z1TRUE.LT.Z1TY5) THEN
            IF(CC2T5.LE.0.0) THEN
              WRITE(STDOUT,*) ' *BUG:XXX FRFT5: CC2T5=',CC2T5
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            IF(CC5.LE.0.0) THEN
              WRITE(STDOUT,*) ' *BUG:XXX FRFT5: CC5=',CC5
              STOP 'Abnormal stop. Errors found.'
            ENDIF
C           No discharge coefficient used.   Compute the contraction
C           coefficient by linear interpolation.
            CC =  CC2T5 + (Z1TRUE - Z1TY52)*(CC5 - CC2T5)
     A                                           /(Z1TY5 - Z1TY52)
            BETAF =  BT2T5 + (Z1TRUE - Z1TY52)*(BT3AT5 - BT2T5)
     A                                           /(Z1TY5 - Z1TY52)
            ALPHAF =  AP2T5 + (Z1TRUE - Z1TY52)*(AP3AT5 - AP2T5)
     A                                           /(Z1TY5 - Z1TY52)
            Y3LIM = Y3LTY2 + (Z1TRUE - Z1TY52)*(Y3LTY5 - Y3LTY2)
     A                                           /(Z1TY5 - Z1TY52)
          ELSE
            WRITE(STDOUT,*) ' *BUG:XXX* WRONG ELEVATION FOR TYPE 52',
     A                ' IN FRFT5.'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
        ELSE
          TYPE = 51
          LFTYPE = 51
          LSTYPE = -1
C          WRITE(STDOUT,*) ' Z1TY51=',Z1TY51,' Z1TY5=',Z1TY51
C          WRITE(STDOUT,*) ' Z1TRUE=',Z1TRUE

          IF(Z1TRUE.GE.Z1TY51.AND.Z1TRUE.LT.Z1TY5) THEN
            IF(CC1T5.LE.0.0) THEN
              WRITE(STDOUT,*) ' *BUG:XXX FRFT5: CC1T5=',CC1T5
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            IF(CC5.LE.0.0) THEN
              WRITE(STDOUT,*) ' *BUG:XXX FRFT5: CC5=',CC5
              STOP 'Abnormal stop. Errors found.'
            ENDIF
C           No discharge coefficient used.   Compute the contraction
C           coefficient by linear interpolation.
            CC =  CC1T5 + (Z1TRUE - Z1TY51)*(CC5 - CC1T5)
     A                                           /(Z1TY5 - Z1TY51)
            BETAF =  BT1T5 + (Z1TRUE - Z1TY51)*(BT3AT5 - BT1T5)
     A                                           /(Z1TY5 - Z1TY51)
            ALPHAF =  AP1T5 + (Z1TRUE - Z1TY51)*(AP3AT5 - AP1T5)
     A                                           /(Z1TY5 - Z1TY51)
            Y3LIM = Y3LTY1 + (Z1TRUE - Z1TY51)*(Y3LTY5 - Y3LTY1)
     A                                           /(Z1TY5 - Z1TY51)
C            WRITE(STDOUT,*) ' FRFT5: BT1T5=',BT1T5,' BT3AT5=',BT3AT5,
C     A            ' Z1TY5=',Z1TY5,' Z1TY51=',Z1TY51, ' Z1TRUE=',Z1TRUE,
C     B            ' Y3LTY5=',Y3LTY5,' Y3LTY1=',Y3LTY1,' Y3LIM=',Y3LIM
          ELSE
            WRITE(STDOUT,*) ' *BUG:XXX* WRONG ELEVATION FOR TYPE 51',
     A                ' IN FRFT5.'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
        ENDIF
      ENDIF
C     Now estimate the flow using the contraction coefficient as
C     the estimate of depth at the vena contracta.
 
      YVC = DUP*YOVERD(NSEC(IAT3D), DUP, A2, CC)
 
 
      HEADVC = HUP + HDATUM - ZBVEC(IAT3D)
      QRF = WFRDF
C      WRITE(STDOUT,*) ' YVC=',YVC,' YRATIO=',YVC/DUP,' QRF=',QRF
      AVC = CC*A2
      DEN = 1.0 + (AVC)**2*(APPLEN*GRAV2/(K1*K2) - (ALP1 - APPLOS)/
     A                                                   A1**2)
      NUM = GRAV2*(HEADVC - YVC) + (ALP1 - APPLOS)*(QRF/A1)**2
 
C     Starting estimate  ignores the cross product terms. This estimate
C     is correct if QRF = 0.
      Q = AVC*SQRT(NUM/DEN)
      IF(QRF.GT.0.0) THEN
        KNT = 0
        NUMFAC = QRF*((ALP1 - APPLOS)/A1**2 - GRAV2*APPLEN/(K1*K2))
 100    CONTINUE
          QT = AVC*SQRT((NUM + Q*NUMFAC)/DEN)
          IF(ABS(Q - QT)/QT.GE.EPSF) THEN
            Q = QT
            KNT = KNT + 1
            IF(KNT.GT.100) THEN
              WRITE(STDOUT,54)
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            GOTO 100
          ENDIF
          Q = QT
      ENDIF
 
      QFREE = Q
 
C      WRITE(STDOUT,*) ' FRFT5: TYPE 5 FLOW=',Q,' with CC=',CC
 
      Q2 = Q
      Q1 = Q2 + QRF
      Q3 = Q
 
      TY6LSS = ((1.0/C46**2 - 1.0))*(Q/A2)**2/GRAV2
C      WRITE(STDOUT,*) ' TYPE 6 LOSS =',TY6LSS
 
C      EUP = ZB1 + Y1 + ALP1*(Q/A1)**2/GRAV2
 
C     Section 2 is always full.  Therefore Y2 is always DUP.
C     The piezometric pressure is atmospheric at the water surface
C     just downstream of the entrance soffit.  This assumes perfect
C     venilation.
 
      Z2 = DUP + ZB2
      Y2 = DUP
C      WRITE(STDOUT,*) ' TOTAL ENERGY ELEV AT SEC. 1=',EUP
      CALL LKTQC
     I          (NSEC(IAT3D),
     M           YVC,
     O           QC)
C      WRITE(STDOUT,*) ' FRFT5: DEPTH AT VC=',YVC,' FROUDE=',Q/QC
C      VC = Q/(CC*A2)
C      WRITE(STDOUT,*) ' VELOCITY AT VC =',VC
C     Find critical depth at exit.
      Y3 = 0.5*DDN
      CALL FNDCDE
     I           (STDOUT, ADRXS3, Q,
     M            Y3)
      YC = Y3
 
C
C     Find the depth at the exit from the culvert barrel.  A hydraulic
C     jump  may exist in the barrel, part way in the barrel, or
C     in some problem cases not at all.
 
      CALL SUPSUB
     I           (STDOUT, IAT3D, ID, YVC, Y3LIM, Q, DUP, YC,
     M            TY6LSS,
     O            Z3, Z43T, JMPLOC)
 
C     Compute values at section 3.  Not done in SUPSUB
      Y3 = Z3 - ZB3
      Z3P = Z43T
      TP = MIN(Y3,DDN)
      CALL XLKT22
     I           (ADRXS3,
     M            TP,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3,
     O            QC3)
      CALL LKTJ
     I         (ADRXS3,
     M          Y3,
     O          J3)
C     Type 5 flow at submergence limit passes to full flow
C     immediately upon submergence.  Y3FREE is used to remember the
C     area defining depth that is to be associated with the
C     submergence limit tailwater at section 43.  Thus for type 5
C     flow the proper value for Y3FREE is DDN.
      Y3FREE = DDN
      Q3FREE = Q3
      Z43T = Y3LIM + ZB3
C      WRITE(STDOUT,*) ' FRFT5: Y3LIM=',Y3LIM,' Z43T=',Z43T,
C     A                ' BETAF=',BETAF
      CALL DPM26
     I          (STDOUT, WFRDF, MFRDF, Z43T,
     O           EXPFLG)
      WRITE(STDOUT,60) JMPLOC
      IF(JMPLOC.EQ.'Jump in the barrel exit') THEN
        WRITE(STDOUT,62)
      ENDIF
C      WRITE(STDOUT,*) ' Z3=',Z3,' Z43=',Z43,' Z4=',Z4,' Z44=',Z44
      IF(EXPFLG.EQ.1) THEN
C       Expansion of flow in departure reach.
        FREED = Z1TRUE - Z4
C        WRITE(STDOUT,*) ' FREED=',FREED
C       Compute the discharge coefficient that is needed for
C       full flow to match the type 5, 51, or 52 flow.  Tailwater
C       is at Z43OLD, set in DPM26.
 
        TP = ((ALP1 - APPLOS)*(Q1/A1)**2/GRAV2 + Z1 -
     A      ((Q3/A3FULL)**2/GRAV2 + Z43OLD) - L23*(Q2/K2FULL)**2
     B     - APPLEN*(Q1*Q2)/(K1*K2FULL))*(GRAV2*(A2FULL/Q2)**2)
        CD5T4 = SQRT(1.0/(1.0 + TP))
        CDF = CD5T4
        Y3PF =  Z43OLD - ZB3
      ELSE
        FREED = 0.0
        Y3PF =  Z43OLD - ZB3
      ENDIF
 
 
 
C     Now make sure that the M43 as computed in DPM26 for type 5 will
C     be closely approximated in FRFT7.  In DPM26 the momentum flux
C     was at part full and super critical.  Also the pressure force
C     terms include the force at the culvert exit headwall at some
C     point in the jump.   That is, Z3 and Z43 may differ even if
C     the barrel is part full.
 
      IF(BETAF.GT.0.0) THEN
        BETA = BETAF
      ELSE
        BETA = 1.0
      ENDIF
C      WRITE(STDOUT,*) ' A3FULL=',A3FULL,' BETA=',BETA,' Q3=',Q3,
C     A         ' A3=',A3,' GRAV=',GRAV,' J3=',J3,' J3Z43=',J3Z43
 
      BETAF = A3FULL*(BETA*Q3**2/A3 + GRAV*(J3 - J3Z43))/Q3**2
      IF(BETAF.LT.0.0) THEN
        WRITE(STDOUT,56)
        BETAF = 0.0
      ENDIF
C     No direct way to estimate ALPHAF.  Make a rough estimate
      IF(BETAF.GT.1.0) THEN
        ALPHAF = 1.0 + 3.0*(BETAF - 1.0)
      ELSE
        ALPHAF = 1.0
      ENDIF
      WRITE(STDOUT,58) TYPE, CD5T4, Y3PF, BETAF, ALPHAF
 
      RETURN
      END
C
C
C
      SUBROUTINE   DOTY6
     I                  (STDOUT, CULCLS, A1TRUE, ALP1T, K1TRUE, Z1TRUE,
     I                   A2FULL, K2FULL, A3FULL, DDN, ZB, CDIS,
     O                   Z3P, Q, ZAT2)
 
C     + + + PURPOSE + + +
C     Compute type 6 flow in the culvert.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      REAL A1TRUE, A2FULL, A3FULL, ALP1T, CDIS, DDN, K1TRUE, K2FULL, Q,
     A     Z1TRUE, Z3P, ZAT2, ZB
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     A1TRUE - area at section 1 for the known water surface elevation
C     ALP1T  - value of energy flux correction coefficient at section 1
C     K1TRUE - conveyance at section 1 for the known elevation there
C     Z1TRUE - known elevation at section 1
C     A2FULL - area at section 2 with full flow in barrel
C     K2FULL - full barrel conveyance at section 2
C     A3FULL - area at section 3 with full flow in barrel
C     DDN    - vertical diameter of culvert barrel at downstream end
C     ZB     - bottom elevation at section 3
C     CDIS   - discharge coefficient
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
      REAL CCON, DIV, NUM, QT, Y3P
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL TY6RAT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL TY6RAT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* No convergence on Q in subroutine DOTY6.')
C***********************************************************************
C     Make the first estimate as if roadflow were zero.  If it is
C     zero, the result is valid.  Otherwise, iterate
C     to include the effect of velocity head of approach from the
C     flow over the road.
C     Estimate piezometric level at section 3 to start the iterations.
 
      Z3P = ZB + 0.75*DDN
      NUM = GRAV2*(Z1TRUE - Z3P)
C      IF(FRCFAC.EQ.0.0) THEN
C       Prismatic barrel.
C        DIV = (1.0 + GRAV2*(CDIS*A2FULL)**2*(APPLEN/(K1TRUE*K2FULL)
C     B       + L23/K2FULL**2 - (ALP1T - APPLOS)/(GRAV2*A1TRUE**2)))
C      ELSE
        DIV = (1.0  + CDIS**2*((A2FULL/A3FULL)**2 - 1.0 +
     B    GRAV2*A2FULL**2*(APPLEN/(K1TRUE*K2FULL)
     C     + FRCFAC - (ALP1T - APPLOS)/(GRAV2*A1TRUE**2))))
C      ENDIF
      KNT = 0
 100  CONTINUE
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          WRITE(STDOUT,50)
          STOP 'Abnormal stop. Errors found.'
        ENDIF
 
 
        Q = CDIS*A2FULL*SQRT(NUM/DIV)
 
C       Adjust the piezometric level at 3 for the latest flow.
        Y3P = DDN*TY6RAT(Q, CULCLS, DDN, A3FULL)
        Z3P = ZB + Y3P
        NUM = GRAV2*(Z1TRUE - Z3P)
 
        IF(WFRDF.GT.0.0) THEN
C         Add terms to the numerator to represent values affected
C         by flow over the road.
          NUM = NUM + WFRDF*((ALP1T - APPLOS)*
     A                         (WFRDF + Q + Q)/(GRAV2*A1TRUE**2)
     B                          - APPLEN*Q/(K1TRUE*K2FULL))
        ENDIF
        QT = CDIS*A2FULL*SQRT(NUM/DIV)
 
C        WRITE(STDOUT,*) ' DOTY6: KNT=',KNT,' Q=',Q,' QT=',QT
        IF(ABS(QT - Q)/Q.GT.EPSF) THEN
          Q = QT
          GOTO 100
        ENDIF
      Q = QT
C     Compute estimate of the piezometric level at section 2.   Estimate
C     the coefficient of contraction from the discharge coefficient.
C     Includes effect of vena contracta near entrance.
 
      CCON = 1.0/(SQRT(1.0/CDIS**2 -1.0) + 1.0)
      ZAT2 = ALP1T*((Q + WFRDF)/A1TRUE)**2/GRAV2 + Z1TRUE -
     A                (Q/(CCON*A2FULL))**2/GRAV2
C      WRITE(STDOUT,*) ' DOTY6: ZAT2=',ZAT2,' CCON=',CCON
      RETURN
      END
C
C
C
      SUBROUTINE   DOTY61
     I                   (STDOUT, A1TRUE, ALP1T, K1TRUE, Z1TRUE, IU, ID,
     I                    ZB, CDIS, AVH,
     M                    ZP,
     O                    Q, ZAT2)
 
C     + + + PURPOSE + + +
C     Compute type 61 flow in the culvert. On entry ZP has
C     estimated part-full piezometric level at section 3.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, IU, STDOUT
      REAL A1TRUE, ALP1T, AVH, CDIS, K1TRUE, Q, Z1TRUE, ZAT2, ZB, ZP
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     A1TRUE - area at section 1 for the known water surface elevation
C     ALP1T  - value of energy flux correction coefficient at section 1
C     K1TRUE - conveyance at section 1 for the known elevation there
C     Z1TRUE - known elevation at section 1
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     ZB     - bottom elevation at section 3
C     CDIS   - discharge coefficient
C     AVH    - area for computing velocity head
C     ZP     - elevation of the piezometric surface at section 3
C     Q      - Flowrate
C     ZAT2   - elevation of water surface at section 2
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'appcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IS, KNT, SFLAG
      REAL AL, AR, CCON, DH, DIV, KL, KR, NUM, QT, ZT
      DOUBLE PRECISION DX, SUM
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FNDCDE, LKTA, LKTQC, SFPSBE
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* No convergence on Q in subroutine DOTY61.')
C***********************************************************************
      Y3 = ZP - ZB
      CALL LKTQC
     I          (ADRXS3,
     M           Y3,
     O           Q)
      KNT = 0
 100  CONTINUE
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          WRITE(STDOUT,50)
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        CALL LKTA
     I           (ADRXS3,
     M            Y3,
     O            A3)
 
        DH = (1.0/CDIS**2 - 1.0)*(Q/AVH)**2/GRAV2
C        WRITE(STDOUT,*) ' DOTY61: DH=',DH,' CDIS=',CDIS,' AVH=',AVH,
C     A                  ' ZP=',ZP
        CALL SFPSBE
     I             (STDOUT, IU, ID, DH, Q, ZP,
     O              IS, ZT, SFLAG)
C        IF(ZT.LT.ZB2 + DUP) THEN
C          WRITE(STDOUT,*) ' PROBLEM IN DOTY61'
C          WRITE(STDOUT,*) ' Invalid values from SFPSBE.  ZT < ZB2 + DUP'
C          WRITE(STDOUT,*) ' ZT=',ZT,' ZB2 + DUP=',ZB2 + DUP
C        ENDIF
        SUM = 0.D0
        AL = AVEC(IU)
        KL = KVEC(IU)
        DO 200 I=IU+1, ID
          AR = AVEC(I)
          KR = KVEC(I)
          DX = ABS(XVEC(I) - XVEC(I-1))
          IF(AL.LE.AR) THEN
C           The flow is expanding.
            SUM = SUM + DX/(KL*KR) + KD(I)*(1.0/AL**2 - 1.0/AR**2)/
     A                                                         GRAV2
          ELSE
C           The flow is contracting.
            SUM = SUM + DX/(KL*KR) + KA(I)*(1.0/AR**2 - 1.0/AL**2)/
     A                                                         GRAV2
          ENDIF
          KL = KR
          AL = AR
 200    CONTINUE
        NUM = GRAV2*(Z1TRUE - ZP)
        DIV = (1.0  + CDIS**2*((AVH/A3)**2 - 1.0 +
     B    GRAV2*AVH**2*(APPLEN/(K1TRUE*K2FULL)
     C     + SUM - (ALP1T - APPLOS)/(GRAV2*A1TRUE**2))))
C        WRITE(STDOUT,*) ' DIV=',DIV,' SUM=',FRCFAC,' ZT=',ZT
C        WRITE(STDOUT,*) 'DOTY61: AVH=',AVH,' NUM=',NUM
C        WRITE(STDOUT,*) ' Z1TRUE=',Z1TRUE,' ZP=',ZP,' WFRDF=',WFRDF,
C     A       ' A1TRUE=',A1TRUE,' K1TRUE=',K1TRUE,' K2FULL=',K2FULL,
C     B       ' ALP1T=',ALP1T,' APPLOS=',APPLOS,' GRAV2=',GRAV2,
C     C       ' APPLEN=',APPLEN
C        WRITE(STDOUT,*) ' A3=',A3
        IF(WFRDF.GT.0.0) THEN
C         Add terms to the numerator to represent values affected
C         by flow over the road.
          NUM = NUM + WFRDF*((ALP1T - APPLOS)*
     A                         (WFRDF + Q + Q)/(GRAV2*A1TRUE**2)
     B                          - APPLEN*Q/(K1TRUE*K2FULL))
        ENDIF
        QT = CDIS*AVH*SQRT(NUM/DIV)
 
C       Adjust the piezometric level at 3 for the latest flow.
        CALL FNDCDE
     I             (STDOUT, ADRXS3, QT,
     M              Y3)
 
        ZP = ZB + Y3
 
C        WRITE(STDOUT,*) ' DOTY61: KNT=',KNT,' Q=',Q,' QT=',QT
        IF(ABS(QT - Q)/Q.GT.EPSF) THEN
          Q = QT
          GOTO 100
        ENDIF
      Q = QT
C     Compute estimate of the piezometric level at section 2.   Estimate
C     the coefficient of contraction from the discharge coefficient.
C     Includes effect of vena contracta near entrance.
 
      CCON = 1.0/(SQRT(1.0/CDIS**2 -1.0) + 1.0)
      ZAT2 = ALP1T*((Q + WFRDF)/A1TRUE)**2/GRAV2 + Z1TRUE -
     A                (Q/(CCON*A2FULL))**2/GRAV2
C      WRITE(STDOUT,*) ' DOTY61: ZAT2=',ZAT2,' CCON=',CCON
      RETURN
      END
C
C
C
      SUBROUTINE   F6BDY
     I                  (STDOUT, HDATUM, HUP, IU, ID, CULCLS,
     O                   EFLAG, ZP, Q6, ZAT2)
 
C     + + + PURPOSE + + +
C     Find the flow and piezometric level at the type 6 flow boundary.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, ID, IU, STDOUT
      REAL HDATUM, HUP, Q6, ZAT2, ZP
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     ZP     - elevation of the piezometric surface at section 3
C     Q6     - free flow of type 6
C     ZAT2   - elevation of water surface at section 2
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'culcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL DDN, DUP, FDRDW, Z, ZSBRDF
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DOTY6, GETFRF, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:689* Flow type 6 limit not defined. Approach',
     A  ' and/or',/,11X,'departure elevation unrealistic.')
C***********************************************************************
C     The type 6 flow boundary is set at an approach head ratio
C     of 1.5.  We need to know the piezometric level at the outlet
C     for some of the transitions.
      DUP = DVEC(IU)
      DDN = DVEC(ID)
      IF(ZB2.GE.ZB3) THEN
C       Barrel slope is >= 0.
        Z = ZB2 + 1.5*DUP
      ELSE
C       Barrel slope is adverse, that is, < 0.
        Z = ZB3 + 1.5*DDN
      ENDIF
      IF(Z.LE.HDATUM) THEN
        WRITE(STDOUT,50)
        EFLAG = 1
        RETURN
      ENDIF
      CALL GETFRF
     I           (Z,
     O            ZSBRDF, FDRDW)
      Y1L = Z - ZB1
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1L,
     O            A1L, T1L, DT1L, J1L, K1L, DK1L, BET1L, DBET1L, ALP1L,
     O            DALP1L)
      CD6 = C46
      CALL DOTY6
     I          (STDOUT, CULCLS, A1L, ALP1L, K1L, Z, A2FULL, K2FULL,
     I           A3FULL, DDN, ZB3, CD6,
     O           ZP, Q6, ZAT2)
 
      Z1TY6 = Z
 
C     Compute the free flow values for the current upstream
C     level.
      Z = HUP + HDATUM
      CALL GETFRF
     I           (Z,
     O            ZSBRDF, FDRDW)
 
      RETURN
      END
C
C
C
      SUBROUTINE   F61BDY
     I                   (STDOUT, TYPE, HDATUM, ZDATUM, HUP, IU, ID,
     I                    CULCLS, Y2VAL,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Find the elevation at section 1 such that the computed water
C     level at section 2 matches the value, Y2VAL.  The flow is
C     always part full and at critical depth at section 3.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, ID, IU, STDOUT, TYPE
      REAL HDATUM, HUP, Y2VAL, ZDATUM
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     TYPE   - culvert flow type
C     HDATUM - Datum for measuring head
C     ZDATUM - datum for local elevation
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     Y2VAL  - value of depth at section 2
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'epscom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'rty6c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, KNT, KNT2, MESS_OUT
      REAL DDN, DUP, F, FDRDW, FHIGH, FLOW, Y, YMAX, Z, ZHIGH, ZLOW,
     A     ZMIN, ZSBRDF, ZT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FMXARG, RTY61
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FMXARG, GETFRF, RGF, RTY61
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* Maximum argument=',F10.3,' reached and no',
     A       ' sign change in F61BDY')
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF CLAIMS NONE IN',
     A       ' F61BDY.')
 60   FORMAT(' *BUG:XXX* RGF: MORE THAN 100 ITERATIONS IN F61BDY.')
C***********************************************************************
C     Allow five residual problem mesages
      MESS_OUT = 5

C     Clear the error flag
      EFLAG = 0
      IUP = IU
      IDN = ID
      OUTUN = STDOUT
      EF = 0
      FTYPE = TYPE
      CLASS = CULCLS
      YAT2 = Y2VAL
      ZDAT = ZDATUM
      DDN = DVEC(IDN)
      DUP = DVEC(IUP)
      ZMIN = YAT2 + ZB2
C      WRITE(STDOUT,*) ' F61BDY: ZMIN=',ZMIN
C     Search for a change in sign of the residual.
      IF(TYPE.EQ.2) THEN
        IF(CD1.GT.0.0) THEN
          Z = Z1TY1 + 0.01*DUP
          ZMIN = Z1TY1
C          WRITE(STDOUT,*)' F61BDY: Z SET BY Z1TY1=',Z
          IF(Z.LT.ZB2 + 1.1*DUP) THEN
            Z = ZB2 + ZB2 + 1.1*DUP
C            WRITE(STDOUT,*) ' F61BDY: Z RESET TO:',Z
          ENDIF
        ELSE
          IF(ZB2.LT.ZB3) THEN
C           Adverse barrel slope
            IF(ZB3.GE.ZB2+DUP) THEN
C             Type 2 flow does not exist.  Zero head has entrance
C             soffit at or below water surface.
              EFLAG = -1
              RETURN
            ELSE
              Z = ZB2 + DUP + 0.1*(ZB2 + DUP - ZB3)
            ENDIF
          ELSE
            Z = 1.2*DUP +  ZB2
          ENDIF
        ENDIF
        Y = Z - ZB1
        IF(Y.LT.0.0) THEN
          Y = 0.1*DUP
          Z = ZB1 + Y
          IF(Z.GT.1.5*DUP + ZB2) THEN
C           Type 2 limit does not exist.
            EFLAG = -1
            RETURN
          ENDIF
        ENDIF
      ELSE
        Z = 1.2*DDN + ZB3
        Y = Z - ZB1
        IF(Y.LT.0.0) THEN
          Y = 0.1*DDN
          Z = ZB1 + Y
          IF(Z.GT.1.5*DDN + ZB3) THEN
C           Type 61 limit does not exist.
            EFLAG = -1
            RETURN
          ENDIF
        ENDIF
      ENDIF
 
      YMAX = FMXARG(ADRXS1)
      ZLOW = 0.0
      ZHIGH = 0.0
      KNT = 0
      KNT2 = 0

 100  CONTINUE
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          WRITE(OUTUN,*) ' No sign change in F61BDY after 100 tries.'
          EFLAG = -1
          GOTO 999
        ENDIF
        F = RTY61(Z)
C        WRITE(OUTUN,*) ' F61BDY: Z=',Z,' Y1=',Z-ZB1,' F=',F
C        WRITE(OUTUN,*) ' ZMIN=',ZMIN,' SBFLAG=',SBFLAG
 
        IF(SBFLAG.EQ.-1) THEN
C         Special action required.  Steady flow profile from critical
C         depth failed to complete or sqrt of negative number.
C         Increase Z and try again.
          KNT2 = KNT2 + 1
          IF(KNT2.GT.100) THEN
C           Apparently type 2 limit does not exist.
            EFLAG = -1
            GOTO 999
          ENDIF
          IF(MESS_OUT.GT.0) THEN
            WRITE(OUTUN,*) ' RTY61: Residual problem. Z=',Z
            MESS_OUT = MESS_OUT - 1
          ENDIF
C         Upstream elevation and flow are related.  If the
C         profile failed at this elevation it will fail at any
C         lower elevation.  Upgrade the minimum elevation.
          IF(Z.GT.ZMIN) ZMIN = Z
          Z = Z + 0.01*DUP
          GOTO 100
        ELSEIF(SBFLAG.EQ.-2) THEN
C         Assume solution does not exist.
          EFLAG = -1
          GOTO 999
        ENDIF
        IF(ABS(F).GT.EPSF) THEN
C         Continue looking for change in sign of residual.
          IF(F.GT.0.0) THEN
            ZHIGH = Z
            FHIGH = F
            IF(ZLOW.EQ.0.0) THEN
              Y = Z - ZB1
              Y = 0.99*Y
              IF(Y.LT.EPSABS) THEN
C               Negative residual not found.  Take to mean that
C               the limit does not exist.
                EFLAG = -1
                GOTO 999
              ENDIF
              ZT = ZB1 + Y
              IF(ZT.LE.ZMIN) THEN
                Z = 0.5*(ZMIN + Z)
                IF((Z-ZMIN)/ZMIN.LT.1.E-7) THEN
C                 Conclude that the limit does not exist.
                  EFLAG = -1
                  GOTO 999
                ENDIF
              ELSE
                Z = ZT
              ENDIF
              GOTO 100
            ENDIF
          ELSE
            ZLOW = Z
            FLOW = F
            IF(ZHIGH.EQ.0.0) THEN
              Y = Z - ZB1
              Y = 1.05*Y
              IF(Y.GT.YMAX) THEN
C               Should not happen but flag just in case.
                WRITE(OUTUN,50) YMAX
                STOP 'Abnormal stop.  Errors found.' 
              ENDIF
              Z = ZB1 + Y
              GOTO 100
            ENDIF
          ENDIF
C         Sign change here.
 
C          WRITE(STDOUT,*) ' F61BDY: ZLOW=',ZLOW,' FLOW=',FLOW,
C     A                   ' ZHIGH=',ZHIGH,' FHIGH=',FHIGH
          CALL RGF
     I            (1.E-6, EPSF, RTY61,
     M             ZLOW, ZHIGH, FLOW, FHIGH,
     O             Z, FLAG)
          IF(FLAG.EQ.1) THEN
            WRITE(STDOUT, 54)
            STOP 'Abnormal stop.  Errors found.' 
          ELSEIF(FLAG.EQ.2) THEN
            WRITE(STDOUT,60)
            STOP 'Abnormal stop.  Errors found.' 
          ENDIF
        ENDIF
 
      IF(TYPE.EQ.2) THEN
        Z1TY2 = Z
        CD2 = CD
      ELSE
        Z1TY61 = Z
        CD61 = CD
      ENDIF
 
 999  CONTINUE
C     Compute the free flow values for the current upstream
C     level.
      Z = HUP + HDATUM
      CALL GETFRF
     I           (Z,
     O            ZSBRDF, FDRDW)
 
      RETURN
      END
C
C
C
      SUBROUTINE   FRFT6
     I                  (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O                   EFLAG, TYPE, EXPFLG, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute free flow of type 6 and its relatives.
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL FREED, HDATUM, HUP, QFREE, ZDATUM
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     ZDATUM - datum for local elevation
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TYPE   - Culvert flow type
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'culcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IFLAG
      REAL DDN, DUP, Q6, TP, TP2, YAT2, ZAT2
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN, SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DOTY6, DOTY61, DPM26, F61BDY, F6BDY, FULBAR, LKTJ,
     A         XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 76   FORMAT(/,' *WRN:574* Unable to force type 6 Cd to match type 61',
     A         ' flow at its limit.',/,11X,' Using type 6 Cd. ',
     B     'Manual adjustment of 2-D table may be needed.')
 78   FORMAT(/,' Type 6 limits:',/,
     A        4X,'Depth at section 1=',F8.3,/,
     B        4X,'Piezometric depth at section 3=',F8.3)
C***********************************************************************
      DDN = DVEC(ID)
      DUP = DVEC(IU) 
      A2 = A2FULL
      IF(CD6.EQ.0.0) THEN
C       The limit for type 6 flow has not been computed.  Find it.
        CALL F6BDY
     I            (STDOUT, HDATUM, HUP, IU, ID, CULCLS,
     O             EFLAG, Z3PTY6, Q6, ZAT2)
        WRITE(STDOUT,78) Z1TY6 - ZB1, Z3PTY6 - ZB3
      ENDIF
 
      IF(CD61.EQ.0.0) THEN
C       Type 61 limit has not been computed.  See if limit exists
C       and if it does, find it.
        IF(ZB2.GE.ZB3) THEN
C         Type 61 flow does not exist.  Therefore limit does not exist.
          CD61 = -1.0
        ELSE
C         Type 2 flow exists and the barrel slope is adverse.
C         Therefore, type 61 exists!
          YAT2 = ZB3 + DDN - ZB2
          CALL F61BDY
     I               (STDOUT, 61, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     I                YAT2,
     O                IFLAG)
          IF(IFLAG.EQ.-1) THEN
            WRITE(STDOUT,*) ' Type 61 limit not found when it must',
     A                  ' exist.'
            STOP 'Abnormal stop. Errors found.'
          ELSE
            CD61 = CD
            AVH61 = AVH
C           Type 61 exists.  Therefore transition is between 61 and 6.
 
            TP = ((ALP1L - APPLOS)*(Q1L/A1L)**2/GRAV2 + Z1L -
     A          ((Q3L/A3FULL)**2/GRAV2 + Z3PART) - L23*(Q2L/K2FULL)**2
     B         - APPLEN*(Q1L*Q2L)/(K1L*K2FULL))*(GRAV2*(A2FULL/Q2L)**2)
            IF(TP.LT.0.0) THEN
              WRITE(STDOUT,76)
              CD61T6 = C46
            ELSE
              CD61T6 = SQRT(1.0/(1.0 + TP))
            ENDIF
            BT61T6 = A3FULL/A3PART
            AP61T6 = (BT61T6)**2
        WRITE(STDOUT,*) ' CD61T6=',CD61T6,' BT61T6=',BT61T6,' AP61T6=',
     A                   AP61T6
          WRITE(STDOUT,*) ' Z1TY61=',Z1TY61,' CD61=',CD61
          ENDIF
        ENDIF
      ENDIF
 
      Z1 = HDATUM + HUP
      Y1 = Z1 - ZB1
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1,
     O            A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1)
 
      IF(Z1.GE.Z1TY6) THEN
C       This is type 6 flow
        TYPE = 6
        CALL DOTY6
     I            (STDOUT, CULCLS, A1, ALP1, K1, Z1, A2FULL, K2FULL,
     I             A3FULL, DDN, ZB3, CD6,
     O             Z3P, Q3, Z2P)
        Y2P = Z2P - ZB2
        Y3 = DDN
        Z2 = Z2P
        Y2 = DUP
        BETAF = 0.0
        BETA3 = 0.0
        ALPHA3 = 0.0
      ELSEIF(CD1T6.GT.0.0) THEN
C       We are between type 1 and type 6.  This is type 62. Full.
C       Interpolate the CD and for Z3P.
        TYPE = 62
        IF(Z1.GT.Z1TY1.AND.Z1.LT.Z1TY6) THEN
C         Proper interval.
          CD = CD1T6 + (Z1 - Z1TY1)*(CD6 - CD1T6)/(Z1TY6 - Z1TY1)
          Z3P= Z3PART + (Z1 - Z1TY1)*(Z3PTY6 - Z3PART)/(Z1TY6 - Z1TY1)
          CALL FULBAR
     I               (STDOUT, A1, ALP1, K1, Z1, A2FULL, K2FULL, A3FULL,
     I                CD, A2FULL, Z3P,
     O                Q3, Z2P)
          Z2 = Z2P
          Y2 = DUP
          BETAF= BT1T6 + (Z1 - Z1TY1)*(BT3ATD - BT1T6)/
     A                                         (Z1TY6 - Z1TY1)
          ALPHAF= AP1T6 + (Z1 - Z1TY1)*(AP3ATD - AP1T6)/
     A                                         (Z1TY6 - Z1TY1)
          CDF = CD
          AVH = A2FULL
          Y3PF = Z3P - ZB3
        ELSE
          WRITE(STDOUT,*) ' PROBLEM WHEN CD1T6 > 0 IN FRFT6'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
      ELSE
        IF(CD2T6.GT.0.0) THEN
C         No type 61 in this case.  Between type 2 and type 6. Full.
          TYPE = 62
          IF(Z1.GT.Z1TY1.AND.Z1.LT.Z1TY6) THEN
            CD = CD2T6 + (Z1 - Z1TY2)*(CD6 - CD2T6)/(Z1TY6 - Z1TY2)
            Z3P= Z3PART + (Z1 - Z1TY2)*(Z3PTY6 - Z3PART)/(Z1TY6 - Z1TY2)
 
            CALL FULBAR
     I                 (STDOUT, A1, ALP1, K1, Z1, A2FULL, K2FULL,
     I                  A3FULL, CD, A2FULL, Z3P,
     O                  Q3, Z2P)
            Z2 = Z2P
            Y2 = DUP
            BETAF= BT2T6 + (Z1 - Z1TY2)*(BT3ATD - BT2T6)/
     A                                           (Z1TY6 - Z1TY2)
            ALPHAF= AP2T6 + (Z1 - Z1TY2)*(AP3ATD - AP2T6)/
     A                                           (Z1TY6 - Z1TY2)
            CDF = CD
            AVH = A2FULL
            Y3PF = Z3P - ZB3
          ELSE
            WRITE(STDOUT,*) ' PROBLEM WHEN CD2T6 > 0 IN FRFT6.'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
        ELSEIF(CD61T6.GT.0.0) THEN
          IF(Z1.GT.Z1TY61.AND.Z1.LT.Z1TY6) THEN
C           Between type 61 and 6.  Full.
            TYPE = 62
            CD = CD61T6 + (Z1 - Z1TY61)*(CD6 - CD61T6)/
     A                                         (Z1TY6 - Z1TY61)
            Z3P= Z3PART + (Z1 - Z1TY61)*(Z3PTY6 - Z3PART)/
     A                                          (Z1TY6 - Z1TY61)
            CALL FULBAR
     I                 (STDOUT, A1, ALP1, K1, Z1, A2FULL, K2FULL,
     I                  A3FULL, CD, A2FULL, Z3P,
     O                  Q3, Z2P)
            Z2 = Z2P
            Y2 = DUP
            BETAF= BT61T6 + (Z1 - Z1TY61)*(BT3ATD - BT61T6)/
     A                                          (Z1TY6 - Z1TY61)
            ALPHAF= AP61T6 + (Z1 - Z1TY61)*(AP3ATD - AP61T6)/
     A                                          (Z1TY6 - Z1TY61)
            CDF = CD
            AVH = A2FULL
            Y3PF = Z3P - ZB3
          ELSE
C           Z1 < Z1TY61 here.
 
            IF(CD2.LT.0.0) THEN
C             There is no type 2.  Therefore, type 61 flow has
C             a constant CD.
              CD = CD6
              Z3P = ZB3 + 0.5*(Z1 - ZB3)
              AVH = A2FULL
              BETAF = -2.0
            ELSE
C             There is type 2.  Type 61 flow has a CD that varies
C             from type 2 value at its limit to the CD at the
C             Type 61 limit.
 
              CD = CD2 + (Z1 - Z1TY2)*(CD61 - CD2)/(Z1TY61 - Z1TY2)
              AVH = AVH2 + (Z1 - Z1TY2)*(AVH61 - AVH2)/(Z1TY61 - Z1TY2)
              IF(Z3PEST.EQ.ZB3) THEN
C               Type 2 is possible but the head range given by the
C               user did not invoke it.
                Z3P = ZB3 + 0.5*(Z1 - ZB3)
              ELSE
                Z3P = Z3PEST
              ENDIF
              BETAF = -1.0
            ENDIF
            TYPE = 61
C            WRITE(STDOUT,*) ' TYPE 61 CD=',CD,' Z3P=',Z3P,' Z1=',Z1,
C     A                        ' Z3PEST=',Z3PEST
            CALL DOTY61
     I                 (STDOUT, A1, ALP1, K1, Z1, IU, ID, ZB3, CD, AVH,
     M                  Z3P,
     O                  Q3, Z2P)
            Z3PEST = Z3P
            Z2 = Z2P
            Y2 = DUP
            CDF = CD
            AVHF = AVH
            Y3PF = Z3P - ZB3
          ENDIF
        ENDIF
      ENDIF
 
      Q1 = WFRDF + Q3
      Q2 = Q3
 
      QFREE = Q3
 
C     MAKE SURE THE PROPER AREA IS USED FOR DEPARTURE REACH
      IF(TYPE.EQ.61) THEN
        Y3 = Z3P - ZB3
      ELSE
       Y3 = DDN
      ENDIF
      Z3 = Y3 + ZB3
      Y3FREE = Y3
      Y2FREE = Y2
      Q3FREE = Q3
      LFTYPE = TYPE
      LSTYPE = -1
      TP2 = MIN(Y3,DDN)
      CALL XLKTAL
     I           (ADRXS3,
     M            TP2,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3)
C     First moment of area should be at piezometric level-not the water
C     surface level.   These levels differ for type 6 flow.
      TP2 = Z3P-ZB3
      CALL LKTJ
     I         (ADRXS3,
     M          TP2,
     O          J3)
 
 
C     DETERMINE DEPARTURE SECTION VALUES AND FREED
 
C      WRITE(STDOUT,*) 'FRFT6: Z3P=',Z3P,' Y3=',Y3,' Z3=',Z3
 
      CALL DPM26
     I          (STDOUT, WFRDF, MFRDF, Z3P,
     O           EXPFLG)
      IF(EXPFLG.EQ.0) THEN
        RETURN
      ENDIF
 
C     DEPARTURE SECTION VALUES ARE IN XS4COM
 
      FREED = Z1 - Z4
 
      RETURN
      END
C
C
C
      SUBROUTINE   FRFT56
     I                   (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O                    EFLAG, TYPE, EXPFLG, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute a flow of type 5 or type 6.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL FREED, HDATUM, HUP, QFREE, ZDATUM
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     ZDATUM - datum for local elevation
C     HUP    - Head upstream
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TYPE   - Culvert flow type
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'typlim.cmn'
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FRFT5, FRFT6
C***********************************************************************
      IF(HHTYPE.EQ.5) THEN
C       Culvert is steep enough to have type 5 flow.
        TYPE = 5
        LFTYPE = 5
        LSTYPE = -1
 
        CALL FRFT5
     I            (STDOUT, HDATUM, HUP, IU, ID, CULCLS,
     O             TYPE, EXPFLG, QFREE, FREED)
      ELSE
        TYPE = 6
        CALL FRFT6
     I            (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O             EFLAG, TYPE, EXPFLG, QFREE, FREED)
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   FRFT12
     I                   (STDOUT, HDATUM, ZDATUM, HUP, DUP, DDN, IU, ID,
     I                    CULCLS, TRUEA1,
     M                    CD1,
     O                    EFLAG, TYPE, CONFLG, EXPFLG, QFREE, FREED)
 
C     + + + PURPOSE + + +
C     Compute flow for type 1 or 2 in a culvert.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CONFLG, EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL CD1, DDN, DUP, FREED, HDATUM, HUP, QFREE, TRUEA1, ZDATUM
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HDATUM - Datum for measuring head
C     ZDATUM - datum for local elevation
C     HUP    - Head upstream
C     DUP    - vertical diameter of culvert barrel at upstream end
C     DDN    - vertical diameter of culvert barrel at downstream end
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     TRUEA1 - area at section 1
C     CD1 - Type 1 discharge coef.  -1 if type 1 not possible. 0 if
c              not tested yet, > 0 if known.
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TYPE   - Culvert flow type
C     CONFLG - CONFLG=0: flow contracts as it enters the culvert and
C              CONFLG=1: flow expands as it enters the culver
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'culcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL FG, SB
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FRFT1, FRFT2
C***********************************************************************
C     CHECK FOR BOTTOM SLOPE AT ENTRANCE AND ELIMINATE OBVIOUS
C     NON-TYPE 1 FLOW
 
      SB = (ZBVEC(IU) - ZBVEC(IU+1))/ABS(XVEC(IU) - XVEC(IU+1))
      IF(SB.GT.0.0 .AND. CD1 >= 0.0) THEN
C       TRY TYPE 1 FLOW
        TYPE = 1
        CALL FRFT1
     I            (STDOUT, HDATUM, HUP, DUP, DDN, IU, ID, CULCLS,
     I             TRUEA1,
     O             EFLAG, FG, TYPE, CONFLG, EXPFLG, QFREE, FREED)
 
        IF(TYPE.NE.1) THEN
          WRITE(STDOUT,*) ' Rejecting Type 1. Trying Type 2.'
C         TRY FLOW TYPE 2
          IF(TYPE.EQ.2) THEN
            CALL FRFT2
     I                (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O                 EFLAG, FG, TYPE, EXPFLG, QFREE, FREED)
          ELSE
            WRITE(STDOUT,*) ' Rejecting Type 2. Trying Type 5.'
C           FORCE COMPUTATION OF TYPE 5 OR 6 AND MAYBE 7.
            TYPE = 5
          ENDIF
        ENDIF
      ELSE
C       TRY FLOW TYPE 2
        TYPE = 2
        FG = 0.0
        CALL FRFT2
     I            (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O             EFLAG, FG, TYPE, EXPFLG, QFREE, FREED)
        IF(TYPE.NE.2) THEN
          WRITE(STDOUT,*) ' Rejecting Type 2. Trying Type 5.'
C         FORCE COMPUTATION OF TYPE 5 OR 6 AND MAYBE 7.
          TYPE = 5
        ENDIF
      ENDIF
 
      RETURN
      END
