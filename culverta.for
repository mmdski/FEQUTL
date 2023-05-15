C
C
C
      REAL FUNCTION   RTY7RF
     I                      (ZTAIL)
 
C     + + + PURPOSE + + +
C     Compute the residual for all cases when the tail water location
C     is at section 43 for submergence of flow over the road.  
C     Submergence for flow over the road and for flow through the 
C     culvert conduit is taken to be at section 43.  Flapgate losses
C     are also estimated here using what paltry information is available. 
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL ZTAIL
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ZTAIL  - elevation of water surface at section 43
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'x44com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'embcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'rty7c.cmn'
      INCLUDE 'embwrq.cmn'
      INCLUDE 'flapgate.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER NSFLAG, NTAB
      REAL DQ3, DROP, FLAP_FORCE, KFLAP, M43, M44, PDV, YC, YT,
     A     ZTLOW, ZTHIGH
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F4TO44, FNDCDE, LKTAB, SBFEBC, XLKT22, XLKTAL
C***********************************************************************
C      WRITE(OUTUN,*) 'ENTRY to RTY7RF: ZTAIL=',ZTAIL
 
C     Find the cross section elements at section 43.
      Z43 = ZTAIL
      Y43 = ZTAIL - ZB43
      CALL XLKTAL
     I           (ADRS43,
     M            Y43,
     O            A43, T43, DT43, J43, K43, DK43, BET43, DBET43, ALP43,
     O            DALP43)
 
C      WRITE(OUTUN,*) ' J43=',J43

C     Establish the flow over the road if any.  Any submergence is
C     included in the final values.
 
      IF(ZTAIL.GT.ZSUB) THEN
        CALL SBFEBC
     I             (OUTUN, Z1T, ZTAIL, PLCWTB, GLCWTB, PHCWTB, GHCWTB,
     I              PSUBTB, GSUBTB, NOFF, SURF, HLCRIT, XRDFL, XRDFR,
     I              HRDFL, HRDFM, HRDFR, QRDFL, QRDFM, QRDFR, TOTHL,
     I              TOTHM, TOTHR, YFL, YFM, YFR, APPL, APPM, APPR, WL,
     I              WM, WR, AELL, AELM, AELR, RMFFAC,
     O              WFRD, MFRD, EFRD)
      ELSE
        WFRD = WFRDF
        MFRD = MFRDF
        EFRD = EFRDF
      ENDIF
 
C     Now find the flow in the culvert for the given upstream level,
C     Z1T and the given tailwater, ZTAIL.
 
      DROP = Z1T - ZTAIL
C      WRITE(OUTUN,*) ' RTY7RF: DROP=',DROP
      DROP = SQRT(DROP)
C      WRITE(OUTUN,*) ' RTY7RF: SQRT(DROP)=',DROP
      CALL LKTAB
     I          (Q3ADR, DROP, 0,
     O           Q3, NTAB, DQ3)
 
C     Find the depth at section 3.
 
      CALL LKTAB
     I          (Y3ADR, ZTAIL, 0,
     O           Y3, NTAB, PDV)
      CALL XLKT22
     I           (ADRXS3,
     M            Y3,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3,
     O            Q3C)
      Z3 = Y3 + ZB3
C      WRITE(OUTUN,*) ' RTY7RF: Y3=',Y3,' A3=',A3,' J3=',J3,' Z3=',Z3,
C     A                ' Q3=',Q3
      Q4 = Q3 + WFRD
 
C     Compute the momentum flux plus impulse at section 43.  If
C     BETAF > 0, then we use a special value computed from ZTAIL
      YT = ZTAIL - ZB3
 
      IF(BETAF.GT.0.0)THEN
        IF(Y3PF.LT.DDN) THEN
          IF(YT.LT.DDN) THEN
            BETA3 = BETAF + (YT - Y3PF)*(BT3ATD - BETAF)/(DDN - Y3PF)
            ALPHA3 = ALPHAF + (YT - Y3PF)*(AP3ATD - ALPHAF)/(DDN - Y3PF)
          ELSE
            BETA3 = BT3ATD
            ALPHA3 = AP3ATD
          ENDIF
        ELSE
C         Special case in some cases in which the tailwater inducing
C         value is above the exit soffit.  Interpolate on tailwater
C         over 0.25 of the distance from the free flow limit elevation
C         of tailwater at section 3, Y3PF + ZB3, and the tailwater
C         level at zero flow, Z1T.
          ZTLOW = Y3PF + ZB3
          ZTHIGH = ZTLOW + 0.25*(Z1T - ZTLOW)
          IF(ZTAIL.GE.ZTHIGH) THEN
            BETA3 = BT3ATD
            ALPHA3 = AP3ATD
          ELSE
            BETA3 = BETAF + 
     A               (ZTAIL - ZTLOW)*(BT3ATD - BETAF)/(ZTHIGH - ZTLOW)
            ALPHA3 = ALPHAF + 
     A               (ZTAIL - ZTLOW)*(AP3ATD - ALPHAF)/(ZTHIGH - ZTLOW)
          ENDIF
        ENDIF
      ELSE
        BETA3 = BET3
        ALPHA3 = ALP3
      ENDIF
C      WRITE(OUTUN,*) ' BETAF=',BETAF,' YT=',YT,' DDN=',DDN,
C     A               ' BETA3=',BETA3
 
      M43 = MFRD + BETA3*Q3*Q3/A3 + GRAV*J43
 
      IF(T7FLAG.EQ.1) THEN
C       Finding the free flow for type 7 with critical flow at
C       section 4.  Find critical depth at section 4
 
        YC = Y4
        IF(YC.EQ.0.0) THEN
          YC = 0.5
        ENDIF
        CALL FNDCDE
     I             (OUTUN, ADRXS4, Q4,
     M              YC)
        Y4 = YC
        Z4 = Y4 + ZB4
        FLAP_FORCE = 0.0
      ELSE
C       Finding submerged flow with tailwater for road flow and
C       culvert flow at section 43.  Z4T gives the depth at
C       section 4 in this case.
        Z4 = Z4T
        Y4 = Z4 - ZB4

C       Include estimated flap gate losses if a flapgate is present.
C       Flapgate losses are not well established.  Only passing mention
C       of such losses are made in various trade literature put out
C       by those whose interest is in minimizing the size of the
C       loss.  The following equation has appeared in several locations
C       and is used until such time as something better is found. 
        IF(FLAP.GT.0.0) THEN
C         Compute the velocity head to use as if the conduits are 
C         flowing full even if they are not.  It appears that the loss
C         equations assumed full-conduit flow with the exit completely
C         submerged.  However, we cannot limit the losses to those
C         conditions without introducing nonsense abrupt changes in
C         the flow.  Therefore, we must include losses for submerged
C         part-full flow also. 
          FLAP_FORCE = FLAP*A3*4.0*(Q3/A3FULL)**2
     A                        *EXP(-1.15*Q3/(A3FULL*SQRT(DDN)))
        ELSE
          FLAP_FORCE = 0.0
        ENDIF            

      ENDIF
 
      CALL XLKTAL
     I           (ADRXS4,
     M            Y4,
     O            A4, T4, DT4, J4, K4, DK4, BET4, DBET4, ALP4, DALP4)
 
C     Find the values at section 44 from those at section 4
      CALL F4TO44
     I           (OUTUN,
     O            NSFLAG)
 
C      WRITE(OUTUN,*) ' NSFLAG=',NSFLAG,' Y4=',Y4,' Y44=',Y44,
C     A               ' J44=',J44,' A44=',A44

C     Compute the momentum flux plus impulse at section 44.
      M44 = BET44*Q4*Q4/A44 + GRAV*J44
 
      RTY7RF = (M43 - FLAP_FORCE)/M44 - 1.0
C      WRITE(OUTUN,*) ' AT EXIT RTY7RF=',RTY7RF,' M43=',M43,' M44=',M44
C      WRITE(OUTUN,*) ' Y43=',Y43,' J43=',J43,' A3=',A3,' BETA3=',BETA3,
C     A               ' MFRD=',MFRD
      RETURN
      END
C
C
C
      REAL FUNCTION   RY1GY2
     I                      (Y)
 
C     + + + PURPOSE + + +
C     Residual function for finding type 1 conditions when the
C     depth at section 2 is given.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - value of the unknown being sought
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'y1gy2.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL APPFAC, ARATIO, CDIN, FDRDW, VH1, VH2
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL DEGCON, FCD123
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL DEGCON, FCD123, GETFRF, XLKTAL
C***********************************************************************
C     Get the conditions at section 1
      Y1L = Y
      Z1L = ZB1 + Y1L
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1L,
     O            A1L, T1L, DT1L, J1L, K1L, DK1L, BET1L, DBET1L, ALP1L,
     O            DALP1L)
 
C     Compute the free flow over the roadway.  WFRDF in embcom.cmn.
 
      CALL GETFRF
     I           (Z1L,
     O            ZSBRDF, FDRDW)
      Q1L = Q2L + WFRDF
 
      VH1 = (Q1L/A1L)**2/GRAV2
      VH2 = (Q2L/A2L)**2/GRAV2
 
 
C     CHECK FOR EXPANSION INSTEAD OF CONTRACTION.
      IF(A1L.GT.A2L) THEN
C       WE HAVE A CONTRACTION(I.E. ACCELERATION OF FLOW)
C       DEFINE THE COEF OF DISCHARGE
        C123 = FCD123(OUTUN, 1, CLASS, D, Z1L)
 
C       MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C       CHANNEL CONTRACTION
 
        CD = DEGCON(C123, A1L, A2L)
 
        RY1GY2 = (ALP1L - APPLOS)*VH1 - (ALP2L + (1./CD**2 -1.))*VH2
     A         + Z1L - Z2L - APPLEN*Q1L*Q2L/(K1L*K2L)
      ELSE
C       WE HAVE AN EXPANSION(I.E. DECELERATION OF FLOW)
        ARATIO = A1L/A2L
        IF(ARATIO.GT.0.95) THEN
C         Interpolate coefficients to make the transistion between
C         the two cases smooth.  Define the coefficient for standard
C         type 1 case.
 
          C123 = FCD123(OUTUN, 1, CLASS, D, Z1L)
 
C         MAKE ADJUSTMENTS TO THE COEF. OF DISCHARGE FOR THE DEGREE OF
C         CHANNEL CONTRACTION
 
          CDIN = DEGCON(C123, A1L, A2L)
 
          CD = CDIN + 20.0*(1.0 - ARATIO)*(0.98 - CDIN)
          APPFAC = 20.0*(1.0 - ARATIO)*APPEXP
        ELSE
          CD = 0.98
          APPFAC = APPEXP
        ENDIF
 
        RY1GY2 = (ALP1L - APPLOS)*VH1 - (ALP2L + (1./CD**2 -1.))*VH2
     A         + Z1L - Z2L - APPLEN*Q1L*Q2L/(K1L*K2L) -
     B              APPFAC*(ALP1L*VH1 - ALP2L*VH2)
 
      ENDIF
 
      END
C
C
C
      SUBROUTINE   FRFT7
     I                  (STDOUT, Z1TRUE, ID, ZSBRDF, OPTION, DZ4,
     I                   Q3VSRD, Y3VSTW,
     M                   Z4TRUE,
     O                   QCLV, FDROP)
 
C     + + + PURPOSE + + +
C     Find the flows with the tailwater for roadway flow given
C     at section 43.  Used for both free flow of type 7,
C     OPTION=FREE, and submerged flow, OPTION=SUBMERGE.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ID, Q3VSRD, STDOUT, Y3VSTW
      REAL DZ4, FDROP, QCLV, Z1TRUE, Z4TRUE, ZSBRDF
      CHARACTER OPTION*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     Z1TRUE - known elevation at section 1
C     ID     - Index for the downstream node for the culvert barrel
C     ZSBRDF - Tailwater elevation to cause submergence of flows
C               over the roadway
C     OPTION - specifies if the flow is FREE or SUBMERGE
C     DZ4    - estimated change in elevation at section 4
C     Q3VSRD - address of function table giving flow at section 3 versus
C              tailwater at section 43 for a given head at section 1
C     Y3VSTW - address of function table giving depth at section 3
C              versus the tailwater at section 43
C     Z4TRUE - known elevation of water surface at section 4
C     QCLV   - Flow in the culvert
C     FDROP  - Free drop value
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'rty7c.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, KNT
      REAL F, FL, FR, Y43MIN, YT, YTAIL, YTMAX, YTMIN, ZHIGH, ZLOW, ZT,
     A     ZTAIL
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL RTY7RF
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL RGF, RTY7RF
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT RGF CLAIMS NONE IN',
     A       ' FRFT7.')
 60   FORMAT(' *BUG:XXX* RGF: MORE THAN 100 ITERATIONS IN FRFT7.')
 62   FORMAT(' Low arg=',F10.4,' High arg=',F10.4,' Residual=',F10.4)
 70   FORMAT(' *BUG:XXX* RGF: INVALID FLAG RETURNED IN FRFT7.')
 74   FORMAT(' No solution found in FRFT7 after 100 tries.')
 75   FORMAT(/,' FRFT7: Minimum depth=',F8.3,' at section 43 found',
     A        ' seeking a negative residual.')
 76   FORMAT(/,' Accepting result with relative residual in impulse+',
     A    'momentum=',1PE10.3)
C***********************************************************************
      IF(Q3VSRD.EQ.0.OR.Y3VSTW.EQ.0) THEN
        WRITE(STDOUT,*) ' *BUG:XXX* Tailwater tables unknown in FRFT7.'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C     SET COMMON BLOCK VALUES FOR THE RESIDUAL FUNCTION
 
      OUTUN = STDOUT
      Z1T = Z1TRUE
      Z4T = Z4TRUE
      ZSUB = ZSBRDF
      Q3ADR = Q3VSRD
      Y3ADR = Y3VSTW
      IF(OPTION.EQ.'FREE    ') THEN
        T7FLAG = 1
      ELSE
        T7FLAG = 0
      ENDIF
      DDN = DVEC(ID)
      YTAIL = Z43OLD - ZB43
 
      IF(T7FLAG.EQ.1) THEN
        ZTAIL = Z43OLD
      ELSE
C       Use change in section 4 elevation as an estimate for the change
C       in the elevation at section 43.
        ZTAIL = Z43OLD + DZ4
      ENDIF
      YTAIL = ZTAIL - ZB43
 
C      WRITE(STDOUT,*) ' FRFT7: Z43OLD=',Z43OLD,' YTAIL=',YTAIL,
C     A              ' ZTAIL=',ZTAIL,' Z43MIN=',Z43MIN
 
C     Make a final check to avoid going above Z1T
      YTMAX = Z1T - ZB43
      YTMIN = Z43OLD - ZB43
      IF(YTAIL.GE.YTMAX) THEN
        YTAIL = 0.5*(YTMAX + YTMIN)
      ENDIF
 
      Y43MIN = Z43MIN - ZB43
      ZHIGH = 0.0
      ZLOW = 0.0
      KNT = 0
 100  CONTINUE
        ZTAIL = YTAIL + ZB43
        F = RTY7RF(ZTAIL)
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          WRITE(STDOUT,74)
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        IF(F.GT.0.0) THEN
C         Positive residual found.
          ZHIGH = ZTAIL
          FR = F
          IF(ZLOW.GT.0.0) GOTO 110
C         Negative residual not yet known.  Continue to search.
          YT = 0.95*YTAIL
          IF(YT.LT.Y43MIN) THEN
            YT = 0.5*(YTAIL + Y43MIN)
            IF((YT-Y43MIN)/Y43MIN.LE.1.E-6) THEN
              WRITE(STDOUT,75) Y43MIN
              IF(ABS(F).LE.0.005) THEN
                WRITE(STDOUT, 76) F
C               Values have been set in RTY7RF
                GOTO 120
              ENDIF
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
          YTAIL = YT
          GOTO 100
        ELSE
C         Negative residual found.
          ZLOW = ZTAIL
          FL = F
          IF(ZHIGH.GT.0.0) GOTO 110
C         Positive residual not yet known.  Continue to search.
          YT = 1.1*YTAIL
          ZT = YT + ZB43
          IF(ZT.GE.Z1T) THEN
            ZT = 0.5*(ZTAIL + Z1T)
            IF(Z1T - ZT.LE.0.0) THEN
              WRITE(STDOUT,*) ' FRFT7: Section 1 elevation reached',
     A              ' seeking a positive residual.'
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ENDIF
          ZTAIL = ZT
          YTAIL = ZT - ZB43
          GOTO 100
        ENDIF
 110  CONTINUE
 
C      WRITE(STDOUT,*) ' FRFT7 BEFORE RGF: ZLOW=',ZLOW,' ZHIGH=',ZHIGH
C      WRITE(STDOUT,*) ' FL=',FL,' FR=',FR
 
      CALL RGF
     I        (1.E-6, EPSF, RTY7RF,
     M         ZLOW, ZHIGH, FL, FR,
     O         Z43, FLAG)
      IF(FLAG.EQ.1) THEN
      	  WRITE(STDOUT, 50)
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        WRITE(STDOUT,62) ZLOW, ZHIGH, FL
        STOP 'Abnormal stop. Errors found.'
      ELSEIF(FLAG.EQ.3) THEN
        WRITE(STDOUT,70)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
 120  CONTINUE
      QCLV = Q3
      Q2 = Q3
      Q1 = Q4
      IF(T7FLAG.EQ.1) THEN
        FDROP = Z1T - Z4
        Z4TRUE = Z4
C        Q3FREE = Q3
C        Y3FREE = Y3
C        Y2FREE = Y2
C        WRITE(STDOUT,*) ' FRFT7: Y4=',Y4,' QCLV=',QCLV,' FDROP=',FDROP
C        WRITE(STDOUT,*) ' TYPE=',TYPE
C        WRITE(STDOUT,*) ' Q1=',Q1,' Q2=',Q2,' Q3=',Q3,' Q4=',Q4
      ENDIF
 
      Z3P = Z43
      Y3P = Z3P - ZB3
      Z43OLD = Z43
      Q43OLD = Q3
 
      RETURN
      END
C
C
C
      SUBROUTINE   QVSTW
     I                  (STDOUT, Z1TRUE, IU, ID, NFRAC, POWER, ZSBRDF,
     I                   CULCLS, AB, ALP1T, K1TRUE, HDATUM, ZDATUM,
     I                   FRTYPE,
     M                   EFLAG, NEXT,
     O                   Q3VSRD, Y3VSTW, Y2VSD)
 
C     + + + PURPOSE + + +
C     Compute the flow for various tailwater elevations, that is,
C     at section 43, for a given upstream head.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FRTYPE, ID, IU, NEXT, NFRAC, Q3VSRD, STDOUT,
     A        Y3VSTW, Y2VSD
      REAL AB, ALP1T, HDATUM, K1TRUE, POWER, Z1TRUE, ZDATUM, ZSBRDF
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     Z1TRUE - known elevation at section 1
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     NFRAC  - Number of fractions for defining partial free drop
C     POWER  - power to use to define break points for function fit
C     ZSBRDF - water surface elevation at section 43 that begins
C              submergence of flow over the roadway
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     AB     - true area at section 1 used to adjust CD for flow
C              types 1, 2, and 3
C     ALP1T  - value of energy flux correction coefficient at section 1
C     K1TRUE - conveyance at section 1 for the known elevation there
C     HDATUM - Datum for measuring head
C     ZDATUM - datum for local elevation
C     FRTYPE - free flow type
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NEXT   - pointer into scratch portion of ITAB/FTAB
C     Q3VSRD - address of function table giving flow at section 3 versus
C              tailwater at section 43 for a given head at section 1
C     Y3VSTW - address of function table giving depth at section 3
C              versus the tailwater at section 43
C     Y2VSD - address of function table giving depth at section 2
C              versus the tailwater at section 43
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'embcom.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'rqvstw.cmn'
      INCLUDE 'embwrq.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER MAXN
      PARAMETER (MAXN=150)
 
C     + + + LOCAL VARIABLES + + +
      INTEGER FLAG, I, IFLAG, J, KNT, KNT2, KNT3
      REAL ARGVEC(MAXN), DDN, DIR, DUP, DZ, F, F1(MAXN), F2(MAXN), FAC,
     A     FL, FOLD, FR, QCULV, QHIGH, QLOW, QOLD, Y2PMIN, Y3VEC(MAXN),
     B      ZTAIL, ZTVEC(MAXN)
      DOUBLE PRECISION DROP(MAXN), MCLVT(MAXN), QCLVT(MAXN), SLOPE,
     A                 SQRTDP(MAXN), TWVEC(MAXN), Y2VEC(MAXN)
      CHARACTER  CQCLV*8, CQROAD*8, ADJLOC(MAXN)*1

C     + + + INTRINSICS + + +
      INTRINSIC ABS, FLOAT, MIN, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL RQVSTW
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FULBAR, PUT1D, RGF3, RQVSTW, SBFEBC, VLCHPP, XLKT22,
     A         VAR_DECIMAL
 
C     + + + OUTPUT FORMATS + + +
 76   FORMAT(/,' Flows for water levels at section 43.',
     1   /,'    Road flow submergence at section 43.',
     A  //,' Upstream head=',F9.4,' Elevation=',F10.4,/,
     B /,1X,' Head at Depth at Head at Depth at Drop    Flow Flow   ',
     B    ' Flow   ',
     C /,1X,' section section  section section  section type through',
     C    ' over   ',
     D /,1X,'   2       2       3 & 43    3     1 to 43      barrel ',
     D    ' roadway ',
     E /,1X,' ------- -------- ------- -------- ------- ----  ------',
     E    ' -------')
 77   FORMAT(1X,F8.3,F9.3,F8.3,F9.3,F8.4,I5,2A8)
 78   FORMAT(
     A  1X, '-------------------------------------------------------',
     A      '---------',
     A    /, 1X,'  Notes: Barrel flows full at section 3 when depth',
     A ' there=',F8.3,
     B /,1X,    '         Barrel flows full at section 2 when depth',
     C ' there=>',F8.3)
 80   FORMAT(/,' Fit of flow through barrel versus square root of drop',
     A  ' is a piecewise',/,' cubic Hermite polynomial.  Cubic spline',
     B  ' fit was not variation',/,' limited at the following drop ',
     C  'values:')
 82   FORMAT(5X,F10.4)
 84   FORMAT(/,' Fit of flow through barrel versus square root of drop',
     A  ' is a '/,' variation-limited cubic spline.')
 86   FORMAT(/,' *WRN:567* QVSTW: Type 1 Profile problems not ',
     A  'resolved by flow adjustment.',/,11X,' Making tailwater',
     B  ' adjustments.')
 88   FORMAT(/,' *WRN:568* QVSTW: Making tailwater adjustment. Old',
     A   ' tailwater level=',F8.3)
 90   FORMAT(/,' *BUG:XXX* Too many tailwater-level adjustments.')
 92   FORMAT(/,' *ERR:695* Drop from section 1 to section 43=',F7.3,
     A  ' invalid.',/,11X,'Drop must be > 0.')
 94   FORMAT(/,' *WRN:569* Subatmospheric head=',F5.1,' at section 2',
     A    ' may cause cavitation',/,'  and loss of performance.')
 96   FORMAT(/,' *WRN:570* Subatmospheric head=',F5.1,' at section 2',
     A    ' is impossible.',/,'  Culvert will NOT perform as computed.')
 98   FORMAT(/,' BUG:XXX* Drop <= 0.0.  Cannot continue.')
C***********************************************************************
C      PVEC(1) = 0.D0
C      PVEC(2) = 0.25D0
C      PVEC(3) = 0.5D0
C      PVEC(4) = 0.75D0
C      PVEC(5) = 1.D0
      Y2PMIN = 1.E30
      DDN = DVEC(ID)
      FQTYPE = FRTYPE
 
      WRITE(STDOUT,76) Z1TRUE - HDATUM, Z1TRUE + ZDATUM
 
      CALL VAR_DECIMAL(Q3,
     O                 CQCLV)
      CALL VAR_DECIMAL(WFRD,
     O                 CQROAD)
      WRITE(STDOUT,77) Z2 - HDATUM, Y2, Z43 - HDATUM, MIN(Y3, DDN),
     A                   Z1TRUE - Z43, FRTYPE, CQCLV, CQROAD
 
C      WRITE(STDOUT,*) ' ENTRY TO QVSTW: NFRAC=',NFRAC,' POWER=',POWER
C      WRITE(STDOUT,*) ' Z1TRUE=',Z1TRUE
C      WRITE(STDOUT,*) ' Z43OLD=',Z43OLD,' Q43OLD=',Q43OLD
      DUP = DVEC(IU)
      IUP = IU
      IDN = ID
      OUTUN = STDOUT
      ZSUB = ZSBRDF
      CLASS = CULCLS
      ABASE = AB
      Z1T = Z1TRUE
 
      TWVEC(1) = Z43OLD
      QCLVT(1) = Q43OLD
      Y3 = Z43OLD - ZB3
      IF(Y3.GT.DDN) Y3 = DDN
      CALL XLKT22
     I           (ADRXS3,
     M            Y3,
     O            A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3,
     O            Q3C)
      Y3VEC(1) = Y3FREE
      Y2VEC(1) = Y2
C      WRITE(STDOUT,*) 'QVSTW: Y3VEC(1)=',Y3VEC(1)

      DROP(1) = Z1TRUE - Z43OLD
      SQRTDP(1) = SQRT(Z1TRUE - Z43OLD)
      QCULV = Q43OLD
      DZ = Z1TRUE - Z43OLD
      IF(DZ.LE.0.0) THEN
        WRITE(STDOUT,92)  DZ
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      ZTVEC(1) = Z43OLD
      DO 490 I=2,NFRAC-1
        ZTVEC(I) = Z43OLD + DZ*(1.0 -
     A               (FLOAT(NFRAC - I)/(NFRAC - 1))**POWER)
 490  CONTINUE
      ZTVEC(NFRAC) = Z1TRUE
      DO 500 I=2,NFRAC-1
C       Select the sequence of tail water values so that they are
C       more closely spaced as the upstream water level is
C       approached.  This puts more values in the region of rapid
C       variation of flow.
 
        ZTAIL = ZTVEC(I)
        Z3P = ZTAIL
        Y3 = Z3P - ZB3
        Y3P = Y3
        IF(Y3.GT.DDN) Y3 = DDN
C       Establish the flow over the road if any.  Any submergence is
C       included in the final values.
 
        IF(ZTAIL.GT.ZSUB) THEN
          CALL SBFEBC
     I             (OUTUN, Z1TRUE, ZTAIL, PLCWTB, GLCWTB, PHCWTB,
     I              GHCWTB, PSUBTB, GSUBTB, NOFF, SURF, HLCRIT, XRDFL,
     I              XRDFR, HRDFL, HRDFM, HRDFR, QRDFL, QRDFM, QRDFR,
     I              TOTHL, TOTHM, TOTHR, YFL, YFM, YFR, APPL, APPM,
     I              APPR, WL, WM, WR, AELL, AELM, AELR, RMFFAC,
     O              WFRD, MFRD, EFRD)
        ELSE
          WFRD = WFRDF
          MFRD = MFRDF
          EFRD = EFRDF
        ENDIF
 
 
        IF(LSTYPE.EQ.4) THEN
C         Full flow with complete downstream submergence.
          CD = C46
          AVH = A2FULL
          CALL FULBAR
     I               (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL, K2FULL,
     I                A3FULL, CD, AVH, Z3P,
     O                QCULV, Z2)
          Y2 = DUP
          Y3 = DDN
          Y2P = Z2 - ZB2
        ELSEIF(LSTYPE.EQ.42) THEN
C         Full flow but with possible variable CD,  and
C         perhaps BETA, and  ALPHA.  AVH is constant.
C         BETA and ALPHA are set and used when computing the
C         departure reach values.   The tailwater is below the
C         culvert soffit but about the piezometric level at the
C         culvert exit.
C         The free flow computations must create the values
C         in the TYPTRN common needed at the free flow end of the
C         tailwater range.  Type 42 ends at a tailwater equal to
C         the elevation of the exit soffit of the culvert barrel.
          AVH = A2FULL
 
          IF(Y3P.GE.DDN) THEN
C           Shift to type 4.
            CD = C46
            LSTYPE = 4
            CQTYPE = 4
          ELSE
            IF(LFTYPE.EQ.6) THEN
              CD = C46
            ELSE
              CD =  CDF + (Y3P - Y3PF)*(C46 - CDF)/(DDN - Y3PF)
            ENDIF
            LSTYPE = 42
            CQTYPE = 42
          ENDIF
          CALL FULBAR
     I               (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL, K2FULL,
     I                A3FULL, CD, AVH, Z3P,
     O                QCULV, Z2)
          Y2 = DUP
          Y3 = DDN
          Y2P = Z2 - ZB2
        ELSEIF(LFTYPE.EQ.6) THEN
C         Full flow with perhaps partial submergence but CD and
C         AVH are constant.
          CD = C46
          AVH = A2FULL
          IF(Y3P.GE.DDN) THEN
            LSTYPE = 4
            CQTYPE = 4
          ELSE
            LSTYPE = 42
            CQTYPE = 42
          ENDIF
          CALL FULBAR
     I               (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL, K2FULL,
     I                A3FULL, CD, AVH, Z3P,
     O                QCULV, Z2)
          Y2 = DUP
          Y3 = DDN
          Y2P = Z2 - ZB2
        ELSEIF(LFTYPE.EQ.62)  THEN
C         Full flow with variable values of CD, and perhaps BETA
C         and ALPHA.  AVH is constant.   See section for LSTYPE=42
C         above.
          AVH = A2FULL
          IF(Y3P.LT.DDN) THEN
            CD =  CDF + (Y3P - Y3PF)*(C46 - CDF)/(DDN - Y3PF)
            LSTYPE = 42
            CQTYPE = 42
          ELSE
            CD = C46
            LSTYPE = 4
            CQTYPE = 4
          ENDIF
          CALL FULBAR
     I               (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL, K2FULL,
     I                A3FULL, CD, AVH, Z3P,
     O                QCULV, Z2)
          Y2 = DUP
          Y3 = DDN
          Y2P = Z2 - ZB2
        ELSEIF(LFTYPE.EQ.61.AND.Y3P.GE.DDN) THEN
C         Exit is flowing full.
          CD = C46
          AVH = A2FULL
          LSTYPE = 4
          CQTYPE = 4
          CALL FULBAR
     I               (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL, K2FULL,
     I                A3FULL, CD, AVH, Z3P,
     O                QCULV, Z2)
          Y2 = DUP
          Y3 = DDN
          Y2P = Z2 - ZB2
        ELSEIF(LFTYPE.EQ.5.OR.LFTYPE.EQ.51.OR.LFTYPE.EQ.52) THEN
C         Transition values computed in FRFT5.  Transition
C         to type 4 from type 42.
          AVH = A2FULL
          IF(Y3P.LT.DDN) THEN
C           Transition region for Cd.
            CD =  CDF + (Y3P - Y3PF)*(C46 - CDF)/(DDN - Y3PF)
            LSTYPE = 42
            CQTYPE = 42
          ELSE
C           Pure type 4.
            CD = C46
            LSTYPE = 4
            CQTYPE = 4
          ENDIF
          CALL FULBAR
     I               (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL, K2FULL,
     I                A3FULL, CD, AVH, Z3P,
     O                QCULV, Z2)
          Y2 = DUP
          Y3 = DDN
          Y2P = Z2 - ZB2
        ELSE
C         Do all cases that have part full flow here.
          KNT3 = 0
 
 499      CONTINUE
C         Set the fixed tailwater into the common block for the
C         residual function, RQVSTW.
          Z43FIX = ZTAIL
          Y3P = ZTAIL - ZB3
 
C         Now find the flow in the culvert for the given upstream level,
C         Z1T and the given tailwater, ZTAIL.  The flow in the culvert
C         should always be in the submerged state.  Also at this
C         point it should be part full at least part of the way.
          KNT = 0
          KNT2 = 0
C          WRITE(OUTUN,*) ' STARTING SEARCH WITH QCULV=',QCULV
C         Search for an interval containing a root.  Set the argument
C         values at the limits of the root-containing interval to zero
C         to serve as flags for the finding of the residuals.
 
          QHIGH = 0.0
          QLOW = 0.0
          FAC = 1.05
          QOLD = 0.0
          FOLD = ZTAIL - Z1T
 100      CONTINUE
            KNT2 = KNT2 + 1
            IF(KNT2.GT.100) THEN
              WRITE(OUTUN,*) ' QVSTW: More than 100 tries for root.'
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            F = RQVSTW(QCULV)
C            WRITE(OUTUN,*) ' KNT=',KNT,' KNT2=',KNT2,' SBFLAG=',SBFLAG
C            WRITE(OUTUN,*) ' QVSTW: QCULV=',QCULV,' F=',F,' FOLD=',FOLD,
C     A                    ' QOLD=',QOLD
 
            IF(KNT.GT.6.AND.SBFLAG.EQ.-1) THEN
              SBFLAG = -2
            ENDIF
            IF(SBFLAG.EQ.-1.AND.LFTYPE.NE.1.AND.LFTYPE.NE.5) THEN
               SBFLAG = -2
            ENDIF
            IF(SBFLAG.EQ.-1) THEN
C             Problems computing the steady flow
C             profile.  Make adjustment based on the free flow
C             type.
              IF(LFTYPE.EQ.1) THEN
C               Free flow type was 1.  Likely problem with too large
C               an increase in flow.  Problem with consistent treatment
C               of energy losses with type 1 flow and its submergence.
C               Backup from the assumed increase and then increase
C               by a smaller amount.
                IF(FAC.GT.1.0) THEN
                  QCULV = QCULV/FAC
                ENDIF
                IF(KNT.EQ.0) THEN
                  FAC = 1.01
                ELSE
                  FAC = 1. + 0.5*(FAC - 1.0)
                ENDIF
                QCULV = FAC*QCULV
                WRITE(OUTUN,*) ' QVSTW: Making Type 1 adjustment.'
                WRITE(OUTUN,*) ' KNT=',KNT,' QCULV=',QCULV,' FAC=',FAC
                KNT = KNT + 1
                IF(KNT.GT.6) THEN
                  WRITE(OUTUN,86)
                  GOTO 100
                ENDIF
              ELSE
                WRITE(OUTUN,*) ' QVSTW: PROFILE PROBLEMS. LFTYPE=',
     A                       LFTYPE
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              GOTO 100
            ELSEIF(SBFLAG.EQ.-2) THEN
C             Problem in subroutine APPRO or in profile computation.
C             Adjust the tailwater  level to find a condition that
C             will permit a solution.
 
              KNT3 = KNT3 + 1
              IF(KNT3.GT.20) THEN
                WRITE(OUTUN,90)
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              WRITE(OUTUN,88) ZTAIL
              DO 102 J=I,NFRAC-1
C               Move each tailwater to the mid-point of the interval
C               above it.
                ZTVEC(J) = 0.5*(ZTVEC(J) + ZTVEC(J+1))
 102          CONTINUE
              ZTAIL = ZTVEC(I)
              WRITE(OUTUN,*) ' New tailwater level=',ZTAIL
              GOTO 499
            ENDIF
            IF(ABS(F).LE.EPSABS) THEN
C             Already close enough.  No need to search further.
            ELSE
C             Check on direction of movement.
              DIR = (F - FOLD)/(QCULV - QOLD)
              IF(DIR.GT.0.0) THEN
                IF(FAC.LT.1.0) THEN
                  FAC = 1.0/FAC
                ENDIF
              ELSE
                WRITE(OUTUN,*) ' QVSTW: Unexpected direction seeking',
     A                         ' a root.'
                IF(QLOW.EQ.0.0) THEN
C                 Set the known negative value.
                  FL = ZTAIL - Z1T
                  IF(FL.GT.-EPSABS) THEN
                    FL = -2.*EPSABS
                  ENDIF
                  QLOW = ABS(FL)
                  GOTO 110
                ENDIF
                IF(QHIGH.EQ.0.0) THEN
                  WRITE(OUTUN,*) ' POSSIBLE SEARCH FAILURE IN QVSTW'
C                  STOP 'Abnormal stop. Errors found.'
                  QCULV = FAC*QCULV
                  GOTO 100
                ENDIF
                IF(FAC.GT.1.0) THEN
                  FAC = 1.0/FAC
                ENDIF
              ENDIF
C              WRITE(OUTUN,*) ' KNT2=',KNT2,' DIR=',DIR,' FAC=',FAC
              FOLD = F
              QOLD = QCULV
              IF(F.GE.0.0) THEN
C               Positive residual found.
                QHIGH = QCULV
                FR = F
 
                IF(QLOW.GT.0.0) GOTO 110
                FL = ZTAIL - Z1T
                IF(FL.GT.-EPSABS) THEN
                  FL = -2.*EPSABS
                ENDIF
                GOTO 110
              ELSE
C               Negative residual found.
                QLOW = QCULV
                FL = F
                IF(QHIGH.GT.0.0) GOTO 110
C               Positive residual not yet known.  Continue to search.
                QCULV = FAC*QCULV
                GOTO 100
              ENDIF
 
 110          CONTINUE
 
              CALL RGF3
     I                 (0.0, EPSABS, RQVSTW,
     M                  QLOW, QHIGH, FL, FR,
     O                  QCULV, FLAG)
 
              IF(ABS(FL).GT.EPSDIF) THEN
                WRITE(OUTUN,*) ' QVSTW: Residual at convergence=',
     A                         ABS(FL),' > ',EPSDIF
              ENDIF
            ENDIF
 
C           Assign the flow type.
            IF(Y3P.LT.DDN) THEN
C             Exit soffit free of water
              IF(Y2.LT.DUP) THEN
                CQTYPE = 3
              ELSE
                CQTYPE = 41
              ENDIF
            ELSE
C             Exit soffit is under water
              IF(Y2.LT.DUP) THEN
                CQTYPE = 31
              ELSE
                CQTYPE = 4
              ENDIF
            ENDIF
 
C            WRITE(OUTUN,*) ' RQVSTW: CQTYPE=',CQTYPE
            LSTYPE = CQTYPE
 
C           Free surface at section 2.  Set false value in Y2P to
C           prevent false warning about subatmospheric pressures.
            Y2P = 1.E30
            IF(CQTYPE.EQ.4) THEN
C             Flow is not part full after all.  Compute with full flow
C             in the culvert barrel.
              CALL FULBAR
     I                   (STDOUT, AB, ALP1T, K1TRUE, Z1TRUE, A2FULL,
     I                    K2FULL, A3FULL, C46, A2FULL, Z3P,
     O                    QCULV, Z2)
C              WRITE(STDOUT,*) ' FULBAR AFTER RGF3 RETURNS QCULV=',QCULV
              LSTYPE = 4
              CQTYPE = 4
              Y2 = DUP
              Y3 = DDN
              Y2P = Z2 - ZB2
           ENDIF
        ENDIF
 
        Q1 = QCULV + WFRD
        Q2 = QCULV
        Q3 = QCULV
 
C       Find the cross section elements at section 3.
        CALL XLKT22
     I             (ADRXS3,
     M              Y3,
     O              A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3, DALP3,
     O              Q3C)
        Z3 = Y3 + ZB3
 
 
C       Store the values for later processing
        TWVEC(I) = ZTAIL
        QCLVT(I) = QCULV
        DROP(I) = Z1TRUE - ZTAIL
        SQRTDP(I) = SQRT(Z1TRUE - ZTAIL)
 
        Y3VEC(I) = Y3
        Y2VEC(I) = Y2
        CALL VAR_DECIMAL(QCULV,
     O                   CQCLV)
        CALL VAR_DECIMAL(WFRD,
     O                   CQROAD)
        WRITE(STDOUT,77) Z2 - HDATUM, Y2, ZTAIL - HDATUM, Y3,
     A                   Z1TRUE - ZTAIL, CQTYPE, CQCLV, CQROAD
        IF(DROP(I).LE.0.0) THEN
          WRITE(STDOUT,98)
          EFLAG = 1
          RETURN
        ENDIF
 
C       Get minimum piezometric level at section 2.
        Y2PMIN = MIN(Y2PMIN, Y2P)
 500  CONTINUE
 
      CALL VAR_DECIMAL(0.0,
     O                 CQCLV)
      CALL VAR_DECIMAL(0.0,
     O                 CQROAD)
      WRITE(STDOUT,77) Z1TRUE - HDATUM, MIN(DUP, Z1TRUE - ZB2),
     A                 Z1TRUE - HDATUM, MIN(Z1TRUE - ZB3, DDN),
     B                 0.0, CQTYPE, CQCLV, CQROAD
 
      WRITE(STDOUT,78) DDN, DUP
 
      IF(DUP - Y2PMIN.GT.(22.0 -0.001*(ZB2 + ZDATUM))) THEN
C       Pressure is sufficiently subatmospheric to permit cavitation
C       near the entrance of the conduit.
        WRITE(STDOUT,94) DUP - Y2PMIN
        IF(DUP - Y2PMIN.GT.(34.0 - 0.001*(ZB2 + ZDATUM))) THEN
          WRITE(STDOUT,96) DUP - Y2PMIN
        ENDIF
      ENDIF
C     Store the zero flow point
      TWVEC(NFRAC) = Z1TRUE
      QCLVT(NFRAC) = 0.0
      DROP(NFRAC) = 0.0
      SQRTDP(NFRAC) = 0.0
      IF(Z1T-ZB3.LT.DDN) THEN
        Y3VEC(NFRAC) = Z1T - ZB3
      ELSE
        Y3VEC(NFRAC) = DDN
      ENDIF
      IF(Z1T-ZB2.LE.DUP) THEN
        Y2VEC(NFRAC) = Z1T - ZB2
      ELSE
        Y2VEC(NFRAC) = DUP
      ENDIF
 
C     Try fitting a cubic spline with square root of drop as
C     the argument and flow in the culvert as the function.  Use
C     end conditions of zero second moment at the left end and
C     force linear slope at the  right end.  Note: the flow at zero
C     drop is always zero.  Note: The "left end" is the end 
C     with index in the vectors of 1.  In this case it is the
C     free-flow limit.  The "right end" is the end with 
C     an index of NFRAC and is the submerged flow limit of
C     zero flow. 
 
      SLOPE = QCLVT(NFRAC-1)/SQRTDP(NFRAC-1)
C      CALL SPLINE(STDOUT, SQRTDP, QCLVT, NFRAC, 2, 0.D0, 1,
C     A                      SLOPE, MCLVT)
      CALL VLCHPP
     I           (STDOUT, NFRAC, SQRTDP, QCLVT, 2, 0.D0, 1, SLOPE,
     O            MCLVT, ADJLOC)
 
      IFLAG = 0
      DO 502 I=1,NFRAC
        IF(ADJLOC(I).ne.' ') THEN
C         At least one slope has been adjusted.
          IF(IFLAG.EQ.0) THEN
C           Write the heading.
            WRITE(STDOUT,80)
            IFLAG = 1
          ENDIF
          WRITE(STDOUT,82) SQRTDP(I)**2
        ENDIF
 502  CONTINUE
      IF(IFLAG.EQ.0) THEN
C       No point was adjusted.
        WRITE(STDOUT,84)
      ENDIF
 
C     Store the results into a temporary table of type 4.  Give a table
C     number of -1. Put values in single precision vectors and reverse
C     order to ascending value of argument.
 
      DO 510 I=1,NFRAC
        ARGVEC(NFRAC - I + 1) = SQRTDP(I)
        F1(NFRAC - I + 1) = QCLVT(I)
        F2(NFRAC - I + 1) = MCLVT(I)
 510  CONTINUE
C      WRITE(STDOUT,*) ' QVSTW BEFORE PUT1D: NEXT=',NEXT
      CALL PUT1D
     I          (STDOUT, -1, 4, NFRAC, ARGVEC, F1, F2,
     M           NEXT,
     O           Q3VSRD)
 
C     Store the results for Y3 in a temporary table of type 2
      DO 520 I=1,NFRAC
        ARGVEC(I) = TWVEC(I)
        F1(I) = Y3VEC(I)
        F2(I) = 0.0
C        WRITE(STDOUT,*) 'I=',I,' TWVEC(I)=',TWVEC(I),
C     A   ' Y3VEC(I)=',Y3VEC(I)
 520  CONTINUE
 
      CALL PUT1D
     I          (STDOUT, -1, 2, NFRAC, ARGVEC, F1, F2,
     M           NEXT,
     O           Y3VSTW)
 
C     Process the values for Y2 in same way as the flow in the 
C     culvert except that the argument is the drop instead of
C     the square root of drop 
      SLOPE = (Y2VEC(NFRAC-1) - Y2VEC(NFRAC)) /DROP(NFRAC-1)
      CALL VLCHPP
     I           (STDOUT, NFRAC, DROP, Y2VEC, 2, 0.D0, 1, SLOPE,
     O            MCLVT, ADJLOC)

C      WRITE(STDOUT,*) ' Spline for Y2:'
C      DO 523 I=1,NFRAC
C        WRITE(STDOUT,9569) I, DROP(I), Y2VEC(I), MCLVT(I)
C9569  FORMAT(I5, F10.4,F10.4,1PE12.5)
C523   CONTINUE
 
      DO 525 I=1,NFRAC
        ARGVEC(NFRAC - I + 1) = DROP(I)
        F1(NFRAC - I + 1) = Y2VEC(I)
        F2(NFRAC - I + 1) = MCLVT(I)
 525  CONTINUE
      CALL PUT1D
     I          (STDOUT, -1, 4, NFRAC, ARGVEC, F1, F2,
     M           NEXT,
     O           Y2VSD)


C     Set the minimum value for tailwater elevation
      Z43MIN = TWVEC(1)
 
C      RETURN
 
C      WRITE(STDOUT,70)
C      WRITE(STDOUT,74)
C      WRITE(STDOUT,72) TWVEC(1), QCLVT(1),
C     A                   MCLVT(1), Y3VEC(1), DROP(1), SQRTDP(1)
C      WRITE(STDOUT,*) ' '
C      DO 600 I=2,NFRAC
C        H = SQRTDP(I) - SQRTDP(I-1)
C        H = TWVEC(I) - TWVEC(I-1)
C        DO 590 J=2,4
C          P = PVEC(J)
C          TW = TWVEC(I-1) + P*H
C          DRP = Z1T - TW
C          CALL LKTAB(Q3VSRD, SNGL(SQRT(DRP)), F, NTAB, DF, 1)
C          CALL LKTAB(Y3VSTW, SNGL(TW), Y3VAL, NTAB, PDV, 1)
C
C          PC = 1.D0 - P
C          FP = H*P*PC*(MCLVT(I-1)*PC - MCLVT(I)*P)
C     A         + (P + P + 1.D0)*PC**2*QCLVT(I-1)
C     B         + (3.D0 - P - P)*P**2*QCLVT(I)
C          FPP = PC*(1.D0 - 3.D0*P)*MCLVT(I-1)
C     A          - P*(2.D0 - 3.D0*P)*MCLVT(I)
C     B          + 6.D0*P*PC*(QCLVT(I) - QCLVT(I-1))/H
C          WRITE(STDOUT,73)    TW, F, DF, Y3VAL
C590     CONTINUE
 
C        WRITE(STDOUT,*) ' '
C        WRITE(STDOUT,72) TWVEC(I), QCLVT(I),
C     A                   MCLVT(I), Y3VEC(I), DROP(I), SQRTDP(I)
C        WRITE(STDOUT,*) ' '
C600   CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   CHKTEL
     I                   (STDOUT, CHK23,
     M                    DE34,
     O                    FLAG, CL34, DE14)
 
C     + + + PURPOSE + + +
C     Check the change in the elevation of the total energy line
C     and write warning messages if there is an increase in the
C     downstream direction.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CHK23, FLAG, STDOUT
      REAL CL34, DE14, DE34
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CHK23  - flag for checking energy balance between sections 2 and 3
C     DE34   - Drop in elevation of total energy line between section
C               3 and 4 for a culvert
C     FLAG   - Result flag
C     CL34   - Head loss factor between sections 3 and 4 for culverts
C     DE14   - Drop in elevation of total energy line between section
C               1 and 4 for a culvert
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'grvcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'x44com.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL DV, KEF3, VH3, VH4, VH44
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *WRN:530* Estimated level of total energy line at',
     A ' section 2=',F8.3,' > ',/,11X,'Estimated level of total',
     B ' energy line at section 1=', F8.3)
 52   FORMAT(' *WRN:531* Estimated level of total energy line at',
     A ' section 3=',F8.3,' > ',/,11X,'Estimated level of total',
     B ' energy line at section 2=', F10.2)
 54   FORMAT(' *WRN:532* Estimated level of total energy line at',
     A ' section 4=',F8.3,' > ',/,11X,'Estimated level of total',
     B ' energy line at section 3=', F8.3)
 58   FORMAT(' *WRN:535* Estimated level of total energy line at',
     A ' section 4=',F8.3,' > ',/,11X,'Estimated level of total',
     B ' energy line at section 1=', F8.3)
 70   FORMAT(11X,'BETA4=',F6.3,' ALPHA4=',F6.3,' BETA3=',F6.3,
     A           ' ALPHA3=',F6.3)
 72   FORMAT(11X,'BETA4=',F6.3,' ALPHA4=',F6.3,' BETA1=',F6.3,
     A    ' ALPHA1=',F6.3)
 74   FORMAT(11X,' ALPHA3=',F6.3,' ALPHA2=',F6.3)
C***********************************************************************
      IF(A1*A2*A3*A4.LE.0.0) THEN
        WRITE(STDOUT,*) ' PROBLEM IN CHKTEL: A ZERO AREA'
        WRITE(STDOUT,*) ' A1=',A1,' A2=',A2,' A3=',A3,' A4=',A4
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      FLAG = 1
 
      IF(Q2.LT.Q1) THEN
C       FLOW OVER THE ROAD IS PRESENT.
 
        IF(DE34.EQ.0.0) THEN
          KEF3 = ALP3
        ELSE
          KEF3 = DE34
        ENDIF
        VH3 = KEF3*(Q3/A3)**2/GRAV2
        VH4 = ALP4*(Q4/A4)**2/GRAV2
 
        ZTEL1 = ZB1 +  Y1 + ALP1*(Q1/A1)**2/GRAV2
        ZTEL4 = ZB4 +  Y4 + VH4
 
        DE14 = ZTEL1 - ZTEL4
 
        IF(ZTEL4.GT.ZTEL1+EPSDIF) THEN
          WRITE(STDOUT,58) ZTEL4, ZTEL1
          WRITE(STDOUT,72) BET4, ALP4, BET1,ALP1
        ENDIF
 
      ELSE
        IF(DE34.EQ.0.0) THEN
          KEF3 = ALP3
        ELSE
          KEF3 = DE34
        ENDIF
        VH3 = KEF3*(Q3/A3)**2/GRAV2
        VH4 = ALP4*(Q4/A4)**2/GRAV2
        ZTEL1 = ZB1 +  Y1 + ALP1*(Q1/A1)**2/GRAV2
        IF(CHK23.EQ.1) THEN
          ZTEL2 = ZB2 + Y2 + ALP2*(Q2/A2)**2/GRAV2
        ENDIF
        ZTEL3 = Z3P +       VH3
        ZTEL4 = ZB4 +  Y4 + VH4
 
 
C       APPLY CHECKS
 
        IF(CHK23.EQ.1) THEN
          IF(ZTEL2.GT.ZTEL1+EPSDIF) THEN
            FLAG = 0
            WRITE(STDOUT,50) ZTEL2, ZTEL1
          ENDIF
          IF(ZTEL3.GT.ZTEL2+EPSDIF) THEN
            FLAG = 0
            WRITE(STDOUT,52) ZTEL3, ZTEL2
            WRITE(STDOUT,74) KEF3, ALP2
          ENDIF
        ENDIF
        IF(ZTEL4.GT.ZTEL3+EPSDIF) THEN
          FLAG = 0
          WRITE(STDOUT,54) ZTEL4, ZTEL3
          WRITE(STDOUT,70) BET4, ALP4, BET3, KEF3
        ENDIF
 
        DE14 = ZTEL1 - ZTEL4
C       COMPUTE LOSS OF VELOCITY HEAD DIFFERENCE BETWEEN SECTION 3 AND 44
 
        VH44 = ALP44*(Q44/A44)**2/GRAV2
        DV = ABS(VH3 - VH44)
        DE34 = ZTEL3 - ZTEL4
        IF(DE34.LT.0.0.AND.ABS(DE34).LT.0.5*EPSDIF) THEN
          DE34 = 0.0
        ENDIF
        IF(DV.EQ.0.0) THEN
          CL34 = 0.0
        ELSE
          CL34 = DE34/DV
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   FRFCLV
     I                   (STDOUT, HDATUM, ZDATUM, HUP, DUP, DDN, IU, ID,
     I                    CULCLS,
     O                    EFLAG, TYPE, CONFLG, EXPFLG, ZSBRDF, QFREE,
     O                    FREED)
 
C     + + + PURPOSE + + +
C     Find the free flow rate in the culvert given by (IU, ID).
C     The upstream head is HUP relative to HDATUM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CONFLG, EFLAG, EXPFLG, ID, IU, STDOUT, TYPE
      REAL DDN, DUP, FREED, HDATUM, HUP, QFREE, ZDATUM, ZSBRDF
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
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TYPE   - Culvert flow type
C     CONFLG - CONFLG=0: flow contracts as it enters the culvert and
C              CONFLG=1: flow expands as it enters the culver
C     EXPFLG - Expansion flag.  EXPFLG=1 if flow expands on exit
C              and EXPFLG=0 if flow does not expand on exit
C     ZSBRDF - water surface elevation at section 43 that begins
C              submergence of flow over the roadway
C     QFREE  - Free flow
C     FREED  - Free drop
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'culcom.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs3com.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL TRUEA1, ZT12UP, ZUP
 
C     + + + INTRINSICS + + +
      INTRINSIC MAX
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FRFT0, FRFT12, FRFT56, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' *WRN:542* At upstream head=',F10.2, ' free flow type',
     A       ' unclear.',/,11X,' Unable to continue.')
C***********************************************************************
C     CLEAR FREE DROP FOR THE CULVERT TO PREVENT CONFUSION
 
      FREED = 0
      ZUP = HUP + HDATUM
 
C     CHECK FOR CONTROL AT SECTION 0.
      CALL FRFT0
     I          (STDOUT, HDATUM, HUP, IU, ID, CULCLS,
     O           EFLAG, TYPE, EXPFLG, ZSBRDF, QFREE, FREED)
      IF(EFLAG.NE.0) RETURN
      IF(TYPE.EQ.0) THEN
        RETURN
      ENDIF
 
C     RESET SECTION 1 VARIABLES. FRFT0 MAY HAVE CHANGED THEM
 
      Z1 = HDATUM + HUP
 
C     FIND THE ELEMENTS AT SECTION 1
 
      Y1 = HDATUM + HUP - ZB1
      CALL XLKTAL
     I           (ADRXS1,
     M            Y1,
     O            A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1)
C     SAVE VALUE OF TRUE AREA AT SECTION 1
 
      TRUEA1 = A1
C     The limits for type 1 and type 2 flow may not be computed
C     at entry to FRFCLV.  CD1=0.0 and CD2=0.0 if the respective
C     limit has not been computed.  If the limit has been sought and
C     does not exist then the CD value is negative.  If the limit
C     has been sought and does exist the CD value gives the discharge
C     coefficient at that limit.  The limits are computed within
C     FRFT1 and FRFT2 if they have not been sought earlier.
 
C      WRITE(STDOUT,*) ' CD1=',CD1,' CD2=',CD2
C      WRITE(STDOUT,*) ' Z1TY1=',Z1TY1,' Z1TY2=',Z1TY2
C      IF(CD1.EQ.0.0.OR.CD2.EQ.0.0) THEN
      IF(CD1.EQ.0.0.AND.CD2.EQ.0.0) THEN
C       Compute an estimated upper limit for the limit of type
C       1 and 2 flow
        ZT12UP = ZBVEC(IU) + 1.50*DUP
      ELSE
C       Take the maximum of the two limits as the upper limit
C       that must be exceeded before we no longer try to find
C       type 1 or type 2 flow.
        ZT12UP = MAX(Z1TY1, Z1TY2)
      ENDIF
 
C      WRITE(STDOUT,*) ' ZUP=',ZUP,' ZT12UP=',ZT12UP,' HUP12=',
C     A               ZT12UP - HDATUM
 
      IF(ZUP.LT.ZT12UP) THEN
C       ESTIMATED FLOW TYPE IS 1 OR 2
 
        CALL FRFT12
     I             (STDOUT, HDATUM, ZDATUM, HUP, DUP, DDN, IU, ID,
     I              CULCLS, TRUEA1,
     M              CD1,
     O              EFLAG, TYPE, CONFLG, EXPFLG, QFREE, FREED)
        IF(TYPE.GT.2) THEN
C         FLOW TYPE MUST BE 5 OR 6 OR 7
          IF(ZUP.GE.ZBVEC(IU) + DUP) THEN
            CALL FRFT56
     I                 (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O                  EFLAG, TYPE, EXPFLG, QFREE, FREED)
          ELSE
C           Free flow type unclear.  Set EFLAG and return.
            WRITE(STDOUT, 51) HUP
            EFLAG = 1
            RETURN
          ENDIF
        ELSE
C          WRITE(STDOUT,52) HUP, Z1, Z2, Z3P, Z4, TYPE, QFREE
C            WRITE(STDOUT,60) ZTEL1, ZTEL2, ZTEL3, ZTEL4
        ENDIF
      ELSE
C       ESTIMATED FLOW TYPE IS 5 OR 6
 
        CALL FRFT56
     I             (STDOUT, HDATUM, ZDATUM, HUP, IU, ID, CULCLS,
     O              EFLAG, TYPE, EXPFLG, QFREE, FREED)
 
        IF(TYPE.LT.5.AND.TYPE.NE.2) THEN
C          WRITE(STDOUT,*) ' REJECTING 6. TRY TYPE 1 AND 2.'
C         FLOW TYPE MUST BE 1 OR 2
          CALL FRFT12
     I               (STDOUT, HDATUM, ZDATUM, HUP, DUP, DDN, IU, ID,
     I                CULCLS, TRUEA1,
     M                CD1,
     O                EFLAG, TYPE, CONFLG, EXPFLG, QFREE, FREED)
          IF(TYPE.GT.2) THEN
            WRITE(STDOUT, 51) HUP
            EFLAG = 1
          ENDIF
        ENDIF
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   FHHTYP
     I                   (STDOUT, CULCLS, NFAC, HUP, HDATUM, DUP, IU,
     I                    ID, QRF, Y3LIMU,
     M                    Y3LIM,
     O                    HHTYPE, CCT5, YEXIT)
 
C     + + + PURPOSE + + +
C     Find the type of high head flow for the current upstream head.
C     Current method for box culverts does not depend on head but
C     done here anyway in anticipation of a better rule that will
C     depend on head.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER HHTYPE, ID, IU, STDOUT
      REAL CCT5, DUP, HDATUM, HUP, NFAC, QRF, Y3LIM, Y3LIMU, YEXIT
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     HUP    - Head upstream
C     HDATUM - Datum for measuring head
C     DUP    - vertical diameter of culvert barrel at upstream end
C     IU     - Index for upstream node for the culvert barrel
C     ID     - Index for the downstream node for the culvert barrel
C     QRF    - Flow over the roadway
C     Y3LIMU - full-flow-inducing depth at culvert exit for type 5 flow
C     Y3LIM  - full-flow-inducing depth at culvert exit for type 5 flow
C     HHTYPE - high-head type
C     CCT5   - contraction coefficient for type 5 flow
C     YEXIT  - depth at exit of the culvert
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'culcom.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'grvcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER IERR, IJUMPL, JMPLOC, ITEMP
      REAL CD6L, CD6LM, CONCO, FKNT, HEAD, Q5, QLM, RBMAX, SBLIM,
     A     TY6LSS, TY6TY5, Y3MAX, Y3T, YCAT3, YEND, YENDL, YHIGH,
     B     YHIGHL, YHLIM, YVC, YVCL
      CHARACTER PROTYP*8, PTYPEL*8
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FNDSZL
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FDCD46, FNDCDE, FNDQ5, FNDSZL, FROVD, LOCJMP
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *WRN:544* Rounding value=',F5.3,' > 0.03 taken as',
     A      ' 0.03 for Figure 16.')
 54   FORMAT(/,' Culvert barrel is smooth.')
 56   FORMAT(/,' Culvert barrel is rough.  Factor=',F8.2)
 58   FORMAT(/,' *WRN:546* Length to diameter ratio=',F8.3,' > 35',
     A  ' when finding',/,11X,'limiting slope for type 6/5 boundary.',
     B  ' Value is inaccurate.')
 60   FORMAT(/,' *WRN:563* Roughness factor for limiting slope for',
     A      ' rough pipes',/,11X,'< 0.10 or > 0.30.  Value is',
     B      ' inaccurate.')
 61   FORMAT(/,' *WRN:596* Roughness factor for limiting slope for',
     A      ' rough pipes > 0.60.  Set to 0.60.')
 62   FORMAT(/,' *WRN:564* Rounding value=',F5.3,' > 0.06 taken as',
     A      ' 0.06 for Figure 15.')
 64   FORMAT(/,' High-head flow type=6. Forced for FLARED inlet.')
 66   FORMAT(/,' High-head flow is type',I3,'.  Boundary So=',F8.4)
 68   FORMAT(' *BUG:XXX* H over D=',F8.3,' too small in FHHTYP.')
 70   FORMAT(/,' Checking flow type 5 in the culvert.')
 72   FORMAT(/,' Flow type 5 rejected.  Subcritical profile prevails',
     A         ' in the barrel.')
 74   FORMAT(/,' Flow type 5 rejected.  High side of jump > 0.8*D')
 76   FORMAT(/,' High-head flow type=6 after rejection of type 5.')
 78   FORMAT(/,' Seeking full-flow-inducing value for type 5 flow.')
 80   FORMAT(/,' Flow type 5 rejected.  Profile in barrel undefined.')
 82   FORMAT(/,' Flow type 5 rejected seeking full-flow-inducing',
     A          ' value.  High-head flow=6.')
 84   FORMAT(/,' User limit=',F8.3,' taken for full-flow-inducing',
     A          ' value. R/B at maximum.')
 85   FORMAT(/,' User limit=',F8.3,' taken for full-flow-inducing',
     A          ' value.  Type 5 profile undefined.')
 86   FORMAT('  Super critical end depth=',F8.3,' taken for',
     A         ' full-flow-inducing value.')
 87   FORMAT(/,' Enhanced R/B yields super critical profile and ',
     A         ' user R/B yields a jump.',/,' Taking end depth=',F8.3,
     B         ' from enhanced R/B profile for the',/,' full-flow-',
     C         'inducing value.')
 88   FORMAT(/,' Enhanced R/B yields super critical profile and ',
     A         ' user R/B yields a jump.',/,' Enhanced R/B profile',
     B         ' has end depth=',F8.3,' < ',F8.3,', the end depth',/,
     C         ' for the user R/B profile.  User limit=',F8.3,
     D         ' taken for full-flow-inducing value.')
 89   FORMAT(/,' Flow type 5 rejected. User R/B profile has jump high',
     A     ' side=',F8.3,' > ',F8.3,/,', the enhanced profile jump',
     B     ' high side.')
 90   FORMAT(/,' Full-flow-inducing value=',F8.3,' defined by',
     A        ' forcing high side of',/,' jump in user R/B profile to',
     B        ' match high side of jump in enhanced R/B profile.')
 98   FORMAT(/,' *BUG:XXX* Iteration count exceeded in FHHTYP while',
     A      ' seeking jump height.')
 99   FORMAT(/,' *BUG:XXX* Problem in FHHTYP when varying jump height',
     A        ' seeking full-flow-inducing value.')
C***********************************************************************
C      WRITE(STDOUT,*) ' FHHTYP: on entry: Y3LIM=',Y3LIM,
C     A               ' Y3LIMU=',Y3LIMU
      IF(CULCLS.EQ.'FLARED'.AND.RBVAL.GT.0.0) THEN
C       Assume the user wants type 6 flow.
        HHTYPE = 6
        WRITE(STDOUT,64)
      ELSE
        HEAD =  HUP + HDATUM - ZB2
        HOVERD = HEAD/DUP
        IF(HOVERD.LT.1.5) THEN
C         Should not happen.
          WRITE(STDOUT,68) HOVERD
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        IF(CULCLS.EQ.'BOX') THEN
C         Assign values that will not be used.
          NBAR = 0.001
          RGHFAC = 0.01
          WRITE(STDOUT,54)
        ELSE
C         Estimate the composite n value for the conduit or conduits.
          NBAR = NFAC*A2FULL*(DUP/4.0)**0.666667/K2FULL
C          WRITE(STDOUT,*) ' NBAR=',NBAR
          IF(NBAR.LT.0.019) THEN
C           Treat the conduit as smooth.
            RGHFAC = 0.01
            WRITE(STDOUT,54)
          ELSE
C           Treat the conduit as rough.  Compute the roughness factor.
            RGHFAC = 29.0*HEAD*(A2FULL*NFAC/K2FULL)**2
            WRITE(STDOUT,56) RGHFAC
            IF(RBVAL.GT.0.03) THEN
              WRITE(STDOUT,50) RBVAL
            ENDIF
            IF(RGHFAC.LT.0.1.OR.RGHFAC.GT.0.3) THEN
              WRITE(STDOUT,60)
            ENDIF
            IF(RGHFAC.GT.0.60) THEN
              WRITE(STDOUT, 61)
              RGHFAC = 0.60
            ENDIF
          ENDIF
        ENDIF
        IF(NBAR.LT.0.019.AND.RBVAL.GT.0.06) THEN
          WRITE(STDOUT,62) RBVAL
        ENDIF
        IF(LOVERD.GT.35) THEN
          WRITE(STDOUT,58) LOVERD
        ENDIF
        SBLIM = FNDSZL(STDOUT, CULCLS, LOVERD, NBAR, RGHFAC,
     A                     RBVAL, TB15AD, TB16AD)
        IF(SZERO.GT.SBLIM) THEN
          HHTYPE = 5
          WRITE(STDOUT,70)
C         Verify type 5 and compute the full-flow-inducing value for
C         it.  Type 5 flow is in a sense unstable and a transition to
C         type 6 may be likely.   Therefore type 6 is selected if there
C         is any reason to doubt type 5.
          CALL FNDQ5
     I              (STDOUT, HDATUM, HUP, IU, CULCLS, RBVAL, QRF,
     O               YVC, Q5, CCT5)
          CALL FDCD46
     I               (STDOUT, CULCLS, RBVAL,
     O                CD6L)
          TY6TY5 = ((1.0/CD6L**2 - 1.0))*(Q5/A2FULL)**2/GRAV2
          YCAT3 = 0.5*DUP
          CALL FNDCDE
     I               (STDOUT, ADRXS3, Q5,
     M                YCAT3)
C          WRITE(STDOUT,*) ' Critical depth at exit=',YCAT3
C         Find the profile in the barrel.
          CALL LOCJMP
     I               (STDOUT, IAT3D, ID, YVC, YCAT3, Q5, DUP,
     M                TY6TY5,
     O                YHIGH, JMPLOC, YEND, IERR, PROTYP)
          YEXIT = YEND
          IF(PROTYP.EQ.'SUB') THEN
C           Reject type 5 because critical depth at the exit at the
C           type 5 flow yields a subcritical profile that drowns
C           any jump created by the supercritical flow at the vena
C           contracta.   Negative sign used to distinguish type 6
C           after rejection of type 5 from selection of type 6
C           from the figures.
            HHTYPE = -6
            WRITE(STDOUT,72)
          ELSEIF(PROTYP.EQ.'MIXED') THEN
C           A jump exists in the barrel maintained by critical depth
C           at the barrel exit.  Is the jump too high?  The Froude
C           numbers encountered will tend to produce an undular jump.
C           Therefore the water surface downstream of the jump will
C           have waves that are a significant proportion of the depth.
C           These waves will tend to seal the barrel so that the
C           air cannot easily enter.
            IF(YHIGH.GT.0.8*DVEC(JMPLOC)) THEN
C             Reject type 5 flow.  Unlikely to prevail with a jump
C             this high.
              HHTYPE = -6
              WRITE(STDOUT,74)
            ENDIF
          ELSEIF(IERR.NE.0) THEN
C           Some problem found in computing the profile.  Could be
C           a nonsense culvert barrel.  However, try to continue
C           by making high-head flow type 6.
            HHTYPE = -6
            WRITE(STDOUT,80)
          ENDIF
        ELSE
          HHTYPE = 6
        ENDIF
        IF(HHTYPE.EQ.6) THEN
          WRITE(STDOUT,66) HHTYPE, SBLIM
        ELSEIF(HHTYPE.EQ.-6) THEN
          HHTYPE = 6
          WRITE(STDOUT,76)
        ELSE
          WRITE(STDOUT,66) HHTYPE, SBLIM
        ENDIF
        IF(HHTYPE.EQ.5) THEN
C         Seek the depth at the culvert exit that will induce full
C         flow in the barrel.
          WRITE(STDOUT,78)
 
C         Find the value of relative rounding/beveling that causes
C         the culvert slope at the boundary between type 5 and type 6 to
C         be the same as the culvert slope.
          CALL FROVD
     I              (STDOUT, CULCLS,
     M               ROVDLM)
          IF(ROVDLM.LE.RBVAL) THEN
C           The culvert already has a rounding/beveling value at the
C           limit for the figures.  Assign the full-flow-inducing value
C           to the user value.
            WRITE(STDOUT,84) Y3LIMU
            Y3LIM = Y3LIMU
          ELSE
            RBMAX = ROVDLM
            FKNT = 8
 100        CONTINUE
C           Find the type 5 flow for the enhanced rounding/beveling
            CALL FNDQ5
     I                (STDOUT, HDATUM, HUP, IU, CULCLS, ROVDLM, QRF,
     O                 YVCL, QLM, CONCO)
            YCAT3 = 0.5*DVEC(ID)
            CALL FNDCDE
     I                 (STDOUT, ADRXS3, QLM,
     M                  YCAT3)
            CALL FDCD46
     I                 (STDOUT, CULCLS, ROVDLM,
     O                  CD6LM)
            TY6LSS = ((1.0/CD6LM**2 - 1.0))*(QLM/A2FULL)**2/GRAV2
            CALL LOCJMP
     I                 (STDOUT, IAT3D, ID, YVCL, YCAT3, QLM, DUP,
     M                  TY6LSS,
     O                  YHIGHL, IJUMPL, YENDL, IERR, PTYPEL)
            IF(IERR.EQ.1) THEN
C             Take user value.
              Y3LIM = Y3LIMU
              WRITE(STDOUT,85) Y3LIM
            ELSE
C             See if enhanced rounding/beveling creates a meaningful
C             profile.
              IF(PTYPEL.EQ.'MIXED'.AND.YHIGHL.GT.0.8*DVEC(IJUMPL)) THEN
                ITEMP = 1
              ELSE
                ITEMP = 0
              ENDIF
              IF(PTYPEL.EQ.'SUB'.OR.ITEMP.EQ.1) THEN
C               Reduce the rounding and beveling value and try again.
C               The profile is not meaningful at the enhanced level.
                FKNT = FKNT - 1.0
                IF(FKNT.GT.0.0) THEN
                  ROVDLM = RBVAL + FKNT*(RBMAX - RBVAL)/8.0
                  GOTO 100
                ELSE
C                 No meaningful profile.  Reject type 5
                  HHTYPE = 6
                  WRITE(STDOUT,82)
                ENDIF
              ELSE
                IF(PTYPEL.EQ.'SUP') THEN
                  IF(PROTYP.EQ.'SUP') THEN
C                   The profiles are of like character.   Therefore
C                   take the end depth from the enhanced rounding/
C                   beveling as the full-flow-inducing value.
                    Y3LIM = YENDL
                    WRITE(STDOUT,86) Y3LIM
                  ELSE
C                   The profiles differ with the enhanced R/B being
C                   super critical and the user R/B having a jump.
C                   This case may not happen.
                    IF(YENDL.GT.YEND) THEN
C                     Use the enhanced R/B ending depth as the limit.
                      Y3LIM = YENDL
                      WRITE(STDOUT,87) Y3LIM
                    ELSE
C                     Result cannot be used.  Take user value.
                      Y3LIM = Y3LIMU
                      WRITE(STDOUT,88) YENDL, YEND, Y3LIM
                    ENDIF
                  ENDIF
                ELSE
C                 The enhanced R/B profile has a jump.  The user R/B
C                 profile may or may not have a jump.  If the user R/B
C                 profile has a jump and it is higher than the enhanced
C                 profile jump, then reject type 5 flow.  Otherwise
C                 find the exit depth, greater than critical depth,
C                 that yields a close match between the high side of
C                 the enhanced R/B profile jump and the high side of
C                 the user R/B profile jump.  Note that YHIGH is
C                 zero if no jump exists.
                  IF(YHIGH.GT.YHIGHL) THEN
C                   Reject type 5.
                    HHTYPE = 6
                    WRITE(STDOUT,89) YHIGH, YHIGHL
                  ELSE
C                   Both remaining cases here because YHIGH= 0 if there
C                   is no jump.  Find critical depth at culvert exit.
                    YCAT3 = 0.5*DUP
                    CALL FNDCDE
     I                         (STDOUT, ADRXS3, Q5,
     M                          YCAT3)
                    YHLIM = YHIGHL
 
C                   Establish the initial values for jump high side,
C                   YHIGHL, end depth, YENDL, and profile type, PTYPEL.
                    CALL LOCJMP
     I                         (STDOUT, IAT3D, ID, YVC, YCAT3, Q5, DUP,
     M                          TY6TY5,
     O                          YHIGHL, JMPLOC, YENDL, IERR, PTYPEL)
                    Y3MAX = 0.8*DVEC(ID)
                    FKNT = 1.0
 200                CONTINUE
                      Y3T =  YCAT3 + FKNT*(Y3MAX - YCAT3)/8.0
                      CALL LOCJMP
     I                           (STDOUT, IAT3D, ID, YVC, Y3T, Q5, DUP,
     M                            TY6TY5,
     O                            YHIGH, JMPLOC, YEND, IERR, PROTYP)
                      IF(IERR.NE.0) THEN
                        WRITE(STDOUT,99)
                        STOP 'Abnormal stop. Errors found.'
                      ENDIF
                      IF(PROTYP.EQ.'MIXED'.AND.PTYPEL.EQ.'MIXED') THEN
                        IF(YHLIM.GE.YHIGHL.AND.YHLIM.LE.YHIGH) THEN
C                         Containing interval found.  Interpolate
C                         for the end depth and use for the
C                         full-flow-inducing value.
                          Y3LIM = YENDL + (YHLIM - YHIGHL)*
     A                             (YEND - YENDL)/(YHIGH - YHIGHL)
                          WRITE(STDOUT,90) Y3LIM
                          GOTO 210
                        ENDIF
                      ENDIF
                      YHIGHL = YHIGH
                      PTYPEL = PROTYP
                      YENDL = YEND
                      FKNT = FKNT + 1.0
                      IF(FKNT.GT.8.0) THEN
                        WRITE(STDOUT,98)
                        STOP 'Abnormal stop. Errors found.'
                      ENDIF
                      GOTO 200
 210                CONTINUE
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C      WRITE(STDOUT,*) ' FHHTYP: on exit: Y3LIM=',Y3LIM,
C     A               ' Y3LIMU=',Y3LIMU
      RETURN
      END
C
C
C
      SUBROUTINE   CHKDEP
     I                   (STDOUT,
     M                    EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Check the cross section at section 43 against the cross section
C     at the culvert exit.  Section 43 must be at least a factor
C     of WIDFAC greater than the culvert exit section.  If it is
C     not, construct one that is and give the user a graphic warning
C     message.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDOUT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'depcom.cmn'
      INCLUDE 'ftable.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER XOFF22
      PARAMETER(XOFF22=8)
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADR, FLAG, I, IE, IS, SPACE, TYPE, XOFF3, XOFF43
      REAL MXA43, QC, T, TMAX, TOLD, Y, YOLD
      DOUBLE PRECISION AOLD, ASUM, HH, JSUM
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE, MAX, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FMXARG
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FMXARG, KIL, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,1X,'Analysis of the cross section at start of',
     A            ' departure reach')
 51   FORMAT(5X,' Width Factor=',F8.2)
 52   FORMAT(/,
     A  1X,' Depth in Depth in  Width factor x  Sect. 43  Difference',
     A/,1X,' sect. 3  sect. 43  culvert width   width',
     C/,1X,' -------- --------  --------------  --------  ----------')
 54   FORMAT(1X,F9.3,F9.3,6X,F10.3,F10.3,F12.3)
 56   FORMAT(/,1X,'Constructing a cross section that will be wide',
     A  ' enough.',/,' THIS CROSS SECTION SHOULD NOT BE USED',
     B  ' FOR FINAL RESULTS.',/,' You must review the departure',
     C  ' cross section you have supplied',/,' and make informed',
     D  ' adjustments to it to satisfy CULVERT.',/,' DO NOT',
     E  ' RELIE ON THE AUTOMATIC ADJUSTMENT MADE BY CULVERT.',
     F  /,' CULVERT forces a rectangular section with a width',
     G  ' of WIDFAC times the',/,' maximum horizontal barrel',
     H  ' opening.  This section is extended',/,' vertically',
     I  ' until the supplied departure cross section has',/,
     J  ' a width that exceeds the width of the rectangular',
     K  ' section.',//,' The adjusted table has alpha and beta',
     L  ' set to 1.0,'/,' the critical flow is computed from',
     M  ' A*SQRT(gA/T), the conveyance is ignored,'/,' and the',
     N  ' area and first moment of area are recomputed.')
C***********************************************************************
      WRITE(STDOUT,50)
      WRITE(STDOUT,51) WIDFAC
      WRITE(STDOUT,52)
 
C     FIND MAX ARGUMENT IN SECTION 43
 
      MXA43 = FMXARG(ADRS43)
      TYPE = ITAB(ADRS43+2)
      XOFF43 = OFFVEC(TYPE)
      XOFF3 = OFFVEC(ITAB(ADRXS3+2))
 
C     FOR EACH DEPTH TABULATED AT SECTION 3, CHECK THE TOP WIDTH
C     IN SECTION 43. DELETE THE LAST DEPTH SINCE THAT REPRESENTS THE
C     TOP OF THE SLOT IN THE CURRENT FEQUTL.  INCLUDE CHECK AGAINST
C     MAX ARGUMENT TO PREVENT ERRORS USING OLD FORM OF CLOSED
C     CONDUIT TABLES WHICH HAVE MORE THAN 2 ARGUMENTS IN THE SLOT.
 
C     Get the address of the first argument of the table at section 3
      IS = ADRXS3 + XTIOFF
C     Get the address of the last argument of the table at section 3
      ITMP = ITAB(ADRXS3)
      IE = ITMP - XOFF3
      FLAG = 0
      DO 100 I=IS,IE,XOFF3
        Y3 = FTAB(I)
        T3 = WIDFAC*FTAB(I+1)
        Y43 = Y3 + ZB3 - ZB43
        IF(Y43.LT.0.0) Y43 = 0.0
        IF(Y43.GT.MXA43)  GOTO 100
        CALL XLKTAL
     I             (ADRS43,
     M              Y43,
     O              A43, T43, DT43, J43, K43, DK43, BET43, DBET43,
     O              ALP43, DALP43)
        WRITE(STDOUT,54) Y3, Y43, T3, T43, T3 - T43
 
        IF(0.9999*T3.GT.T43) THEN
C         CROSS SECTION 43 IS TOO NARROW at this elevation
          FLAG = 1
        ENDIF
 100  CONTINUE
 
      IF(FLAG.EQ.0) THEN
C       TABLE 43 IS OK
        WRITE(STDOUT,'(/,A)') ' BEGTAB accepted.'
        RETURN
      ELSE
        WRITE(STDOUT,'(/,A)') ' BEGTAB rejected as too narrow.'
        WRITE(STDOUT,56)
 
C       BEGTAB REJECTED
C       FIND MAX WIDTH AT THE EXIT OF THE CULVERT
 
        TMAX = 0.0
        DO 120 I=IS,IE,XOFF3
          TMAX = MAX(TMAX, FTAB(I+1))
 120    CONTINUE
 
C       INCREASE TMAX
 
        TMAX = TMAX*WIDFAC
 
C       Construct a table and store in temporary part of the
C       function table system.
 
        ADR = FTP
 
C       Construct the table by copying the depths and top widths
C       from the user supplied table.  Replace each width that
C       is not at least as wide as TMAX.
 
        IS = ADRS43 + XTIOFF
        IE = ITAB(ADRS43)
        SPACE = I + XOFF43 - ADRS43
 
        IF(ADR+SPACE.GT.MRFTAB) CALL KIL
     I                                   (10,
     M                                    ADR, EFLAG)
 
C       Use the same table number again.
        FTAB(ADR+1) = FTAB(ADRS43 + 1)
        ITAB(ADR+2) = 22
        ITAB(ADR+3) = ADR + XTIOFF
        FTAB(ADR+4) = 0.0
        FTAB(ADR+5) = BEGELV
        FTAB(ADR+6) = FTAB(ADRS43+5)
        ITAB(ADR+21) = 0
        ADR = ADR + XTIOFF
 
C       Do the zero depth row.
        YOLD = FTAB(IS)
        TOLD = FTAB(IS + 1)
        IF(TOLD.LT.TMAX) TOLD = TMAX
        FTAB(ADR) = YOLD
        FTAB(ADR+1) = TOLD
        DO 190 I=2,7
         FTAB(ADR+I) = 0.0
 190    CONTINUE
        FTAB(ADR+4) = 1.0
        FTAB(ADR+6) = 1.0
        ASUM = 0.D0
        JSUM = 0.D0
        AOLD = 0.D0
 
        ADR = ADR + XOFF22
        IS = IS + XOFF43
        DO 200 I=IS,IE,XOFF43
          Y = FTAB(I)
          T = FTAB(I+1)
          IF(T.LT.TMAX) T = TMAX
          HH = 0.5D0*DBLE(Y - YOLD)
          ASUM = ASUM + HH*DBLE(T + TOLD)
          JSUM = JSUM + HH*((ASUM + AOLD) -
     A                HH*DBLE(T - TOLD)/3.D0)
 
          QC = ASUM*SQRT(GRAV*ASUM/T)
 
C         Now store the values
          FTAB(ADR) = Y
          FTAB(ADR+1) = T
          FTAB(ADR+2) = ASUM
          FTAB(ADR+3) = 100.
          FTAB(ADR+4) = 1.0
          FTAB(ADR+5) = JSUM
          FTAB(ADR+6) = 1.0
          FTAB(ADR+7) = QC
 
          AOLD = ASUM
          TOLD = T
          YOLD = Y
          ADR = ADR + XOFF22
 200    CONTINUE
C       Put address of last argument in the first element.
        ITAB(FTP) = ADR - XOFF22
        ADR = ADR + XOFF22
 
        ADRS43 = FTP
        FTP = ADR
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   FDCD46
     I                   (STDOUT, CULCLS, RBV,
     O                    CDIS)
 
C     + + + PURPOSE + + +
C     Find the discharge coefficient for flow types 4 and 6
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT
      REAL CDIS, RBV
      CHARACTER CULCLS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CULCLS - Class for culvert shape, BOX, PIPE, ..
C     RBV    - rounding/beveling value
C     CDIS   - discharge coefficient
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'cdcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER NTAB
      REAL DF
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTAB
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *BUG:XXX* Culvert class=',A8,' unknown in FDCD46.')
C***********************************************************************
      IF(CULCLS.EQ.'BOX') THEN
        IF(WWANGL.EQ.0.0) THEN
          CALL LKTAB
     I              (TB5ADR, RBV, 1,
     O               CDIS, NTAB, DF)
        ELSE
          IF(RBV.EQ.0.0) THEN
            IF(WWANGL.LE.30.0) THEN
C             Interpolate between .84 and .87
              CDIS = 0.84 + 0.001*WWANGL
            ELSEIF(WWANGL.LE.75.0) THEN
              CDIS = 0.87
            ELSE
              CDIS = .87 - 0.008*(WWANGL - 75.0)
            ENDIF
          ELSE
            IF(WWANGL.LE.30.0) THEN
              CDIS = .84 + 0.001*WWANGL
              CDIS = KRB*CDIS
            ELSEIF(WWANGL.LE.75.0) THEN
              CALL LKTAB
     I                  (TB5ADR, RBV, 1,
     O                   CDIS, NTAB, DF)
              IF(CDIS.LT.0.87) CDIS = 0.87
            ELSE
              CDIS = .87 - 0.008*(WWANGL - 75.0)
              CDIS = KRB*CDIS
            ENDIF
          ENDIF
        ENDIF
      ELSEIF(CULCLS.EQ.'PIPE') THEN
        CALL LKTAB
     I            (TB5ADR, RBV, 1,
     O             CDIS, NTAB, DF)
        CDIS = KPROJ*CDIS
      ELSEIF(CULCLS.EQ.'MITER') THEN
        CDIS = KPROJ*0.74
      ELSEIF(CULCLS.EQ.'FLARED') THEN
        CDIS = 0.90
      ELSEIF(CULCLS.EQ.'RCPTG') THEN
        CALL LKTAB
     I            (TB5ADR, RBV, 1,
     O             CDIS, NTAB, DF)
C        CDIS = 0.95
      ELSE
        WRITE(STDOUT,50) CULCLS
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
      SUBROUTINE   CLVIN
     I                  (SFAC, STDIN, STDOUT, MNBN, NBRA, NODEID, 
     I                   MFTNUM, MRFTAB, FTPNT, TYPFLG,
     M                   EFLAG, FTKNT, FTP,
     O                   NBN, BRPT, NSEC, XVEC, ZBVEC, KA, KD, HLTAB,
     O                   BNODID, NEGTAB, CHKBAR, IAT3D, IAT6D, BSHAPE,
     O                   sbkind)
 
C     + + + PURPOSE + + +
C     Input the culvert conduit description.  Similar to a regular
C     branch in FEQ but with fewer options
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER BSHAPE, CHKBAR, EFLAG, FTKNT, FTP, IAT3D, IAT6D, sbkind,
     a   STDIN, STDOUT, MFTNUM, MNBN, MRFTAB, NBN, NBRA, TYPFLG
      INTEGER BRPT(8,NBRA), FTPNT(MFTNUM), HLTAB(MNBN), NEGTAB(MNBN),
     A        NSEC(MNBN)
      REAL KA(MNBN), KD(MNBN), SFAC, XVEC(MNBN), ZBVEC(MNBN)
      CHARACTER BNODID(MNBN)*8, NODEID*4
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     SFAC   - Scale factor for stations
C     STDIN     - Fortran unit number for user input file
C     STDOUT   - Fortran unit number for user output and messages
C     MNBN   - Maximum number of nodes on branches
C     NBRA   - Number of branches
C     NODEID - Descriptor string for the culvert barrel nodes
C     MFTNUM - Maximum allowed table number
C     MRFTAB - maximum length of FTAB/ITAB
C     FTPNT  - function table pointer giving the table address for each
C              table number.  If the address is zero the table does not
C              exist.
C     TYPFLG - If TYPFLG=0 then any cross section table type is
C               adequate; else at least equivalent of old type 12
C               (new type 22) must be supplied.  TYPFLG = 1 when called
C               from CULVRT and TYPFLG=0 when called from from XSTMAK.
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTKNT  - function table counter
C     FTP    - next open location in the function table storage
C     NBN    - Number of nodes on branches
C     BRPT   - branch pointer table.  Values for each branch are:
C              ROW       Meaning
C              1         upstream user node number
C              2         downstream user node number
C              3         pointer into branch vector for upstream node
C              4         pointer into branch vector for downstream node
C              5         upstream exterior node number
C              6         downstream exterior node number
C              7         pointer to address in EMC for the branch
C              8         number of unknowns at a node for the branch
C     NSEC   - Table number for cross section at the node
C     XVEC   - Stations for nodes on culvert barrel
C     ZBVEC  - Bottom elevations at nodes for culvert barrel
C     KA     - Acceleration loss.
C     KD     - Deceleration loss.
C     HLTAB  - Point head loss table number
C     BNODID - Id string for nodes on culvert barrel
C     NEGTAB - Memory for negative table ids given explicitly by the
C              user.  Interpolation requests given with a single - 
C              generate a special negative internal number that is NOT
C              counted as a negative table id for NEGTAB. 
C     CHKBAR - flag, if 1, requests checking user node assignments
C     IAT3D  - index to vena contracta location
C     IAT6D  - index to vena contracta location
C     BSHAPE - barrel shape: 0-prismatic, 1-nonprismatic
 
C     + + + LOCAL VARIABLES + + +
      INTEGER MAXN
      PARAMETER (MAXN=8)
      INTEGER BNUM, FIRST, HL, I, IT, J, JZERO, K, LAST, LNODE, LXTAB,
     A        NHAT, NODE, XTAB, N, ITEM_START(MAXN), ITEM_END(MAXN),
     B        IDLEN, HLIDLEN, PFLAG

      REAL CA, CD, CLVLEN, D, DX, DXFAC, DZ, ELEV, FAC, FINT, LOVERD,
     A     SDXOVL, STAT, STATL, STATR, STDL, STDR, X, Z,
     B     XL, XR
      CHARACTER HEAD*120, LINE*120, NAME*8, XCHAR*10, XTABN*5, ZCHAR*10,
     A          JUST*5, TABID*16, HLTABID*16, BLANK*16, TEMPC*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, EXP, FLOAT, LOG
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LENSTR
      REAL GETD
      CHARACTER GETTOK*10, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKTAB, FNDELV, FNDSTA, GETD, GETTOK, inline, KIL,
     A         MKXVEC, SCAN, GET_ITEM_LIMITS, LENSTR, 
     B         READ_CULVERT_ITEMS
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(A80)
 13   FORMAT(I5,A5,2A10,2F5.0,I5)
 16   FORMAT(I5,1X,A8,1X,A5,2A10,2F5.0,I5)
 
C     + + + OUTPUT FORMATS + + +
 11   FORMAT(/,' ',A80)
 14   FORMAT(' ',I5,A5,F10.3,F10.3,2F5.2,I5)
 18   FORMAT(' ',I6,1X,A8,1X,A,F10.3,F10.3,2F5.2,A)
 50   FORMAT(' *ERR:108* CROSS SECTION TABLE NUMBER MUST BE > 0',
     A       ' AT FIRST NODE ON BRANCH')
 52   FORMAT(/,' Culvert barrel length <= 3 times its vertical',
     A   ' diameter.  Culvert losses',/,'  will be overestimated.')
 54   FORMAT(/,' Culvert barrel length <= 6 times its vertical',
     A   ' diameter.  Culvert losses',/,'  may be overestimated.')
 56   FORMAT(/,' FEQUTL will use',I4,' nodes for the culvert barrel.')
84    FORMAT('   Node  Node Id',A,'TabId    Station Elevation',
     A       '   KA   KD',A,'HL Id')
C***********************************************************************
      DXFAC = 0.1
      JUST = 'RIGHT'
      IDLEN = 0
      HLIDLEN = 0
      BLANK = ' '

C     SET UP INITIAL VALUES OF POINTERS, COUNTERS, and flags.

C     Clear the flag used to note one or more occurrences of the 
C     the special flag denoting location of an interpolated cross sections
C     using the proportion of the distance between adjacent actual cross sections.
      PFLAG = 0

C     If TYPFLG is 0 we must clear HLTAB to keep a record of occurances of
C     relative location given explicitly by the user.
      IF(TYPFLG.EQ.0) THEN
        DO 95 I=1,MNBN
          HLTAB(I) = 0
95      CONTINUE
      ENDIF
      NBN = 0
      LNODE = 0
      BNUM = 1
      CALL inline
     I           (STDIN, STDOUT,
     O            LINE)
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     N, ITEM_START, ITEM_END)

      READ(LINE,1) HEAD
C      WRITE(STDOUT,11) HEAD
 100  CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        IF(NODEID.EQ.'NO') THEN
          WRITE(STDOUT,*) ' CULVERT requires using NODEID=YES'
          STOP 'Abnormal stop.  Errors found.' 
        ELSE
          CALL READ_CULVERT_ITEMS(
     I                           STDOUT, LINE, N, ITEM_START,
     I                           ITEM_END,
     M                           EFLAG, IDLEN, HLIDLEN, 
     O                           NODE, NAME, XTABN, XCHAR, 
     O                           ZCHAR, CA, CD, HL)
C          READ(LINE,16,ERR=991) NODE, NAME, XTABN, XCHAR, ZCHAR, CA,
C     A                          CD, HL
        ENDIF
        IF(NODE.EQ.-1) GOTO 120
        IF(LNODE.NE.0) NODE = LNODE + 1
 
        IF(NODE.LE.0) CALL KIL
     I                         (2,
     M                          NODE, EFLAG)
        IF(NODE.LE.LNODE) CALL KIL
     I                             (3,
     M                              NODE, EFLAG)
        IF(XTABN.EQ.'-') THEN
          XTAB = -1000000
        ELSE
          READ(XTABN,'(I5)') XTAB
        ENDIF
 
C       CATCH THE FIRST NODE NUMBER AND SET THE POINTER
 
        IF(LNODE.NE.0) GOTO 102
          BRPT(1,BNUM) = NODE
          BRPT(3,BNUM) = NBN+1
          IF(XTAB.LE.0) THEN
            WRITE(STDOUT,50)
            EFLAG = 1
            XTAB = 1
          ENDIF
          LXTAB = XTAB
 102    CONTINUE
        LNODE = NODE
 
 
C       STORE VALUES FOR THIS NODE
 
        NBN = NBN + 1
        IF(NBN.GT.MNBN) CALL KIL
     I                           (27,
     M                            NBN, EFLAG)
        IF(XTAB.EQ.0) THEN
          XTAB = LXTAB
        ELSE
          LXTAB = XTAB
        ENDIF
        NSEC(NBN) = XTAB
 
C       Remember the locations of negative table ids given
C       by the user.
        IF(XTAB.LT.0.AND.XTAB.NE.-1000000) THEN
          NEGTAB(NBN) = 1
        ELSE
          NEGTAB(NBN) = 0
        ENDIF
 
C       PROCESS THE STATION AND ELEVATION VALUES
        XCHAR = GETTOK(XCHAR)
        ZCHAR = GETTOK(ZCHAR)
        IF(XCHAR(1:3).EQ.'TAB'.OR.XCHAR(1:3).EQ.'tab') THEN
C         FIND THE STATION FROM THE CROSS SECTION TABLE XTAB
C          X = FNDSTA(XTAB, STDOUT, EFLAG)
          CALL FNDSTA
     I               (XTAB, STDOUT,
     O                EFLAG, X)
        ELSEIF(XCHAR.EQ.'          ') THEN
          X = 0.0
        ELSEIF(XCHAR(1:2).EQ.'P='.OR.
     A         XCHAR(1:2).EQ.'p=') THEN
C         Read the p value. 
          READ(XCHAR(3:10),'(F8.0)',ERR=991) X
          PFLAG = 1
          HLTAB(NBN) = 1
        ELSE
          READ(XCHAR,'(F10.0)',ERR=991) X
        ENDIF
        IF(ZCHAR(1:3).EQ.'TAB'.OR.ZCHAR(1:3).EQ.'tab') THEN
C         FIND THE ELEVATION FROM THE CROSS SECTION TABLE XTAB
C          Z = FNDELV(XTAB, STDOUT, EFLAG)
          CALL FNDELV
     I               (XTAB, STDOUT,
     O                EFLAG, Z)
        ELSEIF(ZCHAR.EQ.'          ') THEN
          Z = 0.0
        ELSE
          READ(ZCHAR,'(F10.0)',ERR=991) Z
        ENDIF
        XVEC(NBN) = X
        ZBVEC(NBN) = Z
        KA(NBN) = CA
        KD(NBN) = CD
        IF(TYPFLG.EQ.1) THEN
          HLTAB(NBN) = HL
        ENDIF
        IF(NODEID.NE.'NO') BNODID(NBN) = NAME
 
        GOTO 100
 120  CONTINUE
C     Check if either endpoint invert elevation is exactly zero. 
C     If so, replace with a small positive number to prevent 
C     interpolation problems.
      IF(ZBVEC(1).EQ.0.0) THEN
        ZBVEC(1) = 1.E-5
      ENDIF
      IF(ZBVEC(NBN).EQ.0.0) THEN
        ZBVEC(NBN) = 1.E-5
      ENDIF
 
      IF(TYPFLG.EQ.1) THEN
        CHKBAR = 0
        IF(NBN.GT.2) THEN
C         Check if we can replace the user's definition with our own.
          BSHAPE = 0
          SBKIND = 0
          DO 122 I=2,NBN-1
            IF(NSEC(I).GT.0) THEN
C             Does it differ from the first table given?
              IF(NSEC(I).NE.NSEC(1)) THEN
C               The barrel is non-prismatic with a cross section
C               table given that differs from the first table and
C               precedes the last table.
                BSHAPE = 1
              ENDIF
            ENDIF
            IF(ZBVEC(I).NE.0.0) THEN
C             An invert elevation has been given before the
C             last table.  Bottom slope may have changed.
              SBKIND = 1
            ENDIF
 122      CONTINUE
          IF(SBKIND.GT.0.OR.BSHAPE.GT.0) THEN
C           The conduit appears to have a break in slope and/or is
C           non-prismatic even though the final section has not been
C           checked.  This is too complex for the simple code below
C           for automatic cross section assignment.  We will have to
C           check later for the proper stations for type 5 flow.
 
            CHKBAR = 1
 
          ELSE
C           The conduit has constant slope and is prismatic before
C           we check the last section.  Therefore, it can be
C           described by two nodes.
            NSEC(2) = NSEC(NBN)
            XVEC(2) = XVEC(NBN)
            ZBVEC(2) = ZBVEC(NBN)
            BNODID(2) = BNODID(NBN)
            NBN = 2
          ENDIF
        ENDIF
C       Check for the case of internally generated spacing of
C       cross sections.
        IF(NBN.EQ.2) THEN
C         Only two nodes given.  Signals that the stations are to
C         be internally generated.
c         Also the slope is constant
          SBKIND = 0
          CLVLEN = ABS(XVEC(1) - XVEC(2))
          IT = NSEC(1)
          CALL CHKTAB
     I               (12, STDOUT, FTPNT, MFTNUM,
     M                IT,
     O                EFLAG)
          IF(EFLAG.NE.0) THEN
            STOP 'Abnormal stop.  Errors found.' 
          ENDIF

          D = GETD(IT, STDOUT)/SFAC
          LOVERD = CLVLEN/D
          SDXOVL = DXFAC/LOVERD
C         Compute the number of Gauss points needed to give the
C         minimum distance increment relative to the culvert opening
C         at its upstream end.
 
          NHAT = EXP(3.63759 - 0.505495*(LOG(SDXOVL) + 6.93304)) + .4
 
          CALL MKXVEC
     I               (LOVERD, D,
     M                NHAT,
     O                XVEC, IAT3D, IAT6D)
          WRITE(STDOUT,56) NHAT
 
C         Check for length of the culvert.
          IF(IAT3D.EQ.0.OR.IAT3D.EQ.NHAT) THEN
C           Culvert is 3 D or less long.  Write message and set location
C           to entrance.
            WRITE(STDOUT,52)
            IAT3D = 1
          ENDIF
          IF(IAT6D.EQ.0) THEN
C           Culvert is less than 6 D long.  Leave IAT6D=0 to allow
C           for under representation of losses in short culverts.
            WRITE(STDOUT,54)
          ENDIF
 
C         Create new version of the culvert barrel description
 
          IF(NSEC(1).EQ.NSEC(2)) THEN
C            NEGTAB(1) = 0
C           Prismatic barrel.
            BSHAPE = 0
            DO 125 I=2,NHAT
              NSEC(I) = NSEC(1)
C              NEGTAB(I) = 0
 125        CONTINUE
          ELSE
C           Non-prismatic barrel.
            BSHAPE = 1
            NSEC(NHAT) = NSEC(2)
C            NEGTAB(1) = 0
C            NEGTAB(NHAT) = 0
            DO 130 I=2,NHAT-1
              NSEC(I) = -1000000
C              NEGTAB(I) = 1
 130        CONTINUE
          ENDIF
          ZBVEC(NHAT) = ZBVEC(2)
          KA(NHAT) = KA(2)
          KD(NHAT) = KD(2)
          BNODID(NHAT) = BNODID(2)
          HLTAB(NHAT) = 0
          DO 135 I=2,NHAT-1
            ZBVEC(I) = 0.0
            KA(I) = KA(NHAT)
            KD(I) = KD(NHAT)
            BNODID(I) = ' '
            HLTAB(I) = 0
 135      CONTINUE
 
          NBN = NHAT
          LNODE = LNODE + NHAT -2
        ENDIF
      ENDIF
      BRPT(2,BNUM) = LNODE
      BRPT(4,BNUM) = NBN
      LAST = NBN
      FIRST = BRPT(3,BNUM)
 
      KA(FIRST) = 0.0
      KD(FIRST) = 0.0
 
      IF(TYPFLG.EQ.0) THEN
C       Check for use of a p value for location 
        IF(PFLAG.EQ.1) THEN
C         Scan and evaluate p values.  They appear in sequence between
C         known cross sections. 
C         Clear the memory location of the first p value seen 
C         in a potential sequence of such values.
          IT = 0

          DO 2200 I=FIRST+1, LAST
c            WRITE(STDOUT,*) ' IT=',IT, ' HLTAB(I)=',HLTAB(I)
            IF(IT.EQ.0) THEN
              IF(HLTAB(I).EQ.1) THEN 
C               Set a pointer to remember location of the first p value 
C               that needs to be evaluated.  Note: If p values are used
C               in an interval then all locations must be in those terms.
C               This makes the detection easier and does not reduce the
C               usefulness in our current application.
         
c                WRITE(STDOUT,*) ' Found p-value. IT=',IT
                IT = I
              ENDIF
            ELSE
              IF(HLTAB(I).EQ.0) THEN
C               We have seen at least one p value.  The next non-p value
C               should be a valid station value for interpolation. Interpolate
C               using XVEC(IT-1) as the first X value and XVEC(I) as the last X value.
                XL = XVEC(IT-1)
                XR = XVEC(I)
c                WRITE(STDOUT,*) ' Interpolate with: XL=',XL,' XR=',XR 
                DO 2190 J=IT,I-1
                  XVEC(J) = XL + XVEC(J)*(XR - XL)
                  HLTAB(J) = 0
2190            CONTINUE
C               Clear the memory flag
                IT = 0
              ENDIF
            ENDIF
2200      CONTINUE
        ENDIF
      ENDIF

       
C     CHECK FOR REQUESTS FOR INTERPOLATION OF BOTTOM PROFILE
C     BEFORE CHECKING FOR REQUESTS FOR STATION AND BOTTOM
C     PROFILE.
C     NOTE THAT A GIVEN PROFILE LEVEL OF 0.0 IS TAKEN TO BE A
C     REQUEST FOR INTERPOLATION.  A TRUE VALUE OF ZERO CANNOT
C     BE GIVEN BUT IF AN ESSENTIAL ZERO IS DESIRED INPUT A
C     SMALL NON-ZERO VALUE TO OBTAIN IT.
 
C     AT A MINIMUM THE FIRST AND LAST NODE ON THE BRANCH MUST HAVE
C     NON-ZERO PROFILE VALUES.  THESE ARE NOT CHECKED.
 
      IF(ZBVEC(FIRST).NE.0.0) THEN
        JZERO = 0
C       SCAN ALL NODES BUT THE FIRST FOR REQUESTS
        DO 2400 J=FIRST+1,LAST
          IF(ZBVEC(J).EQ.0.0) THEN
C           REQUEST FOUND- SAVE INDEX OF FIRST REQUEST
            IF(JZERO.EQ.0) THEN
              JZERO = J
            ENDIF
          ELSE
C           NO REQUEST AT THIS POINT.  CHECK FOR PREVIOUS REQUESTS
            IF(JZERO.GT.0) THEN
C             ONE OR MORE INTERPOLATION REQUESTS ARE OUTSTANDING.
C             JZERO POINTS TO THE FIRST REQUEST. JZERO-1 POINTS
C             TO THE PREVIOUS VALID STATION-PROFILE PAIR AND
C             J POINTS TO THE CURRENT VALID PAIR.
 
C             GET THE STATION AND PROFILE FROM THE PREVIOUS VALID
C             PAIR
 
              STATL = XVEC(JZERO-1)
              STDL = ZBVEC(JZERO-1)
 
C             GET CURRENT VALID VALUE
 
              STATR = XVEC(J)
              STDR = ZBVEC(J)
 
              FAC = (STDR - STDL)/(STATR - STATL)
              DO 2300 K=JZERO, J-1
                IF(XVEC(K).NE.0.0) THEN
C                 EXACT ZERO FOR STATION PRECLUDES COMPUTING
C                 PROFILE. WILL BE DONE BY STATION AND PROFILE
C                 INTERPOLATION PROCESS.
 
                  ZBVEC(K) = STDL + (XVEC(K) - STATL)*FAC
                ENDIF
 2300         CONTINUE
              JZERO = 0
            ENDIF
          ENDIF
 2400   CONTINUE
      ENDIF
 
 
C     CHECK FOR REQUESTS FOR INTERPOLATION OF STATION AND
C     ELEVATION AND DO THE REQUESTS. NOTE: THE FIRST AND LAST
C     NODES ON A BRANCH ARE ASSUMED TO HAVE VALID VALUES
C     OF STATION AND ELEVATION.  THESE CANNOT BE CHECKED!
 
      JZERO = 0
C     SCAN ALL NODES BUT THE FIRST FOR REQUESTS
      DO 2000 J=FIRST+1,LAST
        IF(ZBVEC(J).EQ.0.0.AND.XVEC(J).EQ.0.0) THEN
C         REQUEST FOUND- SAVE INDEX OF FIRST REQUEST
          IF(JZERO.EQ.0) THEN
            JZERO = J
          ENDIF
        ELSE
C         NO REQUEST AT THIS POINT.  CHECK FOR PREVIOUS REQUESTS
          IF(JZERO.GT.0) THEN
C           ONE OR MORE INTERPOLATION REQUESTS ARE OUTSTANDING.
C           JZERO POINTS TO THE FIRST REQUEST. JZERO-1 POINTS
C           TO THE PREVIOUS VALID STATION-ELEVATION PAIR AND
C           J POINTS TO THE CURRENT VALID PAIR.
 
C           GET THE FLOAT PNT FORM OF THE NUMBER OF INTERVALS
 
            FINT = FLOAT(J - JZERO + 1)
 
C           GET THE STATION AND ELEVATION FROM THE PREVIOUS VALID
C           STATION-ELEVATION PAIR
 
            STAT = XVEC(JZERO-1)
            ELEV = ZBVEC(JZERO-1)
 
C           FIND THE CHANGE IN STATION AND ELEVATION
 
            DX = XVEC(J) - STAT
            DZ = ZBVEC(J) - ELEV
 
            DO 1000 K=JZERO, J-1
              FAC = FLOAT(K-JZERO+1)/FINT
              XVEC(K) = STAT + DX*FAC
              ZBVEC(K) = ELEV + DZ*FAC
 1000       CONTINUE
            JZERO = 0
          ENDIF
        ENDIF
 2000 CONTINUE
 
      NODE = BRPT(1,BNUM)
      IF(IDLEN.LT.6) IDLEN = 6
      IF(HLIDLEN.LT.6) HLIDLEN = 6
      WRITE(STDOUT,84) BLANK(1:IDLEN-5), BLANK(1:HLIDLEN-5)
      DO 4000 J=FIRST,LAST
        IF(NSEC(J).EQ.-1000000) THEN
          TABID = ' '
          TABID(IDLEN:IDLEN) = '-'
        ELSE
          IF(NSEC(J).GT.0) THEN
            TEMPC = GET_TABID(NSEC(J))
          ELSE
            TEMPC(2:) = GET_TABID(-NSEC(J))
            TEMPC(1:1) = '-'
          ENDIF

          IT = LENSTR(TEMPC)
          TABID = ' '
          TABID(IDLEN-IT+1:IDLEN) = TEMPC(1:IT)
        ENDIF
        IF(HLTAB(J).GT.0) THEN
          TEMPC = GET_TABID(HLTAB(J))
          IT = LENSTR(TEMPC)
          HLTABID = ' '
          HLTABID(HLIDLEN-IT+1:HLIDLEN) = TEMPC(1:IT)
        ELSE
          HLTABID = ' '
        ENDIF
        IF(NODEID.EQ.'NO') THEN
          WRITE(STDOUT,*) ' *BUG*: Should not get here: CLVIN'
        ELSE
          WRITE(STDOUT,18) NODE, BNODID(J), TABID(1:IDLEN), XVEC(J),
     A                       ZBVEC(J), KA(J), KD(J), HLTABID(1:HLIDLEN)
C          WRITE(STDOUT,18) NODE, BNODID(J), XTABN, XVEC(J), ZBVEC(J),
C     A                   KA(J), KD(J), HLTAB(J)
        ENDIF
        XVEC(J) = SFAC*XVEC(J)
        NODE = NODE + 1
 4000 CONTINUE
 
C     CHECK TABLE NUMBERS AND COMPLETE ANY CROSS SECTION INTERPOLATION
C     REQUESTS
      DO 5000 J=FIRST,LAST
        XTAB = NSEC(J)
        IF(XTAB.GT.0) THEN
          IF(TYPFLG.EQ.1) THEN
C           We must have at least the values in old type 12(new type 22)
            CALL CHKTAB
     I                 (12, STDOUT, FTPNT, MFTNUM,
     M                  XTAB,
     O                  EFLAG)
C           ON RETURN XTAB IS REPLACED BY THE ADDRESS OF THE TABLE
            NSEC(J) = XTAB
          ELSE
C           Any cross section type will do.
            CALL CHKTAB
     I                 (20, STDOUT, FTPNT, MFTNUM,
     M                  XTAB,
     O                  EFLAG)
C           ON RETURN XTAB IS REPLACED BY THE ADDRESS OF THE TABLE
            NSEC(J) = XTAB
 
          ENDIF
        ENDIF
 5000 CONTINUE
 
C     DO ANY CROSS SECTION INTERPOLATION REQUESTS
 
      IF(EFLAG.EQ.0) THEN
        CALL SCAN
     I           (STDOUT, NBRA, NBN, MFTNUM, MRFTAB, BRPT, NSEC, XVEC,
     I            ZBVEC, SFAC,
     M            EFLAG, FTP, FTKNT)
      ENDIF
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop.  Errors found.' 
      END
C
C
C
      SUBROUTINE   CULVRT
     I                   (STDIN, STDOUT, STDTAB, NFAC,
     M                    TABDIR, EFLAG, FTP, FTKNT)
 
C     + + + PURPOSE + + +
C     Compute a two-d table for flow through one or more culverts
C     including possible flow over the roadway.
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTKNT, FTP, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL NFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
C     FTKNT  - function table counter
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'ftable.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'cdcom.cmn'
      INCLUDE 'appcom.cmn'
      INCLUDE 'xs0com.cmn'
      INCLUDE 'xs1com.cmn'
      INCLUDE 'xs2com.cmn'
      INCLUDE 'xs3com.cmn'
      INCLUDE 'xs4com.cmn'
      INCLUDE 'x43com.cmn'
      INCLUDE 'x44com.cmn'
      INCLUDE 'typlim.cmn'
      INCLUDE 'typtrn.cmn'
      INCLUDE 'rdfcom.cmn'
      INCLUDE 'culcom.cmn'
      INCLUDE 'embcom.cmn'
      INCLUDE 'depcom.cmn'
      INCLUDE 'embwrq.cmn'
      INCLUDE 'flapgate.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'putget.cmn'
      INCLUDE 'tabupgrade.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CULVERT 
      PARAMETER (CULVERT=1)
      
      INTEGER BSHAPE, CHK23, CHKBAR, CONFLG, EXPFLG, FRTYPE, FTPTMP, I,
     A        ID, IFLAG, IU, IT, J, NBN, NBRA, NFRAC, NHU, NT, NTAB, 
     B        Q3VSRD, TAB, TABTYP, TYPE, TYPFLG, WFLAG, Y3VSTW, Y2VSD,
     C        FTPBASE, ftpup, verbose
      INTEGER BRPT(8,1), NEGTAB(MNBN), DEFAULT_TABVEC(13)
      REAL A1TRUE, ALP1T, CCT5, CL34, DDN, DDROP, DE14, DE34, DEPCON,
     A     DF, DFDC, DFDR, DQ, DROP, DUP, DY2, DZ4, ED, EU, FDRDW, 
     B     FDROP, FDVEC(PMXNHU), FREED, H, HCREST, HDATUM, HHLIM, HP, 
     C     HRDFD, HRDFU, HUOLD, HUP, HUVEC(PMXNHU), K1TRUE, 
     D     PFDVEC(PMXFRC), POWER, QC3, QCLV, QFREE, QMAT(PMXNHU,PMXFRC),
     E     QOLD, RATIO, RDFLOW, SFAC, X3D, X6D, Y3LIMU, YEXIT, 
     F     Z1TRUE, Z3PSAV, Z3SAV, Z3SOF, Z44SAV, Z4OLD, Z4SAV, 
     G     ZDATUM, ZSBRDF, ZT, Y2MAT(PMXNHU,PMXFRC), zrhufd

      REAL*8 FD, easting, northing
      CHARACTER BNODID(MNBN)*8, C5VEC(3)*5, CHAR10*10, CHAR4*4, CHAR5*5,
     A          CHAR6*6, CHAR7*7, CL34C*5, CULCLS*8, HEAD*80, LABEL*50,
     B          LINE*80, LINE120*120, LOSOPT*8, NODEID*4, SAVEIT*4, 
     C          CQFREE*7, CWFRD*7, CBOTH*7, CQCLV*7, TABID*16,
     D          APPTABID*16, DEPTABID*16, BEGTABID*16, CHAR16*16,
     e          zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, DBLE, FLOAT, MAX, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FDFRFC, GETD
      INTEGER LENSTR
      CHARACTER GETTOK*8, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKCFC, CHKDEP, CHKTAB, CHKTEL, CLVIN, EMBSUB, ETABIN,
     A         FDCD46, FDFRFC, FHHTYP, FRFCLV, FRFEMB, FRFT7, GETD,
     B         GETTOK, inline, INPRDP, INSERT, LKTA, LKTAB, LKTK, QVSTW,
     C         TABCHK, TDLK10, TWDOUT, XLKT22, XLKTAL, VAR_DECIMAL,
     D         LENSTR, TWOD13PUT, TWOD13GET, GET_INTERNAL_TAB_NUMBER,
     E         STRIP_L_BLANKS, READ_TABID_PLUS, GET_TABID, READ_TABID
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,1X,A80)
 2    FORMAT(A5,1X,A)
 4    FORMAT(A6,2X,I5,2A5)
 6    FORMAT(A6,1X,F10.0)
 7    FORMAT(A6,1X,F10.0,A10,2A5)
 10   FORMAT(A6,1X,A8)
 12   FORMAT(A6,1X,A8)
 14   FORMAT(3X,1X,F10.0,1X,A4)
 15   FORMAT(A7,1X,F10.0)
 16   FORMAT(A5,1X,F10.0)
 17   FORMAT(8X,F10.0)
 18   FORMAT(A5,1X,F10.0)
 19   FORMAT(9X,F10.0)
 20   FORMAT(A6,1X,A4)
 21   FORMAT(7X,F10.0)
 22   FORMAT(A4,1X,F10.0)
 23   FORMAT(7X,F10.3)
 24   FORMAT(A80)
 26   FORMAT(A5,1X,I5)
 27   FORMAT(A4,1X,I5,1X,9X,F5.0)
1000  FORMAT(6X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 3    FORMAT(/,' Adding upstream head=',F8.3,' at the nominal',
     A        ' lower bound for high head flow.')
 5    FORMAT(/,' Adding upstream head=',F8.3,' at the minimum',
     A   ' crest elevation of roadway.')
 8    FORMAT(/,' *WRN:584* KWING=',F5.3,' ignored.  Applies to BOX',
     A       ' culverts only.')
 9    FORMAT(/,' *WRN:585* KPROJ=',F5.3,' ignored.  Does not apply',
     A       ' to BOX culverts.')
 11   FORMAT(/,' Power for distributing partial free drops=',F5.1)
 13   FORMAT(/,' *ERR:721* Relative projection < 0.')
 25   FORMAT(/,' Angle of inlet bevel=',F5.1,' degrees.')
 28   FORMAT(/,' Relative projection=',F8.3)
 29   FORMAT(/,' *ERR:722* Angle of inlet bevel < 0 or > 90 degrees.')
 30   FORMAT(/,' Rounding/beveling value=',F5.3)
 31   FORMAT(/,' Checking the critical flow variation',
     A     ' at section 2.')
 32   FORMAT(/,' This cross section may cause computational problems',
     B  ' if critical',/,' flow decreases with depth.  CULVERT',
     C  ' may become confused',/,' when the relationship between',
     D  ' critical flow',/,' and depth is not unique.')
 33   FORMAT(/,' Checking the critical flow variation',
     A    ' at section 3.')
 34   FORMAT(/,' *ERR:703* Number of upstream heads=',I5,' too',
     A       ' small to add',/,11X,' high-head flow limit.')
 35   FORMAT(/,' *ERR:704* Number of upstream heads=',I5,' too',
     A       ' small to add',/,11X,' head at roadway crest.')
 36   FORMAT(/,' Wingwall angle=',F5.0,' degrees.')
 37   FORMAT(/,' Ratio of flow depth to vertical diameter at culvert',
     A      ' exit',/,' that causes submergence=',F5.2)
 38   FORMAT(/,'  Maximum ratio of depth to vertical diameter',
     A      ' for Type 1 flow=',F5.3)
 39   FORMAT(/,' Maximum ratio of head to vertical diameter for',
     A        ' type 1 and 2 flow=',F5.3)
 40   FORMAT(
     A /,' Processing cannot continue on the current culvert because',
     B /,' it has either a break in slope or a non-uniform',
     C /,' non-prismatic shape.  The automatic procedure to insert',
     D /,' nodes will not function on culverts like this.  The CULVERT',
     E /,' command expects the culvert barrel to have two nodes',
     F /,' placed at special locations close to its inlet.  The first',
     G /,' special node must be at a distance of 3 times the vertical',
     H /,' diameter of the barrel from its inlet.   The second special',
     I /,' node must be at a distance of 6 times the  vertical',
     J /,' diameter from the inlet.  These special nodes are used',
     K /,' in computing Type 5 flow and must be present.  Culverts too',
     L /,' short for the placement of these special nodes cannot be',
     M /,' analyzed if their barrels have a break in slope or if the',
     N /,' barrel is non-uniform non-prismatic.  Breaks in slope must',
     O /,' be used with caution because the analysis procedure does',
     P /,' not find control points within the culvert barrel.')
 41   FORMAT(/,' Flows for water levels at section 4.',
     A       /,'    Road submergence at section 43.',
     B     //,' Upstream Head=',F9.4,'  Free Drop=',F9.4)
 42   FORMAT('  Processing CULVERT TabId= ',A)
 43   FORMAT(/,
     A1X,  'Partial Depth Head  Head  Head  Drop  Flow   Flow   Total ',
     B  ' Loss Energy Energy',
     C/,1X,'free    at    at    at    at    sec   in     over   flow  ',
     D  ' coef loss   loss  ',
     E/,1X,'drop    sec 3 sec 3 sec44 sec 4 1->4  barrel road         ',
     F  ' 3->4 3->4   1->4  ',
     G/,1X,'------- ----- ----- ----- ----- ----- ------ ------ ------',
     H  ' ---- ------ ------')
 44   FORMAT(
     A1X,  '----------------------------------------------------------',
     B  '--------------------')
 45   FORMAT(1X,F7.4,4F6.2,F6.3,3A7,A5,A7,F7.3)
 46   FORMAT(/,' ERR:705* Rounding/beveling value < 0.0 or > 0.14')
 47   FORMAT(/,' ERR:706* Wingwall angle < 0.0 or > 90.0 degrees.')
 48   FORMAT(1X,'Note: Barrel flows full at section 3 when depth',
     A  ' there=',F8.3)
 49   FORMAT(/,1X,'*WRN:**** Flow Type 7 indicates that the downstream',
     A ' section',/,11X,'of the departure reach drowns all other',
     B ' control points. This is',/,11X,'non-standard but not',
     C ' unusual.  However, Type 7 flow requires',/,11X,'several',
     D ' unverified assumptions. Please check results carefully.')
 50   FORMAT(/,' TABID= ',A,2X,A)
 51   FORMAT(/,' *ERR:707* Type 5 submergence ratio < 0.0 or > 1.0')
 52   FORMAT(/,' ',A5,'=',A)
 53   FORMAT(/,' *ERR:708* Ratio of depth to vertical diameter <= 0.5',
     A       ' or >= 1.0')
 54   FORMAT(/,' ',' APPTAB= ',A)
 55   FORMAT(/,' *ERR:709* Ratio of head to vertical diameter < 1.0 or',
     A        ' >= 1.5')
 56   FORMAT(' ',A6,'=',F10.2)
 57   FORMAT(' ',A6,'=',F10.2,A10,2A5)
 58   FORMAT(' ',A6,'=',F10.1)
 59   FORMAT(' BEGTAB=',I5)
 60   FORMAT(' ',A6,'=',A8)
 61   FORMAT(/,' ',A6,'=',A8)
 62   FORMAT(/,' Culvert barrel length <= 3 times its vertical',
     A   ' diameter.  Culvert losses',/,'  will be overestimated.')
 63   FORMAT(' *ERR:603* ROADWAY MOMENTUM FLUX FACTOR=',F6.3,
     A         ' < 0 OR > 1.')
 64   FORMAT(/,' Discharge coefficient for flow types 4 and 6=',F5.3,/,
     A       12X, ' Type 6 jet assumed supported.')
 65   FORMAT(/,1X,'Departure reach start elev.=',F10.2,/,
     A         1X,'Departure reach end elev.=',F10.2)
 66   FORMAT(/,' Wingwall angle adjustment factor=',F5.3)
 67   FORMAT(/,' Discharge coefficient for flow types 4 and 6=',F5.3,/,
     A       12X, ' Type 6 jet assumed unsupported.')
 68   FORMAT(/,' Projecting-pipe adjustment factor=',F5.3)
 69   FORMAT(/,' Rounding/beveling adjustment factor=',F5.3)
 70   FORMAT(' *ERR:710* At least one submergence table is unknown.',
     A       ' Both must be known.')
 71   FORMAT(/,'  Table type 5 replaced by type 13.')
 72   FORMAT(' *ERR:578* APPROACH REACH LENGTH=',F10.2,' < 0.')
 73   FORMAT(' *ERR:579* APPROACH LOSS COEF. < 0 OR > 1.')
 74   FORMAT(/,' *WRN:582* The culvert barrel has slope breaks.  The',
     A    ' CULVERT command may fail.')
 75   FORMAT(' *ERR:581* APPROACH EXPANSION COEF. < 0 OR > 1.')
 76   FORMAT(' *ERR:582* ROUNDING/BEVELING FACTOR < 1 OR > 1.5.')
 77   FORMAT(' *ERR:583* DEPARTURE CONTRACTION COEF. < 0 OR > 1.')
 78   FORMAT(' *ERR:584* WINGWALL FACTOR < 1 OR > 1.25.')
 79   FORMAT(' *ERR:604* WIDTH FACTOR FOR DEPARTURE REACH <= 1.0')
 80   FORMAT(' *ERR:585* PROJECTION ADJUSTMENT < 0.90 OR > 1.')
 81   FORMAT(' *ERR:605* BEGINNING ELEV. FOR DEP. RCH.=',F10.2,
     A     ' >',' CLVRT EXIT ELEV.=',F10.2)
 82   FORMAT(' *ERR:586* LOSS OPTION: ',A8,' INVALID.')
 83   FORMAT(/,' Culvert barrel length <= 6 times its vertical',
     A   ' diameter.  Culvert losses',/,'  may be overestimated.')
 84   FORMAT(' *ERR:587* CULVERT CLASS:', A8,' INVALID.')
 85   FORMAT(' ROADWAY MOMENTUM FLUX FACTOR=',F6.3)
 86   FORMAT(' *ERR:588* TYPE 4-6 COEF.:',F5.2,'  < 0.60 OR > 1.')
 87   FORMAT(' *ERR:607* TABLE# <= 0')
 88   FORMAT(' ',A6,'=',A4)
 89   FORMAT(' *ERR:608* SUBMERGED ROAD FLOW WITH FREE CULVERT FLOW',
     A       ' NOT YET SUPPORTED.')
 90   FORMAT(' ',A4,'=',F10.2)
 91   FORMAT(' *WRN:547* Road is lower than a culvert soffit',
     A      ' elevation.')
 92   FORMAT(/,' ',A80)
 93   FORMAT(' *ERR:609* Road crest elev. < head datum invalid.')
 94   FORMAT(/,' Number of partial free drops=',I5)
 95   FORMAT(/,' ','Input complete. Begin computations.',//,
     A       ' ', ' Upstream opening=',F8.4,' Downstream opening=',
     B       F8.4)
 96   FORMAT(/,' ','Datum for heads is:',F10.4)
 97   FORMAT(2X,A4,'=',I5,' Flap gate factor=',F5.2)
 98   FORMAT(' *ERR:589* TABLE TYPE NOT 6 or 13.')
 99   FORMAT(/,' *ERR:630* Upstream heads non-increasing at:',F10.4)
 101  FORMAT(/,' *ERR:610* Free drop=',F8.4,' < ',F8.4,' not',
     A       ' yet supported.')
102   FORMAT(' ','DEPTAB= ',A,'  BEGTAB= ',A,'  RMFFAC=',F5.0)
1050  FORMAT(/,' FHWA a coefficient=',F5.2)
1051  FORMAT(/,' FHWA c coefficient=',F8.4)
1052  FORMAT(/,' FHWA Y coefficient=',F8.4)
1053  FORMAT(/,' FHWA barrel slope=',F10.4)
C***********************************************************************
      HLCRIT = 0.15
      HLMAX = 0.32
      HLFLAG = 0
 
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE120)
      CALL READ_TABID_PLUS
     I                   (STDOUT, LINE120,
     O                    EFLAG, TABID, TAB, LINE)

      IT = LENSTR(LINE)
      IF(IT.EQ.0) THEN 
        IT = 1
        LINE = ' '
      ENDIF      
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), LINE(1:IT)
      WRITE(*,42) TABID
 
      CALL TABCHK
     I           (STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)
      TABU = TAB


C     Process get/put actions.
      CALL GET_PUTGET_OPTIONS
     I                (STDIN, STDOUT, LINE, CULVERT,
     M                  EFLAG)

      IF(PUTQ.GT.0) THEN
        CALL TABCHK
     I             (STDOUT, PMXTAB,
     M              PUTQ, TABDIR, EFLAG)
      ENDIF
      IF(PUTY2.GT.0) THEN
        CALL TABCHK
     I             (STDOUT, PMXTAB,
     M              PUTY2, TABDIR, EFLAG)
      ENDIF

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,27,ERR=991) CHAR4, TABTYP, FLAP
      WRITE(STDOUT,97) CHAR4, TABTYP, FLAP
      IF(TABTYP.EQ.5) THEN
        WRITE(STDOUT,71)
        TABTYP = 13
      ENDIF
      IF(TABTYP.NE.13.AND.TABTYP.NE.6) THEN
        WRITE(STDOUT,98)
        EFLAG = 1
      ENDIF
 
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)


      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) CHAR5, LABEL
      WRITE(STDOUT,52) CHAR5, LABEL
 
C     INPUT THE APPROACH SECTION DESCRIPTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE, 24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'APPTAB',
     O                EFLAG, APPTABID, APPTAB)
      WRITE(STDOUT,54)  APPTABID(1:LENSTR(APPTABID))
 
      IF(APPTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              APPTAB,
     O              EFLAG)
      ENDIF
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CHAR6, APPELV
 
      WRITE(STDOUT,56) CHAR6, APPELV
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CHAR6, APPLEN
      WRITE(STDOUT,58) CHAR6, APPLEN
 
      IF(APPLEN.LT.0.0)  THEN
        WRITE(STDOUT,72) APPLEN
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CHAR6, APPLOS
      WRITE(STDOUT,56) CHAR6, APPLOS
      IF(APPLOS.LT.0.0.OR.APPLOS.GT.1.0) THEN
        WRITE(STDOUT, 73)
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CHAR6, APPEXP
      WRITE(STDOUT,56) CHAR6, APPEXP
      IF(APPEXP.LT.0.0.OR.APPEXP.GT.1.0) THEN
        WRITE(STDOUT, 75)
        EFLAG = 1
      ENDIF
 
 
C     INPUT THE CULVERT DESCRIPTION.  TREAT LIKE A BRANCH IN FEQ AND
C     MAKE AS CLOSE TO A BRANCH AS POSSIBLE.
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,20,ERR=991) CHAR6,NODEID
      NODEID = GETTOK(NODEID)
      WRITE(STDOUT,88) CHAR6,NODEID
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,22,ERR=991) CHAR4, SFAC
      WRITE(STDOUT,90) CHAR4, SFAC
 
      NBRA = 1
C     TYPFLG at 1 requests automatic interpolation for cross sections
C     along the barrel if the user only gives 2 crossections: one at the
C     barrel entrance and one at the barrel exit.  NEGTAB is not used
C     in culvert but it is used in XSINTERP where CLVIN is called as
C     well. 
      TYPFLG = 1
      CALL CLVIN
     I          (SFAC, STDIN, STDOUT, MNBN, NBRA, NODEID, MFTNUM,
     I           MRFTAB, FTPNT, TYPFLG,
     M           EFLAG, FTKNT, FTP,
     O           NBN, BRPT, NSEC, XVEC, ZBVEC, KA, KD, HLTAB, BNODID,
     O           NEGTAB, CHKBAR, IAT3D, IAT6D, BSHAPE, sbkind)
 
 
C     INPUT THE CULVERT CLASS.  USED TO DECIDE WHAT RELATIONSHIPS TO USE.
C     CLASS IS ASSIGNED BY THE USER AND IS INDEPENDENT OF THE ACTUAL SHAPE
C     OF THE CULVERT OPENING.  NORMALLY HOWEVER, THE USER SHOULD SELECT
C     A CLASS WHICH IS DESCRIPTIVE OF THE SHAPE.
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,10,ERR=991) CHAR6, CULCLS
      WRITE(STDOUT,61) CHAR6, CULCLS
 
C     SKIP LEADING BLANKS
      CULCLS = GETTOK(CULCLS)
      IF(CULCLS.NE.'BOX'.AND.CULCLS.NE.'PIPE'.AND.CULCLS.NE.'MITER'.
     A    AND.CULCLS.NE.'RCPTG'.AND.CULCLS.NE.'FLARED') THEN
        WRITE(STDOUT,84) CULCLS
        EFLAG = 1
      ENDIF
 
C     INPUT THE DEPARTURE SECTION DATA
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
      CALL GET_DEPARTURE_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
C      CALL inline
C     I          (STDIN, STDOUT,
C     O           LINE)
C      READ(LINE,4,ERR=991)  CHAR6, DEPTAB, CHAR5, C5VEC(1)

      CALL  SET_DEPITM(
     O                 DEPTAB, BEGTAB, RMFFAC,
     M                 EFLAG)

      IF(DEPTAB.GT.0) THEN
        DEPTABID = GET_TABID(DEPTAB)
      ELSE
        DEPTABID = '   '
      ENDIF
      IF(BEGTAB.GT.0) THEN      
        BEGTABID = GET_TABID(BEGTAB)
      ELSE
        BEGTABID = '   '
      ENDIF

      WRITE(STDOUT,102) DEPTABID(1:LENSTR(DEPTABID)), 
     A                 BEGTABID(1:LENSTR(BEGTABID)), RMFFAC
      IF(DEPTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              DEPTAB,
     O              EFLAG)
      ENDIF
      IF(BEGTAB.GT.0) THEN
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              BEGTAB,
     O              EFLAG)
      ELSE
C       MAKE TABLES THE SAME
        BEGTAB = DEPTAB
      ENDIF
 
C     CHECK USER SUPPLIED ROADWAY MOMENTUM FLUX FACTOR
      IF(RMFFAC.LT.0.0.OR.RMFFAC.GT.1.0) THEN
        WRITE(STDOUT,63) RMFFAC
        EFLAG = 1
      ENDIF
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,7,ERR=991) CHAR6, DEPELV, CHAR10, C5VEC(1), C5VEC(2)
      WRITE(STDOUT,57) CHAR6, DEPELV, CHAR10, C5VEC(1), C5VEC(2)
 
      IF(CHAR10.NE.'          ') THEN
C       INPUT THE USER SPECIFIED BEGINNING ELEVATION FOR THE
C       DEPARTURE REACH.
        READ(CHAR10,'(F10.0)') BEGELV
        IF(BEGELV.GT.ZBVEC(NBN)) THEN
          WRITE(STDOUT,81) BEGELV, ZBVEC(NBN)
          EFLAG = 1
        ENDIF
      ELSE
C       IF THE DEPELV > CULVERT EXIT ELEVATION THEN MAKE BEGELV
C       SAME AS CULVERT ELEVATION. OTHERWISE MAKE IT THE SAME AS
C       DEPELV
        IF(DEPELV.GT.ZBVEC(NBN)) THEN
          BEGELV = ZBVEC(NBN)
        ELSE
          BEGELV = DEPELV
        ENDIF
      ENDIF
 
      IF(C5VEC(2).NE.'     ') THEN
C       INPUT USER SPECIFIED CULVERT WIDTH FACTOR
        READ(C5VEC(2),'(F5.0)') WIDFAC
        IF(WIDFAC.LT.1.0) THEN
          WRITE(STDOUT,79)
          EFLAG = 1
          WIDFAC = 1.02
        ENDIF
      ELSE
        WIDFAC = 1.02
      ENDIF
 
      WRITE(STDOUT,65) BEGELV, DEPELV
 
C     Find the elevation datum to use to improve the precision
C     of residual computation.
 
      ZDATUM = MIN(APPELV, ZBVEC(1), ZBVEC(NBN), DEPELV)
      APPELV = DBLE(APPELV) - DBLE(ZDATUM)
      DEPELV = DBLE(DEPELV) - DBLE(ZDATUM)
      BEGELV = DBLE(BEGELV) - DBLE(ZDATUM)
      DO 100 I=1,NBN
        ZBVEC(I) = DBLE(ZBVEC(I)) - DBLE(ZDATUM)
 100  CONTINUE
 
C     INPUT THE EXIT LOSS OPTIONS
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,12,ERR=991) CHAR6, LOSOPT
      IF(LOSOPT.NE.'MOMENTUM') THEN
        WRITE(STDOUT,82) LOSOPT
        EFLAG = 1
      ENDIF
      WRITE(STDOUT,60) CHAR6, LOSOPT
 
      IF(LOSOPT.EQ.'ENERGY') THEN
 
        CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
        READ(LINE,6,ERR=991) CHAR6, DEPCON
        WRITE(STDOUT,56) CHAR6, DEPCON
 
        IF(DEPCON.LT.0.0.OR.DEPCON.GT.1.0) THEN
          WRITE(STDOUT, 77)
          EFLAG = 1
        ENDIF
      ENDIF
C     INPUT THE DISCHARGE COEFFICIENT RELATED FACTORS WHICH ARE
C     FUNCTIONS OF GEOMETRY ONLY.
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE, 24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE, 14,ERR=991)  KRB
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) CHAR5, KWING
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,18,ERR=991) CHAR5, KPROJ
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,14,ERR=991)  C46, CHAR4
      IF(CHAR4.EQ.'    ') THEN
        TY6SUP = 1
      ELSEIF(CHAR4.EQ.'SJET') THEN
        TY6SUP = 1
      ELSE
        TY6SUP = 0
      ENDIF
 
C     Set the default values for type 5 stuff and the
C     internal lookup stuff.
      HHTYPE = 0
      RBVAL = 0.0
      TY5SBF = 0.75
      WWANGL = 0.0
      IT = 9979
      CHAR16 = ' '
C     Clear all the FHWA values.
      AFHWA = 0.0
      CFHWA = 0.0
      YFHWA = 0.0
      SFHWA = 0.0
      DO 250 J=1,13
        IT = IT + 1
        WRITE(CHAR16(1:5),'(I5)') IT
        CALL STRIP_L_BLANKS(
     M                      CHAR16)
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, CHAR16,
     M                               EFLAG,
     O                               DEFAULT_TABVEC(J))
250   CONTINUE

      TB6ADR = DEFAULT_TABVEC(01)
      TB7ADR = DEFAULT_TABVEC(02)
      TB8ADR = DEFAULT_TABVEC(03)
      TB15AD = DEFAULT_TABVEC(04)
      TB16AD(1) = DEFAULT_TABVEC(05)
      TB16AD(2) = DEFAULT_TABVEC(06)
      TB16AD(3) = DEFAULT_TABVEC(07)
      TB16AD(4) = DEFAULT_TABVEC(08)
      TB5ADR = DEFAULT_TABVEC(09)
      TB21ADR = DEFAULT_TABVEC(10)
      TB22ADR = DEFAULT_TABVEC(11)
      TB24ADR = DEFAULT_TABVEC(12)
      TBLPDADR = DEFAULT_TABVEC(13)
      BVANGL = 0.0
      LPOVD = 0.0
 
C     Optional type 5 input may appear here.  If not the
C     roadway heading appears.
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
 
      IF(LINE(1:6).EQ.'TYPE 5'.OR.LINE(1:6).EQ.'Type 5') THEN
C       Type 5 flow parameters appear.
        WRITE(STDOUT,92) HEAD
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,15,ERR=991) CHAR7, RBVAL
        IF(RBVAL.LT.0.0.OR.RBVAL.GT.0.14) THEN
          WRITE(STDOUT,46)
          EFLAG = 1
        ENDIF
        WRITE(STDOUT,30) RBVAL
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,17,ERR=991) BVANGL
        IF(BVANGL.LT.0.0.OR.BVANGL.GT.90.0) THEN
          WRITE(STDOUT,29)
          EFLAG = 1
        ENDIF
        WRITE(STDOUT,25) BVANGL
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,17,ERR=991) WWANGL
        WRITE(STDOUT,36) WWANGL
        IF(WWANGL.LT.0.0.OR.WWANGL.GT.90.0) THEN
          WRITE(STDOUT,47)
          EFLAG = 1
        ENDIF
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,17,ERR=991) LPOVD
        WRITE(STDOUT,28) LPOVD
        IF(LPOVD.LT.0.0) THEN
          WRITE(STDOUT,13)
          EFLAG = 1
        ENDIF
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,19,ERR=991) TY5SBF
        WRITE(STDOUT,37) TY5SBF
        IF(TY5SBF.LE.0.0.OR.TY5SBF.GT.1.0) THEN
          WRITE(STDOUT,51)
          EFLAG = 1
        ENDIF

        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
      ENDIF

C     Insert input of FHWA information after the Type 5 block
      IF(LINE(1:4).EQ.'FHWA') THEN
C       Process the FHWA parameters.
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,1000,ERR=991) AFHWA
C       If AFHWA is zero we will assign a value based on 
C       culvert class.  Otherwise retain the given value. 
        IF(AFHWA.EQ.0.0) THEN
          IF(CULCLS.EQ.'MITER') THEN
            AFHWA = -0.7
          ELSE
            AFHWA = 0.5
          ENDIF
        ENDIF
        WRITE(STDOUT,1050) AFHWA
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,1000,ERR=991) CFHWA
        WRITE(STDOUT,1051) CFHWA
C       Apply factor of 2g to CFHWA
        CFHWA = GRAV2*CFHWA
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,1000,ERR=991) YFHWA
        WRITE(STDOUT,1052) YFHWA


C       Compute the bottom slope of the culvert near the 
C       entrance. Slope > 0 for decline.
        SFHWA = (ZBVEC(1) - ZBVEC(2))/(XVEC(2) - XVEC(1))
        WRITE(STDOUT,1053) SFHWA
C       Read next line to prepare for what follows 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)

      ENDIF
      IF(LINE(1:5).EQ.'TABLE'.OR.LINE(1:5).EQ.'Table') THEN
C       Get the next line and read user assigned table numbers from it.
        WRITE(STDOUT,92) LINE(1:80)
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,*,ERR=991)  TB6ADR, TB7ADR, TB8ADR, TB15AD,
     A      TB16AD(1), TB16AD(2), TB16AD(3), TB16AD(4),
     B      TB5ADR, TB21ADR, TB22ADR, TB24ADR, TBLPDADR
 
C       Transform to internal table number

C       Read the next line.  Should be a heading of
C       some sort.
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
      ENDIF
 
C     Get the function table addresses.
      CALL CHKTAB
     I           (10, STDOUT, FTPNT, MFTNUM,
     M            TB6ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (10, STDOUT, FTPNT, MFTNUM,
     M            TB7ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TB8ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (10, STDOUT, FTPNT, MFTNUM,
     M            TB15AD,
     O            EFLAG)
      DO 120 I=1,4
        CALL CHKTAB
     I             (10, STDOUT, FTPNT, MFTNUM,
     M              TB16AD(I),
     O              EFLAG)
 120  CONTINUE
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TB5ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TB21ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (10, STDOUT, FTPNT, MFTNUM,
     M            TB22ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TB24ADR,
     O            EFLAG)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TBLPDADR,
     O            EFLAG)
 
C     Set defaults for low-head flow limits.
      TY1YTD = 0.95
      TY1HTD = 1.4
 
      IF(LINE(1:6).EQ.'TYPE 1'.OR.LINE(1:6).EQ.'Type 1') THEN
C       Process optional type 1 parameters.
        WRITE(STDOUT,92) LINE(1:80)
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,21,ERR=991) TY1YTD
        WRITE(STDOUT,38) TY1YTD
        IF(TY1YTD.LE.0.5.OR.TY1YTD.GT.1.0) THEN
          WRITE(STDOUT,53)
          EFLAG = 1
        ENDIF
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,23,ERR=991) TY1HTD
        WRITE(STDOUT,39) TY1HTD
        IF(TY1HTD.LT.1.0.OR.TY1HTD.GE.1.6) THEN
          WRITE(STDOUT,55)
          EFLAG = 1
        ENDIF
C       Read the next heading
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
      ENDIF
 
      IF(EFLAG.EQ.0) THEN
C       Check for request for lookup.  Assumes that type 5 section
C       has been given with the parameter values.  The lookup is
C       done anyway because all values in that section have defaults.
        IF(KRB.EQ.0.0) THEN
          IF(BVANGL.EQ.0.0) THEN
C           RBVAL is for rounding.
            CALL LKTAB
     I                (TB21ADR, RBVAL, 1,
     O                 KRB, NTAB, DF)
          ELSE
C           RBVAL is for beveling.  Two-D table.
            CALL TDLK10
     I                 (STDOUT, TB22ADR, 10, RBVAL, BVANGL,
     O                  KRB, DFDR, DFDC)
        
          ENDIF
        ENDIF
        WRITE(STDOUT,69) KRB
        IF(KRB.LT.1.0.OR.KRB.GT.1.5) THEN
          WRITE(STDOUT,76)
          EFLAG = 1
        ENDIF
        
        
        IF(KWING.EQ.0.0) THEN
          IF(CULCLS.EQ.'BOX') THEN
            CALL LKTAB
     I                (TB24ADR, WWANGL, 1,
     O                 KWING, NTAB, DF)
          ELSE
C           Only box culverts are affected by wingwalls.  All others
C           are not changed.
            KWING = 1.0
          ENDIF
        ENDIF
        WRITE(STDOUT,66) KWING
        IF(KWING.LT.1.0.OR.KWING.GT.1.25) THEN
          WRITE(STDOUT,78)
          EFLAG = 1
        ENDIF
        
        IF(KWING.NE.1.0.AND.CULCLS.NE.'BOX') THEN
          WRITE(STDOUT,8) KWING
        ENDIF
        IF(KPROJ.EQ.0.0) THEN
          IF(CULCLS.NE.'BOX'.AND.CULCLS.NE.'RCPTG') THEN
            CALL LKTAB
     I                (TBLPDADR, LPOVD, 1,
     O                 KPROJ, NTAB, DF)
          ELSE
            KPROJ = 1.0
          ENDIF
        ENDIF
        WRITE(STDOUT,68) KPROJ
        IF(KPROJ.LT.0.90.OR.KPROJ.GT.1.0) THEN
          WRITE(STDOUT, 80) KPROJ
          EFLAG = 1
        ENDIF
        IF(KPROJ.NE.1.0.AND.
     A                 (CULCLS.EQ.'BOX'.OR.CULCLS.EQ.'RCPTG')) THEN
          WRITE(STDOUT,9) KPROJ
        ENDIF
        IF(C46.EQ.0.0) THEN
          CALL FDCD46
     I               (STDOUT, CULCLS, RBVAL,
     O                C46)
        ENDIF
        IF(TY6SUP.EQ.1) THEN
          WRITE(STDOUT,64) C46
        ELSE
          WRITE(STDOUT,67) C46
        ENDIF
        IF(C46.LT.0.60.OR.C46.GT.1.0) THEN
          WRITE(STDOUT,86) C46
          EFLAG = 1
        ENDIF
      ENDIF  
C     INPUT THE ROADWAY DESCRIPTION HERE.  STORE IN SAME MANNER AS
C     FOR EMBANKQ.
 
C     This LINE has been read above in process various optional 
C     input blocks. 
      WRITE(STDOUT,92) LINE(1:80)
 
      CALL ETABIN
     I           (STDIN, STDOUT, MFTNUM, FTPNT,
     O            EFLAG, PLCWTB, GLCWTB, PHCWTB, GHCWTB, PSUBTB, GSUBTB)
 
      IF(PSUBTB.EQ.0.OR.GSUBTB.EQ.0) THEN
        WRITE(STDOUT,70)
        EFLAG = 1
      ENDIF
      CALL INPRDP
     I           (STDIN, STDOUT, ZDATUM,
     O            EFLAG, NOFF, MINCRS, MINLOC, OFF, CREST, WIDTH,
     O            APPROC, SURF)
 
 
C     INPUT THE FACTORS CONTROLLING THE UPSTREAM HEADS AND THE
C     DISTRIBUTION OF DROPS
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,26,ERR=991) CHAR5, NFRAC
      WRITE(STDOUT,94) NFRAC
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) CHAR5, POWER
      WRITE(STDOUT,11) POWER
 
C     COMPUTE THE PROPORTIONS OF FREE DROP
 
      DO 200 I=1,NFRAC
        PFDVEC(I) = (FLOAT(I-1)/FLOAT(NFRAC-1))**POWER
 200  CONTINUE
 
C     INPUT THE HEAD SEQUENCE
 
      I = 1
      HUOLD = -1.0
 300  CONTINUE
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
        READ(LINE,'(F10.0)',ERR=991) HUVEC(I)
        WRITE(STDOUT,'(1X,F10.2)') HUVEC(I)
        IF(HUVEC(I).LE.0.0) THEN
          NHU = I - 1
          GOTO 310
        ELSE
          IF(HUVEC(I).LE.HUOLD) THEN
            WRITE(STDOUT,99) HUVEC(I)
            EFLAG = 1
          ENDIF
          HUOLD = HUVEC(I)
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
 
C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.
 
C     ASSIGN LABELED COMMON COPIES OF TABLE ADDRESSES AND
C     BOTTOM ELEVATIONS
 
C     FOR CURRENT VERSION SECTION 0 AND SECTION 1 ARE THE SAME.
C     THIS WILL BE CHANGED WHEN DROP INLETS ARE SUPPORTED.
 
      ADRXS0 = APPTAB
      ADRXS1 = APPTAB
      ADRXS2 = NSEC(1)
      ADRXS3 = NSEC(NBN)
      ADRXS4 = DEPTAB
      ADRS43 = BEGTAB
      ADRS44 = BEGTAB
 
 
      ZB0 = APPELV
      ZB1 = APPELV
      ZB2 = ZBVEC(1)
      ZB3 = ZBVEC(NBN)
      ZB4 = DEPELV
      ZB43 = BEGELV
      ZB44 = BEGELV
 
C     VFAC WILL BE CHANGED WHEN DROP INLETS ARE SUPPORTED
      VFAC = 1.0
 
C     FIND THE MAXIMUM VERTICAL DIAMETER ALONG THE CULVERT
 
      DO 320 I=1,NBN
        DVEC(I) = GETD(NSEC(I), STDOUT)
        IF(DVEC(I).LT.0) THEN
C         Signal an error.
          EFLAG = 1
        ENDIF
 320  CONTINUE
      DUP = DVEC(1)
      DDN = DVEC(NBN)
 
      WRITE(STDOUT,95) DUP, DDN
 
C     Check the barrel description for non-increasing critical
C     flow.
      WRITE(STDOUT, 31)
      CALL CHKCFC
     I           (GRAV, STDOUT, ADRXS2,
     O            WFLAG)
      IF(WFLAG.NE.0) THEN
        WRITE(STDOUT,32)
      ENDIF
      IF(ADRXS3.NE.ADRXS2) THEN
        WRITE(STDOUT,33)
        CALL CHKCFC
     I             (GRAV, STDOUT, ADRXS3,
     O              WFLAG)
        IF(WFLAG.NE.0) THEN
          WRITE(STDOUT,32)
        ENDIF
      ENDIF
C     Set the elevation of the soffit at culvert exit(section 3)
      Z3SOF = ZB3 + DDN
 
 
C     SELECT THE POINT DEFINING HEAD FOR THE CULVERT.
 
      HDATUM = MAX(APPELV, ZBVEC(1), ZBVEC(NBN), DEPELV)
 
      WRITE(STDOUT,96) HDATUM + ZDATUM
 
      IF(MINCRS.LE.HDATUM) THEN
C       MINIMUM ROADWAY ELEVATION MUST BE GREATER THAN THE ELEVATION
C       OF THE DATUM FOR HEAD
 
        WRITE(STDOUT,93)
        EFLAG = 1
        RETURN
      ENDIF
 
      IF(MINCRS.LE.ZBVEC(1) + DUP.OR.MINCRS.LE.ZBVEC(NBN) + DDN) THEN
C       SEND WARNING IF MINIMUM ELEVATION OF ROADWAY IS BELOW
C       THE CULVERT SOFFITS.
        WRITE(STDOUT,91)
      ENDIF
 
      CALL XLKTAL
     I           (ADRXS3,
     M            DDN,
     O            A3FULL, T3, DT3, J3, K3, DK3, BT3ATD, DBET3, AP3ATD,
     O            DALP3)
      CALL LKTA
     I         (ADRXS2,
     M          DUP,
     O          A2FULL)
      CALL LKTK
     I         (ADRXS2,
     M          DUP,
     O          K2FULL)
 
C     Compute the length of the barrel and length related values.
      L23 = ABS(XVEC(1) - XVEC(NBN))
      LOVERD = L23/DUP
 
C     Check user-specified node spacing for inclusion of the
C     critical locations for type 5 flow.
      IF(CHKBAR.EQ.1) THEN
        X3D = 3.*DUP
        X6D = 6.*DUP
 
        IAT3D = -1
        IAT6D = -1
        DO 322 I=2,NBN
          RATIO = ABS(XVEC(I) - XVEC(1))/DUP
          IF(RATIO.GE.2.9.AND.RATIO.LE.3.1) THEN
            IAT3D = I
          ELSEIF(RATIO.GE.5.9.AND.RATIO.LE.6.1) THEN
            IAT6D = I
          ENDIF
 322    CONTINUE
 
        IF(L23.LE.X3D) THEN
C         Culvert length is less than 3 times its vertical diameter.
          IAT3D = 1
          IAT6D = 0
          WRITE(STDOUT,62)
        ELSEIF(L23.LE.X6D) THEN
          WRITE(STDOUT,83)
          IF(IAT3D.GT.0) THEN
            IAT6D = 0
          ENDIF
        ENDIF
        IF(IAT3D.LT.0.OR.IAT6D.LT.0) THEN
C         Cannot continue
          EFLAG = 1
          WRITE(STDOUT,40)
        ENDIF
      ENDIF
 
      IF(SBKIND.EQ.0) THEN
        SZERO = (ZB2 - ZB3)/L23
      ELSE
C       There are breaks in slope.  Take the bottom slope
C       at the entrance.
        SZERO = (ZBVEC(1) - ZBVEC(2))/ABS(XVEC(1) - XVEC(2))
        WRITE(STDOUT,74)
      ENDIF
 
C     Is the barrel prismatic?
      IF(BSHAPE.EQ.1) THEN
C       Non-prismatic.  Compute the factor on square of flow
C       to estimate the friction and eddy losses in the
C       culvert barrel.
        FRCFAC = FDFRFC(1, NBN)
      ELSE
        FRCFAC = L23/K2FULL**2
      ENDIF
 
      IF(EFLAG.NE.0) RETURN
 
C     Add special head values to the head vector.  One is at the
C     nominal boundary between high-head and low-head flow.  The
C     other is at the crest of the roadway embankment if that crest
C     is within the existing head range.  This latter rule is used
C     because false high crests are used to disable computation of
C     the flow over the embankment.
 
      HHLIM = 1.5*DUP + ZB2 - HDATUM
      IF(HHLIM.LE.0.0) HHLIM = -1.0
      IF(HHLIM.GE.HUVEC(NHU)) THEN
C       Do not add the high-head flow limit
        HHLIM = -1.0
      ENDIF
      HCREST = MINCRS - HDATUM
      IF(HCREST.GE.HUVEC(NHU)) THEN
C       Do not add the road crest head
        HCREST = -1.0
      ENDIF
      IF(HHLIM.GT.0.0) THEN
        NT = NHU
        CALL INSERT
     I             (HHLIM, EPSDIF, PMXNHU,
     M              NHU, HUVEC,
     O              EFLAG)
        IF(EFLAG.GT.0) THEN
          WRITE(STDOUT,34) PMXNHU
          RETURN
        ENDIF
        IF(NHU.GT.NT) THEN
          WRITE(STDOUT,3) HHLIM
        ENDIF
      ENDIF
      IF(HCREST.GT.0.0) THEN
        NT = NHU
        CALL INSERT
     I             (HCREST, EPSDIF, PMXNHU,
     M              NHU, HUVEC,
     O              EFLAG)
        IF(EFLAG.GT.0) THEN
          WRITE(STDOUT,35) PMXNHU
          RETURN
        ENDIF
        IF(NHU.GT.NT) THEN
          WRITE(STDOUT,5) HCREST
        ENDIF
      ENDIF
C     Establish the function table storage pointer for any temporary
C     tables to be constructed during the CULVERT command.
 
      FTPTMP = FTP
      
 
C     CHECK TO MAKE SURE THAT THE CROSS SECTION AT THE BEGINNING OF THE
C     DEPARTURE REACH IS LARGER THAN THE CULVERT EXIT SECTION BY
C     AT LEAST THE FACTOR WIDFAC WHEN THE WATER SURFACE IN THE TWO
C     SECTIONS IS AT THE SAME ELEVATION.  If not true, construct
C     a temporary table that does satisfy the requirements.
C     Put out graphic warning message to user.
 
      CALL CHKDEP
     I           (STDOUT,
     M            EFLAG, FTPTMP)
 
      IF(EFLAG.NE.0) RETURN
 
C     SET THE CHARACTER OF THE SECTION AT LOCATION 44.  Always the
C     same as the section at 43.    Section 4 may differ from these
C     two, however,
 
      ADRS44 = ADRS43
 
C     Signal that the boundaries of types 1, 2, 5, 61, and 6 have not
C     been computed.
 
      CD1 = 0.0
      Z1TY1 = ZB1
      Z1TY1L = ZB1
      CD2 = 0.0
      Z1TY2 = ZB1
      CC5 = 0.0
      Z1TY5 = ZB1
      CD61 = 0.0
      Z1TY61 = ZB1
      CD6 = 0.0
      Z1TY6 = ZB1
 
      Z1TY51 = ZB1
      Z1TY52 = ZB1
 
      Z3PEST = ZB3

C     Clear the transitional free flow coefficients of discharge.
      CD1T6 = 0.0
      CD2T6 = 0.0
      CD61T6 = 0.0
      CC2T5 = 0.0
      CC1T5 = 0.0
 
C     Set the user full-flow-inducing value at section 3 for flow type 5
C     and its relatives.
      Y3LIMU = DDN*TY5SBF
 
C     Initialize the flow type memories: last free type and last
C     submerged type.  -1 means undefined.
      LFTYPE = -1
      LSTYPE = -1
 
C     Signal that no value of free flow is known.
 
      Q3FREE = 0.0
      Y3FREE = 0.0
      Y2FREE = 0.0
 
      IU = 1
      ID = NBN
C      IF(TABTYP.EQ.14) THEN
C       Compute a 2-D table of type 14.
C        CALL DOTY14(STDOUT, EFLAG, TAB, TABTYP, HDATUM, ZDATUM,
C     A              DUP, DDN, NFAC, IU, ID, CULCLS, FTPTMP, NHU, HUVEC)
C        RETURN
C      ENDIF
C     FOR EACH UPSTREAM HEAD IN THE LIST, COMPUTE THE FREE FLOW AND
C     THE SUBMERGED FLOWS AT THE INDICATED PARTIAL FREE FLOW DROPS.
 
      FTPBASE = FTPTMP
      DO 2000 I=1,NHU
 
        FTPTMP = FTPBASE
C       CLEAR THE FLAG USED TO LIMIT THE NUMBER OF ERROR MESSAGES
C       WRITTEN FOR INVALID HEADWATER RATIO TO ONE
 
        RATFLG = 0
 
C       CLEAR FLAG FOR HEAD TO EMBANKMENT WIDTH RATIO WARNING
C       MESSAGES
 
        HLFLAG = 0
 
C       CLEAR FLAG FOR PIEZOMETRIC HEAD TO WEIR HEIGHT RATIO
 
        HPFLAG = 0
 
        HUP = HUVEC(I)
 
C       COMPUTE THE ELEVATION IN SECTION 1(THE APPROACH SECTION)
 
        Z1 = HDATUM + HUP
 
C       FIND THE ELEMENTS AT SECTION 1
 
        Y1 = HDATUM + HUP - ZB1
        CALL XLKTAL
     I             (ADRXS1,
     M              Y1,
     O              A1, T1, DT1, J1, K1, DK1, BET1, DBET1, ALP1, DALP1)
        A1TRUE = A1
        ALP1T = ALP1
        K1TRUE = K1
        Z1 = ZB1 + Y1
        EU = Z1
        Z1TRUE = Z1
        YUPTRU = Y1
 
 
 
C        WRITE(STDOUT,*) ' CULVERT: Z1=',Z1
C       COMPUTE FREE FLOW OVER THE ROADWAY IF ANY.  NOTE THAT WE
C       ASSUME THAT CHANGES TO THE  VELOCITY HEAD OF APPROACH
C       FOR THE ROADWAY CAUSED BY FLOW THROUGH THE CULVERTS CAN
C       BE NEGLECTED.  THIS CAN BE MODIFIED BUT MUST WAIT DEVELOPMENT
C       OF A METHOD FOR DISTRIBUTING THE CULVERT FLOW ACROSS THE
C       APPROACH SECTION TO BE CONSISTENT WITH THE METHODS USED
C       TO COMPUTE THE FLOW OVER THE ROADWAY.
 
        IF(Z1TRUE.GT.MINCRS) THEN
C         THERE IS FLOW OVER THE ROADWAY. FIND ITS FREE FLOW VALUE AND
C         THE DOWNSTREAM ELEVATION AT THE FREE FLOW BOUNDARY
 
C         FIND THE HEAD ON THE ROADWAY
 
          HRDFU = HDATUM + HUP - MINCRS
 
          CALL FRFEMB
     I               (HRDFU, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB,
     I                NOFF, OFF, CREST, WIDTH, APPROC, SURF, RMFFAC,
     I                HLCRIT, HLMAX,
     M                HLFLAG, HPFLAG,
     O                ZT, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O                QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O                APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O                WFRDF, MFRDF, EFRDF)
 
          CALL EMBSUB
     I               (MINLOC, MINCRS, NOFF, SURF, TOTHL, TOTHR, HRDFU,
     O                HRDFD, FD, ZSBRDF)
          FDRDW = FD 
C         ZSBRDF GIVES THE TAILWATER ELEVATION WHICH MUST BE
C         REACHED FOR SUBMERGENCE OF THE FLOW OVER THE ROADWAY TO OCCUR.
 
 
        ELSE
          HRDFU = 0.0
          WFRDF = 0.0
          MFRDF = 0.0
          EFRDF = 0.0
          ZSBRDF = MINCRS
          FDRDW = 0.0
        ENDIF
 
C       SET THE GENERIC FLUX VALUES FOR THE ROAD
 
        WFRD = WFRDF
        MFRD = MFRDF
        EFRD = EFRDF
 
C        WRITE(STDOUT,*) '  '
C        WRITE(STDOUT,*) ' HUP=',HUP,' HRDFU=',HRDFU
C        WRITE(STDOUT,*) ' WFRDF=',WFRDF,' FDRDW=',FDRDW
C        WRITE(STDOUT,*) ' MFRDF=',MFRDF
C        WRITE(STDOUT,*) ' ZSBRDF=',ZSBRDF
C       Determine what the high head flow type will be.  If the
C       ratio of head on the culvert opening is less than 1.5
C       the high-head flow type is computed at a head ratio
C       of 1.5.  Above 1.5, for rough pipe culverts, the head
C       can change the flow type.  For smooth pipe and for box
C       culverts(assumed to be smooth) the high-head type is
C       independent of the elevation of water at section 1.
 
        IF(HHTYPE.EQ.0) THEN
C         Establish the high-head flow type at the high-head flow
C         limit
C          HP = 1.5001*DUP + ZB2 - HDATUM
          HP = (TY1HTD+.1001)*DUP + ZB2 - HDATUM
          CALL FHHTYP
     I               (STDOUT, CULCLS, NFAC, HP, HDATUM, DUP, IU, ID,
     I                WFRD, Y3LIMU,
     M                Y3LTY5,
     O                HHTYPE, CCT5, YEXIT)
          IF(HHTYPE.EQ.5) THEN
C           Set the elevation at section 1 at the lower limit of
C           type 5 flow.
            Z1TY5 = HP + HDATUM
            CC5 = CCT5
            CALL XLKT22
     I                 (ADRXS3,
     M                  YEXIT,
     O                  A3, T3, DT3, J3, K3, DK3, BET3, DBET3, ALP3,
     O                  DALP3, QC3)
            BT3AT5 = BET3
            AP3AT5 = ALP3
          ENDIF
        ENDIF
 
        IF((HUP+HDATUM-ZB2)/DUP.GE.1.5.AND.HHTYPE.EQ.5) THEN
C         Review high-head type.
          CALL FHHTYP
     I               (STDOUT, CULCLS, NFAC, HUP, HDATUM, DUP, IU, ID,
     I                WFRD, Y3LIMU,
     M                Y3LTY5,
     O                HHTYPE, CCT5, YEXIT)
        ENDIF
 
 
C       COMPUTE FREE FLOW THROUGH THE CULVERT.  THE FREE DROP FOR
C       THE CULVERT MAY BE  INFLUENCED BY THE FLOW OVER THE ROADWAY.
 
C        WRITE(STDOUT,*) ' CULVERT: BEFORE CALL TO FRFCLV Z1=',Z1
 
        CALL FRFCLV
     I             (STDOUT, HDATUM, ZDATUM, HUP, DUP, DDN, IU, ID,
     I              CULCLS,
     O              EFLAG, TYPE, CONFLG, EXPFLG, ZSBRDF, QFREE, FREED)
 
        IF(FREED.LT.EPSDIF.AND.EXPFLG.NE.0) THEN
          WRITE(STDOUT,101) FREED, EPSDIF
          EFLAG = 1
          WRITE(*,*) '   Free drop too small'
          RETURN
        ENDIF
C       Check for drop from section 1 to section 43.
        IF(Z1TRUE - Z43OLD.LT.EPSDIF) THEN
          WRITE(STDOUT,*) ' *ERR:711* Drop from section 1 to ',
     A     ' section 43 < ',EPSDIF,' not yet supported.'
          WRITE(*,*) '  Drop from section 1 to 43 too small.'
          EFLAG = 1
          RETURN
        ENDIF
        IF(EFLAG.NE.0) GOTO 2001
 
        IF(TYPE.NE.0) THEN
          FDROP = MAX(FREED, FDRDW)
        ELSE
          FDROP = FREED
        ENDIF
 
C       Do the energy check on the free flow.
        IF(EXPFLG.NE.0) THEN
C         Enable checking of energy losses involving section 2
          IF(TYPE.EQ.1.OR.TYPE.EQ.2) THEN
            CHK23 = 1
          ELSE
            CHK23 = 0
          ENDIF
          IF(ALPHA3.EQ.0.0) THEN
            DE34 = 0.0
          ELSE
            DE34 = ALPHA3
          ENDIF
          CALL CHKTEL
     I               (STDOUT, CHK23,
     M                DE34,
     O                IFLAG, CL34, DE14)
        ENDIF
 
C       SAVE THE TYPE OF THE FREE FLOW THROUGH THE CULVERT,
C       THE CRITICAL DEPTH AT THE EXIT OF THE CULVERT, AND
C       THE PIEZOMETRIC DEPTH AT THE EXIT OF THE CULVERT.
C       NEEDED TO PROPERLY COMPUTE SUBMERGED FLOW.
 
        FRTYPE = TYPE
        Z3SAV = Z3
        Z3PSAV = Z3P
        Z44SAV = Z44
        Z4SAV = Z4
        RDFLOW = WFRD
C       Save free flow condition at section 2.
C        Y2MAT(I,NFRAC) = Z2 - ZB2
        Y2MAT(I,NFRAC) = Y2

 
C       Clear the addresses for the temporary tables created
C       by QVSTW to make sure that FRFT7 cannot function
C       if they are not computed for the same upstream head.
        Q3VSRD = 0
        Y3VSTW = 0
        Y2VSD = 0
C       Compute the flow in the culvert as a function
C       of tailwater at section 43.  This must also be
C       done for free flow type 7.
        CALL QVSTW
     I            (STDOUT, Z1TRUE, IU, ID, NFRAC, POWER, ZSBRDF, CULCLS,
     I             A1TRUE, ALP1T, K1TRUE, HDATUM, ZDATUM, FRTYPE,
     M             EFLAG, FTPTMP,
     O             Q3VSRD, Y3VSTW, Y2VSD)
 
        IF(EXPFLG.EQ.0) THEN
 
C         CLEAR THE FLAG USED TO LIMIT THE NUMBER OF ERROR MESSAGES
C         WRITTEN FOR INVALID HEADWATER RATIO TO ONE
 
          RATFLG = 0
 
C         COMPUTE TYPE 7 FREE FLOW HERE
 
C         OUTPUT SPECIAL WARNING MESSAGE.
 
          WRITE(STDOUT,49)
 
          CALL FRFT7
     I              (STDOUT, Z1TRUE, ID, ZSBRDF, 'FREE    ', 0.0,
     I               Q3VSRD, Y3VSTW,
     M               Z4,
     O               QFREE, FDROP)
          IF(EFLAG.NE.0) RETURN
          FRTYPE = 7
C         Do the energy check on the free flow.
C         Disable checking of energy losses involving section 2
          CHK23 = 0
          IF(ALPHA3.EQ.0.0) THEN
            DE34 = 0.0
          ELSE
            DE34 = ALPHA3
          ENDIF
          CALL CHKTEL
     I               (STDOUT, CHK23,
     M                DE34,
     O                IFLAG, CL34, DE14)
 
        ELSE
C         MAKE SURE THAT FREE FLOW IN THE CULVERT HAS NOT SUBMERGED
C         THE FLOW OVER THE ROADWAY FOR OTHER THAN TYPE = 0.
 
          IF(TYPE.NE.0) THEN
            IF(FREED.LT.FDRDW) THEN
              WRITE(STDOUT,89)
              WRITE(STDOUT,*) ' FREED=',FREED,' FDRDW=',FDRDW
              EFLAG = 1
              RETURN
            ENDIF
          ENDIF
        ENDIF
 
        IF(FRTYPE.EQ.7) THEN
          WRITE(STDOUT,41) HUP, FDROP
          WRITE(STDOUT,43)
          H = HDATUM
          IF(WFRD.EQ.0.0) THEN
            WRITE(CL34C,'(F5.2)') CL34
            WRITE(CHAR7,'(F7.3)') DE34
          ELSE
            CL34C = ' ----'
            CHAR7 = ' ------'
          ENDIF
C         Limit Z3 to the soffit at the outlet
          IF(Z3.GT.Z3SOF) THEN
            ZT = Z3SOF
          ELSE
             ZT = Z3
          ENDIF
          CALL VAR_DECIMAL(QFREE,
     O                       CQFREE)
          CALL VAR_DECIMAL(WFRD,
     O                       CWFRD)
          CALL VAR_DECIMAL(QFREE+WFRD,
     O                       CBOTH)

          WRITE(STDOUT,45) 1.0, ZT-ZB3, Z3P-H, Z44-H, Z4-H,
     A                 FDROP, CQFREE, CWFRD, CBOTH,
     B                  CL34C, CHAR7, DE14
        ELSE
          WRITE(STDOUT,41) HUP, FDROP
          WRITE(STDOUT,43)
          H = HDATUM
 
          IF(WFRD.EQ.0.0) THEN
            WRITE(CL34C,'(F5.2)') CL34
            WRITE(CHAR7,'(F7.3)') DE34
          ELSE
            CL34C = ' ----'
            CHAR7 = ' ------'
          ENDIF
C         Limit Z3 to the soffit at the outlet
          IF(Z3SAV.GT.Z3SOF) THEN
            ZT = Z3SOF
          ELSE
             ZT = Z3SAV
          ENDIF
          CALL VAR_DECIMAL(QFREE,
     O                       CQFREE)
          CALL VAR_DECIMAL(RDFLOW,
     O                       CWFRD)
          CALL VAR_DECIMAL(QFREE+RDFLOW,
     O                       CBOTH)
          WRITE(STDOUT,45) 1.0, ZT-ZB3, Z3PSAV-H, Z44SAV-H, Z4SAV-H,
     A                 FDROP, CQFREE, CWFRD, CBOTH,
     B                 CL34C, CHAR7, DE14
        ENDIF
C       Set the values for zero partial free drop.  Known to be zero
C       therefore not computed. 
        QMAT(I,1) = 0.0

C        Y2MAT(I,1) = Z1TRUE - ZB2
        Y2MAT(I,1) = Y2

C       Store the free flow value. Note there is a row for each 
C       upstream head in QMAT.  Each column then contains the 
C       flows at a constant free drop as the upstream head
C       varies.  Also only the non-zero upstream heads are 
C       present in QMAT.  The flow at zero upstream head is
C       always zero!
        QMAT(I,NFRAC) = QFREE + RDFLOW
        FDVEC(I) = FDROP
        QCLV = QFREE
        DDROP = 0.0
        Z4OLD = Z4
        QOLD = QFREE
 
C       FOR EACH OF THE PARTIAL FREE DROPS(EXCLUDING 0.0 AND 1.0)
C       COMPUTE THE FLOW OVER THE ROADWAY AND THE FLOW THROUGH THE
C       CULVERT.
 
        DQ = 0.0
C        WRITE(STDOUT,*) ' CULVRT: NFRAC=',NFRAC
        DO 500 J=NFRAC-1,2,-1
 
C         CLEAR THE FLAG USED TO LIMIT THE NUMBER OF ERROR MESSAGES
C         WRITTEN FOR INVALID HEADWATER RATIO TO ONE
 
          RATFLG = 0
 
          DROP = FDROP*PFDVEC(J)
          ED = EU - FDROP*PFDVEC(J)
 
C         Compute change in section 4 elevation.
          DZ4 = ED - Z4OLD
          Z4OLD = ED
 
C          WRITE(STDOUT,*) ' DZ4=',DZ4
 
C          WRITE(STDOUT,*) ' EU=',EU,' ED=',ED,' FDROP=',FDROP
 
C          WRITE(STDOUT,*) ' FRTYPE=',FRTYPE
 
 
C         Disable checking of energy losses involving section 2
          CHK23 = 0
          CALL FRFT7
     I              (STDOUT, EU, ID, ZSBRDF, 'SUBMERGE', DZ4, Q3VSRD,
     I               Y3VSTW,
     M               ED,
     O               QCLV, DDROP)
 
 
          DQ = QCLV - QOLD
          QOLD = QCLV
C          WRITE(STDOUT,*) ' DQ=',DQ
          IF(DQ.GT.0.0) DQ = 0.0
          IF(ALPHA3.EQ.0.0) THEN
            DE34 = 0.0
          ELSE
            DE34 = ALPHA3
          ENDIF
          CALL CHKTEL
     I               (STDOUT, CHK23,
     M                DE34,
     O                IFLAG, CL34, DE14)
 
          QMAT(I,J) = QCLV + WFRD
C         Find the depth at section 2.  Given in a 1-D table
C         as function of drop from section 1 to section 43. 
          CALL LKTAB
     I              (Y2VSD, Z1TRUE-Z43, 0,
     O                Y2, NTAB, DY2)

          Y2MAT(I,J) = Y2
          IF(Z3.GT.Z3SOF) THEN
            ZT = Z3SOF
          ELSE
            ZT = Z3
          ENDIF
          IF(WFRD.EQ.0.0) THEN
            WRITE(CL34C,'(F5.2)') CL34
            WRITE(CHAR7,'(F7.3)') DE34
          ELSE
            CL34C = ' ----'
            CHAR7 = ' ------'
          ENDIF
          CALL VAR_DECIMAL(QCLV,
     O                       CQCLV)
          CALL VAR_DECIMAL(WFRD,
     O                       CWFRD)
          CALL VAR_DECIMAL(QCLV+WFRD,
     O                       CBOTH)
          WRITE(STDOUT,45) PFDVEC(J), ZT-ZB3, Z3P-H, Z44-H,
     A                     Z4-H, DROP, CQCLV, CWFRD, CBOTH,
     B                     CL34C, CHAR7, DE14
C          ENDIF
 
 500    CONTINUE
        WRITE(STDOUT,44)
        WRITE(STDOUT,48) DDN
 2000 CONTINUE
 2001 CONTINUE
C     OUTPUT THE TABLE
      IF(EFLAG.EQ.0) THEN
        zrhufd = 0.0
        CALL TWDOUT
     I             (STDOUT, STDTAB, TAB, LABEL, NHU, NFRAC, HUVEC,
     I              FDVEC, PFDVEC, QMAT, HDATUM+ZDATUM,
     I              TABTYP, ' CULVERT', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)

        if(twod_cubic_out.eq.'YES') then
          verbose = 1
          CALL twodfit
     I             (STDOUT, TAB, NHU, NFRAC, HUVEC,
     I              FDVEC, PFDVEC, QMAT, HDATUM+ZDATUM,
     I              TABTYP, ' CULVERT', zrhufd, verbose,
     M              FTP,
     O              EFLAG, ftpup)
        endif

        IF(PUTQ.GT.0) THEN
          CALL TWOD13PUT
     I                  (STDOUT, PUTQ, NHU, NFRAC, HUVEC, FDVEC,
     I                   PFDVEC, QMAT, HDATUM+ZDATUM, TABTYP, 
     M                   FTP, 
     O                   EFLAG)
        ENDIF
        IF(PUTY2.GT.0) THEN

C          DO 5000 I=1,NHU
C            WRITE(STDOUT,'(I5,25F10.3)') I, (Y2MAT(I,J), J=1,NFRAC)
C5000      CONTINUE
          CALL TWOD13PUT
     I                  (STDOUT, PUTY2, NHU, NFRAC, HUVEC, FDVEC,
     I                   PFDVEC, Y2MAT, HDATUM+ZDATUM, TABTYP, 
     M                   FTP, 
     O                   EFLAG)
C          NHU = 0
C          NFRAC =0
C          HUVEC = 0.0
C          FDROP = 0.0
C          PFDVEC = 0.0
C          Y2MAT = 0.0
C          HDATUM = 0.0
C          TABTYP = 0
C          CALL TWOD13GET
C     I                  (STDOUT, 
C     O                   PUTY2, NHU, NFRAC, HUVEC, FDVEC,
C     O                   PFDVEC, Y2MAT, HDATUM, TABTYP, 
C     O                   EFLAG)
C        CALL TWDOUT
C     I             (STDOUT, STDOUT, PUTY2, LABEL, NHU, NFRAC, HUVEC,
C     I              FDVEC, PFDVEC, Y2MAT, HDATUM,
C     I              TABTYP,' CULVERT',
C     O              EFLAG)
        ENDIF

      ENDIF
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END

