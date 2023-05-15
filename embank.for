C
C
C
      REAL FUNCTION   FSBRAT
     I                      (ADR)
 
C     + + + PURPOSE + + +
C     Find the submergence ratio for the submergence table
C     stored at address ADR
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADR
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADR    - Address of submergence table
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER BADR, HADR, I, IOFF, TYPE
      CHARACTER TABID*16 
C     + + + INTRINSICS + + +
      INTRINSIC ABS

C     External names
      INTEGER LENSTR
      CHARACTER GET_TABID*16
      EXTERNAL LENSTR, GET_TABID
      
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:600* TYPE=',I4,' invalid type for submergence',
     A     ' table.')
 52   FORMAT(' *ERR:601* TABID=',A,' invalid submergence values at',
     A     ' table start.')
 54   FORMAT(' *ERR:602* TABID=',A,' no submergence in submergence',
     A     ' table.')
C***********************************************************************
      TYPE = ITAB(ADR+2)
      BADR = ADR + OFF234
      HADR = ITAB(ADR)
      IF(TYPE.LT.2.OR.TYPE.GT.4) THEN
C       INVALID TYPE FOR SUBMERGENCE
        WRITE(STD6,50) TYPE
        STOP 'Abnormal stop. Errors found.' 
      ENDIF
 
      IF(TYPE.EQ.2) THEN
          IOFF = 2
      ELSE
          IOFF = 3
      ENDIF
 
C     NOW SEARCH THE TABLE IN ORDER OF ADDRESSES(ASCENDING ARGUMENTS)
C     FOR THE HIGHEST ARGUMENT  WITH A FUNCTION VALUE OF 1.0
C     CHECK FOR CONSISTENCY
 
      IF(ABS(FTAB(BADR+1) - 1.0).GE.1.E-6.OR.FTAB(BADR).NE.0.0) THEN
C       FIRST ARGUMENT IN TABLE MUST BE 0.0 AND FIRST FUNCTION VALUE
C       MUST BE 1.0.
        TABID = GET_TABID(ITAB(ADR+1))
        WRITE(STD6,52) TABID(1:LENSTR(TABID))
        STOP 'Abnormal stop. Errors found.' 
      ENDIF
 
      DO 100 I=BADR,HADR,IOFF
        IF(1.0 - FTAB(I+1).GT.1.E-6) THEN
C         FOUND FIRST FUNCTION VALUE < 1.0.  BACKUP ARGUMENT.
 
          FSBRAT = FTAB(I-IOFF)
          RETURN
        ENDIF
 100  CONTINUE
 
C     FALL THROUGH INDICATES THAT THERE IS NO SUBMERGENCE!
      TABID = GET_TABID(ITAB(ADR+1))
      WRITE(STD6,54) TABID(1:LENSTR(TABID))
      STOP 'Abnormal stop. Errors found.' 
      END
C
C
C
      SUBROUTINE   FTOTHQ
     I                   (HCWTAB, LCWTAB, HLCRIT, HLMAX, SRAT, L, H,
     I                    DEPTH,
     M                    HLFLAG, HPFLAG,
     O                    Q, HTOT, YCREST, AELRAT)
 
C     + + + PURPOSE + + +
C     Compute the free flow using the total head and coefficients
C     defined as a function of head and weir width.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER HCWTAB, HLFLAG, HPFLAG, LCWTAB
      REAL AELRAT, DEPTH, H, HLCRIT, HLMAX, HTOT, L, Q, SRAT, YCREST
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     HCWTAB - High head weir coefficient table
C     LCWTAB - Low head weir coefficient table
C     HLCRIT - Ratio of piezometric head to crest breadth at boundary
C              between low head and high head flow
C     HLMAX  - Maximum value of piezometric head to crest breadth ratio
C               above which a warning message is issued
C     SRAT   - Submergence ratio at the free flow limit
C     L      - Crest breadth for the embankment
C     H      - Piezometric head
C     DEPTH  - Depth of approaching flow
C     HLFLAG - Warning message suppression flag
C     HPFLAG - Warning message suppression flag for invalid weir flow
C     Q      - Flowrate
C     HTOT   - Total head
C     YCREST - Estimated depth on crest of the embankment
C     AELRAT - Estimated ratio of head loss to head loss at incipient
C              submergence
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT, NTAB
      REAL AEL, CW, EL, HEAD, PDV, QMAX, QW, QWOLD, RATIO
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTAB
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *WRN:545* Head to width ratio=',F8.3,' >',
     A  ' maximum ratio=',F8.3)
 52   FORMAT(/,' *WRN:553* (Piezometric head)/(weir height)=',F8.2,
     a       ' > 4.','  Weir',/,10X,'flow may be INVALID.')
C***********************************************************************
      KNT = 0
C     USE LINEAR ITERATION TO INCLUDE THE EFFECT OF VELOCITY HEAD
 
      QWOLD = 0.475812*SQRT_GRAV*H*SQRT(H)
      IF(L.EQ.0.0) THEN
        WRITE(STD6,*) ' *BUG:XXX* L=0 IN FTOTHQ'
        STOP 'Abnormal stop. Errors found.' 
      ENDIF
      RATIO = H/L
 
 100  CONTINUE
 
C       COMPUTE THE VELOCITY HEAD AND ADD TO THE PIEZOMETRIC HEAD
 
        HEAD = H + (QWOLD/DEPTH)**2/GRAV2
 
C       FIND THE WEIR COEFFICIENT
 
 
        IF(RATIO.GT.HLCRIT) THEN
          IF(RATIO.GT.HLMAX) THEN
            IF(HLFLAG.EQ.0) THEN
              WRITE(STD6,50) RATIO, HLMAX
              HLFLAG = 1
            ENDIF
          ENDIF
          CALL LKTAB
     I              (HCWTAB, RATIO, 0,
     O               CW, NTAB, PDV)
        ELSE
          CALL LKTAB
     I              (LCWTAB, HEAD, 0,
     O               CW, NTAB, PDV)
        ENDIF
 
        QW = CW*HEAD*SQRT(HEAD)
 
        IF(ABS(QW - QWOLD)/(0.5*(QWOLD + QW)).LE.5.E-4) THEN
          Q = QW
          HTOT = HEAD
 
C         MAKE AN ESTIMATE OF THE CREST HEIGHT TO ESTIMATE
C         THE MOMENTUM FLUX
          YCREST = (CW**2*HTOT**3/GRAV)**0.3333333
C         COMPUTE THE ENERGY LOSS IMPLIED BY THE FLOW OVER
C         THE WEIR IN THOSE CASES WHEN CRITICAL DEPTH OCCURS AT OR
C         NEAR THE CREST.
 
          AEL = HEAD - 1.5*YCREST
          IF(AEL.LT.0.0) AEL = 0
C         COMPUTE THE HEAD LOSS AT INCIPIENT SUBMERGENCE
 
          EL = HEAD*(1.0 - SRAT) -
     A            (QW/(DEPTH - H + SRAT*HEAD))**2/GRAV2
          IF(EL.LT.0.0) THEN
C           FORCE EFFECTIVE ZERO LOSS FOR COMPUTING CREST
C           ELEVATION
            EL = 1.0
            AEL = 0.0
          ENDIF
C         NOW COMPUTE THE RATIO OF THE APPROACH HEAD LOSS TO THE
C         THE HEAD LOSS AT INCIPIENT SUBMERGENCE.
 
          AELRAT = AEL/EL
          IF(AELRAT.LE.0.0) THEN
            AELRAT = 0.005
          ENDIF
          IF(1.25*H.GT.DEPTH) THEN
C           ISSUE WARNING: H/(WEIR HEIGHT) > 4.  WEIR FLOW MAY BE
C           INVALID.
            IF(HPFLAG.EQ.0) THEN
              WRITE(STD6,52)  H/(DEPTH - H)
              HPFLAG = 1
            ENDIF
          ENDIF
C         CHECK FLOW AGAINST CRITICAL FLOW IN APPROACH SECTION
          QMAX = DEPTH*SQRT(GRAV*DEPTH)
          IF(Q.GT.QMAX) Q = QMAX
          RETURN
        ELSE
          QWOLD = QW
          KNT = KNT + 1
          IF(KNT.GT.100) THEN
            WRITE(STD6,*) ' KNT > 100 IN FTOTHQ.'
            WRITE(STD6,*) ' FTOTHQ: HEAD=',HEAD,' QW=',QW,
     A             ' QWOLD=',QWOLD,' DEPTH=',DEPTH,' CW=',CW,' H=',H
            STOP 'Abnormal stop. Errors found.' 
          ENDIF
 
          GOTO 100
        ENDIF
 
      END
C
C
C
      SUBROUTINE   STOTHQ
     I                   (HCWTAB, LCWTAB, SUBTAB, HLCRIT, L, H, HTAIL,
     I                    DEPTH,
     M                    HTOT,
     O                    Q)
 
C     + + + PURPOSE + + +
C     Compute the submerged flow using the total head and coefficients
C     defined as a function of head and weir width. On entry Q and HTOT
C     should have their free flow values.  Submerged values are
C     returned.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER HCWTAB, LCWTAB, SUBTAB
      REAL DEPTH, H, HLCRIT, HTAIL, HTOT, L, Q
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     HCWTAB - High head weir coefficient table
C     LCWTAB - Low head weir coefficient table
C     SUBTAB - Address for submergence table
C     HLCRIT - Ratio of piezometric head to crest breadth at boundary
C              between low head and high head flow
C     L      - Crest breadth for the embankment
C     H      - Piezometric head
C     HTAIL  - Tailwater head
C     DEPTH  - Depth of approaching flow
C     HTOT   - Total head
C     Q      - Flowrate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT, NTAB
      REAL CW, FRAC, HEAD, PDV, QMAX, QW, QWOLD, RAT, RATIO
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL LKTAB
C***********************************************************************
      KNT = 0
C      WRITE(STD6,*) ' STOTHQ:',HCWTAB, LCWTAB, SUBTAB, HLCRIT, L, H,
C     A            HTAIL, DEPTH
 
C     USE LINEAR ITERATION TO INCLUDE THE EFFECT OF VELOCITY HEAD
 
      RAT = HTAIL/HTOT
      CALL LKTAB
     I          (SUBTAB, RAT, 0,
     O           FRAC, NTAB, PDV)
 
C      WRITE(STD6,*) ' STOTHQ: RAT=',RAT,' FRAC=',FRAC,' HTAIL=',HTAIL,
C     A           ' HTOT=',HTOT
 
      IF(FRAC.EQ.0.0) THEN
C        HTOT = H
        Q = 0.0
        RETURN
      ELSEIF(FRAC.GT.0.9999.AND.H.GT.HTAIL) THEN
C       RETAIN FREE FLOW VALUES
        RETURN
      ELSE
C       USE FRAC ON FREE FLOW VALUE AS AN ESTIMATE OF THE STARTING
C       FLOW FOR LINEAR ITERATION
        QWOLD = FRAC*Q
      ENDIF
      IF(L.EQ.0.0) THEN
        WRITE(STD6,*) ' L=0 IN STOTHQ'
        STOP 'Abnormal stop. Errors found.' 
      ENDIF
      RATIO = H/L
 
C     FIND THE WEIR COEFFICIENT
 
      IF(RATIO.GT.HLCRIT) THEN
        CALL LKTAB
     I            (HCWTAB, RATIO, 0,
     O             CW, NTAB, PDV)
      ELSE
        CALL LKTAB
     I            (LCWTAB, HTOT, 0,
     O             CW, NTAB, PDV)
      ENDIF

 100  CONTINUE
 
C       COMPUTE THE VELOCITY HEAD AND ADD TO THE PIEZOMETRIC HEAD
 
        HEAD = H + (QWOLD/DEPTH)**2/GRAV2
 
 
C       FIND THE SUBMERGENCE EFFECT
 
        RAT = HTAIL/HEAD
        CALL LKTAB
     I            (SUBTAB, RAT, 0,
     O             FRAC, NTAB, PDV)
 
        QW = CW*HEAD*SQRT(HEAD)*FRAC
 
C        WRITE(STD6,*) ' '
C        WRITE(STD6,*) ' KNT=',KNT
C        WRITE(STD6,*) ' STOTHQ: TOTH=',HEAD,' FRAC=',FRAC,' QW=',QW,
C     A             ' QWOLD=',QWOLD,' DEPTH=',DEPTH,' CW=',CW,
C     B             ' HTAIL=',HTAIL
 
        IF(ABS(QW - QWOLD)/(QWOLD + QW + 1.E-3).LE.5.E-4) THEN
          Q = QW
          HTOT = HEAD
        ELSE
          KNT = KNT + 1
          IF(KNT.GT.100) THEN
            WRITE(STD6,*) ' KNT > 100 IN STOTHQ.'
            WRITE(STD6,*) ' STOTHQ: HEAD=',HEAD,' FRAC=',FRAC,' QW=',QW,
     A' QWOLD=',QWOLD,' DEPTH=',DEPTH,' CW=',CW,' H=',H,' HTAIL=',HTAIL
            STOP 'Abnormal stop. Errors found.' 
          ENDIF
          QWOLD = QW
 
          GOTO 100
        ENDIF
 
C     CHECK FLOW AGAINST CRITICAL FLOW IN APPROACH SECTION
      QMAX = DEPTH*SQRT(GRAV*DEPTH)
      IF(Q.GT.QMAX) Q = QMAX
 
C      WRITE(STD6,*) ' EXITING STOTHQ. Q=',Q
      RETURN
 
      END
C
C
C
      SUBROUTINE   EMBSUB
     I                   (MINLOC, MINCRS, NOFF, SURF, TOTHL, TOTHR, HU,
     O                    HD, FREED, ED)
 
C     + + + PURPOSE + + +
C     Find the submergence ratio to use for flow over an embankment
C     and the elevation of the tailwater at which submergence
C     begins. HU is the piezometric head on the embankment,
C     ED is tailwater elevation at which submergence begins.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER MINLOC, NOFF
      REAL ED, HD, HU, MINCRS, TOTHL(*), TOTHR(*)
      REAL*8 FREED
      CHARACTER SURF(*)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     MINLOC - Offset of the minimum crest elevation
C     MINCRS - Minimum crest elevation
C     NOFF   - Number of offsets
C     SURF   - Nature of the embankment surface
C     TOTHL  - Total head at left hand end of segment
C     TOTHR  - Total head at right hand end of segment
C     HU     - Head upstream
C     HD     - Tailwater head at free flow limit
C     FREED  - Free drop
C     ED     - Downstream water surface elevation
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'subcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ISEG
      REAL HTOT, RATIO, RATIOL, RATIOR
 
C     + + + INTRINSICS + + +
      INTRINSIC MAX, MIN
C***********************************************************************
C     FIND THE RATIO TO USE FOR DETERMNING THE FREE DROP VALUE
 
      IF(MINLOC.EQ.1) THEN
        IF(SURF(MINLOC).EQ.'PAVED') THEN
          RATIO = PRAT
        ELSE
          RATIO = GRAT
        ENDIF
      ELSE IF(MINLOC.EQ.NOFF) THEN
        IF(SURF(MINLOC-1).EQ.'PAVED') THEN
          RATIO = PRAT
        ELSE
          RATIO = GRAT
        ENDIF
      ELSE
C       CHECK BOTH SEGMENTS
        IF(SURF(MINLOC).EQ.'PAVED') THEN
          RATIOR = PRAT
        ELSE
          RATIOR = GRAT
        ENDIF
        IF(SURF(MINLOC-1).EQ.'PAVED') THEN
          RATIOL = PRAT
        ELSE
          RATIOL = GRAT
        ENDIF
        RATIO = MIN(RATIOL, RATIOR)
      ENDIF
 
C     COMPUTE THE DOWNSTREAM PIEZOMETRIC HEAD AT WHICH SUBMERGENCE
C     EFFECT BEGINS FOR AT LEAST ONE POINT ON THE WEIR.
 
C     FIND THE UPSTREAM TOTAL HEAD AT THE MINIMUM POINT
 
      ISEG = MINLOC
      IF(MINLOC.EQ.NOFF) ISEG = MINLOC - 1
      HTOT = MAX(TOTHL(ISEG), TOTHR(ISEG))
 
C      WRITE(STD6,*) ' HTOT AT MINCRS=',HTOT,' MINLOC=',MINLOC
C      WRITE(STD6,*) ' RATIO AT MINCRS=',RATIO
 
      HD = RATIO*HTOT
 
      FREED = DBLE(HU) - DBLE(RATIO)*DBLE(HTOT)
 
      ED = HD + MINCRS
 
      RETURN
      END
C
C
C
      SUBROUTINE   SBFEMB
     I                   (NFRAC, EU, HU, MINCRS, MINLOC, PLCWTB, GLCWTB,
     I                    PHCWTB, GHCWTB, PSUBTB, GSUBTB, NOFF, SURF,
     I                    HLCRIT, XL, XR, HL, HM, HR, QL, QM, QR, TOTHL,
     I                    TOTHM, TOTHR, APPL, APPM, APPR, WL, WM, WR,
     I                    QFREE, PFDVEC,
     O                    FREED, QSBVEC)
 
C     + + + PURPOSE + + +
C     SuBmerged Flow EMBankment- Find the submerged flows for
C        the given proportions of free drop and the free flow
C        parameters computed by FRFEMB.
 
      IMPLICIT NONE

      INCLUDE 'stdun.cmn'

C     + + + DUMMY ARGUMENTS + + +
      INTEGER GHCWTB, GLCWTB, GSUBTB, MINLOC, NFRAC, NOFF, PHCWTB,
     A        PLCWTB, PSUBTB
      REAL APPL(*), APPM(*), APPR(*), EU, FREED, HL(*), HLCRIT, HM(*),
     A     HR(*), HU, MINCRS, PFDVEC(*), QFREE, QL(*), QM(*), QR(*),
     B     QSBVEC(*), TOTHL(*), TOTHM(*), TOTHR(*), WL(*), WM(*), WR(*),
     C     XL(*), XR(*)
      CHARACTER SURF(*)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     NFRAC  - Number of fractions for defining partial free drop
C     EU     - Water surface elevation upstream.
C     HU     - Head upstream
C     MINCRS - Minimum crest elevation
C     MINLOC - Offset of the minimum crest elevation
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
C     APPL   - Elevation of approach at left end of line segment
C     APPM   - Elevation of approach at middle of line segment
C     APPR   - Elevation of approach at right end of line segment
C     WL     - Breadth of the crest at left hand end of segment
C     WM     - Breadth of the crest at middle of segment
C     WR     - Breadth of the crest at right hand end of segment
C     QFREE  - Free flow
C     PFDVEC - Partial free drop vector
C     FREED  - Free drop
C     QSBVEC - Submerged flow for each partial free drop
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HCWTAB, IP, ISEG, LCWTAB, SUBTAB
      REAL DEPTH, DX, ED, FRAC, HD, QLEFT, QMID, QRIGHT, QSEG, SHTOT
      DOUBLE PRECISION QSUB, FDROP
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL EMBSUB, STOTHQ
C     *******************************FORMATS****************************
50    FORMAT(/,' *ERR:740* Drop to free flow=',F10.4,
     A        ' < 0.0 in EMBANKQ.')
C***********************************************************************
C         COMPUTE THE SUBMERGENCE EFFECT OVER THE RANGE OF TAIL WATER
C         ELEVATIONS WHICH CAN AFFECT THE FLOW OVER THE EMBANKMENT.
 
C         FOR THE GIVEN UPSTREAM PIEZOMETRIC HEAD, FIND THE
C         DROP BETWEEN THE UPSTREAM AND DOWNSTREAM PIEZOMETRIC HEADS
C         REQUIRED TO ENABLE FREE FLOW.
C         THE DOWNSTREAM EFFECT IS FIRST FELT AT THE POINT OF MAXIMUM
C         HEAD.  RATIO OF TAILWATER PIEZOMETRIC HEAD TO THE UPSTREAM
C         TOTAL HEAD IS THE ARGUMENT USED BY THE USGS.
 
 
          CALL EMBSUB
     I               (MINLOC, MINCRS, NOFF, SURF, TOTHL, TOTHR, HU,
     O                HD, FDROP, ED)
          FREED = FDROP
 
          IF(FREED.LT.0.0) THEN
            WRITE(STD6,50) FREED
            STOP 'Abnormal stop. Error found.'
          ENDIF
C          WRITE(STD6,*) ' '
C          WRITE(STD6,*) ' AT FREE DROP. EU=',EU,' ED=',ED
C          WRITE(STD6,'('' FREE DROP='',F10.3)') FREED
 
C         NOW FOR A SERIES OF PROPORTIONS OF THE FREE DROP COMPUTE
C         THE FLOW OVER THE EMBANKMENT. AT ZERO PROPORTION OF FREE DROP
C         THE FLOW IS ZERO AND AT FREE DROP THE FLOW IS THE FREE FLOW.
C         THE ZERO FLOW AND THE FREE FLOW CONDITION ARE KNOWN SO THEY
C         ARE NOT RECOMPUTED.
 
          QSBVEC(1) = 0.0
          QSBVEC(NFRAC) = QFREE
          DO 4500 IP=NFRAC-1,2,-1
            QSUB = 0.D0
            FRAC = PFDVEC(IP)
 
C            WRITE(STDOUT,*) ' '
C            WRITE(STDOUT,*) ' SUBMERGED FLOWS'
C            WRITE(STDOUT,*) ' FRAC=',FRAC
 
C           COMPUTE THE DOWNSTREAM WATER SURFACE ELEVATION FOR THIS
C           FRACTION OF THE FREE DROP
 
            ED = EU - DBLE(FRAC)*FDROP
 
C            WRITE(STDOUT,'('' EU='',F10.4,'' ED='',F10.4)') EU, ED
 
C           FOR EACH SEGMENT OF THE WEIR FIND THE EFFECT OF THIS
C           VALUE OF DOWNSTREAM PIEZOMETRIC HEAD
 
            DO 4100 ISEG=1,NOFF-1
 
C              WRITE(STDOUT,*) ' '
C              WRITE(STDOUT,*) ' ISEG=',ISEG
 
C             SELECT THE TABLE FOR SUBMERGENCE
              IF(SURF(ISEG).EQ.'PAVED') THEN
                SUBTAB = PSUBTB
              ELSE
                SUBTAB = GSUBTB
              ENDIF
C             SELECT THE WEIR COEF. TABLES BASED ON THE SURFACE CONDITION.
 
              IF(SURF(ISEG).EQ.'PAVED') THEN
                HCWTAB = PHCWTB
                LCWTAB = PLCWTB
              ELSE
                HCWTAB = GHCWTB
                LCWTAB = GLCWTB
              ENDIF
              IF(HL(ISEG).GT.0.0.OR.HR(ISEG).GT.0.0) THEN
C               SEGMENT HAS NON-ZERO FREE FLOW
 
                QLEFT = QL(ISEG)
                SHTOT = TOTHL(ISEG)
                IF(QLEFT.GT.0.0) THEN
C                 NOTE THAT THE UPSTREAM WATER SURFACE ELEVATION
C                 LESS THE PIEZOMETRIC HEAD GIVES THE WEIR CREST
C                 ELEVATION
 
                  HD = ED - (EU - HL(ISEG))
                  IF(HD.LT.0.0) HD = 0.0
 
                  DEPTH = EU - APPL(ISEG)
                  CALL STOTHQ
     I                       (HCWTAB, LCWTAB, SUBTAB, HLCRIT, WL(ISEG),
     I                        HL(ISEG), HD, DEPTH,
     M                        SHTOT,
     O                        QLEFT)
                ELSE
                  QLEFT = 0.0
                ENDIF
 
C               WRITE(STDOUT,*) ' LEFT: DEPTH=',DEPTH,' HD=',HD,' QLEFT=',
C    A                     QLEFT,' SHTOT=',SHTOT, ' HEAD=',HL(ISEG)
 
                QMID = QM(ISEG)
                SHTOT = TOTHM(ISEG)
                IF(QMID.GT.0.0) THEN
                  HD = ED - (EU - HM(ISEG))
                  IF(HD.LT.0.0) HD = 0.0
                  DEPTH = EU - APPM(ISEG)
                  CALL STOTHQ
     I                       (HCWTAB, LCWTAB, SUBTAB, HLCRIT, WM(ISEG),
     I                        HM(ISEG), HD, DEPTH,
     M                        SHTOT,
     O                        QMID)
 
                ELSE
                  QMID = 0.0
                ENDIF
 
C               WRITE(STDOUT,*) ' MID: DEPTH=',DEPTH,' HD=',HD,' QMID=',
C    A                     QMID,' SHTOT=',SHTOT, ' HEAD=',HM(ISEG)
 
                QRIGHT = QR(ISEG)
                SHTOT = TOTHR(ISEG)
                IF(QRIGHT.GT.0.0) THEN
                  HD = ED - (EU - HR(ISEG))
                  IF(HD.LT.0.0) HD = 0.0
                  DEPTH = EU - APPR(ISEG)
                  CALL STOTHQ
     I                       (HCWTAB, LCWTAB, SUBTAB, HLCRIT, WR(ISEG),
     I                        HR(ISEG), HD, DEPTH,
     M                        SHTOT,
     O                        QRIGHT)
                ELSE
                  QRIGHT = 0.0
                ENDIF
 
C               WRITE(STDOUT,*) ' RIGHT: DEPTH=',DEPTH,' HD=',HD,' QRIGHT=',
C    A                     QRIGHT,' SHTOT=',SHTOT, ' HEAD=',HR(ISEG)
 
C               NOW COMPUTE THE FLOW AS AFFECTED BY SUBMERGENCE
 
                DX = ABS(XR(ISEG) - XL(ISEG))
 
C               INTEGRATE OVER THE WETTED LENGTH USING SIMPSON'S RULE
 
                QSEG = DX*(QLEFT + 4.*QMID + QRIGHT)/6.
C          WRITE(STDOUT,*) ' '
C          WRITE(STDOUT,*) ' DX =',DX
C          WRITE(STDOUT,'('' SEG FLOWS:'',3F12.4)') QLEFT, QMID,
C     A                                  QRIGHT
 
C                WRITE(STDOUT,'('' QSEG='',F12.4)')
C     A                          QSEG
 
                QSUB = QSUB + QSEG
              ENDIF
 4100       CONTINUE
 
C           STORE THE VALUE AT THIS FRACTION OF THE FREE DROP
 
            QSBVEC(IP) = QSUB
 4500     CONTINUE
 
 
      RETURN
      END
C
C
C
      SUBROUTINE   FRFEMB
     I                   (HU, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB,
     I                    NOFF, OFF, CREST, WIDTH, APPROC, SURF, RMFFAC,
     I                    HLCRIT, HLMAX,
     M                    HLFLAG, HPFLAG,
     O                    EU, XL, XR, HL, HM, HR, QL, QM, QR, TOTHL,
     O                    TOTHM, TOTHR, YFL, YFM, YFR, APPL, APPM, APPR,
     O                    WL, WM, WR, AELL, AELM, AELR, QFREE, MFREE,
     O                    EFREE)
 
C     + + + PURPOSE + + +
C     FRee Flow EMBankment- Find free flow over an embankment
C     for a given upstream water surface elevation.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER GHCWTB, GLCWTB, HLFLAG, HPFLAG, NOFF, PHCWTB, PLCWTB
      REAL AELL(*), AELM(*), AELR(*), APPL(*), APPM(*), APPR(*),
     A     APPROC(*), CREST(*), EFREE, EU, HL(*), HLCRIT, HLMAX, HM(*),
     B     HR(*), HU, MFREE, MINCRS, OFF(*), QFREE, QL(*), QM(*), QR(*),
     C     RMFFAC, TOTHL(*), TOTHM(*), TOTHR(*), WIDTH(*), WL(*), WM(*),
     D     WR(*), XL(*), XR(*), YFL(*), YFM(*), YFR(*)
      CHARACTER SURF(*)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     HU     - Head upstream
C     MINCRS - Minimum crest elevation
C     PLCWTB - Paved low-head weir coefficient table
C     GLCWTB - Gravel surface low head weir coefficient table
C     PHCWTB - Paved high-head weir coefficient table
C     GHCWTB - Gravel surface high head weir coefficient table
C     NOFF   - Number of offsets
C     OFF    - Offsets for the crest of the embankment
C     CREST  - Crest profile for the high point on the embankment
C     WIDTH  - Breadth of the embankment crest
C     APPROC - Elevation of approach for the embankment
C     SURF   - Nature of the embankment surface
C     RMFFAC - Adjustment factor for roadway momentum flux
C     HLCRIT - Ratio of piezometric head to crest breadth at boundary
C              between low head and high head flow
C     HLMAX  - Maximum value of piezometric head to crest breadth ratio
C               above which a warning message is issued
C     HLFLAG - Warning message suppression flag
C     HPFLAG - Warning message suppression flag for invalid weir flow
C     EU     - Water surface elevation upstream.
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
C     QFREE  - Free flow
C     MFREE  - Momentum flux over roadway for free flow
C     EFREE  - Estimated energy flux for free flow over the embankment
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'subcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HCWTAB, IL, IR, ISEG, LCWTAB
      REAL DEPTH, DX, ESEG, MSEG, QSEG, RATIO, SRAT
      DOUBLE PRECISION EF, MF, QF
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FTOTHQ
C***********************************************************************
C       COMPUTE THE ELEVATION OF THE UPSTREAM WATER SURFACE
 
        EU = HU + MINCRS
 
        QF = 0.D0
        MF = 0.D0
        EF = 0.D0
        DO 4000 ISEG=1,NOFF-1
          IL = ISEG
          IR = ISEG + 1
 
C         Compute free flow for each segment to get the total.
C         retain the intermediate values to use for subsequent
C         submergence computations.
 
C         For each straight line segment of the weir crest find the
C         left and right limits for the wetted length. Four cases:
C         entire crest is above water, entire crest is below water,
C         left end is below water and right end is above water,
C         and right end is below water and left end is above water.
C         Water level is taken at piezometric head level and not at
C         actual wetted length because actual wetted length is unknown
C         and limited tests show that using the piezometric head
C         yields values close to experiment for triangular weirs.
 
          IF(EU.GT.CREST(IL).AND.EU.GT.CREST(IR)) THEN
C           ENTIRE CREST BELOW WATER
            XL(ISEG) = OFF(IL)
            XR(ISEG) = OFF(IR)
            HL(ISEG) = EU - CREST(IL)
            HR(ISEG) = EU - CREST(IR)
            APPL(ISEG) = APPROC(IL)
            APPR(ISEG) = APPROC(IR)
            WL(ISEG) = WIDTH(IL)
            WR(ISEG) = WIDTH(IR)
          ELSE IF(EU.GT.CREST(IL).AND.EU.LE.CREST(IR)) THEN
C           LEFT END WET AND RIGHT END DRY
            XL(ISEG) = OFF(IL)
            WL(ISEG) = WIDTH(IL)
            APPL(ISEG) = APPROC(IL)
            HL(ISEG) = EU - CREST(IL)
 
            XR(ISEG) = OFF(IL) + (OFF(IR) - OFF(IL))*
     A                      (EU - CREST(IL))/(CREST(IR) - CREST(IL))
 
            RATIO = (XR(ISEG) - OFF(IL))/(OFF(IR) - OFF(IL))
            WR(ISEG) = WIDTH(IL) + RATIO*(WIDTH(IR) - WIDTH(IL))
            APPR(ISEG) = APPROC(IL) + RATIO*(APPROC(IR) - APPROC(IL))
            HR(ISEG) = 0.0
          ELSE IF(EU.LE.CREST(IL).AND.EU.GT.CREST(IR)) THEN
C           LEFT END DRY AND RIGHT END WET
            XL(ISEG) = OFF(IR) + (EU - CREST(IR))*(OFF(IL) - OFF(IR))/
     A                                      (CREST(IL) - CREST(IR))
 
            RATIO = (XL(ISEG) - OFF(IR))/(OFF(IL) - OFF(IR))
            WL(ISEG) = WIDTH(IR) + RATIO*(WIDTH(IL) - WIDTH(IR))
            APPL(ISEG) = APPROC(IR) + RATIO*(APPROC(IL) - APPROC(IR))
            HL(ISEG) = 0.0
 
            XR(ISEG) = OFF(IR)
            WR(ISEG) = WIDTH(IR)
            APPR(ISEG) = APPROC(IR)
            HR(ISEG) = EU - CREST(IR)
          ELSE
C           LEFT END DRY AND RIGHT END DRY
 
            HL(ISEG) = -1.0
            HR(ISEG) = -1.0
            GOTO 4000
          ENDIF
 
 
C         COMPUTE THE MIDPOINT VALUES
 
          HM(ISEG) = 0.5*(HR(ISEG) + HL(ISEG))
          WM(ISEG) = 0.5*(WR(ISEG) + WL(ISEG))
          APPM(ISEG) = 0.5*(APPL(ISEG) + APPR(ISEG))
 
 
C          WRITE(STDOUT,*)  ' CHECK OF THE WEIR SEGMENT VALUES'
C          WRITE(STDOUT,*) ' ISEG=',ISEG
C          WRITE(STDOUT,'(''   X:'',2F10.2)') XL(ISEG), XR(ISEG)
C          WRITE(STDOUT,'(''   W:'',2F10.2)') WL(ISEG), WR(ISEG)
C          WRITE(STDOUT,'(''   H:'',2F10.4)') HL(ISEG), HR(ISEG)
C          WRITE(STDOUT,'('' APP:'',2F10.2)') APPL(ISEG), APPR(ISEG)
 
 
C         FOR EACH OF THE THREE LOCATIONS ALONG THIS SEGMENT COMPUTE
C         THE FLOW PER UNIT LENTH OF THE WEIR.  THE FLOW AT ZERO
C         PIEZOMETRIC HEAD IS TAKEN AS ZERO.
 
C         SELECT THE TABLES BASED ON THE SURFACE CONDITION.
 
          IF(SURF(ISEG).EQ.'PAVED') THEN
            HCWTAB = PHCWTB
            LCWTAB = PLCWTB
            SRAT = PRAT
          ELSE
            HCWTAB = GHCWTB
            LCWTAB = GLCWTB
            SRAT = GRAT
          ENDIF
 
 
          IF(HL(ISEG).LE.0.0) THEN
            QL(ISEG) = 0.0
            TOTHL(ISEG) = 0.0
            YFL(ISEG) = 1.0
          ELSE
            DEPTH = EU - APPL(ISEG)
            CALL FTOTHQ
     I                 (HCWTAB, LCWTAB, HLCRIT, HLMAX, SRAT, WL(ISEG),
     I                  HL(ISEG), DEPTH,
     M                  HLFLAG, HPFLAG,
     O                  QL(ISEG), TOTHL(ISEG), YFL(ISEG), AELL(ISEG))
          ENDIF
 
 
          IF(HM(ISEG).LE.0.0) THEN
            QM(ISEG) = 0.0
            TOTHM(ISEG) = 0.0
            YFM(ISEG) = 1.0
          ELSE
            DEPTH = EU - APPM(ISEG)
            CALL FTOTHQ
     I                 (HCWTAB, LCWTAB, HLCRIT, HLMAX, SRAT, WM(ISEG),
     I                  HM(ISEG), DEPTH,
     M                  HLFLAG, HPFLAG,
     O                  QM(ISEG), TOTHM(ISEG), YFM(ISEG), AELM(ISEG))
          ENDIF
 
 
          IF(HR(ISEG).LE.0.0) THEN
            QR(ISEG) = 0.0
            TOTHR(ISEG) = 0.0
            YFR(ISEG) = 1.0
          ELSE
            DEPTH = EU  - APPR(ISEG)
            CALL FTOTHQ
     I                 (HCWTAB, LCWTAB, HLCRIT, HLMAX, SRAT, WR(ISEG),
     I                  HR(ISEG), DEPTH,
     M                  HLFLAG, HPFLAG,
     O                  QR(ISEG), TOTHR(ISEG), YFR(ISEG), AELR(ISEG))
          ENDIF
 
C          WRITE(STDOUT,'('' TOTH:'',1P3E10.3)') TOTHL(ISEG),
C    A         TOTHM(ISEG), TOTHR(ISEG)
C          WRITE(STDOUT,'('' YF:'',1P3E10.3)') YFL(ISEG), YFM(ISEG),
C    A                   YFR(ISEG)
 
 
C         NOW COMPUTE THE FREE FLOW, FREE FLOW
C         MOMENTUM FLUX, AND FREE FLOW ENERGY FLUX FOR THIS SEGMENT
 
          DX = ABS(XR(ISEG) - XL(ISEG))
 
C         INTEGRATE OVER THE WETTED LENGTH USING SIMPSON'S RULE
 
          QSEG  =  DX*(QL(ISEG) + 4.*QM(ISEG) + QR(ISEG))/6.
          MSEG =   DX*(QL(ISEG)**2/YFL(ISEG) + 4.*QM(ISEG)**2/
     A                  YFM(ISEG) + QR(ISEG)**2/YFR(ISEG))/6.
 
          ESEG =   DX*(QL(ISEG)*((QL(ISEG)/YFL(ISEG))**2/GRAV2 +
     A                 (YFL(ISEG) + EU - HL(ISEG)))  +
     B                 4.*QM(ISEG)*((QM(ISEG)/YFM(ISEG))**2/GRAV2 +
     C                    (YFM(ISEG) + EU - HM(ISEG))) +
     D                 QR(ISEG)*((QR(ISEG)/YFR(ISEG))**2/GRAV2 +
     E                 (YFR(ISEG) + EU - HR(ISEG))) )/6.
 
 
C          WRITE(STDOUT,*) ' '
C          WRITE(STDOUT,*) ' DX =',DX
C          WRITE(STDOUT,'('' SEG FLOWS:'',3F12.4)') QL(ISEG), QM(ISEG),
C     A                                  QR(ISEG)
C          WRITE(STDOUT,'('' ISEG='',I5,'' QSEG='',F12.4)') ISEG, QSEG
C          WRITE(STDOUT,'('' MSEG='',F12.4,'' ESEG='',F12.4)')
C     A                             MSEG, ESEG
 
C          WRITE(STDOUT,*) ' AELL=',AELL(ISEG),' AELM=',AELM(ISEG),
C     A                     ' AELR=',AELR(ISEG)
 
C          WRITE(STDOUT,*) ' YFL=',YFL(ISEG),' YFM=',YFM(ISEG),
C     A                    ' YFR=',YFR(ISEG)
 
          QF = QF + QSEG
          MF = MF + MSEG
          EF = EF + ESEG
 4000   CONTINUE
 
 
      QFREE = QF
      MFREE = RMFFAC*MF
      EFREE = EF
C      WRITE(STDOUT,*) ' QFREE=',QFREE,' MFREE=',MFREE,' EFREE=',EFREE
      RETURN
      END
C
C
C
      SUBROUTINE   INPRDP
     I                   (STDIN, STDOUT, ZDATUM,
     O                    EFLAG, NOFF, MINCRS, MINLOC, OFF, CREST,
     O                    WIDTH, APPROC, SURF)
 
C     + + + PURPOSE + + +
C     INPut RoaDway Profile for computing flow over an embankment,
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, MINLOC, NOFF, STDIN, STDOUT
      REAL APPROC(PMXOFF), CREST(PMXOFF), MINCRS, OFF(PMXOFF),
     A     WIDTH(PMXOFF), ZDATUM
      CHARACTER SURF(PMXOFF)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     ZDATUM - datum for local elevation
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NOFF   - Number of offsets
C     MINCRS - Minimum crest elevation
C     MINLOC - Offset of the minimum crest elevation
C     OFF    - Offsets for the crest of the embankment
C     CREST  - Crest profile for the high point on the embankment
C     WIDTH  - Breadth of the embankment crest
C     APPROC - Elevation of approach for the embankment
C     SURF   - Nature of the embankment surface
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL OLDOFF
      CHARACTER HEAD*80, LINE*80
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL inline
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(A80)
 2    FORMAT(4F10.0,1X,A8)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' ',A80)
 52   FORMAT(' ',F10.1,F10.2,2F10.1,1X,A8)
 54   FORMAT(' *ERR:543* Crest below approach invalid.')
 56   FORMAT(' *ERR:544* Weir width must be positive.')
 58   FORMAT(' *ERR:545* More than ',I5,' offsets for embankment.')
 60   FORMAT(' *ERR:546* Invalid surface option: need PAVED or GRAVEL.')
 64   FORMAT(/,' Minimum crest elevation=',F10.2)
 66   FORMAT(' Minimum crest location=',F10.1)
 68   FORMAT(' *ERR:702* Offset on crest non-increasing at offset=',
     A      F10.2)
C***********************************************************************
C     INPUT THE WEIR DEFINITION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) HEAD
      WRITE(STDOUT,50) HEAD
 
      I = 1
      OLDOFF = -1.E30
      MINCRS = 1.E30
 200  CONTINUE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
        READ(LINE,2,ERR=991) OFF(I), CREST(I), WIDTH(I),
     A                       APPROC(I), SURF(I)
 
        IF(OFF(I).LE.OLDOFF) THEN
          WRITE(STDOUT,68) OFF(I)
          EFLAG = 1
        ENDIF
        OLDOFF = OFF(I)
        IF(APPROC(I).NE.0.0) THEN
          APPROC(I) = APPROC(I) - ZDATUM
        ENDIF
C       PROPAGATE VALUES FROM ABOVE TO SAVE TYPING ON INPUT
        IF(WIDTH(I).EQ.0.0.AND.I.GT.1) WIDTH(I) = WIDTH(I-1)
        IF(APPROC(I).EQ.0.0.AND.I.GT.1) APPROC(I) = APPROC(I-1)
        IF(SURF(I).EQ.'       '.AND.I.GT.1) SURF(I) = SURF(I-1)
        WRITE(STDOUT,52) OFF(I), CREST(I), WIDTH(I), APPROC(I) +
     A                   ZDATUM, SURF(I)
 
C       Adjust for the elevation datum
        CREST(I) = CREST(I) - ZDATUM
 
        IF(CREST(I).LE.APPROC(I)) THEN
          WRITE(STDOUT,54)
          EFLAG = 1
        ENDIF
        IF(WIDTH(I).LE.0.0) THEN
          WRITE(STDOUT,56)
          EFLAG =  1
        ENDIF
        IF(CREST(I).LT.MINCRS) THEN
          MINCRS = CREST(I)
          MINLOC = I
        ENDIF
 
        IF(SURF(I).EQ.'END') THEN
          NOFF = I
          GOTO 210
        ELSE IF(SURF(I).EQ.'PAVED'.OR.SURF(I).EQ.'GRAVEL') THEN
          I = I + 1
 
          IF(I.GT.PMXOFF) THEN
            WRITE(STDOUT,58) PMXOFF
            I = PMXOFF
            EFLAG = 1
          ENDIF
          GOTO 200
        ELSE
          WRITE(STDOUT,60)
          EFLAG = 1
          SURF(I) = 'GRAVEL'
          I = I + 1
          IF(I.GT.PMXOFF) THEN
            WRITE(STDOUT,58) PMXOFF
            I = PMXOFF
            EFLAG = 1
          ENDIF
          GOTO 200
        ENDIF
 
 210  CONTINUE
 
      WRITE(STDOUT,64) MINCRS + ZDATUM
      WRITE(STDOUT,66)  OFF(MINLOC)
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.' 
      END
C
C
C
      SUBROUTINE   ETABIN
     I                   (STDIN, STDOUT, MFTNUM, FTPNT,
     O                    EFLAG, PLCWTB, GLCWTB, PHCWTB, GHCWTB, PSUBTB,
     O                    GSUBTB)
 
C     + + + PURPOSE + + +
C     Input standard table numbers for weir coef. variation
C     for flow over embankment shaped weirs.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, GHCWTB, GLCWTB, GSUBTB, MFTNUM, PHCWTB, PLCWTB,
     A        PSUBTB, STDIN, STDOUT
      INTEGER FTPNT(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     MFTNUM - Maximum allowed table number
C     FTPNT  - function table pointer giving the table address for each
C              table number.  If the address is zero the table does not
C              exist.
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     PLCWTB - Paved low-head weir coefficient table
C     GLCWTB - Gravel surface low head weir coefficient table
C     PHCWTB - Paved high-head weir coefficient table
C     GHCWTB - Gravel surface high head weir coefficient table
C     PSUBTB - Paved submergence table
C     GSUBTB - Gravel surface submergence table
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'subcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      INTEGER TAB(6)
      CHARACTER LINE*80, NAME*7, TABID(6)*16
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FSBRAT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKTAB, FSBRAT, inline, STRIP_L_BLANKS,
     A         GET_INTERNAL_TAB_NUMBER
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT('0*ERR:593* Weir Coef. Table Id must be given.')
C***********************************************************************
      DO 100 I=1,6
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,'(A7,A)',ERR=991) NAME, TABID(I)
        WRITE(STDOUT,'(1X,A7,A)') NAME, TABID(I)
        CALL STRIP_L_BLANKS(
     M                      TABID(I))
        IF(TABID(I).NE.' ') THEN
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, TABID(I),
     M                                 EFLAG,
     O                                 TAB(I))
        ELSE
          TAB(I) = 0
        ENDIF
        CALL CHKTAB
     I             (2, STDOUT, FTPNT, MFTNUM,
     M              TAB(I),
     O              EFLAG)
        IF(I.LE.4.AND.TAB(I).EQ.0) THEN
          WRITE(STDOUT,52)
          EFLAG = 1
        ENDIF
 100  CONTINUE
      IF(EFLAG.GT.0) RETURN
 
C     TRANSFER TO THE INDIVIDUAL VARIABLES
 
      PLCWTB = TAB(1)
      GLCWTB = TAB(2)
      PHCWTB = TAB(3)
      GHCWTB = TAB(4)
      PSUBTB = TAB(5)
      GSUBTB = TAB(6)
 
C     ESTABLISH THE TAILWATER RATIO AT WHICH SUBMERGENCE BEGINS.
      IF(PSUBTB.GT.0) THEN
        PRAT = FSBRAT(PSUBTB)
      ENDIF
      IF(GSUBTB.GT.0) THEN
        GRAT = FSBRAT(GSUBTB)
      ENDIF
 
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.' 
      END
C
C
C
      SUBROUTINE   EMBANK
     I                   (STDIN, STDOUT, STDTAB, MINQ,
     M                    EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute flow over embankment shaped weirs using the USGS
C     procedure.  Tables required by this routine must have
C     been input using FTABIN before requesting the computation
C     of flow using this routine.

      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      REAL MINQ
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     MINQ   - minimum target Q to define the minimum for auto-arguments
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'embcom.cmn'
      INCLUDE 'embwrq.cmn'
      INCLUDE 'tabupgrade.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, K, N, NFRAC, NFRAC2, NHU, SFLAG, TABLE, TYPE, HFLAG,
     A        LOCATION_HU, LOCATION_PFD, N_GT, N_GT_TWICE, ftpup,
     b        verbose
      REAL EFREE, EU, FDROP(PMXNHU), FREED, HUP, HUVEC(PMXNHU), MFREE,
     A     PFDVEC(PMXFRC), POWER, Q(PMXNHU,PMXFRC), QFREE,
     B     QSBVEC(PMXFRC), ZDATUM, LIPREC, MINPFD, XMID, FMID,
     C     GLOBAL_ERROR, RERR, RMS_ERROR, QD_SPAN, DROP, VSCALE,
     D     HSCALE, CSHIFT, LOCAL_MINQ, zrhufd

      REAL*8 SIESQR, easting, northing

      CHARACTER CHAR5*5, CHAR5A*5, HEAD*80, LABEL*50, LINE*120,
     A          CHAR6*6, CQ*8, TABID*16, RLINE*80,
     b  zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC FLOAT
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      CHARACTER GET_TABID*16
      EXTERNAL ETABIN, FRFEMB, inline, INPRDP, KIL, SBFEMB, TWDOUT,
     A         VAR_DECIMAL, LENSTR, GET_EMBK_ITEMS, SET_EMBK_ITEMS,
     B          GET_TABID
 
C     + + + INPUT FORMATS + + +
 8    FORMAT(A6,1X,F10.0)
 16   FORMAT(A5,1X,F10.0)
 26   FORMAT(A5,1X,I5)
 
C     + + + OUTPUT FORMATS + + +
 22   FORMAT('TABID=',A)
 28   FORMAT(F10.5,F10.3)
 29   FORMAT('TYPE=   -2')
30    FORMAT(/,' Datum for heads (minimum crest elevation)=',F10.3,
     A         ' after scaling and possible crest shift.')
 50   FORMAT(' TabId= ',A,' TYPE=',I5,' H/L Crit ratio=',F7.2,
     A   ' H/L max ratio=',F7.2,/,4X,'Hor. scale factor=',F10.3,
     B   ' Vert. scale factor=',F10.3,' Crest shift=',F9.3,/, 4X,
     C   ' Minimum flow target=',F10.3 )
 53   FORMAT(1X,' PartialFD      Drop  Dns Head      Flow  RelError')
 55   FORMAT(1X,F10.6,F10.4,F10.4,2X,A8,F10.3)
 60   FORMAT('   Free flow interpolation between head=',F8.3,
     A       ' and ',F8.3,/,'   has estimated maximum relative',
     B       ' error of',F6.3)
 61   FORMAT(/,' Maximum estimated relative error=',F8.3,' is in free',
     A         ' flow',/,5X,' between heads=',F8.3,' and ',F8.3)
 62   FORMAT(/, ' Estimated root-mean-square error=',F8.3)
 63   FORMAT(/,' Maximum estimated relative error=',F8.3,' is in ',
     A         'submerged flow',/,5X,'at ups. head=',F8.3,
     B         ' between PFD=',F8.4,' and ',F8.4)
 64   FORMAT(/,' ',F5.2,' of checked points had error > LIPREC and',
     A         F5.2,' had error > 2*LIPREC.')
 66   FORMAT(' ',A5,'=',F10.2)
 71   FORMAT('0 Table type 5 replaced by type 13.')
 79   FORMAT(/,' Requested linear interpolation precision=',F5.3)
 80   FORMAT(/,' Minimum partial free drop=',F8.3)
 81   FORMAT(/,' Upstream Head=',F9.4,'  Free Flow=',F10.3)
 90   FORMAT('  Processing EMBANKQ TabId= ',A)
 92   format('; Flow defining minimum head=',f10.3)
 94   FORMAT(/,1X,A5,'=',I5)
C***********************************************************************
      CALL GET_EMBK_ITEMS(STDIN, STDOUT,
     M                    EFLAG)
      CALL SET_EMBK_ITEMS(
     M             EFLAG,
     O             TABLE, TYPE, HLCRIT, HLMAX, HSCALE, VSCALE, CSHIFT,
     O             LOCAL_MINQ)

      IF(LOCAL_MINQ.EQ.0.0) THEN
C       Use the global default value because a minimum Q of zero is 
C       invalid.
        LOCAL_MINQ = MINQ
      ENDIF
      IF(TYPE.EQ.5) THEN
        WRITE(STDOUT,71)
        TYPE = 13
      ENDIF
        IF(HLCRIT.LT.0.0) THEN
          WRITE(STDOUT,*) ' *ERR:611* H/L CRIT RATIO < 0.0 IN EMBANK.'
          HLCRIT = 0.15
        ENDIF
        IF(HLMAX.LT.0.0) THEN
          WRITE(STDOUT,*) ' *ERR:612* H/L MAX RATIO < 0.0 IN EMBANK.'
          HLMAX = 0.32
        ENDIF
      TABID = GET_TABID(TABLE)
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), 
     A                  TYPE, HLCRIT, HLMAX, HSCALE, VSCALE, CSHIFT,
     B                  LOCAL_MINQ
      write(stdtab,92) local_minq
 
      WRITE(*,90) TABID
 
C     MAKE SURE TABLE NUMBER IS NOT ALREADY USED IN THIS INPUT
 
      IF(FTPNT(TABLE).NE.0) CALL TAB_IN_USE
     I                                     (TABID,
     M                                       EFLAG)
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)
 
 
C     INPUT THE STANDARD TABLE NUMBERS
 
      CALL ETABIN
     I           (STDIN, STDOUT, MFTNUM, FTPNT,
     O            EFLAG, PLCWTB, GLCWTB, PHCWTB, GHCWTB, PSUBTB, GSUBTB)
 
      IF(PSUBTB.EQ.0.OR.GSUBTB.EQ.0) THEN
        SFLAG = 0
      ELSE
        SFLAG = 1
      ENDIF
 
C     INPUT THE LABEL FOR THE TABLE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(6X,A50)',ERR=991) LABEL
      WRITE(STDOUT,'('' LABEL='',A50)') LABEL
 
C     INPUT THE WEIR DEFINITION
      ZDATUM = 0.0
      CALL INPRDP
     I           (STDIN, STDOUT, ZDATUM,
     O            EFLAG, NOFF, MINCRS, MINLOC, OFF, CREST, WIDTH,
     O            APPROC, SURF)
 

C     Adjust the crest and approach elevation to force local datum
C     to MINCRS.  
C     Adjust for shift in crest
      MINCRS = MINCRS + CSHIFT      
      DO 100 I=1,NOFF
        CREST(I) = CREST(I) + CSHIFT
        CREST(I) = CREST(I) - MINCRS
        APPROC(I) = APPROC(I) - MINCRS
C       Now apply the scale factors
        CREST(I) = CREST(I)*VSCALE
        APPROC(I) = APPROC(I)*VSCALE
        OFF(I) = OFF(I)*HSCALE
        WIDTH(I) = WIDTH(I)*HSCALE
100   CONTINUE
      ZDATUM = MINCRS
      ZDATUM = VSCALE*ZDATUM
      MINCRS = 0.0

      WRITE(STDOUT,30) ZDATUM

C     Input the upstream head sequence to use - piezometric head
 
C     Clear the flag to signal that the line with the first head
C     is already in hand.
      HFLAG = 0 

C     Set defaults for NFRAC and POWER
      NFRAC = 21
      POWER = 2.0

C     Clear the linear interpolation precision to use as a flag
      LIPREC = 0.0
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(A80)',ERR=991) HEAD
      WRITE(STDOUT,'(1X,A80)') HEAD
 
C     Check for the input option. Next line may contain the first
C     of the heads or it may contain a definition of the number of
C     fractions of free drop TO USE or it may contain linear
C     interpolation precision.
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IF(LINE(1:5).EQ.'NFRAC') THEN
C       We have an explicit specification of the number of fractions
C       of free drop to use.
        READ(LINE,26,ERR=991) CHAR5, NFRAC
        WRITE(STDOUT,94) CHAR5, NFRAC
 
C       Get the power to use
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,16,ERR=991) CHAR5, POWER
        WRITE(STDOUT,66) CHAR5,POWER
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        IF(LINE(1:6).EQ.'LIPREC') THEN
C         The user wants table optimization.
          READ(LINE,8, ERR=991) CHAR6, LIPREC
          WRITE(STDOUT,79) LIPREC
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
          READ(LINE,8, ERR=991) CHAR6, MINPFD
          WRITE(STDOUT,80) MINPFD
        ELSE
C         We have the first line for head.
          HFLAG = 1
        ENDIF        
      ELSEIF(LINE(1:6).EQ.'LIPREC') THEN
C       The user  wants table optimization.
        READ(LINE,8, ERR=991) CHAR6, LIPREC
        WRITE(STDOUT,79) LIPREC
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,8, ERR=991) CHAR6, MINPFD
        WRITE(STDOUT,80) MINPFD
      ELSE
C       We have the first line for head.
        HFLAG = 1
      ENDIF

      I = 1 
 300  CONTINUE
        IF(HFLAG.EQ.0) THEN
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
        ENDIF
C       Only skip reading the first time if HFLAG = 1
        HFLAG = 0
        READ(LINE,'(F10.0)',ERR=991) HUVEC(I)
        WRITE(STDOUT,'('' '',F10.2)') HUVEC(I)
        IF(HUVEC(I).LE.0.0) THEN
          NHU = I - 1
          GOTO 310
        ELSE
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
 
C     Apply scale factor to heads.
      DO 312 I=1,NHU
        HUVEC(I) = HUVEC(I)*VSCALE
312   CONTINUE

C     All items have been read-  start the computations.
C     For each upstream head compute the discharge and also

C     SET THE ROADWAY MOMENTUM FLUX FACTOR TO 1.0. NEEDED BY
C     FRFEMB BUT MOMENTUM FLUX NOT USED FOR COMMAND EMBANK.
 
      RMFFAC = 1.0


C     Compute the submerged flow if requested
 
 
      IF(NFRAC.GT.PMXFRC-1) THEN
        WRITE(STDOUT,
     A'('' *ERR:549* MORE THAN '',I5,'' FRACTIONS OF FREE DROP'')')
     B       PMXFRC
 
        EFLAG = 1
      ENDIF
      IF(EFLAG.GT.0) RETURN


      IF(LIPREC.EQ.0.0) THEN
C       COMPUTE THE PROPORTIONS OF FREE DROP TO USE.
 
        DO 500 I=1,NFRAC
          PFDVEC(I) =(FLOAT(I-1)/FLOAT(NFRAC-1))**POWER
 500    CONTINUE
      ELSE
C       Redefine the head sequence and the sequence of
C       partial free drops. 
        CALL EMBANKQ_OPT
     I                  (STDOUT, LIPREC, MINPFD, SFLAG, LOCAL_MINQ,
     M                   NHU, HUVEC, NFRAC, 
     O                   PFDVEC, EFLAG)
      ENDIF

C     Expand the PFD sequence to include the points needed
C     for estimationg interpolation error. 
      IF(SFLAG.GT.0) THEN
        NFRAC2 = NFRAC + NFRAC - 2
        IF(NFRAC2.GT.PMXFRC-1) THEN
          WRITE(STDOUT,
     A'('' *ERR:549* MORE THAN '',I5,'' FRACTIONS OF FREE DROP'')')
     B         PMXFRC
         
          EFLAG = 1
          RETURN
        ENDIF
        K = NFRAC2
        DO 510 I=NFRAC,3,-1
C         Move the I-th point to its new location.
          PFDVEC(K) = PFDVEC(I)
C         Insert the new point
          PFDVEC(K-1) = 0.5*(PFDVEC(K) + PFDVEC(I-1))
C         Point to the next new location
          K = K - 2
510     CONTINUE
      ELSE
        NFRAC2 = NFRAC
      ENDIF
      GLOBAL_ERROR = 0.0
      LOCATION_HU = 0
      LOCATION_PFD = 0            
      N_GT = 0
      N_GT_TWICE = 0 
      QD_SPAN = 0.0
      SIESQR = 0.D0
 
      DO 5000 I=1, NHU
C       CLEAR FLAGS TO RESTRICT WARNINGS TO ONE MESSAGE PER UPSTREAM
C       HEAD
 
        HLFLAG = 0
        HPFLAG = 0
 
C       COMPUTE A FLOW FOR EACH OF THE GIVEN HEAD VALUES
 
        QFREE = 0.0
        HUP = HUVEC(I)

C       Set the flow at zero partial free drop to 0.0
        Q(I,1) = 0.0

        IF(I.EQ.1) THEN
 
          CALL FRFEMB
     I              (HUP, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I               OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I               HLMAX,
     M               HLFLAG, HPFLAG,
     O               EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O               QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O               APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O               QFREE, MFREE, EFREE)
 
          Q(I,NFRAC) = QFREE
        ELSE
C         Compute an intermediate value for error estimation.
          XMID = 0.5*(HUP + HUVEC(I-1))
          CALL FRFEMB
     I           (XMID, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I            OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I            HLMAX,
     M            HLFLAG, HPFLAG,
     O            EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O            QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O            APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O            FMID, MFREE, EFREE)
 
          CALL FRFEMB
     I              (HUP, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I               OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I               HLMAX,
     M               HLFLAG, HPFLAG,
     O               EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O               QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O               APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O               QFREE, MFREE, EFREE)
          Q(I,NFRAC) = QFREE

          RERR = (0.5*(QFREE + Q(I-1,NFRAC)) - FMID)/FMID
          QD_SPAN = QD_SPAN + HUP - HUVEC(I-1) 
          SIESQR = SIESQR + 0.6666667*RERR**2*(HUP - HUVEC(I-1))
          WRITE(STDOUT,60) HUVEC(I-1), HUP, RERR
          IF(ABS(RERR).GT.LIPREC.AND.LIPREC.GT.0.0) THEN
            N_GT = N_GT + 1
            IF(ABS(RERR).GT.2.*LIPREC) THEN
              N_GT_TWICE = N_GT_TWICE + 1
            ENDIF
          ENDIF
          IF(ABS(RERR).GT.GLOBAL_ERROR) THEN
            GLOBAL_ERROR = ABS(RERR)
            LOCATION_HU = I-1
            LOCATION_PFD = NFRAC2
          ENDIF
        ENDIF

        WRITE(STDOUT,81) HUP, QFREE


        IF(SFLAG.GT.0) THEN 
          WRITE(STDOUT,53)
 
          CALL SBFEMB
     I               (NFRAC2, EU, HUP, MINCRS, MINLOC, PLCWTB, GLCWTB,
     I                PHCWTB, GHCWTB, PSUBTB, GSUBTB, NOFF, SURF,
     I                HLCRIT, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     I                QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, APPL, APPM,
     I                APPR, WL, WM, WR, QFREE, PFDVEC,
     O                FREED, QSBVEC)

          FDROP(I) = FREED
          CALL VAR_DECIMAL(QSBVEC(NFRAC2),
     O                       CQ)
          WRITE(STDOUT,55) PFDVEC(NFRAC2), FREED, HUP - FREED, 
     A                     CQ
C        Note that NFRAC2 is always an even integer!
          Q(I,NFRAC2/2+1) = QSBVEC(NFRAC2)
          DO 600 J=NFRAC2-2,2,-2
            Q(I,J/2+1) = QSBVEC(J)

            RERR = (0.5*(QSBVEC(J) + QSBVEC(J+2))/QSBVEC(J+1) - 1.0) 
C           Use XMID as a temp variable.
            XMID = FREED*(PFDVEC(J+2) - PFDVEC(J)) 
            QD_SPAN = QD_SPAN + XMID
            SIESQR = SIESQR + 0.6666667*RERR**2*XMID
            DROP = FREED*PFDVEC(J)
            CALL VAR_DECIMAL(QSBVEC(J),
     O                         CQ)
            WRITE(STDOUT,55) PFDVEC(J), DROP, HUP - DROP, 
     A                        CQ, RERR
            IF(ABS(RERR).GT.LIPREC.AND.LIPREC.GT.0.0) THEN
              N_GT = N_GT + 1
              IF(ABS(RERR).GT.2.*LIPREC) THEN
                N_GT_TWICE = N_GT_TWICE + 1
              ENDIF
            ENDIF
            IF(ABS(RERR).GT.GLOBAL_ERROR) THEN
              GLOBAL_ERROR = ABS(RERR)
              LOCATION_HU = I
              LOCATION_PFD = J
            ENDIF
 600      CONTINUE
        ENDIF
 5000 CONTINUE

      IF(LOCATION_PFD.EQ.NFRAC2) THEN
        WRITE(STDOUT,61) GLOBAL_ERROR, HUVEC(LOCATION_HU),
     A     HUVEC(LOCATION_HU+1)
      ELSE
        WRITE(STDOUT,63) GLOBAL_ERROR, HUVEC(LOCATION_HU),
     A     PFDVEC(LOCATION_PFD), PFDVEC(LOCATION_PFD+2)
      ENDIF
      RMS_ERROR = SQRT(SIESQR/QD_SPAN)
      WRITE(STDOUT,62) RMS_ERROR
      IF(LIPREC.GT.0.0) THEN
        IF(SFLAG.GT.0) THEN
          N = NHU - 1 + NHU*(NFRAC - 1)
        ELSE
          N = NHU - 1
        ENDIF
    
        WRITE(STDOUT,64) FLOAT(N_GT)/N, FLOAT(N_GT_TWICE)/N
      ENDIF
        
      IF(SFLAG.EQ.1) THEN
C       Squeeze out the extra PFD values 
        DO 5010 I=4,NFRAC2,2
          PFDVEC(I/2+1) = PFDVEC(I)
5010    CONTINUE

 
C       NOW OUTPUT THE TABLE TO STDTAB
 
        zrhufd = 0.0
        CALL TWDOUT
     I             (STDOUT, STDTAB, TABLE, LABEL, NHU, NFRAC, HUVEC,
     I              FDROP, PFDVEC, Q, MINCRS+ZDATUM,
     I              TYPE, ' EMBANKQ', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)

        if(twod_cubic_out.eq.'YES') then
          verbose = 1
          CALL twodfit
     I             (STDOUT, TABLE, NHU, NFRAC, HUVEC,
     I              FDROP, PFDVEC, Q, MINCRS+ZDATUM,
     I              TYPE, ' EMBANKQ', zrhufd, verbose,
     M              FTP,
     O              EFLAG,ftpup)
        endif
      ELSE
        WRITE(STDTAB,22)  TABID(1:LENSTR(TABID))
        WRITE(STDTAB,29)
        WRITE(STDTAB,'(''REFL= 0.0'')')
        WRITE(STDTAB,'(''      HEAD DISCHARGE'')')
        WRITE(STDTAB,28) 0.0, 0.0
        DO 2300 I=1,NHU
          WRITE(STDTAB,28) HUVEC(I), Q(I,NFRAC)
 2300   CONTINUE
        WRITE(STDTAB,28) -1.0, 0.0
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
      SUBROUTINE EMBANKQ_OPT
     I                      (STDOUT, LIPREC, MINPFD, SFLAG, MINQ,
     M                       NHU, HUVEC, NFRAC, 
     O                       PFDVEC, EFLAG)

C     Compute a sequence for upstream head and partial free drops
C     to produce an optimized table for EMBANKQ.

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NFRAC,  NHU, STDOUT, SFLAG
      REAL  LIPREC, MINPFD, MINQ, HUVEC(PMXNHU), PFDVEC(PMXFRC)

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     LIPREC - requested linear interpolation precision
C     MINPFD - minimum value of partial free drop in the 2-D table
C     SFLAG - =1 if submerged flows are to be computed, =0 otherwise
C     MINQ  - target minimum flow to re-define minimum head
C     NHU - number of upstream heads.
C     HUVEC - upstream head sequence
C     NFRAC - number of partial free drops
C     PFDVEC - sequence of partial free drops
C     EFLAG - error flag
      

C     + + + COMMON BLOCKS + + +
      INCLUDE 'ftable.cmn'
      INCLUDE 'embcom.cmn'
      INCLUDE 'embwrq.cmn'

C     + + + LOCAL VARIABLES + + +

      INTEGER  I, J, TABLT, TABGT, N, NFRAC_OLD
      REAL ARGRAT, DFROW, DFCOL, HUP, EFREE, QFREE,  MFREE, 
     A     EU, FREED, QSBVEC(PMXNHU), H1, H2, Q1, Q2, B

      REAL*8  XVEC(PMXNHU), FVEC(PMXNHU), BREAK_POINTS(PMXNHU), 
     A        MAX_RERR

      CHARACTER*16 TEMP

C     Call sub programs
      EXTERNAL GET_INTERNAL_TAB_NUMBER

C     + + + OUTPUT FORMATS + + +
 70   FORMAT(/,' EMBANKQ will use ',I5,' upstream heads for the 2-D',
     A         ' table.',/,'  The maximum estimated relative ',
     B         'interpolation error is:',F10.3)
 71   FORMAT(/,' EMBANKQ will use ',I5,' partial free drops for ',
     A         'the 2-D table.',/,'  The maximum estimated',
     B         ' relative interpolation error is:',F10.3)
 72   FORMAT(/,' Using upstream head=',F8.3,' to define partial',
     A         ' free drops.')
73    FORMAT(/,' Minimum head revised to: ',F10.3,
     A         ' for target flow=',F10.3)
98    FORMAT(1PE14.6,1PE14.6)
 99   FORMAT(/,' *ERR:635* Interpolation precision tables missing.',
     A   ' Check for version',/,11X,'of file: TYPE5.TAB.')
C***********************************************************************
C     Define the power-function interpolation precision tables.

      TEMP = '10001'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TEMP,
     M                             EFLAG,
     O                             TABLT)
      TEMP = '10002'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TEMP,
     M                             EFLAG,
     O                             TABGT)
      TABLT = FTPNT(TABLT)
      TABGT = FTPNT(TABGT)
      IF(TABLT.LT.1.OR.TABGT.LT.1) THEN
        WRITE(STDOUT,99) 
        EFLAG = 1
        RETURN
      ENDIF

C     Save entry value of NFRAC for lower limit of point set.
      NFRAC_OLD = NFRAC

C     Define a better minimum head.  Compute free flow at .5 and 1.0
C     times the user-given minimum head value.  Find the power in
C     a simple power-function fit to these two points.    Then compute
C     the min-head value that would give a desired minimum flow.
C     The desired minimum flow is under user control but has a 
C     default value. 
      H2 = HUVEC(1)
      H1 = H2/2.0
      CALL FRFEMB
     I           (H1, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I            OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I            HLMAX,
     M            HLFLAG, HPFLAG,
     O            EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O            QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O            APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O            Q1, MFREE, EFREE)
      CALL FRFEMB
     I           (H2, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I            OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I            HLMAX,
     M            HLFLAG, HPFLAG,
     O            EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O            QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O            APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O            Q2, MFREE, EFREE)

C     Compute the power of the simple power function that fits
C     these two points (and (0,0) as well).
      B = LOG(Q1/Q2)/LOG(H1/H2)
      H1 = H2*(MINQ/Q2)**(1.0/B)
      WRITE(STDOUT,73) H1, MINQ
      HUVEC(1) = H1      

C     Define the point set for computing the free flows to fit with
C     a cubic spline. Find the argument ratio, ARGRAT, for a power
C     of 2.5, close to the maximum for EMBANKQ.

      CALL TDLK10
     I           (STDOUT, TABGT, 10, 2.5, LIPREC,
     O            ARGRAT, DFROW, DFCOL)

C     Use the maximum and minimum head from the user and ARGRAT to
C     compute the number of heads to use and compute a new head
C     sequence.
      N = INT(LOG(HUVEC(NHU)/HUVEC(1))/LOG(ARGRAT) + 1.0) + 1
      IF(N.LT.NFRAC_OLD) N = NFRAC_OLD
      CALL RATIOPNT
     I             (N, DBLE(HUVEC(1)), DBLE(HUVEC(NHU)), 
     O                    XVEC)
      NHU = N
C     Now compute the free flows for this sequence of heads. Set flags
C     to suppress warning messages. 

      HLFLAG = 1
      HPFLAG = 1
      RMFFAC = 1.0
      DO 110 I=1,NHU

        HUP = REAL(XVEC(I))

         CALL FRFEMB
     I              (HUP, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I               OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I               HLMAX,
     M               HLFLAG, HPFLAG,
     O               EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O               QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O               APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O               QFREE, MFREE, EFREE)

        IF(EFLAG.NE.0) THEN
          RETURN
        ENDIF
        FVEC(I) = DBLE(QFREE)
110   CONTINUE

C      WRITE(STDOUT,*) ' DUMP OF POINTS FOR TESTING'
C      WRITE(STDOUT,*)' NHU=',NHU
C      DO 111 I=1,NHU
C        WRITE(STDOUT,98) XVEC(I), FVEC(I)
C111   CONTINUE

C     Try to find improved breakpoints. 

      CALL FINDBRK
     I            (STDOUT, NHU, XVEC, FVEC, DBLE(LIPREC),
     I             1, 1, 
     O             N, BREAK_POINTS, EFLAG, MAX_RERR)

      IF(EFLAG.NE.0) THEN
        RETURN
      ENDIF
      NHU = N
      DO 120 I=1,N
        HUVEC(I) = REAL(BREAK_POINTS(I))
120   CONTINUE

      WRITE(STDOUT,70) NHU, MAX_RERR
      IF(SFLAG.EQ.0) RETURN


C     Select a head near the middle of the vector and compute the 
C     free and submerged flows for that head to define the 
C     basis for finding good breakpoints for the partial free 
C     drops.  
      I = FLOAT(2*NHU)/3. + 1.
      IF(I.GT.NHU) I = NHU
      HUP = HUVEC(I)
      CALL FRFEMB
     I           (HUP, MINCRS, PLCWTB, GLCWTB, PHCWTB, GHCWTB, NOFF,
     I            OFF, CREST, WIDTH, APPROC, SURF, RMFFAC, HLCRIT,
     I            HLMAX,
     M            HLFLAG, HPFLAG,
     O            EU, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     O            QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, YFL, YFM, YFR,
     O            APPL, APPM, APPR, WL, WM, WR, AELL, AELM, AELR,
     O            QFREE, MFREE, EFREE)
      WRITE(STDOUT,72) HUP

      
C     Define the sequence of partial free drops to use.
      
      CALL TDLK10
     I           (STDOUT, TABLT, 10, 0.5, LIPREC,
     O            ARGRAT, DFROW, DFCOL)

C     Use the maximum PFD, 1.0, and the MINPFD to compute
C     the number of PFD's.

      N = INT(LOG(1.0/MINPFD)/LOG(ARGRAT) + 1.0) + 1
      IF(N.LT.NFRAC_OLD) N = NFRAC_OLD
      CALL RATIOPNT
     I             (N, DBLE(MINPFD), 1.D0, 
     O                    XVEC)
C     Note: XVEC does not contain the zero point; PFDVEC does!
C     History has come back to haunt us!
      NFRAC = N + 1
      PFDVEC(1) = 0.0
      DO 115 I=2,NFRAC
        PFDVEC(I) = REAL(XVEC(I-1))
115   CONTINUE
      CALL SBFEMB
     I           (NFRAC, EU, HUP, MINCRS, MINLOC, PLCWTB, GLCWTB,
     I            PHCWTB, GHCWTB, PSUBTB, GSUBTB, NOFF, SURF,
     I            HLCRIT, XRDFL, XRDFR, HRDFL, HRDFM, HRDFR, QRDFL,
     I            QRDFM, QRDFR, TOTHL, TOTHM, TOTHR, APPL, APPM,
     I            APPR, WL, WM, WR, QFREE, PFDVEC,
     O            FREED, QSBVEC)
      DO 116 I=NFRAC,2,-1
        FVEC(I-1) = QSBVEC(I)
116   CONTINUE

C     Try to find improved breakpoints. 

      CALL FINDBRK
     I            (STDOUT, N, XVEC, FVEC, DBLE(LIPREC),
     I             1, 1, 
     O             NFRAC, BREAK_POINTS, EFLAG, MAX_RERR)
      IF(EFLAG.NE.0) THEN
        WRITE(STDOUT,*) ' FINDBRK failure'
        RETURN
      ENDIF
      NFRAC = NFRAC + 1
      PFDVEC(1) = 0.0
      DO 150 J=2,NFRAC
        PFDVEC(J) = BREAK_POINTS(J-1)
150   CONTINUE
      WRITE(STDOUT,71) NFRAC, MAX_RERR

      RETURN
      END
C
C
C 
      SUBROUTINE SET_EMBK_ITEM_DEFAULTS()

C     Set the default values in the vectors used to
      IMPLICIT NONE
      INCLUDE 'embkitm.cmn'

C***********************************************************************
C     Default for: TABID 
      EMBKITMCTAB(  1) = '    '           
C     Default for: TABLE - note # is ignored in the standard scanner 
      EMBKITMCTAB(  2) = '    '           
C     Default for: TYPE
      EMBKITMITAB(  1) = 13
C     Default for: HLCRIT
      EMBKITMFTAB(  2) = 0.15
C     Default for: HLMAX
      EMBKITMFTAB(  3) = 0.32
C     Default for: HSCALE
      EMBKITMFTAB(  4) = 1.0
C     Default for: VSCALE
      EMBKITMFTAB(  5) = 1.0
C     Default for: CSHIFT
      EMBKITMFTAB(  6) = 0.0
C     Default for: LOCAL_MINQ
      EMBKITMFTAB(  7) = 0.0

      RETURN
      END
C
C
C
      SUBROUTINE  SET_EMBK_ITEMS(
     M             EFLAG,
     O             TAB, TYPE, HLCRIT, HLMAX, HSCALE, VSCALE, CSHIFT,
     O             LOCAL_MINQ)

C     Set items in EMBANK
C     All values not set explicitly by user are at their default value.

      IMPLICIT NONE

      INTEGER TAB, TYPE, GETQ, GETY2, EFLAG
      REAL HLCRIT, HLMAX, HSCALE, VSCALE, CSHIFT, LOCAL_MINQ
    

C     Local
      CHARACTER*16 KEY1, KEY2

      INCLUDE 'embkitm.cmn'
      INCLUDE 'stdun.cmn'
C***********************************************************************
C     Set the table id.  There  are two strings allowed as the variable name:
C     TABID or TABLE
      KEY1 = EMBKITMCTAB(  1)
      KEY2 = EMBKITMCTAB(  2)
      IF(KEY1.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY1,
     M                                EFLAG,
     O                                TAB)
      ELSEIF(KEY2.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY2,
     M                                EFLAG,
     O                                TAB)
      ELSE
        TAB = 0
      ENDIF
    
C     Set the value for TYPE
      TYPE =   EMBKITMITAB( 1)
C     Set the value for HLCRIT
      HLCRIT = EMBKITMFTAB( 2)
C     Set the value for HLMAX
      HLMAX =  EMBKITMFTAB( 3)   
C     Set the value for HSCALE
      HSCALE = EMBKITMFTAB( 4)
C     Set the value for VSCALE
      VSCALE = EMBKITMFTAB( 5)
C     Set the value for CSHIFT
      CSHIFT = EMBKITMFTAB( 6)
C     Set the value for the local MINQ
      LOCAL_MINQ = EMBKITMFTAB(7)
    
      RETURN
      END
C
C
C
      SUBROUTINE GET_EMBK_ITEMS(STDIN, STDOUT,
     M                   EFLAG)

C     Get the table id and various options for EMBANKQ command

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG

      INCLUDE 'embkitm.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  INTVAL, REAVAL, 
     A         CHRVAL, DPRVAL, EXACT, LOWER,
     B         NUMERIC, CHAR, N_SYMBOL
      PARAMETER(N_SYMBOL=9, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, 
     B          EXACT=0,LOWER=1, NUMERIC=0, CHAR=1)

      INTEGER MAX_LINE
     A        

      EXTERNAL GET_NAMED_ITEMS, SET_EMBK_ITEM_DEFAULTS

C     + + + SAVED VALUES + + +
      INTEGER GROUP(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL), GROUP_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE  SYMBOL_TABLE, GROUP, RESPONSE_TYPE, CONVERT_RULE,
     A      GROUP_INDEX

      DATA  SYMBOL_TABLE /
     *'TABID','TABLE','TYPE','HLCRIT','HLMAX','HSCALE','VSCALE',
     A 'CSHIFT','MINQ'/
                                                                      
      DATA GROUP  /
     *CHAR, CHAR, 7*NUMERIC/          

      DATA GROUP_INDEX /
     *1, 2, 1, 2, 3, 4, 5, 6, 7/
                                       
      DATA RESPONSE_TYPE  /
     *CHRVAL,CHRVAL,INTVAL, 6*REAVAL/
    
      DATA CONVERT_RULE /
     *LOWER,LOWER,EXACT,6*LOWER/
   
C***********************************************************************
C     Set Defaults
 
      CALL SET_EMBK_ITEM_DEFAULTS()

      MAX_LINE = 1
      CALL GET_NAMED_ITEMS(
     I                     STDIN, STDOUT,  MAX_LINE, N_SYMBOL,
     I  GROUP, RESPONSE_TYPE, CONVERT_RULE, GROUP_INDEX, SYMBOL_TABLE,
     I  MAXR_EMBKITM, MAXDP_EMBKITM, MAXC_EMBKITM, 'EMBANKQ items',
     O  EMBKITMITAB, EMBKITMFTAB, EMBKITMDTAB, EMBKITMCTAB,
     O  EFLAG)
      
      RETURN

      END
