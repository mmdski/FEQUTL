C
C
C
      DOUBLE PRECISION FUNCTION   FFPRES
     I                                  (Y)
 
C     + + + PURPOSE + + +
C     Compute the free flow profile residual.  Flow at downstream end of
C     the channel is at critical depth.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      DOUBLE PRECISION Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'chncom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL QC
      DOUBLE PRECISION RESULT, TEMP
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE
 
C     + + + EXTERNAL FUNCTIONS + + +
      DOUBLE PRECISION FDXDY
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FDXDY, LKTQC, QFUN
C     ***********************************FORMAT*************************
50    FORMAT(' FFPRES:YL=',F10.4,' YR=',F10.4,' RESULT=',1PE12.4,
     A       ' FFPRES=',1PE12.4,' QERR=',1PE12.4,' Q=',1PE12.4)
C***********************************************************************
C      WRITE(STD6,*) ' ENTERING FFPRES: Y=',Y
 
C     FIND CRITICAL FLOW AT DOWNSTREAM END OF THE CHANNEL
 
      YR = Y
      CALL LKTQC
     I          (XSADR,
     M           YR,
     O           QC)
      Q = QC
 
C     COMPUTE THE LENGTH OF THE PROFILE FROM DEPTH YL TO DEPTH YR
C     AND SUBTRACT THE CHANNEL LENGTH, L(FROM CHANCOM).
C     NOTE: YR < YL IN ALL CASES.
 
      CALL QFUN
     I         (YL, YR, FDXDY,
     O          QERR, RESULT)
      TEMP = RESULT - DBLE(L)  
      FFPRES = TEMP
      LSTRES = TEMP
C      WRITE(STD6,50) YL, YR, RESULT, FFPRES, QERR, Q
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION   SBPRES
     I                                  (FLOW)
 
C     + + + PURPOSE + + +
C     Compute the submerged flow profile residual.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      DOUBLE PRECISION FLOW
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     FLOW   - Flowrate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'chncom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      DOUBLE PRECISION TEMP
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE
 
C     + + + EXTERNAL FUNCTIONS + + +
      DOUBLE PRECISION FDXDY
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FDXDY, QFUN
C***********************************************************************
C      WRITE(STD6,*) ' ENTERING SBPRES: FLOW=',FLOW,' YL=',YL,' YR=',YR
 
C     YL AND YR ARE FIXED AND IN CHNCOM.
 
      Q = FLOW
 
C     COMPUTE THE LENGTH OF THE PROFILE FROM DEPTH YL TO DEPTH YR
C     AND SUBTRACT THE CHANNEL LENGTH, L(FROM CHNCOM).
 
      CALL QFUN
     I         (YL, YR, FDXDY,
     O          QERR, TEMP)
      TEMP = TEMP - DBLE(L)
      SBPRES = TEMP
      LSTRES = TEMP
C      WRITE(STD6,*) ' SBPRES=',SBPRES
 
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION   FDXDY
     I                                 (Y)
 
C     + + + PURPOSE + + +
C     Compute the inverse water surface slope for steady flow in a
C     prismatic channel.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      DOUBLE PRECISION Y
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     Y      - maximum depth in a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'chncom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL A, ALPHA, BETA, DALPHA, DBETA, DK, DT, J, K, QC, T, TP
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE, SNGL
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKT22
C***********************************************************************
      TP = SNGL(Y)
      CALL XLKT22
     I           (XSADR,
     M            TP,
     O            A, T, DT, J, K, DK, BETA, DBETA, ALPHA, DALPHA, QC)
 
C      WRITE(STD6,*) 'FDXDY: Y=',Y,' Q=',Q,' QC=',QC,' SBOT=',SBOT,
C     A               ' K=',K     
 
      FDXDY = (1.D0 - (DBLE(Q)/DBLE(QC))**2)/
     A          (DBLE(SBOT) - (DBLE(Q)/DBLE(K))**2)
C     WRITE(STD6,*) ' FDXDY=', FDXDY
 
      RETURN
      END
C
C
C
      SUBROUTINE FIND_ROOT_INTERVAL(STDOUT, F, XMIN, XMAX)

C     Do a detailed search for intervals containing a root. 

C     Use for SBFCHN first and then extend.

      IMPLICIT NONE
      INTEGER STDOUT
      REAL XMIN, XMAX
      REAL*8 F

      EXTERNAL F

C     Local
      INTEGER I, J, N, KNT, ROOT_INTERVAL(500)


      REAL*8 FN, FI, XVEC(5000), FVEC(5000), XTMIN, XTMAX

C     *****************************FORMAT*******************************
50    FORMAT(' NO ROOT INTERVAL FOUND!')
52    FORMAT(' MORE THAN ONE ROOT INTERVAL!',/,
     A  '  KNT=',I5,(8F10.4))
C***********************************************************************
      XTMIN = XMIN
      XTMAX = XMAX*1.01
      IF(XTMIN.EQ.0.D0) THEN
        XTMIN = 0.1*XTMAX
      ENDIF
C     Argument is in cfs.  Assign a minimum of 100 points to any range
C     and a maximum of 2 points per  cfs.
      N = 2.*(XTMAX - XTMIN)
      IF(N.LT.100) N = 100
      IF(N.GT.5000) N = 5000
      FN = DBLE(N-1)
      DO 100 I=1,N
        FI = DBLE(I-1)
        XVEC(I) = XTMIN +   FI*(XTMAX - XTMIN)/FN
        FVEC(I) = F(XVEC(I))
100   CONTINUE

      KNT = 0
      DO 200 I= 2,N
        IF(FVEC(I)*FVEC(I-1).LT.0.D0) THEN
          KNT = KNT + 1
          ROOT_INTERVAL(KNT) = I
        ENDIF
200   CONTINUE

      IF(KNT.EQ.0) THEN
        WRITE(STDOUT,50) 
      ELSEIF(KNT.GT.1) THEN
        WRITE(STDOUT,52) KNT, (XVEC(ROOT_INTERVAL(J)),J=1,KNT)
      ENDIF
      RETURN
      END
      
C
C
C
      SUBROUTINE   SBFCHN
     I                   (STDOUT,
     M                    QRATIO, EFLAG)
 
C     + + + PURPOSE + + +
C     Compute submerged flow through a prismatic channel.
C     most values are in CHNCOM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER STDOUT, EFLAG
      REAL QRATIO 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     QRATIO - ratio of flows used to provide initial estimate.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'chncom.cmn'
      INCLUDE 'epscom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER RFLAG
      REAL  QMAX, QMIN, QSAVE, YNEAR
      DOUBLE PRECISION LHAT, FL, FR, QMID
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE, MAX, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      DOUBLE PRECISION FDXDY, SBPRES
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FDXDY, QFUN, SBPRES, DBLRGF
C     *******************************FORMATS****************************
 51   FORMAT(/,' *BUG in SBFCHN: No sign change on entry to DBLRGF.')
 52   FORMAT(/,' *BUG in SBFFCHN: More than 100 iterations in DBLRGF.')
 53   FORMAT(/,
     A ' *ERR/WRN:762* in SBFCHN: Argument collapse with residual=',
     B   1PE10.3,' more than twice EPSINT in DBLRGF.')
C***********************************************************************
      EFLAG = 0
C      WRITE(STDOUT,*) ' SBFCHN: YL=',YL,' YR=',YR
C     INITIAL FLOWS MUST BE SELECTED WITH CARE IF SBOT > 0.0.  WE
C     MUST MAKE SURE THAT NO INTEGRAL CROSSES THE NORMAL FLOW LEVEL.
C     SINCE THE CHANNEL IS PRISMATIC, NORMAL FLOW EXISTS WHENEVER
C     YL = YR.  IF YR < YL, AND SBOT > 0.0 THEN THE MIMIMUM FLOW IS
C     NORMAL FLOW PLUS SOME SMALL TOLERANCE.  IF YR > YL AND
C     SBOT > 0.0 THEN THE MAXIMUM FLOW IS NORMAL FLOW LESS SOME
C     SMALL TOLERANCE.
 
C     IF SBOT IS <=0.0 NORMAL FLOW DOES NOT EXIST. THUS THE MIMIMUM FLOW
C     IS ESSENTIALLY ZERO AND THE MAXIMUM FLOW THE PREVIOUS FLOW.
C     December 15, 1999- allow maximum flow to be larger than the 
C     previous flow.  Needed in some cases when the critical flow 
C     decreases with increasing stage.  This case may require 
C     detailed study to solve correctly-at some point in the future. 
      QSAVE = Q
      IF(SBOT.LE.0.0) THEN
        QMIN = 0.0
C       Removed Feb. 4, 2000.  Caused more problems than it solved. 
        QMAX = 1.0*Q
      ELSE
 
        IF(YR.LT.YL) THEN
C         CHECK IF THE CHANNEL IS LONG ENOUGH TO ESTABLISH ESSENTIAL
C         NORMAL FLOW AT THE UPSTREAM END FOR THE GIVEN DOWNSTREAN DEPTH.
          YNEAR = MAX((1.0 - NDDREL)*YL, YL - NDDABS)
          Q = QN
          CALL QFUN
     I             (YNEAR, YR, FDXDY,
     O              QERR, LHAT)
 
C          WRITE(STDOUT,*) ' SBFCHN: YNEAR=',YNEAR
C          WRITE(STDOUT,*) ' SBFCHN: LHAT=',LHAT,' L=',L
          IF(LHAT.LT.L) THEN
C           FLOW IS AT NORMAL RATE
            RETURN
          ENDIF
          Q = QSAVE
          QMIN = QN + EPSARG
          QMAX = Q
        ELSE
C         CHECK IF THE CHANNEL IS LONG ENOUGH TO ESTABLISH ESSENTIAL
C         NORMAL FLOW AT THE UPSTREAM END FOR THE GIVEN DOWNSTREAN DEPTH.
          YNEAR = MIN((1.0 + NDDREL)*YL, YL + NDDABS)
          Q = QN
          CALL QFUN
     I             (YNEAR, YR, FDXDY,
     O              QERR, LHAT)
 
C          WRITE(STDOUT,*) ' SBFCHN: YNEAR=',YNEAR
C          WRITE(STDOUT,*) ' SBFCHN: LHAT=',LHAT,' L=',L
          IF(LHAT.LT.L) THEN
C           FLOW IS AT NORMAL RATE
            RETURN
          ENDIF
          Q = QSAVE
          IF(Q.GT.QN) Q = QN - EPSARG
          QMIN = 0.0
          QMAX = Q
        ENDIF
      ENDIF


C      IF(YR.LT.YL) THEN
C        FL = 9.0*DBLE(L)
C        FR = -0.5*DBLE(L)
C      ELSE
C        FL = -0.5*DBLE(L)
C        FR = 9.0*DBLE(L)
C      ENDIF

C     Estimate the root based on past ratio with the previous
C     flow.  
      QMID = QRATIO*QMAX
      IF(QMID.LE.QMIN) THEN
        QMID = 0.125*QMAX + 0.875*QMIN
      ELSEIF(QMID.GE.QMAX) THEN
        QMID = 0.125*QMIN + 0.875*QMAX
      ENDIF
C     Establish sign change for regula falsi.  
C     Does detailed analysis for checking
C      CALL FIND_ROOT_INTERVAL(STDOUT, SBPRES,
C     M                        QMIN, QMAX, QMID,
C     O                        FL, FR)
C      CALL FIND_ROOT_INTERVAL(STDOUT, SBPRES, QMIN, QMAX)
     
C     Establish sign change for regula falsi.  

      IF(YR.LT.YL) THEN
        FL = 9.0*DBLE(L)
        FR = -0.5*DBLE(L)
      ELSE
        FL = -0.5*DBLE(L)
        FR = 9.0*DBLE(L)
      ENDIF

      CALL DBLRGF
     I           (DBLE(EPSARG), EPSINT, SBPRES,
     M            QMIN, QMAX, FL, FR, QMID,
     O            RFLAG)
      IF(RFLAG.GT.0) THEN
        EFLAG = 1
        IF(RFLAG.EQ.1) THEN
          WRITE(STDOUT,51)
        ELSEIF(RFLAG.EQ.2) THEN
          WRITE(STDOUT,52)
        ELSE
          WRITE(STDOUT,53) FL
        ENDIF
      ENDIF
C     Compute the solution ratio to use in estimating the starting
C     point for the next upstream head in the sequence. 
      QRATIO = QMID/QMAX
C      IF(QRATIO.GT.0.98) QRATIO = 0.98
     
      Q = QMID

      RETURN
      END
C
C
C
      SUBROUTINE   FRFCHN
     I                   (STDOUT, HUP, HDATUM, ZBOTL, ZBOTR,
     M                    RATIO,
     O                    EFLAG, QFREE, FDROP)
 
C     + + + PURPOSE + + +
C     Find free flow in a prismatic channel for the CHANRAT command.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDOUT
      REAL FDROP, HDATUM, HUP, QFREE, RATIO, ZBOTL, ZBOTR
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     HUP    - Head upstream
C     HDATUM - Datum for measuring head
C     ZBOTL  - Bottom elevation at left section
C     ZBOTR  - Bottom elevation at right section
C     RATIO  - Estimated ratio of downstream depth to upstream
C              depth.
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     QFREE  - Free flow
C     FDROP  - Free drop value
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'chncom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER RFLAG
      REAL A, ALPHA, BETA, DALPHA, DBETA, DK, DT, J, K, T,
     A     YMAX, YMIN, YNEAR, ZL, ZR
      DOUBLE PRECISION LHAT, FL, FR, YMID
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE, MAX, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      DOUBLE PRECISION FDXDY, FFPRES
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FDXDY, FFPRES, FNDCDE, QFUN, DBLRGF, XLKTAL
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:627* SLOPE OF CHANNEL IS STEEP.')
 51   FORMAT(/,' *BUG in FRFCHN: No sign change on entry to DBLRGF.')
 52   FORMAT(/,' *BUG in FRFCHN: More than 100 iterations in DBLRGF.')
 53   FORMAT(/,' *BUG in FRFCHN: Argument collapse with residual more',
     A         ' than twice EPSINT in DBLRGF.')
C***********************************************************************
C     WRITE(STDOUT,*) ' ENTERING FRFCHN. HUP=',HUP
C     FIND DEPTH AT UPSTREAM END(LEFT END).
      ZL = HUP + HDATUM
      YL = ZL - ZBOTL
C     WRITE(STDOUT,*) ' FRFCHN: YL=',YL
C     DOES NORMAL FLOW EXIST?
      IF(SBOT.GT.0.0) THEN
C       YES.
 
        CALL XLKTAL
     I             (XSADR,
     M              YL,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA, DALPHA)
 
        QN = K*SQRT(SBOT)
 
C       FIND CRITICAL DEPTH AT QN
        YMIN = YL
        CALL FNDCDE
     I             (STDOUT, XSADR, QN,
     M              YMIN)
 
        IF(YMIN.GT.YL) THEN
          WRITE(STDOUT,50)
          EFLAG = 1
          RETURN
        ENDIF
C       WRITE(STDOUT,*) ' FRFCHN. QN=',QN,' CRIT. DEPTH=',YMIN
 
C       DETERMINE IF THE CHANNEL IS LONG ENOUGH TO CAUSE THE
C       FLOW PROFILE WHEN STARTING AT YMIN(THAT IS AT CRITICAL
C       DEPTH FOR THE NORMAL FLOW)  TO ATTAIN CLOSELY TO
C       NORMAL DEPTH( YL IN THIS CASE) BEFORE THE END OF THE
C       CHANNEL IS REACHED.
 
        YNEAR = MAX((1.0 - NDDREL)*YL, YL - NDDABS)
        Q = QN
        CALL QFUN
     I           (YNEAR, YMIN, FDXDY,
     O            QERR, LHAT)
 
C        WRITE(STDOUT,*) ' FRFCHN: YNEAR=',YNEAR
C        WRITE(STDOUT,*) ' FRFCHN: LHAT=',LHAT,' L=',L,' QERR=',QERR
        IF(LHAT.LT.L) THEN
C         THE FLOW IS THE SAME AS NORMAL FLOW.
          YR = YMIN
          ZR = YR + ZBOTR
 
          FDROP = ZL - ZR
          QFREE = Q
          RETURN
        ENDIF
 
C       THE ROOT IS SOMEWHERE IN YMIN < Y < YL.  IF YR = YL THEN
C       THE PROFILE LENGTH, LHAT, IS 0 AND IF YR = YMIN THEN
C       LHAT = INFINITY.
      ELSE
C       IN THIS CASE THERE IS NO NORMAL DEPTH AND THE MINIMUM DEPTH
C       CAN APPROACH ZERO.
        YMIN = EPSARG
      ENDIF
C     MAXIMUM DEPTH AT THE DOWNSTREAM END IS YL.
 
      YMAX = YL
C     Start with rough estimates of the two extremes.
      FL = -DBLE(L)
      FR = 1.5*DBLE(L)
      YMID = RATIO*YL
      IF(YMID.LT.YMIN) THEN 
        YMID = 0.125*YMAX + 0.875*YMIN
      ELSEIF(YMID.GT.YMAX) THEN
        YMID = 0.125*YMIN + 0.875*YMAX
      ENDIF
      CALL DBLRGF
     I           (DBLE(EPSARG), EPSINT, FFPRES,
     M            YMAX, YMIN, FL, FR, YMID, 
     O            RFLAG)
      IF(RFLAG.GT.0) THEN
        EFLAG = 1
        IF(RFLAG.EQ.1) THEN
          WRITE(STDOUT,51)
        ELSEIF(RFLAG.EQ.2) THEN
          WRITE(STDOUT,52)
        ELSE
          WRITE(STDOUT,53)
        ENDIF
      ENDIF
C     Compute the solution ratio to use in estimating the starting
C     point for the next upstream head in the sequence. 
      RATIO = YMID/YL 
      YR = YMID
      ZR = YR + ZBOTR
 
      FDROP = ZL - ZR
      QFREE = Q
 
      RETURN
      END
C
C
C
      SUBROUTINE   CHNTAB
     I                   (STDIN, STDOUT, STDTAB, GRV, MINQ,
     M                    TABDIR, EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute a 2D table for flow through prismatic channel reach,
C     that is, for command CHANRAT
 
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL GRV, MINQ
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     GRV    - value of acceleration due to gravity
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'epscom.cmn'
      INCLUDE 'chncom.cmn'
      INCLUDE 'ftable.cmn'
      INCLUDE 'tabupgrade.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, N, NFRAC, NHU, TAB, TABTYP, WFLAG, XSTAB,
     A        LOCATION_HU, LOCATION_PFD, N_GT, N_GT_TWICE,
     B        XSTAB_NUMBER, ftpup, verbose
      REAL BOTSLP, DROP, ELEV, FDROP, FDVEC(PMXNHU), HDATUM, HUOLD, HUP,
     A     HUVEC(PMXNHU), LENGTH, PFDVEC(PMXFRC), POWER, QFREE,
     B     QMAT(PMXNHU,PMXFRC), QRATIO, ZBOTL, YRATIO, ZBOTR,
     C     GLOBAL_ERROR, QD_SPAN, RMS_ERROR, LIPREC, MINPFD, FMID, XMID,
     D     RERR, LOCAL_MINQ, zrhufd, epsint_sngl

      REAL*8 SIESQR, easting, northing, dnull
      CHARACTER CHAR4*4, CHAR5*5, CHAR6*6, CHRVEC(5)*5, HEAD*80,
     A          LABEL*50, LINE*80, SAVEIT*4, CQ*8, TABID*16,
     B          XSTABID*16, ELEV_STRING*10,
     c          zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, FLOAT, MAX, MIN
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      CHARACTER GET_TABID*16
      EXTERNAL CHKCFC, CHKTAB, FRFCHN, inline, SBFCHN, TABCHK, TWDOUT,
     A         CHANRAT_OPT, LENSTR, READ_TABID_PLUS, READ_TABID,
     B         STRIP_L_BLANKS, FNDELV, GET_TABID, get_east_north
 
      data dnull/-33d6/
C     + + + INPUT FORMATS + + +
C 1    FORMAT(7X,I5,1X,A4)
 2    FORMAT(A5,1X,A)
 4    FORMAT(A6,1X,I5,2A5)
 6    FORMAT(A6,1X,F10.0,1X,8X,A10)
 8    FORMAT(A6,1X,F10.0)
 16   FORMAT(A5,1X,F10.0)
 24   FORMAT(A80)
 26   FORMAT(A5,1X,I5)
 27   FORMAT(A4,1X,I5,5A5)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' TABID= ',A,'  TYPE=',i5)
 51   FORMAT(/,' Upstream Head=',F9.4,'  Free Drop=',F9.4,
     A       ' Normal Flow=',F15.3)
 52   FORMAT(/,' ',A5,'=',A)
 53   FORMAT(1X,' PartialFD  Dns Head      Flow    QdrErr  RelError')
 54   FORMAT(/,' ','XSTAB= ',A)
 55   FORMAT(1X,F10.6,F10.4,2X,A8,F10.4)
 56   FORMAT(' ',A6,'=',F10.4)
 57   FORMAT(1X,F10.6,F10.4,2X,A8,F10.4,F10.3)
 58   FORMAT(' ',A6,'=',F10.1,' Middle elevation=',F10.4)
 60   FORMAT('   Free flow interpolation between head=',F9.4,
     A       ' and ',F9.4,/,'   has estimated maximum relative',
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
 68   format('; Flow defining minimum head=',f10.3)
 71   FORMAT(/,'  Table type 5 replaced by type 13.')
 72   FORMAT(/,' *ERR:628* Channel length <= 0.0. must be > 0.0.')
 74   FORMAT(/,' *ERR:629* ERRKND=',I5,' invalid. range 0 through 1.')
 76   FORMAT(/,' Adaptive Simpson''s rule for steady flow.')
 77   FORMAT(/,' Romberg rule deleted. Using Adaptive Simpson''s rule.')
 78   FORMAT(/,' Integration error tolerance=',F7.3)
 79   FORMAT(/,' Requested linear interpolation precision=',F5.3)
 80   FORMAT(/,' Minimum partial free drop=',F8.3)
 81   FORMAT(/,' Upstream Head=',F9.4,'  Free Drop=',F9.4)
 82   FORMAT(3X,' Ups WS Elevation=',F9.4)
83    format(/,' Minimum flow target=',F10.3 )
 87   FORMAT(/' *ERR:607* TABLE# <= 0')
 88   FORMAT(/,' Checking cross section table for possible critical',
     A       ' flow',/,5X,' problems.')
 89   FORMAT(/,' *WRN:551* CHANRAT command may not converge because',
     A       ' critical',/,10X,' flow decreases.')
 90   FORMAT('  Processing CHANRAT TabId= ',A)
 91   FORMAT(/,' Normal depth deviation-absolute=',F8.4,
     a       /,1X,'Normal depth deviation-relative=',F8.4)
 92   FORMAT(/,' ',A80)
 94   FORMAT(/,' ',A5,'=',I5)
 95   FORMAT(/,' ','Input complete. begin computations')
 96   FORMAT(/,' ','Datum for heads is:',F10.2)
 97   FORMAT(' ',A4,'=',I5,5A5)
 98   FORMAT(' *ERR:589* table type not 6 or 13.')
 99   FORMAT(/,' *ERR:630* Upstream heads non-increasing at:',F10.2)
C***********************************************************************
c      write(stdout,*) ' chntab: ftp=',ftp,' on entry.'

C     Clear the global value of error.
      GLOBAL_ERROR = 0.0

C     Clear the linear interpolation precision 
      LIPREC = 0.0

      call GET_chnrt_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_chnrt_ITEMS(
     M             EFLAG,
     O             TAB, tabtyp, errknd, inthow, epsint_sngl, nddabs, 
     o             nddrel, LOCAL_MINQ,
     O             zone, hgrid, vdatum, unitsys, basis, 
     O             easting, northing)

      epsint = dble(epsint_sngl)

c     Check local min flow.  Default is 0.0 to force use of global 
c     value.
      if(local_minq.eq.0.0) then
        local_minq = minq
      endif

      IF(TABTYP.EQ.5) THEN
        WRITE(STDOUT,71)
        TABTYP = 13
      ENDIF
      IF(TABTYP.NE.13.AND.TABTYP.NE.6) THEN
        WRITE(STDOUT,98)
        EFLAG = 1
      ENDIF
 
      IF(ERRKND.LT.0.OR.ERRKND.GT.1) THEN
        WRITE(STDOUT,74) ERRKND
        EFLAG = 1
      ENDIF
     
       
      if(epsint.le.0.0) then
        IF(GRV.GT.15.0) THEN
C         US standard unit 
          EPSINT = 0.1D0
        ELSE
C         SI or metric
          EPSINT = 0.030480D0
        ENDIF
      ENDIF

      if(nddabs.le.0.0) then
        IF(GRV.GT.15.0) THEN
          NDDABS = 0.005
        ELSE
          NDDABS = 0.001524
        ENDIF
      ENDIF
      

      TABID = GET_TABID(TAB)
     
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), tabtyp

C     MAKE SURE TABLE NUMBER IS NOT ALREADY USED IN THIS INPUT
 
      IF(FTPNT(TAB).NE.0) CALL TAB_IN_USE
     I                                   (TABID,
     M                                     EFLAG)

      IF(INTHOW.EQ.1) THEN
        WRITE(STDOUT,76)
      ELSE
        WRITE(STDOUT,77)
        inthow = 1
      ENDIF
      WRITE(STDOUT,78) EPSINT
      WRITE(STDOUT,91) NDDABS, NDDREL
      write(stdout,83) local_minq
      write(stdtab,68) local_minq
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) CHAR5, LABEL
      WRITE(STDOUT,52) CHAR5, LABEL
 
C     INPUT THE CROSS SECTION TABLE FOR THE CHANNEL
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'XSTAB',
     O                EFLAG, XSTABID, XSTAB)
      WRITE(STDOUT,54)  XSTABID(1:LENSTR(XSTABID))

C     Save the table number
      XSTAB_NUMBER = XSTAB 
      IF(XSTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              XSTAB,
     O              EFLAG)
        if(easting <= dnull) then
c         Get values from the cross-section table. 
          call get_east_north(
     i                        stdout, xstab,
     o                        easting, northing)
        endif
      ENDIF
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CHAR6, BOTSLP
 
      WRITE(STDOUT,56) CHAR6, BOTSLP
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CHAR6, LENGTH, ELEV_STRING
      CALL STRIP_L_BLANKS(
     M                    ELEV_STRING)
      IF(ELEV_STRING.EQ.'TAB'.OR.
     A   ELEV_STRING.EQ.'tab') THEN
C       Get the elevation from the function table.
        CALL FNDELV
     I             (XSTAB_NUMBER, STDOUT,
     O              EFLAG, ELEV)
      ELSE
        READ(ELEV_STRING,'(F10.0)',ERR=991) ELEV
      ENDIF

      WRITE(STDOUT,58) CHAR6, LENGTH, ELEV
 
      IF(LENGTH.LE.0.0)  THEN
        WRITE(STDOUT,72)
        EFLAG = 1
      ENDIF
 
C     INITIALIZE THE CONSTANT PART OF CHNCOM
 
      SBOT = BOTSLP
      L = LENGTH
      XSADR = XSTAB
 
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
      WRITE(STDOUT,94) CHAR5, NFRAC
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) CHAR5, POWER
      WRITE(STDOUT,66) CHAR5,POWER
 
 
C     INPUT THE HEAD SEQUENCE
 
      I = 1
      HUOLD = -1.0
 300  CONTINUE
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
C       Check for optional parameters for definition of 
C       table optimization. 
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
 
          ENDIF
        ENDIF
          GOTO 300
 310  CONTINUE


      IF(NFRAC.GT.PMXFRC-1) THEN
        WRITE(STDOUT,
     A'('' *ERR:549* MORE THAN '',I5,'' FRACTIONS OF FREE DROP'')')
     B       PMXFRC
 
        EFLAG = 1
      ENDIF
 
      IF(EFLAG.NE.0) RETURN

C     Adjust default value of EPSINT if table optimization
C     was selected.
      IF(LIPREC.GT.0.0) THEN
        IF(ABS(EPSINT - 0.1D0)/0.1D0.LE.1.D-6.OR.
     A     ABS(EPSINT - 0.03048D0)/0.03048D0.LE.1.D-6) THEN
C         Default value.  Therefore make smaller.
          EPSINT = EPSINT/2.D0
        ENDIF
      ENDIF
C     CHECK FOR THE CRITICAL FLOW AND CELERITY IN THE TABLE.
      WRITE(STDOUT,88)
      CALL CHKCFC
     I           (GRV, STDOUT, XSTAB,
     O            WFLAG)
      IF(WFLAG.NE.0) THEN
        WRITE(STDOUT,89)
      ENDIF
 
C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.
 
      WRITE(STDOUT,95)
 
C     SELECT THE POINT DEFINING HEAD FOR THE CHANNEL.
 
      ZBOTL = ELEV + 0.5*BOTSLP*LENGTH
      ZBOTR = ELEV - 0.5*BOTSLP*LENGTH
 
c     Drop to free flow at zero flow is equal to the drop in 
c     elevation on the defining plane but only if there is a drop
c     in elevation from ups to dns.  Otherwise it is zero.

      zrhufd = zbotl - zbotr
      if(zrhufd.lt.0.0) zrhufd = 0.0

      HDATUM = MAX(ZBOTL, ZBOTR)
 
      WRITE(STDOUT,96) HDATUM

      IF(LIPREC.GT.0.0) THEN
C       Recompute the head and partial free drop sequence
        CALL CHANRAT_OPT
     I                  (STDOUT, LIPREC, MINPFD, HDATUM, 
     I                   ZBOTL, ZBOTR, LOCAL_MINQ,
     M                   NHU, HUVEC, NFRAC, 
     O                   PFDVEC, EFLAG)
        IF(EFLAG.GT.0) THEN
          STOP 'Abnormal stop. Error(s) found.'
        ENDIF
      ELSE
C       COMPUTE THE PROPORTIONS OF FREE DROP
 
        DO 200 I=1,NFRAC
          PFDVEC(I) = (FLOAT(I-1)/FLOAT(NFRAC-1))**POWER
 
 200    CONTINUE

      ENDIF

C     Now compute the table using the head and PFD sequences. 
C     Set the depth ratio to initial value for assistance
C     in finding root.
      YRATIO = 0.5
C     Clear the sum of interpolation error squared and the quadrature
C     span used to compute the root-mean-squared error.
      SIESQR = 0.0D0
      QD_SPAN = 0.0

      N_GT = 0
      N_GT_TWICE = 0
      DO 2000 I=1,NHU
 
        HUP = HUVEC(I)
 
        IF(I.EQ.1) THEN
          CALL FRFCHN
     I               (STDOUT, HUP, HDATUM, ZBOTL, ZBOTR,
     M                YRATIO,
     O                EFLAG, QFREE, FDROP)
          QMAT(I,1) = 0.0
          QMAT(I,NFRAC) = QFREE
          FDVEC(I) = FDROP
        ELSE
C         Compute the flow between adjacent heads for error 
C         estimation.
          XMID = 0.5*(HUVEC(I-1) + HUP)
          CALL FRFCHN
     I               (STDOUT, XMID, HDATUM, ZBOTL, ZBOTR,
     M                YRATIO,
     O                EFLAG, FMID, FDROP)
          CALL FRFCHN
     I               (STDOUT, HUP, HDATUM, ZBOTL, ZBOTR,
     M                YRATIO,
     O                EFLAG, QFREE, FDROP)
          QMAT(I,1) = 0.0
          QMAT(I,NFRAC) = QFREE
          FDVEC(I) = FDROP

          IF(EFLAG.NE.0) THEN
            RETURN
          ENDIF
          RERR = (0.5*(QFREE + QMAT(I-1,NFRAC)) - FMID)/FMID
          QD_SPAN = QD_SPAN + HUP - HUVEC(I-1) 
          SIESQR = SIESQR + 0.6666667*RERR**2*(HUP - HUVEC(I-1))
          WRITE(STDOUT,60) HUVEC(I-1), HUP, RERR
          IF(ABS(RERR).GT.GLOBAL_ERROR) THEN
            GLOBAL_ERROR = ABS(RERR)
            LOCATION_HU = I-1
            LOCATION_PFD = NFRAC
          ENDIF
          IF(ABS(RERR).GT.LIPREC.AND.LIPREC.GT.0.0) THEN
            N_GT = N_GT + 1
            IF(ABS(RERR).GT.2.*LIPREC) THEN
              N_GT_TWICE = N_GT_TWICE + 1
            ENDIF
          ENDIF
        ENDIF
 
        IF(SBOT.GT.0.0) THEN
          WRITE(STDOUT,51) HUP, FDROP, QN
          WRITE(STDOUT,82) HDATUM + HUP
        ELSE
          WRITE(STDOUT,81) HUP, FDROP
          WRITE(STDOUT,82) HDATUM + HUP
        ENDIF
 
        WRITE(STDOUT,53)
 
        CALL VAR_DECIMAL(QFREE,
     O                   CQ)
        WRITE(STDOUT,55) PFDVEC(NFRAC), HUP - FDROP, CQ, QERR
 
C       FOR EACH OF THE PARTIAL FREE DROPS(EXCLUDING 0.00 AND 1.0)
C       COMPUTE THE FLOW THROUGH THE CHANNEL.

C       Initialize the ratio of flow rates used to provide an initial
C       estimate for the root. 
        QRATIO = 0.98        
        DO 500 J=NFRAC-1,2,-1

C         Compute a value midway between PFD's for error estimation.
          
          DROP = FDROP*0.5*(PFDVEC(J) + PFDVEC(J+1))
          YR = YL + ZBOTL - DROP - ZBOTR
          IF(SBOT.GT.0.0.AND.ABS(YL - YR).LT.
     A           MIN(NDDABS+EPSARG,(NDDREL+EPSARG)*YL)) THEN
C           USE NORMAL FLOW.
            FMID = QN
          ELSE
C            WRITE(STDOUT,*) 'Compute midway between PFD values'
            CALL SBFCHN
     I                 (STDOUT,
     M                  QRATIO, EFLAG)
            FMID = Q
          ENDIF

C         Compute the next tabulated value. 
          DROP = FDROP*PFDVEC(J)
          YR = YL + ZBOTL - DROP - ZBOTR
          IF(SBOT.GT.0.0.AND.ABS(YL - YR).LT.
     A           MIN(NDDABS+EPSARG,(NDDREL+EPSARG)*YL)) THEN
C           USE NORMAL FLOW.
            Q = QN
          ELSE
C           Compute submerged flow.  YL and Q in CHNCOM are at the
C           values they last had.  YL remains constant once set by
C           free flow computation.  Q will change as we submerge
C           the flow. Use the last value as the starting value.
C            WRITE(STDOUT,*) ' Compute at tabulated value.'
            CALL SBFCHN
     I                 (STDOUT,
     M                  QRATIO, EFLAG)
          ENDIF
          QMAT(I,J) = Q
          RERR = (0.5*(Q + QMAT(I,J+1)) - FMID)/FMID
          IF(ABS(RERR).GT.GLOBAL_ERROR) THEN
            GLOBAL_ERROR = ABS(RERR)
            LOCATION_HU = I
            LOCATION_PFD = J 
          ENDIF
          IF(ABS(RERR).GT.LIPREC.AND.LIPREC.GT.0.0) THEN
            N_GT = N_GT + 1
            IF(ABS(RERR).GT.2.*LIPREC) THEN
              N_GT_TWICE = N_GT_TWICE + 1
            ENDIF
          ENDIF
          SIESQR = SIESQR + 0.666667*RERR**2*FDROP*
     A                               (PFDVEC(J+1) - PFDVEC(J))
          QD_SPAN = QD_SPAN + FDROP*(PFDVEC(J+1) - PFDVEC(J))

          CALL VAR_DECIMAL(Q,
     O                     CQ)
          WRITE(STDOUT,57) PFDVEC(J), HUP - DROP, CQ, QERR,
     A                     RERR
 
 500    CONTINUE
 2000 CONTINUE

      IF(LOCATION_PFD.EQ.NFRAC) THEN

        WRITE(STDOUT,61) GLOBAL_ERROR, HUVEC(LOCATION_HU),
     A     HUVEC(LOCATION_HU+1)

      ELSE
        WRITE(STDOUT,63) GLOBAL_ERROR, HUVEC(LOCATION_HU),
     A     PFDVEC(LOCATION_PFD), PFDVEC(LOCATION_PFD+1)
      ENDIF

      RMS_ERROR = SQRT(SIESQR/QD_SPAN)
      WRITE(STDOUT,62) RMS_ERROR
      IF(LIPREC.GT.0.0) THEN
        N = NHU - 1 + NHU*(NFRAC - 1)
        WRITE(STDOUT,64) FLOAT(N_GT)/N, FLOAT(N_GT_TWICE)/N
      ENDIF

C     OUTPUT THE TABLE
      CALL TWDOUT
     I           (STDOUT, STDTAB, TAB, LABEL, NHU, NFRAC, HUVEC, FDVEC,
     I            PFDVEC, QMAT, HDATUM,
     I            TABTYP, ' CHANRAT', zrhufd,
     i            zone, hgrid, vdatum, unitsys, basis,
     i            easting, northing,
     O            EFLAG)

      if(twod_cubic_out.eq.'YES') then
        verbose = 1

        CALL twodfit
     I           (STDOUT,  TAB, NHU, NFRAC, HUVEC, FDVEC,
     I            PFDVEC, QMAT, HDATUM,
     I            TABTYP, ' CHANRAT', zrhufd, verbose,
     M            FTP,
     O            EFLAG, ftpup)
      endif


C      DO 890 I=1,NHU
C        WRITE(STDOUT,68) HUVEC(I), QMAT(I,NFRAC)
C68    FORMAT(F10.6,1PE15.6)
890   CONTINUE
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE CHANRAT_OPT
     I                      (STDOUT, LIPREC, MINPFD, HDATUM,
     I                       ZBOTL, ZBOTR, MINQ,
     M                       NHU, HUVEC, NFRAC, 
     O                       PFDVEC, EFLAG)

C     Compute a sequence for upstream head and partial free drops
C     to produce an optimized table for CHANRAT.

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NFRAC,  NHU, STDOUT
      REAL  HDATUM, LIPREC, MINPFD, ZBOTL, ZBOTR, MINQ,
     A     HUVEC(PMXNHU), PFDVEC(PMXFRC)

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     LIPREC - requested linear interpolation precision
C     MINPFD - minimum value of partial free drop in the 2-D table
C     HDATUM - head datum for the table
C     ZBOTL - bottom elev. at left end of channel
C     ZBOTR - bottom elev. at right end of channel
C     NHU - number of upstream heads.
C     HUVEC - upstream head sequence
C     NFRAC - number of partial free drops
C     PFDVEC - sequence of partial free drops
C     EFLAG - error flag
      

C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'chncom.cmn'
      INCLUDE 'ftable.cmn'

C     + + + LOCAL VARIABLES + + +

      INTEGER  I, J, TABLT, TABGT, N, NFRAC_OLD
      REAL ARGRAT, B, DFROW, DFCOL, H1, H2, HUP, QFREE, 
     A  FDROP, YRATIO, Q1, Q2,  QRATIO, DROP

      REAL*8  XVEC(PMXNHU), FVEC(PMXNHU), BREAK_POINTS(PMXNHU),
     A        MAX_RERR
      CHARACTER CHAR16*16

C     Called program units
      EXTERNAL GET_INTERNAL_TAB_NUMBER

C     + + + OUTPUT FORMATS + + +
 70   FORMAT(/,' CHANRAT will use ',I5,' upstream heads for the 2-D',
     A         ' table.',/,'  The maximum estimated relative ',
     B         'interpolation error is:',F10.3)
 71   FORMAT(/,' CHANRAT will use ',I5,' partial free drops for ',
     A         'the 2-D table.',/,'  The maximum estimated',
     B         ' relative interpolation error is:',F10.3)
 72   FORMAT(/,' Using upstream head=',F8.3,' to define partial',
     A         ' free drops.')
73    FORMAT(/,' Minimum head revised to: ',F10.3,
     A         ' for target flow=',F10.3)
 99   FORMAT(/,' *ERR:635* Interpolation precision tables missing.',
     A   ' Check for version',/,11X,'of file: TYPE5.TAB.')
C***********************************************************************
C     Define the power-function interpolation precision tables.
      CHAR16 = '10001'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, CHAR16,
     M                             EFLAG,
     O                             TABLT)
      TABLT = FTPNT(TABLT)
      CHAR16 = '10002'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, CHAR16,
     M                             EFLAG,
     O                             TABGT)

      TABGT = FTPNT(TABGT)

      IF(TABLT.LT.1.OR.TABGT.LT.1) THEN
        WRITE(STDOUT,99) 
        EFLAG = 1
        RETURN
      ENDIF

C     Define a better minimum head.  Compute free flow at .5 and 1.0
C     times the user-given minimum head value.  Find the power in
C     a simple power-function fit to these two points.    Then compute
C     the min-head value that would give a desired minimum flow.
C     The desired minimum flow is under user control but has a 
C     default value. 
      H2 = HUVEC(1)
      H1 = 0.5*H2
      YRATIO = 0.5
      CALL FRFCHN
     I           (STDOUT, H1, HDATUM, ZBOTL, ZBOTR,
     M            YRATIO,
     O            EFLAG, Q1, FDROP)
      CALL FRFCHN
     I           (STDOUT, H2, HDATUM, ZBOTL, ZBOTR,
     M            YRATIO,
     O            EFLAG, Q2, FDROP)

C     Compute the power of the simple power function that fits
C     these two points (and (0,0) as well).
      B = LOG(Q1/Q2)/LOG(H1/H2)
      H1 = H2*(MINQ/Q2)**(1.0/B)
      WRITE(STDOUT,73) H1, MINQ
      HUVEC(1) = H1      


C     Save entry value of NFRAC for lower limit of point set.
      NFRAC_OLD = NFRAC

C     Define the point set for computing the free flows to fit with
C     a cubic spline. Find the argument ratio, ARGRAT, for a power
C     of 3.5, close to the maximum for CHANRAT.

      CALL TDLK10
     I           (STDOUT, TABGT, 10, 3.5, LIPREC,
     O            ARGRAT, DFROW, DFCOL)

C     Use the maximum and minimum head  and ARGRAT to
C     compute the number of heads to use and compute a new head
C     sequence.
      N = INT(LOG(HUVEC(NHU)/HUVEC(1))/LOG(ARGRAT) + 1.0) + 1
      IF(N.LT.NFRAC_OLD) N = NFRAC_OLD
      CALL RATIOPNT
     I             (N, DBLE(HUVEC(1)), DBLE(HUVEC(NHU)), 
     O                    XVEC)
      NHU = N
C     Now compute the free flows for this sequence of heads.
      YRATIO = 0.5
      DO 110 I=1,NHU

        HUP = REAL(XVEC(I))

        CALL FRFCHN
     I             (STDOUT, HUP, HDATUM, ZBOTL, ZBOTR,
     M              YRATIO,
     O              EFLAG, QFREE, FDROP)

        IF(EFLAG.NE.0) THEN
          RETURN
        ENDIF
        FVEC(I) = DBLE(QFREE)
110   CONTINUE

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
C     Select a head near the middle of the vector and compute the 
C     free and submerged flows for that head to define the 
C     basis for finding good breakpoints for the partial free 
C     drops.  
      I = INT(REAL(2*NHU)/3. + 1.)
      IF(I.GT.NHU) I = NHU
      HUP = HUVEC(I)
      WRITE(STDOUT,72) HUP
      CALL FRFCHN
     I           (STDOUT, HUP, HDATUM, ZBOTL, ZBOTR,
     M            YRATIO,
     O            EFLAG, QFREE, FDROP)
      IF(EFLAG.NE.0) THEN
        RETURN
      ENDIF

      
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
      QRATIO = 0.98        
      DO 140 J=N-1,1,-1
        DROP = FDROP*REAL(XVEC(J))
        YR = YL + ZBOTL - DROP - ZBOTR
        IF(SBOT.GT.0.0.AND.ABS(YL - YR).LT.
     A         MIN(NDDABS+EPSARG,(NDDREL+EPSARG)*YL)) THEN
C         USE NORMAL FLOW.
          FVEC(J) = DBLE(QN)
        ELSE
C         COMPUTE SUBMERGED FLOW.  YL AND Q IN CHNCOM ARE AT THE
C         VALUES THEY LAST HAD.  YL REMAINS CONSTANT ONCE SET BY
C         FREE FLOW COMPUTATION.  Q WILL CHANGE AS WE SUBMERGE
C         THE FLOW. USE THE LAST VALUE AS THE STARTING VALUE.
 
          CALL SBFCHN
     I               (STDOUT,
     M                QRATIO, EFLAG)
          FVEC(J) = DBLE(Q)
        ENDIF
 
 140  CONTINUE
      FVEC(N) = QFREE

C     Try to find improved breakpoints. 

      CALL FINDBRK
     I            (STDOUT, N, XVEC, FVEC, DBLE(LIPREC),
     I             2, 2, 
     O             NFRAC, BREAK_POINTS, EFLAG, MAX_RERR)
      IF(EFLAG.NE.0) THEN
        RETURN
      ENDIF
      PFDVEC(1) = 0.0
      NFRAC = NFRAC + 1
      DO 150 J=2,NFRAC
        PFDVEC(J) = BREAK_POINTS(J-1)
150   CONTINUE
      WRITE(STDOUT,71) NFRAC, MAX_RERR

      RETURN
      END



C
C
C 
      SUBROUTINE SET_chnrt_ITEM_DEFAULTS()

C     Set the default values in the vectors used to
      IMPLICIT NONE
      INCLUDE 'chnrtitm.cmn'

C***********************************************************************
C     Default for: TABID 
      chnrtITMCTAB(  1) = '    '           
C     Default for: TABLE - note # is ignored in the standard scanner 
      chnrtITMCTAB(  2) = '    '           
C     Default for: TYPE
      chnrtITMITAB(  1) = 13
C     Default for: ERRKND-absolute error
      chnrtITMiTAB(  2) = 0
C     Default for: INTHOW - simpson adaptive-only option
      chnrtITMiTAB(  3) = 1
C     Default for: EPSINT - special value to signal need for default assignment later
      chnrtITMFTAB(  4) = -1.0
C     Default for: NDDABS - special value to signal need for default assignment later
      chnrtITMFTAB(  5) = -1.0
C     Default for: NDDREL 
      chnrtITMFTAB(  6) = 0.005
C     Default for: LOCAL_MINQ 
      chnrtITMFTAB(  7) = 0.0
C     Default for: ZONE
      chnrtITMCTAB(  3) = 'NONE'           
C     Default for: HGRID
      chnrtITMCTAB(  4) = 'NONE'           
C     Default for: VDATUM
      chnrtITMCTAB(  5) = 'NONE'           
C     Default for: UNITSYS
      chnrtITMCTAB(  6) = 'NONE' 
C     Default for: BASIS
      chnrtITMCTAB(  7) = 'NONE' 
c     Default for: EASTING
      chnrtITMDTAB(  5) = -33d6
c     Default for: NORTHING
      chnrtITMDTAB(  6) = -33d6
      

      RETURN
      END
C
C
C
      SUBROUTINE  SET_chnrt_ITEMS(
     M             EFLAG,
     O             TAB, TYPE, errknd, inthow, epsint, nddabs, nddrel,
     O             LOCAL_MINQ, zone, hgrid, vdatum, unitsys, basis,
     O             easting, northing)

C     Set items in CHANRAT
C     All values not set explicitly by user are at their default value.

      IMPLICIT NONE

      INTEGER TAB, TYPE, GETQ, GETY2, EFLAG, errknd, inthow
      REAL epsint, nddabs, nddrel, LOCAL_MINQ
      real*8 easting, northing
      character*8 zone, hgrid, vdatum, unitsys, basis
    

C     Local
      CHARACTER*16 KEY1, KEY2

      INCLUDE 'chnrtitm.cmn'
      INCLUDE 'stdun.cmn'
C***********************************************************************
C     Set the table id.  There  are two strings allowed as the variable name:
C     TABID or TABLE
      KEY1 = chnrtITMCTAB(  1)
      KEY2 = chnrtITMCTAB(  2)
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
      TYPE =   chnrtITMITAB( 1)
C     Set the value for errknd
      errknd = chnrtITMITAB( 2)
C     Set the value for inthow
      inthow =  chnrtITMITAB( 3)   
C     Set the value for epsint
      epsint = chnrtITMFTAB( 4)
C     Set the value for nddabs
      nddabs = chnrtITMFTAB( 5)
C     Set the value for nddrel
      nddrel = chnrtITMFTAB( 6)
C     Set the value for the local MINQ
      LOCAL_MINQ = chnrtITMFTAB(7)
C     Set the value for ZONE
      ZONE = chnrtITMCTAB( 3)
c     Set the value for HGRID
      HGRID = chnrtITMCTAB( 4)
c     Set the value for VDATUM
      VDATUM = chnrtITMCTAB( 5)
c     Set the value for UNITSYS
      UNITSYS = chnrtITMCTAB( 6)             
c     Set the value for BASIS
      basis = chnrtITMCTAB( 7)             
c     Set the value for EASTING
      EASTING =  chnrtITMDTAB( 5)
c     Set the value for NORTHING
      NORTHING = chnrtITMDTAB( 6)
    
      RETURN
      END

C
C
C
      SUBROUTINE GET_chnrt_ITEMS(STDIN, STDOUT,
     M                   EFLAG)

C     Get the table id and various options for CHANRAT command

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG

      INCLUDE 'chnrtitm.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  INTVAL, REAVAL, NONE,
     A         CHRVAL, DPRVAL, EXACT, LOWER,
     B         NUMERIC, CHAR, N_SYMBOL, NXTBLK
      PARAMETER(N_SYMBOL=17, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, NXTBLK=2, NONE=0,
     B          EXACT=0,LOWER=1, NUMERIC=0, CHAR=1)

      INTEGER MAX_LINE
     A        

      EXTERNAL GET_NAMED_ITEMS, SET_chnrt_ITEM_DEFAULTS

C     + + + SAVED VALUES + + +
      INTEGER GROUP(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL), GROUP_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE  SYMBOL_TABLE, GROUP, RESPONSE_TYPE, CONVERT_RULE,
     A      GROUP_INDEX

      DATA  SYMBOL_TABLE /
     *'TABID','TABLE','TYPE','ERRKND','INTHOW','EPSINT','NDDABS',
     A 'NDDREL','MINQ','LABEL','ZONE','HGRID','VDATUM','UNITSYS',
     b 'EASTING','NORTHING','BASIS'/
                                                                      
      DATA GROUP  /
     *CHAR, CHAR, 7*NUMERIC,NXTBLK,4*CHAR,2*NUMERIC,CHAR/          

      DATA GROUP_INDEX /
     *1, 2, 1, 2, 3, 4, 5, 6, 7, 0, 3, 4, 5, 6, 9, 11, 7/
                                       
      DATA RESPONSE_TYPE  /
     *CHRVAL,CHRVAL,3*INTVAL, 4*REAVAL, NONE, 4*CHRVAL, 2*DPRVAL,
     a  CHRVAL/
    
      DATA CONVERT_RULE /
     *LOWER,LOWER,3*EXACT,4*LOWER,EXACT,4*LOWER,2*LOWER,LOWER/
   
C***********************************************************************
C     Set Defaults
 
      CALL SET_chnrt_ITEM_DEFAULTS()

      MAX_LINE = 4
      CALL GET_NAMED_ITEMS(
     I                     STDIN, STDOUT,  MAX_LINE, N_SYMBOL,
     I  GROUP, RESPONSE_TYPE, CONVERT_RULE, GROUP_INDEX, SYMBOL_TABLE,
     I  MAXR_chnrtITM, MAXDP_chnrtITM, MAXC_chnrtITM, 'CHANRAT items',
     O  chnrtITMITAB, chnrtITMFTAB, chnrtITMDTAB, chnrtITMCTAB,
     O  EFLAG)
      
      RETURN

      END

