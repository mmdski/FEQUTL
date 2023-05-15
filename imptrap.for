C
C
C
      SUBROUTINE   TD13_FDROP
     I                   (STDOUT, IP, PTYPE, HU,
     I                    HBASE,
     O                    FDROP)
 
C     + + + PURPOSE + + +
C     Find drop to free flow for a given upstream head.
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER  IP, STDOUT, PTYPE
      REAL  HU, HBASE, FDROP
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT   - Fortran unit number for user output and messages
C     IP     - address of the table
C     PTYPE  - expected table type
C     HU     - upstream head
C     HBASE  - datum for heads
C     DROP   - drop to free flow
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER INC
      PARAMETER(INC=4)
 
C     + + + LOCAL VARIABLES + + +
      INTEGER LHU, LP, LPFD,  TAB, TYPE, L
      REAL DDROPU, HDROPL, HDROPR, HMAX, HUL, HUR
      CHARACTER TABID*16
 
C     External names
      INTEGER LENSTR
      CHARACTER GET_TABID*16      
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT('*WRN:04* HU > HMAX in two-D TABID=',A,' HU=',
     A       F10.2,' HMAX=',F10.2)
 54   FORMAT('*ERR:73* Unexpected type in two-D TABID=',A,
     A       ' Type=',I5,' expected Type=',I5)
C***********************************************************************
C     GET VALUES FROM HEADER-  MUST BE INITIALIZED ON INPUT SO THAT
C     POINTERS ARE SET TO THE PROPER VALUES AT START
 
      TAB = ITAB(IP+1)
      TABID = GET_TABID(TAB)
      L = LENSTR(TABID)
      TYPE = ITAB(IP+2)
      LHU = ITAB(IP+4)
      HMAX = FTAB(IP+9)

      IF(TYPE.NE.PTYPE) THEN
        WRITE(STDOUT,54) TABID(1:L), TYPE, PTYPE
        STOP 'Abnormal stop: errors found.' 
      ENDIF
 
      IF(HU.GT.HMAX) THEN
        WRITE(STDOUT,52) TABID(1:L), HU, HMAX
        HU = HMAX
      ENDIF
 
C     DETERMINE IF THERE IS FLOW AND IF IT IS KNOWN TO BE FREE
 
      IF(HU.LE.0.0) THEN
C       Return a zero for FDROP
        FDROP = 0.0
        RETURN
      ENDIF
 
C     UPSTREAM LEVEL IS ABOVE BASE.  FIND THE UPSTREAM LEVEL INTERVAL
C     CONTAINING HU
 
      IF(HU.GE.FTAB(LHU)) THEN
 100    CONTINUE
          IF(HU.LE.FTAB(LHU+INC)) GOTO 120
          LHU = LHU + INC
          GOTO 100
      ELSE
 110    CONTINUE
          LHU = LHU - INC
          IF(HU.GE.FTAB(LHU)) GOTO 120
          GOTO 110
      ENDIF
 120  CONTINUE
      ITAB(IP+4) = LHU

C      WRITE(STDOUT,*) ' LHU AFTER SEARCH=',LHU
 
C     FREE FLOW MAY RESULT IF THE DOWNSTREAM HEAD IS
C     SMALL ENOUGH
 
C     FIND THE HEAD DROP CORRESPONDING TO HU
 
      HUL = FTAB(LHU)
      HUR = FTAB(LHU+INC)
      HDROPL = FTAB(LHU+1)
      HDROPR = FTAB(LHU+INC+1)
      DDROPU = (HDROPR - HDROPL)/(HUR - HUL)
      FDROP = HDROPL + DDROPU*(HU - HUL)

 
      RETURN
 
      END
C
C
C
      SUBROUTINE LKT_ZB(
     I                  SET_INTERVAL, X,
     O                  ZB)

C     Find the bottom elevation at position X


      IMPLICIT NONE

      INTEGER SET_INTERVAL

      REAL*8 X, ZB

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER LOC
      REAL*8 FRAC

C     Called functions

      INTEGER FIND_BARREL_INTERVAL

      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL
C***********************************************************************
      IF(SET_INTERVAL.EQ.1) THEN

        LOC = FIND_BARREL_INTERVAL(
     I                             X)
        CALL  SET_BARREL_INTERVAL(
     I                            LOC)
      ENDIF

      FRAC = (X - X_L)/(X_R - X_L)
      ZB = Z_L + FRAC*(Z_R - Z_L)
      RETURN
      END
C
C
C
      SUBROUTINE LKT_AD(
     I                  SET_INTERVAL, X, Y,
     O                  AD, ALPHAD)

C     Find the area and alpha in the barrel at stations X and depth Y.


      IMPLICIT NONE

      INTEGER SET_INTERVAL

      REAL*8 X, Y, AD, ALPHAD

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER LOC
      REAL  AL, AR, ALPHAL, ALPHAR ,  TP,
     A      A, T, DT, J, K, DK, BETA, DBETA, ALPHA, DALPHA, QC
      REAL*8 FRAC

C     Called functions

      INTEGER FIND_BARREL_INTERVAL

      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL
C***********************************************************************
      IF(PRIS_FLAG.EQ.1) THEN
C       Prismatic channel.  Look up elements at left end only

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        AD = DBLE(A)
        ALPHAD = DBLE(ALPHA)
      ELSE
        IF(SET_INTERVAL.EQ.1) THEN

          LOC = FIND_BARREL_INTERVAL(
     I                               X)
          CALL  SET_BARREL_INTERVAL(
     I                              LOC)
        ENDIF

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              AL, T, DT, J, K, DK, BETA, DBETA, ALPHAL,
     O              DALPHA, QC)
        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_R,
     M              TP,
     O              AR, T, DT, J, K, DK, BETA, DBETA, ALPHAR,
     O              DALPHA, QC)
        FRAC = (X - X_L)/(X_R - X_L)
        AD = DBLE(AL) + FRAC*(DBLE(AR) - DBLE(AL))
        ALPHAD = DBLE(ALPHAL) + FRAC*(DBLE(ALPHAR) - DBLE(ALPHAL))
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE LKT_JDA(
     I                   SET_INTERVAL, X, Y,
     O                   AD, JD)

C     Find the first moment in the barrel at stations X and depth Y.

      IMPLICIT NONE

      INTEGER SET_INTERVAL

      REAL*8 X, Y, AD, JD

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER LOC
      REAL  AL, AR, JL, JR, TP,
     A     A, T, DT, J, K, DK, BETA, DBETA, ALPHA, DALPHA, QC

      REAL*8 FRAC
C     Called functions

      INTEGER FIND_BARREL_INTERVAL

      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL
C***********************************************************************
      IF(PRIS_FLAG.EQ.1) THEN
C       Prismatic channel.  Look up elements at left end only

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        AD = DBLE(A)
        JD = DBLE(J)

      ELSE
        IF(SET_INTERVAL.EQ.1) THEN

          LOC = FIND_BARREL_INTERVAL(
     I                               X)
          CALL  SET_BARREL_INTERVAL(
     I                              LOC)
        ENDIF

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              AL, T, DT, JL, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_R,
     M              TP,
     O              AR, T, DT, JR, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        FRAC = (X - X_L)/(X_R - X_L)
        JD = DBLE(JL) + FRAC*(DBLE(JR) - DBLE(JL))
        AD = DBLE(AL) + FRAC*(DBLE(AR) - DBLE(AL))
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE LKT_QCD(
     I                  SET_INTERVAL, X, Y,
     O                  QCD)

C     Find critical flow in the barrel at stations X and depth Y.


      IMPLICIT NONE

      INTEGER SET_INTERVAL

      REAL*8 X, Y, QCD

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER LOC
      REAL  QCL, QCR, TP,
     A     A, T, DT, J, K, DK, BETA, DBETA, ALPHA, DALPHA, QC
      REAL*8 FRAC

C     Called functions

      INTEGER FIND_BARREL_INTERVAL

      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL
C***********************************************************************
      IF(PRIS_FLAG.EQ.1) THEN
C       Prismatic channel.  Look up elements at left end only

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        QCD = DBLE(QC)
      ELSE
        IF(SET_INTERVAL.EQ.1) THEN

          LOC = FIND_BARREL_INTERVAL(
     I                               X)
          CALL  SET_BARREL_INTERVAL(
     I                              LOC)
        ENDIF

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QCL)
        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_R,
     M              TP,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QCR)
        FRAC = (X - X_L)/(X_R - X_L)
        QCD = DBLE(QCL) + FRAC*(DBLE(QCR) - DBLE(QCL))
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE LKT_KD(
     I                  SET_INTERVAL, X, Y,
     O                  KD)

C     Find the conveyance in the barrel at stations X and depth Y.


      IMPLICIT NONE

      INTEGER SET_INTERVAL

      REAL*8 X, Y, KD

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER LOC
      REAL  KL, KR, TP,
     A     A, T, DT, J, K, DK, BETA, DBETA, ALPHA, DALPHA, QC
      REAL*8 FRAC

C     Called functions

      INTEGER FIND_BARREL_INTERVAL

      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL
C***********************************************************************
      IF(PRIS_FLAG.EQ.1) THEN
C       Prismatic channel.  Look up elements at left end only

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              A, T, DT, J, K, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        KD = DBLE(K)
      ELSE
        IF(SET_INTERVAL.EQ.1) THEN

          LOC = FIND_BARREL_INTERVAL(
     I                               X)
          CALL  SET_BARREL_INTERVAL(
     I                              LOC)
        ENDIF

        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_L,
     M              TP,
     O              A, T, DT, J, KL, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        TP = SNGL(Y)
        CALL XLKT22
     I             (ADRS_R,
     M              TP,
     O              A, T, DT, J, KR, DK, BETA, DBETA, ALPHA,
     O              DALPHA, QC)
        FRAC = (X - X_L)/(X_R - X_L)
        KD = DBLE(KL) + FRAC*(DBLE(KR) - DBLE(KL))
      ENDIF
      RETURN
      END


C
C
C
      SUBROUTINE GET_YCYNYM
     I                    (STDOUT, X, EPSARG, EPSF, EPSABS, 
     M                     YCRIT, YNORM,
     O                     YMAX, RFLAG)
C     Find the values of critical depth, normal depth (if it exists),
C     and the maximum argument.  Sets the channel interval 
C     before calling FIND_YCYNYM.  On entry, YCRIT and YNORM should
C     be estimates of the critical depth and normal depth respectively.
C     If normal depth does not exist, the input value of YNORM is ignored
C     and a value of -1.0 is returned for the normal depth. 


      IMPLICIT NONE
      
      INTEGER STDOUT, RFLAG

      REAL*8 X, YCRIT, YNORM, YMAX, EPSARG, EPSF, EPSABS


      INCLUDE 'barrel.cmn'

C     External function

      INTEGER FIND_BARREL_INTERVAL

C     Called routines

      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL, 
     A         FIND_YCYNYM
C     Local

      INTEGER LOC

      
C***********************************************************************
        LOC = FIND_BARREL_INTERVAL(
     I                             X)
        CALL SET_BARREL_INTERVAL(
     I                           LOC)
        CALL FIND_YCYNYM
     I                  (STDOUT, X, EPSARG, EPSF, EPSABS,
     M                   YCRIT, YNORM,
     O                   YMAX, RFLAG)

      RETURN
      END
c
C
C
C
      SUBROUTINE   DFNDCD
     I                   (STDOUT, EPSARG, EPSF, EPSABS, YMAX,  
     M                    YE,
     O                    FLAG)
 
C     + + + PURPOSE + + +
C     Find critial depth for the given bottom slope. Full double
C     precsion version for UFGCULV.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG, STDOUT
      REAL*8  EPSARG, EPSF, EPSABS, YE, YMAX
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     YE     - Normal depth with initial estimate on entry 
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'barrel.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL*8 FL, FR, Y, YL, YR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL*8 YC_RESID
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL  GETTBN, YC_RESID, GET_TABID, LENSTR
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:743* TABID=',A,' Table overflow seeeking critical',
     A       ' depth for flow=',F10.2,11X,' D=',F10.2)
 52   FORMAT(' *ERR:626* TABID=',A,' Underflow seeking norm. depth for',
     A       ' flow=',F10.3)
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT FDBLRGF CLAIMS NONE IN',
     A       ' DFNDCD.')
 60   FORMAT(' *BUG:XXX* FDBLRGF ITERATION>100 IN DFNDCD')
C***********************************************************************
      FLAG = 0

C      WRITE(STDOUT,*) ' Entering DFNDCD with YE=',YE


C     SEARCH FOR A positive RESIDUAL
      YR = YE
      YL = -1.D0
 100  CONTINUE
        FR = YC_RESID(YR)
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
        IF(FR.GE.0.0D0) THEN
          GOTO 110
        ELSE
          YL = YR
          FL = FR
          YR = 0.7*YR + 0.3*VERT_D
          IF(ABS(YR - YMAX).LE.EPSARG) THEN
            TABID = GET_TABID(GETTBN(ADRS_L))
            WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), QD, VERT_D
            FLAG = 1
            GOTO 1000
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
C     POSITIVE RESIDUAL FOUND- SEARCH FOR NEGATIVE RESIDUAL
 
      IF(YL.LT.0.0D0) THEN
        YL = 0.6D0*YR
 120    CONTINUE
          FL = YC_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
          IF(FL.LE.0.0D0) THEN
            GOTO 130
          ELSE
            FR = FL
            YR = YL
            YL = 0.6D0*YL
            IF(ABS(YL).LT.EPSABS) THEN
              TABID = GET_TABID(GETTBN(ADRS_L))
              WRITE(STDOUT, 52) TABID(1:LENSTR(TABID)), QD
              FLAG = 1
              GOTO 1000
            ENDIF
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF
 
C     WE HAVE A SIGN CHANGE OR ONE OR BOTH POINTS HAVE ZERO RESIDUAL

C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG, EPSF, YC_RESID,
     M             YL, YR, FL, FR,
     O             Y, FLAG)
      YE = Y
      IF(FLAG.EQ.3) FLAG = 0

C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' Y =',Y,' FL=',FL
C      WRITE(STDOUT,*) ' EPSARG=',EPSARG,' EPSF=',EPSF
 
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
        GOTO 1000
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        GOTO 1000
      ENDIF
 
1000  CONTINUE 
      RETURN
      END
C
C
C
      SUBROUTINE   DFNDND
     I                   (STDOUT, EPSARG, EPSF, EPSABS, YMAX,  
     M                    YE,
     O                    FLAG)
 
C     + + + PURPOSE + + +
C     Find normal depth for the given bottom slope. Full double
C     precsion version for UFGCULV.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG, STDOUT
      REAL*8  EPSARG, EPSF, EPSABS, YE, YMAX
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     YE     - Normal depth with initial estimate on entry 
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'barrel.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL*8 FL, FR,  RES_AT_D, Y, YL, YR
      CHARACTER TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER GETTBN, LENSTR
      REAL*8 YN_RESID
      CHARACTER GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL  GETTBN, YN_RESID, GET_TABID, LENSTR
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:744* TABID=',A,' Soffit found seeking normal',
     A       ' depth for flow=',F10.2,11X,' D=',F10.2)
 52   FORMAT(' *ERR:626* TABID=',A,' Underflow seeking norm. depth for',
     A       ' flow=',F10.3)
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT FDBLRGF CLAIMS NONE IN',
     A       ' DFNDND.')
 60   FORMAT(' *BUG:XXX* FDBLRGF ITERATION>100 IN DFNDND')
C***********************************************************************
      FLAG = 0

C      WRITE(STDOUT,*) ' Entering DFNDND with YE=',YE

C     Compute the residual at the vertical diameter in the culvert. 
C     The residual function is defined such that if > 0 then
C     QD > the normal flow at the given depth. 
 
      RES_AT_D = YN_RESID(VERT_D)
C      WRITE(STDOUT,*) ' RES_AT_D=',RES_AT_D

      IF(RES_AT_D.GT.0.D0) THEN
C       The normal depth is above the culvert soffit. 
C       We assume that the capacity of the barrel 
C       is not increased when we are in the slot.  Thus
C       the normal depth should be made as large as 
C       possible and still remain within the range of the
C       cross-section function tables involved. 
 
        YE = YMAX
        GOTO 1000
      ENDIF
 
      IF(YE.GT.VERT_D) THEN
C       Adjust to be less than VERT_D
        YE = 0.5D0*VERT_D
      ENDIF

C     SEARCH FOR A negative RESIDUAL
      YR = YE
      YL = -1.D0
 100  CONTINUE
        FR = YN_RESID(YR)
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
        IF(FR.LE.0.0D0) THEN
          GOTO 110
        ELSE
          YL = YR
          FL = FR
          YR = 0.7*YR + 0.3*VERT_D
          IF(ABS(YR - VERT_D).LE.EPSARG) THEN
            TABID = GET_TABID(GETTBN(ADRS_L))
            WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), QD, VERT_D
            FLAG = 1
            GOTO 1000
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
C     Negative RESIDUAL FOUND- SEARCH FOR positive RESIDUAL
 
      IF(YL.LT.0.0D0) THEN
        YL = 0.6D0*YR
 120    CONTINUE
          FL = YN_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
          IF(FL.GE.0.0D0) THEN
            GOTO 130
          ELSE
            FR = FL
            YR = YL
            YL = 0.6D0*YL
            IF(ABS(YL).LT.EPSABS) THEN
              TABID = GET_TABID(GETTBN(ADRS_L))
              WRITE(STDOUT, 52) TABID(1:LENSTR(TABID)), QD
              FLAG = 1
              GOTO 1000
            ENDIF
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF
 
C     WE HAVE A SIGN CHANGE OR ONE OR BOTH POINTS HAVE ZERO RESIDUAL

C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG, EPSF, YN_RESID,
     M             YL, YR, FL, FR,
     O             Y, FLAG)
      YE = Y
      IF(FLAG.EQ.3) FLAG = 0

C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' Y =',Y,' FL=',FL
C      WRITE(STDOUT,*) ' EPSARG=',EPSARG,' EPSF=',EPSF
 
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
        GOTO 1000
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
        GOTO 1000
      ENDIF
 
1000  CONTINUE 
      RETURN
      END
C
C
C
      REAL*8 FUNCTION E(X,Y)

C     Compute the conservation quantity.

      IMPLICIT NONE
      REAL*8 X, Y

      INCLUDE 'barrel.cmn'
      INCLUDE 'grvcom.cmn'

      REAL*8 AD, ALPHAD
C***********************************************************************
      CALL LKT_AD(
     I            0, X, Y,
     O            AD, ALPHAD)

      E = ALPHAD*(QD/AD)**2/DBLE(GRAV2) + Y*COS_THETA

      RETURN
      END
C
C
C
      REAL*8 FUNCTION RHSE(X,Y)

C     Right-hand-side function.
      IMPLICIT NONE

      REAL*8 X, Y

      INCLUDE 'barrel.cmn'

C     Local variables
      
      REAL*8 KD
C***********************************************************************
      CALL LKT_KD(
     I            0, X, Y,
     O            KD)
     
      RHSE = SIN_THETA  - ABS(QD)*QD/KD**2
      RETURN
      END
C
C
C
      REAL*8 FUNCTION YC_RESID(Y)

C     Residual function for finding critical depth.

      IMPLICIT NONE
      REAL*8  Y

      INCLUDE 'barrel.cmn'

C     Local

      REAL*8 QCD
C***********************************************************************
      CALL LKT_QCD(
     I             0, XLOC, Y,
     O             QCD)

      YC_RESID = (COS_THETA - (QD/QCD)**2)
      RETURN
      END
C
C
C
      REAL*8 FUNCTION YN_RESID(Y)

C     Residual function for finding normal depth.

      IMPLICIT NONE
      REAL*8 Y

      INCLUDE 'barrel.cmn'

C     Local
      
      REAL*8 KD
C***********************************************************************
      CALL LKT_KD(
     I            0, XLOC, Y,
     O            KD)

      YN_RESID = 1.D0 - SIN_THETA*(KD/QD)**2

      RETURN
      END
C
C
C
      SUBROUTINE FIND_YCYNYM
     I                    (STDOUT, X, EPSARG, EPSF, EPSABS, 
     M                     YCRIT, YNORM,
     O                     YMAX, RFLAG)

C     Find the values of critical depth, normal depth (if it exists),
C     and the maximum argument.  Assumes the channel interval
C     has been set in barrel.cmn.

      IMPLICIT NONE
      
      INTEGER STDOUT, RFLAG
      REAL*8 X, YCRIT, YNORM, YMAX, EPSARG, EPSF, EPSABS


      INCLUDE 'barrel.cmn'

C     Local

      INTEGER FLAG

C***********************************************************************
      RFLAG = 1

C      WRITE(STDOUT,*) ' FIND_YCYNYM: YMAX_L=',YMAX_L
C      WRITE(STDOUT,*) ' YMAX_R=',YMAX_R

     
      YMAX = MIN(YMAX_L, YMAX_R)

C     Set the location for the residual functions
      XLOC = X

C     Compute critical depth.   YCRIT is an estimate.

      CALL DFNDCD
     I           (STDOUT, EPSARG, EPSF, EPSABS, YMAX,  
     M            YCRIT,
     O            FLAG)
      IF(FLAG.NE.0) THEN
        RFLAG = -6
C        RETURN
      ENDIF

C      WRITE(STDOUT,*) ' SIN_THETA=',SIN_THETA
C     Compute normal depth.
      IF(SIN_THETA.GT.0.D0) THEN
        CALL DFNDND
     I             (STDOUT, EPSARG, EPSF, EPSABS, YMAX,
     M              YNORM,
     O              FLAG)
        IF(FLAG.NE.0) THEN
          RFLAG = -7
          RETURN
        ENDIF
      ELSE
        YNORM = -1.D0
      ENDIF
        
      RETURN
      END
C
C
C
      SUBROUTINE SETUP_FOR_SECANT(
     I                            Y, YUP, YDN,
     O                            YA, YB)

C     Set an increment for starting the secant method so that the two 
C     points are close and both within the solution limits.  The solution
C     limits can be close together and we do not want to fall outside
C     those limits. 

      IMPLICIT NONE

      REAL*8 Y, YUP, YDN, YA, YB

C     Local variables
      
      REAL*8 DYTRY, DYLIM, YTRY
C***********************************************************************
      DYTRY = 0.005D0*Y
      DYLIM = 0.5D0*(YUP - YDN)
      IF(DYTRY.GT.DYLIM) DYTRY = DYLIM
      YTRY = Y + DYTRY
      IF(YTRY.GT.YUP) YTRY = Y - DYTRY
      YA = Y
      YB = YTRY
      RETURN
      END
C
C
C
      SUBROUTINE SET_LIMITS_AND_ESTIMATE(
     I                              STATE, SUB, SUP, YCRIT, 
     I                              YNORM, YMAX,
     M                              Y,
     O                              YUP, YDN, YCRIT_YES)

C     Set the upper and lower limits for a profile based on the 
C     state of the solution and the location of the current 
C     result.  Also checks the current estimate, Y, for the root,
C     and adjusts to be valid.

      IMPLICIT NONE
      INTEGER STATE, SUB, SUP, YCRIT_YES
      REAL*8 Y, YCRIT, YNORM, YMAX, YUP, YDN

C***********************************************************************
C     YCRIT_YES is a flag that signals that the computations can 
C     encounter critical flow as they progress. 

C     Clear the flag that signals that the computations can approach 
C     critical depth.
      YCRIT_YES = 0 
      IF(YNORM.LT.0.D0) THEN
C       Normal depth does not exist.  
        IF(STATE.EQ.SUB) THEN
          YUP = YMAX
          YDN = YCRIT
          IF(Y.LT.YCRIT) THEN
C           Move to the closest boundary.
            Y = YCRIT
          ELSEIF(Y.GT.YMAX) THEN
            Y = YMAX
          ENDIF
        ELSE
          YUP = YCRIT
          YDN = 0.D0
          YCRIT_YES = 1
          IF(Y.GT.YCRIT) THEN
            Y = YCRIT
          ENDIF
        ENDIF
      ELSE
C       Normal depth exists. Three cases: local slope is mild, steep, or
C       critical.
        IF(YCRIT.LT.YNORM) THEN
C         Local channel slope is mild.
          IF(STATE.EQ.SUB) THEN
            IF(Y.GE.YNORM) THEN
              YUP = YMAX
              YDN = YNORM
            ELSEIF(Y.GE.YCRIT) THEN
              YUP = YNORM
              YDN = YCRIT
            ELSE
C             Invalid solution: Expected subcritical result on a mild
C             slope but found a super critical estimate.  Move estimate
C             to closest boundary.
              Y = YCRIT
            ENDIF
          ELSE
C           Solution is super critical here.
            IF(Y.LE.YCRIT) THEN
              YUP = YCRIT
              YDN = 0.D0
              YCRIT_YES = 1
            ELSE
C             Invalid solution: Expected super critical result on a
C             mild slope but found a subcritcal estimate. Move 
C             estimate to closest boundary.
              Y = YCRIT
            ENDIF
          ENDIF
        ELSEIF(YCRIT.GT.YNORM) THEN
C         Local channel slope is steep.
          IF(STATE.EQ.SUB) THEN
            IF(Y.GE.YCRIT) THEN
              YUP = YMAX
              YDN = YCRIT
              YCRIT_YES = 1
            ELSE
C             Invalid solution: Expected subcritical result on a steep
C             slope but found a super critical estimate.  Move to 
C             closest boundary.
              Y = YCRIT
            ENDIF
          ELSE
C           Solution is super critical here.
            IF(Y.GT.YCRIT) THEN
C             Invalid solution: Expected a super critical result on a
C             steep slope but found a subcritical estimate.  Move to
C             closet boundary.
              Y = YCRIT
            ELSEIF(Y.GE.YNORM) THEN
              YUP = YCRIT
              YDN = YNORM
            ELSE
              YUP = YNORM
              YDN = 0.D0
            ENDIF
          ENDIF
        ELSE
C         Local channel slope is critical.
          IF(STATE.EQ.SUB) THEN
            IF(Y.GE.YCRIT) THEN
              YUP = YMAX
              YDN = YCRIT
            ELSE
C             Invalid solution: Expected a subcritical result on a
C             critical slope but found a super critical estimate.
C             Move to closest boundary.
              Y = YCRIT              
            ENDIF
          ELSE
C           Solution super critical here.
            IF(Y.GT.YCRIT) THEN
C             Invalid solution: Expected a super critical result on a 
C             critical slope but found a subcritical result. Move to
C             closest boundary.
              Y = YCRIT
            ELSE
              YUP = YCRIT
              YDN = 0.D0
            ENDIF
          ENDIF
        ENDIF
      ENDIF    
      RETURN
      END

C
C
C
      SUBROUTINE GET_ROOT
     I                   (H, F1, P1, X2, Y2, F, P, EPSF, EPSARG, EPSABS, 
     I                    YMAX, YMIN,
     M                    ROOT, 
     O                    RFLAG)

C     Find a root for subroutine IMPTRAP.  Use secant method as 
C     implemented in subroutine SECANT. 

      IMPLICIT NONE

      INTEGER RFLAG
      REAL*8 H, F1, P1, X2, Y2, EPSF, EPSARG, EPSABS, YMAX, YMIN,
     A       ROOT

      REAL*8 F, P
      EXTERNAL F, P


C     Local variables

      INTEGER I, MAXIT
      REAL*8 YL, YR, DY, OLDVAL, FL, FR, DF, DYTEMP

      DATA MAXIT/30/
     
C***********************************************************************
      RFLAG = 1
C     Set the initial starting points
      YL = Y2
      YR = ROOT
      OLDVAL = 1.D0
C      WRITE(STD6,*) ' ENTERING GET_ROOT. YL=',YL,' YR=',YR
 
      FL = (P(X2,YL) - P1 - 0.5D0*H*(F1 + F(X2,YL)))/P(X2,YL)
      IF(ABS(FL).LE.EPSF) THEN
        ROOT = YL
        RETURN
      ENDIF
C      WRITE(STD6,*) ' FL=',FL
      DY = YR - YL
      DO 100 I=1,MAXIT
        FR = (P(X2,YR) - P1 - 0.5D0*H*(F1 + F(X2,YR)))/P(X2,YR)
C       WRITE(STD6,*) ' I=',I,' FR=',FR
        IF(ABS(FR).LE.EPSF) THEN
C          WRITE(6,*) ' Convergence on FR with EPSF'
          ROOT = YR
          RETURN
        ENDIF
        DF = FL - FR
C        WRITE(STD6,50) I, FL, FR, DY
C50    FORMAT(' I=1',I5,' FL=',F10.5,' FR=',F10.5,' DY=',F15.8)
        IF(DF.NE.0.0) THEN
          DYTEMP = 0.0D0
          OLDVAL = DY/DF
          DY = FR*OLDVAL
          YL = YR
          FL = FR
          YR = YR + DY
          IF(YR.GT.YMAX) THEN
            YR = 0.5*(YL + YMAX)
            DYTEMP = YR - YL
          ELSEIF(YR.LT.YMIN) THEN
            YR = 0.5*(YL + YMIN)
            DYTEMP = YR - YL
          ENDIF
C          WRITE(STD6,*) ' DY=',DY
          IF(ABS(DY/YR).LE.EPSARG.OR.ABS(DY).LE.EPSABS) THEN
C            WRITE(6,*) ' Convergence on EPSARG or EPSABS'
            ROOT = YR
            RETURN
          ENDIF
          IF(DYTEMP.NE.0.0) THEN
            DY = DYTEMP
          ENDIF
        ELSE
C         DF = 0. SIGNAL POSSIBLE ERROR
          ROOT = YR
          IF(ABS(FR).LE.EPSF) THEN
C            WRITE(6,*) ' Convergence on EPSF when DF=0.0'
            RETURN
          ELSE
            RFLAG = 0
            RETURN
          ENDIF
        ENDIF
 
 100  CONTINUE
 
C     DROP THROUGH INDICATES NO CONVERGENCE WITHIN MAXIT ITERATIONS
      RFLAG = 0
      ROOT = YR
      RETURN
      END
C
C
C
      SUBROUTINE  INVERT_P(
     I                     PVALUE, X, Y, P, YMAX, YMIN, EPSF, 
     I                     EPSABS, EPSARG, 
     M                     ROOT,
     O                     RFLAG)

C     Invert P for IMPTRAP.  Use secant method based on SECANT.


      IMPLICIT NONE
      INTEGER RFLAG
      REAL*8  PVALUE, X, Y, EPSF, EPSARG, EPSABS, ROOT,
     A        YMAX, YMIN

      REAL*8  P
      EXTERNAL  P


C     Local variables

      INTEGER I, MAXIT
      REAL*8 YL, YR, DY, OLDVAL, FL, FR, DF, DYTEMP

      DATA MAXIT/30/
     
C***********************************************************************
      RFLAG = 1
C     Set the initial starting points
      YL = Y
      YR = ROOT
      OLDVAL = 1.D0
C      WRITE(STD6,*) ' ENTERING INVERT_P. YL=',YL,' YR=',YR
 
      FL = (PVALUE - P(X,YL))/PVALUE

C      WRITE(STD6,*) ' FL=',FL
      DY = YR - YL
      DO 100 I=1,MAXIT
        FR = (PVALUE - P(X,YR))/PVALUE
C       WRITE(STD6,*) ' I=',I,' FR=',FR
        IF(ABS(FR).LE.EPSF) THEN
          ROOT = YR
          RETURN
        ENDIF
        DF = FL - FR
C        WRITE(STD6,*) ' INVERT_P:I=',I,' FL=',FL,' FR=',FR,' DY=',DY
        IF(DF.NE.0.0) THEN
          DYTEMP = 0.0D0
          OLDVAL = DY/DF
          DY = FR*OLDVAL
          YL = YR
          FL = FR
          YR = YR + DY
          IF(YR.GT.YMAX) THEN
            YR = 0.5*(YL + YMAX)
            DYTEMP = YR - YL
          ELSEIF(YR.LT.YMIN) THEN
            YR = 0.5*(YL + YMIN)
            DYTEMP = YR - YL
          ENDIF
C          WRITE(STD6,*) ' DY=',DY
          IF(ABS(DY/YR).LE.EPSARG.OR.ABS(DY).LE.EPSABS) THEN
            ROOT = YR
            RETURN
          ENDIF
          IF(DYTEMP.NE.0.0) THEN
            DY = DYTEMP
          ENDIF
        ELSE
C         DF = 0. SIGNAL POSSIBLE ERROR
          ROOT = YR
          IF(ABS(FR).LE.EPSF) THEN
            RETURN
          ELSE
            RFLAG = 0
            RETURN
          ENDIF
        ENDIF
 
 100  CONTINUE
 
C     DROP THROUGH INDICATES NO CONVERGENCE WITHIN MAXIT ITERATIONS
      RFLAG = 0
      ROOT = YR
      RETURN
      END
C
C
C
      SUBROUTINE   ERRORS_IN_IMPTRAP(
     I                               STDOUT, RFLAG, N , X, NMAX, H, H2)

C     Report errors in subroutine IMPTRAP

      IMPLICIT NONE

      INTEGER STDOUT, RFLAG, N, NMAX

      REAL*8 X, H, H2


C     ****************************FORMATS*******************************
50    FORMAT(/,' *ERR:745*: Computational problem(s) in IMPTRAP on',
     A /, 11X, ' step after station=',F12.4,' and after index=',I10,'.',
     B /, 11X,  ' The problem was identified as:')
51    FORMAT(' Maximum number of steps=',I5,' exceeded.')
52    FORMAT(' Minimum step length reached seeking solution at full'  ,
     a      ' step=',1PD10.3)
53    FORMAT(' Minimum step length reached seeking solution at first',
     A       ' half-step=',1PD10.3)
54    FORMAT(' Minimum step length reached seeking solution at second',
     A       ' half-step=',1PD10.3)
55    FORMAT(' Minimum step length reached seeking solution after',
     A       ' extrapolation.')
56    FORMAT(' Critical depth not found.')
57    FORMAT(' Normal depth not found when it should exist.')
58    FORMAT(' Starting depth exactly matches critical depth. Desired ',
     A       'profile unknown.')
59    FORMAT(' Minimum step length reached seeking to meet error',
     A       ' tolerance.')
C***********************************************************************

      WRITE(STDOUT,50) X, N
      IF(RFLAG.EQ.-1) THEN
        WRITE(STDOUT,51) NMAX
      ELSEIF(RFLAG.EQ.-2) THEN
        WRITE(STDOUT,52) H
      ELSEIF(RFLAG.EQ.-3) THEN
        WRITE(STDOUT,53) H2
      ELSEIF(RFLAG.EQ.-4) THEN
        WRITE(STDOUT,54) H2
      ELSEIF(RFLAG.EQ.-5) THEN
        WRITE(STDOUT,55)
      ELSEIF(RFLAG.EQ.-6) THEN
        WRITE(STDOUT,56) 
      ELSEIF(RFLAG.EQ.-7) THEN
        WRITE(STDOUT,57)
      ELSEIF(RFLAG.EQ.-8) THEN
        WRITE(STDOUT,58)
      ELSEIF(RFLAG.EQ.-9) THEN
        WRITE(STDOUT,59)
      ENDIF
      RETURN
      END 
C
C
C
      SUBROUTINE IMPTRAP                  
     I                  (STDOUT, PRISMATIC, XS, XE, YS, EPS, EPSARG, 
     I                   EPSF, EPSABS, HMINA, HMAXA, YCRITA, YNORMA, 
     I                   YMAXA, P, F, FIND_YCYNYM, EXTRAP, YSFAC, 
     I                   DHLIM, NMAX,
     M                   N,
     O                   XVEC, YVEC, YCVEC, YNVEC, RFLAG)

C     Solve for the steady 1-D water surface profile using the 
C     conservation  form of the equations,

C                dP(x,y)
C                -------  = F(x,y)
C                dx

C     using the implicit trapezoidal rule, varying the step length to
C     keep the estimated local truncation error, in relative terms,
C     less than EPS.  If EXTRAP = 1, apply one Richardson extrapolation
C     correction to the result before starting the next step. 

C     Assumes that the proper channel interval has been set in the
C     common blocks internal to the functions P, F, and FIND_YCYNYM.

      IMPLICIT NONE

      INTEGER EXTRAP, N, NMAX, PRISMATIC, RFLAG, STDOUT
      REAL*8 XS, XE, YS, EPS, EPSARG, EPSF, EPSABS, HMINA, HMAXA, 
     A       YCRITA, YNORMA, YMAXA, YSFAC, DHLIM
      REAL*8 XVEC(NMAX), YVEC(NMAX), YCVEC(NMAX), YNVEC(NMAX)

      REAL*8 F, P
      EXTERNAL F, P, FIND_YCYNYM

C         Definition of dummy arguments.

C     STDOUT-  unit for output of messages.

C     PRISMATIC- if 1, then channel is primatic, if 0 non-prismatic.

C     XS-     starting value of x, distance along the channel.
C     XE-     ending value of x
C     YS-     starting value of y, the maximum depth in a cross section
C             The starting value should not exactly equal critical depth
C             to give a clear indication of which of two solutions is being
C             sought.  The rate of change of depth with distance at critical
C             depth is infinite.  Thus the distance between critical depth
C             and the depth that differs by 1 percent is negligible in any
C             engineering application. 
C     EPS-    limit for the maximum local truncation error in estimates
C             of E relative to E.

C     The following three values relate to the non-linear solution required
C     in the implicit trapezoidal method:
C     EPSARG- relative tolerance on the change in an argument when finding
C             a root. 
C     EPSF-   toleracne on the value of the function when finding a root.
C     EPSABS- absolute tolerance on the change in an argument when 
C             finding a root.

C     HMINA-   minimum step size permitted.
C     HMAXA-   upper limit for step size.

C     The following three values apply at XS only if the channel is
C     non-prismatic, and at all points if the channel is prismatic. 
C     YCRITA-  critical depth
C     YNORMA-  normal depth if it exists; otherwise, -1.0 
C     YMAXA-   maximum argument supported in the cross-section function
C              table defining the channel shape at XS.
C     P-      function of x and y, either the specific energy or the 
C             impulse-momentum function. 
C     F-      function of x and y, giving the correct right-hand side for
C             conservation law being used.
C     FIND_YCYNYM- subroutine to return values of critical depth, normal
C              depth, and maximum allowed depth.
C     EXTRAP- if 1 apply an extrapolation correction, otherwise, no
C             correction is applied.
C     YSFAC-  factor on YS giving the initial step.
C     DHLIM-  maximum increase ratio for a step.  DHLIM = 2 limits increases
C             to twice the previous value.
C     NMAX-   maximum number of steps permitted.
C     N-      on entry gives the index into the result vectors at which
C             to start storing values.  On exit gives the last index 
C             filed.  We will store values  such that the X location 
C             will always increase with increase in index.  
C             Also the initial point is only stored on the first 
C             entry to the IMPTRAP in a series of entries.
C     XVEC(*)-values of x at which a solution has been computed.
C     RFLAG-  result flag: 1-normal completion.  See ERRORS_IN_IMPTRAP.


C     External routines

      EXTERNAL SET_LIMITS_AND_ESTIMATE, SETUP_FOR_SECANT, GET_ROOT
C     Local variables

      INTEGER I, STATE, SUB, SUP, YCRIT_YES, INC
      PARAMETER (SUB=0, SUP=2)

      REAL*8 H, H2, HMAX, HMIN, X1, X2, XMID, F1, FMID, P1, P2, P2END,
     A       Y1, Y2, YT, YMID, Y2END, YCRIT2, YNORM2, YMAX2, 
     B       YUP2, YDN2, YCRIT_MID, YNORM_MID, YMAX_MID, YUP_MID, 
     C       YDN_MID, ROOT, PMID, RELERR, Y2BEST, P2BEST,
     D       HFAC

      

C***********************************************************************
      RFLAG = 1
      I = N

C     Set the initial step length and set the sign of the maximum step.
      IF(XE.GT.XS) THEN
C        H = YS*YSFAC
        H = HMAXA
        HMAX = HMAXA
        HMIN = HMINA
        INC = 1
      ELSE
C        H = -YS*YSFAC
        H = -HMAXA
        HMAX = -HMAXA
        HMIN = -HMINA
        INC = -1
      ENDIF

      IF((I.EQ.1.AND.INC.EQ.1).OR.
     A   (I.EQ.NMAX.AND.INC.EQ.-1)) THEN
C       Store the starting values
        XVEC(I) = XS
        YVEC(I) = YS
      ENDIF

C      WRITE(STDOUT,*) ' Initial H=',H,' YS=',YS

C     Initialize the local values
      X1 = XS
      Y1 = YS
      YCRIT2 = YCRITA
      YNORM2 = YNORMA
      YMAX2 = YMAXA
      YCVEC(I) = YCRIT2
      YNVEC(I) = YNORM2

C      WRITE(STDOUT,*) ' YCRIT2=',YCRIT2,' YNORM2=',YNORM2,
C     A                 ' YMAX2=',YMAX2
      YCRIT_MID = YCRIT2
      YNORM_MID = YNORM2
      YMAX_MID = YMAX2

C     Establish the state for the solution
      IF(YS.LT.YCRIT2) THEN
        STATE = SUP
      ELSEIF(YS.GT.YCRIT2) THEN
        STATE = SUB
      ELSE
C       Error-not known which profile is desired.
        RFLAG = -8
        GOTO 1000
      ENDIF

C      WRITE(STDOUT,*) ' STATE=',STATE

C     Establish the upper and lower bounds for the solution.
C     Applies at each point for a prismatic channel but will
C     be over ridden at each point for a non-prismatic channel.
        
      Y2 = Y1 
      CALL SET_LIMITS_AND_ESTIMATE(
     I                              STATE, SUB, SUP, YCRIT2, 
     I                              YNORM2, YMAX2,
     M                              Y2,
     O                              YUP2, YDN2, YCRIT_YES)

      YUP_MID = YUP2
      YDN_MID = YDN2

C      WRITE(STDOUT,*) ' YUP2=',YUP2,' YDN2=',YDN2

C     Start loop over steps until XE is reached, critical depth
C     is approached, minimum step is reached, or non-linear 
C     equation solution convergence failure occurs.

100   CONTINUE
        F1 = F(X1,Y1)
        P1 = P(X1,Y1)

C       Start loop within the step until the local error tolerance
C       is met or some computational failure occurs.
200     CONTINUE

          X2 = X1 + H

C         Check for the last step.
          IF(H.GT.0) THEN
            IF(X2.GE.XE) THEN
              X2 = XE
              H = X2 - X1
            ENDIF
          ELSE
            IF(X2.LE.XE) THEN
              X2 = XE
              H = X2 - X1
            ENDIF
          ENDIF

C         If the channel is non-prismatic, we must compute
C         new values of critical and normal depths.
          IF(PRISMATIC.EQ.0) THEN
          CALL FIND_YCYNYM
     I                    (STDOUT, X2, EPSARG, EPSF, EPSABS, 
     M                     YCRIT2, YNORM2,
     O                     YMAX2, RFLAG)
            IF(RFLAG.NE.1) THEN
C             Computational failure 
              GOTO 1000
            ENDIF
            IF(RFLAG.NE.1) THEN
C             Computational failure when normal depth is known to 
C             exist.
              GOTO 1000
            ENDIF
C           Compute the new profile limits and adjust the estimate.
            CALL SET_LIMITS_AND_ESTIMATE(
     I                                    STATE, SUB, SUP, YCRIT2, 
     I                                    YNORM2, YMAX2,
     M                                    Y2,
     O                                    YUP2, YDN2, YCRIT_YES)
          ENDIF
            
C         Make an estimate of the root to start the root finding 
C         process. Y2 is a good estimate (we hope!) but we need to set a 
C         nearby point that is valid for starting the secant
C         method for finding the root. Any new root must 
C         be between YUP and YDN.  At this point we know that
C         Y2 is in the proper range.  The upper and lower limits
C         could be close therefore we must choose the increment
C         to stay within the solution interval. 
          YT = Y2
          CALL SETUP_FOR_SECANT(
     I                          YT, YUP2, YDN2,
     O                          Y2, ROOT)


C         Find the solution at X2 with a step of H.
          CALL GET_ROOT
     I                 (H, F1, P1, X2, Y2, F, P, EPSF, EPSARG, EPSABS, 
     I                  YUP2, YDN2,
     M                  ROOT, 
     O                  RFLAG)
          Y2 = ROOT
C          WRITE(STDOUT,*) ' Y2=',Y2,' RFLAG=',RFLAG

          IF(RFLAG.EQ.0) THEN
C           No solution found or invalid solution found.  
C           Reduce the step and try again.
            H = 0.9D0*H
            IF(ABS(H).GT.ABS(HMIN)) THEN
              GOTO 200
            ELSE
              RFLAG = -2
              GOTO 1000
            ENDIF
          ENDIF      
          
          P2 = P(X2,Y2)

C         Now half the step and take two steps to X2.
          
          H2 = 0.5D0*H
          XMID = X1 + H2

C         Make estimate of the root at the new location. 

          YMID = 0.5D0*(Y1 + Y2)

          IF(PRISMATIC.EQ.0) THEN
            YCRIT_MID = YCRIT2
            YNORM_MID = YNORM2
            CALL FIND_YCYNYM
     I                      (STDOUT, XMID, EPSARG, EPSF, EPSABS, 
     M                       YCRIT_MID, YNORM_MID,
     O                       YMAX_MID, RFLAG)
            IF(RFLAG.NE.1) THEN
C             Computational failure 
              GOTO 1000
            ENDIF
            IF(RFLAG.NE.1) THEN
C             Computational failure when normal depth is known to 
C             exist.
              GOTO 1000
            ENDIF
C           Compute the new profile limits and adjust the estimate.
            CALL SET_LIMITS_AND_ESTIMATE(
     I                                    STATE, SUB, SUP, YCRIT_MID, 
     I                                    YNORM_MID, YMAX_MID,
     M                                    YMID,
     O                                    YUP_MID, YDN_MID, YCRIT_YES)
          ENDIF

          YT = YMID
          CALL SETUP_FOR_SECANT(
     I                          YT, YUP_MID, YDN_MID,
     O                          YMID, ROOT)

C         Find solution at XMID with step of H/2
          CALL GET_ROOT
     I                 (H2, F1, P1, XMID, YMID, F, P, EPSF, 
     I                  EPSARG, EPSABS,  YUP_MID, YDN_MID,
     M                  ROOT, 
     O                  RFLAG)
          YMID = ROOT
C          WRITE(STDOUT,*) ' YMID=',YMID,' RFLAG=',RFLAG

          IF(RFLAG.EQ.0) THEN
C           No solution found.  Reduce the step and try again.
            H = 0.9D0*H
            IF(ABS(H).GT.ABS(HMIN)) THEN
              GOTO 200
            ELSE
              RFLAG = -3
              GOTO 1000
            ENDIF
          ENDIF      

          PMID = P(XMID, YMID)
          FMID = F(XMID, YMID)

C         Find solution at X2 with a step of H/2 from XMID                         
C         We already have the limits for this point.  And we 
C         know that Y2 is within the proper limits. 
          YT = Y2 
          CALL SETUP_FOR_SECANT(
     I                          YT, YUP2, YDN2,
     O                          Y2END, ROOT)

C         Find solution at X2 with second step of H/2
          CALL GET_ROOT
     I                 (H2, FMID, PMID, X2, Y2END, F, P, EPSF, 
     I                  EPSARG, EPSABS,  YUP2, YDN2,
     M                  ROOT, 
     O                  RFLAG)
          Y2END = ROOT
C          WRITE(STDOUT,*) ' Y2END=',Y2END,' RFLAG=',RFLAG

          IF(RFLAG.EQ.0) THEN
C           No solution found.  Reduce the step and try again.
            H = 0.9D0*H
            IF(ABS(H).GT.ABS(HMIN)) THEN
              GOTO 200
            ELSE
              RFLAG = -4
              GOTO 1000
            ENDIF
          ENDIF  

          P2END = P(X2, Y2END)

C         We now have the two estimates at X2: Y2 and Y2 end.              

          RELERR = ABS((P2 - P2END)/(3.D0*P2END))
C          WRITE(STDOUT,*) ' RELERR=',RELERR,' EPS=',EPS
C          WRITE(STDOUT,*) ' H=',H
C          WRITE(STDOUT,*) ' '


          IF(RELERR.LE.EPS) GOTO  300

C         Error is too large.  Reduce the step and try again. 

          H = (0.9D0*(EPS/RELERR)**0.33333333333333D0)*H

          IF(ABS(H).LT.ABS(HMIN)) THEN
            RFLAG = -9
            WRITE(STDOUT,*) ' Y2=',Y2,' Y2END=',Y2END
            WRITE(STDOUT,*) ' YCRIT2=',YCRIT2,' YNORM2=', YNORM2
            GOTO 1000
          ENDIF
          GOTO 200

300     CONTINUE


        Y2BEST = Y2END

C       We have completed a step of H with an estimated value
C       of P, P2END, that satisfies the local tolerance for
C       truncation error.  The value of Y, Y2END is the 
C       current best estimate.  Do optional extrapolation correction.
        IF(EXTRAP.EQ.1) THEN
          P2BEST = P2END + (P2END - P2)/3

C         Solve for the value of Y that matches P2BEST and that 
C         also falls within the valid limits for a solution. 

          YT = Y2BEST
          CALL SETUP_FOR_SECANT(
     I                          YT, YUP2, YDN2,
     O                          Y2BEST, ROOT)
          CALL INVERT_P(
     I                  P2BEST, X2, Y2BEST, P, YUP2, YDN2, EPSF, 
     I                  EPSABS, EPSARG, 
     M                  ROOT,
     O                  RFLAG)
          IF(RFLAG.EQ.0) THEN
C           No solution found.  Reduce the step and try again.
            H = 0.9D0*H
            IF(ABS(H).GT.ABS(HMIN)) THEN
             
              GOTO 200
            ELSE
              RFLAG = -4
              GOTO 1000
            ENDIF
          ENDIF  
          Y2BEST = ROOT
        ENDIF

C       At this point we have found a solution at X2.  Store the
C       point and prepare for the next step, if there is one.
        
        I = I + INC
        IF(I.GT.NMAX) THEN
          RFLAG = -1
          N = I - 1
          GOTO 1000
        ELSEIF(I.LT.1) THEN
          RFLAG = -1
          N = I + 1
          GOTO 1000
        ENDIF 

        N = I
        XVEC(I) = X2
        YVEC(I) = Y2BEST
        YCVEC(I) = YCRIT2
        YNVEC(I) = YNORM2

        IF(H.GT.0) THEN
          IF(X2.GE.XE) THEN
            GOTO 1000
          ENDIF
        ELSE
          IF(X2.LE.XE) THEN
            GOTO 1000
          ENDIF
        ENDIF

C       Check to see if we are close to critical depth but only 
C       when the course of the computations can attain critical depth.

        IF(YCRIT_YES.EQ.1) THEN
          IF(ABS(Y2BEST - YCRIT2)/YCRIT2.LT.0.01D0) THEN
            RFLAG = 2
            GOTO 1000
          ENDIF  
        ENDIF
C       Adjust the step
        IF(RELERR.GT.0.D0) THEN
          HFAC = 0.9*(EPS/RELERR)**0.3333333333D0
          IF(HFAC.LE.DHLIM) THEN
            H = HFAC*H
          ELSE
            H = DHLIM*H
          ENDIF
        ELSE
          H = DHLIM*H
        ENDIF

        IF(ABS(H).GT.ABS(HMAX)) THEN
          H = HMAX
        ENDIF
        X1 = X2
        Y1 = Y2BEST
        GOTO 100

1000  CONTINUE
      IF(RFLAG.LT.1) THEN
        CALL ERRORS_IN_IMPTRAP(
     I                         STDOUT, RFLAG, N , X1, NMAX, H, H2)
      ENDIF

      RETURN
      END                

C
C
C
