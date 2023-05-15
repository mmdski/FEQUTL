C
C
C
      SUBROUTINE   FDROOT
     I                   (B, FUN, EPSF,
     M                    A,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Find a sign change in the function FUN starting with
C     the limits A and B.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG
      REAL A, B, EPSF
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL FUN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     B      - upper limit for finding a root
C     FUN    - function involved in root search
C     EPSF   - Convergence tolerance for function values
C     A      - initial lower limit and final lower limit
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL FA, FB, FM, XM
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, FLOAT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FUN
C***********************************************************************
      FA = FUN(A)
      FB = FUN(B)
      IF(FA*FB.LE.0.0) RETURN
      IF(ABS(FA).LE.EPSF.OR.ABS(FB).LE.EPSF) RETURN
 
      DO 200 I=-8,15
        XM = A + (B - A)*FLOAT(I)/16.0
        IF(XM.LE.0.0) GOTO 200
        FM = FUN(XM)
C        WRITE(6,*) ' FDROOT: XM=',XM,' FM=',FM
        IF(FM*FB.LE.0.0) THEN
          A = XM
          RETURN
        ENDIF
 200  CONTINUE
      EFLAG = 1
      RETURN
      END
C
C
C
      SUBROUTINE   REGFAL
     I                   (EPSX, EPSF, F,
     M                    A, B,
     O                    XM, FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the sign of F() must differ
C     at the ends of the interval. The root is returned as XM.
C     Uses modified regula falsi.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL A, B, EPSF, EPSX, XM
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL F
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag. If FLAG=0, solution claimed; if FLAG=1, 
C         no solution possible; if FLAG=2, more than 100 iterations.
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL FL, FM, FMOLD, FR, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C***********************************************************************
      FLAG = 0
      XL = A
      XR = B
 
      FL = F(XL)
      FR = F(XR)
      IF(ABS(FL).LE.EPSF) THEN
        XM = XL
        RETURN
      ENDIF
      IF(ABS(FR).LE.EPSF) THEN
        XM = XR
        RETURN
      ENDIF
      FMOLD = FL
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 
C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)
        FM = F(XM)
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LT.EPSX) THEN
          A = XL
          B = XR
          RETURN
        ENDIF
 
        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(FL.EQ.0.0.OR.FR.EQ.0.0) THEN
          A = XL
          B = XR
          RETURN
        ENDIF
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          FLAG = 2
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.5*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0) THEN
C           YES
            FR = 0.5*FR
          ENDIF
          FMOLD = FM
        ENDIF
        GOTO 100
      END
C
C
C
      SUBROUTINE   REGFLT
     I                   (EPSX, EPSF, F,
     M                    A, B, FL, FR,
     O                    XM, FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the signs of FL and FR must differ
C     because they are the values of F() at the ends of the interval.
C     The root is returned as XM. Uses modified regula falsi.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL A, B, EPSF, EPSX, FL, FR, XM
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL F
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     FL     - Function value on the left end of the interval
C     FR     - Value of function on the right end of the interval
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL FM, FMOLD, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C***********************************************************************
      FLAG = 0
      XL = A
      XR = B
 
      FMOLD = FL
 
C      WRITE(STDOUT,53)
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 
C        WRITE(STDOUT,52) XL, XR, FL, FR
 
C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)
        FM = F(XM)
        IF(FM.LT.-0.99E30) THEN
          FLAG = 3
          RETURN
        ENDIF
 
        IF(FM.EQ.0.0) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LT.EPSX) THEN
C         When argument collapse causes termination, return
C         the argument value that gives the minimum absolute value
C         of the residual function.
          A = XL
          B = XR
          IF(ABS(FL).LT.ABS(FR)) THEN
C           Take the left point as the answer.
            XM = XL
          ELSE
C           Take the right point as the answer.  Note that the
C           calling routine may check for the residual at convergence.
C           This is always taken to be in FL.
            XM = XR
            FL = FR
          ENDIF
          RETURN
        ENDIF
 
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          FLAG = 2
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.9*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0) THEN
C           YES
            FR = 0.9*FR
          ENDIF
          FMOLD = FM
        ENDIF
        GOTO 100
      END
C
C
C
      SUBROUTINE   RGF
     I                (EPSX, EPSF, F,
     M                 A, B, FL, FR,
     O                 XM, FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the signs of FL and FR must differ
C     because they are the values of F() at the ends of the interval.
C     The root is returned as XM. Uses modified regula falsi.
 
C     same as REGFLT- needed to avoid indirect recursion.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL A, B, EPSF, EPSX, FL, FR, XM
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL F
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     FL     - Function value on the left end of the interval
C     FR     - Value of function on the right end of the interval
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL FM, FMOLD, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C***********************************************************************
      FLAG = 0
      XL = A
      XR = B
 
      FMOLD = FL
 
C      WRITE(STDOUT,53)
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 
C        WRITE(STD6,52) XL, XR, FL, FR
C52    FORMAT(' RGF: XL=',F10.4,' XR=',F10.4,
C     A       ' FL=',F10.4,' FR=',F10.4)
 
C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)
        FM = F(XM)
 
        IF(FM.LT.-0.99E30) THEN
          FLAG = 3
          RETURN
        ENDIF
        IF(FM.EQ.0.0) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LT.EPSX) THEN
C         Argument collapse .
 
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          FL = FM
          FLAG = 2
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.9*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0) THEN
C           YES
            FR = 0.9*FR
          ENDIF
          FMOLD = FM
        ENDIF
        GOTO 100
      END
C
C
C
      SUBROUTINE   FDBLRGF
     I                   (EPSX, EPSF, F,
     M                    A, B, FL, FR,
     O                    ROOT, FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the signs of FL and FR must differ
C     because they are the values of F() at the ends of the interval.
C     The root is returned as ROOT. Uses modified regula falsi.
C     Version of RGF with full double precision.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL*8 A, B
      REAL*8 EPSF, EPSX, FL, FR, ROOT
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL*8 F
 

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     FL     - Function value on the left end of the interval
C     FR     - Value of function on the right end of the interval
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      
      REAL*8 FM, FMOLD, XL, XR, XM
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C***********************************************************************
      FLAG = 0
      XL = A
      XR = B
 
      FMOLD = FL
 
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 

C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)

C       Get the function at the current intermediate argument. 
        FM = F(XM)


        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          ROOT = XM
          FR = DBLE(KNT)
          RETURN
        ENDIF
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LT.EPSX) THEN
C         Argument collapse .
 
          A = XL
          B = XR
          FL = FM
          ROOT = XM
          FR = DBLE(KNT)
          IF(ABS(FM).GT.2.D0*EPSF) THEN
            FLAG = 3
          ENDIF
          RETURN
        ENDIF
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          ROOT = XM
          FLAG = 2
          FR = DBLE(KNT)
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.75D0*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0D0) THEN
C           YES
            FR = 0.75D0*FR
          ENDIF
          FMOLD = FM
        ENDIF

        GOTO 100
      END
C
C
C
      SUBROUTINE   DBLRGF
     I                   (EPSX, EPSF, F,
     M                    A, B, FL, FR, XM,
     O                    FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the signs of FL and FR must differ
C     because they are the values of F() at the ends of the interval.
C     The root is returned as XM. Uses modified regula falsi.
C     Version of RGF with double precision function evaluation. 
C     XM enters with value in this version.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL A, B
      REAL*8 EPSF, EPSX, FL, FR, XM
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL*8 F

      INCLUDE 'stdun.cmn'
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     FL     - Function value on the left end of the interval
C     FR     - Value of function on the right end of the interval
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      
      REAL*8 FM, FMOLD, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, DBLE
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C     *****************************FORMATS******************************
50    FORMAT(' DBLRGF: XL=',1PE12.5,' FL=',1PE12.5,' XR=',
     A        1PE12.5,' FR=',1PE12.5)
52      FORMAT(' DBLRGF: XM=',1PE12.5,' FM=',1PE12.5) 
C***********************************************************************
      FLAG = 0
      XL = DBLE(A)
      XR = DBLE(B)
 
      FMOLD = FL
 
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 
C        WRITE(STD6,50) XL, FL, XR, FR

C       Get the function at the current intermediate argument. 

        FM = F(XM)

        IF(FM.EQ.0.0D0) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LT.EPSX) THEN
C         Argument collapse .
 
          A = XL
          B = XR
          FL = FM
          IF(ABS(FM).GT.2.*EPSF) THEN
            FLAG = 3
          ENDIF
          RETURN
        ENDIF
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          FLAG = 2
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.5*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0) THEN
C           YES
            FR = 0.5*FR
          ENDIF
          FMOLD = FM
        ENDIF

C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)
        GOTO 100
      END
C
C
C
      SUBROUTINE   RGF3
     I                 (EPSX, EPSF, F,
     M                  A, B, FL, FR,
     O                  XM, FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the signs of FL and FR must differ
C     because they are the values of F() at the ends of the interval.
C     The root is returned as XM. Uses modified regula falsi.
 
C     Same as RGF- needed to avoid indirect recursion.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL A, B, EPSF, EPSX, FL, FR, XM
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL F
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     FL     - Function value on the left end of the interval
C     FR     - Value of function on the right end of the interval
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL FM, FMOLD, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C***********************************************************************
      FLAG = 0
      XL = A
      XR = B
 
      IF(ABS(FL).LE.EPSF) THEN
        XM = A
        RETURN
      ELSEIF(ABS(FR).LE.EPSF) THEN
        XM = B
        FL = FR
        RETURN
      ENDIF
 
      FMOLD = FL
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 
C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)
        FM = F(XM)
 
        IF(FM.LT.-0.99E30) THEN
          FLAG = 3
          RETURN
        ENDIF
        IF(FM.EQ.0.0) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LT.EPSX) THEN
C         When argument collapse causes termination, return
C         the argument value that gives the minimum absolute value
C         of the residual function.
          A = XL
          B = XR
          IF(ABS(FL).LT.ABS(FR)) THEN
C           Take the left point as the answer.
            XM = XL
          ELSE
C           Take the right point as the answer.  Note that the
C           calling routine may check for the residual at convergence.
C           This is always taken to be in FL.
            XM = XR
            FL = FR
          ENDIF
          RETURN
        ENDIF
 
        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          FLAG = 2
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.5*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0) THEN
C           YES
            FR = 0.5*FR
          ENDIF
          FMOLD = FM
        ENDIF
        GOTO 100
      END
C
C
C
      SUBROUTINE   RGF5
     I                 (EPSX, EPSF, F,
     M                  A, B, FL, FR,
     O                  XM, FLAG)
 
C     + + + PURPOSE + + +
C     Find at least one root of the function, F(), in the
C     interval (A, B).  On entry the signs of FL and FR must differ
C     because they are the values of F() at the ends of the interval.
C     The root is returned as XM. Uses modified regula falsi.
 
C     same as RGF- needed to avoid indirect recursion.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG
      REAL A, B, EPSF, EPSX, FL, FR, XM
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      REAL F
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     EPSX   - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     F      - Function defining the root being sought
C     A      - Left end of interval containing at least one root
C     B      - Right end of interval containing at least one root
C     FL     - Function value on the left end of the interval
C     FR     - Value of function on the right end of the interval
C     XM     - interpolated point for next residual function evaluation
C     FLAG   - Result flag
 
C      INCLUDE 'stdun.cmn'

C     + + + LOCAL VARIABLES + + +
      INTEGER KNT
      REAL LOCAL_EPSX, FM, FMOLD, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL F
C***********************************************************************
      LOCAL_EPSX = EPSX
      IF(LOCAL_EPSX.LT.1.E-6) THEN
C       Allow argument collapse to roundoff even though the calling
C       program requests a smaller tolerance.  Used to avoid 
C       non-convergence on residual value when the residual 
C       function is changing rapidly near the root.  In these cases
C       roundoff differences in argument value can cause residual
C       differences that exceed the residual tolerance.
        LOCAL_EPSX = 1.E-6
      ENDIF
      FLAG = 0
      XL = A
      XR = B
 
      FMOLD = FL
 
C      WRITE(STDOUT,53)
 
      IF(FL*FR.GT.0.0) THEN
C       NO SOLUTION POSSIBLE
        FLAG = 1
        RETURN
      ENDIF
 
      KNT = 0
 100  CONTINUE
 
C        WRITE(STD6,52) KNT, XL, XR, FL, FR
C52    FORMAT(' RGF5: KNT=',I5,' XL=',F10.5,' XR=',F10.5,' FL=',F10.5,
C     A             ' FR=',F10.5)
 
C       COMPUTE THE INTERMEDIATE POINT LOCATION
 
        XM = (FR*XL - FL*XR)/(FR - FL)
        FM = F(XM)
C        WRITE(STD6,54) FM, EPSF
C54    FORMAT(' FM=',F10.6,' EPSF=',F10.6)
        IF(FM.LT.-0.99E30) THEN
          FLAG = 3
          RETURN
        ENDIF
        IF(FM.EQ.0.0) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        IF(ABS(XL - XR)/(ABS(XL) + ABS(XR)).LE.LOCAL_EPSX) THEN
C         When argument collapse causes termination, return
C         the argument value that gives the minimum absolute value
C         of the residual function.
          A = XL
          B = XR
          IF(ABS(FL).LT.ABS(FR)) THEN
C           Take the left point as the answer.
            XM = XL
          ELSE
C           Take the right point as the answer.  Note that the
C           calling routine may check for the residual at convergence.
C           This is always taken to be in FL.
            XM = XR
            FL = FR
          ENDIF
          RETURN
        ENDIF
 
        IF(ABS(FM).LT.EPSF) THEN
          A = XL
          B = XR
          FL = FM
          RETURN
        ENDIF
 
        KNT = KNT + 1
        IF(KNT.GT.100) THEN
          A = XL
          B = XR
          FLAG = 2
          RETURN
        ENDIF
        IF(FL*FM.LE.0.0) THEN
C         RETAIN LEFT POINT TO RETAIN SIGN CHANGE
 
          XR = XM
          FR = FM
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FMOLD*FM.GT.0.0) THEN
C           YES
            FL = 0.5*FL
          ENDIF
          FMOLD = FM
        ELSE
C         RETAIN RIGHT POINT TO RETAIN SIGN CHANGE
 
          XL = XM
          FL = FM
 
C         HAVE INTERMEDIATE POINTS BEEN ON THE SAME SIDE OF THE ROOT
C         TWICE IN SUCCESSION?
          IF(FM*FMOLD.GT.0.0) THEN
C           YES
            FR = 0.5*FR
          ENDIF
          FMOLD = FM
        ENDIF
        GOTO 100
      END
C
C
C
      SUBROUTINE   SECANT
     I                   (XLA, EPSARG, EPSF, EPSABS, MAXIT, XMIN, XMAX, 
     I                    FUN,
     M                    XRA,
     O                    FLAG)
 
C     + + + PURPOSE + + +
C     Find a root of the function, FUN, using the secant method.
C     XLA and XRA give the starting points.  XMAX and XMIN give the
C     maximum and minimum values allowed for the root.  The final
C     root is returned in XRA.  FUN is a double precision function
C     of one argument.  All other values must be passed to FUN using
C     a common block.  FLAG = 0: solution found. FLAG=1: MAXIT
C     interations and no solution. FLAG=2: computed derivative is
C     zero and the function is still larger than the tolerance.
C     Indicates a local max or min of the function.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER FLAG, MAXIT
      REAL EPSARG, EPSABS, EPSF, XLA, XMAX, XMIN, XRA
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      DOUBLE PRECISION FUN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     XLA    - Left hand end of starting interval
C     EPSARG - Convergence tolerance on arguments
C     EPSF   - Convergence tolerance for function values
C     EPSABS - Absolute convergence tolerance on arguments
C     MAXIT  - Maximum number of iterations
C     XMIN   - Lower limit for root
C     XMAX   - Upper limit for root
C     FUN    - Integrand function
C     XRA    - Right hand end of starting interval
C     FLAG   - Result flag
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      DOUBLE PRECISION DF, DX, DXTEMP, FL, FR, OLDVAL, XL, XR
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FUN
C***********************************************************************
      XL = XLA
      XR = XRA
      OLDVAL = 1.D0
C      WRITE(STD6,*) ' ENTERING SECANT. XL=',XL,' XR=',XR
 
      FLAG = 0
      FL = FUN(XL)
C      WRITE(STD6,*) ' FL=',FL
      DX = XR - XL
      DO 100 I=1,MAXIT
        FR = FUN(XR)
C       WRITE(STD6,*) ' I=',I,' FR=',FR
        IF(ABS(FR).LE.EPSF) THEN
          XRA = XR
          RETURN
        ENDIF
        DF = FL - FR
C        WRITE(STD6,*) ' SECANT:I=',I,' FL=',FL,' FR=',FR,' DX=',DX
        IF(DF.NE.0.0) THEN
          DXTEMP = 0.0
          OLDVAL = DX/DF
          DX = FR*OLDVAL
          XL = XR
          FL = FR
          XR = XR + DX
          IF(XR.GE.XMAX) THEN
            XR = 0.5*(XL + XMAX)
            DXTEMP = XR - XL
          ELSEIF(XR.LE.XMIN) THEN
            XR = 0.5*(XL + XMIN)
            DXTEMP = XR - XL
          ENDIF
C          WRITE(STD6,*) ' DX=',DX
          IF(ABS(DX/XR).LE.EPSARG.OR.ABS(DX).LE.EPSABS) THEN
            XRA = XR
            RETURN
          ENDIF
          IF(DXTEMP.NE.0.0) THEN
            DX = DXTEMP
          ENDIF
        ELSE
C         DF = 0. SIGNAL POSSIBLE ERROR
          XRA = XR
          IF(ABS(FR).LE.EPSF) THEN
            RETURN
          ELSE
            FLAG = 2
            RETURN
          ENDIF
        ENDIF
 
 100  CONTINUE
 
C     DROP THROUGH INDICATES NO CONVERGENCE WITHIN MAXIT ITERATIONS
      FLAG = 1
      XRA = XR
      RETURN
      END
