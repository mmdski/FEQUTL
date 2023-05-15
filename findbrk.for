C     ***********
C     *         *
C     * RATIOPNT
C     *         *
C     ***********

      SUBROUTINE RATIOPNT(N, A, B, 
     O                    TAU)

C     + + + PURPOSE + + +
C     Compute distribution of points that have a constant ratio.

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N

      REAL*8 A, B, TAU(N)

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N - number of points
C     A - lower limit of point range
C     B - upper limit of point range
C     TAU- sequence of points having a constant ratio.

C     + + + LOCAL VARIABLES + + +
      INTEGER J

      REAL*8  R

C***********************************************************************
      R = EXP(LOG(B/A)/DBLE(N-1))
      TAU(1) = A
      DO 100 J=2,N-1
        TAU(J) = TAU(J-1)*R
100   CONTINUE
      TAU(N) = B
      RETURN
      END 

C     ***********
C     *         *
C     * CUBIC_XPLINE_LOOKUP
C     *         *
C     ***********

      SUBROUTINE CUBIC_SPLINE_LOOKUP(ARG, N, XVEC, FVEC, FPVEC, FPPVEC,
     M                               LAST,
     O                               F, FP, FPP)

C     + + + PURPOSE + + +
C     Lookup an argument, ARG, in a cubic spline and return the 
C     function value, first derivative, and second derivative. 
C     We assume that FPVEC and FPPVEC have been computed properly.

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, LAST

      REAL*8 ARG, XVEC(N), FVEC(N), FPVEC(N), FPPVEC(N), F, FP, FPP

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARG - argument at which values are to be found
C     N - number of breakpoints for cubic spline
C     XVEC - breakpoint argument sequence
C     FVEC - function value at breakpoints
C     FPVEC - first derivative values at breakpoints
C     FPPVEC - second derivative values at breakpoints
C     LAST - index to the breakpoint at the left of the 
C            last breakpoint interval, panel, found on
C            previous call. 
C     F - function value at ARG
C     FP - derivative value at ARG
C     FPP - second derivative value at ARG


C     + + + LOCAL VARIABLES + + +
      INTEGER L
      REAL*8 DX, HALFDX, FPL, FPPL, P
C***********************************************************************
      L = LAST
C     Find the interval that contains the argument.  Assume well-behaved
C     calls, that is no checking for being outside of range. 

      IF(ARG.GE.XVEC(L)) THEN
100     CONTINUE
          IF(ARG.GT.XVEC(L+1)) THEN
            L = L + 1
            GOTO 100
          ENDIF
      ELSE
200     CONTINUE
          L = L - 1
          IF(ARG.LT.XVEC(L)) GOTO 200
      ENDIF

      LAST = L
C     L should point to the left end of the interval containing ARG. 
      
      DX = ARG - XVEC(L)
      HALFDX = 0.5*DX
      P = DX/(XVEC(L+1) - XVEC(L))
      FPL = FPVEC(L)
      FPPL = FPPVEC(L) 

      FPP = FPPL + P*(FPPVEC(L+1) - FPPL)
      FP = FPL + HALFDX*(FPPL + FPP)
      F = FVEC(L) + HALFDX*( FPL + FP - HALFDX*(FPP - FPPL)/3.D0)
      RETURN
      END
            
C     ***********
C     *         *
C     * PARABOLIC_SPLINE_LOOKUP
C     *         *
C     ***********

      SUBROUTINE PARABOLIC_SPLINE_LOOKUP(ARG, N, XVEC, FVEC, FPVEC,
     M                               LAST,
     O                               F, FP)

C     + + + PURPOSE + + +
C     Lookup an argument, ARG, in a parabolic spline and return the 
C     function value, and first derivative
C     We assume that FPVEC  have been computed properly.

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, LAST

      REAL*8 ARG, XVEC(N), FVEC(N), FPVEC(N), F, FP

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARG - argument at which values are to be found
C     N - number of breakpoints for cubic spline
C     XVEC - breakpoint argument sequence
C     FVEC - function value at breakpoints
C     FPVEC - first derivative values at breakpoints
C     LAST - index to the breakpoint at the left of the 
C            last breakpoint interval, panel, found on
C            previous call. 
C     F - function value at ARG
C     FP - derivative value at ARG

C     + + + LOCAL VARIABLES + + +
      INTEGER L
      REAL*8 DX, HALFDX, FPL, P
C***********************************************************************
      L = LAST
C     Find the interval that contains the argument.  Assume well-behaved
C     calls.

      IF(ARG.GE.XVEC(L)) THEN
100     CONTINUE
          IF(ARG.GT.XVEC(L+1)) THEN
            L = L + 1
            GOTO 100
          ENDIF
      ELSE
200     CONTINUE
          L = L - 1
          IF(ARG.LT.XVEC(L)) GOTO 200
      ENDIF

      LAST = L
C     L should point to the left end of the interval containing ARG. 
      
      DX = ARG - XVEC(L)
      HALFDX = 0.5*DX
      P = DX/(XVEC(L+1) - XVEC(L))
      FPL = FPVEC(L)

      FP = FPL + P*(FPVEC(L+1) - FPL)
      F = FVEC(L) + HALFDX*(FPL + FP)
      RETURN
      END


C     ***********
C     *         *
C     * LINEAR_SPLINE_LOOKUP
C     *         *
C     ***********

      SUBROUTINE LINEAR_SPLINE_LOOKUP(ARG, N, XVEC, FVEC,
     M                               LAST,
     O                               F)

C     + + + PURPOSE + + +
C     Lookup an argument, ARG, in a linear spline and return the 
C     function value.

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, LAST

      REAL*8 ARG, XVEC(N), FVEC(N), F

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARG - argument at which values are to be found
C     N - number of breakpoints for cubic spline
C     XVEC - breakpoint argument sequence
C     FVEC - function value at breakpoints
C     LAST - index to the breakpoint at the left of the 
C            last breakpoint interval, panel, found on
C            previous call. 
C     F - function value at ARG

C     + + + LOCAL VARIABLES + + +
      INTEGER L
      REAL*8 DX, FL, P
C***********************************************************************
      L = LAST
C     Find the interval that contains the argument.  Assume well-behaved
C     calls.

      IF(ARG.GE.XVEC(L)) THEN
100     CONTINUE
          IF(ARG.GT.XVEC(L+1)) THEN
            L = L + 1
            GOTO 100
          ENDIF
      ELSE
200     CONTINUE
          L = L - 1
          IF(ARG.LT.XVEC(L)) GOTO 200
      ENDIF

      LAST = L
C     L should point to the left end of the interval containing ARG. 
      
      DX = ARG - XVEC(L)
      P = DX/(XVEC(L+1) - XVEC(L))
      FL = FVEC(L)

      F = FL + P*(FVEC(L+1) - FL)
      RETURN
      END
C     
C     ***********
C     *         *
C     * FIND_FPP
C     *         *
C     ***********

      SUBROUTINE FIND_FPP(N, X, F, FP,
     O              FPP)

C     + + + PURPOSE + + +
C     Compute the second derivative at the breakpoints of a cubic
C     spline given the first derivatives.

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      REAL*8 X(N), F(N), FP(N), FPP(N)

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N - number breakpoints in the cubic spline
C     X - breakpoint sequence
C     F - function sequence
C     FP - derivative sequence
C     FPP - second derivative sequence

C     + + + LOCAL VARIABLES + + +
      INTEGER J
      REAL*8 H
C*******************************************************************
C     Do all but last point using the panel to the right of 
C     each point. 
      DO 100 J=1,N-1
        H = X(J+1) - X(J)
        FPP(J) = (6.*(F(J+1) - F(J))/H - 4.*FP(J) - 2.*FP(J+1))/H
100   CONTINUE
C     Do last point using the panel to the left of that point.
C     Note that the value of H is valid here.
      FPP(N) = (2.*FP(N-1) +4.*FP(N) -6.*(F(N) - F(N-1))/H)/H
      RETURN
      END

C     ***********
C     *         *
C     * LOOKUP_IG
C     *         *
C     ***********

      REAL*8 FUNCTION LOOKUP_IG(ARG) 

C     + + + PURPOSE + + +
C     Lookup the indefinite integral computed using the mid-point
C     rule on the integrand, G, where G is computed from the
C     cubic-spline fit to a series of points describing the 
C     function values.  

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL*8 ARG

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARG - argument for the lookup
      
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'goodbrk.cmn'

C     + + + LOCAL VARIABLES + + +
      INTEGER L, K
      REAL*8 DX_CS, HALFDX_CS, FPL, FPPL, P, F, FP, FPP, ARG_QD, 
     A       DX_QD, G
C***********************************************************************
      L = IG_LAST
C     Find the interval that contains the argument.  Assume well-behaved
C     calls.

      IF(ARG.GE.XVEC_QD(L)) THEN
100     CONTINUE
          IF(ARG.GT.XVEC_QD(L+1)) THEN
            L = L + 1
            GOTO 100
          ENDIF
      ELSE
200     CONTINUE
          L = L - 1
          IF(ARG.LT.XVEC_QD(L)) GOTO 200
      ENDIF

      IG_LAST = L
C     L should point to the left end of the interval containing ARG. This is
C     in the argument series for the integral, NOT the argument series for
C     the cubic spline.  Compute the argument value for the evaluation of
C     the integrand. 
      DX_QD = ARG - XVEC_QD(L)
      ARG_QD = XVEC_QD(L) + 0.5*DX_QD

C     Find the integrand value at ARG_QD.  Get the cubic-spline pointer.
      
      K = QD_PNT(L)
      DX_CS = ARG_QD - XVEC_CS(K)
      HALFDX_CS = 0.5*DX_CS
      P = DX_CS/(XVEC_CS(K+1) - XVEC_CS(K))
      FPL = FPVEC(K)
      FPPL = FPPVEC(K) 

      FPP = FPPL + P*(FPPVEC(K+1) - FPPL)
      FP = FPL + HALFDX_CS*(FPPL + FPP)
      F = FVEC_CS(K) + HALFDX_CS*( FPL + FP - HALFDX_CS*
     A            (FPP - FPPL)/3.D0)
      G = SQRT(ABS(FPP)/F)
      
      LOOKUP_IG = IGVEC(L) + G*DX_QD
      
      RETURN
      END
C     ***********
C     *         *
C     * IG_RES
C     *         *
C     ***********

      REAL*8 FUNCTION IG_RES(ARG)

C     + + + PURPOSE + + +
C     Residual function for finding good linear-spline breakpoints

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL*8 ARG

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARG - argument for the residual computation

C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'goodbrk.cmn'

C     + + +FUNCTIONS+ + +
      REAL*8 LOOKUP_IG
      EXTERNAL LOOKUP_IG
C**********************************************************************
      IG_RES = (IGVAL - LOOKUP_IG(ARG))/IGVAL
      RETURN
      END

C     ***********
C     *         *
C     * FIND_GOOD_POINTS
C     *         *
C     ***********

      SUBROUTINE FIND_GOOD_POINTS(STDOUT, NPNTS, XVEC, FVEC, 
     I       RELATIVE_ERROR, LEFT_SLOPE, RIGHT_SLOPE,
     M       N_GOOD,
     O       GOOD_POINTS, NO_GOOD, MAX_RERR)

C     + + + PURPOSE + + +
C     Given a series of arguments and function values, XVEC and FVEC,
C     find a sequence of good breakpoints for interpolating 
C     the function represented by the given values using a linear
C     spline and with the target of having an interpolation error
C     of no more than RELATIVE_ERROR. 

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNTS, N_GOOD, NO_GOOD, STDOUT
      REAL*8 XVEC(NPNTS), FVEC(NPNTS), RELATIVE_ERROR, LEFT_SLOPE,
     A       RIGHT_SLOPE, GOOD_POINTS(N_GOOD), MAX_RERR
      
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - standard output unit
C     NPNTS - number of values in XVEC and FVEC
C     XVEC - series of argument values
C     FVEC - function values at arguments in XVEC
C     RELATIVE_ERROR - target relative error for linear spline
C                      interpolation in the function tabulated
C                      in XVEC and FVEC
C     LEFT_SLOPE - slope at left end of range
C     RIGHT_SLOPE - slope at right end of range
C     N_GOOD - number of points in the good sequence
C     GOOD_POINTS - the sequence of good breakpoints
C     NO_GOOD - error flag. =1 if failure, =0 otherwise
C     MAX_RERR - absolute value of the maximum relative error
C               in linear interpolation using the GOOD_POINTS

C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'goodbrk.cmn'


C     + + + LOCAL VARIABLES + + +
C     Define some patterns because we need many different
C     but closely related values:

C     Values defining the cubic spline fit to the data
C     points

C       XVEC  - the argument values
C       FVEC  - the function values
C       FPVEC - the first derivative values
C       FPPVEC- the second derivative values.

C     FPVEC, and FPPVEC are in common block: goodbrk.cmn.
C     FVEC_CS and XVEC_CS are also in goodbrk.cmn and are the same
C     as FVEC and XVEC.

C     The indefinite integral variables are:

C     IGVEC - indefinite integral.  g(x) is the integrand
C             and so IG denotes 'integral of g'
C     XVEC_QD-  argument points for integral.  Contains
C               all points in XVEC plus more.
C     QD_PNT- pointer to the left panel point of the
C             cubic spline of the panel that contains
C             the point.  If XL is the left panel point
C             for a panel in the cubic spline, and XR is
C             the right panel point for the same panel,
C             then all points  XL <= x < XR will have
C             a pointer to the left panel point, XL. 
C       These are in goodbrk.cmn.

C     We will also have a linear spline approximation.
C     GOOD_POINTS  - argument values
C     FVEC_LS  - Function values. 

      INTEGER I, J, K, FLAG, LAST, LAST2, N_GOOD_MAX, NIN, N_TO_ADD,
     A        L, NBIG, IFLAG
      PARAMETER (NBIG=2000)

      REAL*8  ARG, F, FP, FPP, G, XL, XR, RFL, RFR, FPPL, XROOT,
     A        DELTA_IG, EPSF, EPSX, ROOT, FVEC_LS(PMXNHU),
     B        RERR, FHAT, DX, HALFDX, P, FPL,
     C        INTERVAL_LENGTH, OFFSET(3), NEWX(NBIG), S(NBIG,6),
     D        COEF(4,NBIG),GAMMA, TESTFPP(NBIG)

C     + + + FUNCTIONS + + +
      REAL*8 IG_RES
      EXTERNAL IG_RES

      DATA OFFSET/0.4D0, 0.5D0, 0.6D0/, EPSX/0.000001D0/,EPSF/0.00001D0/
C     ******************************FORMATS*****************************
 30   FORMAT(I5,F10.3,I5)
 31   FORMAT(/,' *BUG in FIND_GOOD_POINTS: No sign change on entry',
     A         ' to FDBLRGF.')
 32   FORMAT(/,' *BUG in FIND_GOOD_POINTS: More than 100 iterations',
     A         ' in FDBLRGF.')
 33   FORMAT(/,' *BUG in FIND_GOOD_POINTS: Argument collapse with ',
     A         'residual more than twice EPSF in FDBLRGF.')
 50   FORMAT(' ',I5,F12.7,F10.4, 1PE12.4,1PE12.4,1PE12.4,4(1PE12.4))
 52   FORMAT(' ', F12.7, 1PE12.4, 1PE12.4)
 54   FORMAT(F10.6,F10.4,F10.4,F10.4)
 58   FORMAT(/,' Results from cubic-spline fit.')
 60   FORMAT(' Index    Argument  Function    1stDeriv    2ndDeriv',
     A      '   Integrand')
 62   FORMAT(/,' Computation of indefinite integral of breakpoint',
     A         ' density.')
 64   FORMAT(' Index    Argument    Integral   Integrand    Function',
     A    '    1stDeriv    2ndDeriv')
 66   FORMAT(' ',I5,F12.6,F12.6)
 68   FORMAT(30X,F12.5,F12.5,1PE12.4,1PE12.4)
 70   FORMAT(/,' *BUG in FIND_GOOD_POINTS: Space for integral of the',
     A         ' breakpoint density exhausted.')
 72   FORMAT(/,' Number of breakpoints=',I5, 
     A  ' > the maximum allowed=',I5,'.',
     A  ' This should not happen and probably indicates a',/,
     B  '  subtle error at some point.  Make sure that LIPREC is',
     C  '  greater than 0.005 and NFRAC=60 or more.')
C***********************************************************************
      IG_LAST = 1
      LAST = 1
      N_GOOD_MAX = N_GOOD
C     Fit the points with a cubic spline using the supplied derivative
C     end conditions. 

      CALL SPLINE(STDOUT, XVEC, FVEC, NPNTS, 1, LEFT_SLOPE, 1,
     I            RIGHT_SLOPE,
     O            FPVEC)

C     Compute the second derivatives. 
      CALL FIND_FPP(NPNTS, XVEC, FVEC, FPVEC,
     O              FPPVEC)

C      GAMMA = 2.5D0
C      CALL TAUTSP ( XVEC, FVEC, NPNTS, GAMMA, NBIG,
C     O                    S,
C     O                    NEWX, COEF, L, K, IFLAG )
C
C      CALL SPLINE
C     I           (STDOUT, XVEC, FPVEC, NPNTS, 3, 0.5D0, 3, 0.5D0,
C     O            TESTFPP)
C
C
C      WRITE(STDOUT,58)
C      WRITE(STDOUT,*) ' NPNTS=',NPNTS,' L=',L
C      WRITE(STDOUT,60)
C      DO 90 I=1,NPNTS
C        WRITE(STDOUT,50) I, XVEC(I), FVEC(I), FPVEC(I), 
C     A                   FPPVEC(I),
C     B                   SQRT(ABS(FPPVEC(I))/FVEC(I)),
C     C      S(I,4), TESTFPP(I)
C90    CONTINUE 

C     Transfer FVEC and XVEC to the common-block values
C     FVEC_CS and XVEC_CS
      DO 100 I=1,NPNTS
        FVEC_CS(I) = FVEC(I)
        XVEC_CS(I) = XVEC(I)
100   CONTINUE

C     Now define the argument points for the integral.  Expand on those
C     in XVEC. 
      N_TO_ADD = 2
      K = 0
      DO 105 I=1,NPNTS-1
C       Is there a sign change in FPP in the panel between
C       I and I + 1?
        IF(FPPVEC(I)*FPPVEC(I+1).LT.0.0) THEN
C         There is a sign change.  This means that the integrand
C         has a weak singularity, that is, it has an infinite 
C         derivative as its goes to zero.  This comes about because
C         of the sqrt of the absolute value involved in the
C         integrand.  Therefore we must introduce a point at the
C         zero of the integrand and also insert points in each
C         of the sub-intervals so formed. 

          XROOT = XVEC(I) - FPPVEC(I)*(XVEC(I+1) - XVEC(I))/
     A                                       (FPPVEC(I+1) - FPPVEC(I))
C         Process the left-hand subinterval.
          K = K + 1
          IF(K.GT.PMXNIG) THEN
            WRITE(STDOUT,70) 
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          XVEC_QD(K) = XVEC(I)
          QD_PNT(K) = I
C         Compute number of points to insert.  Minimum is one. 
          NIN = INT(DBLE(N_TO_ADD)*(XROOT - XVEC(I))/
     A                                     (XVEC(I+1) - XVEC(I)) + 1.0)
          DO 101 J=1,NIN
            K = K + 1
            IF(K.GT.PMXNIG) THEN
              WRITE(STDOUT,70) 
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            XVEC_QD(K) = XVEC(I) + DBLE(J)*(XROOT - XVEC(I))/
     A                                             DBLE(NIN + 1)
            QD_PNT(K) = I
101       CONTINUE
C         Process the right-hand subinterval.
          K = K + 1
          IF(K.GT.PMXNIG) THEN
            WRITE(STDOUT,70) 
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          XVEC_QD(K) = XROOT
          QD_PNT(K) = I
C         Compute number of points to insert.  Minimum is one. 
          NIN = INT(DBLE(N_TO_ADD)*(XVEC(I+1) - XROOT)/
     A                                     (XVEC(I+1) - XVEC(I)) + 1.0)
          DO 102 J=1,NIN
            K = K + 1
            IF(K.GT.PMXNIG) THEN
              WRITE(STDOUT,70) 
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            XVEC_QD(K) = XROOT + DBLE(J)*(XVEC(I+1) - XROOT)/
     A                                             DBLE(NIN + 1)
            QD_PNT(K) = I
102       CONTINUE
        ELSE
C         There is no sign change.  Transfer the point at the beginning
C         of the current panel.
          K = K + 1
          IF(K.GT.PMXNIG) THEN
            WRITE(STDOUT,70) 
            STOP 'Abnormal stop. Errors found.'
          ENDIF
          XVEC_QD(K) = XVEC(I)
          QD_PNT(K) = I
C         Add the intermediate points.
          DO 104 J=1,N_TO_ADD
            K = K + 1
            IF(K.GT.PMXNIG) THEN
              WRITE(STDOUT,70) 
              STOP 'Abnormal stop. Errors found.'
            ENDIF
            XVEC_QD(K) = XVEC(I) + DBLE(J)*(XVEC(I+1) - XVEC(I))/
     A                                             DBLE(N_TO_ADD + 1)
            QD_PNT(K) = I

104       CONTINUE
        ENDIF
105   CONTINUE
C     Transfer the last point. 
      K = K + 1
      IF(K.GT.PMXNIG) THEN
        WRITE(STDOUT,70) 
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      XVEC_QD(K) = XVEC(NPNTS)
      QD_PNT(K) = NPNTS 
      N_QD = K
      
C     Compute the indefinite integral using the mid-point rule.
C     This rule is more accurate in many cases than the trapezoidal 
C     rule, especially when weak singularities are present in the
C     integrand.  

C      WRITE(STDOUT,62) 
C      WRITE(STDOUT,64)
      IGVEC(1) = 0.D0
C      WRITE(STDOUT,66) 1, XVEC_QD(1), IGVEC(1)
      DO 110 I=1,N_QD -1
C       Find mid-point of current interval.
        INTERVAL_LENGTH = XVEC_QD(I+1) - XVEC_QD(I)
        ARG = XVEC_QD(I) + 0.5*INTERVAL_LENGTH

C       Compute the values for the cubic spline at argument
C       ARG. Get pointer into the cubic-spline description.
        K = QD_PNT(I)
        DX = ARG - XVEC(K)
        HALFDX = 0.5*DX
        P = DX/(XVEC(K+1) - XVEC(K))
        FPL = FPVEC(K)
        FPPL = FPPVEC(K) 

        FPP = FPPL + P*(FPPVEC(K+1) - FPPL)
        FP = FPL + HALFDX*(FPPL + FPP)
        F = FVEC(K) + HALFDX*( FPL + FP - HALFDX*(FPP - FPPL)/3.D0)
      
C       Compute the integrand value and the integral increment.
        G = SQRT(ABS(FPP)/(F))

C        WRITE(STDOUT,68) G, F, FP, FPP
        IGVEC(I+1) = IGVEC(I) + G*INTERVAL_LENGTH
C        WRITE(STDOUT,66) I+1, XVEC_QD(I+1), IGVEC(I+1)

110   CONTINUE

C     Estimate the number of points needed.
      N_GOOD = INT(IGVEC(N_QD)/SQRT(8.0*RELATIVE_ERROR) + 1) + 1
      IF(N_GOOD.GT.N_GOOD_MAX) THEN
C       Limit reached.
        WRITE(STDOUT,72) N_GOOD, N_GOOD_MAX
        STOP 'Abnormal stop.  Errors Found.'
      ENDIF
200   CONTINUE
C        WRITE(STDOUT,*) ' '
C        WRITE(STDOUT,*) ' N_GOOD=', N_GOOD

C       Now assign the good breakpoints so that each panel of the 
C       linear spline gets the same amount of the cumulative 
C       error density.  
        DELTA_IG = IGVEC(N_QD)/(N_GOOD - 1)
        GOOD_POINTS(1) = XVEC(1)
        GOOD_POINTS(N_GOOD) = XVEC(NPNTS)
        FLAG = 0
C        WRITE(STDOUT,30) 1, GOOD_POINTS(1), FLAG
      
        IGVAL = 0.0
        DO 150 I=2,N_GOOD-1
          IGVAL = IGVAL + DELTA_IG
          XL = GOOD_POINTS(I-1)
          RFL = IG_RES(XL)
          XR = GOOD_POINTS(N_GOOD)
          RFR = IGVAL - IGVEC(N_QD)              

          CALL FDBLRGF
     I                (EPSX, EPSF, IG_RES,
     M                 XL, XR, RFL, RFR,
     O                 ROOT, FLAG)
          IF(FLAG.GT.0) THEN
            NO_GOOD = 1
            IF(FLAG.EQ.1) THEN
              WRITE(STDOUT,31)
              STOP 'Abnormal stop. Errors found.'
            ELSEIF(FLAG.EQ.2) THEN
              WRITE(STDOUT,32)
              STOP 'Abnormal stop. Errors found.'
            ELSE
              WRITE(STDOUT,33)
            ENDIF
          ENDIF
        
          GOOD_POINTS(I) = ROOT
C          WRITE(STDOUT,30) I, ROOT, FLAG
150     CONTINUE
C        WRITE(STDOUT,30) N_GOOD, GOOD_POINTS(N_GOOD), 0

C       Do a check computation
C       Get function value at the breakpoints for the linear spline.
        DO 170 I=1, N_GOOD
          CALL CUBIC_SPLINE_LOOKUP(GOOD_POINTS(I), NPNTS, XVEC, FVEC,
     I                              FPVEC, FPPVEC,
     M                              LAST,
     O                              FVEC_LS(I), FP, FPP)
170     CONTINUE

        LAST = 1
        LAST2 = 1
        MAX_RERR = 0.D0
        XL = GOOD_POINTS(1)
        DO 180 I=2,N_GOOD
          XR = GOOD_POINTS(I)
          DX = XR - XL
          DO 175 J=1,3
            ARG = XL + DX*OFFSET(J)
            CALL CUBIC_SPLINE_LOOKUP(ARG, NPNTS, XVEC, FVEC,
     I                               FPVEC, FPPVEC,
     M                               LAST,
     O                               F, FP, FPP)
            CALL LINEAR_SPLINE_LOOKUP(ARG, N_GOOD, GOOD_POINTS, FVEC_LS,
     M                                LAST2,
     O                                FHAT)
            RERR = (FHAT - F)/F
            MAX_RERR = MAX(ABS(RERR), MAX_RERR)                      
175       CONTINUE
          XL = XR
180     CONTINUE
        IF(MAX_RERR.GT.RELATIVE_ERROR) THEN
          IF(N_GOOD.LT.N_GOOD_MAX) THEN
C           Add another breakpoint and try again.
            N_GOOD = N_GOOD + 1
            GOTO 200
          ENDIF
        ENDIF        
      
      RETURN
      END

C     ***********
C     *         *
C     * FIND_POWER_FUNCTION_DERIVATIVE  *
C     *         *
C     ***********

      SUBROUTINE FIND_POWER_FUNCTION_DERIVATIVE(N, K, X, F,
     O                                          FP)


C     + + + PURPOSE + + +
C     Find derivative at index K of the function tabulated in X and F
C     by computing the slope in log-log space and using the derivative
C     at that point of the simple power function defined by the 
C     log-log slope. 

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, K
      REAL*8 X(N), F(N), FP

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N - number of points in X and F
C     K - index of point in X at which derivative is needed
C     X - sequence of arguments
C     F - sequence of function values; all > 0
C     FP - computed derivative of power function at K

C     + + + LOCAL VARIABLES + + +
      REAL*8 B, H1, H2, HSUM

C***********************************************************************
C     The three points involved in computing the slope of a 
C     parabola are by:    X0, X0 + H1, and X0 + H1 + H2

      IF(K.EQ.1) THEN
C       WE ARE AT X0.
        H1 = LOG(X(K+1)/X(K))
        H2 = LOG(X(K+2)/X(K+1))
        HSUM = H1 + H2
        B =-(2.*H1 + H2)*LOG(F(K))/(H1*HSUM)
     A               + HSUM*LOG(F(K+1))/(H1*H2)
     B               - H1*LOG(F(K+2))/(HSUM*H2)

      ELSEIF(K.EQ.N) THEN
C       WE ARE AT X0+H1+H2
        H1 = LOG(X(K-1)/X(K-2))
        H2 = LOG(X(K)/X(K-1))
        HSUM = H1 + H2
        B = +H2*LOG(F(K-2))/(H1*HSUM)
     A                -HSUM*LOG(F(K-1))/(H1*H2)
     B                +(2.*H2 + H1)*LOG(F(K))/(HSUM*H2)
      ELSE
C       WE ARE AT X0+H1
        H1 = LOG(X(K)/X(K-1))
        H2 = LOG(X(K+1)/X(K))
        HSUM = H1 + H2
        B = -H2*LOG(F(K-1))/(H1*HSUM)
     A                +(H2 - H1)*LOG(F(K))/(H1*H2)
     B                + H1*LOG(F(K+1))/(HSUM*H2)
      ENDIF

      FP = B*F(K)/X(K)
      RETURN
      END
C     ***********
C     *         *
C     * FIND_PARABOLIC_DERIVATIVE  *
C     *         *
C     ***********

      SUBROUTINE FIND_PARABOLIC_DERIVATIVE(N, K, X, F,
     O                                          FP)

C     + + + PURPOSE + + +
C     Find derivative at index K of the function tabulated in X and F
C     by computing the slope of the parabola fitted to the point
C     at K and the two points nearest to the point at index K.

      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N, K
      REAL*8 X(N), F(N), FP

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N - number of points in X and F
C     K - index of point in X at which derivative is needed
C     X - sequence of arguments
C     F - sequence of function values; all > 0
C     FP - computed derivative at K

C     + + + LOCAL VARIABLES + + +

      REAL*8 H1, H2, HSUM
C***********************************************************************
C     The three points involved in computing the slope of a 
C     parabola are denoted by:  X0, X0 + H1, and X0 + H1 + H2

      IF(K.EQ.1) THEN
C       WE ARE AT X0.
        H1 = X(K+1) - X(K)
        H2 = X(K+2) - X(K+1)
        HSUM = H1 + H2
        FP =-(2.*H1 + H2)*F(K)/(H1*HSUM)
     A               + HSUM*F(K+1)/(H1*H2)
     B               - H1*F(K+2)/(HSUM*H2)

      ELSEIF(K.EQ.N) THEN
C       WE ARE AT X0+H1+H2
        H1 = X(K-1) - X(K-2)
        H2 = X(K) - X(K-1)
        HSUM = H1 + H2
        FP = +H2*F(K-2)/(H1*HSUM)
     A                -HSUM*F(K-1)/(H1*H2)
     B                +(2.*H2 + H1)*F(K)/(HSUM*H2)
      ELSE
C       WE ARE AT X0+H1
        H1 = X(K) - X(K-1)
        H2 = X(K+1) - X(K)
        HSUM = H1 + H2
        FP = -H2*F(K-1)/(H1*HSUM)
     A                +(H2 - H1)*F(K)/(H1*H2)
     B                + H1*F(K+1)/(HSUM*H2)
      ENDIF
      RETURN
      END

C
C
C

      SUBROUTINE FINDBRK
     I                   (STDOUT, NPNTS, XVEC, FVEC, RELATIVE_ERROR,
     I                    LEFT_OPTION, RIGHT_OPTION, 
     O                    N_GOOD, GOOD_POINTS, NO_GOOD, MAX_RERR)

C     + + + PURPOSE + + +
C     Find a series of breakpoints that will hopefully result
C     in a relative error of linear interpolation of no greater than
C     RELATIVE_ERROR.  The function is assumed to be defined
C     by the point set in (XVEC, FVEC) as fitted by a cubic 
C     spline.  The only requirement is that a cubic spline 
C     interpolate the function with accuracy much better than
C     the requested linear-interpolation error.

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

C     + + + DUMMY ARGUMENTS + + +

      INTEGER NPNTS, STDOUT, N_GOOD, NO_GOOD, LEFT_OPTION,
     A        RIGHT_OPTION
      REAL*8 XVEC(NPNTS), FVEC(NPNTS), GOOD_POINTS(NPNTS),
     A       MAX_RERR, RELATIVE_ERROR

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - standard output unit
C     NPNTS - number of values in XVEC and FVEC
C     XVEC - series of argument values
C     FVEC - function values at arguments in XVEC
C     RELATIVE_ERROR - target relative error for linear spline
C                      interpolation in the function tabulated
C                      in XVEC and FVEC
C     LEFT_OPTION - code for the slope at left end of range
C     RIGHT_OPTION - code for the slope at right end of range
C     N_GOOD - number of points in the good sequence
C     GOOD_POINTS - the sequence of good breakpoints
C     NO_GOOD - error flag. =1 if failure, =0 otherwise
C     MAX_RERR - absolute value of the maximum relative error
C               in linear interpolation using the GOOD_POINTS

C     + + + LOCAL VARIABLES + + +
      REAL*8 LEFT_SLOPE, RIGHT_SLOPE
C***********************************************************************
C     Define the end conditions.  Use derivative of simple power
C     function fitted to three points. 
      IF(LEFT_OPTION.EQ.1) THEN
C       Compute derivative using values as is.
        CALL FIND_PARABOLIC_DERIVATIVE(NPNTS, 1, XVEC, FVEC,
     O                                            LEFT_SLOPE)
      ELSE
C       Compute derivative using a fitted power function.
        CALL FIND_POWER_FUNCTION_DERIVATIVE(NPNTS, 1, XVEC, FVEC,
     O                                          LEFT_SLOPE)
      ENDIF
      IF(RIGHT_OPTION.EQ.1) THEN            
        CALL FIND_PARABOLIC_DERIVATIVE(NPNTS, NPNTS, XVEC, FVEC,
     O                                            RIGHT_SLOPE)
      ELSE
        CALL FIND_POWER_FUNCTION_DERIVATIVE(NPNTS, NPNTS, XVEC, FVEC,
     O                                            RIGHT_SLOPE)
      ENDIF

      N_GOOD = NPNTS
      CALL FIND_GOOD_POINTS(STDOUT, NPNTS, XVEC, FVEC, 
     I       RELATIVE_ERROR, LEFT_SLOPE, RIGHT_SLOPE,
     M       N_GOOD,
     O       GOOD_POINTS, NO_GOOD, MAX_RERR)

      RETURN
      END      
