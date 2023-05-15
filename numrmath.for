C     Contains stuff for some standard numerical methods.

C
C
C
      REAL FUNCTION   FDCNMN
     I                      (N, X, Y, XL, XR)
 
C     + + + PURPOSE + + +
C     Find a constrained minimum.  For X between XL and XR find
C     minimum value of Y.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      REAL X(N), XL, XR, Y(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     X      - Offsets of points on cross section boundary
C     Y      - vector of values defining minimum
C     XL     - Offset at left hand end of segment
C     XR     - Offset at right hand end of segment
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL YMIN
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
C***********************************************************************
      YMIN = 1.E30
      DO 100 I=1,N
        IF(X(I).GE.XL.AND.X(I).LE.XR) THEN
          YMIN = MIN(YMIN, Y(I))
        ENDIF
 100  CONTINUE
      FDCNMN = YMIN
      RETURN
      END
C
C
C
      SUBROUTINE   VLCHPP
     I                   (STDOUT, N, X, Y, LCODE, LVAL, RCODE, RVAL,
     O                    M, ADJLOC)
 
C     + + + PURPOSE + + +
C     Find a variation-limited cubic Hermite piecewise polynomial
C     that interpolates Y at each of the break points in X.  The
C     derivative of the function, y, with respect to x is stored
C     in M.   Such a function will have montone variation between
C     successive breakpoints and interpolation between breakpoints
C     does not introduce any extreme points.  That is, all extreme
C     values are already present in Y.
 
C     Method: Fit a cubic spline to the N points in X and Y and
C     then adjust the derivatives in M, as required to make the
C     variation monotone between break points.  The end conditions
C     for the spline are given in LCODE, LVAL for the left end
C     and RCODE, RVAL for the right end.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER LCODE, N, RCODE, STDOUT
      DOUBLE PRECISION LVAL, M(N), RVAL, X(N), Y(N)
      character*1 ADJLOC(N)

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     N      - number of points
C     X      - breakpoints for the spline
C     Y      - function values for spline
C     LCODE  - Code for left hand end condition
C     LVAL   - Value of the left hand end condition
C     RCODE  - Code for the right-hand end condition for the spline
C     RVAL   - Value of the right-hand end condition
C     M      - Vector of first derivatives
C     ADJLOC - values at which slope was adjusted to limit variation
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL SPLINE, VRLIM

c     Local
      integer nadj
C***********************************************************************
      CALL SPLINE
     I           (STDOUT, X, Y, N, LCODE, LVAL, RCODE, RVAL,
     O            M)
 
C     Check for monotone variation
 
      CALL VRLIM
     I          (N, X, Y,
     M           M,
     O           nadj, ADJLOC)
 
      RETURN
      END
C
C
C
      SUBROUTINE   IMTQL2
     I                   (N,
     M                    D, Z, E,
     O                    IERR)
 
C     + + + PURPOSE + + +

C     This subroutine finds the eigenvalues and first components of the
C     eigenvectors of a symmetric tridiagonal matrix by the implicit QL
C     method.

C     This subroutine is a translation of the algol procedure IMTQL2,
C     NUM. MATH. 12, 377-383(1968) by Martin and Wilkinson,
C     as modified in NUM. MATH. 15, 450(1970) by Dubrulle,
C     Handbook for Auto. Comp., Vol.II-Linear Algebra, 241-248(1971).
C     This is a modified version of the 'EISPACK' routine IMTQL2.
C     Called by GRULE and adapted from the same source. 
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER IERR, N
      DOUBLE PRECISION D(N), E(N), Z(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N      - the order of the matrix
C     D      - diagonal elements of the input matrix on input,
C            - eigenvalues in ascending order on output.  If an
C              error exit is made, the eigenvalues are correct but
C              unordered for indices 1, 2, ..., IERR-1.
C     Z      - contains the first row of the identity matrix on input.
C              contains the first components of the orthonormal 
C              eigenvectors of the symmetric tridiagonal matrix on
C              output.  If an error exit is  made, z contains the 
C              eigenvectors  associated with the stored eigenvalues.
C     E      - contains the subdiagonal elements of the input matrix
C              in its first N-1 positions.  E(N) is arbitrary. Values
C              are destroyed in processing. 
C     IERR   - is zero for normal return; is J if the J-th eigenvalue
C              has not been  determined after 30 iterations.
 
C     + + + SAVED VALUES + + +
      DOUBLE PRECISION MACHEP
      SAVE MACHEP
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, II, J, K, L, M, MML
      DOUBLE PRECISION B, C, F, G, P, R, S
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SIGN, SQRT
 
C     + + + DATA INITIALIZATIONS + + +
      DATA MACHEP/2.2D-16/
C***********************************************************************
      IERR = 0
      IF (N .EQ. 1) GO TO 1001

      E(N) = 0.0D0
      DO 240 L = 1, N
        J = 0
C     LOOK FOR SMALL SUB-DIAGONAL ELEMENT 
  105   DO 110 M = L, N
          IF (M .EQ. N) GO TO 120
          IF (ABS(E(M)) .LE. MACHEP * (ABS(D(M)) + ABS(D(M+1))))
     X       GO TO 120
  110   CONTINUE

  120   P = D(L)
        IF (M .EQ. L) GO TO 240
        IF (J .EQ. 30) GO TO 1000
        J = J + 1
C        FORM SHIFT 
        G = (D(L+1) - P) / (2.0D0 * E(L))
        R = SQRT(G*G+1.0D0)
        G = D(M) - P + E(L) / (G + SIGN(R, G))
        S = 1.0D0
        C = 1.0D0
        P = 0.0D0
        MML = M - L

C       FOR I=M-1 STEP -1 UNTIL L DO 
        DO 200 II = 1, MML
          I = M - II
          F = S * E(I)
          B = C * E(I)
          IF (ABS(F) .LT. ABS(G)) GO TO 150
          C = G / F
          R = SQRT(C*C+1.0D0)
          E(I+1) = F * R
          S = 1.0D0 / R
          C = C * S
          GO TO 160
  150     S = F / G
          R = SQRT(S*S+1.0D0)
          E(I+1) = G * R
          C = 1.0D0 / R
          S = S * C
  160     G = D(I+1) - P
          R = (D(I) - G) * S + 2.0D0 * C * B
          P = S * R
          D(I+1) = G + P
          G = C * R - B
C         FORM FIRST COMPONENT OF VECTOR
          F = Z(I+1)
          Z(I+1) = S * Z(I) + C * F
  200     Z(I) = C * Z(I) - S * F

         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE

C     ORDER EIGENVALUES AND EIGENVECTORS 
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)

         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE

         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
  300 CONTINUE

      GO TO 1001
C     SET ERROR -- NO CONVERGENCE TO AN
C     EIGENVALUE AFTER 30 ITERATIONS 
 1000 IERR = L
 1001 RETURN
      END
C
C
C
      SUBROUTINE   GRULE
     I                  (N,
     O                   X, W)
 
C     + + + PURPOSE + + +
C     Compute the N  abscissas X(I) and the corresponding weights
C     W(I), of the n-point Gauss-Legendre integration rule, normalizing
C     to the interval [-1,1]. The abscissas appear in asending order.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      DOUBLE PRECISION W(N), X(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N      - Number of points in the Gauss rule
C     X      - Abscissas for Gaussian quadrature
C     W      - Weights for Gaussian quadrature
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IERR, NM1
      DOUBLE PRECISION ABI, B(256), MUZERO
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE, SQRT
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL IMTQL2
C***********************************************************************
 
C     Method
 
C           The coefficients of the three-term recurrence relation
C        for the corresponding set of orthogonal polynomials are
C        used to form a symmetric tridiagonal matrix, whose
C        eigenvalues (determined by the implicit QL-method with
C        shifts) are just the desired nodes.  The first components of
C        the orthonormalized eigenvectors, when properly scaled,
C        yield the weights.  This technique is much faster than using a
C        root-finder to locate the zeroes of the orthogonal polynomial.
C        for further details, see Ref. 1. 
 
C     References
 
C        1.  Golub, G. H., and Welsch, J. H., "Calculation of Gaussian
C            Quadrature Rules," Mathematics of Computation 23 (April,
C            1969), pp. 221-230.
 
C       Adapted from subroutine GAUSSQ from a public domain library
C       discontinued years ago. 
 
      NM1 = N - 1
      MUZERO = 2.0D0
      DO 11 I = 1, NM1
        X(I) = 0.0D0
        ABI = DBLE(I)
        B(I) = ABI/SQRT(4*ABI*ABI - 1.0D0)
   11 CONTINUE
      X(N) = 0.0D0
 
C     THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
C     THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
C     B THE OFF-DIAGONAL ELEMENTS.
 
      W(1) = 1.0D0
      DO 105 I = 2, N
        W(I) = 0.0D0
  105 CONTINUE
C
      CALL IMTQL2
     I           (N,
     M            X, W, B,
     O            IERR)
      DO 110 I = 1, N
        W(I) = MUZERO * W(I) * W(I)
  110 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE   SQUANK
     I                   (A, B, EPS, FUN,
     O                    ERR, ERR5, NO, RESULT)
 
C     + + + PURPOSE + + +
C     SQUANK -- Simpson quadrature used adaptively, noise killed.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NO
      DOUBLE PRECISION A, B, EPS, ERR, ERR5, RESULT
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      DOUBLE PRECISION FUN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     A      - Lower limit for integration
C     B      - Upper limit for integration
C     EPS    - Convergence tolerance
C     FUN    - Integrand function
C     ERR    - Error indicator
C     ERR5   - Error indicator
C     NO     - Number of function evaluations
C     RESULT - Value of the integral
 
C     + + + LOCAL VARIABLES + + +
      INTEGER LEV, LEVTAG, NIM, NOM, NUM
      DOUBLE PRECISION ADIFF, ADIFF1, CEPS, CEPSF, CEPST, CRIT, DIFF,
     A                 DIM, DUM, EFACT, EPMACH, EST, EST1, EST2,
     B                 ESTST(30), FACERR, FX1, FX2, FX3, FX3ST(30), FX4,
     C                 FX5, FX5ST(30), PREDIF(30), QCEPS, SIM, SUM, T0,
     D                 T1, T2, T3, THIRD, X1, X2, X3, X3ST(30), X4, X5,
     E                 X5ST(30), XZERO
 
C     + + + INTRINSICS + + +
      INTRINSIC DABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FUN
C***********************************************************************

C       Adapted from subroutine SQUANK from public domain library
C       which once existed at Stanford University. 
C       This subroutine is basically the same as  the  one  published  
C       by J. N. Lyness  as  ACM  Algorithm  379,  Comm.  A.C.M.  13, 
C       1970, pp. 260-262.

C          ****    STAGE ONE    ****
C      INITIALISE QUANTITIES REQUIRED FOR CENTRAL CALCULATION (STAGE 3).
 
      SUM = 0.0
      SIM = 0.0
      DUM = 0.0
      DIM = 0.0
      CEPSF = 180.0*EPS/(B - A)
      CEPS = CEPSF
      ADIFF  = 0.0
      LEVTAG = -1
      FACERR = 1.0
      XZERO = A
      EFACT = 0.0
      NIM = 1
      LEV = 0
C     FIRST INTERVAL
      X1 = A
      X5 = B
      X3 = 0.5*(A+B)
      FX1= FUN(X1)
      FX3= FUN(X3)
      FX5= FUN(X5)
      NO = 3
      EST = FX1 + FX5 + 4.0*FX3
 
C          ****    STAGE TWO    ****
C      SET A STARTING VALUE FOR TOLERANCE IN CASE THAT CEPSF = 0.0
 
      IF(CEPSF) 295,205,295
  205 LEVTAG = 0
      EPMACH = 1.
 208  EPMACH = .5*EPMACH
      EST1 = EPMACH+1.
      IF (EST1.GT.1.) GO TO 208
      FACERR = 15.0
      CEPS = EPMACH*DABS (FX1)
      IF(FX1) 295,210,295
  210 CEPS = EPMACH*DABS (FX3)
      LEVTAG = 3
      IF(FX3) 295,215,295
  215 CEPS = EPMACH*DABS (FX5)
      IF(FX5) 295,220,295
  220 CEPS = EPMACH
  295 QCEPS = 0.25*CEPS
C     INITIALISING COMPLETE
 
C          ****    STAGE THREE    ****
C      CENTRAL CALCULATION.
C                   REQUIRES X1,X3,X5,FX1,FX3,FX5,EST,ADIFF.
 
  300 CONTINUE
      X2 = 0.5*(X1 + X3)
      X4 = 0.5*(X3 + X5)
      FX2= FUN(X2)
      FX4= FUN(X4)
      NO = NO + 2
      EST1 = FX1 + 4.0*FX2 + FX3
      EST2 = FX3 + 4.0*FX4 + FX5
      ADIFF1 = ADIFF
      DIFF = EST + EST - EST1 - EST2
      IF(LEV - 30) 305,800,800
  305 ADIFF= DABS (DIFF)
      CRIT = ADIFF - CEPS
      IF(CRIT)  700,700,400
C     END OF CENTRAL LOOP
C      NEXT STAGE IS STAGE FOUR  IN CASE OF NO NATURAL CONVERGENCE
C      NEXT STAGE IS STAGE SEVEN IN CASE OF    NATURAL CONVERGENCE
 
C          ****    STAGE FOUR    ****
C      NO NATURAL CONVERGENCE.  A COMPLEX SEQUENCE OF INSTRUCTIONS
C      FOLLOWS WHICH ASSIGNS CONVERGENCE AND / OR ALTERS TOLERANCE
C      LEVEL IN UPWARD DIRECTION IF THERE ARE INDICATIONS OF ROUND OFF
C      ERROR.
 
  400 CONTINUE
      IF(ADIFF1 - ADIFF)  410, 410, 500
C      IN A NORMAL RUN WITH NO ROUND OFF ERROR PROBLEM,  ADIFF1  IS
C      GREATER THAN ADIFF  AND THE REST OF STAGE FOUR IS OMITTED.
  410 IF(LEV - 5) 500,415,415
  415 EFACT = EFACT + CEPS *(X1 - XZERO)*FACERR
      XZERO = X1
      FACERR = 15.0
C      THE REST OF STAGE FOUR DEALS WITH UPWARD ADJUSTMENT OF
C      TOLERANCE (CEPS) BECAUSE OF SUSPECTED ROUND OFF ERROR TROUBLE.
      IF(ADIFF-2.0*CEPS)  420,420,425
C      SMALL JUMP IN CEPS. ASSIGN CONVERGENCE
  420 CEPS =     ADIFF
      LEVTAG = 0
      GO TO 780
  425 IF( ADIFF1 - ADIFF)  435,430,435
C      LARGE JUMP IN CEPS
  430 CEPS = ADIFF
      LEVTAG = 3
      GO TO 445
C      FACTOR TWO  JUMP IN CEPS
  435 CEPS = 2.0*CEPS
      IF(LEVTAG - 3) 440,445,445
  440 LEVTAG = 2
  445 QCEPS = 0.25*CEPS
 
C          ****    STAGE FIVE    ****
C       NO ACTUAL CONVERGENCE.
C     STORE RIGHT HAND ELEMENTS
 
  500 CONTINUE
      NIM = 2*NIM
      LEV = LEV + 1
      ESTST(LEV) = EST2
      X3ST(LEV)= X4
      X5ST(LEV)= X5
      FX3ST(LEV)=FX4
      FX5ST(LEV)=FX5
      PREDIF(LEV) = ADIFF
 
C          ****    STAGE SIX    ****
C      SET UP QUANTITIES FOR CENTRAL CALCULATION.
 
C     READY TO GO AHEAD AT LEVEL LOWER WITH LEFT HAND ELEMENTS
C     X1 AND FX1 ARE THE SAME AS BEFORE
      X5 = X3
      X3 = X2
      FX5 = FX3
      FX3 = FX2
      EST = EST1
      GO TO 300
 
C          ****    STAGE SEVEN    ****
C      NATURAL CONVERGENCE IN PREVIOUS INTERVAL. THE  FOLLOWING  COMPLEX
C      SEQUENCE  CHECKS  PRIMARILY THAT TOLERANCE LEVEL IS NOT TOO HIGH.
C      UNDER CERTAIN CIRCUMSTANCES  NON CONVERGENCE IS ASSIGNED AND / OR
C      TOLERANCE LEVEL IS RE-SET.
 
  700 CONTINUE
C      CHECK THAT IT WAS NOT LEVEL ZERO INTERVAL. IF SO ASSIGN
C      NON CONVERGENCE
      IF( LEV ) 400,400, 705
 
C      LEVTAG =-1   CEPS = CEPSF, ITS ORIGINAL VALUE.
C      LEVTAG = 0   CEPS IS GREATER THAN CEPSF. REGULAR SITUATION.
C      LEVTAG = 2   CEPS IS GREATER THAN CEPSF. CEPS PREVIOUSLY
C                           ASKED FOR A BIG JUMP, BUT DID NOT GET ONE.
C      LEVTAG = 3   CEPS IS GREATER THAN CEPSF. CEPS PREVIOUSLY HAD
C                           A BIG JUMP.
  705 IF(LEVTAG) 800,710,710
C      IN A NORMAL RUN WITH NO ROUND OFF ERROR PROBLEM,  LEVTAG = -1
C      AND THE REST OF STAGE SEVEN IS OMITTED.
  710 CEPST = 15.0*CEPS
C      CEPST HERE IS FACERR*CURRENT VALUE OF CEPS
      IF(CRIT) 715,800,800
  715 IF(LEVTAG - 2)  720,740,750
C      LEVTAG = 0
  720 IF(ADIFF) 800,800,725
  725 IF(ADIFF - QCEPS) 730,800,800
  730 IF(ADIFF - CEPSF) 770,770,735
  735 LEVTAG = 0
      CEPS = ADIFF
      EFACT = EFACT + CEPST*(X1 - XZERO)
      XZERO = X1
      GO TO 445
C      LEVTAG = 2
  740 LEVTAG = 0
      IF(ADIFF) 765,765,725
C      LEVTAG = 3
  750 LEVTAG = 0
      IF(ADIFF) 775,775,730
  765 CEPS = ADIFF1
      GO TO 775
  770 LEVTAG = -1
      FACERR = 1.0
      CEPS = CEPSF
  775 EFACT = EFACT + CEPST*(X1 - XZERO)
      XZERO = X1
  780 CONTINUE
      QCEPS = 0.25*CEPS
C          ****    STAGE EIGHT    ****
C      ACTUAL CONVERGENCE IN PREVIOUS INTERVAL.  INCREMENTS ADDED INTO
C      RUNNING SUMS
 
C     ADD INTO SUM AND SIM
  800 CONTINUE
C  THE FOLLOWING CODING AVOIDS LOSS OF PRECISION IN ACCUMULATING THE
C  SUM.  IT CORRESPONDS TO
C     SUM = SUM + T1
 
      T1 = (EST1+EST2)*(X5-X1)
      T0 = SUM+T1
      T2 = DABS(T1)
      T3 = DABS(SUM)
      IF (T2.GT.T3) GO TO 801
      DUM = DUM + T1 - (T0-SUM)
      GO TO 802
  801 DUM = DUM + SUM - (T0-T1)
  802 SUM = T0
      IF (LEVTAG) 805,810,810
C      WE ADD INTO SIM ONLY IF WE ARE CLEAR OF ROUND OFF LEVEL.
  805 T1 = DIFF*(X5-X1)
      T0 = SUM+T1
      T2 = DABS(T1)
      T3 = DABS(SIM)
      IF (T2.GT.T3) GO TO 806
      DIM = DIM + T1 - (T0-SIM)
      GO TO 807
  806 DIM = DIM + SIM - (T0-T1)
  807 SIM = T0
  810 CONTINUE
 
C          ****    STAGE NINE    ****
C      SORT OUT WHICH LEVEL TO GO TO.  THIS INVOLVES NIM NUMBERING
C      SYSTEM DESCRIBED BEFORE STAGE ONE.
 
  905 NUM = NIM/2
      NOM = NIM - 2*NUM
      IF(NOM) 910,915,910
  910 NIM = NUM
      LEV = LEV - 1
      GO TO 905
  915 NIM = NIM + 1
C     NEW LEVEL IS SET. IF LEV=0 WE HAVE FINISHED
      IF( LEV ) 1100,1100,1000
 
C          ****    STAGE TEN    ****
C      SET UP QUANTITIES FOR CENTRAL CALCULATION.
 
 1000 CONTINUE
      X1 = X5
      FX1= FX5
      X3 = X3ST(LEV)
      X5 = X5ST(LEV)
      FX3= FX3ST(LEV)
      FX5= FX5ST(LEV)
      EST= ESTST(LEV)
      ADIFF = PREDIF(LEV)
      GO TO 300
 
C          ****    STAGE ELEVEN    ****
C      CALCULATION NOW COMPLETE. FINALISE.
 
 1100 CONTINUE
      EFACT = EFACT + CEPS *(B- XZERO)*FACERR
      ERR     = EFACT/180.0
      THIRD= (SUM+DUM)/12.0
      ERR5 =-(SIM+DIM)/180.0
      RESULT  = THIRD + ERR5
      RETURN
      END
C
C
C
      SUBROUTINE   QFUN
     I                 (A, B, FUN,
     O                  QERR, RESULT)
 
C     + + + PURPOSE + + +
C     Integrate the double precision function, FUN, between the limits
C     of A and B.  Standard error limits are obtained from EPSCOM.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL A, B, QERR
      DOUBLE PRECISION RESULT
 
C     + + + DUMMY ARGUMENT FUNCTIONS + + +
      DOUBLE PRECISION FUN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     A      - Lower limit for integration
C     B      - Upper limit for integration
C     FUN    - Integrand function
C     QERR   - Quadrature error
C     RESULT - Value of the integral
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'epscom.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER LIMFLG, N
      DOUBLE PRECISION ABSERR, ADBL, BDBL, EPS, ERR, ERR5, TP
 
C     + + + INTRINSICS + + +
      INTRINSIC DBLE
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FUN, SQUANK
C     *********************************FORMATS**************************
C50    FORMAT(' QFUN: ERR=',1PE12.4,' ERR5=',1PE12.4,' N=',I7)
C***********************************************************************
C     SET THE DOUBLE PRECISION LIMITS
 
      ADBL = DBLE(A)
      BDBL = DBLE(B)
 
      IF(ADBL.GT.BDBL) THEN
C       SWITCH LIMITS BECAUSE SOME QUADRATURE ROUTINES FAIL TO FUNCTION
C       IF A > B.
 
        LIMFLG = 1
        TP = BDBL
        BDBL = ADBL
        ADBL = TP
      ELSE
        LIMFLG = 0
      ENDIF
 
C     SET THE ERROR TOLERANCES
      EPS = EPSINT
 
 
      IF(ADBL.NE.BDBL) THEN
        IF(INTHOW.EQ.1) THEN
C         USE ADAPTIVE SIMPSON'S RULE
          CALL SQUANK
     I               (ADBL, BDBL, EPS, FUN,
     O                ERR, ERR5, N, RESULT)
          QERR = ERR
C          WRITE(STD6,50) ERR, ERR5, N
        ELSE
          CALL SQUANK
     I               (ADBL, BDBL, EPS, FUN,
     O                ERR, ERR5, N, RESULT)
          QERR = ERR
C          WRITE(STD6,50) ERR, ERR5, N
        ENDIF
      ELSE
        RESULT = 0.0D0
      ENDIF
      IF(LIMFLG.EQ.1) THEN
        RESULT = -RESULT
      ENDIF
 
      RETURN
      END
