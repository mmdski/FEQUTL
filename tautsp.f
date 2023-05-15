C
C
C
      SUBROUTINE TAUTSP ( TAU, GTAU, NTAU, GAMMA, NBIG,
     O                    S,
     O                    BREAK, COEF, L, K, IFLAG )
C  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR    
CONSTRUCTS CUBIC SPLINE INTERPOLANT TO GIVEN DATA
C         TAU(I), GTAU(I), I=1,...,NTAU.
C  IF  GAMMA .GT. 0., ADDITIONAL KNOTS ARE INTRODUCED WHERE NEEDED TO
C  MAKE THE INTERPOLANT MORE FLEXIBLE LOCALLY. THIS AVOIDS EXTRANEOUS
C  INFLECTION POINTS TYPICAL OF CUBIC SPLINE INTERPOLATION AT KNOTS TO
C  RAPIDLY CHANGING DATA.
C
C  PARAMETERS
C            INPUT
C  TAU      SEQUENCE OF DATA POINTS. MUST BE STRICTLY INCREASING.
C  GTAU     CORRESPONDING SEQUENCE OF FUNCTION VALUES.
C  NTAU     NUMBER OF DATA POINTS. MUST BE AT LEAST  4 .
C  GAMMA    INDICATES WHETHER ADDITIONAL FLEXIBILITY IS DESIRED.
C          = 0., NO ADDITIONAL KNOTS
C          IN (0.,3.), UNDER CERTAIN CONDITIONS ON THE GIVEN DATA AT
C                POINTS I-1, I, I+1, AND I+2, A KNOT IS ADDED IN THE
C                I-TH INTERVAL, I=2,...,NTAU-2. SEE DESCRIPTION OF METH-
C                OD BELOW. THE INTERPOLANT GETS ROUNDED WITH INCREASING
C                GAMMA. A VALUE OF  2.5  FOR GAMMA IS TYPICAL.
C          IN (3.,6.), SAME , EXCEPT THAT KNOTS MIGHT ALSO BE ADDED IN
C                INTERVALS IN WHICH AN INFLECTION POINT WOULD BE PERMIT-
C                TED.  A VALUE OF  5.5  FOR GAMMA IS TYPICAL.
C            OUTPUT
C  BREAK, COEF, L, K  GIVE THE PP-REPRESENTATION OF THE INTERPOLANT.
C          SPECIFICALLY, FOR BREAK(I) .LE. X .LE. BREAK(I+1), THE
C        INTERPOLANT HAS THE FORM
C  F(X) = COEF(1,I) +DX(COEF(2,I) +(DX/2)(COEF(3,I) +(DX/3)COEF(4,I)))
C        WITH  DX = X - BREAK(I) AND I=1,...,L .
C  IFLAG   = 1, OK
C          = 2, INPUT WAS INCORRECT. A PRINTOUT SPECIFYING THE MISTAKE
C            WAS MADE.
C            WORKSPACE
C  S     IS REQUIRED, OF SIZE (NTAU,6). THE INDIVIDUAL COLUMNS OF THIS
C        ARRAY CONTAIN THE FOLLOWING QUANTITIES MENTIONED IN THE WRITE-
C        UP AND BELOW.
C     S(.,1) = DTAU = TAU(.+1) - TAU
C     S(.,2) = DIAG = DIAGONAL IN LINEAR SYSTEM
C     S(.,3) = U = UPPER DIAGONAL IN LINEAR SYSTEM
C     S(.,4) = R = RIGHT SIDE FOR LINEAR SYSTEM (INITIALLY)
C            = FSECND = SOLUTION OF LINEAR SYSTEM , NAMELY THE SECOND
C                       DERIVATIVES OF INTERPOLANT AT  TAU
C     S(.,5) = Z = INDICATOR OF ADDITIONAL KNOTS
C     S(.,6) = 1/HSECND(1,X) WITH X = Z OR = 1-Z. SEE BELOW.
C
C  ------  M E T H O D  ------
C  ON THE I-TH INTERVAL, (TAU(I), TAU(I+1)), THE INTERPOLANT IS OF THE
C  FORM
C  (*)  F(U(X)) = A + B*U + C*H(U,Z) + D*H(1-U,1-Z) ,
C  WITH  U = U(X) = (X - TAU(I))/DTAU(I). HERE,
C       Z = Z(I) = ADDG(I+1)/(ADDG(I) + ADDG(I+1))
C  (= .5, IN CASE THE DENOMINATOR VANISHES). WITH
C       ADDG(J) = ABS(DDG(J)), DDG(J) = DG(J+1) - DG(J),
C       DG(J) = DIVDIF(J) = (GTAU(J+1) - GTAU(J))/DTAU(J)
C  AND
C       H(U,Z) = ALPHA*U**3 + (1 - ALPHA)*(MAX(((U-ZETA)/(1-ZETA)),0)**3
C  WITH
C       ALPHA(Z) = (1-GAMMA/3)/ZETA
C       ZETA(Z) = 1 - GAMMA*MIN((1 - Z), 1/3)
C  THUS, FOR 1/3 .LE. Z .LE. 2/3,  F  IS JUST A CUBIC POLYNOMIAL ON
C  THE INTERVAL I. OTHERWISE, IT HAS ONE ADDITIONAL KNOT, AT
C         TAU(I) + ZETA*DTAU(I) .
C  AS  Z  APPROACHES  1, H(.,Z) HAS AN INCREASINGLY SHARP BEND  NEAR 1,
C  THUS ALLOWING  F  TO TURN RAPIDLY NEAR THE ADDITIONAL KNOT.
C     IN TERMS OF F(J) = GTAU(J) AND
C       FSECND(J) = 2.DERIVATIVE OF  F  AT  TAU(J),
C  THE COEFFICIENTS FOR (*) ARE GIVEN AS
C       A = F(I) - D
C       B = (F(I+1) - F(I)) - (C - D)
C       C = FSECND(I+1)*DTAU(I)**2/HSECND(1,Z)
C       D = FSECND(I)*DTAU(I)**2/HSECND(1,1-Z)
C  HENCE CAN BE COMPUTED ONCE FSECND(I),I=1,...,NTAU, IS FIXED.
C   F  IS AUTOMATICALLY CONTINUOUS AND HAS A CONTINUOUS SECOND DERIVAT-
C  IVE (EXCEPT WHEN Z = 0 OR 1 FOR SOME I). WE DETERMINE FSCND(.) FROM
C  THE REQUIREMENT THAT ALSO THE FIRST DERIVATIVE OF  F  BE CONTIN-
C  UOUS. IN ADDITION, WE REQUIRE THAT THE THIRD DERIVATIVE BE CONTINUOUS
C  ACROSS  TAU(2) AND ACROSS  TAU(NTAU-1) . THIS LEADS TO A STRICTLY
C  DIAGONALLY DOMINANT TRIDIAGONAL LINEAR SYSTEM FOR THE FSECND(I)
C  WHICH WE SOLVE BY GAUSS ELIMINATION WITHOUT PIVOTING.

      IMPLICIT NONE
      INTEGER IFLAG,K,L,NTAU,   I,METHOD,NTAUM1
      INTEGER NBIG
      REAL*8 BREAK(NBIG),COEF(4,NBIG),GAMMA,GTAU(NTAU),S(NBIG,6),
     A      TAU(NTAU)
     *    ,ALPHA,C,D,DEL,DENOM,DIVDIF,ENTRY,ENTRY3,FACTOR,FACTR2,GAM
     *    ,ONEMG3,ONEMZT,RATIO,SIXTH,TEMP,X,Z,ZETA,ZT2
      REAL*8 ALPH
C***********************************************************************
      ALPH(X) = MIN(1.,ONEMG3/X)
C
C  THERE MUST BE AT LEAST  4  INTERPOLATION POINTS.
      IF (NTAU .GE. 4)                  GO TO 5
      PRINT 600,NTAU
  600 FORMAT(8H NTAU = ,I4,20H. SHOULD BE .GE. 4 .)
                                        GO TO 999
C
CONSTRUCT DELTA TAU AND FIRST AND SECOND (DIVIDED) DIFFERENCES OF DATA
C
    5 NTAUM1 = NTAU - 1
      DO 6 I=1,NTAUM1
         S(I,1) = TAU(I+1) - TAU(I)
         IF (S(I,1) .GT. 0.)            GO TO 6
         PRINT 610,I,TAU(I),TAU(I+1)
  610    FORMAT(7H POINT ,I3,13H AND THE NEXT,2E15.6,15H ARE DISORDERED)
                                        GO TO 999
    6    S(I+1,4) = (GTAU(I+1)-GTAU(I))/S(I,1)
      DO 7 I=2,NTAUM1
    7    S(I,4) = S(I+1,4) - S(I,4)
C
CONSTRUCT SYSTEM OF EQUATIONS FOR SECOND DERIVATIVES AT  TAU. AT EACH
C  INTERIOR DATA POINT, THERE IS ONE CONTINUITY EQUATION, AT THE FIRST
C  AND THE LAST INTERIOR DATA POINT THERE IS AN ADDITIONAL ONE FOR A
C  TOTAL OF  NTAU  EQUATIONS IN  NTAU  UNKNOWNS.
C
      I = 2
      S(2,2) = S(1,1)/3.
      SIXTH = 1./6.
      METHOD = 2
      GAM = GAMMA
      IF (GAM .LE. 0.)   METHOD = 1
      IF ( GAM .LE. 3.)                 GO TO 9
      METHOD = 3
      GAM = GAM - 3.
    9 ONEMG3 = 1. - GAM/3.
C                 ------ LOOP OVER I ------
   10 CONTINUE
C          CONSTRUCT Z(I) AND ZETA(I)
      Z = .5
                                        GO TO (19,11,12),METHOD
   11 IF (S(I,4)*S(I+1,4) .LT. 0.)      GO TO 19
   12 TEMP = ABS(S(I+1,4))
      DENOM = ABS(S(I,4)) + TEMP
      IF (DENOM .EQ. 0.)                GO TO 19
      Z = TEMP/DENOM
      IF (ABS(Z - .5) .LE. SIXTH)  Z = .5
   19 S(I,5) = Z
C   ******SET UP PART OF THE I-TH EQUATION WHICH DEPENDS ON
C         THE I-TH INTERVAL
      IF (Z - .5)                       21,22,23
   21 ZETA = GAM*Z
      ONEMZT = 1. - ZETA
      ZT2 = ZETA**2
      ALPHA = ALPH(ONEMZT)
      FACTOR = ZETA/(ALPHA*(ZT2-1.) + 1.)
      S(I,6) = ZETA*FACTOR/6.
      S(I,2) = S(I,2) + S(I,1)*((1.-ALPHA*ONEMZT)*FACTOR/2. - S(I,6))
C     IF Z = 0 AND THE PREVIOUS Z = 1, THEN D(I) = 0. SINCE THEN
C     ALSO U(I-1) = L(I+1) = 0, ITS VALUE DOES NOT MATTER. RESET
C     D(I) = 1 TO INSURE NONZERO PIVOT IN ELIMINATION.
      IF (S(I,2) .LE. 0.) S(I,2) = 1.
      S(I,3) = S(I,1)/6.
                                        GO TO 25
   22 S(I,2) = S(I,2) + S(I,1)/3.
      S(I,3) = S(I,1)/6.
                                        GO TO 25
   23 ONEMZT = GAM*(1. - Z)
      ZETA = 1. - ONEMZT
      ALPHA = ALPH(ZETA)
      FACTOR = ONEMZT/(1. - ALPHA*ZETA*(1.+ONEMZT))
      S(I,6) = ONEMZT*FACTOR/6.
      S(I,2) = S(I,2) + S(I,1)/3.
      S(I,3) = S(I,6)*S(I,1)
   25 IF (I .GT. 2)                     GO TO 30
      S(1,5) = .5
C  ******THE FIRST TWO EQUATIONS ENFORCE CONTINUITY OF THE FIRST AND OF
C        THE THIRD DERIVATIVE ACROSS TAU(2).
      S(1,2) = S(1,1)/6.
      S(1,3) = S(2,2)
      ENTRY3 = S(2,3)
      IF (Z - .5)                       26,27,28
   26 FACTR2 = ZETA*(ALPHA*(ZT2-1.) + 1.)/(ALPHA*(ZETA*ZT2-1.)+1.)
      RATIO = FACTR2*S(2,1)/S(1,2)
      S(2,2) = FACTR2*S(2,1) + S(1,1)
      S(2,3) = -FACTR2*S(1,1)
                                        GO TO 29
   27 RATIO = S(2,1)/S(1,2)
      S(2,2) = S(2,1) + S(1,1)
      S(2,3) = -S(1,1)
                                        GO TO 29
   28 RATIO = S(2,1)/S(1,2)
      S(2,2) = S(2,1) + S(1,1)
      S(2,3) = -S(1,1)*6.*ALPHA*S(2,6)
C       AT THIS POINT, THE FIRST TWO EQUATIONS READ
C              DIAG(1)*X1 + U(1)*X2 + ENTRY3*X3 = R(2)
C       -RATIO*DIAG(1)*X1 + DIAG(2)*X2 + U(2)*X3 = 0.
C       ELIMINATE FIRST UNKNOWN FROM SECOND EQUATION
   29 S(2,2) = RATIO*S(1,3) + S(2,2)
      S(2,3) = RATIO*ENTRY3 + S(2,3)
      S(1,4) = S(2,4)
      S(2,4) = RATIO*S(1,4)
                                        GO TO 35
   30 CONTINUE
C  ******THE I-TH EQUATION ENFORCES CONTINUITY OF THE FIRST DERIVATIVE
C        ACROSS TAU(I). IT HAS BEEN SET UP IN STATEMENTS 35 UP TO 40
C        AND 21 UP TO 25 AND READS NOW
C         -RATIO*DIAG(I-1)*XI-1 + DIAG(I)*XI + U(I)*XI+1 = R(I) .
C        ELIMINATE (I-1)ST UNKNOWN FROM THIS EQUATION
      S(I,2) = RATIO*S(I-1,3) + S(I,2)
      S(I,4) = RATIO*S(I-1,4) + S(I,4)
C
C  ******SET UP THE PART OF THE NEXT EQUATION WHICH DEPENDS ON THE
C        I-TH INTERVAL.
   35 IF (Z - .5)                       36,37,38
   36 RATIO = -S(I,6)*S(I,1)/S(I,2)
      S(I+1,2) = S(I,1)/3.
                                        GO TO 40
   37 RATIO = -(S(I,1)/6.)/S(I,2)
      S(I+1,2) = S(I,1)/3.
                                        GO TO 40
   38 RATIO = -(S(I,1)/6.)/S(I,2)
      S(I+1,2) = S(I,1)*((1.-ZETA*ALPHA)*FACTOR/2. - S(I,6))
C         ------  END OF I LOOP ------
   40 I = I+1
      IF (I .LT. NTAUM1)                GO TO 10
      S(I,5) = .5
C
C        ------  LAST TWO EQUATIONS  ------
C  THE LAST TWO EQUATIONS ENFORCE CONTINUITY OF THIRD DERIVATIVE AND
C  OF FIRST DERIVATIVE ACROSS  TAU(NTAU-1).
      ENTRY = RATIO*S(I-1,3) + S(I,2) + S(I,1)/3.
      S(I+1,2) = S(I,1)/6.
      S(I+1,4) = RATIO*S(I-1,4) + S(I,4)
      IF (Z - .5)                       41,42,43
   41 RATIO = S(I,1)*6.*S(I-1,6)*ALPHA/S(I-1,2)
      S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
      S(I,3) = -S(I-1,1)
                                        GO TO 45
   42 RATIO = S(I,1)/S(I-1,2)
      S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
      S(I,3) = -S(I-1,1)
                                        GO TO 45
   43 FACTR2 = ONEMZT*(ALPHA*(ONEMZT**2-1.)+1.)/
     *               (ALPHA*(ONEMZT**3-1.)+1.)
      RATIO = FACTR2*S(I,1)/S(I-1,2)
      S(I,2) = RATIO*S(I-1,3) + FACTR2*S(I-1,1) + S(I,1)
      S(I,3) = -FACTR2*S(I-1,1)
C     AT THIS POINT, THE LAST TWO EQUATIONS READ
C             DIAG(I)*XI + U(I)*XI+1 = R(I)
C      -RATIO*DIAG(I)*XI + DIAG(I+1)*XI+1 = R(I+1)
C     ELIMINATE XI FROM LAST EQUATION
   45 S(I,4) = RATIO*S(I-1,4)
      RATIO = -ENTRY/S(I,2)
      S(I+1,2) = RATIO*S(I,3) + S(I+1,2)
      S(I+1,4) = RATIO*S(I,4) + S(I+1,4)
C
C        ------ BACK SUBSTITUTION ------
C
      S(NTAU,4) = S(NTAU,4)/S(NTAU,2)
   50    S(I,4) = (S(I,4) - S(I,3)*S(I+1,4))/S(I,2)
         I = I - 1
         IF (I .GT. 1)                  GO TO 50
      S(1,4) = (S(1,4)-S(1,3)*S(2,4)-ENTRY3*S(3,4))/S(1,2)
C
C        ------ CONSTRUCT POLYNOMIAL PIECES ------
C
      BREAK(1) = TAU(1)
      L = 1
      DO 70 I=1,NTAUM1
         COEF(1,L) = GTAU(I)
         COEF(3,L) = S(I,4)
         DIVDIF = (GTAU(I+1)-GTAU(I))/S(I,1)
         Z = S(I,5)
         IF (Z - .5)                    61,62,63
   61    IF (Z .EQ. 0.)                 GO TO 65
         ZETA = GAM*Z
         ONEMZT = 1. - ZETA
         C = S(I+1,4)/6.
         D = S(I,4)*S(I,6)
         L = L + 1
         DEL = ZETA*S(I,1)
         BREAK(L) = TAU(I) + DEL
         ZT2 = ZETA**2
         ALPHA = ALPH(ONEMZT)
         FACTOR = ONEMZT**2*ALPHA
         COEF(1,L) = GTAU(I) + DIVDIF*DEL
     *             + S(I,1)**2*(D*ONEMZT*(FACTOR-1.)+C*ZETA*(ZT2-1.))
         COEF(2,L) = DIVDIF + S(I,1)*(D*(1.-3.*FACTOR)+C*(3.*ZT2-1.))
         COEF(3,L) = 6.*(D*ALPHA*ONEMZT + C*ZETA)
         COEF(4,L) = 6.*(C - D*ALPHA)/S(I,1)
         COEF(4,L-1) = COEF(4,L) - 6.*D*(1.-ALPHA)/(DEL*ZT2)
         COEF(2,L-1) = COEF(2,L) - DEL*(COEF(3,L) -(DEL/2.)*COEF(4,L-1))
                                        GO TO 68
   62    COEF(2,L) = DIVDIF - S(I,1)*(2.*S(I,4) + S(I+1,4))/6.
         COEF(4,L) = (S(I+1,4)-S(I,4))/S(I,1)
                                        GO TO 68
   63    ONEMZT = GAM*(1. - Z)
         IF (ONEMZT .EQ. 0.)            GO TO 65
         ZETA = 1. - ONEMZT
         ALPHA = ALPH(ZETA)
         C = S(I+1,4)*S(I,6)
         D = S(I,4)/6.
         DEL = ZETA*S(I,1)
         BREAK(L+1) = TAU(I) + DEL
         COEF(2,L) = DIVDIF - S(I,1)*(2.*D + C)
         COEF(4,L) = 6.*(C*ALPHA - D)/S(I,1)
         L = L + 1
         COEF(4,L) = COEF(4,L-1) + 6.*(1.-ALPHA)*C/(S(I,1)*ONEMZT**3)
         COEF(3,L) = COEF(3,L-1) + DEL*COEF(4,L-1)
         COEF(2,L) = COEF(2,L-1)+DEL*(COEF(3,L-1)+(DEL/2.)*COEF(4,L-1))
         COEF(1,L) = COEF(1,L-1)+DEL*(COEF(2,L-1)+(DEL/2.)*(COEF(3,L-1)
     *                  +(DEL/3.)*COEF(4,L-1)))
                                        GO TO 68
   65    COEF(2,L) = DIVDIF
         COEF(3,L) = 0.
         COEF(4,L) = 0.
   68    L = L + 1
   70    BREAK(L) = TAU(I+1)
      L = L - 1
      K = 4
      IFLAG = 1
                                        RETURN
  999 IFLAG = 2
                                        RETURN
      END
