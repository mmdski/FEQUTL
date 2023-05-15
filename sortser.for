C     Stuff for sorting and searching of all kinds.

C
C
C
      SUBROUTINE   INSERT
     I                   (VALUE, EPS, MAXN,
     M                    N, VEC,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Insert the real value, VALUE, in the proper location in
C     the ascending-ordered vector of distinct non-negative real
C     values, VEC.  On entry N gives the number of items in the
C     list and on exit gives the new number of items in the list.
C     MAXN gives the maximum extent permitted for VEC.  EPS is the
C     absolute tolerance for VALUE to be taken as matching
C     an item in the list.  If a match occurs VALUE is not
C     added to the list. EFLAG flags that there is no room to
C     add the item.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, MAXN, N
      REAL EPS, VALUE, VEC(MAXN)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     VALUE  - item to insert
C     EPS    - tolerance for a match
C     MAXN   - maximum extent of arrays
C     N      - number of items in list
C     VEC    - vector getting the inserted value
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
C***********************************************************************
      IF(N+1.GT.MAXN) THEN
        EFLAG = 1
      ELSE
C       There is room to add to the list.  Do the limiting items
C       first.
        IF(VALUE.GT.VEC(N) + EPS) THEN
C         Add to the top of the list.
          N = N + 1
          VEC(N) = VALUE
        ELSEIF(VALUE.LT.VEC(1) - EPS) THEN
C         Add to the beginning of the list.
          DO 100 I=N,1,-1
            VEC(I+1) = VEC(I)
 100      CONTINUE
          VEC(1) = VALUE
          N = N + 1
        ELSE
C         Insert at some intermediate point.  Find the interval
C         that contains VALUE.  Check for a match and if not
C         found insert it.
          DO 120 I=2,N
            IF(VALUE.LE.VEC(I)) THEN
C             Found upper end of interval containing VALUE.
              IF(VALUE.GT.VEC(I-1)+EPS.AND.VALUE.LT.VEC(I)-EPS) THEN
C               Insert it.
                DO 110 J=N,I,-1
                  VEC(J+1) = VEC(J)
 110            CONTINUE
                VEC(I) = VALUE
                N = N + 1
                GOTO 130
              ENDIF
            ENDIF
 120      CONTINUE
 130      CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   SORT
     I                 (N,
     M                  X)
 
C     + + + PURPOSE + + +
C     Sort the real array, X(*), of length N into ascending
C     numerical order.
 
      IMPLICIT NONE

C     DUMMY ARGUMENTS
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      REAL X(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N      - Number of items to sort
C     X      - Values to be sorted
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL T
C***********************************************************************
      DO 200 I=2,N
        J = I
        T = X(J)
 100    CONTINUE
          IF(J.GT.1.AND.X(J-1).GT.T) THEN
            X(J) = X(J-1)
            J = J-1
            GOTO 100
          ENDIF
        X(J) = T
 200  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   SORT2
     I                  (N,
     M                   RVAL, RVALA)
 
C     + + + PURPOSE + + +
C     Sort RVAL into ascending order and move other values at
C     same time.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      INTEGER RVALA(N)
      CHARACTER RVAL(N)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N      - Number of items to sort
C     RVAL   - Values to sort
C     RVALA  - Value to move as RVAL is sorted
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, TA
      CHARACTER T*8
C***********************************************************************
C     USE INSERTION SORT FROM PROGRAMMING PEARLS
C     COMM. OF ACM APRIL 1984
 
      DO 200 I=2,N
         J = I
         T = RVAL(J)
         TA = RVALA(J)
 100     CONTINUE
 
            IF(J.GT.1.AND.RVAL(J-1).GT.T) THEN
               RVAL(J) = RVAL(J-1)
               RVALA(J) = RVALA(J-1)
 
               J = J-1
 
               GOTO 100
            ENDIF
         RVAL(J) = T
         RVALA(J) = TA
 
 200  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   SORT2R
     I                   (N,
     M                    RVAL, RVALA)
 
C     + + + PURPOSE + + +
C     Sort RVAL into ascending order and move other values at same time
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      REAL RVAL(N), RVALA(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N      - Number of items to sort
C     RVAL   - Values to sort
C     RVALA  - Value to move as RVAL is sorted
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL T, TA
C***********************************************************************
C     USE INSERTION SORT FROM PROGRAMMING PEARLS
C     COMM. OF ACM APRIL 1984
 
      DO 200 I=2,N
         J = I
         T = RVAL(J)
         TA = RVALA(J)
 100     CONTINUE
 
            IF(J.GT.1.AND.RVAL(J-1).GT.T) THEN
               RVAL(J) = RVAL(J-1)
               RVALA(J) = RVALA(J-1)
 
               J = J-1
 
               GOTO 100
            ENDIF
         RVAL(J) = T
         RVALA(J) = TA
 
 200  CONTINUE
 
      RETURN
      END
