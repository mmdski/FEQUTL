C
C
C
      SUBROUTINE REVERSE_IN_PLACE_D(X, N)

C     Reverse the elements of a vector in place.

      IMPLICIT NONE
      INTEGER N
      REAL*8 X(N)

C     Local
      INTEGER I, IS, IE
      REAL*8 T
C***********************************************************************
      DO 100 I=1,N/2
        IS = I
        IE = N + 1 - I
        T = X(IS)
        X(IS) = X(IE)
        X(IE) = T
100   CONTINUE
      RETURN
      END

C
C
C
      SUBROUTINE INTERSECT(AI, BI, CI, AJ, BJ, CJ,
     O                     STATUS, XS, YS)

C     Find the intersection between two lines in a plane

      IMPLICIT NONE
      CHARACTER STATUS*5

      REAL*8 AI, BI, CI, AJ, BJ, CJ, XS, YS

C     Local

      REAL*8 D, TPA, TPB, DIV, DIST, PDI, PDJ
C***********************************************************************
C     The lines will intersect, be parallel, or be coincident 
C     (a degenerate case of being parallel).

      TPA = AI*BJ
      TPB = AJ*BI
      DIV = MAX(ABS(TPA), ABS(TPB), 1.0d0)
      D = TPA - TPB

      IF(ABS(D)/DIV.LE.1.D-12) THEN
C       Treat the lines as parallel.   Check for near coincidence.
C       Compute the perpendicular distance from the origin to 
C       each line and compare.  If the distances are nearly equal,
C       the lines are coincident. 

        STATUS = 'PAR'
        PDI = ABS(CI)/SQRT(AI**2 + BI**2)
        PDJ = ABS(CJ)/SQRT(AJ**2 + BJ**2) 
        DIST = ABS(PDI - PDJ)

        IF(DIST/PDI.LE.1.E-10) THEN
          STATUS = 'SAME'
        ENDIF
      ELSE
C       Compute the intersection point

        STATUS = 'CROSS'
        XS = (-CI*BJ + CJ*BI)/D
        YS = (-AI*CJ + AJ*CI)/D
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE FIND_INTERSECT(XA, YA, ZA, AA, BA, CA,  NA, 
     I                          A, B, C, XLL, YLL, XUR, YUR, 
     O                          XS, YS, ZS, FLAG)

C     Find the point of intersection between the approach trace polyline
C     and the line given  by (A, B, C)
      IMPLICIT NONE

      INCLUDE 'stdun.cmn'
      INTEGER NA, FLAG
      REAL*8 XA(NA), YA(NA), ZA(NA), AA(NA), BA(NA), CA(NA), 
     A       XLL(NA), YLL(NA), XUR(NA), YUR(NA),
     A        A, B, C, XS, YS, ZS

C     Called program units
      EXTERNAL INTERSECT

C     Local
      INTEGER I

      REAL*8 D, DIST

      CHARACTER STATUS*5
C***********************************************************************
C     Check each line segment until we find an intersection within one.

      DO 100 I=1, NA-1
C        WRITE(STD6,*) ' '
C        WRITE(STD6,*) ' I=',I
C        WRITE(STD6,50) A, B, C
C50    FORMAT(' A=',F15.3,' B=',F15.3,' C=',F15.3)
C        WRITE(STD6,52) AA(I), BA(I), CA(I)
C52    FORMAT(' AA=',F15.3,' BA=',F15.3,' CA=',F15.3)

        CALL  INTERSECT(A, B, C, AA(I), BA(I), CA(I), 
     O                     STATUS, XS, YS)
        IF(STATUS.EQ.'CROSS') THEN
C         Check if it crosses within the line segment. 

C          WRITE(STD6,54) XS, YS
C54    FORMAT(' XS=',F15.2,' YS=',F15.2)

          IF(XS.LE.XUR(I).AND.XS.GE.XLL(I).AND.
     A       YS.LE.YUR(I).AND.YS.GE.YLL(I)) THEN
C           Find the elevation at this point. 
            D = SQRT((XS - XA(I))**2 + (YS - YA(I))**2)
            DIST = SQRT((XA(I+1) - XA(I))**2 + (YA(I+1) - YA(I))**2)
            ZS = ZA(I) + D*(ZA(I+1) - ZA(I))/DIST
            FLAG = 1
            RETURN
          ENDIF            
        ENDIF
100   CONTINUE
      FLAG = 0
      RETURN
          
      END      

       
C
C
C
      SUBROUTINE  MKEMBANK(STDIN, STDOUT, EFLAG)

C     Reads an X, Y, Z trace for the crest of an overflow weir 
C     and an optional trace for the approach surface upstream of the
C     weir and produces the properly formated embankment crest 
C     definition for the EMBANKQ and CULVERT commands.  The weir
C     table references are not supplied and must be added by the
C     user. 

      IMPLICIT NONE
      INTEGER EFLAG, STDIN, STDOUT

      INCLUDE 'arsize.prm'

C     Called program units
      EXTERNAL inline, STRIP_L_BLANKS

C     Local variables.

      INTEGER NC, NA, I, FLAG, REVERSE

      character*8 zone, hgrid, vdatum, unitsys, basis
      REAL WIDTH
      real*8 easting, northing, shift, zmin, xatzmin, yatzmin

C     XC, YC, and ZC give the trace of the crest and XA, YA, and ZA give
C     the trace for the approach.  The X and Y are in the horizontal plane 
C     and Z is in the vertical. 

C     We need to designate 2 different lines:
C     5. The line that is perpendicular to a line at a point. AP, BP, CP
C     6. A generic line.  A, B, C
      REAL*8 XC(PMXOFF), YC(PMXOFF), ZC(PMXOFF),
     A       XA(PMXOFF), YA(PMXOFF), ZA(PMXOFF),
     B       XUR(PMXOFF), YUR(PMXOFF), XLL(PMXOFF), YLL(PMXOFF),
     C       AA(PMXOFF), BA(PMXOFF), CA(PMXOFF),
     D       X, Y, Z, AP, BP, CP, A, B, C, EPS, XS, YS, ZS, ZL, ZR, S,
     E       ZSORG, MINWH, WH

      CHARACTER LINE*80, SURFACE*8, NOTE*80

      DATA EPS/1.D-3/
C     *****************************FORMATS******************************
50    FORMAT(/,' WIDTH=',F10.2)
51    FORMAT(/,' Minimum weir height=',F10.3)
52    FORMAT(/,' SURFACE=',A8)
54    FORMAT(/,' *ERR:772* Line:',A,/,11X,' is unknown in MKEMBANK.')
55    format(/,' SHIFT=',f10.3)
56    FORMAT(/,' *ERR:773* No approach trace intersection at X=',
     A       F12.2,' Y=',F12.2)
58    FORMAT('    OFFSET     CREST     WIDTH  APPROACH SURFACE')
60    FORMAT(F10.1, F10.3, F10.2, F10.2,1X,A8,5X,'''',A)
61    FORMAT(
     A  'App. surface <',F5.2,' below crest. Shifted ',F5.2,' down.')
62    FORMAT(F10.1, F10.3, F10.1,' -99999.0 ' ,1X,A8)
64    FORMAT('  Processing MKEMBANK')
 65   format('ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a        ' BASIS=',a8,/, 'EASTING=',0PF15.3,' NORTHING=',F15.3)

C***********************************************************************

      WRITE(*,64)

      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)


      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      shift = 0.0
      if(line(1:5) == 'SHIFT') then
        read(line,'(6x,f10.0)') shift
      else
        backspace(stdin)
      endif
      
      write(stdout,55) shift

C     Input the constant width for the weir crest.  May add option for 
C     varying width later
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(6X,F10.0)') WIDTH
      WRITE(STDOUT,50) WIDTH

C     Input the minimum weir height.  This will be imposed 
C     with a message as to how much the approach surface was depressed to
C     attain it. 

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      IF(LINE(1:7).EQ.'SURFACE') THEN
C       User has not given a minimum weir height. 
        MINWH = 0.0
        BACKSPACE(STDIN)
      ELSE
        READ(LINE,'(10X,F10.0)') MINWH
        IF(MINWH.LT.0.0) MINWH = 0.0
      ENDIF
      WRITE(STDOUT,51) MINWH


C     Input the constant surface type.  Two types: PAVED, or  GRAVEL
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(8X,A8)') SURFACE
      CALL STRIP_L_BLANKS(
     M                    SURFACE)
      WRITE(STDOUT,52) SURFACE


C     Now get the points on the crest of the weir.  
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
    
      WRITE(STDOUT,'(/,A)') LINE
C     Test for need to reverse
      REVERSE = INDEX(LINE, 'REVERSE')
      
      NC = 0
      NA = 0
      zmin = 1.d50
100   CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        IF(LINE(1:3).NE.'END') THEN
          READ(LINE,*) X, Y, Z
          NC = NC + 1
          XC(NC) = X
          YC(NC) = Y
          ZC(NC) = Z + shift
          if(zc(nc) < zmin) then
            zmin = zc(nc)
            xatzmin = xc(nc)
            yatzmin = yc(nc)
          endif
          GOTO 100
        ENDIF
        easting = xatzmin
        northing = yatzmin
C       Check on reversing the order of the points on the crest trace
        IF(REVERSE.GT.0) THEN
          WRITE(STDOUT,*) ' REVERSING CREST TRACE. NC=', NC
          CALL REVERSE_IN_PLACE_D(XC, NC)
          CALL REVERSE_IN_PLACE_D(YC, NC)
          CALL REVERSE_IN_PLACE_D(ZC, NC)
        ENDIF

C     Get the next line of input.  Check for action to take.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,'(/,A)') LINE
      IF(LINE(1:8).EQ.'APPROACH') THEN
C       Test for need to reverse
        REVERSE = INDEX(LINE, 'REVERSE')
      
C       Read the approach trace.  We assume it extends beyond the end of the 
C       crest trace so that we can find the proper point for computing the
C       approach elevation. 
110     CONTINUE
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
          IF(LINE(1:3).NE.'END') THEN
            READ(LINE,*) X, Y, Z
            NA = NA + 1
            XA(NA) = X
            YA(NA) = Y
            ZA(NA) = Z + shift
            GOTO 110
          ENDIF

C       Check on reversing the order of the points on the approach trace
        IF(REVERSE.GT.0) THEN
          WRITE(STDOUT,*) ' REVERSING APPROACH TRACE. NA=', NA
          CALL REVERSE_IN_PLACE_D(XA, NA)
          CALL REVERSE_IN_PLACE_D(YA, NA)
          CALL REVERSE_IN_PLACE_D(ZA, NA)
        ENDIF
      
C       Read the final END
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        
C       Set the windows for each line segment and find the 
C       coefficients for the line through each line segment.
        DO 120 I=1,NA-1
          XLL(I) = MIN(XA(I), XA(I+1)) - EPS
          YLL(I) = MIN(YA(I), YA(I+1)) - EPS
          XUR(I) = MAX(XA(I), XA(I+1)) + EPS 
          YUR(I) = MAX(YA(I), YA(I+1)) + EPS
          CALL FIND_LINE_COEF(XA(I), YA(I), XA(I+1), YA(I+1),
     A                        A, B , C)
          AA(I) = A
          BA(I) = B
          CA(I) = C
120     CONTINUE

      ELSEIF(LINE(1:3).EQ.'END') THEN
C       No approach trace given.  Approach elevation will be set to 
C       a large negative value to produce a neglibible velocity of
C       approach
      ELSE
        WRITE(STDOUT,54) LINE
        STOP 'Abnormal stop.  Error found.'
      ENDIF

      WRITE(STDOUT,*) ' '
      if(zone /= 'NONE') then
        write(stdout,65) zone, hgrid, vdatum, unitsys, basis,
     a          easting, northing
      endif

      WRITE(STDOUT,58)

C     The basic data is now available.   
      IF(NA.GT.0) THEN
C       The approach trace is given. 

        ZL = -1.D0
        ZR = -1.D0
        ZS = -1.D0
        S = 0.D0

C       At the end points of the crest trace find the elevation 
C       in the approach trace where the perpendicular to the
C       crest trace intersects the  approach trace. 
C       Find the coef for the line through the first two points on 
C       the crest trace.
        CALL FIND_LINE_COEF(XC(1), YC(1), XC(2), YC(2),
     O                 A, B, C)

C       Find the line perpendicular to the line (A, B, C) and passing 
C       through point ( XC(1), YC(1) )
        CALL FIND_PERPENDICULAR( XC(1), YC(1), A, B, C,
     O                          AP, BP, CP)
        CALL FIND_INTERSECT(XA, YA, ZA, AA, BA, CA, NA, AP, BP, CP,
     I                      XLL, YLL, XUR, YUR, 
     O                      XS, YS, ZS, FLAG)
        IF(FLAG.EQ.0) THEN
          WRITE(STDOUT,56) XC(1), YC(1)
          EFLAG = 1
        ELSE
          WH = ZC(1) - ZS
          IF(WH.LT.MINWH) THEN
            WRITE(NOTE,61) MINWH, MINWH - WH
            ZS = ZC(1) - MINWH
          ELSE
            NOTE = ' '
          ENDIF
          WRITE(STDOUT,60) S, ZC(1), WIDTH, ZS, SURFACE, NOTE
C          WRITE(STDOUT,60) S, ZC(1), WIDTH, ZS, SURFACE
        ENDIF
C       Do for the interior points. Find the elevation for the perpendicular
C       for each line segment and take the average.
        DO 200 I=2,NC-1
          S = S + SQRT((XC(I) - XC(I-1))**2 + (YC(I) - YC(I-1))**2)

C         Left-hand line segment
          CALL FIND_LINE_COEF(XC(I-1), YC(I-1), XC(I), YC(I),
     O                   A, B, C)

C         Find the line perpendicular to the line (A, B, C) and passing 
C         through point ( XC(I), YC(I) )
          CALL FIND_PERPENDICULAR( XC(I), YC(I), A, B, C,
     O                            AP, BP, CP)
          CALL FIND_INTERSECT(XA, YA, ZA, AA, BA, CA,  NA, AP, BP, CP,
     I                        XLL, YLL, XUR, YUR, 
     O                        XS, YS, ZL, FLAG)
          IF(FLAG.EQ.0) THEN
            WRITE(STDOUT,56) XC(I), YC(I)
            EFLAG = 1
          ENDIF

C         Right-hand line segment
          CALL FIND_LINE_COEF(XC(I), YC(I), XC(I+1), YC(I+1),
     O                   A, B, C)

C         Find the line perpendicular to the line (A, B, C) and passing 
C         through point ( XC(I), YC(I) )
          CALL FIND_PERPENDICULAR( XC(I), YC(I), A, B, C,
     O                            AP, BP, CP)
          CALL FIND_INTERSECT(XA, YA, ZA, AA, BA, CA,  NA, AP, BP, CP,
     I                        XLL, YLL, XUR, YUR, 
     O                        XS, YS, ZR, FLAG)
          IF(FLAG.EQ.0) THEN
            WRITE(STDOUT,56) XC(I), YC(I)
            EFLAG = 1
          ENDIF
          
C         Output the line
          ZS = 0.5*(ZL + ZR)
          WH = ZC(I) - ZS
          IF(WH.LT.MINWH) THEN
            WRITE(NOTE,61) MINWH, MINWH - WH
            ZS = ZC(I) - MINWH
          ELSE
            NOTE = ' '
          ENDIF
          WRITE(STDOUT,60) S, ZC(I), WIDTH, ZS, SURFACE, NOTE
200     CONTINUE
C       Do the last point
        I = NC
        S = S + SQRT((XC(I) - XC(I-1))**2 + (YC(I) - YC(I-1))**2)
        CALL FIND_LINE_COEF(XC(NC-1), YC(NC-1), XC(NC), YC(NC),
     O                 A, B, C)

C       Find the line perpendicular to the line (A, B, C) and passing 
C       through point ( XC(NC), YC(NC) )
        CALL FIND_PERPENDICULAR( XC(NC), YC(NC), A, B, C,
     O                          AP, BP, CP)
        CALL FIND_INTERSECT(XA, YA, ZA, AA, BA, CA,  NA, AP, BP, CP,
     I                      XLL, YLL, XUR, YUR, 
     O                      XS, YS, ZS, FLAG)
        IF(FLAG.EQ.0) THEN
          WRITE(STDOUT,56) XC(NC), YC(NC)
        ELSE
          WH = ZC(NC) - ZS
          IF(WH.LT.MINWH) THEN
            WRITE(NOTE,61) MINWH, MINWH - WH
            ZS = ZC(NC) - MINWH
          ELSE
            NOTE = ' '
          ENDIF
          WRITE(STDOUT,60) S, ZC(NC), WIDTH, ZS,'END     ' , NOTE
        ENDIF

      ELSE
C      The approach trace is not given. 
        S = 0.D0
        DO 300 I=1,NC-1
          WRITE(STDOUT,62) S, ZC(I), WIDTH, SURFACE
          S = S + SQRT((XC(I) - XC(I+1))**2 + (YC(I) - YC(I+1))**2)
300     CONTINUE         
        WRITE(STDOUT,62) S, ZC(NC), WIDTH, 'END     '
      ENDIF

      RETURN
      END
              



      