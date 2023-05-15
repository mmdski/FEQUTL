C     Code to help make estimated cross sections below a 
C     reservoir when preconstruction topography is lost


      SUBROUTINE  MAKELAKEXS(STDIN, STDOUT, EFLAG)

C     Uses reservoir capacity table, invert description, and 
C     above-water elevation offset data to estimate below-water
C     cross sections. 

      IMPLICIT NONE
      INTEGER EFLAG, STDIN, STDOUT

      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'

C     Called program units
      EXTERNAL inline, STRIP_L_BLANKS

C     Local variables.

      INTEGER I, J, TABZS, TABSZ, TABLPR, NOFF, NXS, NTAB


      REAL ZATOFF(6), OFFSET(100,6), STATCONV, OFFCONV,
     A     STATION(100), TS(100), ZBOT(100), DF,
     B     YS(100), TP(100), Y(100)

      REAL*8 SUM, SUMT, P

      CHARACTER LINE*80, TABIDZS*16, TABIDSZ*16, TABIDLPR*16,
     A          XSID(100)*2, TABID*16

C     *****************************FORMATS******************************
50    FORMAT(A2,2X,F10.1,6F6.2,F10.0,2F10.1)
52    FORMAT('FEQX',/,'TABID=',A6,' OUT22 NEWBETAM',/,
     A  'STATION=',F10.4,/,'NAVM=   0',/,'NSUB 0.080 0.05 0.08',
     B /, '    OFFSET ELEVATION SUBS')
54    FORMAT(F10.2,F10.3,I5)
C***********************************************************************
C     Get the tabids for the invert description.

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)

      READ(LINE,'(3A16)') TABIDZS, TABIDSZ, TABIDLPR

      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TABIDZS,
     M                             EFLAG,
     O                             TABZS)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TABZS,
     O            EFLAG)

      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TABIDSZ,
     M                             EFLAG,
     O                             TABSZ)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TABSZ,
     O            EFLAG)

      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TABIDLPR,
     M                             EFLAG,
     O                             TABLPR)
      CALL CHKTAB
     I           (2, STDOUT, FTPNT, MFTNUM,
     M            TABLPR,
     O            EFLAG)


C     Hard code the elevations for now- on a tight schedule

      NOFF = 6
      ZATOFF(1) = 800.
      ZATOFF(2) = 600.0
      ZATOFF(3) = 435.
      ZATOFF(4) = 435.
      ZATOFF(5) = 600.
      ZATOFF(6) = 800.0

C     Hardcode station conversion and offset conversion.
      STATCONV = 1.0
      OFFCONV = 2000.0

C     input the station and offset data in advancing station order. 

C     Input the heading line
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      
      I = 0
100   CONTINUE

        CALL inline
     I    (STDIN, STDOUT,
     O     LINE)

        IF(LINE(1:3).NE.'END') THEN
          I = I + 1
          READ(LINE,'(A2,F8.0,6F5.0)') 
     A         XSID(I), STATION(I), (OFFSET(I,J), J=1,6)
          GOTO 100
        ENDIF
      NXS = I
C     Dump values for checking
      P = 1.0D0

      DO 110 I=1,NXS
C       Compute top width at water surface
        TS(I) = OFFCONV*(OFFSET(I,4) - OFFSET(I,3))
C       Compute the elevation of the bottom
        CALL LKTAB
     I    (TABZS, STATION(I), 1,
     O     ZBOT(I), NTAB, DF)
C       Compute the depth.
        YS(I) = ZATOFF(3) - ZBOT(I)

C       Compute the parameter 
        TP(I) = TS(I)/(YS(I)**P)

        WRITE(STDOUT,50) XSID(I), STATION(I), (OFFSET(I,J), J=1,6), 
     A         TS(I), ZBOT(I), YS(I)
110   CONTINUE


C     Compute the surface area. 
      SUM = 0.D0
      DO 120 I=2,NXS
        SUM = SUM + 0.5D0*(TS(I) + TS(I-1))*
     A    (STATION(I) - STATION(I-1))
120   CONTINUE

      WRITE(STDOUT,*) ' Approx surface area at 435=',SUM/43560.D0

      RETURN

C     Compute volume assuming a generalized parabola
      SUM =0.D0
      DO 130 I=2,NXS
        SUM = SUM + 0.5D0*(TS(I)*YS(I) + TS(I-1)*YS(I-1))/(1.D0 + P)*
     A      (STATION(I) - STATION(I-1))
130   CONTINUE

      WRITE(STDOUT,*) 'Approx volume at 435=',SUM/43560.D0

C     Compute volume at 421 feet.
      DO 140 I=1,NXS
        Y(I) = 421.0 - ZBOT(I)
140   CONTINUE

      SUM = 0.D0
      SUMT = 0.D0
      DO 150 I=2,NXS
        SUM = SUM + 0.5D0*(TP(I)*Y(I)**(P+1.D0) +
     A                      TP(I-1)*Y(I-1)**(P+1.D0))/(P +1.D0)*
     A      (STATION(I) - STATION(I-1))
        SUMT = SUMT + 0.5D0*(TP(I)*Y(I)**P +
     A                      TP(I-1)*Y(I-1)**P)*
     A      (STATION(I) - STATION(I-1))
150   CONTINUE

      WRITE(STDOUT,*) 'Approx area at 421=',SUMT/43560.D0
      WRITE(STDOUT,*) 'Approx volume at 421=',SUM/43560.D0
      


C     Compute volume at 384. feet.
      DO 160 I=1,17
        Y(I) = 384.0 - ZBOT(I)
160   CONTINUE

      SUM = 0.D0
      SUMT = 0.D0
      DO 170 I=2,17
        SUM = SUM + 0.5D0*(TP(I)*Y(I)**(P+1.D0) +
     A                      TP(I-1)*Y(I-1)**(P+1.D0))/(P +1.D0)*
     A      (STATION(I) - STATION(I-1))
        SUMT = SUMT + 0.5D0*(TP(I)*Y(I)**P +
     A                      TP(I-1)*Y(I-1)**P)*
     A      (STATION(I) - STATION(I-1))
170   CONTINUE

      WRITE(STDOUT,*) 'Approx area at 384=',SUMT/43560.D0
      WRITE(STDOUT,*) 'Approx volume at 384=',SUM/43560.D0 


C     Output FEQX format using a triangular channel for the 
C     below water portion.

      DO 180 I=1,NXS
        TABID = 'RDXS'//XSID(I)
        WRITE(STDOUT,52) TABID, STATION(I)/3280.84

        WRITE(STDOUT,54) OFFSET(I,1)*609.60, ZATOFF(1)*0.3048, 1
        WRITE(STDOUT,54) OFFSET(I,2)*609.60, ZATOFF(2)*0.3048, 1
        WRITE(STDOUT,54) OFFSET(I,3)*609.60, ZATOFF(3)*0.3048, 2
        WRITE(STDOUT,54) 
     A    609.60*(OFFSET(I,3) + OFFSET(I,4))/2.0, ZBOT(I)*0.3048, 2
        WRITE(STDOUT,54) OFFSET(I,4)*609.60, ZATOFF(4)*0.3048, 3
        WRITE(STDOUT,54) OFFSET(I,5)*609.60, ZATOFF(5)*0.3048, 3
        WRITE(STDOUT,54) OFFSET(I,6)*609.60, ZATOFF(6)*0.3048, -1
        WRITE(STDOUT,*) ' '
        WRITE(STDOUT,*) ' '
180   CONTINUE


      RETURN
    
      END          