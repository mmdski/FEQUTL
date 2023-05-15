C
C
C
      SUBROUTINE  MKWSPRO(STDIN, STDOUT, TABDIR, 
     O                    EFLAG)

C     Make a WSPRO compatiable cross section description from 
C     an FEQX, FEQXLST, or FEQXEXT cross section description.
C     Initial version does not cover all options. 


      IMPLICIT NONE
      INTEGER EFLAG, STDIN, STDOUT
      INTEGER TABDIR(*)

      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'fldway.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'sincom.cmn'


C     Called program units
      EXTERNAL INFEQX, INFQXE

C     Local
      INTEGER I, J, MODE, M3, R, OLDSUB, M7, M
      REAL LEFT, RIGHT, ZMAX, NWSPRO(PMXSUB), SA(PMXSUB)
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8, zone*8, hgrid*8,
     a          vdatum*8, unitsys*8, basis*8


      CHARACTER NAME*7, FORMAT*8
C     **************************FORMATS*********************************
2     FORMAT(5X,A)
4     FORMAT(A)

52    FORMAT(/,' Cross section name=',A)
54    FORMAT(/,' Processing cross section in ',A,' format.')
55    FORMAT(/,' *ERR:751* Format= ',A,' unknown in MKWSPRO.'/,
     A 11X,' Must be: FEQX, FEQXLST, or FEQXEXT')
56    FORMAT('XS   ',A5,F10.2)
57    FORMAT(/,' The cross section in WSPRO format is:')
58    FORMAT(
     A 'GR',8X,F10.2,',',F9.2,5X,F10.2,',',F9.2,5X,F10.2,',',F9.2)
60    FORMAT('N ',8X,7F10.3)
62    FORMAT('SA',8X,7F10.2)
C***********************************************************************
C     Read the name to place with the cross section in the WSPRO format.
      READ(STDIN,2) NAME
      WRITE(STDOUT,52) NAME


C     Get the format of the FEQUTL cross section
      READ(STDIN,4) FORMAT
      WRITE(STDOUT,54) FORMAT

      IF(FORMAT.EQ.'FEQXEXT') THEN
        
        CALL INFQXE
     I             (STDIN, STDOUT,
     M              TABDIR, EFLAG,
     O              TABU, STATU, NPNTU, NSUBU, XU, ZU, SBU, NU, LEFT,
     O              RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX, LSNU, NVARU,
     O              NATYU, YATNU, NNYU, zone, hgrid, vdatum, unitsys,
     O              basis)
 
      ELSEIF(FORMAT.EQ.'FEQX') THEN
        MODE = 1
        CALL INFEQX
     I             (STDIN, STDOUT, MXPNTU, MODE,
     M              TABDIR, EFLAG,
     O              TABU, STATU, NPNTU, NSUBU, NAVMU, XU, ZU, SBU, NU,
     O              LEFT, RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX, 
     O              zone, hgrid, vdatum, unitsys, basis)
      ELSEIF(FORMAT.EQ.'FEQXLST') THEN
        MODE = 2
        CALL INFEQX
     I             (STDIN, STDOUT, MXPNTU, MODE,
     M              TABDIR, EFLAG,
     O              TABU, STATU, NPNTU, NSUBU, NAVMU, XU, ZU, SBU, NU,
     O              LEFT, RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX, zone, 
     O              hgrid, vdatum, unitsys, basis)
      ELSE
        WRITE(STDOUT,55) FORMAT
        STOP 'Abnormal stop.  Error(s) found.'
      ENDIF

C     The FEQXEXT option for varying Manning's n and for using
C     a weighted Manning's n within a subsection are not yet
C     supported. 
      
      WRITE(STDOUT,57)
C     Output the name and the station of the cross section.  
      WRITE(STDOUT,56) NAME, STATU

C     Output the GR cards  with three points per card. 
      M3 = NPNTU/3
C     Get the even multiple of 3 that is less or equal to NPNTU, the number
C     of points on the boundary of the cross section.
      M3 = 3*M3

C     Get the remainder. May be 0, 1, or 2
      R = NPNTU - M3

      DO 100 I=1,M3,3
        WRITE(STDOUT,58) (XU(I+J), ZU(I+J), J=0,2)
100   CONTINUE
      IF(R.GT.0) THEN
        WRITE(STDOUT,58) (XU(M3+J), ZU(M3+J), J=1,R)
      ENDIF

C     Output the values of Manning's n for the subsections.  Find the
C     change in subsection to get the value of n and of the offset at
C     the boundary between adjacent subsections.   The final value of 
C     SA will be the final point on the cross section.  This last
C     point of SA will not be output. 
      OLDSUB = SBU(1)
      M = 0
      DO 200 I=2,NPNTU
        IF(SBU(I).NE.OLDSUB) THEN
C         Save the old value because it applies to the previous 
C         subsection.
          M = M + 1
          NWSPRO(M) = NU(OLDSUB)
          SA(M) = XU(I)
          OLDSUB = SBU(I)
        ENDIF
200   CONTINUE

      M7 = M/7
      M7 = 7*M7
      R = M - M7
      DO 300 I=1,M7,7
        WRITE(STDOUT,60) (NWSPRO(I+J),J=0,6)
300   CONTINUE
      IF(R.GT.0) THEN
        WRITE(STDOUT,60) (NWSPRO(M7+J),J=1,R)
      ENDIF

C     Reduce the count by 1 to clip the last value.
      M = M - 1
      M7 = M/7
      M7 = 7*M7
      R = M - M7
      DO 400 I=1,M7,7
        WRITE(STDOUT,62) (SA(I+J), J=0,6)
400   CONTINUE
      IF(R.GT.0) THEN
        WRITE(STDOUT,62) (SA(M7+J),J=1,R)
      ENDIF
     



      RETURN
      END


