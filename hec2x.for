C     Stuff for processing HEC-2 cross sections


C
C
C
      SUBROUTINE   HECSTB
     I                   (IN, OUT, STDOUT, STADIR, SFAC, STATTB, 
     I                     UFAC)
 
C     + + + PURPOSE + + +
C     Read the HEC-2 input file and establish the sinuosity table
C     for the CHANNEL command.  Convert to station form in the
C     same way as in the HEC2X command.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER IN, OUT, STDOUT
      REAL SFAC, STADIR, STATTB, UFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IN     - Fortran unit number for user input file
C     OUT    - Unit number for the output of the sinuosity table
C     STDOUT - Fortran unit number for user output and messages
C     STADIR - Direction of stationing
C     SFAC   - Scale factor for stations
C     STATTB - Initial station value for all flow lines
C     UFAC   - Factor for conversion of units.
 
C     + + + LOCAL VARIABLES + + +
      REAL AXIS, STCHL, STCHLT, STCHR, STCHRT, STLOB, STROB, TPA, TPL,
     A     TPR, XLCH, XLOBL, XLOBR
      CHARACTER CHR5*5, LINE*80, STAOUT*36
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL MKFMT
 
C     + + + INPUT FORMATS + + +
 2    FORMAT(A80)
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT('OFFS',10X,F10.2,10X,F10.2)
C***********************************************************************
C     CREATE THE SKELETON FORMAT
      STAOUT(1:10) = '(/,''STAT'','
      STAOUT(16:20) = ',10X,'
      STAOUT(26:30) = ',10X,'
      STAOUT(36:36) = ')'
 
C     OUTPUT THE COMMAND AND OTHER HEADER INFORMATION TO THE
C     OUT FILE
      WRITE(OUT,'(A,/,A,/,A)') 'CHANNEL','SINDEF=CUBIC',
     A'HEAD      LOB                AXIS                 ROB'
 
C     INITIALIZE THE STATION VALUES
 
      STLOB = STATTB
      AXIS = STATTB
      STROB = STATTB
 
 
 100  CONTINUE
        READ(IN,2) LINE
 
        IF(LINE(1:2).EQ.'X1') THEN
C         PROCESS THE X1 CARD
 
          READ(LINE,'(16X,5F8.0)',ERR=991) STCHLT, STCHRT,
     A              XLOBL, XLOBR, XLCH
          IF(STADIR.GE.0.0) THEN
            STLOB = STLOB + XLOBL
            AXIS = AXIS + XLCH
            STROB = STROB + XLOBR
          ELSE
            STLOB = STLOB - XLOBL
            AXIS = AXIS - XLCH
            STROB = STROB - XLOBR
          ENDIF
 
          IF(STCHLT.NE.0.0) STCHL = STCHLT
          IF(STCHRT.NE.0.0) STCHR = STCHRT
 
          TPL = STLOB/SFAC
          TPA = AXIS/SFAC
          TPR = STROB/SFAC
          CALL MKFMT
     I              (TPA, 10,
     O               CHR5)
 
 
          STAOUT(11:15) = CHR5
          STAOUT(21:25) = CHR5
          STAOUT(31:35) = CHR5
 
          WRITE(OUT,STAOUT) TPL, TPA, TPR
          WRITE(OUT,52) UFAC*STCHL, UFAC*STCHR
          GOTO 100
 
        ELSEIF(LINE(1:7).EQ.'ENDFILE'.OR.LINE(1:2).EQ.'ER') THEN
C         END OF THE INPUT FILE.  END THE SINUOSITY TABLE
          WRITE(OUT,'(A)')  'END'
 
C         REWIND THE INPUT FILE
          REWIND(IN)
          RETURN
        ENDIF
        GOTO 100
 
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
 
      END
C
C
C
      SUBROUTINE   INHECX
     I                   (STDIN, STDOUT,
     O                    SAVOPT, OUTOPT, BETOPT, MONFLG, EFLAG, MODE,
     O                    BEGTAB, TABINC, BEGSTA, STADIR, SFAC,
     O                    CONFLG, zone, hgrid, vdatum, unitsys, basis,
     O                    easting, northing)
 
C     + + + PURPOSE + + +
C     Input controlling information for HEC2X processing
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER BEGTAB, EFLAG, STDIN, STDOUT, TABINC
      REAL BEGSTA, SFAC, STADIR
      real*8 easting, northing
      CHARACTER BETOPT*8, MODE*8, MONFLG*8, OUTOPT*8, SAVOPT*8,
     A          CONFLG*8, zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     SAVOPT - Function table saving option
C     OUTOPT - Output option for the table file for cross section function
C               tables
C     BETOPT - Option for computing flux coefficients and critical flow
C     MONFLG - Monotonicity flag value
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     MODE   - Mode of processing cross section boundary specification:
C               MODE=1: fixed format.  MODE=2: list format.
C     BEGTAB - beginning cross section table number
C     TABINC - Increment for generating table numbers
C     BEGSTA - Beginning station
C     STADIR - Direction of stationing
C     SFAC   - Scale factor for stations
C     CONFLG - conversion flag for units of measure
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      CHARACTER CHAR4*4, CIN*72, INFILE*64, LINE*80, OUTFIL*64
      LOGICAL THERE
 
C     + + + EXTERNAL FUNCTIONS + + +
      CHARACTER GETTOK*8
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL GETTOK, inline, SETOPT, os_file_style
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(A4,1X,A8,1X,A8)
 2    FORMAT(6X,1X,A64)
 3    FORMAT(7X,1X,A64)
 4    FORMAT(8X,A72)
 5    FORMAT(8X,2I5)
 6    FORMAT(8X,F10.0,F5.0)
 7    FORMAT(A4,1X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' ',A4,'=',A8,1X,A8)
 52   FORMAT(' ','HEC2X INPUT FILE NAME:',A64)
 53   FORMAT(' ','FEQX OUTPUT FILE NAME:',A64)
 54   FORMAT(' ','OPTIONS:',A71)
 55   FORMAT(' BEGINNING TABLE NUMBER=',I5,' TABLE NUMBER INCREMENT=',
     A       I5)
 56   FORMAT(' BEGINNING STATION=',F10.2,' STATION DIRECTION=',F5.0)
 57   FORMAT(' STATIONING DIVISOR TO CONVERT FEET TO DESIRED',
     A        ' UNIT=',F10.2)
 66   FORMAT(' Selection of beta option "NEWBETA" implies checking for',
     A      ' monotonicity.')
 90   FORMAT(' *ERR:646* ',A8,' is invalid mode for HEC2X command.')
 92   FORMAT(' *WRN:561 NO OUTPUT FILE NAME FOR MODE=INDIRECT.',
     A     '  USING NAME:INDIRECT')
 94   FORMAT(' *WRN:515 NO OUTPUT FILE NAME FOR MODE=CHANNEL.',
     A     '  USING NAME:CHANNEL')
C***********************************************************************

c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)
 

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) CHAR4, MODE, CONFLG
      WRITE(STDOUT,51) CHAR4, MODE, CONFLG
 
C     STRIP LEADING BLANKS
      MODE = GETTOK(MODE)
 
C     CHECK FOR VALID MODES
      IF(MODE.NE.'DIRECT'.AND.MODE.NE.'INDIRECT'.AND.MODE.NE.'direct'.
     A      AND.MODE.NE.'indirect'.AND.MODE.NE.'CHANNEL'.AND.
     B      MODE.NE.'channel') THEN
        WRITE(STDOUT,90) MODE
        EFLAG = 1
        RETURN
      ENDIF
 
C     INPUT THE INPUT FILE NAME FOR THE HEC2 DATA
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) INFILE
      CALL MAYBE_ADD_HOME(
     M                    infile)

      call os_file_style(
     m                            infile)
      WRITE(STDOUT,52) INFILE
 
C     CHECK IF THE FILE EXISTS
      INQUIRE(FILE=INFILE, EXIST=THERE)
      IF(THERE) THEN
        OPEN(UNIT=STD48, FILE=INFILE, STATUS='OLD')
      ELSE
        WRITE(STDOUT,*) ' FILE NAMED:',INFILE,' NOT FOUND.'
        WRITE(STDOUT,*) ' CHECK SPELLING OF HEC2X INPUT FILE.'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C     INPUT THE OUTPUT FILE NAME. USED ONLY IF MODE IS INDIRECT BUT INPUT
C     IN ANY CASE.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,3,ERR=991)  OUTFIL
      CALL MAYBE_ADD_HOME(
     M                    outfil)
      call os_file_style(
     m                            outfil)

      WRITE(STDOUT,53) OUTFIL
 
      IF(MODE.EQ.'INDIRECT'.OR.MODE.EQ.'indirect') THEN
        IF(OUTFIL.EQ.' ') THEN
          WRITE(STDOUT,92)
          OUTFIL = 'INDIRECT'
        ENDIF
        OPEN(UNIT=STD49, FILE=OUTFIL, STATUS='UNKNOWN')
      ENDIF
 
      IF(MODE.EQ.'CHANNEL'.OR.MODE.EQ.'channel') THEN
        IF(OUTFIL.EQ.' ') THEN
          WRITE(STDOUT,94)
          OUTFIL = 'CHANNEL'
        ENDIF
        OPEN(UNIT=STD49, FILE=OUTFIL, STATUS='UNKNOWN')
      ENDIF
C     GET THE FEQUTL OPTIONS FOR PROCESSING THE CROSS SECTIONS
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) CIN
      WRITE(STDOUT,54) CIN
 
C     SET THE OPTION FLAGS
 
      CALL SETOPT
     I           (STDOUT, CIN,
     O            SAVOPT, OUTOPT, MONFLG, BETOPT)
      IF(BETOPT(1:7).EQ.'NEWBETA') THEN
        MONFLG = 'MONOTONE'
        WRITE(STDOUT,66)
      ENDIF
 
C     INPUT THE BEGINING TABLE NUMBER AND THE TABLE INCREMENT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,5,ERR=991) BEGTAB, TABINC
      WRITE(STDOUT,55) BEGTAB, TABINC
 
C     INPUT THE BEGINNING STATION AND THE STATIONING DIRECTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) BEGSTA, STADIR
      WRITE(STDOUT,56) BEGSTA, STADIR
 
 
C     INPUT THE STATIONING FACTOR TO CONVERT FEET TO THE DESIRED UNITS
C     BY DIVISION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,7,ERR=991) CHAR4, SFAC
      IF(SFAC.EQ.0.0) SFAC = 1.0
      WRITE(STDOUT,57) SFAC
 
      RETURN
 
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   INSFPX
     I                   (STDIN, STDOUT, SFPNAM, SFPN,
     O                    SAVOPT, OUTOPT, BETOPT, MONFLG, EFLAG, MODE,
     O                    BEGTAB, TABINC, BEGSTA, STADIR, SFAC,
     O                    zone, hgrid, vdatum, unitsys, basis,
     O                    easting, northing)
 
C     + + + PURPOSE + + +
C     Input controlling information for processing cross section
C     descriptions from existing steady-flow profile programs.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER BEGTAB, EFLAG, SFPN, STDIN, STDOUT, TABINC
      REAL BEGSTA, SFAC, STADIR
      real*8 easting, northing
      CHARACTER BETOPT*8, MODE*8, MONFLG*8, OUTOPT*8, SAVOPT*8,
     A          SFPNAM*8, zone*8, hgrid*8, vdatum*8, unitsys*8,
     b          basis*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     SFPNAM - steady-flow program name
C     SFPN   - number of characters in the steady-flow program name
C     SAVOPT - Function table saving option
C     OUTOPT - Output option for the table file for cross section function
C               tables
C     MONFLG - Monotonicity flag value
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     MODE   - mode of process cross sections: DIRECT or INDIRECT
C     BEGTAB - beginning cross section table number
C     TABINC - Increment for generating table numbers
C     BEGSTA - Beginning station
C     STADIR - Direction of stationing
C     SFAC   - Scale factor for stations
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      CHARACTER CHAR4*4, CIN*72, INFILE*64, LINE*80, OUTFIL*64
      LOGICAL THERE
 
C     + + + EXTERNAL FUNCTIONS + + +
      CHARACTER GETTOK*8
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL GETTOK, inline, SETOPT, os_file_style
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(A4,1X,A8)
 2    FORMAT(6X,1X,A64)
 3    FORMAT(7X,1X,A64)
 4    FORMAT(8X,A72)
 5    FORMAT(8X,2I5)
 6    FORMAT(8X,F10.0,F5.0)
 7    FORMAT(A4,1X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' ',A4,'=',A8)
 52   FORMAT(' ',A,' INPUT FILE NAME:',A)
 53   FORMAT(' ','FEQX OUTPUT FILE NAME:',A)
 54   FORMAT(' ','OPTIONS:',A)
 55   FORMAT(' BEGINNING TABLE NUMBER=',I5,' TABLE NUMBER INCREMENT=',
     A       I5)
 56   FORMAT(' BEGINNING STATION=',F10.2,' STATION DIRECTION=',F5.0)
 57   FORMAT(' STATIONING DIVISOR TO CONVERT FEET TO DESIRED',
     A        ' UNIT=',F10.2)
 66   FORMAT(' Selection of beta option "NEWBETA" implies checking for',
     A      ' monotonicity.')
 90   FORMAT(' *ERR:646* ',A8,' is invalid mode for ',A,' command.')
 92   FORMAT(' *WRN:561 No output file name for MODE=INDIRECT.',
     A     '  using name:INDIRECT')
 94   FORMAT(' *WRN:515 No output file name for MODE=CHANNEL.',
     A     '  using name:CHANNEL')
 96   FORMAT(/,' File named:',A,' not found.  check spelling of',
     A       ' file name.')
C***********************************************************************
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) CHAR4, MODE
      WRITE(STDOUT,51) CHAR4, MODE
 
C     STRIP LEADING BLANKS
      MODE = GETTOK(MODE)
 
C     CHECK FOR VALID MODES
      IF(MODE.NE.'DIRECT'.AND.MODE.NE.'INDIRECT'.AND.MODE.NE.'direct'.
     A      AND.MODE.NE.'indirect'.AND.MODE.NE.'CHANNEL'.AND.
     B      MODE.NE.'channel') THEN
        WRITE(STDOUT,90) SFPNAM(1:SFPN), MODE
        EFLAG = 1
        RETURN
      ENDIF
 
C     INPUT THE INPUT FILE NAME FOR THE SFP INPUT DATA STREAM
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991)  INFILE
      CALL MAYBE_ADD_HOME(
     M                     infile)
      call os_file_style(
     m                            infile)
      WRITE(STDOUT,52) SFPNAM(1:SFPN), INFILE
 
C     CHECK IF THE FILE EXISTS
      INQUIRE(FILE=INFILE, EXIST=THERE)
      IF(THERE) THEN
        OPEN(UNIT=STD48, FILE=INFILE, STATUS='OLD')
      ELSE
        WRITE(STDOUT,96) INFILE
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C     INPUT THE OUTPUT FILE NAME. USED ONLY IF MODE IS INDIRECT BUT INPUT
C     IN ANY CASE.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,3,ERR=991)  OUTFIL
      CALL MAYBE_ADD_HOME(
     M                    outfil)
      call os_file_style(
     m                            outfil)
      WRITE(STDOUT,53) OUTFIL
 
      IF(MODE.EQ.'INDIRECT'.OR.MODE.EQ.'indirect') THEN
        IF(OUTFIL.EQ.' ') THEN
          WRITE(STDOUT,92)
          OUTFIL = 'INDIRECT'
        ENDIF
        OPEN(UNIT=STD49, FILE=OUTFIL, STATUS='UNKNOWN')
      ENDIF
 
      IF(MODE.EQ.'CHANNEL'.OR.MODE.EQ.'channel') THEN
        IF(OUTFIL.EQ.' ') THEN
          WRITE(STDOUT,94)
          OUTFIL = 'CHANNEL'
        ENDIF
        IF(SFPNAM(1:SFPN).EQ.'WSPRO') THEN
C         Open a scratch file attached to STD50
          OPEN(UNIT=STD50, STATUS='SCRATCH')
        ENDIF
        OPEN(UNIT=STD49, FILE=OUTFIL, STATUS='UNKNOWN')
      ENDIF
C     GET THE FEQUTL OPTIONS FOR PROCESSING THE CROSS SECTIONS
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) CIN
      WRITE(STDOUT,54) CIN
 
C     SET THE OPTION FLAGS
 
      CALL SETOPT
     I           (STDOUT, CIN,
     O            SAVOPT, OUTOPT, MONFLG, BETOPT)
      IF(BETOPT(1:7).EQ.'NEWBETA') THEN
        MONFLG = 'MONOTONE'
        WRITE(STDOUT,66)
      ENDIF
 
C     INPUT THE BEGINING TABLE NUMBER AND THE TABLE INCREMENT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,5,ERR=991)  BEGTAB, TABINC
      WRITE(STDOUT,55) BEGTAB, TABINC
 
C     INPUT THE BEGINNING STATION AND THE STATIONING DIRECTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991)  BEGSTA, STADIR
      WRITE(STDOUT,56) BEGSTA, STADIR
 
 
C     INPUT THE STATIONING FACTOR TO CONVERT FEET TO THE DESIRED UNITS
C     BY DIVISION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,7,ERR=991) CHAR4, SFAC
      IF(SFAC.EQ.0.0) SFAC = 1.0
      WRITE(STDOUT,57) SFAC
 
      RETURN
 
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   FIXNH
     I                  (STDOUT, NUMNH, STN, SIZEXZ,
     M                   NPNT, X, Z)
 
C     + + + PURPOSE + + +
C     Add coordinate points to the cross section boundary to match
C     offsets at which the boundaries for horizontal variation of
C     n are given.  Needed for handling HEC-2 input resulting from
C     translation from the Kansas City District steady water surface
C     profile program.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
C      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNT, NUMNH, SIZEXZ, STDOUT
      REAL STN(NUMNH), X(SIZEXZ), Z(SIZEXZ)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NUMNH  - number of n values given
C     STN    - offsets for n values
C     NPNT   - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, JADD, JS
      REAL DIFF, XNH, ZNH
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL SORT2R
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' *ERR:719* Cross section offset non-increasing at',
     A       ' offset=',F10.2)
 52   FORMAT(/,' *ERR:720* No space left in cross section boundary',
     A   'when adding point for NH card.')
C***********************************************************************
      JADD = NPNT
      JS = 2
      DO 200 I=1,NUMNH
        XNH = STN(I)
        IF(XNH.GE.X(NPNT)) GOTO 110
        DO 100 J=JS,NPNT
          IF(XNH.GE.X(J-1).AND.XNH.LE.X(J)) THEN
C           In the interval.  Is it at one of the end points?
            IF(XNH.NE.X(J-1).AND.XNH.NE.X(J)) THEN
C             Not at the end points.  Add a point to the boundary.
              DIFF = X(J) - X(J-1)
              IF(DIFF.LE.0.0) THEN
                WRITE(STDOUT,50) X(J)
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              ZNH = Z(J-1) + (XNH - X(J-1))*(Z(J) - Z(J-1))/DIFF
              JADD = JADD + 1
              IF(JADD.GT.SIZEXZ) THEN
                WRITE(STDOUT,52) JADD
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              X(JADD) = XNH
              Z(JADD) = ZNH
              JS = J
              GOTO 110
            ENDIF
          ENDIF
 100    CONTINUE
 110    CONTINUE
 200  CONTINUE
      IF(JADD.GT.NPNT) THEN
C       Points were added.
        NPNT = JADD
        CALL SORT2R
     I             (NPNT,
     M              X, Z)
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   SCNHEC
     I                   (IN, STDOUT, MXPNTU, STADIR,
     M                    STATTB, EFLAG, NCFLAG, GRFLAG, LNFLAG,
     O                    NPNTU, NSUBU, NAVMU, SCALE, SHIFT, XU, ZU,
     O                    SBU, NU, LEFT, RIGHT, SECID, NHFLAG, X4FLAG)
 
C     + + + PURPOSE + + +
C     Scan HEC2 input in the file, IN, and get the next
C     cross section and return the values needed for computing
C     a cross section table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, GRFLAG, IN, LNFLAG, MXPNTU, NAVMU, NCFLAG, NHFLAG,
     A        NPNTU, NSUBU, STDOUT, X4FLAG
      INTEGER SBU(MXPNTU)
      REAL LEFT, NU(20), RIGHT, SCALE, SHIFT, STADIR, STATTB,
     A     XU(MXPNTU), ZU(MXPNTU)
      CHARACTER SECID*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IN     - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     MXPNTU - Maximum number of points on a cross section boundary
C     STADIR - Direction of stationing
C     STATTB - Initial station value for all flow lines
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NCFLAG - Flag for the NC card
C     GRFLAG - Flag set when GR cards are encountered
C     LNFLAG - Flag to signal the state of the current line buffer
C     NPNTU  - Number of points on boundary of a cross section
C     NSUBU  - Number of subsections
C     NAVMU  - Flag for averaging roughness
C     SCALE  - Scale factor to apply to the offsets
C     SHIFT  - Value of vertical shift to apply to the points on the
C               cross section boundary
C     XU     - Offsets for points on boundary of cross section
C     ZU     - Elevation of points on boundary of cross section
C     SBU    - Subsection numbers for the line segments-upstream
C               location
C     NU     - Vector for Manning's n values
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
C     SECID  - Section identification value
C     NHFLAG - Flag for the NH card
C     X4FLAG - Flag for X4 card
 
C     + + + SAVED VALUES + + +
      INTEGER NELT, NUMNH, NUMST
      REAL EL(254), ELT(20), STA(254), STAT(20), STCHL, STCHR, STN(25),
     A     VALN(25), XNCH, XNL, XNR
      CHARACTER LINE*80
      SAVE EL, ELT, LINE, NELT, NUMNH, NUMST, STA, STAT, STCHL, STCHR,
     A     STN, VALN, XNCH, XNL, XNR
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, JE, JS, K
      REAL DUMNUM, EL2(254), NUMSTT, PXSECE, PXSECR, STA2(254), STCHLT,
     A     STCHRT, XLCH, XNCHT, XNLT, XNRT
 
C     + + + INTRINSICS + + +
      INTRINSIC MIN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL FIXNH, inline, SORT2R
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *WRN:586* Left overbank n$=$0.0 on first NC card.',
     A      '  Setting n to 1.0.')
 51   FORMAT(' *WRN:587* Channel n$=$0.0 on first NC card.',
     A      '  Setting n to 1.0.')
 52   FORMAT(' *WRN:588* Right overbank n$=$0.0 on first NC card.',
     A      '  Setting n to 1.0.')
 53   FORMAT(' *ERR:637* CARDS OUT OF ORDER. LAST CARD READ:')
 54   FORMAT(' *ERR:638* GR CARDS FOUND BUT NC OR NH CARD MISSING.')
 55   FORMAT(' *ERR:639* X1 REFERS TO PREVIOUS GR DATA BUT NO GR DATA',
     A       ' IS IN HAND.')
 56   FORMAT(' *ERR:640* NUMNH=',I3,' > 20 ON NH CARD.')
 58   FORMAT(' *ERR:641* NELT=',I3,' > 20 ON X4 CARD.')
 60   FORMAT(' *ERR:642* NH CARD DOES NOT COVER CROSS SECTION.')
 61   FORMAT('  Last offset on NH card=',F10.1,' Last offset on GR',
     A       ' card=',F10.1)
 62   FORMAT(' *ERR:643* MANNING N=',F8.3,' <= 0 IN SUBSECTION',I4)
C***********************************************************************
C     GET THE NEXT LINE FROM THE HEC2 INPUT FILE
 
 90   CONTINUE
        IF(LNFLAG.EQ.0) THEN
          CALL inline
     I              (IN, STDOUT,
     O               LINE)
        ELSE
C         CLEAR THE LINE BUFFER FLAG TO SIGNAL THAT THE BUFFER HAS
C         BEEN PROCESSED.
          LNFLAG = 0
        ENDIF
 
C       USE AN IF-THEN-ELSE SEQUENCE TO DECIDE WHAT TO DO.  NOTE THAT
C       NCFLAG AND GRFLAG MUST BE CLEARED BY THE CALLING PROGRAM.
 
        IF(LINE(1:2).EQ.'NC') THEN
C         WE HAVE A CARD WHICH MIGHT CHANGE VALUES OF N
          READ(LINE,'(2X,F6.0,9F8.0)',ERR=991) XNLT, XNRT, XNCHT
          IF(NCFLAG.EQ.0) THEN
            IF(XNLT.EQ.0.0) THEN
              WRITE(STDOUT,50)
C              EFLAG = 1
              XNLT = 1.0
            ENDIF
            IF(XNCHT.EQ.0.0) THEN
              WRITE(STDOUT,51)
C              EFLAG = 1
              XNCHT = 1.0
            ENDIF
            IF(XNRT.EQ.0.0) THEN
C              EFLAG = 1
              XNRT = 1.0
              WRITE(STDOUT,52)
            ENDIF
            IF(EFLAG.EQ.0) THEN
              NCFLAG = 1
              NHFLAG = 0
            ENDIF
          ELSE
            NHFLAG = 0
C           Take value of zero to mean that n is unchanged.
            IF(XNLT.EQ.0.0) THEN
              XNLT = -1.0
            ENDIF
            IF(XNCHT.EQ.0.0) THEN
              XNCHT = -1.0
            ENDIF
            IF(XNRT.EQ.0.0) THEN
             XNRT = -1.0
            ENDIF
            
          ENDIF
 
          IF(XNLT.GT.0.0.AND.XNLT.LE.9.0) THEN
            XNL = XNLT
          ELSE
C            XNL = 9.0
          ENDIF
          IF(XNRT.GT.0.0.AND.XNRT.LE.9.0) THEN
            XNR = XNRT
          ELSE
C            XNR = 9.0
          ENDIF
          IF(XNCHT.GT.0.0.AND.XNCHT.LE.9.0) THEN
            XNCH = XNCHT
          ELSE
C            XNCH = 9.0
          ENDIF
 
        ELSEIF(LINE(1:2).EQ.'X1') THEN
C         PROCESS THE X1 CARD
 
          READ(LINE,'(A8,3F8.0,16X,3F8.0)',ERR=991) SECID, NUMSTT,
     A              STCHLT, STCHRT, XLCH, PXSECR, PXSECE
          SCALE = PXSECR
          SHIFT = PXSECE
          IF(STADIR.GE.0.0) THEN
            STATTB = STATTB + XLCH
          ELSE
            STATTB = STATTB - XLCH
          ENDIF
 
          IF(NUMSTT.GT.0.0) NUMST = NUMSTT + .1
          IF(STCHLT.NE.0.0) STCHL = STCHLT
          IF(STCHRT.NE.0.0) STCHR = STCHRT
 
C         LOOK AHEAD FOR AN X4 CARD.  MUST FOLLOW X1 CARD IMMEDIATELY.
 
          CALL inline
     I              (IN, STDOUT,
     O               LINE)
          IF(LINE(1:2).EQ.'X4') THEN
            X4FLAG = 1
 
C           INPUT THE X4 CARD
 
            READ(LINE,'(2X,F6.0,9F8.0)',ERR=991) DUMNUM, ELT(1),
     A       STAT(1), ELT(2), STAT(2), ELT(3), STAT(3), ELT(4),
     B       STAT(4), ELT(5)
            NELT = DUMNUM + .1
            IF(NELT.GT.20) THEN
              WRITE(STDOUT,58) NELT
              EFLAG = 1
              RETURN
            ENDIF
 
C           SET THE OFFSET
            J = 0
 
 230        CONTINUE
              IF(NELT.GT.4+J) THEN
                CALL inline
     I                    (IN, STDOUT,
     O                     LINE)
                IF(LINE(1:2).NE.'X4') THEN
                  WRITE(STDOUT,53)
                  WRITE(STDOUT,'(1H ,A80)') LINE
                  EFLAG = 1
                  RETURN
                ELSE
                  READ(LINE,'(2X,F6.0,9F8.0)', ERR=991) STAT(J+5),
     A             ELT(J+6), STAT(J+6), ELT(J+7), STAT(J+7), ELT(J+8),
     B             STAT(J+8), ELT(J+9), STAT(J+9), ELT(J+10)
                ENDIF
 
                J = J + 5
                GOTO 230
              ENDIF
          ELSE
            LNFLAG = 1
            X4FLAG = 0
            NELT = 0
          ENDIF
 
          IF(NUMSTT.EQ.0) THEN
C           IS GR DATA IN HAND ALREADY?
            IF(GRFLAG.EQ.0) THEN
              WRITE(STDOUT,55)
              EFLAG = 1
              RETURN
            ENDIF
            IF(NCFLAG.EQ.0.AND.NHFLAG.EQ.0) THEN
C             NO N VALUES IN HAND YET. CARDS OUT OF ORDER
              WRITE(STDOUT,54)
              EFLAG = 1
              RETURN
            ENDIF
 
C           CROSS SECTION IN HAND. PROCESS IT INTO FEQX FORM.
 
            GOTO 500
 
          ENDIF
 
        ELSEIF(LINE(1:2).EQ.'GR') THEN
C         PROCESS GR CARDS.  NUMST > 0 HERE.
 
          GRFLAG = 1
 
          DO 100 I=1,NUMST,5
            JS =  I
            JE = MIN(JS+4, NUMST)
            READ(LINE,'(2X,F6.0,9F8.0)',ERR=991) (EL(J),STA(J),J=JS,JE)
            IF(JE.LT.NUMST) THEN
C             GET THE NEXT CARD
              CALL inline
     I                  (IN, STDOUT,
     O                   LINE)
              IF(LINE(1:2).NE.'GR') THEN
                WRITE(STDOUT,53)
                WRITE(STDOUT,'(1H ,A80)') LINE
                EFLAG = 1
                RETURN
              ENDIF
            ENDIF
 100      CONTINUE
 
          IF(NCFLAG.EQ.0.AND.NHFLAG.EQ.0) THEN
C           NO N VALUES IN HAND YET. CARDS OUT OF ORDER
            WRITE(STDOUT,54)
            EFLAG = 1
            RETURN
          ENDIF
 
C         CROSS SECTION DEFINED.  CONVERT IT TO THE FEQX DESCRIPTION
 
          GOTO 500
        ELSEIF(LINE(1:2).EQ.'ER') THEN
C         END OF HEC2 INPUT
          SECID = 'END'
          RETURN
        ELSEIF(LINE(1:7).EQ.'ENDFILE') THEN
C         END OF FILE IN HEC2 INPUT
          SECID = 'END'
          RETURN
        ELSEIF(LINE(1:2).EQ.'NH') THEN
C         INPUT THE NH CARD AND SET FLAGS TO REFLECT VARIABLE
C         HORIZONTAL N VALUES.
 
          READ(LINE,'(2X,F6.0,9F8.0)',ERR=991) DUMNUM, VALN(1),
     A     STN(1), VALN(2), STN(2), VALN(3), STN(3), VALN(4), STN(4),
     B     VALN(5)
          NUMNH = DUMNUM + .1
          IF(NUMNH.GT.20) THEN
            WRITE(STDOUT,56) NUMNH
            EFLAG = 1
            RETURN
          ENDIF
 
C         SET THE OFFSET
          J = 0
 
 200      CONTINUE
            IF(NUMNH.GT.4+J) THEN
              CALL inline
     I                  (IN, STDOUT,
     O                   LINE)
              IF(LINE(1:2).NE.'NH') THEN
                WRITE(STDOUT,53)
                WRITE(STDOUT,'(1H ,A80)') LINE
                EFLAG = 1
                RETURN
              ELSE
                READ(LINE,'(2X,F6.0,9F8.0)', ERR=991) STN(J+5),
     A            VALN(J+6), STN(J+6), VALN(J+7), STN(J+7), VALN(J+8),
     B            STN(J+8), VALN(J+9), STN(J+9), VALN(J+10)
              ENDIF
 
              J = J + 5
              GOTO 200
            ENDIF
 
          NHFLAG = 1
C          NCFLAG = 0
 
        ENDIF
 
C       READ THE NEXT CARD
        GOTO 90
 
 500  CONTINUE
 
 
C     PROCESS THE EFFECT OF THE X4 CARD IF ANY WAS FOUND.  ASSUME THAT
C     THE X4 CARD MODIFIES BUT DOES NOT CHANGE THE RETAINED GR CARD
C     DATA.  TRANSFER THE RETAINED GR DATA TO THE WORK SPACE IN
C     ANY CASE.
 
      NPNTU = NUMST
      DO 300 J=1,NUMST
        STA2(J) = STA(J)
        EL2(J) = EL(J)
 300  CONTINUE
 
 
      IF(X4FLAG.EQ.1) THEN
C       X4 IS ACTIVE.  APPEND TO THE LIST OF THE CURRENT GR CARD
C       PATTERN AND SORT TO PLACE THE STATIONS IN ASCENDING ORDER
        NPNTU = NUMST + NELT
        DO 310 J=NUMST+1,NPNTU
          STA2(J) = STAT(J-NUMST)
          EL2(J) = ELT(J-NUMST)
 310    CONTINUE
 
C       SORT THE VALUES INTO INCREASING ORDER OF STATION
 
        CALL SORT2R
     I             (NPNTU,
     M              STA2, EL2)
        X4FLAG = 0
      ENDIF
 
C     CONVERT THE CURRENT HEC2 CROSS SECTION TO FEQX FORM.  THEN USE
C     THE USUAL METHOD FOR COMPUTING CROSS SECTION TABLES.
 
      NAVMU = 0
C      WRITE(STDOUT,*) ' NUMST=',NUMST,' SCALE=',SCALE,' SHIFT=',SHIFT
      IF(SCALE.EQ.0.0) SCALE = 1.0
C     ASSIGN SUBSECTIONS
      IF(NCFLAG.EQ.1.AND.NHFLAG.EQ.0) THEN
        NSUBU = 3
        NU(1) = XNL
        NU(2) = XNCH
        NU(3) = XNR
 
       DO 510 J=1,NPNTU
C         WRITE(STDOUT,*) ' J=',J,' STA2(J)=',STA2(J),' EL2(J)=',EL2(J)
          XU(J) = STA2(J)*SCALE
          ZU(J) = EL2(J) + SHIFT
          IF(STA2(J).LT.STCHL) THEN
            SBU(J) = 1
          ELSEIF(STA2(J).GE.STCHL.AND.STA2(J).LT.STCHR) THEN
            SBU(J) = 2
          ELSE
            SBU(J) = 3
          ENDIF
 510    CONTINUE
        SBU(NPNTU) = -1
      ELSE
        IF(NHFLAG.NE.1) THEN
          WRITE(STDOUT,*) ' BUG IN FEQX2. NHFLAG INVALID.'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
C       Check for matching of boundary points.  If match not found
C       add a point to the boundary.
        CALL FIXNH
     I            (STDOUT, NUMNH, STN, 254,
     M             NPNTU, STA2, EL2)
 
        NSUBU = NUMNH
        DO 210 J=1,NUMNH
          IF(VALN(J).LE.0.0) THEN
            WRITE(STDOUT,62) VALN(J), J
            EFLAG = 1
            VALN(J) = 1.0
          ENDIF
          IF(VALN(J).GT.9.0) THEN
            VALN(J) = 9.0
          ENDIF
          NU(J) = VALN(J)
 210    CONTINUE
 
        IF(STN(NUMNH).LT.STA2(NPNTU)) THEN
C         SPECIFICATION OF N DOES NOT COVER THE CROSS SECTION
          WRITE(STDOUT,60)
          WRITE(STDOUT,61) STN(NUMNH), STA2(NPNTU)
          EFLAG = 1
          STN(NUMNH) = STA2(NPNTU)
        ENDIF
C       ASSIGN VALUES TO THE LINE SEGMENTS ON THE CROSS SECTION
C       BOUNDARY
        K = 1
        DO 220 J=1,NPNTU
          XU(J) = STA2(J)*SCALE
          ZU(J) = EL2(J) + SHIFT
 
          SBU(J) = -1
          IF(STA2(J).LT.STN(K)) THEN
            SBU(J) = K
          ELSE
            K = K + 1
            SBU(J) = K
          ENDIF
 220    CONTINUE
        SBU(NPNTU) = -1
      ENDIF
 
 
C     SET SCALE AND SHIFT TO DEFAULT VALUES
 
      SCALE = 1.0
      SHIFT = 0.0
      LEFT = 1.E20
      RIGHT = -1.E20
 
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   HEC2X
     I                  (STDIN, STDOUT, STDTAB, NFAC,
     M                   TABDIR, FTP,
     O                   EFLAG)
 
C     + + + PURPOSE + + +
C     Abstract HEC2 cross sections from a HEC2 input deck.
C     Two options:  Compute the cross section tables directly or
C     compute the FEQX input form for later editing before the
C     cross section tables are computed.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL NFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     TABDIR - Table directory to remember table numbers
C     FTP    - next open location in the function table storage
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'stdun.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xtadd.cmn'
 
C     + + + SAVED VALUES + + +
      CHARACTER STAOUT*18
      SAVE STAOUT
 
C     + + + LOCAL VARIABLES + + +
      INTEGER BEGTAB, GRFLAG, I, IN, J, JJ, LIM, LNFLAG, MESG, NCFLAG, 
     A        NHFLAG, TABINC, X4FLAG, ITEMP, locflag
      REAL BEGSTA, LEFT, RIGHT, SCALE, SFAC, SHIFT, STADIR, STATT, TP,
     A     XOLD, ZMAX, UFAC
      CHARACTER BETOPT*8, F10X*5, MODE*8, MONFLG*8, OUTOPT*8, SAVOPT*8,
     A          SECID*8, CONFLG*8, COFF*10, CELEV*10,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8

 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSTAB, HECSTB, INHECX, MKFMT, SCNHEC, TABCHK, TABOUT,
     A         VAR_DECIMAL, GET_INTERNAL_TAB_NUMBER, STRIP_L_BLANKS
 
C     + + + DATA INITIALIZATIONS + + +
      DATA STAOUT/'(''STATION='','/
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(//)
 51   FORMAT(/,' TABID=',I8)
 52   FORMAT(' NAVM=',I5,'  SCALE=',F10.3,' SHIFT=',F10.3)
 53   FORMAT(' STATION=',F10.3,' LEFT=',F10.1,' RIGHT=',F10.1)
 54   FORMAT(' STATION=',F10.3)
 55   FORMAT(1X,'NSUB',I5,6F6.3)
 56   FORMAT(10X,6F6.3)
 58   FORMAT(' *ERR:504* NUMBER OF SUBSECTIONS=',I5,' > ',I5)
 61   FORMAT(' *ERR:506* ONLY ONE POINT GIVEN ON BOUNDARY OF THE',
     1  ' CROSS SECTION.')
 67   FORMAT('0*WRN:554* Extending left end of cross section',
     A   ' by ',F8.3)
 68   FORMAT(/,' *WRN:555* Extending right end of cross section',
     A   ' by ',F8.3)
 69   FORMAT(/,' *WRN:556* Some point in cross section higher than',
     A  ' either end.',/,10X,'All area above minimum end elevation',
     b  ' is ignored.')
 70   FORMAT('     OFFSET ELEVATION SUBS',3X,A8)
 72   FORMAT(' ',F10.2,F10.2, I5)
 75   FORMAT('TABID=',I8,2X,A8,1X,A8,1X,A8,1X,A8)
 76   FORMAT('TABID=',I8,'  EXTEND',2X,A8,1X,A8,1X,A8,1X,A8)
 77   FORMAT('FEQX')
 79   FORMAT('NAVM=0')
 80   FORMAT('NSUB',I5,6F10.3)
 81   FORMAT('    OFFSET ELEVATION SUBS',2X,A8)
 82   FORMAT(2A10,I5)
 83   FORMAT(9X,6F10.3)
 84   FORMAT(I10)
 85   format(' ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a' BASIS=',a8,/,' EASTING=',f15.3,' NORTHING=',f15.3)
 86   format('ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a' BASIS=',a8,/,'EASTING=',f15.3,' NORTHING=',f15.3)
C***********************************************************************
C     Clear the flags for values not set
      GISID = ' '
      NORTHING = 0.D0
      EASTING = 0.D0
      zone = 'NONE'
      hgrid = 'NONE'
      vdatum = 'NONE'
      unitsys = 'NONE'
      basis = 'NONE'

C     CLEAR THE ERROR FLAG.  USED TO DETECT PROBLEMS THAT REQUIRE
C     EARLY EXITS FROM PROCESSING
 
      EFLAG = 0
 
C     CLEAR THE VALUES NOT NEEDED TO REPRESENT THE CROSS SECTIONS BUT
C     STILL USED IN THE COMPUTATION OF THE TABLE.
 
      DO 95 I=1,PMXSUB
        NNYU(I) = 0
        NVARU(I) = 0
 95   CONTINUE
 
      DO 96 I=1,PMXPNT
        SNU(I) = 1.0
        LSNU(I) = 0.0
 96   CONTINUE
 
      SNFLGU = 0
 
C     INPUT THE CONTROLLING DATA
      CALL INHECX
     I           (STDIN, STDOUT,
     O            SAVOPT, OUTOPT, BETOPT, MONFLG, EFLAG, MODE, BEGTAB,
     O            TABINC, BEGSTA, STADIR, SFAC, CONFLG,
     O            zone, hgrid, vdatum, unitsys, basis,
     O            easting, northing)
c     Location information only supplied if vdatum or unit system is
c     active!
      locflag = 0
      if(vdatum /= 'NONE' .or. vdatum /= 'NA') then
        locflag = 1
      endif
      if(unitsys /= 'NONE' .or. unitsys /= 'NA') then
        locflag =1 
      endif
 
      IF(EFLAG.EQ.1) RETURN

C     Set conversion options
      IF(ABS(NFAC - 1.0).LE.0.0001) THEN
C       Metric
        IF(CONFLG.NE.'        ') THEN
C         Assume units are to be converted from english to metric
          UFAC = 0.3048
        ELSE
          UFAC = 1.0
        ENDIF
      ELSE
        IF(CONFLG.NE.'        ') THEN
C         Assume units are to be converted from metric to english
          UFAC = 3.280840
        ELSE
          UFAC = 1.0
        ENDIF
      ENDIF
 
C     SET THE INPUT UNIT NUMBER
 
      IN = STD48
 
      IF(MODE.EQ.'CHANNEL'.OR.MODE.EQ.'channel') THEN
C       DO A SCAN OF THE INPUT FILE TO ESTABLISH THE
C       SINUOSITY TABLE.
 
        CALL HECSTB
     I             (IN, STD49, STDOUT, STADIR, SFAC, BEGSTA, UFAC)
 
      ENDIF
 
C     GET THE NEXT CROSS SECTION FROM THE HEC2 INPUT.  WHEN A
C     COMPLETE CROSS SECTION HAS BEEN DEFINED, RETURN WITH
C     THE FEQX INPUT AND PROCESS.  RETURN A FLAG WHEN NO
C     CROSS SECTION HAS BEEN FOUND.  THIS TERMINATES THE PROCESSING
C     OF THE FILE.
 
C     START THE TABLE NUMBERS
 
      TABU = BEGTAB
 
C     ENABLE CONVEYANCE WARNING MESSAGES
 
      NOCM = 0
      SLOT = 0.0
 
 
C     START THE STATION AT BEGSTA
 
      STATU = BEGSTA
 
C     CLEAR THE FLAGS FOR CARD SEQUENCE CHECKING
      NCFLAG = 0
      GRFLAG = 0
      NHFLAG = 0
      LNFLAG = 0
      X4FLAG = 0
 
 100  CONTINUE
 
        IF(TABU.GE.0) CALL TABCHK
     I                           (STDOUT, PMXTAB,
     M                            TABU, TABDIR, EFLAG)
 
        CALL SCNHEC
     I             (IN, STDOUT, MXPNTU, STADIR,
     M              STATU, EFLAG, NCFLAG, GRFLAG, LNFLAG,
     O              NPNTU, NSUBU, NAVMU, SCALE, SHIFT, XU, ZU, SBU, NU,
     O              LEFT, RIGHT, SECID, NHFLAG, X4FLAG)
 
        IF(SECID(1:3).EQ.'END') THEN
          CLOSE(STD48)
          IF(MODE.EQ.'indirect'.OR.MODE.EQ.'INDIRECT')THEN
            CLOSE(STD49)
          ELSEIF(MODE.EQ.'CHANNEL'.OR.MODE.EQ.'channel') THEN
            WRITE(STD49,'(A)') 'ENDCHAN'
            CLOSE(STD49)
          ENDIF
          RETURN
        ENDIF
        IF(EFLAG.GT.0) RETURN
 
        WRITE(STDOUT,50)
        WRITE(STDOUT,51) TABU
        LEFT = SCALE*LEFT
        RIGHT = SCALE*RIGHT
        IF(LEFT.GE.RIGHT) THEN
          WRITE(STDOUT,54) STATU/SFAC
        ELSE
          WRITE(STDOUT,53) STATU/SFAC, LEFT, RIGHT
        ENDIF
        write(stdout,85)  zone, hgrid, vdatum, unitsys, basis,
     a                    easting, northing
        WRITE(STDOUT,52) NAVMU, SCALE, SHIFT
 
        IF(NSUBU.GT.PMXSUB) THEN
          WRITE(STDOUT,58) NSUBU, PMXSUB
          EFLAG = 1
          NSUBU=PMXSUB
        ENDIF
 
C        WRITE(STDOUT,55)  NSUBU, (NU(J),J=1,NSUBU)
        WRITE(STDOUT,55) NSUBU, (NU(J),J=1,MIN(6,NSUBU))
        DO 280 JJ=7,NSUBU,6
          LIM = NSUBU - JJ
          IF(LIM.GT.5) LIM = 5
          WRITE(STDOUT,56) (NU(JJ+J), J=0,LIM)
280     CONTINUE
 
C       SHIFT and SCALE applied in SCNHEC already.  Check
C       for monotonicity. 
        ZMAX = -9999999.
        XOLD = -1.E20
 
        WRITE(STDOUT,70) SECID
 
        DO 110 J=1,NPNTU
          WRITE(STDOUT,72) XU(J), ZU(J), SBU(J)
 
 
C         FIND MAXIMUM ELEVATION IN CROSS SECTION FOR LATER CHECKING
          ZMAX = MAX(ZU(J), ZMAX)
 
C         CHECK FOR MONOTONICITY OF TOP WIDTH.  THIS REQUIRES THAT
C         THE OFFSET NEVER DECREASE.
 
          IF(MONFLG.EQ.'MONOTONE') THEN
            IF(XU(J).LT.XOLD) THEN
              WRITE(STDOUT,*)
     A  ' *ERR:508* SECTION VIOLATES MONOTONICITY AT OFFSET=', XU(J)
              EFLAG = EFLAG + 1
            ENDIF
          ENDIF
 
          XOLD = XU(J)
 
 110    CONTINUE
 
        IF(NPNTU.LE.1) THEN
          WRITE(STDOUT,61)
          EFLAG = 1
        ENDIF
 
C       CHECK THE CROSS SECTION FOR NONSENSE BEHAVIOR AT THE END
 
        IF(ZU(1).LE.ZU(2).AND.XU(2).GT.XU(1)) THEN
C         THE LEFT MOST LINE SEGMENT HAS UPWARD SLOPE, THEREFORE HIGH
C         POINT IS NOT AT THE LEFT LIMIT.
 
          WRITE(STDOUT,'(/,A,A)') ' *WRN:502* UNEXPECTED SLOPE',
     A     ' AT LEFT END.'
          WRITE(STDOUT,*) '  SLOPE EXPECTED TO BE < 0  AT LEFT BOUNDARY'
        ENDIF
        IF(ZU(NPNTU).LE.ZU(NPNTU-1).AND.XU(NPNTU).GT.XU(NPNTU-1)) THEN
C         THE RIGHT MOST LINE SEGMENT HAS DOWNWARD SLOPE
 
          WRITE(STDOUT,'(/,A,A)') ' *WRN:503* UNEXPECTED SLOPE',
     A     ' AT RIGHT END.'
          WRITE(STDOUT,*) '  SLOPE EXPECTED TO BE > 0 AT RIGHT BOUNDARY'
        ENDIF
 
        IF(EXTEND.EQ.1) THEN
C         CHECK FOR ONE END BEING HIGHER THAN THE OTHER
 
          IF(ABS(ZMAX - ZU(1)).GT.EPSDIF) THEN
            WRITE(STDOUT,67) ZMAX - ZU(1)
          ENDIF
          IF(ABS(ZMAX - ZU(NPNTU)).GT.EPSDIF) THEN
            WRITE(STDOUT,68) ZMAX - ZU(NPNTU)
          ENDIF
        ELSE
C         CHECK FOR AN INTERMEDIATE POINT BEING HIGHER THAN EITHER
C         END
          IF(ABS(ZMAX - ZU(1)).GT.EPSDIF.AND.
     A       ABS(ZMAX - ZU(NPNTU)).GT.EPSDIF) THEN
            WRITE(STDOUT,69)
          ENDIF
 
          ZMAX = MIN(ZU(1), ZU(NPNTU))
        ENDIF
 
C       ASSIGN THE VALUES OF N FROM THE SUBSECTIONS TO THE LINE SEGMENT
C       LOCATIONS.
        DO 120 I=1,NPNTU-1
          LSNU(I) = NU(SBU(I))
 120    CONTINUE

C       Convert units here.
        DO 122 J=1,NPNTU 
          XU(J) = UFAC*XU(J)
          ZU(J) = UFAC*ZU(J)
122     CONTINUE
        ZMAX = UFAC*ZMAX
 
        IF(MODE.EQ.'DIRECT'.OR.MODE.EQ.'direct') THEN
 
C         COMPUTE ELEMENTS FOR CURRENT CROSS SECTION
 
          IF(NPNTU.GT.1) THEN
C           FIND THE MAXIMUM AND MINIMUM ARGUMENT VALUES
            ZMINU= 9999999.
            ZMAXU = -9999999.
            DO 150 J=1,NPNTU
              IF(ZU(J).GT.ZMAXU) ZMAXU = ZU(J)
              IF(ZU(J).LT.ZMINU) ZMINU = ZU(J)
 150        CONTINUE
 
            IF(EXTEND.EQ.0) ZMAXU = ZMAX
 
            IF(EFLAG.EQ.0) THEN
 
              CALL CXSTAB
     I                   (STDOUT, NSUBU, NAVMU, NFAC, MXPNTU, LEFT,
     I                    RIGHT, BETOPT, SNFLGU, NVARU, NATYU, YATNU,
     I                    NNYU,
     M                    NPNTU, ZMINU, ZMAXU, XU, ZU, SBU, EFLAG, LSNU,
     M                    SNU,
     O                    NU, NDEPU, XSTU)
            ENDIF
          ENDIF
 
C         OUTPUT THE TABLE IF NO ERRORS AND IF OUTPUT IS REQUESTED
 
          IF(EFLAG.NE.0.OR.TABU.EQ.0) GOTO 200
            IF(NOCM.EQ.1) THEN
              MESG = 0
            ELSE
              MESG = 1
            ENDIF
            STATT = STATU/SFAC
C           Change the generated table number to a tabid and get
C           an internal table number.
            WRITE(TABID,84) TABU
            CALL STRIP_L_BLANKS(
     M                          TABID)
            CALL GET_INTERNAL_TAB_NUMBER
     I                                  (STDOUT, TABID,
     M                                   EFLAG,
     O                                   ITEMP)
c           compute derivatives of square root of conveyance, alpha, beta, 
c           da, and dq

            call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

            CALL TABOUT
     I                 (STDOUT, STDTAB, ITEMP, STATT, ZMINU, MESG,
     I                  SAVOPT, OUTOPT, BETOPT, 
     i                  zone, hgrid, vdatum, unitsys, basis,
     i                    khflag, alphaflag, betaflag, maflag, mqflag,
     M                  NDEPU, XSTU, FTP)
 200      CONTINUE
        ELSE
C         OUTPUT TO A FILE ATTACHED TO UNIT STD49
 
          WRITE(STD49,77)
 
          IF(EXTEND.EQ.0) THEN
            WRITE(STD49,75) TABU, MONFLG, BETOPT, SAVOPT, OUTOPT
          ELSE
            WRITE(STD49,76) TABU, MONFLG, BETOPT, SAVOPT, OUTOPT
          ENDIF
 
          TP = STATU/SFAC
          CALL MKFMT
     I              (TP, 10,
     O               F10X)
          STAOUT(13:18) = F10X//')'
          WRITE(STD49,STAOUT) TP
          write(std49,86)  zone, hgrid, vdatum, unitsys, basis,
     a                    easting, northing

          WRITE(STD49,79)
    
          WRITE(STD49,80) NSUBU, (NU(J),J=1,MIN(6,NSUBU))
          DO 290 JJ=7,NSUBU,6
            LIM = NSUBU - JJ
            IF(LIM.GT.5) LIM = 5
            WRITE(STD49,83) (NU(JJ+J), J=0,LIM)
290       CONTINUE

C          WRITE(STD49,80) NSUBU, (NU(J), J=1,NSUBU)
          WRITE(STD49,81) SECID
 
          DO 300 J=1,NPNTU
            CALL VAR_DECIMAL(XU(J),
     O                       COFF)
            CALL VAR_DECIMAL(ZU(J),
     O                       CELEV)
            WRITE(STD49,82) COFF, CELEV, SBU(J)
 300      CONTINUE
          WRITE(STD49,50)
        ENDIF
 
        TABU = TABU + TABINC
        GOTO 100
 
      END
