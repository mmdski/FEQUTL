C  ***********************************************************************
C  *  Warning:  This program is large and complex and  extensive         *  
C  *  knowledge of its design, purpose, and limitations is required      *  
C  *  in order to apply it properly.  Application of this program by an  *
C  *  unqualified user for any other purpose than an educational one is  *
C  *  not only unwise but is also unethical.  The user of this           *
C  *  program is totally responsible for its use and application and for *
C  *  any actions or events which follow therefrom.  Any user of this    *
C  *  program  holds the developer of the program harmless from          *
C  *  damages of any kind.                                               *
C  *                                                                     *       
C  *  The developer has used reasonable care in the construction and     *
C  *  testing of the program.  However, in a program of this size and    *
C  *  complexity, it is impossible to verify more than a minute number of*
C  *  possible options or applications.  The developer is continuing to  *
C  *  modify and use the program and is interested in information on     *
C  *  operational problems encountered in its application.  However, the *
C  *  developer gives no assurance that the problem can or will be       *
C  *  rectified.                                                         *
C  *                                                                     *       
C  *  This program is not to be sold in any form modified or otherwise.  *
C  ***********************************************************************

C
C
C
      SUBROUTINE CLEAR_GHOME()

C     Clear the GHOME part of file names. 

      IMPLICIT NONE

      INCLUDE 'home.cmn'
C***********************************************************************
      GHOME = ' '
      RETURN
      END

C
C
C
      SUBROUTINE MAKE_STANDARD_FILE_NAMES(FNAME,
     O                                    FNAME2, FNAME3)

C     Given the name for the input file to FEQUTL form the standard names for
C     the remaining file names. 

      IMPLICIT NONE

      CHARACTER*64 FNAME, FNAME2, FNAME3


C     Local

      INTEGER I, L, N

C***********************************************************************
C     The input file name may or may not have an extension.  Also
C     this should work with more than one period in the file name. 
C     The period closest to the end of the string will be taken to 
C     be the delimiter for the extension.  If no period is found, 
C     then the whole name is used for the base name of the remaining
C     file names unless the period is in the first position and is
C     the only period present. 

      N = LEN_TRIM(FNAME)

C     Set L for case of no period found
      L = N
      DO I=N,1,-1

        IF(FNAME(I:I).EQ.'.') THEN
          IF(I.GT.1) THEN
            L = I - 1
          ELSE
            L = N
          ENDIF
          EXIT
        ENDIF

      END DO

      FNAME2 = FNAME(1:L)//'.out'
      FNAME3 = FNAME(1:L)//'.tab'
      RETURN
      END
        


C
C
C
      INTEGER FUNCTION   GETUSB
     I                         (IBN)
 
C     + + + PURPOSE + + +
C     Dummy function for calls from fequtl to routines which
C     refer to a branch number.
      IMPLICIT NONE

 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER IBN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IBN    - Internal branch number
C***********************************************************************
      GETUSB = IBN
      RETURN
      END
C
C
C
      SUBROUTINE   INIT()
 
C     + + + PURPOSE + + +
C     Initialize common variables
 
      IMPLICIT NONE

C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomd.cmn'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'xscomb.cmn'
      INCLUDE 'fldway.cmn'
      INCLUDE 'bridge.cmn'
      INCLUDE 'flotab.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
      INTEGER DVEC(XSCOML), UVEC(XSCOML)
 
C     + + + EQUIVALENCES + + +
      EQUIVALENCE (NPNTD,DVEC(1)), (NPNTU,UVEC(1))
C***********************************************************************
c     zero come common blocks
      DO 100 J=1,XSCOML
        DVEC(J)=0
        UVEC(j) = 0 !added  5 Sept 2007 to avoid uninitialized var errors 
 100    CONTINUE
      MAXNFT=20
      MXPNTU=PMXPNT
      MXPNTD=PMXPNT
      MXPNTB=PMXPNT
      NPNTU=0
      MAXNPZ=25
C     SET FLOODWAY FLAG TO OFF. MUST BE TURNED ON BY USER TO
C     BE ACTIVE
      FLOOD = 0
c     Zero lsnu and lsnd to avoid failures when debugging with
c     full checking turned on using LF95
      lsnu = 0.0
      lsnd = 0.0
      RETURN
      END
C
C
C
      SUBROUTINE   MKFMT
     I                  (XA, W,
     O                   FMSTRG)
 
C     + + + PURPOSE + + +
C     Make a fixed format code, excluding () where W is the field
C     width, and X is a defining number.  The goal is to provide
C     as many decimal places as will fit, leaving room for a
C     decimal point and possibly a sign.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER W
      REAL XA
      CHARACTER FMSTRG*(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     XA     - Trial number to define a format
C     W      - Field with in characters for the format
C     FMSTRG - Format string
 
C     + + + LOCAL VARIABLES + + +
      INTEGER D
      REAL X
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, INT, LOG10
C***********************************************************************
      FMSTRG = ' '
      X = ABS(XA)
      IF(X.NE.0.0) THEN
        D = W - (INT(LOG10(X)) + 3)
        IF(D.GT.W-3) THEN
          D = W - 3
        ELSEIF(D.LT.0) THEN
          D = 0
        ENDIF
      ELSE
        D = W - 3
      ENDIF
 
      WRITE(FMSTRG,'(1HF,I2,1H.,I1)') W, D
      RETURN
      END
C
C
C
      CHARACTER*7 FUNCTION   PUT7
     I                           (X)
 
C     + + + PURPOSE + + +
C     Function to convert a real number into a special compact
C     form of output to retain 4 significant figures for
C     numbers in the range -1e12 < x < 1e12.  This range includes
C     all reasonable flows for any river on earth!
C     The smallest non-zero flow is 1e-6, again smaller than any
C     flow of interest in a stream.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      REAL X
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     X      - Value to be put
 
C     + + + LOCAL VARIABLES + + +
      CHARACTER RESULT*7, WORK*10
C***********************************************************************
      WRITE(WORK,'(4PE10.3)') X
      IF(WORK(8:8).EQ.'+'.AND.WORK(9:9).NE.'0') THEN
C       OVERFLOW
        RESULT = ' ******'
      ELSEIF(WORK(8:8).EQ.'-'.AND.WORK(9:9).NE.'0') THEN
C       UNDERFLOW
        RESULT = ' 0000+0'
      ELSE
        RESULT(1:5) = WORK(1:5)
        RESULT(6:6) = WORK(8:8)
        RESULT(7:7) = WORK(10:10)
      ENDIF
      PUT7 = RESULT
      RETURN
      END
C
C
C
      INTEGER FUNCTION   STRLEN
     I                         (STR)
 
C     + + + PURPOSE + + +
C     Return the actual length of the character array,
C     excluding trailing blanks.  A string of all blanks is taken
C     to be empty.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      CHARACTER STR*(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STR    - string for length computation
 
C     + + + LOCAL VARIABLES + + +
      INTEGER DONE, LENT
 
C     + + + INTRINSICS + + +
      INTRINSIC LEN
C***********************************************************************
      LENT = LEN(STR)
      DONE = 0
 10   CONTINUE
        IF (STR(LENT:LENT).EQ.' ') THEN
          LENT = LENT - 1
        ELSE
          DONE = 1
        END IF
      IF (LENT.GT.0.AND.DONE.EQ.0) GO TO 10
 
      STRLEN= LENT
 
      RETURN
      END
C
C
C
      SUBROUTINE   TABCHK
     I                   (STDOUT, PMXTAB,
     M                    TAB, TABDIR, EFLAG)
 
C     + + + PURPOSE + + +
C     Check for valid table number
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, PMXTAB, STDOUT, TAB
      INTEGER TABDIR(PMXTAB)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     PMXTAB - Maximum value of function table number
C     TAB    - Table number
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT('0*ERR:510* Duplicate table id.')
C***********************************************************************
      IF(TABDIR(TAB).EQ.0) GOTO 110
        WRITE(STDOUT,52)
        EFLAG = EFLAG + 1
        RETURN
 110  CONTINUE
        TABDIR(TAB) = TAB
        RETURN
      END
C
C
C
      BLOCKDATA   XOFFIN
 
C     + + + PURPOSE + + +
C     Initialize the offset list for cross sections.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'offcom.cmn'
 
C     + + + DATA INITIALIZATIONS + + +
      DATA OFFVEC/6,10*0,8,7*0,5, 6, 8, 7, 8, 10,
     A      4*0, 7, 8, 11, 11, 12, 15/
C***********************************************************************
      END
C
C
C
      SUBROUTINE   CHKCFC
     I                   (GRAV, STDOUT, ADRS,
     O                    WFLAG)
 
C     + + + PURPOSE + + +
C     Check the critical flow and celerity in the table given by
C     the address, ADRS, and report regions of decrease in
C     critical flow or celerity.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ADRS, STDOUT, WFLAG
      REAL GRAV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     GRAV   - value of acceleration due to gravity
C     STDOUT - Fortran unit number for user output and messages
C     ADRS   - Address of function table
C     WFLAG  - Flag to suppress more than a single warning message
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'offcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HA, IOFF, TYPE, XOFF, doff, dtype
      REAL A, C, P, Q, QOLD, T, Y, MXSLOT
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, SQRT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *BUG:XXX.  Invalid type in CHKCFC') 
 52   FORMAT(' *WRN:521* Critical flow decreases by ',F7.2,
     A  ' per cent at depth=',F10.3)
 54   FORMAT(' Table tests OK.')
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        MXSLOT = 0.07
      ELSE
        MXSLOT = 0.02134
      ENDIF
      TYPE = ITAB(ADRS+2)
      XOFF = OFFVEC(TYPE)
      HA = ITAB(ADRS)
      doff = itab(adrs+21)
      if(doff.gt.0) then
        dtype = 10
      else
        dtype = 0
      endif
 
      type = type - dtype
      IF(TYPE.LT.20.OR.TYPE.GT.25) THEN
        WRITE(STDOUT,50) TYPE + dtype
        STOP 'Abnormal stop. Bug found.'
      ENDIF
      IF(TYPE.EQ.20.OR.TYPE.EQ.21.OR.TYPE.EQ.23.OR.TYPE.EQ.24) THEN
C       IOFF GIVES THE OFFSET FROM THE TABLE ADDRESS TO THE FIRST
C       POSITIVE DEPTH  ENTRY IN THE TABLE.
 
        IOFF = XTIOFF + XOFF
 
C        COLD = 0.0
        QOLD = 0.0
        WFLAG = 0
 
 100    CONTINUE
 
          Y = FTAB(ADRS + IOFF)
          T = FTAB(ADRS + IOFF + 1)
          A = FTAB(ADRS + IOFF + 2)
 
          C = SQRT(GRAV*A/T)
          Q = A*C
 
C          IF(C.LE.COLD) THEN
C            P = 100.*ABS(C - COLD)/COLD
C            WRITE(STDOUT,50) P, Y
C            WFLAG = 1
C          ENDIF
 
C          COLD = C
 
          IF(Q.LE.QOLD.AND.T.GT.MXSLOT) THEN
            P = 100.*ABS(Q - QOLD)/QOLD
            WRITE(STDOUT,52) P, Y
            WFLAG = 1
          ENDIF
 
          QOLD = Q
          IOFF = IOFF + XOFF
 
          IF(ADRS + IOFF.LE.HA) GOTO 100
      ELSE
C       IOFF GIVES THE OFFSET FROM THE TABLE ADDRESS TO THE FIRST
C       POSITIVE DEPTH  ENTRY IN THE TABLE.
 
        IOFF = XTIOFF + XOFF
 
C        COLD = 0.0
        QOLD = 0.0
        WFLAG = 0
 
 200    CONTINUE
 
          Y = FTAB(ADRS + IOFF)
          T = FTAB(ADRS + IOFF + 1)
          A = FTAB(ADRS + IOFF + 2)
          Q = FTAB(ADRS + IOFF + 7)
          C = Q/A
 
C          IF(C.LE.COLD) THEN
C            P = 100.*ABS(C - COLD)/COLD
C            WRITE(STDOUT,50) P, Y
C            WFLAG = 1
C          ENDIF
 
C          COLD = C
 
          IF(Q.LE.QOLD.AND.T.GT.0.07) THEN
            P = 100.*ABS(Q - QOLD)/QOLD
            WRITE(STDOUT,52) P, Y
            WFLAG = 1
          ENDIF
 
          QOLD = Q
          IOFF = IOFF + XOFF
 
          IF(ADRS + IOFF.LE.HA) GOTO 200
      ENDIF
      IF(WFLAG.EQ.0) WRITE(STDOUT,54)
      RETURN
      END
C
C
C
      SUBROUTINE   NXTCMD
     I                   (STDIN, STDOUT, NCMD, CMDTAB, CMDVAL,
     O                    NEXT)
 
C     + + + PURPOSE + + +
C     Read next command, skipping blank commands as needed, and
C     determine action by defining NEXT.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NCMD, NEXT, STDIN, STDOUT
      INTEGER CMDVAL(NCMD)
      CHARACTER CMDTAB(NCMD)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     NCMD   - Number of commands
C     CMDTAB - Vector holding the command names
C     CMDVAL - Vector holding the command value for branching
C     NEXT   - Code for the next command
 
C     + + + SAVED VALUES + + +
      CHARACTER BLANK*8
      SAVE BLANK
 
C     + + + LOCAL VARIABLES + + +
      INTEGER INDEX
      CHARACTER CARD(10)*8, LINE*80
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BINSER, inline
 
C     + + + DATA INITIALIZATIONS + + +
      DATA BLANK/'        '/
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(10A8)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' ',10A8)
 52   FORMAT(/,' End of file before FINISH found.  FINISH supplied.')
C***********************************************************************
 100  CONTINUE
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
        READ(LINE,1) CARD
        WRITE(STDOUT,51) CARD
        IF(CARD(1).EQ.BLANK) GOTO 100
 
      IF(CARD(1).EQ.'ENDFILE') THEN
C       End of file found before next command.  Missing FINISH.
        WRITE(STDOUT,52)
        CARD(1) = 'FINISH'
      ENDIF
C     FIND COMMAND IN CMDTAB
 
      CALL BINSER
     I           (CARD(1), NCMD, CMDTAB,
     O            INDEX)
 
      IF(INDEX.LE.0) THEN
        NEXT = 0
      ELSEIF(INDEX.LE.NCMD) THEN
        NEXT = CMDVAL(INDEX)
      ELSE
        WRITE(STDOUT,*) ' ADDRESS PROBLEM IN NXTCMD'
        WRITE(STDOUT,*) ' INDEX=',INDEX
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
      RETURN
      END
C
C
C
      SUBROUTINE   CHKTAB
     I                   (CTYPE, LOUT, FTPNT, MFT,
     M                    TAB,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Checks for address of a table and returns its address in place of
C     TAB.  Otherwise writes error and sets TAB = 1 on return.
C     IF TAB = 0 on entry return without doing anything.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CTYPE, EFLAG, LOUT, MFT, TAB
      INTEGER FTPNT(MFT)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     CTYPE  - type class number for checking valid table types
C     LOUT   - Fortran unit number for user output and messages
C     FTPNT  - function table pointer giving the table address for each
C              table number.  If the address is zero the table does not
C              exist.
C     MFT    - Maximum allowed table number
C     TAB    - Table number
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     Common blocks


C     + + + LOCAL VARIABLES + + +
      INTEGER ADRS
      CHARACTER CHAR16*16
 

C     External Functions

      INTEGER GETTBN

C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      CHARACTER GET_TABID*16
      EXTERNAL CHKTYP, GETTBN, GET_TABID, LENSTR
 
C     + + + OUTPUT FORMATS + + +
 1    FORMAT(/,' *ERR:81* TABID = ',A,' does not exist.')
 2    FORMAT(/,' *ERR:82* TABLE# = ',I5,' too large. Set to current',
     1       ' maximum of:',I5)
50    FORMAT(/,' *BUG:XXX* Table number mismatch: External=',I5,
     A       ' Internal=',I5,' Address=',I10)
C***********************************************************************
      IF(TAB.EQ.0) RETURN
 
C     MAKE SURE TABLE NUMBER IS IN VALID RANGE
 
      IF(TAB.GT.MFT) THEN
        WRITE(LOUT,2) TAB, MFT
        TAB = MFT
        EFLAG = 1
      ENDIF
 
      ADRS = FTPNT(TAB)
      
      IF(ADRS.GT.0) GOTO 10
        CHAR16 = GET_TABID(TAB)
        WRITE(LOUT,1) CHAR16(1:LENSTR(CHAR16))
        EFLAG = 1
        TAB = 1
        RETURN
 10   CONTINUE
      IF(TAB.NE.GETTBN(ADRS)) THEN
        WRITE(LOUT,50) TAB, GETTBN(ADRS), ADRS
        STOP ' Abnormal stop. Bug found.' 
      ENDIF      
      TAB = ADRS
      CALL CHKTYP
     I           (LOUT, ADRS, CTYPE,
     O            EFLAG)
      RETURN
      END
C
C
C
      SUBROUTINE   INSYS
     I                  (STDSYS, MAXCMD, NARG, STDOUT,
     O                   UNITS, NFAC, GRAV,
     O                   MAXKNT, EPS, NCMD, CMDTAB, CMDVAL)
 
C     + + + PURPOSE + + +
C     Input the default control values for the program.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER MAXCMD, MAXKNT, NARG, NCMD, STDIN, STDOUT, STDSYS, STDTAB
      INTEGER CMDVAL(MAXCMD)
      REAL EPS, GRAV, NFAC
      CHARACTER CMDTAB(MAXCMD)*8, UNITS*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDSYS - Fortran unit number for system table definition
C     MAXCMD - Maximum number of commands
C     NARG   - Number of command line arguments found
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     UNITS  - Definition of the unit system: ENGLISH or METRIC
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     GRAV   - value of acceleration due to gravity
C     MAXKNT - Iteration count limit
C     EPS    - Convergence tolerance
C     NCMD   - Number of commands
C     CMDTAB - Vector holding the command names
C     CMDVAL - Vector holding the command value for branching
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'gnicom.cmn'
 
C     + + + SAVED VALUES + + +
      CHARACTER ENGL*8, METRIC*8
      SAVE ENGL, METRIC
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
      REAL LATITUDE, ELEVATION
      CHARACTER DUMMY*8, NFAC_FLAG*8, LINE*80
 
C     + + + EXTERNAL FUNCTIONS + + +
      CHARACTER GETTOK*8
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL GETTOK, GRULE, SORT2, STRIP_L_BLANKS, inline
 
C     + + + DATA INITIALIZATIONS + + +
      DATA METRIC/'METRIC'/,ENGL/'ENGLISH'/
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(8X,2I5)
 2    FORMAT(A8,I5)
 3    FORMAT(8X,A8)
 4    FORMAT(8X,A8,1X,A8,F10.0,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT('0*ERR:503* TOO MANY COMMANDS IN SYSTEM FILE. LIMIT=',
     1        I5)
 52   FORMAT('0*ERR:613* ',A8,' IS AN UNKNOWN UNITS OPTION.')
54    FORMAT(/,' FEQUTL uses:',
     A     /,5X,' NFAC=',F10.6,' for Manning''s equation.',
     A     /,5X,' GRAV=',F10.4,' for g.')
56    FORMAT(/,' NFAC and GRAV are nominal engineering values')
58    FORMAT(/,' NFAC is exact.  GRAV is for latitude ',
     A      F5.1,' degrees and elevation ',F6.0)
C***********************************************************************
C     SKIP THE NEXT THREE CARDS OF INFORMATION

      CALL inline
     I          (STDSYS, STDOUT,
     O           LINE)

      IF(LINE(1:5).EQ.'STDIN') THEN
C       Read two more lines to skip over obselete information
        CALL inline
     I            (STDSYS, STDOUT,
     O             LINE)
        CALL inline
     I            (STDSYS, STDOUT,
     O             LINE)
      ELSE
         BACKSPACE(STDSYS)
      ENDIF
 
      CALL inline
     I          (STDSYS, STDOUT,
     O           LINE)
      READ(LINE,4) UNITS, NFAC_FLAG, LATITUDE, ELEVATION
      LATITUDE = LATITUDE*1.7453293E-2
      UNITS = GETTOK(UNITS)
 
      CALL STRIP_L_BLANKS(
     M                    NFAC_FLAG)

      IF(UNITS.EQ.ENGL) THEN
        IF(NFAC_FLAG.EQ.' '.OR.NFAC_FLAG.EQ.'NOMINAL') THEN
C         Use the approximate value
          NFAC = 1.49
          GRAV = 32.2
        ELSE
C         Use the "exact" value
          NFAC = 1.485919
          GRAV = 32.1726 - 0.08495*COS(2.*LATITUDE) 
     A                         - 3.09E-6*ELEVATION
        ENDIF
      ELSEIF(UNITS.EQ.METRIC) THEN
        NFAC = 1.
        IF(NFAC_FLAG.EQ.' '.OR.NFAC_FLAG.EQ.'NOMINAL') THEN
C         Make GRAV consistent with 32.2! 
          GRAV = 9.81456
        ELSE
          GRAV = 32.1726 - 0.08495*COS(2.*LATITUDE) 
     A                      - 3.09E-6*ELEVATION/0.3048
          GRAV = GRAV*0.3048
        ENDIF
      ELSE
        WRITE(STDOUT,52) UNITS
        STOP 'Abnormal stop. Errors found.'
      ENDIF

      WRITE(STDOUT,54) NFAC, GRAV
 
      IF(NFAC_FLAG.EQ.' '.OR.NFAC_FLAG.EQ.'NOMINAL') THEN
        WRITE(STDOUT,56)
      ELSE
        WRITE(STDOUT,58) LATITUDE*57.29578, ELEVATION
      ENDIF
      CALL inline
     I          (STDSYS, STDOUT,
     O           LINE)
      READ(LINE,1) NCMD, NGS
      IF(NGS.EQ.0) NGS = 5
 
      IF(NCMD.LE.MAXCMD) GOTO 90
        WRITE(STDOUT,51) MAXCMD
        STOP 'Abnormal stop. Errors found.'
 90   CONTINUE
 
      DO 100 J=1,NCMD
        CALL inline
     I            (STDSYS, STDOUT,
     O             LINE)
        READ(LINE,2) CMDTAB(J),CMDVAL(J)
        CMDTAB(J) = GETTOK(CMDTAB(J))
 100  CONTINUE
 
C     SORT THE COMMANDS IN ASCENDING ORDER
 
      CALL SORT2
     I          (NCMD,
     M           CMDTAB, CMDVAL)
 
      MAXKNT = 10
      EPS = 0.005
 
C     COMPUTE THE RULE FOR GAUSSIAN INTEGRATION
 
      CALL GRULE
     I          (NGS,
     O           XGS, WGS)
 
      RETURN
      END
c
c 
c 
      subroutine get_single_named_item(stdout, line,
     o                                 name, item, eflag)

c     Get a single named item on the line 

      implicit none

      integer stdout, eflag

      character*(*) line, name, item

c     Local
      integer ieq, iend

c     ***************************formats********************************
50    format(/,'*ERR:XXX* Equal sign missing after name=',a)
c***********************************************************************
      ieq = index(line, '=')
      if(ieq == 0) then
c       Missing delimiter
        name = adjustl(line)
        iend = len_trim(name)
        write(stdout,50) name(1:iend)
        eflag = 1
      else
        name = adjustl(line(1:ieq-1))
        item = adjustl(line(ieq+1:))
      endif

      return
      end

C
C
C
      SUBROUTINE SET_VERSION()

C     Sets the current version number and date in the version 
C     common block

      INCLUDE 'version.cmn'
C***********************************************************************
      VERSION_NUMBER = 5.80
      VERSION_DATE = '6 October 2008'
      RETURN
      END
C
C
C
      PROGRAM   FEQUTL
 
C     + + + PURPOSE + + +
 
      IMPLICIT NONE
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'morg.prm'
      INCLUDE 'stdun.cmn'
      INCLUDE 'fldway.cmn'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'ftable.cmn'
      INCLUDE 'xscom.cmn'
      INCLUDE 'grvcom.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'clcom.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'julian.cmn'
      INCLUDE 'tabid.cmn'
      INCLUDE 'version.cmn'
      include 'tabupgrade.cmn'
      include 'home.cmn'
      include 'whatos.cmn'
      include 'datetime.cmn'
      include 'grid_datum.cmn'
      include 'maketabindex.cmn'
      include 'svn.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CFLAG, FTKNT, FTP, GFLAG, IS, ITABA, J, MAXKNT, MODE,
     A        NCMD, NEXT, STDFLD, STDIN, STDSYS, STDTAB, TABLE,
     B        STDEXT_OPTION, N, IOS, conf_flag, ieq, stdscr
   
      INTEGER CMDVAL(MAXCMD), TABDIR(PMXTAB+5), VALUES(8)
      REAL EPS, HSLOT, NFAC, MINQ, velocity
      CHARACTER BOT*8, CMDTAB(MAXCMD)*8, FHEAD*91, FNAME*64, LEFT*8,
     A          LINE*80, NXTNAM*128, REPLY*3, RIGHT*8, UNITS*8,
     B          DATE*8, ZEIT*10, ZONE*5, TABID*16, FNAME2*64,
     C          FNAME3*64, GISID*16, cmd_line_args(0:10)*64,
     d          conf_file*64, item_name*8, current_file*128

      LOGICAL THERE, trk_files
 
C     + + + INTRINSICS + + +
      INTRINSIC IABS, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER IARGC, GET_UNIT, what_os
      CHARACTER*16 GET_TABID, GET_GISID
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHANEL, CHNTAB, CRITQ, CULVRT, EMBANK, EXPCON, FEQX,
     A         FLDIN, FQXE, FTABIN, GETARG, HEC2X, IARGC, INIT, inline,
     B         INSYS, KIL, MULCON, NXTCMD, PIPES, QCLIM, RITTER, SEWER,
     C         SPBRID, UFGATE, WPRO14, WPROQZ, WPROX, XSTMAK,
     D         GET_UNIT, FREE_UNIT, RESET_KOUNT_OF_INTERNAL_TABIDS,
     E         GET_TABID, GET_GISID, os_file_style, what_os
 
C     + + + INPUT FORMATS + + +
 40   FORMAT(7X,F10.0)
 41   FORMAT(7X,A3)
 42   FORMAT(9X,A3)
 44   FORMAT(7X,F10.0)
 45   FORMAT(A64)
 46   FORMAT(5X,F10.0,1X,7X,F10.0)
 49   FORMAT(6X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT('0*ERR:501* Unknown command.')
 52   FORMAT(' *NOTE* Known command found after error.')
 54   FORMAT('0*WRN:518* One or more errors found. Table file is ',
     A       'incomplete or invalid.')
56    format(/,' Global values: ZONE=',A8,' HGRID=',A8,' VDATUM=',a8,
     a        ' UNITSYS=',a8, ' BASIS=',a8)
58    format(16x, 'VDATUM_SHIFT=',f10.4)
 70   FORMAT(' DZLIM=',F10.6,' NRZERO=',F10.6,' USGSBETA=',A3)
 72   FORMAT(' Version: ',F6.2,' Version date: ',A17,5X,
     A          'Date/time of run: ',A4,'/',A2,'/',A2,': ',
     B          A2,'.',A2,'.',A6)
 74   FORMAT(' *WRN:526* INVALID RESPONSE FOR USGSBETA. TAKEN AS:NO')
 76   FORMAT(' EPSARG=',1PE10.4,/,' EPSF=',1PE10.4,' EPSABS=',1PE10.4)
 77   FORMAT('*WRN:558 Vertical extension of a cross section end in',
     A   ' order to match',/,10X,' the other end''s elevation may lead',
     B   ' to nonsense results.')
 78   FORMAT(/,' Automatic cross section extension to match end',
     A  ' elevations selected.')
 79   FORMAT(/,' No automatic cross section extension.  Cross section',
     A  /,5X,' table arguments limited to minimum end elevation.')
 80   FORMAT(' *WRN:559* Invalid response for EXTEND.  Taken as: NO')
 81   FORMAT(/,' EXTEND, the cross section extension option, missing.',
     A  ' Extension assumed.',/,4x,'EXTEND appears after ESPF in the',
     B  ' input.  Options are:',/,4X,'EXTEND=NO- cross section table',
     C  ' arguments limited to minimum',/,4X,'of the two end point',
     D  ' elevations.',/,4X,'EXTEND=YES- lower elevation end point',
     E  ' extended to higher elevation.',/,4X,'Previous versions',
     D  ' used EXTEND=YES by default.')
82    format(/,' Global home name is: ',a)
83    format(/,' FEQUTL configuration from conf file: ',a) 
84    format(/,' FEQUTL configuration from master-input file.')
85    format(/,' Master-input file: ',a,/,
     a         ' Master-output file: ',a,/,
     b         ' Function-table file: ',a)

86    format(/,';Version number is ',i8,' for the Subversion working ',
     a         'copy') 
87    format(/,';*WRN:XXX* The working copy contains local ',
     a 'modifications.',/,';',5x,
     b 'The version number given is therefore',
     b' not current and may not properly',/,';',5x,
     b  'provide the version',
     c' number required to retrieve all the same files ',
     d'at a later time.')
88    format('; Global home name is: ',a)
89    format('; Local home name is: ',a)

 90   FORMAT('  Processing FTABIN')
91    FORMAT(/,' EPSABS is missing or zero.  Value set to EPSF.',/,
     A    ' Check input and add EPSABS after EPSF on the same line.')
92    FORMAT(/,' Minimum flow target for EMBANKQ and CHANRAT=',F10.3)
 93   FORMAT(1X,A16,1X,A16,4X,A4,F13.2,4X,A8,2X,A8,2X,A8,F10.3,
     a f10.1,f10.2)
94    FORMAT(/,' Table Id-------- GISID-----------  FldOpt',
     A ' BF Elevation  FEQ Invert Left----- Right---- Loss-----',
     a ' FldwyArea FldwyVel-') 
95    FORMAT(/,' Invalid number of command-line agruments: FEQUTL',
     A/,' expects exactly one argument or exactly three arguments.',
     B/,' If one argument is given, FEQUTL strips the last extension,',
     C/,' if there is one, and appends .out and .tab respectivly to',
     D/,' create the second and third file names.  Otherwise give',
     E/,' three file names: (1) user-input file, (2) user-output file,',
     D/,' and (3) function-table file.')
96    FORMAT(/,' The operating system cannot open the user-output',
     A/,' file: ',A,'.','  Check for invalid characters in the',
     B/,' name.  If part or all of the path is given with the name,',
     C/,' make sure that all directories exist as spelled.')
97    FORMAT(/,' The operating system cannot open the function-table',
     A/,' file: ',A,'.','  Check for invalid characters in the',
     B/,' name.  If part or all of the path is given with the name,',
     C/,' make sure that all directories exist as spelled.')
C***********************************************************************
      trk_files = .false. !Disable file tracking.  Implementation incomplete
C
C     strings for Unix what command
      line =
     &  '@(#)FEQUTL - Full Equations Flow Routing Model Utility Program'
      line =
     &      '@(#)FEQUTL - Franz, D.D., and Melching, C.S., WRIR 97-4037'
      line = '@(#)FEQUTL - Contact: h2osoft@usgs.gov'
      line = '@(#)FEQUTL - Version: 4.81x 1998/05/22'
C     set string for use with RCS ident command
      line =
     &'$Id: fequtl.f,v 4.4 1996/02/28 19:16:54 rsregan Exp rsregan $'
C
C     Create the mechanism for assignment of I/O unit numbers

      CALL INITIALIZE_UNITS
      
C     The argument to GET_UNIT intends to write to the console should
C     an error arise in assigning unit numbers.  May not work on
C     all systems!

      STD5 = GET_UNIT(0)
      STD6 = GET_UNIT(0)
      STD7 = GET_UNIT(0)
      STD10 = GET_UNIT(0)
      STD48 = GET_UNIT(0)
      STD49 = GET_UNIT(0)
      STD50 = GET_UNIT(0)
      stdscr = get_unit(0) !define a unit for a scratch file that mirrors the master-input file
      STDIN = STD5
      STDOUT = STD6
      STDSYS = STDIN
      STDTAB = STD7
 
      STDFLD = STD10
      OUTPUT = 1
      IN = STDIN
      LOUT = STDOUT
      SLOT = 1.E30
      TIME = 0.0D0
      FTP = 1
      FTKNT = 0
      GFLAG = 0
      NARG = 0
      NXT = 1


c     set options for testing table fitting
      ty13_to_ty43 = 'NO'
      upgrade_xsec_tab = 'NO'

c     Clear the location-status flag for function tables.
      ft_loc_status = ' '


C     Clear the home-name values
      CALL CLEAR_HOME()
      CALL CLEAR_GHOME()

C     Set the code for the carriage-return character
      CALL SETCR()

C     Clear the base julian time for time series tables.  Not needed
C     now but will be needed if time series tables are added to FEQUTL.
      TAB_789_JTBASE = 0.D0

C     Clear the count of internal table ids
      CALL RESET_KOUNT_OF_INTERNAL_TABIDS()

C     Clear the symbol table for processing table id's
      NUM_TABID = 0

c     set the Operating System flag
      osis = what_os()

C     PROCESS COMMAND LINE ARGUMENTS
 
      NARG = IARGC()
c     Bring all command-line arguments into a vector and then scan for special
c     arguments, such as the standard-header file.  Process the special arguments,
c     if any, then set things so that the subsequent code will work as it did 
c     prior to the addition of special arguments. 
      if(narg.gt.10) then
        write(*,*) ' narg > 10. Too many command-line arguments.'
        STOP 'Abnormal stop. Errors found.'
      else
        do j=1,narg - MORG
          CALL GETARG
     I               (j+MORG,
     O                cmd_line_args(j)  )
        end do

c       Now scan for special arguments:  These are:
c        -conf  signals that a file name for the fequtl configuration
c               file is given.  

        conf_flag = 0
        do j=1,narg
          if(cmd_line_args(j)(1:5) == '-conf') then
c           found configuration-file name.  The file name
c           is the following argument.
            conf_file = cmd_line_args(j+1)
            conf_flag = 1
            call os_file_style(
     m                                    conf_file)
c           reset stdsys to a unique unit number
            stdsys = get_unit(0)
c           Adjust narg to match what would have been true had 
c           -conf not been present.
            narg = narg -2
            exit
          endif
        end do
      endif

      IF(NARG.EQ.1+MORG) THEN
C       If only one file argument is given we assume that the 
C       extension, if any, from the file given is stripped and the 
C       remaining two file names are formed by adding .out and .tab
C       to the file name given by the user name stripped of its extension.

        STDEXT_OPTION = 1
      ELSEIF(NARG.LT.3+MORG) THEN
        WRITE(*,95)
        STOP 'Abnormal stop. Errors found.'
      ELSE
        STDEXT_OPTION = 0
      ENDIF
 
 
C     GET THE FIRST FILE ARGUMENT
 
      CALL GETARG
     I           (1+MORG,
     O            FNAME)
c     fix common OS deviations in file names
      call os_file_style(
     m                           fname)
      INQUIRE(FILE=FNAME, EXIST=THERE)
      IF(THERE) THEN
        OPEN(STDIN, FILE = FNAME, STATUS = 'OLD')
      ELSE
        N = LEN_TRIM(FNAME)
        WRITE(*,*) ' '
        WRITE(*,*) ' File named: ',FNAME(1:N),' not found.'
        WRITE(*,*) ' Please check spelling of master-input file.'
        STOP 'Abnormal stop. Errors found.'
      ENDIF

      IF(STDEXT_OPTION.EQ.1) THEN
C       Form the other two file names 
        CALL MAKE_STANDARD_FILE_NAMES(FNAME,
     O                                FNAME2, FNAME3)
      ELSE
 
C       GET THE SECOND FILE ARGUMENT
        
        CALL GETARG
     I             (2+MORG,
     O              FNAME2)
c       fix common OS deviations in file names
        call os_file_style(
     m                             fname2)
        
C       GET THE THIRD FILE ARGUMENT
        
        CALL GETARG
     I             (3+MORG,
     O              FNAME3)
c       fix common OS deviations in file names
        call os_file_style(
     m                              fname3)

      ENDIF

C     Make sure that no file names match. 
      IF(FNAME.EQ.FNAME2.OR.FNAME.EQ.FNAME3.OR.FNAME2.EQ.FNAME3) THEN
C       One or more names match.
          WRITE(*,*) ' '
          WRITE(*,*) ' Two or more names given as command-line'
          WRITE(*,*)
     A ' arguments are the same.  All names must be unique.'
          STOP 'Abnormal stop: errors found.'
      ENDIF  
      OPEN(STDOUT, FILE = FNAME2, STATUS = 'UNKNOWN',IOSTAT=IOS)
      IF(IOS.NE.0) THEN
        N = LEN_TRIM(FNAME2)
        WRITE(*,96) FNAME2(1:N)
        STOP 'Abnormal stop. Errors found.'
      ENDIF

      OPEN (STDTAB, FILE=FNAME3, STATUS = 'UNKNOWN', IOSTAT=IOS)
      IF(IOS.NE.0) THEN
        N = LEN_TRIM(FNAME3)
        WRITE(*,97) FNAME3(1:N)
        STOP 'Abnormal stop. Errors found.'
      ENDIF

c     Make a copy of the master-input file as a scratch file here.  We may 
c     need it if we later want to compute md5 digests of various input blocks. 
c
      if(trk_files) then
        call copy_master_input_file(stdin, stdscr)
      endif

C     Clear the cross section slot values

      CALL CLEAR_SLOT(STDTAB)


      CALL SET_VERSION()
      CALL DATE_AND_TIME(DATE, ZEIT, ZONE, VALUES)
      WRITE(STDOUT,72) VERSION_NUMBER, VERSION_DATE, DATE(1:4), 
     A      DATE(5:6), DATE(7:8), ZEIT(1:2), ZEIT(3:4), ZEIT(5:10)
c     Save the version and run time stamps
      WRITE(version_run_date_time_string,72) 
     a      VERSION_NUMBER, VERSION_DATE, DATE(1:4), 
     b      DATE(5:6), DATE(7:8), ZEIT(1:2), ZEIT(3:4), ZEIT(5:10)
      
 
      CALL WHAT_EXECUTABLE(STDOUT)


      write(stdout,85) fname(1:len_trim(fname)),
     a                 fname2(1:len_trim(fname2)),
     b                 fname3(1:len_trim(fname3))

c     Write the version/run string to  the standard table file.  Make an echoing comment.
      write(stdtab,'(a)') '* Created by program: fequtl'
      write(stdtab,'(a1,a)') '*', version_run_date_time_string
      call svn_report(stdout, stdtab)   !Report the repository loc, etc.

      

C     SET VARIOUS VALUES
 
      EFLAG = 0
 
C     INTIALIZE SELECTED VALUES IN COMMON
 
      CALL INIT
 
C     BE SURE TO SET NPNTD TO 0 IN BLOCK DATA SUBROUTINE
 
C     INITIALIZE TABLE DIRECTORY
 
      DO 100 J=1,PMXTAB
        TABDIR(J) = 0
        FTPNT(J) = 0
 100  CONTINUE
 
      if(conf_flag == 1) then
        INQUIRE(FILE=conf_file, EXIST=THERE)
        N = LEN_TRIM(conf_file)
        IF(THERE) THEN
         
           OPEN(stdsys, FILE = conf_file, STATUS = 'OLD')
           write(*,*) ' FEQUTL configuration from conf file.'
           write(stdout,83) conf_file(1:n)
        ELSE
          WRITE(*,*) ' '
          WRITE(*,*) ' File named: ',conf_file(1:n),' not found.'
          WRITE(*,*) ' Please check spelling of configuration file.'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
      else
        write(*,*) ' FEQUTL configuration from master-input file.'
        write(stdout,84) 
        
      endif

C     INPUT THE CONTROL INFORMATION FROM THE SYSTEM DATASET.
c     Info may be at the head of the master-input file, old 
c     method, or in its own file, new method.  In the former
c     case, stdsys is the same as stdin.  Otherwise, stdsys
c     is unique.
 
      CALL INSYS
     I          (STDSYS, MAXCMD, NARG, STDOUT,
     O           UNITS, NFAC, GRAV, MAXKNT, EPS,
     O           NCMD, CMDTAB, CMDVAL)
 
      GRAV2 = GRAV + GRAV
      SQRT_GRAV = SQRT(GRAV)
 
C     SET THE VALUE FOR THE MAXIMUM WEIR COEF FOR A BROAD CRESTED WEIR
C     WITH CRITICAL DEPTH AT THE CREST.  USED IN EMBANKQ FOR COMPUTING
C     ENERGY LOSSES WHEN FLOW ON THE WEIR CREST IS CLOSE TO CRITICAL.
 
      BCWMAX = SQRT(GRAV)/(1.5*SQRT(1.5))
 
C     Set the special small difference value in epscom
      IF(GRAV.GT.15.0) THEN
        EPSDIF = 0.01
      ELSE
        EPSDIF = 0.003048
      ENDIF

      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      READ(LINE, 49,ERR=991) DZLIM
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      READ(LINE, 40,ERR=991) NRZERO
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      READ(LINE, 42,ERR=991) REPLY
 
 
      WRITE(STDOUT,*) '   '
      WRITE(STDOUT,70) DZLIM, NRZERO, REPLY
      IF(NRZERO.LE.0.0.OR.DZLIM.LE.0.0) THEN
        WRITE(STDOUT,*) '*ERR:502* INVALID VALUES FOR NRZERO OR',
     A                  ' DZLIM'
        WRITE(STDOUT,*) '   BOTH MUST BE > 0.0'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      IF(REPLY.EQ.'YES') THEN
        IUSGS = 1
      ELSE IF(REPLY.EQ.'NO') THEN
        IUSGS = 0
      ELSE
        WRITE(STDOUT,74)
        IUSGS = 0
      ENDIF
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      READ(LINE,45,ERR=991) NXTNAM
      IF(NXTNAM(1:6).NE.'EPSARG') THEN
        WRITE(STDOUT,*) ' CHECK INPUT STREAM. THIS VERSION OF FEQUTL'
        WRITE(STDOUT,*) ' HAS ADDED EPSARG AND EPSF AFTER USGSBETA'
        WRITE(STDOUT,*) ' TRY EPSARG=5.E-4 AND EPSF=1.E-3'
        STOP 'Abnormal stop. Errors found.'
      ELSE
        READ(NXTNAM,44) EPSARG
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
        READ(LINE,46,ERR=991) EPSF, EPSABS
        IF(EPSABS.EQ.0.0) THEN
          WRITE(STDOUT,91) 
          IF(GRAV.GT.15.0) THEN
            EPSABS = EPSF
          ELSE
            EPSABS = EPSF*0.3048
          ENDIF
        ENDIF
      ENDIF
 
      WRITE(STDOUT,76) EPSARG, EPSF, EPSABS
 
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      READ(LINE,45,ERR=991) NXTNAM
      IF(NXTNAM(1:6).NE.'EXTEND') THEN
        WRITE(STDOUT,81)
        GXTEND = 1
        WRITE(STDOUT,77)
      ELSE
        READ(NXTNAM,41) REPLY
        IF(REPLY.EQ.'YES') THEN
          GXTEND = 1
          WRITE(STDOUT,78)
          WRITE(STDOUT,77)
        ELSE IF(REPLY.EQ.'NO') THEN
          GXTEND = 0
          WRITE(STDOUT,79)
        ELSE
          WRITE(STDOUT,80)
          GXTEND = 0
          WRITE(STDOUT,79)
        ENDIF
      ENDIF

C     Seek the minimum flow value.  May not be present.
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      IF(LINE(1:4).NE.'MINQ') THEN
C       Put the line back!
        BACKSPACE(stdsys)
C       Set the default values
        IF(GRAV.GT.15.0) THEN
          MINQ = 0.2
        ELSE
          MINQ = 0.2/35.315
        ENDIF
      ELSE
         READ(LINE(6:15),'(F10.0)') MINQ
      ENDIF

      WRITE(STDOUT,92) MINQ
      WRITE(STDOUT,*) ' '
      
c     seek the control value for processing tables for cubic interpolation. 
C     May not be present.
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      IF(LINE(1:5).NE.'TWOD_') THEN
C       Put the line back!
        BACKSPACE(stdsys)
C       Set the default value
        twod_cubic_out = 'NO'
      ELSE
        is = index(line,'NO')
        if(is.gt.0) then
          twod_cubic_out = 'NO'
        else
          is = index(line,'YES')
          if(is.gt.0) then
            twod_cubic_out = 'YES'
          else
            write(stdout,*) ' Invalid response for TWOD_CUBIC_OUT'
            write(stdout,*) ' Must be : YES or NO.'
            write(stdout,*) ' Setting value to NO.'
            twod_cubic_out = 'NO'
          endif
        endif
      ENDIF

c     seek the global home name. May not be present.
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      IF(LINE(1:5).NE.'GHOME') THEN
C       Put the line back!
        BACKSPACE(stdsys)
C       Default value already set above
      else
c       Skip over the =, read, and then shift left to 
c       strip leading spaces if any. 
        read(line(7:),'(a)') ghome
        ghome = adjustl(ghome)
        write(stdout,82) ghome(1:len_trim(ghome))
        write(stdtab,88) ghome(1:len_trim(ghome))
        call getsvn_rev(stdout, ghome,
     o                      svn_rev, svn_mod)
        if(svn_rev >0 ) then
          write(stdtab,86) svn_rev
          if(svn_mod > 0) then
            write(stdtab,87)
          endif
        endif
      
      endif

c     Seek the global grid/datum/unitsys information. Consists of five lines
c     of input.  If the first line is present, then all five lines must be 
c     present!  Makes life a bit simpler:)
      CALL inline
     I          (stdsys, STDOUT,
     O           LINE)
      IF(LINE(1:6).NE.'G_ZONE') THEN
C       Put the line back!
        BACKSPACE(stdsys)
C       Set the default values
        G_ZONE = 'NONE'
        G_HGRID = 'NONE'
        G_VDATUM = 'NONE'
        G_UNITSYS = 'NONE'
        G_BASIS = 'NONE'
      else
        call get_single_named_item(stdout, line,
     o                             item_name, g_zone, eflag)
        CALL inline
     I            (stdsys, STDOUT,
     O             LINE)
        call get_single_named_item(stdout, line,
     o                             item_name, g_hgrid, eflag)
        CALL inline
     I            (stdsys, STDOUT,
     O             LINE)
        call get_single_named_item(stdout, line,
     o                             item_name, g_vdatum, eflag)
        CALL inline
     I            (stdsys, STDOUT,
     O             LINE)
        call get_single_named_item(stdout, line,
     o                             item_name, g_unitsys, eflag)
        CALL inline
     I            (stdsys, STDOUT,
     O             LINE)
        call get_single_named_item(stdout, line,
     o                             item_name, g_basis, eflag)
        if(eflag == 0) then 
          write(stdout,56) g_zone, g_hgrid, g_vdatum, g_unitsys, g_basis
        else
          STOP 'Abnormal stop. Errors found.'
        endif
      endif


      if(conf_flag == 1) then 
        call skip_header_if_present(
     I                              STDIN, STDOUT, NCMD, CMDTAB)
      endif

c     We have processed all the header values and are about to begin processing commands. 
c     At this point we analyze the master-input file to compute various md5 digests.
c     This will be made optional later by user input in the standard header.  

      if(trk_files) then
        call find_md5_for_fequtl_input(stdscr, stdout, stdsys, stdin, 
     i                               fname, NCMD, CMDTAB, cmdval,
     i                               conf_file, ghome)
      endif



C     READ NEXT COMMAND AND DETERMINE ACTION
 
      CFLAG=0
 500  CONTINUE
 
C       RESET THE EXTEND OPTION TO RECOVER THE GLOBAL VALUE
 
        EXTEND = GXTEND
 
        CALL NXTCMD
     I             (STDIN, STDOUT, NCMD, CMDTAB, CMDVAL,
     O              NEXT)
 
        IF(NEXT.GT.0) GO TO 110
          EFLAG = 1
          IF(CFLAG.EQ.1) GOTO 105
            WRITE(STDOUT,51)
            CFLAG = 1
 105      CONTINUE
        GOTO 500
 
 110    CONTINUE
        IF(CFLAG.EQ.0) GOTO 120
          WRITE(STDOUT,52)
          CFLAG = 0
 120    CONTINUE
 
      GFLAG = GFLAG + EFLAG
      EFLAG = 0
 
      GOTO(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
     A      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
     B      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 240, 241,
     C      242,243, 244), NEXT
 
      WRITE(STDOUT,*) ' *BUG:501* INVALID ADDRESS FOR A COMMAND'
      WRITE(STDOUT,*) ' NEXT=',NEXT 
      STOP 'Abnormal stop. Errors found.'
 
 1    CONTINUE
        MODE = 1
        SNFLGU = 0
 
        CALL FEQX
     I           (STDIN, STDOUT, STDTAB, NFAC, MODE,
     M            TABDIR, FTP,
     O            EFLAG)
        GOTO 1000
 
 2    CONTINUE
 
        CALL FLDIN
     I            (STDIN, STDOUT, STDFLD,
     O             EFLAG, FHEAD)
        GOTO 1000
 
 3    CONTINUE
        CALL SPBRID
     I             (STDIN, STDOUT, STDTAB, NFAC,
     M              TABDIR, EFLAG)
        GOTO 1000
 
 4    CONTINUE
C       FIND TWO-D TABLE FOR CULVERT FLOW WITH POSSIBLE FLOW OVER
C       THE ROADWAY
 
        CALL CULVRT
     I             (STDIN, STDOUT, STDTAB, NFAC,
     M              TABDIR, EFLAG, FTP, FTKNT)
 
C        WRITE(STDOUT,*) ' EFLAG=',EFLAG
        GOTO 1000
 
 5    CONTINUE
C       FINISH COMMAND
        IF(FLOOD.EQ.1) THEN
C         IF FLOODWAY OPTION IS ENABLED OUTPUT SUMMARY TABLE
 
          WRITE(STDOUT,*) ' '
          WRITE(STDOUT,*) ' Summary of Floodway Computation results'

          WRITE(STDOUT,94)
 
          DO 5000 ITABA=1,PMXTAB
 
            IF(FLDOPT(ITABA).NE.'    ') THEN
              IF(FLDRT(ITABA).EQ.-1.E30) THEN
                RIGHT = '  -inf  '
              ELSE
                WRITE(RIGHT,'(F8.1)') FLDRT(ITABA)
              ENDIF
              IF(FLDLT(ITABA).EQ.1.E30) THEN
                LEFT = '  +inf  '
              ELSE
                WRITE(LEFT,'(F8.1)') FLDLT(ITABA)
              ENDIF
              IF(FEQBOT(ITABA).EQ.-1.E30) THEN
                BOT = '  same  '
              ELSE
                WRITE(BOT,'(F8.2)') FEQBOT(ITABA)
              ENDIF
              TABID = GET_TABID(ITABA)
              GISID = GET_GISID(FTPNT(ITABA))
              if(fldflow(itaba).gt.0.0) then
                velocity = fldflow(itaba)/fldarea(itaba)
              else
                velocity = 0.0
              endif
              WRITE(STDOUT, 93 )
     A          TABID, GISID, FLDOPT(ITABA), FLDELV(ITABA), BOT, LEFT,
     B          RIGHT, FLDLOS(ITABA), fldarea(itaba), velocity
            ENDIF
 5000     CONTINUE
 
        ENDIF
        WRITE(*,*) ' FINISH found.'
        IF(GFLAG.GT.0) THEN
          WRITE(STDOUT,54)
          WRITE(*,*) ' Errors reported.'
        ELSE
          WRITE(*,*) ' No errors reported.'
        ENDIF
 
 

 
        STOP 
 
 6    CONTINUE
C        CALL SAME(STDIN,STDOUT,EFLAG)
        GOTO 1000
 
 7    CONTINUE
C       PROCESS CROSS SECTION FROM DIGITIZER
C       DISABLE SINUOSITY ELEMENTS
        SNFLGU = 0
        MODE = 3
        CALL FEQX
     I           (STDIN, STDOUT, STDTAB, NFAC, MODE,
     M            TABDIR, FTP,
     O            EFLAG)
        GOTO 1000
 
 8    CONTINUE
C       PROCESS CROSS SECTION IN LIST FORMAT
 
C       DISABLE SINUOSITY ELEMENTS
        SNFLGU = 0
        MODE = 2
        CALL FEQX
     I           (STDIN, STDOUT, STDTAB, NFAC, MODE,
     M            TABDIR, FTP,
     O            EFLAG)
        GOTO 1000
 9    CONTINUE
C       COMPUTE FLOW OVER THE ROADWAY FROM A CROSS-SECTION TABLE
C       ASSUMING CRITICAL DEPTH AND THAT VELOCITY HEAD OF APPROACH
C       WILL BE USED FOR THE HEAD ON THE RESULTING TABLE.
C       THE USER SUPPLIES INFORMATION ON THE TABLE INTERVAL AS WELL
C       AS A FACTOR TO USE TO REDUCE CRITICAL FLOW.
 
C        CALL FNDRF(STDIN, STDOUT, STDTAB, NFAC, TABDIR,
C     A                 EFLAG)
        GOTO 1000
 
 10   CONTINUE
C       FIND CROSS SECTION TABLE FOR A SEWER PIPE.
        CALL SEWER
     I            (STDIN, STDOUT, STDTAB, NFAC,
     M             TABDIR, EFLAG, FTP)
        GOTO 1000
 
 11   CONTINUE
C       FIND CROSS SECTION TABLE FOR MULTIPLE SEWER PIPES
        CALL PIPES
     I            (STDIN, STDOUT, STDTAB, NFAC,
     M             TABDIR, EFLAG, FTP)
        GOTO 1000
 
 12   CONTINUE
C       INPUT FUNCTION TABLES FOR ACCESS BY SUBSEQUENT COMMANDS.
        WRITE(*,90)
 
        HSLOT = 0.0
        IS = STDIN
 1200   CONTINUE
          make_tab_index = 'NO  '  !Disable making a table index during input of function tables
          current_file = ' '       !Argument needed even if no table index is made
          index_knt = 0
          CALL FTABIN
     I               (IS, STDIN, STDOUT, OUTPUT, FTP, PMXTAB, MRFTAB,
     I                HSLOT, current_file,
     M                EFLAG, FTKNT, FTPNT,
     O                TABLE, NXTNAM)
 
          IF(TABLE.EQ.-1) GOTO 300
          IF(IS.NE.STDIN) THEN
            CALL FREE_UNIT(STDOUT, IS)
          ENDIF
 
          IF(TABLE.EQ.-15) THEN
            IS = GET_UNIT(STDOUT)
            INQUIRE(FILE=NXTNAM, EXIST=THERE)
            IF(THERE) THEN
              OPEN(UNIT=IS, FILE=NXTNAM, STATUS='OLD')
            ELSE
              WRITE(STDOUT,*) ' FILE NAMED:',NXTNAM,' NOT FOUND.'
              WRITE(STDOUT,*) ' CHECK FUNCTION TABLE FILE NAME.'
              STOP 'Abnormal stop. Errors found.'
            ENDIF
          ELSEIF(TABLE.EQ.-16) THEN
C           Set the local value of the home directory
            CALL SET_HOME(NXTNAM)
c           Check for version control and report
            
            call getsvn_rev(stdout, nxtnam,
     o                      svn_rev, svn_mod)
            write(stdtab,89) nxtnam
            if(svn_rev >0 ) then
              write(stdtab,86) svn_rev
              if(svn_mod > 0) then
                write(stdtab,87)
              endif
            endif

            
          ELSE
            IS = STDIN
          ENDIF
 
          GOTO 1200
 300    CONTINUE
        IF(IS.NE.STDIN) THEN
          CALL FREE_UNIT(STDOUT,IS)
        ENDIF
 
        WRITE(STDOUT,*) ' '
        WRITE(STDOUT,'('' There are now'',I7,
     A   '' function tables stored'')') FTKNT
        WRITE(STDOUT,'(3X, I9,'' locations out of '',I9,'' are used'')')
     A   FTP, MRFTAB
        GOTO 1000
 
 13   CONTINUE
C       COMPUTE FLOW OVER EMBANKMENT SHAPED WEIRS USING USGS
C       PROCEDURE
        CALL EMBANK
     I             (STDIN, STDOUT, STDTAB, MINQ,
     M              EFLAG, FTP)
        GOTO 1000
 
 14   CONTINUE
C       COMPUTE HYDRAULIC JUMP TABLE FOR FORCED JUMP LOCATION
C        CALL JUMP(STDIN, STDOUT, STDTAB)
        GOTO 1000
 
 15   CONTINUE
C       COMPUTE CRITICAL FLOW TABLE FOR A CONSTRICTION ASSUMING
C       NO DOWNSTREAM EFFECT
        CALL CRITQ
     I            (GRAV, STDIN, STDOUT, STDTAB,
     M             EFLAG, FTP)
        GOTO 1000
 
 16   CONTINUE
C       FIND THE PEAK FLOW FROM THE GENERALIZED RITTER SOLUTION
 
        CALL RITTER
     I             (GRAV, STDIN, STDOUT,
     M              EFLAG)
        GOTO 1000
 
 17   CONTINUE
C       FIND THE FRACTION OF MAXIMUM CAPACITY AS A FUNCTION OF TIME
C
C        CALL PFIND(GRAV, STDIN, STDOUT, STDTAB, EFLAG)
        GOTO 1000
 
 18   CONTINUE
C       FIND CROSS SECTION TABLE FOR MULTIPLE CONDUITS
 

        CALL MULCON
     I             (STDIN, STDOUT, STDTAB, NFAC,
     M              TABDIR, EFLAG, FTP)

        GOTO 1000
 
 19   CONTINUE
C       FIND TWO-D TABLE FOR FLOW THROUGH A PRISMATIC CHANNEL
C       AT SUBCRTICAL SLOPE.
 
        CALL CHNTAB
     I             (STDIN, STDOUT, STDTAB, GRAV, MINQ,
     M              TABDIR, EFLAG, FTP)
 
C        WRITE(STDOUT,*) ' EFLAG=',EFLAG
        GOTO 1000
 
 20   CONTINUE
C       FIND TWO-D TABLE FOR FLOW THROUGH A TRANSITION.
        CALL EXPCON
     I             (STDIN, STDOUT, STDTAB, GRAV,
     M              TABDIR,
     O              EFLAG)
        GOTO 1000
 
 21   CONTINUE
C       COMPUTE CROSS SECTION TABLES OR PRODUCE FEQX INPUT FROM
C       HEC2 INPUT FILES.
 
        CALL HEC2X
     I            (STDIN, STDOUT, STDTAB, NFAC,
     M             TABDIR, FTP,
     O             EFLAG)
        GOTO 1000
 
 22   CONTINUE
C       COMPUTE A LIMIT TO THE CRITICAL FLOW IN A CLOSED CONDUIT
C       CROSS SECTION
        CALL QCLIM
     I            (STDIN, STDOUT,
     M             EFLAG)
        GOTO 1000
 
 23   CONTINUE
C       INTERPOLATE AND OUTPUT CROSS SECTIONS
 
        CALL XSTMAK
     I             (STDIN, STDOUT, STDTAB,
     M              EFLAG, FTP, FTKNT)
        GOTO 1000
 
 24   CONTINUE
 
C        CALL SDWEIR(GRAV, STDIN, STDOUT, STDTAB, EFLAG, FTP, FTKNT)
          WRITE(STDOUT,*) ' SIDEWEIR NOT YET COMPLETE'
        GOTO 1000
 
 25   CONTINUE
        MODE = 0
C       DISABLE SINUOSITY ELEMENTS
        SNFLGU = 0
        CALL FQXE
     I           (STDIN, STDOUT, STDTAB, NFAC,
     M            TABDIR, FTP,
     O            EFLAG)
        GOTO 1000
 
 26   CONTINUE
        CALL CHANEL
     I             (STDIN, STDOUT, STDTAB, NFAC,
     M              TABDIR, FTP,
     O              EFLAG, MODE)
        GOTO 1000
 
 27   CONTINUE
        CALL WPROX
     I            (STDIN, STDOUT, STDTAB, NFAC,
     M             TABDIR, FTP,
     O             EFLAG)
        GOTO 1000
 
 28   CONTINUE
        CALL WPROQZ
     I             (STDIN, STDOUT,
     O              EFLAG)
        GOTO 1000
 
 29   CONTINUE
        CALL WPRO14
     I             (STDIN, STDOUT, STDTAB,
     M              TABDIR, EFLAG)
        GOTO 1000
 
 30   CONTINUE
        CALL UFGATE
     I             (GRAV, STDIN, STDOUT, STDTAB,
     M              EFLAG, TABDIR)
        GOTO 1000
31    CONTINUE
        CALL RISERCLV(GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR)

        GOTO 1000 
32    CONTINUE
        CALL   ORIFICE(GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR, FTP)

      GOTO 1000

33    CONTINUE
      CALL PUMPITEMS
     I              ('SFWMD   ', GRAV, STDIN, STDOUT, STDTAB,
     M               EFLAG, FTP)
      GOTO 1000

34    CONTINUE
      CALL PUMPITEMS
     I              ('PUMPLOSS', GRAV, STDIN, STDOUT, STDTAB,
     M               EFLAG, FTP)
      GOTO 1000

35    CONTINUE
        CALL SET_SLOT(STDIN, STDOUT, STDTAB)
      GOTO 1000

36    CONTINUE
        CALL CLEAR_SLOT(STDTAB)
      GOTO 1000

37    CONTINUE
        CALL   INV_GATE(GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR)
      GOTO 1000

38    CONTINUE
        CALL   UFGCULV
     I                (GRAV, STDIN, STDOUT, STDTAB, FTP,
     M                 EFLAG, TABDIR)
      GOTO 1000

39    CONTINUE
        CALL MKEMBANK(STDIN, STDOUT, EFLAG)
      GOTO 1000


240   CONTINUE
        CALL  MKWSPRO(STDIN, STDOUT, TABDIR, 
     O                    EFLAG)
      GOTO 1000

241   CONTINUE
        CALL WPRO14_NEW
     I             (STDIN, STDOUT, STDTAB,
     M              TABDIR, EFLAG)
        GOTO 1000

242   CONTINUE
        CALL LPRFIT
     I             (STDIN, STDOUT, STDTAB, FTP,
     M              EFLAG)
        GOTO 1000


243   CONTINUE
        CALL  MAKELAKEXS(STDIN, STDOUT, EFLAG)

        GOTO 1000

244   CONTINUE
        CALL SET_SLOTE(STDIN, STDOUT, STDTAB)

        GOTO 1000

 1000 CONTINUE


      WRITE(STDOUT,*) ' '
      GFLAG = GFLAG + EFLAG
      EFLAG = 0
      GOTO 500
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
 
      END

C
C
C
      subroutine skip_header_if_present(
     I                                  STDIN, STDOUT, NCMD, CMDTAB)
 
C     + + + PURPOSE + + +
C     Read files from stdin, skipping any header block lines
c     and returning with stdin set so that it is at the first 
c     command.  Depends on the first command, at least, not having
c     anything following it on its line.   That has been my practice
c     but it is not strictly required.  ddf 8 June 2004.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NCMD, STDIN, STDOUT
      CHARACTER CMDTAB(NCMD)*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     NCMD   - Number of commands
C     CMDTAB - Vector holding the command names
 
C     + + + LOCAL VARIABLES + + +
      integer i, match
      CHARACTER  LINE*80
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BINSER, inline
 
C***********************************************************************
c     Read first 80 chars of each line, see if the line will match 
c     any of the commands.  If not, go get the next line, otherwise,
c     backspace stdin, and exit. 

100   continue

        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)

        if(line == 'ENDFILE') then
          write(*,*) ' No commands found in master-input file.'
          STOP 'Abnormal stop. Errors found.'
        endif
        match = 0
        do i=1,ncmd
          if(line == cmdtab(i)) then
            match = 1
            exit
          endif
        enddo 

        if(match == 0) then

          goto 100

        else
          backspace(stdin)
          return
        endif

        write(*,*) 'bug in handling configuration file.'
        STOP 'Abnormal stop. Bug found.'
      end


C  ***********************************************************************
C  *  Warning:  This program is large and complex and  extensive         *  
C  *  knowledge of its design, purpose, and limitations is required      *  
C  *  in order to apply it properly.  Application of this program by an  *
C  *  unqualified user for any other purpose than an educational one is  *
C  *  not only unwise but is also unethical.  The user of this           *
C  *  program is totally responsible for its use and application and for *
C  *  any actions or events which follow therefrom.  Any user of this    *
C  *  program  holds the developer of the program harmless from          *
C  *  damages of any kind.                                               *
C  *                                                                     *       
C  *  The developer has used reasonable care in the construction and     *
C  *  testing of the program.  However, in a program of this size and    *
C  *  complexity, it is impossible to verify more than a minute number of*
C  *  possible options or applications.  The developer is continuing to  *
C  *  modify and use the program and is interested in information on     *
C  *  operational problems encountered in its application.  However, the *
C  *  developer gives no assurance that the problem can or will be       *
C  *  rectified.                                                         *
C  *                                                                     *       
C  *  This program is not to be sold in any form modified or otherwise.  *
C  ***********************************************************************

