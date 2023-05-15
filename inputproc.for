C     Routines used for processing input lines from the user.

C
C
C
      SUBROUTINE CLEAR_PUTGET()

C     Clear all put/get flags for storing and accessing function
C     tables stored in FTAB/ITAB.  These flags should be cleared
C     before the processing of any command that used any of them.

      IMPLICIT NONE

      INCLUDE 'putget.cmn'

C     Local 
      INTEGER I      
C***********************************************************************
      DO 100 I=1,MXOVER
        PGOVER(I) = 0
100   CONTINUE
      RETURN
      END
      
C           
C
C
C
      SUBROUTINE   GET_MULTIPLE_REAL_VALUES(
     I      STDIN, STDOUT, INTVAL, REAVAL, CONTINUATION_VALUE, 
     I      OPT, NVAL, 
     M      LINE, ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN,
     M      EFLAG2, TERM, TERML, TERMCLS, ITEM_KNT,
     O      NUMBER_OF_VALUES, THE_VALUES)

C     Get multiple real values following a variable name.


      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER CONTINUATION_VALUE, EFLAG2, INTVAL, NUMBER_OF_VALUES, 
     A        NVAL, OPT, REAVAL, STDIN, STDOUT, ITEM_KNT
      INTEGER CLEN(NVAL), IVAL(NVAL), ITEM_TYPE(NVAL), TERML(NVAL),
     A        TERMCLS(NVAL)
      REAL RVAL(NVAL), THE_VALUES(NVAL)
      REAL*8 DPVAL(NVAL)
      CHARACTER CVAL(NVAL)*(*), LINE*(*), TERM(NVAL)*1


C     Definition of arguments.

C     STDIN- unit number for input 
C     STDOUT- unit number for output of errors and messages.
C     INTVAL- code for an integer value
C     REAVAL- code for a real value
C     CONTINUATION_VALUE- code for continuing to the next line
C     LINE- the character string to scan.
C     NVAL- the maximum number of values to expect.  It is an error if
C           this number is exceeded.
C     OPT-  if 0 then
C       ITEM_TYPE- vector giving the type of value expected for each value in
C              LINE
C       IVAL, RVAL, DPVAL, CVAL- vectors of the correct type for storing
C           the values found in LINE.  
C       CLEN- vector giving the length of the character string if the
C             value is of type CHARACTER.
C       EFLAG- set to 1 if an error was found in LINE.
C     else if 1
C       ITEM_TYPE- vector giving the type of value found in LINE
C       CVAL- gives the value returned
C       CLEN- give the length of the value

C     else if 2
C       like 0 but ITEM_TYPE gives the type as found in LINE and
C       CVAL and CLEN also contain the string form of argument.

C     In either case:

C        TERM- gives the single character value for the terminator
C        TERML- gives the terminator length in case we want to 
C               make it more than one character.
C        TERMCLS- gives the class number of the terminator.  
C        ITEM_KNT- the number of values actually found. A value may be null or
C             defaulted and it is still counted.
C     NUMBER_OF_VALUES- number of values returned
C     THE_VALUES- vector containing the values found


C     External program units

      INTEGER LENSTR

      EXTERNAL LENSTR, GETVAL, inline, CHK_FOR_GIVEN_TYPE
      
C     Local variables

      INTEGER IE, IS, IT, KNT, J
C     ********************************FORMATS***************************
50    FORMAT(/,' *ERR:728* There were ',I3,' invalid real numbers ',
     A         'found in the current line.')
56    FORMAT(/,' Unable to continue due to previous errors.')
58    FORMAT(/,' *ERR:729* Found ',I3,' values on line > limit of ',I3)
60    FORMAT(/,' Processing:',A)
62    FORMAT(/,' Seeking additional values.')
C***********************************************************************

      IF(ITEM_TYPE(2).EQ.INTVAL) THEN
C       The user has given a count of the subsections.  Get the values
C       of roughness on the current line. NSUB gives the total
C       number of values expected.

        NUMBER_OF_VALUES = IVAL(2)
        IS = 3
        IE = ITEM_KNT
        KNT = 0

128     CONTINUE
          DO 130 J=IS, IE
            KNT = KNT + 1
            IF(KNT.GT.NVAL) THEN
              WRITE(STDOUT,58) KNT, NVAL
              STOP 'Abnormal stop.  Errors found.' 
            ENDIF
            THE_VALUES(KNT) = RVAL(J)
130       CONTINUE
          IF(KNT.EQ.NUMBER_OF_VALUES) THEN
C           All the claimed  values have been found. 
C           We are done.
            GOTO 200
          ELSE
C           Get the next line and parse it.  It should only contain
C           numeric data. 
            WRITE(STDOUT,62)
            CALL inline
     I                 (STDIN, STDOUT,
     O                  LINE)
            IT = LENSTR(LINE)
            WRITE(STDOUT,60) LINE(1:IT)
            LINE(IT+1:IT+1) = ''''
            CALL GETVAL
     I                 (STDOUT, LINE, NVAL, OPT,
     O                  ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, 
     O                  CLEN, EFLAG2, TERM, TERML, TERMCLS,
     O                  ITEM_KNT)
            IF(EFLAG2.NE.0) THEN
               WRITE(STDOUT,56)
               STOP 'Abnormal stop.  Errors found.' 
            ENDIF
            EFLAG2 = 0
            CALL CHK_FOR_GIVEN_TYPE(
     I                      1, ITEM_KNT, ITEM_TYPE, REAVAL,
     O                      EFLAG2)
            IF(EFLAG2.NE.0) THEN
              WRITE(STDOUT,50) EFLAG2
              STOP 'Abnormal stop.  Errors found.' 
            ENDIF
132         CONTINUE
            IS = 1
            IE = ITEM_KNT
            GOTO 128
        ENDIF           
      ELSE
        KNT = 0
        IS = 2 
133     CONTINUE
C         The user has not given a count of the items
          IF(ITEM_TYPE(ITEM_KNT).EQ.CONTINUATION_VALUE) THEN
            IE = ITEM_KNT - 1
          ELSE
            IE = ITEM_KNT
          ENDIF
          EFLAG2 = 0
          CALL CHK_FOR_GIVEN_TYPE(IS, IE, ITEM_TYPE, REAVAL,
     O                    EFLAG2)
          IF(EFLAG2.NE.0) THEN
            WRITE(STDOUT,50) EFLAG2
            STOP 'Abnormal stop.  Errors found.' 
          ENDIF
          DO 135 J=IS,IE
            IF(KNT.GT.NVAL) THEN
              WRITE(STDOUT,58) KNT, NVAL
              STOP 'Abnormal stop.  Errors found.' 
            ENDIF
            KNT = KNT + 1
            THE_VALUES(KNT) = RVAL(J)
135       CONTINUE
          IF(ITEM_TYPE(ITEM_KNT).EQ.CONTINUATION_VALUE) THEN
C           Get the next line and process it.
            CALL inline
     I                 (STDIN, STDOUT,
     O                  LINE)
            IT = LENSTR(LINE)
            WRITE(STDOUT,60) LINE(1:IT)
            LINE(IT+1:IT+1) = ''''
            CALL GETVAL
     I                 (STDOUT, LINE, NVAL, OPT,
     O                  ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, 
     O                  CLEN, EFLAG2, TERM, TERML, TERMCLS,
     O                  ITEM_KNT)
            IF(EFLAG2.NE.0) THEN
               WRITE(STDOUT,56) LINE
               STOP 'Abnormal stop.  Errors found.' 
            ENDIF
            IS = 1
            GOTO 133
          ELSE
C            All lines processed
             NUMBER_OF_VALUES = KNT
             GOTO 200     
          ENDIF
      ENDIF

      WRITE(STDOUT,64) 
64    FORMAT(/,' *BUG:XXX* GET_MULTIPLE_REAL_VALUES:',
     A         ' Should not be here.')
      STOP 'Abnormal stop: bug found.'
200   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   GET_TABID_VALUES(
     I      STDIN, STDOUT, INTVAL, REAVAL, CHRVAL, CONTINUATION_VALUE, 
     I      OPT, NVAL, 
     M      ISTART, LINE, ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN,
     M      EFLAG2, TERM, TERML, TERMCLS, ITEM_KNT,
     O      NUMBER_OF_VALUES, THE_VALUES)

C     Get multiple table ids following a variable name.  The line
C     on which the name appears has already been parsed with GETVAL.
C     Continuation lines, if any, are processed here as well.  The
C     variable name MUST be followed by an = to distinguish between
C     the variable and the responses.

      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER CONTINUATION_VALUE, EFLAG2, INTVAL, NUMBER_OF_VALUES, 
     A        NVAL, OPT, REAVAL, CHRVAL, STDIN, STDOUT, ITEM_KNT,
     B        ISTART
      INTEGER CLEN(NVAL), IVAL(NVAL), ITEM_TYPE(NVAL), TERML(NVAL),
     A        TERMCLS(NVAL), THE_VALUES(NVAL)
      REAL RVAL(NVAL)
      REAL*8 DPVAL(NVAL)
      CHARACTER CVAL(NVAL)*(*), LINE*(*), TERM(NVAL)*1


C     Definition of arguments.

C     STDIN- unit number for input 
C     STDOUT- unit number for output of errors and messages.
C     INTVAL- code for an integer value
C     REAVAL- code for a real value
C     CHRVAL- code for a character values(identifier) not a string
C     CONTINUATION_VALUE- code for continuing to the next line
C     LINE- the character string to scan.
C     NVAL- the maximum number of values to expect.  It is an error if
C           this number is exceeded.
C     ISTART- On entry:index into ITEM_TYPE, and related vectors at which to 
C             start processing
C           - On exit: index value  of the terminating identifier if any;
C              0 if no terminating identifier was found.
C     OPT-  if 0 then
C       ITEM_TYPE- vector giving the type of value expected for each value in
C              LINE
C       IVAL, RVAL, DPVAL, CVAL- vectors of the correct type for storing
C           the values found in LINE.  
C       CLEN- vector giving the length of the character string if the
C             value is of type CHARACTER.
C       EFLAG- set to 1 if an error was found in LINE.
C     else if 1
C       ITEM_TYPE- vector giving the type of value found in LINE
C       CVAL- gives the value returned
C       CLEN- give the length of the value

C     else if 2
C       like 0 but ITEM_TYPE gives the type as found in LINE and
C       CVAL and CLEN also contain the string form of argument.

C     In any case:

C        TERM- gives the single character value for the terminator
C        TERML- gives the terminator length in case we want to 
C               make it more than one character.
C        TERMCLS- gives the class number of the terminator.  
C        ITEM_KNT- the number of values actually found. A value may be null or
C             defaulted and it is still counted.
C     NUMBER_OF_VALUES- number of values returned
C     THE_VALUES- vector containing the values found


C     External program units

      INTEGER LENSTR

      EXTERNAL LENSTR, GETVAL, inline, GET_INTERNAL_TAB_NUMBER
      
C     Local variables

      INTEGER EQUAL
      PARAMETER (EQUAL=8)
      INTEGER IS, IT, KNT, ITEMP
      CHARACTER TABID*16
C     ********************************FORMATS***************************
50    FORMAT(/,' *ERR:730* Expected an integer or id but found: ',A,
     A       ' in the line being processed.')
56    FORMAT(/,' Errors found but attempting to continue.')
58    FORMAT(/,' *ERR:729* Found ',I3,' values on line > limit of ',I3)
60    FORMAT(/,' Processing:',A)
62    FORMAT(/,' Seeking additional values.')
C***********************************************************************
      KNT = 0
      IS = ISTART
133   CONTINUE
        
        IF(IS.GT.ITEM_KNT) THEN
C         End of line, no continuation, therefore done!
          ISTART = 0
          NUMBER_OF_VALUES = KNT
          GOTO 200
        ENDIF

C       Classify the next item in the results from GETVAL
        IF(TERMCLS(IS).EQ. EQUAL) THEN
C         End of processing for the current identifier
          ISTART = IS
          NUMBER_OF_VALUES = KNT
          GOTO 200
        ELSE
          IF(ITEM_TYPE(IS).EQ.CONTINUATION_VALUE) THEN
C           Get the next line, parse it with GETVAL and process
            WRITE(STDOUT,62) 
            CALL inline
     I                 (STDIN, STDOUT,
     O                  LINE)
            IT = LENSTR(LINE)
            WRITE(STDOUT,60) LINE(1:IT)
            LINE(IT+1:IT+1) = ''''
            CALL GETVAL
     I                 (STDOUT, LINE, NVAL, OPT,
     O                  ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, 
     O                  CLEN, EFLAG2, TERM, TERML, TERMCLS,
     O                  ITEM_KNT)
            IF(EFLAG2.NE.0) THEN
               WRITE(STDOUT,56) LINE
            ENDIF
            IS = 1
            GOTO 133
          ELSE
C           We have a value.  Count it and save it.
            KNT = KNT + 1
            IF(KNT.GT.NVAL) THEN
              WRITE(STDOUT,58) KNT, NVAL
              KNT = KNT - 1
            ENDIF
            IF(ITEM_TYPE(IS).EQ.CHRVAL.OR.ITEM_TYPE(IS).EQ.INTVAL) THEN
            TABID = CVAL(IS)
            CALL GET_INTERNAL_TAB_NUMBER
     I                                  (STDOUT, TABID,
     M                                   EFLAG2,
     O                                   ITEMP)

              THE_VALUES(KNT) = ITEMP
            ELSE
C             Invalid item found. 
              WRITE(STDOUT,50) CVAL(IS)(1:CLEN(IS))
              EFLAG2 = 1

              THE_VALUES(KNT) = 0
            ENDIF
            IS = IS + 1
            GOTO 133
          ENDIF
        ENDIF


200   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE GET_PUTGET_OPTIONS
     I                (STDIN, STDOUT, LINE, COMMAND,
     M                  EFLAG)

C     Get the standard set of put/get function table options
C     for a variety of commands.  The option values are 
C     all table numbers and are stored in common block: PUTGET.
C     

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER COMMAND, STDIN, STDOUT, EFLAG

      CHARACTER LINE*(*)

      INCLUDE 'putget.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  NVAL, INTVAL, REAVAL, CONTINUATION_VALUE,
     A         CHRVAL, DPRVAL, EXACT_TYPE, LOWER_TYPE,
     B         N_SYMBOL
      PARAMETER(N_SYMBOL=7, NVAL=50, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, CONTINUATION_VALUE=5,
     B          EXACT_TYPE=0,LOWER_TYPE=1)
      INTEGER CULVERT, CHANRAT, EMBANKQ, ORIFICE
      PARAMETER (CULVERT=1, CHANRAT=2, EMBANKQ=3, ORIFICE=4)


      INTEGER EFLAG2, I, J, IS, IE, IT, IB, OPT, IP, 
     A                ITEM_KNT, ISTART, NUMBER_OF_VALUES, N
      INTEGER CLEN(NVAL), IVAL(NVAL), TERML(NVAL), TERMCLS(NVAL),
     A        ITEM_TYPE(NVAL), THE_VALUES(MXGET)
      REAL RVAL(NVAL)
      REAL*8 DPVAL(NVAL)
      CHARACTER CVAL(NVAL)*256, TERM(NVAL)*1,
     A          KEY*16, TABID*16


C     + + + SAVED VALUES + + +
      INTEGER  VALID_FOR(N_SYMBOL,4),
     B        START_INDEX(N_SYMBOL), END_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE SYMBOL_TABLE, VALID_FOR, START_INDEX, END_INDEX


C     External names

      INTEGER LENSTR
      EXTERNAL CLEAR_PUTGET, GET_TABID_VALUES, GETVAL, 
     A         LENSTR, LSTAB
 
      DATA  SYMBOL_TABLE
     1      /'PUTQ    ','PUTMF   ','PUTMF3  ','PUTY2   ','PUTY3   ',
     A       'GETQ    ','GETMF   '/

      INTEGER SI1, EI1, EI2
      PARAMETER( SI1= MXGET+7, EI1= MXGET+6, EI2=2*MXGET+7)
      DATA START_INDEX
     A     /1, 2, 3, 4, 5, 6, SI1/

      DATA END_INDEX
     A     /1, 2, 3, 4, 5, EI1, EI2/

      DATA ((VALID_FOR(I,J),J=1,4),I=1,7)/
     A   1, 1, 1, 1,
     B   0, 1, 1, 1,
     C   1, 0, 0, 0,
     D   1, 0, 0, 0,
     E   1, 0, 0, 0,
     F   1, 0, 0, 0,
     G   1, 0, 0, 0/

C     ****************************FORMATS*******************************
50    FORMAT(/,' *ERR:731* Name=',A,' is unknown as an input value.')
51    FORMAT(/,' *ERR:732* Name=',A,' not valid for current command.')
52    FORMAT(/,' *ERR:733* Expected no more than ',I3,' values for option ',
     A           '/,11X,A,' but found ',I3,' values instead.')
56    FORMAT(/,' Unable to continue due to previous errors.')
60    FORMAT(/,' Processing:',A)

C***********************************************************************
C     Clear the local error flag.
      EFLAG2 = 0

C     Clear any existing put/get option values.  A value of zero is the
C     default for each option.
      CALL CLEAR_PUTGET()

C     Select option for GETVAL to return BOTH the string value
C     and the converted value for numeric responses.  
C     In some cases the value will be recomputed. 
      OPT = 2

C     On entry, line will contain the user input line containing the
C     options, if any.  There may be none given.  The line as 
C     passed may not contain all of the options.  If not, then a 
C     continuation value will appear at the end of the line and the
C     remaining options will be on the subsequent line or lines of
C     the file attached to STDIN. 


      IT = LENSTR(LINE)
      IF(IT.EQ.0) IT = IT + 1
      WRITE(STDOUT,60) LINE(1:IT)

      LINE(IT+1:IT+1) = ''''

      CALL GETVAL
     I           (STDOUT, LINE, NVAL, OPT,
     O            ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O            EFLAG2, TERM, TERML, TERMCLS,
     O            ITEM_KNT)
      
      IF(EFLAG2.NE.0) THEN
        WRITE(STDOUT,56)
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

      IF(ITEM_KNT.EQ.0) THEN
C       No options given.  Default values accepted. 
        GOTO 500
      ELSE
C       One or more items present. The first item in the list should
C       be an option name.  Look it up in the symbol table and see
C       if it is valid.

        IS = 1
200     CONTINUE
          KEY = CVAL(IS)          
          CALL LSTAB
     I             (KEY, SYMBOL_TABLE, N_SYMBOL,
     O              IP)
          IF(IP.EQ.0) THEN
C           Item not found in table.
            WRITE(STDOUT,50) KEY(1:LENSTR(KEY))
            STOP 'Abnormal stop. Errors found.' 
          ELSE
C           Is the item valid for the command?
            IF(VALID_FOR(IP,COMMAND).EQ.0) THEN
C             No.
              WRITE(STDOUT,51) KEY(1:LENSTR(KEY))
              STOP 'Abnormal stop. Errors found.' 
            ELSE
C             We have a valid option here.
              ISTART = IS + 1 
              CALL   GET_TABID_VALUES(
     I                      STDIN, STDOUT, INTVAL, REAVAL, CHRVAL, 
     I                      CONTINUATION_VALUE, OPT, NVAL, 
     M                      ISTART, LINE, ITEM_TYPE, IVAL, RVAL, 
     M                      DPVAL, CVAL, CLEN, EFLAG2, TERM, TERML,
     M                      TERMCLS, ITEM_KNT,
     O                      NUMBER_OF_VALUES, THE_VALUES)
              EFLAG = EFLAG2
              N = END_INDEX(IP) - START_INDEX(IP)+ 1
              IF(NUMBER_OF_VALUES.GT.N) THEN
C               Too many reponses to an item
                WRITE(STDOUT,52)  N, KEY(1:LENSTR(KEY)),
     A                    NUMBER_OF_VALUES
                STOP 'Abnormal stop. Errors found.' 
              ELSE
C               Store the values in the overlay vector using 
C               the indices.  This will set the value of the
C               actual named get/put option.
                IB = START_INDEX(IP)
                IE = IB + NUMBER_OF_VALUES - 1
                DO 210 I=IB,IE
                  J = I - IB + 1
                  PGOVER(I) = THE_VALUES(J)
210             CONTINUE

                IS = ISTART
                IF(IS.GT.0) THEN
C                 Continue processing more options.
                  GOTO 200
                ELSE
C                 All options processed.  
                  GOTO 500
                ENDIF
              ENDIF            
            ENDIF
          ENDIF
      ENDIF
500   CONTINUE
      RETURN
      END

C
C
C
      SUBROUTINE READ_TABID_PLUS
     I                         (STDOUT, LINEA,
     O                          EFLAG, TABID, TABLE, FOLLOW)

C     Get the next table id and any following material

      IMPLICIT NONE
      INTEGER EFLAG, STDOUT, TABLE
      CHARACTER LINEA*(*), TABID*16, FOLLOW*(*)

C     External names
      INTEGER NONBLANK_NONZERO
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO

C     Local

      INTEGER IT, NXT
      CHARACTER LINE*120
C     *****************************FORMATS******************************
 50   FORMAT(/,' *ERR:746* Expected line with a table id but found:',
     A       /,5X,A)
 52   FORMAT(/,' *ERR:747* No = found as expected in input line:',/,A)
C***********************************************************************
      LINE = LINEA
C     Get the name field for the table id and discard.  Do a search for
C     the equal sign because some table ids might have no intervening
C     space after the equal. 
      
      CALL STRIP_L_BLANKS(
     M                      LINE)
      IF(LINE(1:3).NE.'TAB') THEN
        WRITE(STDOUT,50) LINE
        STOP 'Abnormal stop: errors found.' 
      ENDIF
      IT = INDEX(LINE, '=')
      IF(IT.EQ.0) THEN
        WRITE(STDOUT,52) LINE
        STOP 'Abnormal stop: errors found.' 
      ENDIF

      NXT = IT + 1
C     Get the TABID field
      CALL NXTTOK
     I           (LINE,
     M            NXT,
     O            TABID)

C     Get the following material.  Could be multiple items
C     that are processed later. 

      FOLLOW = LINE(NXT:)

      CALL STRIP_L_BLANKS(
     M                    FOLLOW)

C     Convert TABID to an internal number
      IF(TABID(1:1).EQ.'-') THEN
        NXT = -1
        TABID = TABID(2:16)
      ELSE
        NXT = 1
      ENDIF
      IF(NONBLANK_NONZERO(TABID).EQ.0) THEN
        TABLE = 0
      ELSE
        CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TABID,
     M                             EFLAG,
     O                             TABLE)
      ENDIF
      TABLE = NXT*TABLE
       
      RETURN
      END      
C
C
C
      SUBROUTINE READ_TABID
     I                     (STDOUT, LINEA, MATCH,
     O                      EFLAG, TABID, TABLE)

C     Read a table id from the line.  The name used to
C     reference the table id must match MATCH 


      IMPLICIT NONE
      INTEGER EFLAG, STDOUT, TABLE
      CHARACTER LINEA*(*), TABID*16, MATCH*(*)

C     External names
      INTEGER LENSTR, NONBLANK_NONZERO
      EXTERNAL LENSTR, STRIP_L_BLANKS, GET_INTERNAL_TAB_NUMBER,
     A         NONBLANK_NONZERO

C     Local

      INTEGER IT, NXT, L
      CHARACTER LINE*120
C     *****************************FORMATS******************************
 50   FORMAT(/,' *ERR:748* Expected line with ',A,' but found:',
     A       /,5X,A)
 52   FORMAT(/,' *ERR:749* No = found as expected in input line:',/,A)
C***********************************************************************
      LINE = LINEA
      L = LENSTR(MATCH)
C     Get the name field for the table id and discard.  Do a search for
C     the equal sign because some table ids might have no intervening
C     space after the equal. 
      
      CALL STRIP_L_BLANKS(
     M                      LINE)
      IF(LINE(1:L).NE.MATCH) THEN
        WRITE(STDOUT,50) MATCH, LINE
        STOP 'Abnormal stop: errors found.' 
      ENDIF
      IT = INDEX(LINE, '=')
      IF(IT.EQ.0) THEN
        WRITE(STDOUT,52) LINE
        STOP 'Abnormal stop: errors found.' 
      ENDIF

      NXT = IT + 1
C     Get the TABID field
      CALL NXTTOK
     I           (LINE,
     M            NXT,
     O            TABID)

C     Convert TABID to an internal number including the negative
C     sign that is sometimes used to signal some special action.
      IF(TABID(1:1).EQ.'-') THEN
        NXT = -1
        TABID = TABID(2:16)
      ELSE
        NXT = 1
      ENDIF
      IF(NONBLANK_NONZERO(TABID).EQ.0) THEN
        TABLE = 0
      ELSE
        CALL GET_INTERNAL_TAB_NUMBER
     I                              (STDOUT, TABID,
     M                               EFLAG,
     O                               TABLE)
      ENDIF
      TABLE = NXT*TABLE
       
      RETURN
      END      
C
C
C
      SUBROUTINE   TAB_IN_USE
     I                       (TABID,
     M                         EFLAG)

C     Report that a table id is already being used.
      
      IMPLICIT NONE

      INTEGER EFLAG
      CHARACTER TABID*16

      INCLUDE 'stdun.cmn'
C     *****************************FORMAT*******************************
50    FORMAT(/,'*ERR:750 Table id= ',A,
     A            ' is in use and is not available.')
C***********************************************************************
      EFLAG = 1
      WRITE(STD6,50) TABID
      RETURN
      END


