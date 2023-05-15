C
C
C
      SUBROUTINE READ_CULVERT_ITEMS(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     M                            EFLAG, IDLEN, HLIDLEN, 
     O                            NODE, NAME, XTABN, XCHAR, 
     O                            ZCHAR, CA, CD, HL)

C     Get the items of data from a CULVERT barrel input line.  
      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A        NODE, HL, EFLAG, IDLEN, HLIDLEN
      REAL CA, CD
      CHARACTER LINE*120, NAME*8, XTABN*5, XCHAR*10, ZCHAR*10 

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:739* Only ',I3,' items given in ',
     A   'culvert-barrel description line.  Need at least five items.')

C***********************************************************************
      IF(NITEM.LT.5) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

      N = 1
C     Process the NODE
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      IF(TPC(1:1).EQ.' ') THEN
        NODE = 0
      ELSE
        READ(TPC,*) NODE
      ENDIF

C     Process the node id
      N = 2
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      NAME = TPC

C     Process the table id
      N = 3
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      KEY = TPC
      LKEY = LENSTR(KEY)
      IDLEN = MAX(IDLEN, LKEY)
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
        IF(KEY(1:1).NE.'-') THEN
C         We have an id given.
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 ITAB)
          WRITE(XTABN,'(I5)') ITAB
        ELSE
C         We have a string starting with a - here
          IF(LKEY.GT.1) THEN
C           It is a tabid prefixed with a minus sign.
C           Strip the minus, define the internal number, and apply
C           the minus to the internal number
            KEY = KEY(2:16)
            CALL GET_INTERNAL_TAB_NUMBER
     I                                  (STDOUT, KEY,
     M                                   EFLAG,
     O                                   ITAB)
            WRITE(XTABN,'(I5)') -ITAB
          ELSE
            XTABN = TPC
          ENDIF
        ENDIF
      ELSE
        XTABN = ' '   
      ENDIF
      
C     Process the station 
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      XCHAR = TPC

C     Process the invert elevation
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      ZCHAR = TPC

C     Process the local acceleration losses
      N = 6
      IF(N.GT.NITEM) THEN
        CA = 0.0
      ELSE
        IS = ITEM_START(N)
        IE = ITEM_END(N)
        TPC = LINE(IS:IE)
        CALL STRIP_L_BLANKS(
     M                      TPC)
        IF(TPC(1:1).EQ.' ') THEN
          CA = 0.0
        ELSE
          READ(TPC(1:5),'(F5.0)') CA
        ENDIF
      ENDIF

C     Process local decceleration losses
      N = 7
      IF(N.GT.NITEM) THEN
        CD = 0.0
      ELSE
        IS = ITEM_START(N)
        IE = ITEM_END(N)
        TPC = LINE(IS:IE)
        CALL STRIP_L_BLANKS(
     M                      TPC)
        IF(TPC(1:1).EQ.' ') THEN
          CD = 0.0
        ELSE
          READ(TPC(1:5),'(F5.0)') CD
        ENDIF
      ENDIF

C     Process local structure losses
      N = 8
      IF(N.GT.NITEM) THEN
        HL = 0
      ELSE
        IS = ITEM_START(N)
        IE = ITEM_END(N)
        TPC = LINE(IS:IE)
        CALL STRIP_L_BLANKS(
     M                      TPC)
        KEY = TPC
        IF(NONBLANK_NONZERO(KEY).GT.0) THEN
C         We have a reference to a loss table.
          HLIDLEN = MAX(HLIDLEN,LENSTR(KEY))
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 HL)
        ELSE
          HL = 0.0
        ENDIF
      ENDIF

      RETURN
      END
C
C
C
      SUBROUTINE SET_DEPITM_DEFAULTS()

C     Set the default values in the vectors used to process a
C     deparature reach table spec.
      IMPLICIT NONE
      INCLUDE 'depitm.cmn'

C***********************************************************************
C     Default for: DEPTAB
      DEPITMCTAB(  1) = '    '           
C     Default for: BEGTAB
      DEPITMCTAB(  2) = '    '
C     Default for: RMFFAC
      DEPITMFTAB(  1) = 1.0
      RETURN
      END
C
C
C
      SUBROUTINE  SET_DEPITM(
     O                       DEPTAB, BEGTAB, RMFFAC,
     M                       EFLAG)

C     Set items for a departure reach table spec.
C     All values not set explicitly by user are at their default value.
      IMPLICIT NONE

    
      INTEGER BEGTAB, DEPTAB, EFLAG 
      REAL RMFFAC
    

C     Local
      CHARACTER*16 KEY

      INCLUDE 'depitm.cmn'
      INCLUDE 'stdun.cmn'
C***********************************************************************
      KEY = DEPITMCTAB(  1)
      IF(KEY.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY,
     M                                EFLAG,
     O                                DEPTAB)
      ELSE
        DEPTAB = 0
      ENDIF
      KEY = DEPITMCTAB(  2)
      IF(KEY.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY,
     M                                EFLAG,
     O                                BEGTAB)
      ELSE
        BEGTAB = 0
      ENDIF
      RMFFAC = DEPITMFTAB(  1)
      RETURN
      END
C
C
C
      SUBROUTINE GET_DEPARTURE_ITEMS(STDIN, STDOUT,
     M                   EFLAG)

C     Get the values from departure reach table specification

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG

      INCLUDE 'depitm.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  INTVAL, REAVAL, 
     A         CHRVAL, DPRVAL, EXACT, LOWER,
     B         NUMERIC, CHAR, N_SYMBOL
      PARAMETER(N_SYMBOL=3, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, 
     B          EXACT=0,LOWER=1, NUMERIC=0, CHAR=1)

      INTEGER MAX_LINE
     A        

      EXTERNAL GET_NAMED_ITEMS, SET_DEPITM_DEFAULTS

C     + + + SAVED VALUES + + +
      INTEGER GROUP(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL), GROUP_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE  SYMBOL_TABLE, GROUP, RESPONSE_TYPE, CONVERT_RULE,
     A      GROUP_INDEX

      DATA  SYMBOL_TABLE /
     *'DEPTAB','BEGTAB','RMFFAC'/
                                                                      
      DATA GROUP  /
     *CHAR,CHAR,NUMERIC/          

      DATA GROUP_INDEX /
     *1,2,1/
                                       
      DATA RESPONSE_TYPE  /
     *CHRVAL,CHRVAL,REAVAL/
    
      DATA CONVERT_RULE /
     *LOWER,LOWER,LOWER/
   
C***********************************************************************
C     Set Defaults
 
      CALL SET_DEPITM_DEFAULTS()

      MAX_LINE = 1
      CALL GET_NAMED_ITEMS(
     I                     STDIN, STDOUT,  MAX_LINE, N_SYMBOL,
     I  GROUP, RESPONSE_TYPE, CONVERT_RULE, GROUP_INDEX, SYMBOL_TABLE,
     I  MAXR_DEPITM, MAXDP_DEPITM, MAXC_DEPITM, 'Dep. reach items',
     O  DEPITMITAB, DEPITMFTAB, DEPITMDTAB, DEPITMCTAB,
     O  EFLAG)
      
      RETURN

      END
