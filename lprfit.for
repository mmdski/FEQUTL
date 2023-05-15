C
C
C 
      SUBROUTINE SET_LPRFIT_ITEM_DEFAULTS()

C     Set the default values in the vectors used to
      IMPLICIT NONE
      INCLUDE 'lprfititem.cmn'

C***********************************************************************
C     Default for: TABID 
      LPRFITITEMCTAB(  1) = '    '           
C     Default for: TABLE - note # is ignored in the standard scanner 
      LPRFITITEMCTAB(  2) = '    '           
C     Default for: FIT_WITH
      LPRFITITEMCTAB(  3) = 'VLSPLINE'
C     Default for: CHK_OPTION
      LPRFITITEMCTAB(  4) = 'NATURAL'
C     Default for: LEFT_SLOPE
      LPRFITITEMCTAB(  5) = 'LINEAR'
C     Default for: RIGHT_SLOPE
      LPRFITITEMCTAB(  6) = 'LINEAR'

C     Default for:ARGFAC 
      LPRFITITEMFTAB(  1) = 1.0           
C     Default for:INFAC
      LPRFITITEMFTAB(  2) =  27.0
C     Default for:OUTFAC
      LPRFITITEMFTAB(  3) = 43560.0

C     Default for: ZONE
      LPRFITITEMCTAB(  7) = 'NONE'
C     Default for: HGRID
      LPRFITITEMCTAB(  8) = 'NONE'
C     Default for: VDATUM
      LPRFITITEMCTAB(  9) = 'NONE'
C     Default for: UNITSYS
      LPRFITITEMCTAB( 10) = 'NONE'
C     Default for: BASIS
      LPRFITITEMCTAB( 11) = 'NONE'

C     Default for: EASTING
      LPRFITITEMDTAB(  3) = -33d6
C     Default for: NORTHING
      LPRFITITEMDTAB(  4) =  -33d6
      

      RETURN
      END
C
C
C
      SUBROUTINE  SET_LPRFIT_ITEMS(
     M             EFLAG,
     O             TAB, FIT_WITH, CHK_OPTION, LEFT_SLOPE, 
     O             RIGHT_SLOPE, ARGFAC, INFAC, OUTFAC,
     O             zone, hgrid, vdatum, unitsys, basis, 
     O             easting, northing)

C     Set items in LPRFIT
C     All values not set explicitly by user are at their default value.

      IMPLICIT NONE

      INTEGER EFLAG, TAB
      REAL ARGFAC, INFAC, OUTFAC
      CHARACTER*16 FIT_WITH, CHK_OPTION, LEFT_SLOPE, RIGHT_SLOPE
      character*8 zone, hgrid, vdatum, unitsys, basis
      real*8 easting, northing
    

C     Local
      CHARACTER*16 KEY1, KEY2

      INCLUDE 'lprfititem.cmn'
      INCLUDE 'stdun.cmn'
C***********************************************************************
C     Set the table id.  There  are two strings allowed as the variable name:
C     TABID or TABLE
      KEY1 = LPRFITITEMCTAB(  1)
      KEY2 = LPRFITITEMCTAB(  2)
      IF(KEY1.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY1,
     M                                EFLAG,
     O                                TAB)
      ELSEIF(KEY2.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY2,
     M                                EFLAG,
     O                                TAB)
      ELSE
        TAB = 0
      ENDIF
C     Valuet for: FIT_WITH
      FIT_WITH = LPRFITITEMCTAB(  3) 
C     Value for: CHK_OPTION
      CHK_OPTION = LPRFITITEMCTAB(  4) 
C     Value for: LEFT_SLOPE
      LEFT_SLOPE = LPRFITITEMCTAB(  5) 
C     Value for: RIGHT_SLOPE
      RIGHT_SLOPE =  LPRFITITEMCTAB(  6)

C     Value for:ARGFAC 
      ARGFAC = LPRFITITEMFTAB(  1)
C     Value for:INFAC
      INFAC =  LPRFITITEMFTAB(  2) 
C     Value for:OUTFAC
      OUTFAC =  LPRFITITEMFTAB(  3) 
c     value for: ZONE
      zone = LPRFITITEMCTAB(  7)
c     value for: HGRID
      hgrid = LPRFITITEMCTAB(  8)
c     value for: VDATUM
      vdatum = LPRFITITEMCTAB(  9)
c     value for: UNITSYS
      unitsys = LPRFITITEMCTAB(  10)
c     value for: BASIS
      basis = LPRFITITEMCTAB(  11)

c     value for: EASTING
      easting = LPRFITITEMDTAB(  3)
c     value for: northing
      northing = LPRFITITEMDTAB(  4)

    
      RETURN
      END
C
C
C
      SUBROUTINE GET_LPRFIT_ITEMS(STDIN, STDOUT,
     M                   EFLAG)

C     Get the table id and various options for EMBANKQ command

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG

      INCLUDE 'lprfititem.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  INTVAL, REAVAL, CHRVAL, DPRVAL, EXACT, LOWER,
     B         NUMERIC, CHAR, N_SYMBOL, NXTBLK, ENDSIG, NONE
      PARAMETER(N_SYMBOL=23, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, EXACT=0, LOWER=1, 
     B          NUMERIC=0, CHAR=1, NXTBLK=2, ENDSIG=3,
     C          NONE=0)

      INTEGER MAX_LINE
     A        

      EXTERNAL GET_NAMED_ITEMS, SET_LPRFIT_ITEM_DEFAULTS

C     + + + SAVED VALUES + + +
      INTEGER GROUP(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL), GROUP_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE  SYMBOL_TABLE, GROUP, RESPONSE_TYPE, CONVERT_RULE,
     A      GROUP_INDEX

      DATA  SYMBOL_TABLE /
     A'TABID','TABLE','FIT_WITH','CHK_OPTION','LEFT_SLOPE',
     B'RIGHT_SLOPE','ARGFAC','INFAC','OUTFAC','DATA',
     C 'ARGUMENT','Argument','ELEVATION','Elevation','DEPTH',
     D 'Depth','ZONE','HGRID','VDATUM','UNITSYS',
     E 'EASTING','NORTHING','BASIS'/
                                                                      
      DATA GROUP  /
     A   6*CHAR, 3*NUMERIC,ENDSIG,6*NXTBLK,4*CHAR,2*NUMERIC,CHAR/          

      DATA GROUP_INDEX /
     A  1, 2, 3, 4, 5, 6, 1, 2, 3, 7*0, 7, 8, 9, 10,
     b  5, 7,11/
                                       
      DATA RESPONSE_TYPE  /
     A  6*CHRVAL,3*REAVAL,7*NONE,4*CHRVAL,2*DPRVAL,CHRVAL/
    
      DATA CONVERT_RULE /
     A  4*EXACT,2*LOWER,3*LOWER,7*EXACT,4*LOWER,2*LOWER,LOWER/
   
C***********************************************************************
C     Set Defaults
 
      CALL SET_LPRFIT_ITEM_DEFAULTS()

      MAX_LINE = 14
      CALL GET_NAMED_ITEMS(
     I                     STDIN, STDOUT,  MAX_LINE, N_SYMBOL,
     I  GROUP, RESPONSE_TYPE, CONVERT_RULE, GROUP_INDEX, SYMBOL_TABLE,
     I  MAXR_LPRFITITEM, MAXDP_LPRFITITEM, MAXC_LPRFITITEM, 
     I  'LPRFIT items',
     O  LPRFITITEMITAB, LPRFITITEMFTAB, LPRFITITEMDTAB, LPRFITITEMCTAB,
     O  EFLAG)
      
      RETURN

      END

C
C
C
      SUBROUTINE READ_LPRFIT_ITEMS(
     I                             STDOUT, LINE, NITEM, ITEM_START,
     I                             ITEM_END,
     M                             EFLAG,
     O                             ARG, NF, FVEC)

C     Get the items of data from a LPRFIT line
      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM), NF,
     A        EFLAG
      REAL*8 ARG, FVEC(10)
      CHARACTER LINE*120

C     Local

      INTEGER I, IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      INTEGER NONBLANK_NONZERO, LENSTR
      EXTERNAL STRIP_L_BLANKS, NONBLANK_NONZERO, LENSTR,
     A         GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:763* Only ',I3,' items given in ',
     A   'an LPRFIT input line.  Need at least 2 items.')
52    FORMAT(/,' *ERR:753* Conversion error in field ',I2,' in:',/,
     A        A)
C***********************************************************************
      IF(NITEM.LT.2) THEN
        WRITE(STDOUT,50) NITEM
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF

      N = 1
C     Process the argument value
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
       READ(TPC,*,ERR=999) ARG

C     Process the storage values- one or more
      
      NF = NITEM - 1
      DO 100 I=1,NF
        N = N + 1
        IS = ITEM_START(N)
        IE = ITEM_END(N)
        TPC = LINE(IS:IE)
        CALL STRIP_L_BLANKS(
     M                      TPC)
        READ(TPC,'(F18.0)',ERR=999) FVEC(I)
100   CONTINUE

      
      RETURN
999   CONTINUE
      WRITE(STDOUT,52) N, LINE
      STOP 'Abnormal stop.  Errors found.' 
      END


C
C
C
      SUBROUTINE   LPRFIT
     I                   (STDIN, STDOUT, STDTAB, FTP,
     M                    EFLAG)
 
C     + + + PURPOSE + + +
C     Fit a sequence of storage values to estimate a surface area for 
C     the capacity table of a level-pool reservoir.  

      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     FTP - pointer to next open location in the function-table storage system
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      include 'grid_datum.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADRS, I, J, N, NF, NITEM, MAXN, TABLE, LCODE, RCODE,
     A        NEXT, NTAB, LFLAG

      PARAMETER (MAXN=8)

      INTEGER ITEM_START(MAXN), ITEM_END(MAXN)

      REAL ARGFAC, INFAC, OUTFAC, YS(MNDEP), SS(MNDEP), AS(MNDEP),
     A     YOLD, SOLD, AOLD, DY, YNEW, SNEW, ANEW, DA

      REAL*8 Y(MNDEP), S(MNDEP), A(MNDEP), ARG, FVEC(10),
     A       OLDARG, SUM, RATIO, LVAL, RVAL, TERM, H1, H2, DS,
     B       IFAC(4), easting, northing

      CHARACTER*16 FIT_WITH, CHK_OPTION, LEFT_SLOPE, RIGHT_SLOPE,
     A             TABID
      CHARACTER LINE*120, JUST*5, AREA*10, VOLUME*10, ADJLOC(MNDEP)*1,
     a          zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
    
  
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      CHARACTER GET_TABID*16, PUT10D*10
      EXTERNAL  inline,  KIL, LENSTR, VAR_DECIMALD,
     A GET_LPRFIT_ITEMS, SET_LPRFIT_ITEMS, GET_TABID,
     B PUT1D, LKTAB, chk_vdatum_unitsys

      DATA IFAC/4.D0,2.D0,4.D0,1.D0/ 
C     + + + INPUT FORMATS + + +
 8    FORMAT(A6,1X,F10.0)
 16   FORMAT(A5,1X,F10.0)
 26   FORMAT(A5,1X,I5)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' TabId= ',A)
52    FORMAT(' Fitting with a cubic spline.')
54    FORMAT(' Fitting with a variation limited cubic spline.')
56    FORMAT(' Fitting with a piecewise cubic Hermite polynomial.')
58    FORMAT(/,'*ERR:764* ',A,' is an unknown fitting option.',/,
     A     9X,' Using VLSPLINE.')
60    FORMAT(/,' Variation-checking option is: ',A,/,
     A         ' Left-hand slope value or option is: ',A,/,
     B     ' Right-hand slope value or option is: ',A)
62    FORMAT(/,' Conversion factor for argument is:',F10.2,/,
     A         ' Conversion factor for incoming values is:',F10.2,/,
     B         ' Conversion factor for outgoing values is:',F10.2)
64    FORMAT(/,' *ERR:765* Number of tabulated values in LPRFIT=',
     A   I5,' > ',I5,/,' Increase value of MNDEP in file',
     B  ' arsize.prm and recompile.')
 65   format('ZONE=',a8,' HGRID=',a8,' VDATUM=',a8,' UNITSYS=',a8,
     a        ' BASIS=',a8,/, 'EASTING=',0PF15.3,' NORTHING=',F15.3)

66    FORMAT(/,'  Argument        Storage')
68    FORMAT(F10.2,10(1PE15.8))
70    FORMAT(/,'*ERR:766* Unable to read ',A,
     A   ' for left-hand slope value.')
72    FORMAT(/,'*ERR:767* Unable to read ',A,
     A   ' for right-hand slope value.')
74    FORMAT('TABID=',A16,/,'TYPE=   -4')
 76   format( 'REFL=       0.0 FAC=',F10.2,
     B /,' ELEVATION    VOLUME      AREA')
80    FORMAT(/,' *ERR:768* Storage=',F10.3,
     A      ' <  0.0 at argument=',F10.3)
82    FORMAT(/,' *ERR:769* Area=',F10.3,
     A      ' <  0.0 at argument=',F10.3)
84    FORMAT(/,' *ERR:770* At argument=',F10.3,' storage=',F10.3,
     A       ' does not increase.')
86    FORMAT(/,' *ERR:771* At argument=',F10.3,' area=',F10.3,
     A       ' does not increase.')
88    FORMAT(/,
     A' *BUG:XXX* Error in area in interval ending at argument=',F10.3)
 90   FORMAT('  Processing LPRFIT TabId= ',A)
92    FORMAT(/,' Detailed check of tabulated values.',
     A'  Three intermediate values are computed for checking.')
94    FORMAT('  Argument     Storage        Area  Delta Area')
96    FORMAT(F10.3,1PE12.4,1PE12.4,1PE12.3)
97    FORMAT(/,
     A' Capacity table not checked at user request.  May be invalid.')
98    FORMAT(/,
     A' The table for level-pool reservoir capacity appears valid.')
C***********************************************************************
C     Set the value for the temporary storage of function tables
      NEXT = FTP

C     Set the justification option for GET_ITEM_LIMITS
      JUST = 'RIGHT'

      CALL GET_LPRFIT_ITEMS(STDIN, STDOUT,
     M                    EFLAG)
      CALL  SET_LPRFIT_ITEMS(
     M             EFLAG,
     O             TABLE, FIT_WITH, CHK_OPTION, LEFT_SLOPE, 
     O             RIGHT_SLOPE, ARGFAC, INFAC, OUTFAC,
     O             zone, hgrid, vdatum, unitsys, basis, 
     O             easting, northing)

      call chk_vdatum_unitsys(stdout, vdatum, unitsys, 
     a          ' during output of table for LPR capacity')



      TABID = GET_TABID(TABLE)
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)) 

      WRITE(STDOUT,*) ' '
      IF(FIT_WITH.EQ.'CSPLINE') THEN
        WRITE(STDOUT,52)
      ELSEIF(FIT_WITH.EQ.'VLSPLINE') THEN
        WRITE(STDOUT,54)
      ELSEIF(FIT_WITH.EQ.'PCHERMITE') THEN
        WRITE(STDOUT,56)
      ELSE
        WRITE(STDOUT,58) FIT_WITH(1:LENSTR(FIT_WITH))
      ENDIF

      WRITE(STDOUT,60) CHK_OPTION, LEFT_SLOPE, RIGHT_SLOPE

      WRITE(STDOUT,62) ARGFAC, INFAC, OUTFAC
      RATIO = INFAC/OUTFAC
      WRITE(*,90) TABID
 
C     MAKE SURE TABLE NUMBER IS NOT ALREADY USED IN THIS INPUT
 
      IF(FTPNT(TABLE).NE.0) CALL TAB_IN_USE
     I                                     (TABID,
     M                                       EFLAG)
 
C     Get the heading line.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)


      I = 0
      OLDARG = -1.D100
      WRITE(STDOUT,66)

100   CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)

        CALL READ_LPRFIT_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG,
     O                         ARG, NF, FVEC)

        IF(EFLAG.EQ.0)  THEN

          IF(ARG.GT.OLDARG) THEN
            SUM = 0.D0
            DO 110 J=1,NF
              SUM = SUM + FVEC(J)
110         CONTINUE

            I = I + 1
            IF(I.GT.MNDEP) THEN
              WRITE(STDOUT,64) I, MNDEP
              STOP 'Abnormal stop.  Error(s) found.'
            ENDIF
            WRITE(STDOUT,68) ARG, (FVEC(J), J=1,NF) 
            Y(I) = (ARG )*ARGFAC
            S(I) = SUM*RATIO
            OLDARG = ARG
            GOTO 100
          ELSE
            N = I
          ENDIF
        ELSE
          RETURN
        ENDIF


C     Now that we finally have the data in hand, do the computations.

     

C     Compute the end conditions for the fit.
      LCODE = 1
      IF(LEFT_SLOPE.EQ.'LINEAR') THEN
        LVAL = (S(2) - S(1))/(Y(2) -  Y(1))
      ELSEIF(LEFT_SLOPE.EQ.'PARABOLIC') THEN
        H1 = Y(2) - Y(1)
        H2 = Y(3) - Y(2)
        LVAL = -(2*H1 + H2)*S(1)/(H1*(H1+H2))
     A         + (H1 + H2)*S(2)/(H1*H2)
     B         - H1*S(3)/((H1+H2)*H2)

      ELSE
C       Try to read the value to get the slope.
        READ(LEFT_SLOPE,'(F15.0)',ERR=998) LVAL
      ENDIF

      RCODE = 1
      IF(RIGHT_SLOPE.EQ.'LINEAR') THEN
        RVAL = (S(N) - S(N-1))/(Y(N) - Y(N-1))
      ELSEIF(RIGHT_SLOPE.EQ.'PARABOLIC') THEN
        H1 = Y(N-1) - Y(N-2)
        H2 = Y(N) - Y(N-1)
        RVAL = H2*S(N-2)/(H1*(H1+H2)) 
     A        - (H1 + H2)*S(N-1)/(H1*H2)
     B        + (2*H2 + H1)*S(N)/((H1+H2)*H2)

      ELSE
        READ(RIGHT_SLOPE,'(F15.0)',ERR=999) RVAL
      ENDIF

      IF(FIT_WITH.EQ.'VLSPLINE') THEN
        CALL VLCHPP
     I             (STDOUT, N, Y, S, LCODE, LVAL, RCODE, RVAL,
     O              A, ADJLOC)
      ELSEIF(FIT_WITH.EQ.'CSPLINE') THEN
        CALL SPLINE
     I              (STDOUT, Y, S, N, LCODE, LVAL, RCODE, RVAL,
     O               A)
      ELSEIF(FIT_WITH.EQ.'PCHERMITE') THEN
        A(1) = LVAL
        A(N) = RVAL
        DO 120 I=2,N-1
          H1 = Y(I) - Y(I-1)
          H2 = Y(I+1) - Y(I)
          A(I) =  -H2*S(I-1)/(H1*(H1+H2))
     A            - (H1 - H2)*S(I)/(H1*H2)
     B            + H1*S(I+1)/((H1+H2)*H2)
120     CONTINUE
      ENDIF

C     Store the table in the function-table system using temporary storage. 
C     That is, the table is only known in this routine and the storage used
C     is released on exit. 
      DO 125 I=1,N
        YS(I) = Y(I)
        SS(I) = S(I)
        AS(I) = A(I)
125   CONTINUE
      CALL PUT1D
     I          (STDOUT, -1, 4, N, YS, SS, AS,
     M           NEXT,
     O           ADRS)

C     Scan the table values checking validity for three cases:
C     Case 1: NATURAL
C           Assumes a natural reservoir site which implies the following:
C         1.  Storage is always positive (excluding initial point) 
C             and increases with increases in argument
C         2.  Area is always positive (excluding initial point) 
C             and increases with increases in argument 

C     Case 2: CONSTRUCTED
C           Assumes a constructed reservoir which implies the following:
C         1. Storage is always positive and increases with increase in argument
C         2. Area is always positive.
         
C     Case 3: NONE
C           No checking for variation. 

      YOLD = YS(1)
      CALL LKTAB(ADRS, YOLD, 0,
     O                  SOLD , NTAB, AOLD)
 
C     Clear the local error flag
      LFLAG = 0
      IF(CHK_OPTION.EQ.'NATURAL') THEN
        IF(SOLD.LT.0.0) THEN
          WRITE(STDOUT,80) SOLD, YOLD
          LFLAG = 1
        ENDIF
        IF(AOLD.LT.0.0) THEN
          WRITE(STDOUT,82) AOLD, YOLD
          LFLAG = 1
        ENDIF
      ENDIF

      WRITE(STDOUT,92)
      WRITE(STDOUT,94)
      DA = 0.0 
      WRITE(STDOUT,96) YOLD, SOLD, AOLD, DA
      DO 129 I=2,N

        DY = YS(I) - YS(I-1)
        DS = SS(I) - SS(I-1)
        SUM = AOLD
        DO 127 J=1,4
          YNEW = YOLD + DY*FLOAT(J)/4.0

          CALL LKTAB(ADRS, YNEW, 0,
     O                      SNEW , NTAB, ANEW)
          DA = ANEW - AOLD
          WRITE(STDOUT,96) YNEW, SNEW, ANEW, DA
          IF(CHK_OPTION.NE.'NONE') THEN
C           Check for values
            IF(SNEW.LT.0.0) THEN
              WRITE(STDOUT,80) SNEW, YNEW
            LFLAG = 1
            ENDIF
            IF(ANEW.LE.0.0) THEN
              WRITE(STDOUT,82) ANEW, YNEW
              LFLAG = 1
            ENDIF

            IF(CHK_OPTION.EQ.'NATURAL') THEN
C             Check for variation.
              IF(SNEW.LE.SOLD) THEN
                WRITE(STDOUT,84) YNEW, SNEW
                LFLAG = 1
              ENDIF
              IF(ANEW.LE.AOLD) THEN
                WRITE(STDOUT,86) YNEW, ANEW
                LFLAG = 1
              ENDIF
            ELSEIF(CHK_OPTION.EQ.'CONSTRUCTED') THEN
C             Check for variation.
              IF(SNEW.LE.SOLD) THEN
                WRITE(STDOUT,84) YNEW, SNEW
                LFLAG = 1
              ENDIF
            ENDIF
          ENDIF
C         Update the sum for numerical integration
             
          SUM = SUM + IFAC(J)*ANEW
          SOLD = SNEW
          AOLD = ANEW

127     CONTINUE

C       Check area by integration
        SUM = SUM*DY/12.D0
        IF(ABS((SUM - DS)/DS).GT.0.001) THEN
          WRITE(STDOUT,88) YNEW
          STOP 'Abnormal stop.  Bug found.'
        ENDIF

        YOLD = YNEW
129   CONTINUE

      IF(LFLAG.EQ.0) THEN
        IF(CHK_OPTION.NE.'NONE') THEN
          WRITE(STDOUT,98)
        ELSE
          WRITE(STDOUT,97)
        ENDIF
      ELSE
        EFLAG = 1
      ENDIF

C     Output a table of type 4

      WRITE(STDTAB,74) TABID

      if(zone /= 'NONE' .and. zone /= 'NONE')  then
c       Output the location information.
        write(stdtab, 65) zone, hgrid, vdatum, unitsys, basis, easting, 
     a        northing
      endif

      write(stdtab,76) outfac

      DO 140 I=1,N
        CALL VAR_DECIMALD(S(I),
     O                     VOLUME)

        CALL VAR_DECIMALD(A(I),
     O                       AREA)
        WRITE(STDTAB,'(F10.3,A10,A10)') Y(I), VOLUME, AREA
140   CONTINUE
      IF(Y(N).GT.0.D0) THEN
        TERM = -1.D0
      ELSE
        TERM = Y(N) - 1.D0
      ENDIF
      WRITE(STDTAB,'(F10.0)') TERM


      RETURN
998   CONTINUE
       WRITE(STDOUT,70) LEFT_SLOPE
      STOP 'Abnormal stop. Error(s) found.'
999   CONTINUE
      WRITE(STDOUT,72) RIGHT_SLOPE
      STOP 'Abnormal stop. Error(s) found.'

      END



