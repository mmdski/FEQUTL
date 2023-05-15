C
C
C
      SUBROUTINE TYPE10_OUT(
     I                      STDOUT, STDTAB, TABNUM, N_ROW_ARGS, 
     I                      N_COL_ARGS, ROW_ARGS, HDATUM, COL_ARGS,
     I                      BODY, HDIN, HDITEM, HDOUT, ROWIN, ROWLAB,
     I                      ROWOUT, LABEL)


C     Output a function table of type 10.

      IMPLICIT NONE

C     Dummy arguments

      INCLUDE 'arsize.prm'

      INTEGER STDOUT, STDTAB, TABNUM, N_ROW_ARGS, N_COL_ARGS

      REAL HDATUM

      CHARACTER  ROW_ARGS(N_ROW_ARGS)*8, COL_ARGS(N_COL_ARGS)*8, 
     A     BODY(MRDT10,N_COL_ARGS)*8, HDIN*24, HDITEM*24, HDOUT*24,
     A     ROWIN*24, ROWLAB*8, ROWOUT*24, LABEL*50


C     Local

      INTEGER I, IS, J, N, TYPE

      CHARACTER LINE*120, TABID*16

C     Called subprograms

      CHARACTER GET_TABID*16
      INTEGER LENSTR

      EXTERNAL LENSTR
C     ********************************FORMATS***************************
50    FORMAT('TABID=',A)
52    FORMAT('TYPE=',I5)
54    FORMAT('HDATUM=',F10.3)
56    FORMAT('LABEL=',A)
58    FORMAT('  -1.0')
C***********************************************************************      

C     Output the table number
      TABID = GET_TABID(TABNUM)
      WRITE(STDTAB,50) TABID(1:LENSTR(TABID))

C     Construct table type followed by the formats.
      LINE = ' '
      TYPE = -10
      WRITE(LINE(1:10),52) TYPE

      IS = 12
      N = LENSTR(HDIN)
      LINE(IS:IS+N-1) = HDIN(1:N)
      IS = IS + N + 1
      N = LENSTR(HDOUT)
      LINE(IS:IS+N-1) = HDOUT(1:N)
      IS = IS + N + 1
      N = LENSTR(HDITEM)
      LINE(IS:IS+N-1) = HDITEM(1:N)
      IS = IS + N + 1
      N = LENSTR(ROWIN)
      LINE(IS:IS+N-1) = ROWIN(1:N)
      IS = IS + N + 1
      N = LENSTR(ROWOUT)
      LINE(IS:IS+N-1) = ROWOUT(1:N)
      IS = IS + N 
      WRITE(STDTAB,'(A)') LINE(1:IS)

C     Output the datum for heads
      WRITE(STDTAB,54) HDATUM

C     Output the table label
      WRITE(STDTAB,56) LABEL
C     Output the headings and the body of the table.  Note that the 
C     output is in character-string format.  The source routine 
C     has the task of creating the correct strings and formats for
C     the table.

      N = LENSTR(HDIN)
      WRITE(STDTAB,HDIN(2:N-1)) ROWLAB, (COL_ARGS(J), J=1,N_COL_ARGS)
      DO 100 I=1,N_ROW_ARGS
        WRITE(STDTAB,HDIN(2:N-1)) ROW_ARGS(I), 
     A                     (BODY(I,J), J=1,N_COL_ARGS)
100   CONTINUE
      WRITE(STDTAB,58)

      RETURN
      END
C
C
C
      REAL FUNCTION GATE_RES(HG)

C     Compute the gate residual function for an underflow
C     gate. 

      REAL HG

C     Common block needed values

      INCLUDE 'invgate.cmn'

C     Local variables.
      INTEGER MFTAB, FTYPE_CODE
      REAL  DQED, DQEU, NEWHG, Q
      real*8 jtime

      DATA JTIME/0.0d0/, MFTAB/0/

C***********************************************************************
C     Find the flow at the current gate opening 

      CALL TDLK15
     I           (GATE_RES_STDOUT, GATE_RES_TABLE, MFTAB, JTIME, 
     I            GATE_RES_EDN, GATE_RES_EUP, HG,
     I            GATE_RES_HBASE,
     O            Q, DQED, DQEU, NEWHG, GATE_RES_FTYPE, FTYPE_CODE)

      GATE_RES = (GATE_RES_FLOW - Q)/GATE_RES_FLOW
C      WRITE(GATE_RES_STDOUT,*) ' GATE_RES:HG=',HG,' Q=',Q,
C     A    ' RES=',GATE_RES

      RETURN
      END
C
C
C
      SUBROUTINE   FIND_GATE_OPENING(
     I         STDOUT, UD_TABLE, DU_TABLE, EUP, EDN, GATE_FLOW,
     I         HG_MAX, HBASE,
     O         HG, P, FTYPE, RESULT_FLAG, FTYPE_CODE)

C     Find the opening for an underflow gate so that the flow
C     through the gate will match the GATE_FLOW given in the
C     dummy argument list.   We assume that the GATE_FLOW differs
C     from zero.

      IMPLICIT NONE

      INTEGER STDOUT, UD_TABLE, DU_TABLE, RESULT_FLAG, FTYPE_CODE

      REAL EUP, EDN, GATE_FLOW, HG, P, HBASE, HG_MAX

      CHARACTER*8 FTYPE

C     Definition of dummy arguments

C     STDOUT- Fortran unit for standard output
C     UD_TABLE- address of type 15 table describing flow through 
C       the gate from upstream node  to downstream node.
C     DU_TABLE-address of type 15 table describing flow through
C       the gate from downstream node to upstream node.
C     EUP- water-surface elevation at upstream node
C     EDN- water-surface elevation at downstream node
C     GATE_FLOW- target flow for the gate
C     HG_MAX- maximum gate opening in the gate-flow tables.
C     HBASE- datum for head for the gate
C     HG- gate opening to match the target flow.
C     P- gate opening as fraction of maximum opening.
C     RESULT_FLAG- 1 if a solution exists, 0 otherwise. 


C     Common blocks

      INCLUDE 'epscom.cmn'
      INCLUDE 'invgate.cmn'

C     Local variables

      INTEGER MFTAB, FLAG

      REAL  DQED, DQEU, NEWHG, Q, A, B, HG_RESULT
      real*8 jtime


C     Called program units

      REAL GATE_RES
      EXTERNAL GATE_RES

C     *************************FORMATS**********************************
50    FORMAT(/,' *ERR:XXX* No sign change in REGFAL in subroutine',
     A       ' FIND_GATE_OPENING.') 
52    FORMAT(/,' *ERR:XXX* No solution in REGFAL in subroutine',
     A       ' FIND_GATE_OPENING after 100 iteration.')
C***********************************************************************
C     Clear the multiplying-factor table address and the time 
C     value.
      MFTAB = 0
      jTIME = 0.d0

C     Determine if a solution exists.  Set the gate to its maximum
C     opening with the given heads and see if the flow exceeds the
C     target value.  If not, no solution is possible.  


      IF(GATE_FLOW.GT.0.0) THEN
        CALL TDLK15
     I             (STDOUT, UD_TABLE, MFTAB, JTIME, EDN, EUP, HG_MAX,
     I              HBASE,
     O              Q, DQED, DQEU, NEWHG, FTYPE, FTYPE_CODE)
      ELSE
        CALL TDLK15
     I             (STDOUT, DU_TABLE, MFTAB, JTIME, EUP, EDN, HG_MAX,
     I              HBASE,
     O              Q, DQED, DQEU, NEWHG, FTYPE, FTYPE_CODE)
      ENDIF
C      WRITE(STDOUT,*) ' FIND_GATE_OPENING: Q=',Q,' GATE_FLOW=',
C     A                  GATE_FLOW 
      IF(ABS(Q).LT.ABS(GATE_FLOW)) THEN
C       No solution possible.
        HG = HG_MAX
        P = 1.0
        RESULT_FLAG = 0
        GOTO 9000
      ENDIF

C     Solution appears to be possible here.  We will solve for the 
C     root using modified Regula Falsi because we have a bound on
C     each end for a possible root and because the derivatives of
C     the flow from the gate can be discontinuous.  

      A = 0.0
      B = HG_MAX
     
C     Set values in the common block for function GATE_RES
      GATE_RES_FLOW = ABS(GATE_FLOW)
      GATE_RES_HBASE = HBASE
      GATE_RES_STDOUT = STDOUT
      IF(GATE_FLOW.GT.0.0) THEN
        GATE_RES_TABLE = UD_TABLE
        GATE_RES_EUP = EUP
        GATE_RES_EDN = EDN
      ELSE
        GATE_RES_TABLE = DU_TABLE
        GATE_RES_EUP = EDN
        GATE_RES_EDN = EUP
      ENDIF

      CALL   REGFAL
     I              (EPSARG, EPSF, GATE_RES,
     M               A, B,
     O               HG_RESULT, FLAG)

      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT,50) 
        STOP ' Abnormal stop.  Errors found.' 
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,52)
        WRITE(STDOUT,*) ' MIN_HG=',A,' MAX_HG=',B 
        STOP ' Abnormal stop. Errors found.' 
      ENDIF

      RESULT_FLAG = 1
      HG = HG_RESULT
      P = HG/HG_MAX
      FTYPE = GATE_RES_FTYPE

9000  CONTINUE
      RETURN
      END

C
C
C
      SUBROUTINE   GET_INV_GATE(STDIN, STDOUT,
     M   EFLAG,
     O   DU_TABLE, UD_TABLE, N_COL_BDYS, N_ROW_BDYS, CONTROL_TAB,
     O   COL_BDYS, ROW_BDYS, FLOOD_ELEV, FLOOD_FLOW,
     O   CPNT_ELEV_TS, UPS_ELEV_TS, DNS_ELEV_TS, GATE_FLOW_TS, 
     O   DRAIN_LOC, MIN_FLOW, OUTPUT_LEVEL, REVERSE_FLOW)

C     Get the values needed for the INV_GATE command

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG
      INTEGER DU_TABLE, UD_TABLE, N_COL_BDYS, N_ROW_BDYS, CONTROL_TAB

      REAL COL_BDYS(MCDT10), ROW_BDYS(MRDT10), FLOOD_ELEV, FLOOD_FLOW,
     A     MIN_FLOW

      CHARACTER CPNT_ELEV_TS*64, UPS_ELEV_TS*64, DNS_ELEV_TS*64, 
     A          GATE_FLOW_TS*64, DRAIN_LOC*3, OUTPUT_LEVEL*3,
     B          REVERSE_FLOW*3

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  NVAL, INTVAL, REAVAL, CONTINUATION_VALUE,
     A         CHRVAL, DPRVAL, EXACT_TYPE, LOWER_TYPE,
     B         N_SYMBOL
      PARAMETER(N_SYMBOL=15, NVAL=98, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, CONTINUATION_VALUE=5,
     B          EXACT_TYPE=0,LOWER_TYPE=1)

      INTEGER EFLAG2, I, IT, OPT, IP, SELECTION, ITEM_KNT
      INTEGER CLEN(NVAL), IVAL(NVAL), TERML(NVAL), TERMCLS(NVAL),
     A        ITEM_TYPE(NVAL)
      REAL RVAL(NVAL)
      REAL*8 DPVAL(NVAL)
      CHARACTER CVAL(NVAL)*256, TERM(NVAL)*1, LINE*120,
     A          KEY*16

      INTEGER LENSTR

      EXTERNAL LENSTR, GET_MULTIPLE_REAL_VALUES, GETVAL, 
     A         CHK_AND_CONVERT_RESPONSE, LSTAB

C     + + + SAVED VALUES + + +
      INTEGER SYMBOL_VALUE(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE SYMBOL_VALUE, SYMBOL_TABLE


      DATA  SYMBOL_TABLE
     1      /'CPNT_ELEV', 'UPS_ELEV', 'DNS_ELEV', 'GATE_FLOW',
     A       'DRAIN_LOC', 'DU_TABLE', 'UD_TABLE', 'CONTROL_TAB',
     B       'COL_BDYS', 'ROW_BDYS', 'FLOOD_ELEV', 'FLOOD_FLOW',
     C       'MIN_FLOW', 'OUTPUT_LEVEL', 'REVERSE_FLOW'/

      DATA SYMBOL_VALUE
     1      /         1,         2,         3,         4,         5,
     A                6,         7,         8,         9,        10,
     B               11,        12,        13,        14,         15/

      DATA RESPONSE_TYPE
     A     /CHRVAL, CHRVAL, CHRVAL, CHRVAL, CHRVAL, INTVAL,
     B            INTVAL, INTVAL, -1, -1, DPRVAL, DPRVAL, DPRVAL,
     C            CHRVAL, CHRVAL/
      DATA CONVERT_RULE
     A     /5*EXACT_TYPE, 5*EXACT_TYPE, 3*LOWER_TYPE, 2*EXACT_TYPE/
C     *****************************FORMATS******************************
 52   FORMAT(/,' *BUG:XXX* Invalid index=',I5,' for name=',A,' in',
     A       ' subroutine GET_INV_GATE.')
54    FORMAT(/,' *ERR:XXX* Name=',A16,' is unknown INV_GATE.')
56    FORMAT(/,' Unable to continue due to previous errors.')
60    FORMAT(/,' Processing:',A)
62    FORMAT(/,' Seeking additional values.')
C***********************************************************************
C     Clear the local error flag for subroutine GETVAL
      EFLAG2 = 0
C     Read lines of input and process each one until the expected number
C     of lines or an end of header signal is found.  


C     SET DEFAULTS
      DU_TABLE = 0
      UD_TABLE = 0
      N_COL_BDYS = 0
      N_ROW_BDYS = 0
      CONTROL_TAB = 0
      FLOOD_ELEV = 0.0
      FLOOD_FLOW = 0.0
      CPNT_ELEV_TS = ' '
      UPS_ELEV_TS = ' '
      DNS_ELEV_TS = ' '
      GATE_FLOW_TS= ' '
      DRAIN_LOC = 'DNS'
      MIN_FLOW = 0.0
      OUTPUT_LEVEL = 'MIN'
      REVERSE_FLOW ='NO '
      
C     Select option for GETVAL to return BOTH the string value
C     and the converted value for numeric responses.  
C     In some cases the value will be recomputed. 
      OPT = 2
C     Start a loop over input lines
100   CONTINUE


        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)

        IT = LENSTR(LINE)
        WRITE(STDOUT,60) LINE(1:IT)

        LINE(IT+1:IT+1) = ''''

        CALL GETVAL
     I             (STDOUT, LINE, NVAL, OPT,
     O              ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O              EFLAG2, TERM, TERML, TERMCLS,
     O              ITEM_KNT)
      WRITE(STDOUT,*) 
     A  ' Return from GETVAL in GET_INV_GATE: ITEM_KNT=',ITEM_KNT
      WRITE(STDOUT,97)
97    FORMAT(1X,12X,'ITEM','   LEN  TYPE T  TCLS')
      DO 9213 I=1,ITEM_KNT
        WRITE(STDOUT,99) CVAL(I), CLEN(I), ITEM_TYPE(I), TERM(I),
     A                  TERMCLS(I)
99    FORMAT(' ',A16,' ',I5,' ',I5,' ',A1,' ',I5)
9213  CONTINUE

        
        IF(EFLAG2.NE.0) THEN
C         Error in parsing the line of input.
          
          WRITE(STDOUT,56) 
          STOP 'Abnormal stop.  Errors found.' 
        ELSE
C         No errors reported.  Process the items found on the current
C         line. 

C         Check for the end of input
          IF(CVAL(1).EQ.'END') THEN
C           end of input. we are done.
            RETURN
          ENDIF
          
          CALL CHK_AND_CONVERT_RESPONSE(STDOUT,
     I           N_SYMBOL, SYMBOL_TABLE, RESPONSE_TYPE,
     I           CONVERT_RULE, ITEM_KNT, ITEM_TYPE, CVAL, CLEN,
     M           EFLAG,
     O           IVAL, RVAL, DPVAL)

          I = 1
110       CONTINUE          
C           Find the value for the next item from the 
C           symbol table. 
            KEY = CVAL(I)(1:CLEN(I))
            CALL LSTAB
     I                (KEY, SYMBOL_TABLE, N_SYMBOL,
     O                 IP)
            IF(IP.EQ.0) THEN
C             error-symbol not found
              WRITE(STDOUT,54) KEY
              EFLAG = 1

            ELSE
              SELECTION = SYMBOL_VALUE(IP)            
            ENDIF

            GOTO(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
     A          15),SELECTION

              WRITE(STDOUT,52) I, KEY
              STOP 'Abnormal stop. Errors found.'
 
 1          CONTINUE
C             File name for the series at the control point.   
              CPNT_ELEV_TS = CVAL(I+1)
              I = I + 2                                          
              GOTO 149
 2          CONTINUE
C             File name for the series at the ups node for the gate
              UPS_ELEV_TS = CVAL(I+1)
              I = I + 2
              GOTO 149
 3          CONTINUE
C             File name for the series at the dns node for the gate
              DNS_ELEV_TS = CVAL(I+1)
              I = I + 2
              GOTO 149
 4          CONTINUE
C             File name for the series giving the flow through the gate.
              GATE_FLOW_TS = CVAL(I+1)
              I = I + 2
              GOTO 149
 5          CONTINUE
C             Location of node for the reservoir
              DRAIN_LOC = CVAL(I+1)
              I = I + 2
              GOTO 149
 6          CONTINUE
C             Table number for dns to upstream flow through the gate
              DU_TABLE = IVAL(I+1)
              I = I + 2
              GOTO 149
 7          CONTINUE
C             Table number for  ups to dns flow through the gate
              UD_TABLE = IVAL(I+1)
              I = I + 2
              GOTO 149
 8          CONTINUE
C             Id number for the control table.
              CONTROL_TAB = IVAL(I+1)
              I = I + 2
              GOTO 149
 9          CONTINUE
C             Process the column arguments.
              CALL  GET_MULTIPLE_REAL_VALUES(
     I          STDIN, STDOUT, INTVAL, REAVAL, CONTINUATION_VALUE, 
     I          OPT, NVAL, 
     M          LINE, ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN,
     M          EFLAG2, TERM, TERML, TERMCLS, ITEM_KNT,
     O          N_COL_BDYS, COL_BDYS)

C               Force next line
                I = 2048
              GOTO 149
10          CONTINUE
C             Process the row arguments
              CALL  GET_MULTIPLE_REAL_VALUES(
     I          STDIN, STDOUT, INTVAL, REAVAL, CONTINUATION_VALUE, 
     I          OPT, NVAL, 
     M          LINE, ITEM_TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN,
     M          EFLAG2, TERM, TERML, TERMCLS, ITEM_KNT,
     O          N_ROW_BDYS, ROW_BDYS)
              
C               Force next line
                I = 2048
              GOTO 149
11          CONTINUE
C             Elevation at control point at beginning of flooding. 
              FLOOD_ELEV = RVAL(I+1)
              I = I + 2
              GOTO 149
12          CONTINUE
C             Flow at the control point at beginning of flooding
              FLOOD_FLOW = RVAL(I+1)
              I = I + 2
              GOTO 149
13          CONTINUE
C             Minimum flow that is treated as non-zero.  All flows
C             less than this value are set to zero before solving
C             for the gate setting. 
              MIN_FLOW = RVAL(I+1)
              I = I + 2
              GOTO 149
14          CONTINUE
C             Output level
              OUTPUT_LEVEL = CVAL(I+1)
              I = I + 2
              GOTO 149
15          CONTINUE
C             Reverse flow signal
              REVERSE_FLOW = CVAL(I+1)
              I = I + 2
              GOTO 149

 149        CONTINUE
              IF(I.GT.ITEM_KNT) THEN
C               Get the next line from the input
                GOTO 100
              ELSE
C               Get the next item from the current line
                  GOTO 110
              ENDIF            
        ENDIF

      END
C
C
C
      INTEGER FUNCTION GET_TYPE13_ADDRESS(TY15_ADRS)

C     Get the address of the first type 13 table in a table of
C     type 15 given its address.

      INTEGER TY15_ADRS

      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'

C     Local

      INTEGER IP
C***********************************************************************
      IP = TY15_ADRS
      GET_TYPE13_ADDRESS = ITAB(IP+5+1)
      RETURN
      END
C
C
C
      REAL FUNCTION GET_MAX_GATE(TY15_ADRS)

C     Get the maximum gate opening.

      INTEGER TY15_ADRS

      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'

C     Local

      INTEGER IP
C***********************************************************************
      IP = ITAB(TY15_ADRS)
      GET_MAX_GATE = FTAB(IP)
      RETURN
      END

      
C
C
C
      SUBROUTINE   INV_GATE(GRAV, STDIN, STDOUT, STDTAB,
     M                    EFLAG, TABDIR)
 
C     Inverts an underflow gate to create a control table to operate
C     the gate.  Gate inversion computes the gate opening to match
C     the flow produced by some fixed structure so that the gate
C     operation will approximate the fixed structure.  This then
C     creates a first approximation to the gate control tables. 
 
      IMPLICIT NONE

C     Dummy arguments
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL GRAV
 
C     DEFINITIONS 
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TABDIR - Table directory to remember table numbers
 
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'

C     Define key local values:

C     UD_TABLE= table number giving the type 15 rating table for the
C                gate for flow from the upstream node to the downstream
C                node.
C     DU_TABLE= table number giving the type 15 rating table for
C                the gate for flow from the downstream node to the 
C                upstream node. 
C     CPNT_ELEV_TS = file name for the time series of elevation at the
C                 control point for the structure.  The control
C                 is the exterior node giving the elevation 
C                 defining the operation of the gate.
C     UPS_ELEV_TS = file name for the time series of elevation at the
C                node upstream of the gate.  Upstream is defined
C                by the user with the rule that the flow, given
C                in another time series, must always be positive
C                for flow from the upstream node to the downstream
C                node.  Given only when the control-point location
C                is different than the upstream-node location.  
C     DNS_ELEV_TS = file name for the time series of elevation at the 
C                node downstream of the gate. 
C     GATE_FLOW_TS = file name for the time series of flows through the gate.

C     FLOOD_ELEV = the elevation at the control point that defines
C                  the zero point on the sequence of arguments for
C                  rows of the control table.  This is the elevation
C                  at the control point that often signals flood 
C                  hazard at some point downstream.  
C     FLOOD_FLOW = the flow at the control point when the elevation
C                  is at FLOOD_ELEV.  Used to estimate the gate
C                  operation for draining of the reservoir. The
C                  difference between the FLOOD_FLOW and the current
C                  flow at the cotnrol point gives the maximum 
C                  flow release possible. 
C     CONTROL_TAB = the table number for the control table computed
C                    by the GATE_INV command.
C     COL_BDYS = the  sequence of  boundary values defining  the cells  
C                for the columns of the type 10 table.  The midpoint of the
C                cell will become the argument value in the table.     
C                downstream node for the gate.  
C     ROW_BDYS = the  sequence of  boundary values defining  the cells 
C                for the rows of the type 10 table.  The midpoint of the
C                cell will become the argument value in the table. 
C       Note: The combination of ROW_BDYS and COL_BDYS defines an
C             array of relative gate openings at the intersection of 
C             each pair of  arguments.  A gate opening of 0.0
C             denotes a closed gate and a gate opening of 1.0 denotes
C             a gate fully open.  The gate openings for points intermediate
C             to those tabulated are defined by linear interpolation.

C     DRAIN_LOC = the location defining the node that represents the
C                 the reservoir to be drained.  Has two values UPS or
C                 DNS.

C     MIN_FLOW =  all flows less than this value are treated as zero flow.
C                 
C     OUTPUT_LEVEL= user control on output level:  MIN gives the minimum
C                   level with the summary tables only.  MAX gives the 
C                   results for each flow greater than MIN_FLOW.
C     REVERSE_FLOW= if YES, reverse flows are included.  If NO reverse
C                   flows are excluded.  If YES, MIN_FLOW refers to the
C                   absolute value of the flow. 
C     N_PER_CELL(*,*)- number of cases per cell.
C     SUM_PER-CELL(*,*) - sum of the gate openings in a cell.
C     MAX_PER_CELL(*,*) - maximum gate opening in a cell.
C     MIN_PER_CELL(*,*) - minimum gate opening in a cell.

      INTEGER DU_TABLE, UD_TABLE, N_COL_BDYS, N_ROW_BDYS, CONTROL_TAB,
     A        CPNT_UNIT, UPS_UNIT, DNS_UNIT, GATE_UNIT, ICODE,
     B        UPS_GIVEN, RESULT_FLAG, MJD, YR, MN, DY,
     C        N_PER_CELL(MRDT10,MCDT10), I, J, KNT, FTYPE_CODE

      REAL COL_BDYS(MCDT10+1), ROW_BDYS(MRDT10+1), FLOOD_ELEV, 
     A     FLOOD_FLOW, CPNT_ELEV, UPS_ELEV, DNS_ELEV, GATE_FLOW, 
     B     MIN_FLOW, HDATUM, HU, HD, HDIFF, HG_MAX, HG, P, HR,
     C     SUM_PER_CELL(MRDT10,MCDT10),
     D     MAX_PER_CELL(MRDT10,MCDT10),
     E     MIN_PER_CELL(MRDT10,MCDT10)


      CHARACTER CPNT_ELEV_TS*64, UPS_ELEV_TS*64, DNS_ELEV_TS*64, 
     A          GATE_FLOW_TS*64, DRAIN_LOC*3, FTYPE*8,
     B          PC_PER_CELL(MRDT10,MCDT10)*8,
     C          MAXC_PER_CELL(MRDT10,MCDT10)*6,
     D          MINC_PER_CELL(MRDT10,MCDT10)*6,
     E          NC_PER_CELL(MRDT10,MCDT10)*6,
     F          OUTPUT_LEVEL*3, REVERSE_FLOW*3,
     G          ROW_ARGS(MRDT10)*8, COL_ARGS(MCDT10)*8, 
     H          HDIN*24, HDITEM*24, HDOUT*24,
     I          ROWIN*24, ROWLAB*8, ROWOUT*24, LABEL*50,
     J          CHR3*3, CHR3_1*3



      LOGICAL THERE

      REAL*8 CPNT_JT, UPS_JT, DNS_JT, GATE_JT, EPS, FRAC

C     External program units

      INTEGER GET_UNIT, GET_TYPE13_ADDRESS

      REAL GETHDD, GET_MAX_GATE, TP

      EXTERNAL GET_INV_GATE, GET_UNIT, CHKTAB, CKTY15,
     A         GET_MAX_GATE, GET_TYPE13_ADDRESS, FIND_GATE_OPENING

      DATA EPS/1.D-6/

C     ***********************************FORMATS************************
50    FORMAT(/,' Elevation time series file for control point=',/,5X,
     A     A,/,' not found.  Please check spelling of name and path.')
51    FORMAT(/,' Elevation time series file for downstream point=',/,5X,
     A     A,/,' not found.  Please check spelling of name and path.')
52    FORMAT(/,' Flow time series file for the underflow gate=',/,5X,
     A     A,/,' not found.  Please check spelling of name and path.')
53    FORMAT(/,' Elevation time series file for upstream point=',/,5X,
     A     A,/,' not found.  Please check spelling of name and path.')
54    FORMAT(/,' *ERR:XXX* Time mismatch in time series:'/,
     A    11X,' Control point time series JT=',F20.10,/,
     B    11X,' Ups point time series JT=',F20.10,/,
     C    11X,' Dns point time series JT=',F20.10,/,
     D    11X,' Gate flow time series JT=',F20.10)
56    FORMAT(I5,I5,I3,I3,F10.2,F10.3,F10.3,F10.1,F10.3,F10.3,A8,I5)
57    FORMAT(I5,I3,I3,F10.2,F10.3)
58    FORMAT(/,
     A  '  KNT Year Mn Dy      Hour  Ups head  Dns head      Flow',
     B  '   Opening  Fraction FlowTyp Flag' )
60    FORMAT(/,' Datum for heads=',F10.3,' Maximum gate opening=',F8.3)
62    FORMAT(/,' *ERR:XXX* Ups level',F10.3,' falls outside the',
     A         ' min and max values of:',2F10.3)
64    FORMAT(/,' *ERR:XXX Head difference=',F10.3,' falls outside the',
     A         ' min and max values of:',2F10.3)
66    FORMAT(/,' Number of occurrences per cell:')
68    FORMAT(25I6)
70    FORMAT(/,' Mean gate opening per cell:')
72    FORMAT(F7.2,1X,25A6)
73    FORMAT(8X,25A6)
74    FORMAT(/,' Summary of Results')
75    FORMAT(' Cell  ','  Type ',24F6.1)
76    FORMAT(' Bdys  ','       ',24F6.1)
C***********************************************************************


      CALL   GET_INV_GATE(STDIN, STDOUT,
     M   EFLAG,
     O   DU_TABLE, UD_TABLE, N_COL_BDYS, N_ROW_BDYS, CONTROL_TAB,
     O   COL_BDYS, ROW_BDYS, FLOOD_ELEV, FLOOD_FLOW,
     O   CPNT_ELEV_TS, UPS_ELEV_TS, DNS_ELEV_TS, GATE_FLOW_TS, 
     O   DRAIN_LOC, MIN_FLOW, OUTPUT_LEVEL, REVERSE_FLOW)


C     Check the values so far.

      WRITE(STDOUT,61) DU_TABLE, UD_TABLE, CONTROL_TAB,
     A                 FLOOD_ELEV, FLOOD_FLOW, CPNT_ELEV_TS,
     B                 UPS_ELEV_TS, DNS_ELEV_TS, GATE_FLOW_TS,
     C                 DRAIN_LOC, MIN_FLOW, OUTPUT_LEVEL,
     D                 REVERSE_FLOW
61    FORMAT(/,' DU_TABLE=',I6,' UD_TABLE=',I6,' CONTROL_TAB=',I6,
     A         ' FLOOD_ELEV=',F8.3,' FLOOD_FLOW=',F8.1,/,
     B         ' CPNT_ELEV_TS=',A,/,' UPS_ELEV_TS=',A,/,
     C         ' DNS_ELEV_TS=',A,/,
     C         ' GATE_FLOW_TS=',A,/,' DRAIN_LOC=',A,
     D         ' MIN_FLOW=',F10.2,' OUTPUT_LEVEL=',A,
     E         ' REVERSE_FLOW=',A)

      WRITE(STDOUT,63) (COL_BDYS(J), J=1,N_COL_BDYS)
63    FORMAT(/,' COL_BDYS=',(10F10.2))
      WRITE(STDOUT,65) (ROW_BDYS(J), J=1,N_ROW_BDYS)
65    FORMAT(/,' ROW_BDYS=',(10F10.2))

C     Initialize the arrays used for counting and recording.
      DO 110 I=1, N_ROW_BDYS-1
        DO 100 J=1, N_COL_BDYS-1
          N_PER_CELL(I,J) = 0
          SUM_PER_CELL(I,J) = 0.0
          MAX_PER_CELL(I,J) = 0.0
          MIN_PER_CELL(I,J) = 1.0
100     CONTINUE
110   CONTINUE


C     Process the time series references

      INQUIRE(FILE=CPNT_ELEV_TS, EXIST=THERE)                               

      IF(THERE) THEN 
C       Assign a unit number and attempt to open the file.
        CPNT_UNIT = GET_UNIT(STDOUT)
        OPEN(CPNT_UNIT, FILE=CPNT_ELEV_TS, FORM = 'UNFORMATTED',
     A       STATUS = 'OLD')
        READ(CPNT_UNIT) ICODE
      ELSE                                                              
        WRITE(STDOUT,50) CPNT_ELEV_TS
        STOP 'Abnormal stop: errors found.' 
      ENDIF                                                             

      INQUIRE(FILE=DNS_ELEV_TS, EXIST=THERE)                               

      IF(THERE) THEN 
C       Assign a unit number and attempt to open the file.
        DNS_UNIT = GET_UNIT(STDOUT)
        OPEN(DNS_UNIT, FILE=DNS_ELEV_TS, FORM = 'UNFORMATTED',
     A       STATUS = 'OLD')
        READ(DNS_UNIT) ICODE
      ELSE                                                              
        WRITE(STDOUT,51) DNS_ELEV_TS
        STOP 'Abnormal stop: errors found.'
      ENDIF                                                             

      INQUIRE(FILE=GATE_FLOW_TS, EXIST=THERE)                               

      IF(THERE) THEN 
C       Assign a unit number and attempt to open the file.
        GATE_UNIT = GET_UNIT(STDOUT)
        OPEN(GATE_UNIT, FILE=GATE_FLOW_TS, FORM = 'UNFORMATTED',
     A       STATUS = 'OLD')
        READ(GATE_UNIT) ICODE
      ELSE                                                              
        WRITE(STDOUT,52) GATE_FLOW_TS
        STOP 'Abnormal stop: errors found.'
      ENDIF                                                             


      UPS_GIVEN = 0
      IF(UPS_ELEV_TS.NE.' ') THEN
        UPS_GIVEN = 1
C       User did give an upstream location distinct from 
C       the control point location.  Therefore, the upstream
C       location is not the same as the control point.

        INQUIRE(FILE=UPS_ELEV_TS, EXIST=THERE)                               

        IF(THERE) THEN 
C         Assign a unit number and attempt to open the file.
          UPS_UNIT = GET_UNIT(STDOUT)
          OPEN(UPS_UNIT, FILE=UPS_ELEV_TS, FORM = 'UNFORMATTED',
     A         STATUS = 'OLD')
        READ(UPS_UNIT) ICODE
        ELSE                                                              
          WRITE(STDOUT,53) UPS_ELEV_TS
          STOP 'Abnormal stop: errors found.'
        ENDIF                                                             
      ENDIF

C     Check that the function tables are known.
      CALL CHKTAB
     I           (15, STDOUT, FTPNT, PMXTAB,
     M            UD_TABLE,
     O            EFLAG)
      IF(EFLAG.EQ.0) THEN
C       Check the contents of the table to make sure that
C       the tables referenced in the table also exist and
C       are of the proper type.
        CALL CKTY15
     I             (UD_TABLE, STDOUT,
     O              EFLAG)
      ENDIF
      CALL CHKTAB
     I           (15, STDOUT, FTPNT, PMXTAB,
     M            DU_TABLE,
     O            EFLAG)
      IF(EFLAG.EQ.0) THEN
C       Check the contents of the table to make sure that
C       the tables referenced in the table also exist and
C       are of the proper type.
        CALL CKTY15
     I             (DU_TABLE, STDOUT,
     O              EFLAG)
      ENDIF

C     Find the datum for the type 15 tables.  Type 15 tables do 
C     do not store their datum; that value is present in the
C     type 13 tables within the type 15 table. 

       HDATUM = GETHDD(GET_TYPE13_ADDRESS(UD_TABLE))


C     Get the maximum gate opening.
      HG_MAX = GET_MAX_GATE(UD_TABLE)
            
      WRITE(STDOUT,60) HDATUM, HG_MAX

      IF(OUTPUT_LEVEL.EQ.'MAX')  WRITE(STDOUT,58)

      KNT = 0
200   CONTINUE
C       We read each time step from the files and then invert the
C       gate relationship for each time step and tabulate the
C       results. 

        READ(CPNT_UNIT) CPNT_JT, CPNT_ELEV
        IF(CPNT_JT.LT.1.D-10) GOTO 9000
        MJD = INT(CPNT_JT)
        FRAC = CPNT_JT - DBLE(MJD)
        HR = 24.D0*FRAC
        CALL INVMJD
     I             (MJD,
     O              YR, MN, DY)
        IF(UPS_GIVEN.EQ.0) THEN
          UPS_ELEV = CPNT_ELEV
          UPS_JT = CPNT_JT
        ELSE
          READ(UPS_UNIT) UPS_JT, UPS_ELEV
        ENDIF
        READ(DNS_UNIT) DNS_JT, DNS_ELEV
        READ(GATE_UNIT) GATE_JT, GATE_FLOW

        IF(REVERSE_FLOW.EQ.'NO ') THEN
          IF(GATE_FLOW.LE.MIN_FLOW) GATE_FLOW = 0.0
        ELSE
          IF(ABS(GATE_FLOW).LE.MIN_FLOW) GATE_FLOW = 0.0
        ENDIF
C       Check for synchronization
        IF(ABS(CPNT_JT - UPS_JT).GT.EPS.OR.
     A     ABS(CPNT_JT - DNS_JT).GT.EPS.OR.
     B     ABS(CPNT_JT - GATE_JT).GT.EPS) THEN
          WRITE(STDOUT,54) CPNT_JT, UPS_JT, DNS_JT, GATE_JT
          STOP ' Abnormal stop.  Errors found.' 
        ENDIF

        HU = UPS_ELEV - FLOOD_ELEV
        HD = DNS_ELEV - FLOOD_ELEV

        IF(GATE_FLOW.NE.0.0) THEN
C         Solve for the gate setting that would have to exist to match the
C         given flow and water surface elevations.


          CALL   FIND_GATE_OPENING(
     I           STDOUT, UD_TABLE, DU_TABLE, UPS_ELEV, DNS_ELEV, 
     I           GATE_FLOW, HG_MAX, HDATUM,
     O           HG, P, FTYPE, RESULT_FLAG, FTYPE_CODE)



          IF(RESULT_FLAG.EQ.1) THEN
            KNT = KNT + 1
            IF(OUTPUT_LEVEL.EQ.'MAX') THEN
              WRITE(STDOUT,56) KNT, YR, MN, DY, HR, HU, HD, GATE_FLOW, 
     A                       HG, P, FTYPE, RESULT_FLAG
            ENDIF
          ENDIF

        ELSE
          HG = 0.0
          P = 0.0
          RESULT_FLAG = 0

        ENDIF

        IF(RESULT_FLAG.EQ.1) THEN
C         Increment the various matrices.  First find the row and 
C         column of the cell that is to be incremented. 

          HDIFF = HU - HD

          IF(HU.LT.ROW_BDYS(1).OR.HU.GT.ROW_BDYS(N_ROW_BDYS)) THEN
            WRITE(STDOUT,62) HU, ROW_BDYS(1), ROW_BDYS(N_ROW_BDYS)
            EFLAG = 1
            HU = ROW_BDYS(2)
          ENDIF
          IF(HDIFF.LT.COL_BDYS(1).OR.HDIFF.GT.COL_BDYS(N_COL_BDYS)) THEN
            WRITE(STDOUT,64) HDIFF, COL_BDYS(1), COL_BDYS(N_COL_BDYS)
            EFLAG = 1
            HDIFF = COL_BDYS(2)
          ENDIF 
          DO 130 I=1,N_ROW_BDYS-1
            IF(ROW_BDYS(I+1).GE.HU) THEN
C             I gives the row 
              GOTO 132
            ENDIF
130       CONTINUE
          WRITE(STDOUT,*) ' Bug in INV_GATE. Should not get here. 130'
          STOP ' Abnormal stop.  Bug found' 
132       CONTINUE
          DO 135 J=1,N_COL_BDYS-1
            IF(COL_BDYS(J+1).GE.HDIFF) THEN
C             J gives the column
              N_PER_CELL(I,J)  =  N_PER_CELL(I,J) + 1
              SUM_PER_CELL(I,J) = SUM_PER_CELL(I,J) + P
              MAX_PER_CELL(I,J) = MAX(P, MAX_PER_CELL(I,J))
              MIN_PER_CELL(I,J) = MIN(P, MIN_PER_CELL(I,J))
              GOTO 137
            ENDIF
135       CONTINUE
          WRITE(STDOUT,*) ' Bug in INV_GATE.  Should not get here. 135'
          STOP ' Abnormal stop.  Bug found.' 
137       CONTINUE
        ENDIF  
        GOTO 200
9000  CONTINUE

C     Close open files
      CALL FREE_UNIT(STDOUT, CPNT_UNIT)
      CALL FREE_UNIT(STDOUT, DNS_UNIT)
      CALL FREE_UNIT(STDOUT, GATE_UNIT)
      IF(UPS_GIVEN.NE.0) THEN
        CALL FREE_UNIT(STDOUT, UPS_UNIT)
      ENDIF

C     Compute the mean gate opening in each cell.
      
      DO 300 I=1,N_ROW_BDYS-1
        DO 290 J=1,N_COL_BDYS
          IF(N_PER_CELL(I,J).GT.0) THEN
            SUM_PER_CELL(I,J) = SUM_PER_CELL(I,J)/N_PER_CELL(I,J)
            WRITE(PC_PER_CELL(I,J)(1:6),'(F6.3)') SUM_PER_CELL(I,J)
            WRITE(MINC_PER_CELL(I,J),'(F6.2)') MIN_PER_CELL(I,J)
            WRITE(MAXC_PER_CELL(I,J),'(F6.2)') MAX_PER_CELL(I,J)
            WRITE(NC_PER_CELL(I,J),'(I6)') N_PER_CELL(I,J)
          ELSE
            PC_PER_CELL(I,J)(1:6) = '   0.0'
            MAXC_PER_CELL(I,J) = '   0.0'
            MINC_PER_CELL(I,J) = '   0.0'
            NC_PER_CELL(I,J) = '     0'
          ENDIF
290     CONTINUE
300   CONTINUE

      WRITE(STDOUT,66)
      DO 320 I=1,N_ROW_BDYS-1
        WRITE(STDOUT,68) (N_PER_CELL(I,J), J=1,N_COL_BDYS-1)
320   CONTINUE
C      WRITE(STDOUT,70) 
C      DO 330 I=1,N_ROW_BDYS-1
C        WRITE(STDOUT,72) (SUM_PER_CELL(I,J), J=1,N_COL_BDYS-1)
C330   CONTINUE
      WRITE(STDOUT,74)
      WRITE(STDOUT,75) (COL_BDYS(J),J=1,N_COL_BDYS-1)
      WRITE(STDOUT,76) (COL_BDYS(J),J=2,N_COL_BDYS)
      DO 340 I=1,N_ROW_BDYS-1
        WRITE(STDOUT,72) ROW_BDYS(I), 'AvrP:', 
     A             (PC_PER_CELL(I,J), J=1,N_COL_BDYS-1)
        WRITE(STDOUT,72) ROW_BDYS(I+1), 'MinP:', 
     B                  (MINC_PER_CELL(I,J), J=1,N_COL_BDYS-1)
        WRITE(STDOUT,73) 'MaxP:', (MAXC_PER_CELL(I,J), J=1,N_COL_BDYS-1)
        WRITE(STDOUT,73) 'N   :', (NC_PER_CELL(I,J), J=1,N_COL_BDYS-1)
        WRITE(STDOUT,*) ' '
340   CONTINUE
      
C     Create the items needed to output the table of type 10
      DO 350 I=1,N_ROW_BDYS-1
        TP = 0.5*(ROW_BDYS(I) + ROW_BDYS(I+1))
        WRITE(ROW_ARGS(I)(1:6),'(F6.2)') TP
350   CONTINUE
      DO 360 J=1,N_COL_BDYS - 1
        TP = 0.5*(COL_BDYS(J) + COL_BDYS(J+1))
        WRITE(COL_ARGS(J)(1:6),'(F6.2)') TP
360   CONTINUE
      
C     Set the formats for the table.
C     Note: The number of columns in the table is always one more than
C     the number of column arguments because there is a column of the 
C     row arguments.  Therefor the number of column boundaries, being 
C     one more than the number of column arguments, also gives the 
C     number of columns in the table. 
      WRITE(CHR3,'(I3)') N_COL_BDYS
      WRITE(CHR3_1,'(I3)') N_COL_BDYS - 1
      HDIN = '''('//CHR3//'A6)'''
      HDOUT = '''(1X,'//CHR3//'A6)'''
      HDITEM = '''(F6.0)'''
      ROWIN =  '''('//CHR3//'F6.0)'''
      ROWOUT = '''(1X,F6.1,'//CHR3_1//'F6.2)'''
      ROWLAB = 'Fstage'
      lABEL = ' Replace with desired value'

C      WRITE(STDOUT,*) ' HDIN=',HDIN
C      WRITE(STDOUT,*) ' HDOUT=',HDOUT
C      WRITE(STDOUT,*) ' HDITEM=',HDITEM
C      WRITE(STDOUT,*) ' ROWIN=',ROWIN
C      WRITE(STDOUT,*) ' ROWOUT=',ROWOUT
      CALL TYPE10_OUT(
     I                STDOUT, STDTAB, CONTROL_TAB, N_ROW_BDYS-1, 
     I                N_COL_BDYS-1, ROW_ARGS, FLOOD_ELEV, COL_ARGS,
     I                PC_PER_CELL, HDIN, HDITEM, HDOUT, ROWIN, ROWLAB,
     I                ROWOUT, LABEL)





      RETURN
      END
