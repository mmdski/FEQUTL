C
C
C
      SUBROUTINE   SCNPRO
     I                   (IN, STDOUT, MXPNTU, STADIR,
     M                    STATTB, SBU, LSNU, LNFLAG, XFLAG, OLDID,
     M                    FIRST, STATL, STATR,
     O                    NPNTU, NSUBU, XU, ZU, NU, NVARU, NATYU, YATNU,
     O                    NNYU, EFLAG, LEFT, RIGHT, XSNAME, NCON,
     O                    XFLINE)
 
C     + + + PURPOSE + + +
C     Scan WSPRO input in the file, IN, and get the next
C     cross section and return the values needed for computing
C     a cross section table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FIRST, IN, LNFLAG, MXPNTU, NCON, NPNTU, NSUBU,
     A        STDOUT, XFLAG
      INTEGER NNYU(20), NVARU(20), SBU(MXPNTU)
      REAL LEFT, LSNU(MXPNTU), NATYU(9,20), NU(20), RIGHT, STADIR,
     A     STATL, STATR, STATTB, XFLINE(6), XU(MXPNTU), YATNU(9,20),
     B     ZU(MXPNTU)
      CHARACTER OLDID*2, XSNAME*5
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     IN     - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     MXPNTU - Maximum number of points on a cross section boundary
C     STADIR - Direction of stationing
C     STATTB - stationing for the axis
C     SBU    - Subsection numbers for the line segments
C     LSNU   - Line segment Manning's n value
C     LNFLAG - Flag to signal the state of the current line buffer
C     XFLAG  - cross section processing flag
C     OLDID  - identification value for the previous card
C     FIRST  - flag to detect the first header
C     STATL  - stationing for the flow path to the left
C     STATR  - stationing for the flow path to the right
C     NPNTU  - Number of points on boundary of a cross section
C     NSUBU  - Number of subsections
C     XU     - Offsets for points on boundary of cross section
C     ZU     - Elevation of points on boundary of cross section
C     NU     - Vector for Manning's n values
C     NVARU  - Flag for variation of Manning's n in each subsection
C     NATYU  - Mannings's n value at depth in YATNU
C     YATNU  - Depth values for the Manning's n values in NATYU
C     NNYU   - Number of values for Manning's n variation with depth
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
C     XSNAME - WSPRO name for the cross section
C     NCON   - flag for variation of Manning's n in the vertical
C     XFLINE - offsets assigned to flow lines
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'wsproxs.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER MXNVAL, NRID, HEADER, CNTRL, XSD, ICD, BDD, CVD
      REAL DEFALT, DGTORD
      PARAMETER(MXNVAL=40, NRID=35, HEADER=1, CNTRL=2, XSD=3,
     A     ICD=4, BDD=5, CVD=6, DEFALT=-1.E30,DGTORD=1.745329E-2)
 
C     + + + SAVED VALUES + + +
      INTEGER RIDCAT(NRID)
      REAL DSTAT
      CHARACTER LINE*80, RIDTAB(NRID)*2
      SAVE DSTAT, LINE, RIDCAT, RIDTAB
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CAT, I, IP, J, NPI, NUMVAL, OPT, OPTION
      INTEGER CLEN(MXNVAL), IVAL(MXNVAL), TERML(MXNVAL), 
     A        TERMCLS(MXNVAL), TYPE(MXNVAL)
      REAL AXISL, COSSK, DELY, FAC, LEFTL, RIGHTL, RVAL(MXNVAL),
     A     XARG(MXNSA), XE, XEPS, XINVAR, XS
      DOUBLE PRECISION DPVAL(MXNVAL)
      CHARACTER CVAL(MXNVAL)*256, HEADID*2, RECID*2, SECID*5, 
     A          TERM(MXNVAL)*1
 
C     + + + INTRINSICS + + +
      INTRINSIC COS, MOD
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL FDCNMN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BINSER, FDCNMN, GETVAL, inline, INSPT
 
C     + + + DATA INITIALIZATIONS + + +
      DATA
     *RIDTAB/'* ','AB','AS','BD','BL','BP',
     *       'BR','CC','CD','CG','CV',
     *       'ER','EX','FL','GR','GT',
     *       'HP','J1','J3','KD','N',
     *       'ND','PW','PX','Q','SA',
     *       'SD','SK','T1','T2','T3',
     *       'WS','XR','XS','XT'/,
     *RIDCAT/    2,    5,    1,    5,    5,    5,
     *           1,    6,    5,    6,    1,
     *           2,    2,    3,    3,    3,
     *           2,    2,    2,    5,    3,
     *           3,    5,    2,    4,    3,
     *           1,    4,    2,    2,    2,
     *           4,    1,    1,    1/
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:536* Record id=',A2,' unknown.')
 51   FORMAT(' *ERR:541* X and Y value count=',I5,' is odd. Both ',
     A    'values ',/,11X,'must be on the same input record.')
 52   FORMAT(' *ERR:542* Section id= ',A5,' has no station value.')
 53   FORMAT(' *ERR:521* Number of subareas=',I5,' inconsistent with',
     A    ' number of n values=',I5)
 54   FORMAT(' *ERR:522* Number of n values=',I5,' inconsistent with',
     A       ' number of breakpoints on ND.')
 56   FORMAT(' *ERR:523* No GR cards found following XT card.')
 58   FORMAT(' *ERR:524* There are',I3,' items for flow-line data.',
     A        ' Only 1, 3, or 5 items are valid.')
C***********************************************************************
C     GET THE NEXT LINE FROM THE WSPRO INPUT FILE
 
 90   CONTINUE
        IF(LNFLAG.EQ.0) THEN
          CALL inline
     I              (IN, STDOUT,
     O               LINE)
        ELSE
C         CLEAR THE LINE BUFFER FLAG TO SIGNAL THAT THE BUFFER WILL
C         HAVE BEEN PROCESSED UPON EXIT FROM THIS ROUTINE.
          LNFLAG = 0
        ENDIF
 
C       Get the fixed-field part of the record with an internal
C       read from a character string.   The value of the
C       record id and the state of various flags determines the
C       action to take.
 
 
        READ(LINE(1:10),'(A2,1X,I1,1X,A5)') RECID, OPTION, SECID
 
C        WRITE(STDOUT,*)' RECID=',RECID,' SECID=',SECID
 
C       If the record id is a header record, it means that the
C       cross section we are processing is complete, if there is
C       a cross section.  XFLAG is the indicator for cross section
C       processing.  XFLAG=0 means no cross section is
C       underway,  XFLAG=1 means there is a cross section in
C       process.
 
C       Find the category of the record id. HEADER-header record,
C       CNTRL-control record, XSD-cross section data, ICD-initial
C       condition data, BDD-bridge data, and CVD-culvert data.
 
        CALL BINSER
     I             (RECID, NRID, RIDTAB,
     O              IP)
        IF(IP.EQ.0) THEN
          EFLAG = 1
          WRITE(STDOUT,50) RECID
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        IF(RECID.EQ.'ER') THEN
          IF(XFLAG.EQ.0) THEN
C           All information in the file has been processed.
            XFLAG = -1
            RETURN
          ENDIF
        ENDIF
        CAT = RIDCAT(IP)
 
C        WRITE(STDOUT,*) ' FOR RECID=',RECID,' CAT=',CAT,' XFLAG=',XFLAG
 
        IF(CAT.NE.XSD) THEN
          IF(XFLAG.EQ.1) THEN
C           We now have the description of a cross section in
C           hand.  Convert it to FEQXEXT form and return the values.
C           Remember that there are good values in LINE waiting to be
C           processed.
            LNFLAG = 1
 
C           Request that the whole cross section be processed.
 
            LEFT = 0.0
            RIGHT = 0.0
C           Determine the description for the current cross section.
C           There are many choices to make here.  The description
C           could have been input with the current header, part or
C           all of the description could come from the last description
C           read for a previous header, or part of the description
C           could come from a template section that was stored
C           earlier.  Ground points must be input with a template
C           header and SA points may be input with a template header.
C           If the header record is XT then this section stores the
C           current ground and if present, SA information.
C           No cross section data should be returned.  Processing
C           of the WSPRO input should continue.
 
            IF(HEADID.NE.'XT') THEN
C             Compute the skew factor
              COSSK = COS(SKEW*DGTORD)
 
              IF(GTFLAG.EQ.1) THEN
C               We have a current GT record. Do not allow GT values to
C               propagate.
 
                GTFLAG = 0
 
C               Get the limits to use.  Note that the retained cross
C               section does not retain the truncation limits.
 
                LEFT = XLIML
                RIGHT = XLIMR
 
C               It appears that the definition of ground points using
C               GT is treated as if it were input of GR cards.  Thus
C               a subsequent cross section heading without a GT or GR
C               takes its ground points from the most recent previous
C               occurrence of GT or GR.
 
                GRSRD = SRD
 
                IF(YSHIFT.EQ.DEFALT) THEN
C                 No vertical shift given.  Use the value of
C                 VSLOPE and the difference in stations for any
C                 vertical adjustment.
 
                  DELY = (SRD - XTSRD)*VSLOPE
                ELSE
                  DELY = YSHIFT
                ENDIF
                NGP = XTNGP
                DO 210 I=1,NGP
                  Y(I) = XTYGR(I) + DELY
                  YGR(I) = Y(I)
 210            CONTINUE
 
C               Do the horizontal adjustments.  Find the invariant
C               point.  Combine with skew adjustment.
                FAC = COSSK*SCALE
                IF(XORIG.GT.DEFALT) THEN
                  XINVAR = XORIG
                ELSE
                  XINVAR = 0.0
                ENDIF
                DO 220 I=1,NGP
                  X(I) = FAC*(XTXGR(I) - XINVAR)  + XINVAR
                  XGR(I) = X(I)
 220            CONTINUE
 
C               Now find the subarea limits for the ground points
C               taken from the template section.
                IF(SAFLAG.EQ.1) THEN
C                 The user has input the SA under the same header as
C                 the GT record.
                  NSA = NSAT + 1
                  DO 225 I=1,NSAT
                    XSA(I) = COSSK*XSAT(I)
                    XSAT(I) = XSA(I)
 225              CONTINUE
                ELSE
                  IF(XTNSA.GT.0) THEN
C                   User has ommitted the SA record but it exists
C                   in the template record.
                    NSA = XTNSA + 1
                    NSAT = XTNSA
                    DO 230 I=1,XTNSA
                      XSA(I) = FAC*(XTXSA(I) - XINVAR) + XINVAR
                      XSAT(I) = XSA(I)
 230                CONTINUE
                  ELSE
C                   User has ommitted the SA record and it does not
C                   exist in the template record.  Use the most recent
C                   value available.
 
                    NSA = NSAT + 1
                    DO 235 I=1,NSAT
                      XSA(I) = COSSK*XSAT(I)
                      XSAT(I) = XSA(I)
 235                CONTINUE
                  ENDIF
                ENDIF
              ELSE
C               Source of ground points is not the template section.
C               Take from the most recent value of ground points
C               and adjust for any non-zero valley slope.  The adjustment
C               affects the retained ground point values just as if
C               the ground points had been input directly.
 
C                WRITE(STDOUT,*) ' SRD=',SRD,' GRSRD=',GRSRD,
C     A                   ' VSLOPE=',VSLOPE
 
                DELY = (SRD - GRSRD)*VSLOPE
                WRITE(STDOUT,*) ' DELY=',DELY
 
                DO 240 I=1,NGP
                  Y(I) = YGR(I) + DELY
                  X(I) = COSSK*XGR(I)
                  YGR(I) = Y(I)
                  XGR(I) = X(I)
 240            CONTINUE
 
                NSA = NSAT + 1
                DO 245 I=1,NSAT
                  XSA(I) = COSSK*XSAT(I)
                  XSAT(I) = XSA(I)
 245            CONTINUE
              ENDIF
 
C             Update the station of the most recent ground points.
              GRSRD = SRD
 
C             Now the ground points and the subarea limits are
C             defined in all cases.  Decipher the meaning of the
C             the N  and ND records.  If the number of
C             values on the N record is the same as NSA, it means that
C             the values on the N record represent values of n that do
C             not vary with depth.  If the number of values on the N
C             record is 2*NSA, then they represent the break-point values
C             of Manning's n.  These are only valid if ND has set the
C             breakpoints.
              IF(NN.EQ.NSA) THEN
C               Constant values of Manning's n.  Values are already
C               in NVAL(*).
              ELSEIF(NN.EQ.NSA+NSA) THEN
C               Breakpoint variation of n in the vertical
                DO 200 I=1,NSA
                  J = I + I
                  TOPN(I) = NVAL(J)
                  BOTN(I) = NVAL(J-1)
 200            CONTINUE
C               Transfer the breakpoint depths.
                IF(NN.NE.NND) THEN
                  WRITE(STDOUT,54) NN, NND
                  EFLAG = 1
                ELSE
                  DO 205 I=1,NSA
                    J = I +I
                    TOPD(I) = NNDBRK(J)
                    BOTD(I) = NNDBRK(J-1)
 205              CONTINUE
                ENDIF
              ELSE
                WRITE(STDOUT,53) NSA, NN
                EFLAG = 1
              ENDIF
 
C             Now decipher any flow length data.   Friction slope
C             averaging is ignored.
              XS = XGR(1)
              XE = XGR(NGP)
              XEPS = 0.001*(XE - XS)
              IF(XEPS.GT.1.0) XEPS = 1.0
              IF(FLFLAG.EQ.1.AND.NFL.GT.0) THEN
C               There is flow length data present.
                IF(NFL.EQ.1) THEN
                  FLEN(1) = FLDAT(1)
C                 There is only one flow length.  Therefore it is the
C                 lenght for all points.
                  LEFTL = FLEN(1)
                  AXISL = FLEN(1)
                  RIGHTL = FLEN(1)
C                 Set the offsets for sinuousity computation.  There
C                 are always 6 offsets because there will always be
C                 6 flow lines.
                  XFLINE(1) = XS - 2.*XEPS
                  XFLINE(2) = XS - XEPS
                  XFLINE(3) = XS
                  XFLINE(4) = XE
                  XFLINE(5) = XE + XEPS
                  XFLINE(6) = XE + 2.*XEPS
                ELSEIF(NFL.EQ.3) THEN
                  FLEN(1) = FLDAT(1)
                  XFL(1) = COSSK*FLDAT(2)
                  YFL(1) = FDCNMN(NGP, X, Y, DEFALT, XFL(1))
                  FLEN(2) = FLDAT(3)
                  YFL(2) = FDCNMN(NGP, X, Y, XFL(1), -DEFALT)
                  NFL = 2
                  IF(YFL(1).LT.YFL(2)) THEN
C                   Assume axis is for first flow length.
                    LEFTL = FLEN(1)
                    AXISL = FLEN(1)
                    RIGHTL = FLEN(2)
                    XFLINE(1) = XS - 2.*XEPS
                    XFLINE(2) = XS - XEPS
                    XFLINE(3) = XS
                    XFLINE(4) = XFL(1)
                    XFLINE(5) = XFL(1) + XEPS
                    XFLINE(6) = XE
                  ELSE
                    LEFTL = FLEN(1)
                    AXISL = FLEN(2)
                    RIGHTL = FLEN(2)
                    XFLINE(1) = XS
                    XFLINE(2) = XFL(1) - XEPS
                    XFLINE(3) = XFL(1)
                    XFLINE(4) = XE
                    XFLINE(5) = XE + XEPS
                    XFLINE(6) = XE + 2.*XEPS
                  ENDIF
                ELSEIF(NFL.EQ.5) THEN
                  FLEN(1) = FLDAT(1)
                  XFL(1) = COSSK*FLDAT(2)
                  YFL(1) = FDCNMN(NGP, X, Y, DEFALT, XFL(1))
                  FLEN(2) = FLDAT(3)
                  XFL(2) = COSSK*FLDAT(4)
                  YFL(2) = FDCNMN(NGP, X, Y, XFL(1), XFL(2))
                  FLEN(3) = FLDAT(5)
                  YFL(3) = FDCNMN(NGP, X, Y, XFL(2), -DEFALT)
                  NFL = 3
                  LEFTL = FLEN(1)
                  AXISL = FLEN(2)
                  RIGHTL = FLEN(3)
                  XFLINE(1) = XS
                  XFLINE(2) = XFL(1) - XEPS
                  XFLINE(3) = XFL(1)
                  XFLINE(4) = XFL(2)
                  XFLINE(5) = XFL(2) + XEPS
                  XFLINE(6) = XE
                ELSE
                  WRITE(STDOUT,58) NFL
                  EFLAG = 1
                ENDIF
                IF(STADIR.GE.0.0) THEN
                  STATL = STATL + LEFTL
                  STATTB = STATTB + AXISL
                  STATR = STATR + RIGHTL
                ELSE
                  STATL = STATL - LEFTL
                  STATTB = STATTB - AXISL
                  STATR = STATR - RIGHTL
                ENDIF
              ELSE
C               No flow line lengths given.  Thus SRD must give the
C               stations for all flow lengths.
                STATL = SRD + DSTAT
                STATTB = SRD + DSTAT
                STATR = SRD + DSTAT
                XFLINE(1) = XS - 2.*XEPS
                XFLINE(2) = XS - XEPS
                XFLINE(3) = XS
                XFLINE(4) = XE
                XFLINE(5) = XE + XEPS
                XFLINE(6) = XE + 2.*XEPS
              ENDIF
 
C             WSPRO cross section is defined here.  Convert the
C             description to FEQXEXT form and return.
 
C             Set the number of sub areas(subsections)
              NSUBU = NSA
 
C             Establish the variation of n for each subsection.
              IF(NN.EQ.NSA) THEN
C               n is constant in each subsection
                DO 300 I=1,NSA
                  NVARU(I) = 0
                  NU(I) = NVAL(I)
 300            CONTINUE
                NCON = 1
              ELSE
                NCON = 0
C               n varies with hydraulic depth in each subsection
                DO 305 I=1,NSA
                  NU(I) = NVAL(2*I-1)
                  NVARU(I) = 1
                  NNYU(I) = 3
                  YATNU(1,I) = 0.0
                  NATYU(1,I) = BOTN(I)
                  YATNU(2,I) = BOTD(I)
                  NATYU(2,I) = BOTN(I)
                  YATNU(3,I) = TOPD(I)
                  NATYU(3,I) = TOPN(I)
 305            CONTINUE
              ENDIF
C             Transfer the ground points
              DO 307 I=1,NGP
                XU(I) = X(I)
                ZU(I) = Y(I)
 307          CONTINUE
              NPNTU = NGP
 
C             Add the points in XSA to the ground points if needed.
C             Set the sinuousity variation to 1 and leave off the
C             the first and last values in XARG to be consistent
C             with other usage.  Points are not added if they
C             exactly match a point already in the boundary
C             specification.
 
              DO 310 I=1,NSA - 1
                XARG(I+1) = XSA(I)
 310          CONTINUE
              NPI = NSA + 1
              CALL INSPT
     I                  (STDOUT, NPI, XARG, 1,
     M                   NPNTU, XU, ZU, SBU, LSNU,
     O                   EFLAG)
 
C             Assign the subsections to the ground points.  Add a
C             large offset to the end of XSA.  If no SA record
C             ever appears in the input NSAT is 0 and NSA becomes 1.
 
              XSA(NSA) = -DEFALT
              J = 1
              DO 315 I=1,NPNTU-1
                IF(XU(I).LT.XSA(J)) THEN
                  SBU(I) = J
                ELSE
                  J = J + 1
                  SBU(I) = J
                ENDIF
 315          CONTINUE
              SBU(NPNTU) = -1
 
C             Now, mark some cross section data flags as being
C             old input.
              IF(NDFLAG.EQ.1) NDFLAG = 2
              IF(NFLAG.EQ.1) NFLAG = 2
              IF(SAFLAG.EQ.1) SAFLAG = 2
              IF(GRFLAG.EQ.1) GRFLAG = 2
              FLFLAG = 0
              NFL = 0
              RETURN
            ELSE
C             XT header here.  Transfer the user input to the template
C             section for later use.
              IF(GRFLAG.EQ.1) THEN
C               Ground points have been input with the XT header
C               record.   Transfer them to the XT storage.
                XTNGP = NGP
                DO 207 I=1,NGP
                  XTXGR(I) = XGR(I)
                  XTYGR(I) = YGR(I)
 207            CONTINUE
              ELSE
                WRITE(STDOUT,56)
                EFLAG = 1
              ENDIF
              IF(SAFLAG.EQ.1) THEN
C               SA record has been input with XT.
                XTNSA = NSAT
                DO 208 I=1,NSAT
                  XTXSA(I) = XSAT(I)
 208            CONTINUE
              ELSE
                XTNSA = 0
              ENDIF
C             Signal a fresh start for cross sections.
              XFLAG = 0
            ENDIF
          ELSEIF(CAT.EQ.HEADER) THEN
C           Remember the name
            XSNAME = SECID
C           Get the values from the remainder of the record.
            OPT = 0
            DO 101 I=1,MXNVAL
              TYPE(I) = 2
C             Set the default value.
              RVAL(I) = DEFALT
 101        CONTINUE
 
C            CALL GETVAL
C     I                 (STDOUT, LINE(11:80), MXNVAL, OPT,
C     O                  TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, EFLAG,
C     O                  NUMVAL)
             CALL GETVAL
     I                  (STDOUT, LINE(11:80), MXNVAL, OPT,
     O                   TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O                   EFLAG, TERM, TERML, TERMCLS,
     O                   NUMVAL)
 
            IF(NUMVAL.LT.1) THEN
              WRITE(STDOUT,52) SECID
              EFLAG = 1
            ENDIF
 
            HEADID = RECID
            SRD = RVAL(1)
            IF(FIRST.EQ.1) THEN
C             Compute the change in station required.
              DSTAT = STATTB - SRD
              FIRST = 0
            ENDIF
C           Process the header records of interest
            IF(RECID.EQ.'XS'.OR.RECID.EQ.'AS'.OR.RECID.EQ.'BR') THEN
              IF(RVAL(2).GT.DEFALT) THEN
                SKEW = RVAL(2)
              ELSE
                SKEW = 0.0
              ENDIF
              IF(RVAL(5).GT.DEFALT.AND.RECID.NE.'BR') THEN
                VSLOPE = RVAL(5)
              ENDIF
              XFLAG = 1
            ELSEIF(RECID.EQ.'XT') THEN
C             Remember the station for the template.
              XTSRD = SRD
              IF(RVAL(2).GT.DEFALT) THEN
                VSLOPE = RVAL(2)
              ENDIF
              XFLAG = 1
            ELSEIF(RECID.EQ.'XR') THEN
              IF(RVAL(5).GT.DEFALT) THEN
                SKEW = RVAL(5)
              ELSE
                SKEW = 0.0
              ENDIF
              XFLAG = 1
            ENDIF
          ENDIF
        ELSE
          IF(XFLAG.EQ.1) THEN
C           Cross section data records here. Use subroutine GETVAL to
C           get the free-form portion of the record.
C           MXNVAL gives the maximun number of values that can be
C           returned.  GETOPT=1 has GETVAL return the variable
C           type found in the vector, TYPE.  NUMVAL gives the number of
C           values found including those given by dual commas or by
C           the asterisk place holder.  We set a standard value, not
C           expected in the input, into each of the four vectors,
C           IVAL-integer value, RVAL-real value, DPVAL-double precision
C           value, and CVAL-an identifier value.  Here we treat all values
C           as real.   GETVAL will complain if any item in the remainder
C           of LINE is not a number.
 
            OPT = 0
            DO 100 I=1,MXNVAL
              TYPE(I) = 2
C             Set the default value.
              RVAL(I) = DEFALT
 100        CONTINUE
 
C            CALL GETVAL
C     I                 (STDOUT, LINE(11:80), MXNVAL, OPT,
C     O                  TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, EFLAG,
C     O                  NUMVAL)
             CALL GETVAL
     I                  (STDOUT, LINE(11:80), MXNVAL, OPT,
     O                   TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O                   EFLAG, TERM, TERML, TERMCLS,
     O                   NUMVAL)
 
            IF(RECID.EQ.'GR') THEN
C             Number of values should be even.
              IF(MOD(NUMVAL,2).NE.0) THEN
                WRITE(STDOUT,51) NUMVAL
                EFLAG = 1
              ENDIF
 
C             Ground point data.  Is this the first card or is
C             it a continuation card?
              IF(OLDID.NE.'GR') THEN
C               This is the first card.  Set the ground point counter
C               and set the ground point data flag.
                NGP = NUMVAL/2
                GRFLAG = 1
C               Remember the station of the GR data for use with VSLOPE
                GRSRD = SRD
 
                DO 105 I=1,NGP
                  J = I + I
                  YGR(I) = RVAL(J)
                  XGR(I) = RVAL(J-1)
 105            CONTINUE
              ELSE
C               This is a continuation record.  Add
C               to the ground point information.
 
                DO 110 I=1,NUMVAL/2
                  NGP = NGP + 1
                  J = I + I
                  YGR(NGP) = RVAL(J)
                  XGR(NGP) = RVAL(J-1)
 110            CONTINUE
              ENDIF
            ELSEIF(RECID.EQ.'N ') THEN
C             Manning's n record.
              IF(OLDID.NE.'N ') THEN
C               This is the first record.
                NN = NUMVAL
                NFLAG = 1
C               Not yet known what the values on the N record mean.
C               Place them in NVAL(*) for later processing.
                DO 115 I=1,NUMVAL
                  NVAL(I) = RVAL(I)
 115            CONTINUE
              ELSE
C               This is a continuation record.
                DO 120  I=1,NUMVAL
                  NN = NN + 1
                  NVAL(NN) = RVAL(I)
 120            CONTINUE
              ENDIF
            ELSEIF(RECID.EQ.'ND') THEN
C             Depth breakpoints for vertical variation of n
              IF(OLDID.NE.'ND') THEN
C               This is the first record.
                NND = NUMVAL
                NDFLAG = 1
C               Transfer and process later when all have been input.
                DO 125 I=1,NUMVAL
                  NNDBRK(I) = RVAL(I)
 125            CONTINUE
              ELSE
C               This is a continuation record.
                DO 130 I=1,NUMVAL
                  NND = NND + 1
                  NNDBRK(NND) = RVAL(I)
 130            CONTINUE
              ENDIF
            ELSEIF(RECID.EQ.'GT') THEN
C             Adjustment values for converting previous ground points
C             to the current ground points.
              GTFLAG = 1
              YSHIFT = RVAL(1)
              IF(RVAL(2).GT.DEFALT) THEN
                XLIML = RVAL(2)
              ENDIF
              IF(RVAL(3).GT.DEFALT) THEN
                XLIMR = RVAL(3)
              ENDIF
              IF(RVAL(4).GT.DEFALT) THEN
                SCALE = RVAL(4)
              ELSE
                SCALE = 1.0
              ENDIF
              XORIG = RVAL(5)
            ELSEIF(RECID.EQ.'SA') THEN
C             Subarea limits.
              SAFLAG = 1
              IF(OLDID.NE.'SA') THEN
C               This is the first record.  Use NSAT as a counter and
C                 adjust later. Save in a temp location and make
C                 the selection of source later.
                NSAT = NUMVAL
                DO 135 I=1,NUMVAL
                  XSAT(I) = RVAL(I)
 135            CONTINUE
              ELSE
C               This is a continuation record.
                DO 140 I=1,NUMVAL
                  NSAT = NSAT + 1
                  XSAT(NSAT) = RVAL(I)
 140            CONTINUE
              ENDIF
            ELSEIF(RECID.EQ.'FL') THEN
C             Flow line information.  Store information and process
C             when all data for a cross section is in hand.
              FLFLAG = 1
              IF(OLDID.NE.'FL') THEN
C               First record.
                NFL = NUMVAL
                DO 145 I=1,NUMVAL
                  FLDAT(I) = RVAL(I)
 145            CONTINUE
              ELSE
C               This is a continuation record.
                DO 150 I=1,NUMVAL
                  NFL = NFL + 1
                  FLDAT(NFL) = RVAL(I)
 150            CONTINUE
              ENDIF
            ENDIF
          ENDIF
        ENDIF
 
        OLDID = RECID
 
        GOTO 90
 
      END
C
C
C
      SUBROUTINE   PROINT()
 
C     + + + PURPOSE + + +
C     Initialize the values in the WSPRO common block to
C     starting values.
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'wsproxs.cmn'
C***********************************************************************
      ASFLAG = 0
      FLFLAG = 0
      GRFLAG = 0
      GTFLAG = 0
      NFLAG = 0
      NDFLAG = 0
      SAFLAG = 0
      XRFLAG = 0
      XSFLAG = 0
      XTFLAG = 0
      NSAT = 0
      NSA = 0
      NFL = 0
      XTNSA = 0
      RETURN
      END
C
C
C
      SUBROUTINE   WPRO14
     I                   (STDIN, STDOUT, STDTAB,
     M                    TABDIR, EFLAG)
 
C     + + + PURPOSE + + +
C     Construct a table of type 14 from one or more .prt files from
C     WSPRO.  These files will have user specified output that will
C     contain the values needed and in the proper order.  The
C     order is established by the WSPROQZ command implemented in
C     subroutine WPROQZ.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER MAXN
      PARAMETER(MAXN=16)
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IT, KNT, MVAL, NFRAC, NZD, NZD2, NZDF, NZDL, OPT, TAB
      INTEGER CLEN(MAXN), IVAL(MAXN), TERML(MAXN), TERMCLS(MAXN),  
     A        TYPE(MAXN)
      REAL FROUDE, HDATUM, HDVEC(PMXNHU), HUMAT(PMXNHU,PMXFRC), MAXQ,
     A     PFQVEC(PMXFRC), Q(PMXNHU,PMXFRC), QFVEC(PMXNHU), RVAL(MAXN),
     B     WSEL, zrhufd
      real*8 DPVAL(MAXN), easting, northing

      CHARACTER CVAL(MAXN)*256, HEAD*80, LABEL*50, LINE*80, LINE1*81,
     A          NAME*64, TERM(MAXN)*1, XSID*5, 
     b     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8

      LOGICAL THERE
 
      CHARACTER GETTOK*5, TABID*16

C     + + + EXTERNAL NAMES + + +
      EXTERNAL GETTOK, GETVAL, inline, TABCHK, TWDOUT, os_file_style
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5)
 2    FORMAT(6X,A)
 3    FORMAT(A)
 24   FORMAT(A80)
 25   FORMAT(4X,5X,13X,F8.0,F8.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' Table id for 2-D table of type 14= ',A)
 51   FORMAT(/,' Label for 2-D table: ',A)
 54   FORMAT(/,' ',A80)
 56   FORMAT(' Current file name is:',A)
 58   FORMAT(/,' *ERR:519* Expected Q card but found:',A2)
 60   FORMAT(/,' *WRN:510* Froude number=',F7.2,' > 1.0 at XSID=',A5)
62    FORMAT(' Not enough user tables.  One or more profiles',
     A           ' failed in WSPRO.')
 96   FORMAT(/,' FILE NAMED:',A,' NOT FOUND.  CHECK SPELLING OF',
     A       ' FILE NAME.')
 98   FORMAT(/,' Found ',I3,' downstream elevations and ',I3,
     A      ' partial maximum flows.')
 99   FORMAT(/,' All files have been processed.')
C***********************************************************************
C     Get the table number to be used by the two-D table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE,'TAB',
     O                EFLAG, TABID, TAB)
C      READ(LINE,1,ERR=991) TABID
      WRITE(STDOUT,50) TABID
      CALL TABCHK
     I           (STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)

c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

 
C     Get the label for the 2-D table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) LABEL
      WRITE(STDOUT,51) LABEL
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,54) HEAD
 
C     Input each file name and process in order.
C     Each file must have a descriptive header that agrees exactly with
C     its contents.  This means that each file must represent a
C     complete series of flows for a given downstream head.  The downstream
C     heads must be presented in increasing order and the flows
C     must be presented in increasing order.   Zero flow is not included
C     in any of the files because the result for zero flow is trivial.
 
C     Initialize the global number of zero depths counter
      NZD = 1
      NZD2 = 1
 100  CONTINUE
        CALL inlineb
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,3,ERR=991) NAME
 
        IF(NAME.NE.' ') THEN
           CALL MAYBE_ADD_HOME(
     M                         name)

          call os_file_style(
     m                              name)
          INQUIRE(FILE=NAME, EXIST=THERE)
          IF(THERE) THEN
            WRITE(STDOUT,56) NAME
            OPEN(UNIT=STD48, FILE=NAME, STATUS='OLD')
          ELSE
            WRITE(STDOUT,96) NAME
            STOP 'Abnormal stop. Errors found.'
          ENDIF
 
C         Initialize the local number of zero depths counter.
          NZDL = 1
 200      CONTINUE
            READ(STD48,24,END=900) LINE
            IF(LINE(14:20).EQ.'WSPROZQ') THEN
C             Found the standard header line.  Read the number of flows,
C             number of downstream elevations, and the head datum.
              READ(LINE,'(40X,I5,5X,I5,8X,F10.0)') NFRAC, NZDF, HDATUM
C              WRITE(STDOUT,*) ' NFRAC=',NFRAC, ' NZDF=',NZDF,
C     A                ' HDATUM=',HDATUM
 
C             Skip over the two lines defining the user output.
              READ(STD48,24) LINE
              READ(STD48,24) LINE
 
C             Set GETVAL to get the number of values and return the
C             types.  Types must be self defining.  That is, a real
C             value must have a decimal point.
              OPT = 1
              KNT = 0
 300          CONTINUE
                READ(STD48,24) LINE1
C               Search for user quote
                IT = INDEX(LINE1,'''')
                IF(IT.EQ.0) THEN
                  LINE1(81:81) = ''''
                ENDIF
                  
                CALL GETVAL
     I                     (STDOUT, LINE1, MAXN, OPT,
     O                      TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O                      EFLAG, TERM, TERML, TERMCLS,
     O                      MVAL)
                DO 310 I=2,MVAL
                  READ(CVAL(I)(1:CLEN(I)),*) RVAL(I)
 310            CONTINUE
                IF(CVAL(1)(1:1).EQ.'Q') THEN
C                 Flow definition.
                  DO 400 I=2,MVAL
                    KNT = KNT + 1
                    Q(NZD,KNT) = RVAL(I)
                    IF(KNT.EQ.NFRAC) THEN
C                     Current downstream depth is complete.
                      NZD = NZD + 1
                      NZDL = NZDL + 1
                      KNT = 0
                      IF(NZDL.GT.NZDF) THEN
C                       Flows from the current file have been completed.
                        NZDL = 1
                        WRITE(STDOUT,*) ' Flows completed'
                        GOTO 401
                      ENDIF
                    ENDIF
 400              CONTINUE
              ELSE
                WRITE(STDOUT,58) CVAL(1)(1:2)
                EFLAG = 1
              ENDIF
              GOTO 300
 401          CONTINUE
 
C             Now look for the string heading the user output.
 
              KNT = 0
              NZDL = 1
 500          CONTINUE
                READ(STD48,24,END=910) LINE
                IF(LINE(4:13).EQ.'FIRST USER') THEN
C                  WRITE(STDOUT,*) ' FOUND: FIRST USER'
C                 Found the line.  Skip over the blank line and heading.
                  READ(STD48,24) LINE
                  READ(STD48,24) LINE
C                 Read the output lines, and check the Froude number.
C                 The downstream elevation is in the line with cross
C                 section id EXIT and the upstream elevation is in the
C                 line with cross section id APPRO.
 600              CONTINUE
                    READ(STD48,'(A)') LINE
                    XSID = LINE(5:9)
                    XSID = GETTOK(XSID)
                    IF(XSID.EQ.'EXIT'.OR.XSID.EQ.'APPRO') THEN
                      READ(LINE,25) FROUDE, WSEL
                    ELSE
                      FROUDE = 0.0
                      WSEL = 0.0
                    ENDIF
                    IF(FROUDE.GT.1.0) THEN
                      WRITE(STDOUT,60) FROUDE, XSID
                    ENDIF
                    IF(XSID.EQ.'EXIT') THEN
C                     Remember the downstream depth.  Note that this
C                     statement is executed NZDF times with the
C                     same water surface elevation.
 
                      HDVEC(NZD2) = WSEL - HDATUM
                    ENDIF
                    IF(XSID.EQ.'APPRO') THEN
C                     Get the upstream elevation and save it.
                      KNT = KNT + 1
C                     Leave space for zero fraction of free flow
C                     by adding one to KNT when storing the upstream
C                     heads.  Required by TWDOUT.
                      HUMAT(NZD2,KNT+1) = WSEL - HDATUM
                      IF(KNT.EQ.NFRAC) THEN
C                       The current downstream elevation's flows are
C                       complete.
                        NZD2 = NZD2 + 1
                        NZDL = NZDL + 1
                        KNT = 0
                        IF(NZDL.GT.NZDF) THEN
C                         Current file is processed.
                          NZDL = 1
                          GOTO 101
                        ENDIF
                      ENDIF
C                     APPRO is the last line of interest in the current
C                     user output.
                      GOTO 601
                    ENDIF
                    GOTO 600
 601              CONTINUE
              ENDIF
              GOTO 500
          ELSE
            GOTO 200
          ENDIF
 101      CONTINUE
          CLOSE(STD48)
          GOTO 100
        ENDIF
 
C     All the files have been processed.
      WRITE(STDOUT,99)
      NZD = NZD - 1
      NZD2 = NZD2 - 1
      WRITE(STDOUT,98) NZD, NFRAC
C     Extract the maximum flows.
      DO 700 I=1,NZD
        QFVEC(I) = Q(I,NFRAC)
 700  CONTINUE
 
C     Compute the partial free flows.  Add the zero point.
      PFQVEC(1) = 0.0
      MAXQ = Q(NZD,NFRAC)
      DO 800 I=1,NFRAC
        PFQVEC(I+1) = Q(NZD,I)/MAXQ
 800  CONTINUE
      NFRAC = NFRAC + 1
      zrhufd = 0.0
      CALL TWDOUT
     I           (STDOUT, STDTAB, TAB, LABEL, NZD, NFRAC, QFVEC, HDVEC,
     I            PFQVEC, HUMAT, HDATUM,
     I            14,'   WSPRO', zrhufd,
     i            zone, hgrid, vdatum, unitsys, basis,
     i            easting, northing,
     O            EFLAG)
 
      RETURN
 900  CONTINUE
        WRITE(STDOUT,*) ' END OF FILE SEEKING WSPROZQ HEADER'
        STOP 'Abnormal stop. Errors found.'
 
 910  CONTINUE
        WRITE(STDOUT,*) ' End of file seeking "FIRST USER"'
        WRITE(STDOUT,62) 
        STOP 'Abnormal stop. Errors found.'
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
 
      END
C
C
C
      SUBROUTINE   WPRO14_NEW
     I                   (STDIN, STDOUT, STDTAB,
     M                    TABDIR, EFLAG)
 
C     + + + PURPOSE + + +
C     Construct a table of type 14 from one or more .lst files from
C     WSPRO version V061698.  These files will have user specified output that will
C     contain the values needed and in the proper order.  The
C     order is established by the WSPROQZ command implemented in
C     subroutine WPROQZ.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL PARAMETERS + + +
      INTEGER MAXN
      PARAMETER(MAXN=16)
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, IT, KNT, MVAL, NFRAC, NZD, NZD2, NZDF, NZDL, OPT, TAB,
     a        first, hdeflag
      INTEGER CLEN(MAXN), IVAL(MAXN), TERML(MAXN), TERMCLS(MAXN),  
     A        TYPE(MAXN)
      REAL FROUDE, HDATUM, HDVEC(PMXNHU), HUMAT(PMXNHU,PMXFRC), MAXQ,
     A     PFQVEC(PMXFRC), Q(PMXNHU,PMXFRC), QFVEC(PMXNHU), RVAL(MAXN),
     B     WSEL, zrhufd, hd, hdfirst
      real*8 DPVAL(MAXN), easting, northing

      CHARACTER CVAL(MAXN)*64, HEAD*80, LABEL*50, LINE*80, LINE1*81,
     A          NAME*256, TERM(MAXN)*1, XSID*7,
     b   zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
      LOGICAL THERE
 
      CHARACTER GETTOK*7, TABID*16, LINE_IN_ERROR*80

C     + + + EXTERNAL NAMES + + +
      EXTERNAL GETTOK, GETVAL, inline, TABCHK, TWDOUT, os_file_style
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5)
 2    FORMAT(6X,A)
 3    FORMAT(A)
 24   FORMAT(A80)
 25   FORMAT(20X,F10.0,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' Table id for 2-D table of type 14= ',A)
 51   FORMAT(/,' Label for 2-D table: ',A)
 54   FORMAT(/,' ',A80)
 56   FORMAT(' Current file name is:',A)
 58   FORMAT(/,' *ERR:519* Expected Q card but found:',A2)
 60   FORMAT(/,' *WRN:510* Froude number=',F7.2,' > 1.0 at XSID=',A7)
62    FORMAT(' Not enough user tables.  One or more profiles',
     A           ' failed in WSPRO.')
64    FORMAT(/,' *ERR:760* One or more of ERROR or FATAL',
     A      ' found in WSPRO results.',/,4X,'Line is: ',A)
65    format(/,'*ERR:XXX* Tailwater elevations not constant when they',
     a        ' should be.',/,10x,
     b        '  First elevation=',f10.3, ' Differing elevation=',f10.3)
66    format(/,'Errors found prevent output of the type 14 table.')
 96   FORMAT(/,' FILE NAMED:',A,' NOT FOUND.  CHECK SPELLING OF',
     A       ' FILE NAME.')
 98   FORMAT(/,' Found ',I3,' downstream elevations and ',I3,
     A      ' partial maximum flows.')
 99   FORMAT(/,' All files have been processed.')
C***********************************************************************
      hdeflag = 0
C     The latest version of WSPRO finally prints results in a user table
C     to a precision that allows us to use the elevation and the flow
C     in the table.  Previously we could not use the flow because it
C     was printed to such a low precision that the small flows 
C     would not make a valid type 14 table. 

C     Get the table number to be used by the two-D table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE,'TAB',
     O                EFLAG, TABID, TAB)
C      READ(LINE,1,ERR=991) TABID
      WRITE(STDOUT,50) TABID
      CALL TABCHK
     I           (STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)
 
c     Get location items that may be present. If they are not present
c     they will be set to default values.  The default requests FEQUTL
c     to omit the items.  
      call GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)
      call  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

C     Get the label for the 2-D table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) LABEL
      WRITE(STDOUT,51) LABEL
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,54) HEAD
 
C     Input each file name and process in order.
C     Each file must have a descriptive header that agrees exactly with
C     its contents.  This means that each file must represent a
C     complete series of flows for a given downstream head.  The downstream
C     heads must be presented in increasing order and the flows
C     must be presented in increasing order.   Zero flow is not included
C     in any of the files because the result for zero flow is trivial.
 
C     Initialize the global number of zero depths counter
      NZD = 1
      NZD2 = 1
 100  CONTINUE
        CALL inlineb
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,3,ERR=991) NAME
 
        IF(NAME.NE.' ') THEN
           CALL MAYBE_ADD_HOME(
     M                         name)
          call os_file_style(
     m                              name)
          INQUIRE(FILE=NAME, EXIST=THERE)
          IF(THERE) THEN
            WRITE(STDOUT,56) NAME
            OPEN(UNIT=STD48, FILE=NAME, STATUS='OLD')
          ELSE
            WRITE(STDOUT,96) NAME
            STOP 'Abnormal stop. Errors found.'
          ENDIF
C         Check for problems in the WSPRO output
          EFLAG = 0
150       CONTINUE
            READ(STD48,'(A80)',END=160) LINE
            IF(INDEX(LINE,'ERROR').GT.0.OR.
     A         INDEX(LINE,'FATAL').GT.0) THEN
              LINE_IN_ERROR = LINE
              EFLAG = 1
            ENDIF
            GOTO 150
160       CONTINUE
          IF(EFLAG.GT.0) THEN
            WRITE(STDOUT,64) LINE_IN_ERROR 
            STOP 'Abnormal stop.  Error(s) found.'
          ENDIF
          REWIND(STD48)
            
 

C         Initialize the local number of zero depths counter.
          NZDL = 1
 200      CONTINUE
            READ(STD48,24,END=900) LINE
            IF(LINE(14:20).EQ.'WSPROZQ') THEN
C             Found the standard header line.  Read the number of flows,
C             number of downstream elevations, and the head datum.
              READ(LINE,'(40X,I5,5X,I5,8X,F10.0)') NFRAC, NZDF, HDATUM
              WRITE(STDOUT,*) ' NFRAC=',NFRAC, ' NZDF=',NZDF,
     A                ' HDATUM=',HDATUM
 
C             Skip over the two lines defining the user output.
              READ(STD48,24) LINE
              READ(STD48,24) LINE
 
C             Set GETVAL to get the number of values and return the
C             types.  Types must be self defining.  That is, a real
C             value must have a decimal point.
              OPT = 1
              KNT = 0
 300          CONTINUE
                READ(STD48,24) LINE1
C               Search for user quote
                IT = INDEX(LINE1,'''')
                IF(IT.EQ.0) THEN
                  LINE1(81:81) = ''''
                ENDIF
                  
                CALL GETVAL
     I                     (STDOUT, LINE1, MAXN, OPT,
     O                      TYPE, IVAL, RVAL, DPVAL, CVAL, CLEN, 
     O                      EFLAG, TERM, TERML, TERMCLS,
     O                      MVAL)
                DO 310 I=2,MVAL
                  READ(CVAL(I)(1:CLEN(I)),*) RVAL(I)
 310            CONTINUE
                IF(CVAL(1)(1:1).EQ.'Q') THEN
C                 Flow definition.
                  DO 400 I=2,MVAL
                    KNT = KNT + 1
                    Q(NZD,KNT) = RVAL(I)
                    IF(KNT.EQ.NFRAC) THEN
C                     Current downstream depth is complete.
                      NZD = NZD + 1
                      NZDL = NZDL + 1
                      KNT = 0
                      IF(NZDL.GT.NZDF) THEN
C                       Flows from the current file have been completed.
                        NZDL = 1
                        WRITE(STDOUT,*) ' Flows completed'
                        GOTO 401
                      ENDIF
                    ENDIF
 400              CONTINUE
              ELSE
                WRITE(STDOUT,58) CVAL(1)(1:2)
                EFLAG = 1
              ENDIF
              GOTO 300
 401          CONTINUE
 
C             Now look for the string heading the user output.
 
              KNT = 0
              NZDL = 1
c             Clear the flag to track the sequences of downstream depths. 
              first = 0
 500          CONTINUE
                READ(STD48,24,END=910) LINE
                IF(LINE(30:41).EQ.'User Defined') THEN
C                  WRITE(STDOUT,*) ' FOUND: User Defined'
C                 Found the line.  Skip over the blank line and heading.
                  READ(STD48,24) LINE
                  READ(STD48,24) LINE
                  READ(STD48,24) LINE
C                 Read the output lines, and check the Froude number.
C                 The downstream elevation is in the line with cross
C                 section id EXIT and the upstream elevation is in the
C                 line with cross section id APPRO.  Note: in the
C                 latest version there are two APPRO lines: the
C                 first is in the unconstricted case and the 
C                 second is the one we want.  The new version prepends
C                 numbers to the names.  We use those numbers to 
C                 make sure we get what we want. 
 600              CONTINUE
                    READ(STD48,'(A80)') LINE
                    XSID = LINE(2:8)
C                          XSID = GETTOK(XSID)
                    IF(XSID.EQ.'1 EXIT '.OR.XSID.EQ.'5 APPRO'.OR.
     A                 XSID.EQ.'6 APPRO') THEN
                      READ(LINE,25) FROUDE, WSEL
                    ELSE
                      FROUDE = 0.0
                      WSEL = 0.0
                    ENDIF
                    IF(FROUDE.GT.1.0) THEN
                      WRITE(STDOUT,60) FROUDE, XSID
                    ENDIF
                    IF(XSID.EQ.'1 EXIT') THEN
C                     Remember the downstream depth.  Note that this
C                     statement is executed NZDF times with the
C                     same water surface elevation.
 
                      hd = wsel - hdatum
                      HDVEC(NZD2) = hd
                    ENDIF
                    IF(XSID.EQ.'5 APPRO'.OR.XSID.EQ.'6 APPRO') THEN
C                     Get the upstream elevation and save it.
                      KNT = KNT + 1
C                     Leave space for zero fraction of free flow
C                     by adding one to KNT when storing the upstream
C                     heads.  Required by TWDOUT.
                      HUMAT(NZD2,KNT+1) = WSEL - HDATUM
                      if(first == 0) then
c                       First time through here for the current sequence of downstream
c                       heads.  They should all have the same value.
                        first = 1
                        hdfirst = hd
                      else
                        if(hd .ne. hdfirst) then
                          write(stdout,65) hdfirst + hdatum, hd + hdatum
                          hdeflag = hdeflag + 1
                        endif
                      endif
                      IF(KNT.EQ.NFRAC) THEN
C                       The current downstream elevation's flows are
C                       complete.
                        NZD2 = NZD2 + 1
                        NZDL = NZDL + 1
                        KNT = 0
                        first = 0
                        IF(NZDL.GT.NZDF) THEN
C                         Current file is processed.
                          NZDL = 1
                          GOTO 101
                        ENDIF
                      ENDIF
C                     APPRO is the last line of interest in the current
C                     user output.
                      GOTO 601
                    ENDIF
                    GOTO 600
 601              CONTINUE
              ENDIF
              GOTO 500
          ELSE
            GOTO 200
          ENDIF
 101      CONTINUE
          CLOSE(STD48)
          GOTO 100
        ENDIF
 
C     All the files have been processed.
      WRITE(STDOUT,99)
      NZD = NZD - 1
      NZD2 = NZD2 - 1
      WRITE(STDOUT,98) NZD, NFRAC
C     Extract the maximum flows.
      DO 700 I=1,NZD
        QFVEC(I) = Q(I,NFRAC)
 700  CONTINUE
 
C     Compute the partial free flows.  Add the zero point.
      PFQVEC(1) = 0.0
      MAXQ = Q(NZD,NFRAC)
      DO 800 I=1,NFRAC
        PFQVEC(I+1) = Q(NZD,I)/MAXQ
 800  CONTINUE
      NFRAC = NFRAC + 1
      zrhufd = 0.0
      if(hdeflag == 0) then
        CALL TWDOUT
     I             (STDOUT, STDTAB, TAB, LABEL, NZD, NFRAC, QFVEC, 
     I              HDVEC, PFQVEC, HUMAT, HDATUM,
     I              14,'   WSPRO', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)
      else
        write(stdout,66) 
        EFLAG =1
      endif
      RETURN
 900  CONTINUE
        WRITE(STDOUT,*) ' END OF FILE SEEKING WSPROZQ HEADER'
        STOP 'Abnormal stop. Errors found.'
 
 910  CONTINUE
        WRITE(STDOUT,*) ' End of file seeking "FIRST USER"'
        WRITE(STDOUT,62) 
        STOP 'Abnormal stop. Errors found.'
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
 
      END
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
      SUBROUTINE   WPROQZ
     I                   (STDIN, STDOUT,
     O                    EFLAG)
 
C     + + + PURPOSE + + +
C     Compute the flow and elevation cards for input to WSPRO in order
C     to produce the output needed to create a table of type 14 to
C     describe the structure being modeled in WSPRO.  The command
C     WSPROT14 is then used to read the .prt file from WSPRO and
C     create a table of type 14
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, JE, JS, KNT, KNT2, KNTL, LEFLAG, LIM, MAXPRO, N,
     A        NFRAC, NZD, NZDMAX, QCTAB, QNTAB, TAB, TOTAL
      REAL DNSWSE, HDATUM, K, MAXQ, MAXQS, OLDZD, PFQMIN, POWER,
     A     Q(PMXNHU*PMXFRC), QCOFF, QCRIT, QNORM, Y, ZBQC, ZBQN,
     B     ZDMAT(PMXNHU*PMXFRC), ZDVEC(PMXNHU)
      CHARACTER CHAR5*5, HEAD*80, LINE*80, QNTABID*16, QCTABID*16
      character*8 zone, hgrid, vdatum, unitsys, basis
      real*8 easting, northing, shift
 
C     + + + INTRINSICS + + +
      INTRINSIC FLOAT, MAX, MIN, SQRT
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      EXTERNAL CHKTAB, FNDELV, inline, LKTK, LKTQC, LENSTR
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(6X,I5)
 2    FORMAT(7X,F10.0)
 3    FORMAT(4F10.0)
 16   FORMAT(A5,1X,F10.0)
 24   FORMAT(A80)
 25   FORMAT(7X,I5)
 26   FORMAT(A5,1X,I5)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' Cross section table id for normal flow= ',A)
 51   FORMAT(/,' Cross section table id for critical flow=',A)
 52   FORMAT(/,' Elevation of datum for computing heads=',F10.3)
 53   FORMAT(/,' Minimum fraction of maximum flow=',F10.3)
 54   FORMAT(/,' *ERR:520* Negative depth=',F8.2,' at normal-flow',
     A        ' section.')
 55   FORMAT(/,' *ERR:526* Negative depth=',F8.2,' at critical-flow',
     A        ' section.')
 56   FORMAT(/,' *ERR:527* Flow undefined at dsn elevation=',F10.3)
 58   FORMAT(/,'*         WSPROZQ from FEQUTL. NFRAC=',I5,' NZD=',I5,
     A                ' HDATUM=',F10.3)
 60   FORMAT('*                 Q FR#  WSEL')
C 61   FORMAT('J3                5  14    3')
 61   FORMAT('UT                7  27    5')
 62   FORMAT('Q ',8X,F12.2,F12.2,F12.2,F12.2,F12.2)
 63   FORMAT('Q ',8X,F12.2,F12.2,F12.2,F12.2,F12.2,',')
 64   FORMAT('WS',8X,F12.4,F12.4,F12.4,F12.4,F12.4)
 65   FORMAT('WS',8X,F12.4,F12.4,F12.4,F12.4,F12.4,',')
 66   FORMAT(' Power for spacing of partial maximum flows=',F6.2)
 75   FORMAT(/,' Profile limit for WSPRO program=',I5)
 76   FORMAT(/,' *ERR:535* NFRAC=',I5,' > MAXPRO=',I5,'.  Unable',
     A       ' to continue.')
 92   FORMAT(/,' ',A80)
 94   FORMAT(/,' Number of partial maximum flows=',I5)
 95   FORMAT(' ',F10.3,F10.1,F10.7,F10.3)
 96   FORMAT(/,' WSPRO must be run ',I3,' times.  Once for each',
     A    ' input block below.')
97    format(/,' Vertical shift=',f10.3)
C***********************************************************************
C     Get the table number of the cross section that may be used
C     to define the maximum flow using an estimated water-surface slope.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
C      READ(LINE,1,ERR=991) QNTAB
      CALL READ_TABID
     I                (STDOUT, LINE, 'QNTAB',
     O                 EFLAG, QNTABID, QNTAB)
      WRITE(STDOUT,50) QNTABID(1:LENSTR(QNTABID))
C     Remember the internal table number because CHKTAB converts
C     QNTAB to the address and FNDELV needs the table number!
      TAB = QNTAB
      CALL CHKTAB
     I           (20, STDOUT, FTPNT, PMXTAB,
     M            QNTAB,
     O            EFLAG)
      IF(QNTAB.GT.0.AND.EFLAG.EQ.0) THEN
C       Find the bottom elevation for the normal flow cross section
C       table.
        CALL FNDELV
     I             (TAB, STDOUT,
     O              EFLAG, ZBQN)
      ENDIF
 
C     Get the table number of the cross section that defines critical
C     flow in the structure opening using a user-supplied offset
C     in elevation from the downstream water-surface elevation.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
C      READ(LINE,1,ERR=991) QCTAB
      CALL READ_TABID
     I                (STDOUT, LINE, 'QCTAB',
     O                 EFLAG, QCTABID, QCTAB)
      WRITE(STDOUT,51) QCTABID(1:LENSTR(QCTABID))
      TAB = QCTAB
      CALL CHKTAB
     I           (20, STDOUT, FTPNT, PMXTAB,
     M            QCTAB,
     O            EFLAG)
      IF(QCTAB.GT.0.AND.EFLAG.EQ.0) THEN
C       Find the bottom elevation for the critical-flow cross section
C       table.
        CALL FNDELV
     I             (TAB, STDOUT,
     O              EFLAG, ZBQC)
      ENDIF
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
      shift = 0.d0
      if(line(1:5) == 'SHIFT') then
        read(line,'(6x,f10)', err=991) shift
      else
        backspace(stdin)
      endif
      write(stdout,97) shift
 
 
C     Get the elevation of the datum for computing heads
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) HDATUM
      hdatum = hdatum + shift
      WRITE(STDOUT,52) HDATUM
 
C     Get the minimum value of partial maximum flow to use in the
C     code 14 table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) PFQMIN
      WRITE(STDOUT,53) PFQMIN
 
C     Get the values defining the distribution of partial free flows.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,26,ERR=991) CHAR5, NFRAC
      WRITE(STDOUT,94) NFRAC
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) CHAR5, POWER
      WRITE(STDOUT,66) POWER
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,25,ERR=991) MAXPRO
      WRITE(STDOUT,75) MAXPRO
 
      IF(NFRAC.GT.MAXPRO) THEN
        WRITE(STDOUT,76) NFRAC, MAXPRO
        EFLAG = 1
        RETURN
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,24,ERR=991) HEAD
      WRITE(STDOUT,92) HEAD
 
C     Input the defining the downstream heads and the maximum and minimum
C     flows to use for each downstream head.
 
      OLDZD = -1.E30
      NZD = 0
      KNT = 0
 100  CONTINUE
        LEFLAG = 0
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,3,ERR=991) DNSWSE, MAXQ, MAXQS, QCOFF
        dnswse = dnswse + shift
        IF(DNSWSE.GT.OLDZD) THEN
          NZD = NZD + 1
          ZDVEC(NZD) = DNSWSE 
 
C         Compute and output the flow and water-surface elevation
C         lines for WSPRO.  These must be transfered by the user
C         to the WSPRO input prepared by the user.
 
          IF(MAXQ.EQ.0.0) THEN
C           No value supplied by the user.  Try to use the water-
C           surface slope.
            IF(QNTAB.GT.0) THEN
C             Cross section table exists.
              IF(MAXQS.GT.0.0) THEN
C               Slope exists. Make the flow estimate.
                Y = DNSWSE - ZBQN
                IF(Y.LE.0.0) THEN
                  WRITE(STDOUT, 54) Y
                  EFLAG = 1
                  LEFLAG = 1
                  QNORM = 0.0
                ELSE
                  CALL LKTK
     I                     (QNTAB,
     M                      Y,
     O                      K)
                  QNORM = K*SQRT(MAXQS)
                ENDIF
              ELSE
                QNORM = 0.0
              ENDIF
            ELSE
              QNORM = 0.0
            ENDIF
            IF(QCTAB.GT.0) THEN
C             Critical flow table exists.
              Y = DNSWSE + QCOFF - ZBQC
              IF(Y.LE.0.0) THEN
                WRITE(STDOUT,55) Y
                EFLAG = 1
                LEFLAG = 1
                QCRIT = 0.0
              ELSE
                CALL LKTQC
     I                    (QCTAB,
     M                     Y,
     O                     QCRIT)
              ENDIF
            ELSE
              QCRIT = 0.0
            ENDIF
            IF(QCRIT.EQ.0.0.AND.QNORM.EQ.0.0) THEN
              WRITE(STDOUT,56) DNSWSE
              EFLAG = 1
              LEFLAG = 1
            ELSE
              IF(QCRIT.GT.0.0.AND.QNORM.GT.0.0) THEN
                MAXQ = MIN(QCRIT, QNORM)
              ELSE
                MAXQ = MAX(QCRIT, QNORM)
              ENDIF
            ENDIF
          ENDIF
          WRITE(STDOUT,95) DNSWSE, MAXQ, MAXQS, QCOFF
          IF(LEFLAG.EQ.0) THEN
            DO 200 I=1,NFRAC
              KNT = KNT + 1
              Q(KNT) = MAXQ*(PFQMIN + (1.0 - PFQMIN)*
     A                      (FLOAT(I-1)/FLOAT(NFRAC-1))**POWER)
 200        CONTINUE
          ENDIF
 
          OLDZD = DNSWSE
          GOTO 100
        ENDIF
 
C     Create the same pattern in ZDMAT as in Q
      KNT = 0
      DO 500 I=1,NZD
        DO 400 J=1,NFRAC
          KNT = KNT + 1
          ZDMAT(KNT) = ZDVEC(I)
 400    CONTINUE
 500  CONTINUE
 
C     All values have been computed here.  Dump in format required by
C     WSPRO.  Output in groups to accommodate limited profile storage
C     that exists in many WSPRO versions.
 
      NZDMAX = MAXPRO/NFRAC
      N = NZD/NZDMAX
      IF(N*NZDMAX.LT.NZD) THEN
C       Catch the fractional part if it exists.
        N = N + 1
      ENDIF
      WRITE(STDOUT,96) N
      IF(N.EQ.1) NZDMAX = NZD
      KNT = 0
      KNT2 = 0
      DO 1000 I=1,N
        JS = 1 + (I - 1)*NZDMAX
        JE = JS + NZDMAX - 1
        IF(JE.GT.NZD) THEN
C         Last block of values may not have NZDMAX values of
C         downstream elevation.  Adjust value to match.
          NZDMAX = NZD - JS + 1
          JE = NZD
        ENDIF
 
        KNTL = 0
C       Output the standard request for output with a standard heading
        WRITE(STDOUT,58) NFRAC, NZDMAX, HDATUM
        WRITE(STDOUT,60)
        WRITE(STDOUT,61)
 
        TOTAL = NZDMAX*NFRAC
 300    CONTINUE
          IF(KNTL+5.LT.TOTAL) THEN
C           Output five values with a trailing comma
            WRITE(STDOUT,63) (Q(KNT+J), J=1,5)
            KNT = KNT + 5
            KNTL = KNTL + 5
            GOTO 300
          ELSE
C           Output the last record without a trailing comma.
            LIM = MIN(5,TOTAL-KNTL)
            WRITE(STDOUT,62) (Q(KNT+J), J=1,LIM)
            KNT = KNT + LIM
          ENDIF
 
 
        KNTL = 0
 600    CONTINUE
          IF(KNTL+5.LT.TOTAL) THEN
C           Output five values with a trailing comma
            WRITE(STDOUT,65) (ZDMAT(KNT2+J), J=1,5)
            KNT2 = KNT2 + 5
            KNTL = KNTL + 5
            GOTO 600
          ELSE
C           Output the last record without a trailing comma.
            LIM = MIN(5,TOTAL-KNTL)
            WRITE(STDOUT,64) (ZDMAT(KNT2+J), J=1,LIM)
            KNT2 = KNT2 + LIM
          ENDIF
 1000 CONTINUE
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
 
      END
C
C
C
      SUBROUTINE   WPROX
     I                  (STDIN, STDOUT, STDTAB, NFAC,
     M                   TABDIR, FTP,
     O                   EFLAG)
 
C     + + + PURPOSE + + +
C     Abstract WPRO cross sections from a WPRO input deck.
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
      INCLUDE 'wsproxs.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xtadd.cmn'
 
C     + + + SAVED VALUES + + +
      CHARACTER STAOUT*18
      SAVE STAOUT
 
C     + + + LOCAL VARIABLES + + +
      INTEGER BEGTAB, CHFLAG, ENDSTR, FIRST, I, IDSN, IN, ISUB, J, K,
     A        KNT, LNFLAG, MESG, NCON, OLDSUB, SFPN, TABINC, XFLAG,
     B        ITEMP
      REAL BEGSTA, FLOFF(PMXSEC,6), LEFT, RIGHT, SFAC, STADIR,
     A     STAT(PMXSEC,3), STATL, STATR, STATT, TP, XFLINE(6), XOLD,
     B     ZMAX
      CHARACTER BETOPT*8, F10X*5, LINE*120, MODE*8, MONFLG*8, OLDID*2,
     A          OUTOPT*8, SAVOPT*8, SFPNAM*8, XSNAME*5,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8


 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER STRLEN
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSTAB, INSFPX, MKFMT, PROINT, SCNPRO, STRLEN, TABCHK,
     A         TABOUT, GET_INTERNAL_TAB_NUMBER, STRIP_L_BLANKS
 
C     + + + DATA INITIALIZATIONS + + +
      DATA STAOUT/'(''STATION='','/
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(//)
 51   FORMAT(/,' TABID=',I8)
 52   FORMAT(' NAVM=',A5,'  SCALE=',F10.3,' SHIFT=',F10.3)
 53   FORMAT(' STATION=',F10.3,' LEFT=',F10.1,' RIGHT=',F10.1)
 54   FORMAT(' STATION=',F10.3)
 55   FORMAT(1X,'NSUB',I5,6F6.3,1X,/,(10X,6F6.3))
 58   FORMAT(' *ERR:504* NUMBER OF SUBSECTIONS=',I5,' > ',I5)
 61   FORMAT(' *ERR:506* ONLY ONE POINT GIVEN ON BOUNDARY OF THE',
     1  ' CROSS SECTION.')
 67   FORMAT(/,'*WRN:554* Extending left end of cross section',
     A   ' by ',F8.3)
 68   FORMAT(/,' *WRN:555* Extending right end of cross section',
     A   ' by ',F8.3)
 69   FORMAT(/,' *WRN:556* Some point in cross section higher than',
     A  ' either end.',/,10X,'All area above minimum end elevation',
     b  ' is ignored.')
 70   FORMAT('     OFFSET ELEVATION SUBS',3X,A8)
 72   FORMAT(' ',F10.2,F10.2, I5)
 73   FORMAT(' ',F10.2,F10.2, I5,F6.3,8(F5.1,F6.3))
 75   FORMAT('TABID=',I8,2X,A8,1X,A8,1X,A8,1X,A8)
 76   FORMAT('TABID=',I8,'  EXTEND',2X,A8,1X,A8,1X,A8,1X,A8)
 77   FORMAT('FEQXEXT')
 78   FORMAT('VARN=HYDY')
 79   FORMAT('VARN=NCON')
 80   FORMAT('NSUB',I5,6F10.3,/,(9X,6F10.3))
 81   FORMAT('    OFFSET ELEVATION SUBS',2X,A8)
 82   FORMAT(F10.2,F10.2, I5,F6.3,8(F5.1,F6.3))
 84   FORMAT('CHANNEL',/,'SINDEF=LINEAR',/,
     A'HEAD      Left     Leftr      Axis     Axisr    Rightl    ',
     B' Right')
 86   FORMAT('STAT',6F10.4)
 88   FORMAT('OFFS',6F10.2)
 90   FORMAT(I10)
C***********************************************************************
      zone = 'NONE'
      hgrid = 'NONE'
      vdatum = 'NONE'
      unitsys = 'NONE'
      basis = 'NONE'
C     Clear values not set
      GISID = ' '
      NORTHING = 0.D0
      EASTING = 0.D0

C     CLEAR THE ERROR FLAG.  USED TO DETECT PROBLEMS THAT REQUIRE
C     EARLY EXITS FROM PROCESSING
 
      EFLAG = 0
 
C     CLEAR THE VALUES that may not be needed to represent the
C     cross section but are  USED IN THE COMPUTATION OF THE TABLE.
 
      DO 95 I=1,PMXSUB
        NNYU(I) = 0
        NVARU(I) = 0
 95   CONTINUE
 
      DO 96 I=1,PMXPNT
        SNU(I) = 1.0
        LSNU(I) = 0.0
 96   CONTINUE
 
      SNFLGU = 0
 
C     Set the value of the command being processed with its length
      SFPNAM = 'WSPRO'
      SFPN = 5
 
C     INPUT THE CONTROLLING DATA
      CALL INSFPX
     I           (STDIN, STDOUT, SFPNAM, SFPN,
     O            SAVOPT, OUTOPT, BETOPT, MONFLG, EFLAG, MODE, BEGTAB,
     O            TABINC, BEGSTA, STADIR, SFAC,
     O            zone, hgrid, vdatum, unitsys, basis,
     O            easting, northing)
      IF(EFLAG.EQ.1) RETURN
 
C     SET THE INPUT UNIT NUMBER
 
      IN = STD48
 
      IF(MODE.EQ.'CHANNEL'.OR.MODE.EQ.'channel') THEN
C       Set flag for convenience
        CHFLAG = 1
        IDSN = STD50
      ELSE
        CHFLAG = 0
        IDSN = STD49
      ENDIF
 
C     GET THE NEXT CROSS SECTION FROM THE WPRO INPUT.  WHEN A
C     COMPLETE CROSS SECTION HAS BEEN DEFINED, RETURN WITH
C     THE FEQX INPUT AND PROCESS.  RETURN A FLAG WHEN NO
C     CROSS SECTION HAS BEEN FOUND.  THIS TERMINATES THE PROCESSING
C     OF THE FILE.
 
C     START THE TABLE NUMBERS
 
      TABU = BEGTAB
 
C     ENABLE CONVEYANCE WARNING MESSAGES
 
      NOCM = 0
      SLOT = 0.0
 
C     Set the initial station value.  Should be in the same scale
C     as the stations on the WSPRO records.
 
      STATU = BEGSTA
 
C     Initialize the values in common for processing WSPRO input
 
      CALL PROINT
 
C     Clear the flag for the input line being present
      LNFLAG = 0
 
      OLDID = '  '
 
C     Set the flag for the first call
      FIRST = 1
C     Clear the open channel cross section counter.
      KNT = 0
 100  CONTINUE
 
C       Check the table number.
 
        IF(TABU.GE.0) CALL TABCHK
     I                           (STDOUT, PMXTAB,
     M                            TABU, TABDIR, EFLAG)
 
C       Clear the flag for cross section presence.
        XFLAG = 0
 
        CALL SCNPRO
     I             (IN, STDOUT, MXPNTU, STADIR,
     M              STATU, SBU, LSNU, LNFLAG, XFLAG, OLDID, FIRST,
     M              STATL, STATR,
     O              NPNTU, NSUBU, XU, ZU, NU, NVARU, NATYU, YATNU, NNYU,
     O              EFLAG, LEFT, RIGHT, XSNAME, NCON, XFLINE)
 
        IF(XFLAG.EQ.-1) THEN
          CLOSE(STD48)
          IF(MODE.EQ.'indirect'.OR.MODE.EQ.'INDIRECT')THEN
            CLOSE(STD49)
          ELSEIF(MODE.EQ.'CHANNEL'.OR.MODE.EQ.'channel') THEN
C           Output the sinuousity information to the file
C           on STD49 and append to it the material in the
C           file on STD50.  Output the command and the
C           heading information.
            WRITE(STD49,84)
            DO 102 I=1,KNT
              WRITE(STD49,86) STAT(I,1), STAT(I,1), STAT(I,2),
     A           STAT(I,2), STAT(I,3), STAT(I,3)
              WRITE(STD49,88) (FLOFF(I,J),J=1,6)
 102        CONTINUE
            WRITE(STD49,'(A)') 'END'
            REWIND(STD50)
 103        CONTINUE
              READ(STD50,'(A)',END=104) LINE
              ENDSTR = STRLEN(LINE)
              WRITE(STD49,'(A)') LINE(1:ENDSTR)
              GOTO 103
 104        CONTINUE
            WRITE(STD49,'(A)') 'ENDCHAN'
            CLOSE(STD49)
            CLOSE(STD50)
          ENDIF
          RETURN
        ENDIF
        IF(EFLAG.GT.0) RETURN
 
        IF(CHFLAG.EQ.1.AND.XU(1).LT.XU(NPNTU)) THEN
C         Making a CHANNEL command.  Save the sinuousity-defining
C         data.
          KNT = KNT + 1
          STAT(KNT,1) = STATL/SFAC
          STAT(KNT,2) = STATU/SFAC
          STAT(KNT,3) = STATR/SFAC
          DO 105 I=1,6
            FLOFF(KNT,I) = XFLINE(I)
 105      CONTINUE
        ENDIF
        WRITE(STDOUT,50)
        WRITE(STDOUT,51) TABU
        IF(LEFT.GE.RIGHT) THEN
          WRITE(STDOUT,54) STATU/SFAC
        ELSE
          WRITE(STDOUT,53) STATU/SFAC, LEFT, RIGHT
        ENDIF
        IF(NCON.EQ.1) THEN
          WRITE(STDOUT,52) '    0', SCALE, YSHIFT
        ELSE
          WRITE(STDOUT,52) 'HYDY ',SCALE, YSHIFT
        ENDIF
 
        IF(NSUBU.GT.PMXSUB) THEN
          WRITE(STDOUT,58) NSUBU, PMXSUB
          EFLAG = 1
          NSUBU=PMXSUB
        ENDIF
 
        WRITE(STDOUT,55)  NSUBU, (NU(J),J=1,NSUBU)
 
C       CHECK FOR MONOTONICITY
        ZMAX = -9999999.
        XOLD = -1.E20
 
        WRITE(STDOUT,70) XSNAME
 
        OLDSUB = 0
        DO 110 J=1,NPNTU
 
          ISUB = SBU(J)
          IF(ISUB.LT.0) THEN
            WRITE(STDOUT,72) XU(J), ZU(J), ISUB
          ELSE
            IF(ISUB.EQ.OLDSUB) THEN
              WRITE(STDOUT,72) XU(J), ZU(J), ISUB
            ELSE
 
              IF(NVARU(ISUB).GT.0) THEN
C               n varies with depth.  Output the details.
                WRITE(STDOUT,73) XU(J), ZU(J), SBU(J), NATYU(1,ISUB),
     A           (YATNU(K,ISUB), NATYU(K,ISUB),K=2,NNYU(ISUB))
              ELSE
                WRITE(STDOUT,72) XU(J), ZU(J), ISUB
              ENDIF
            ENDIF
          ENDIF
          OLDSUB = ISUB
 
 
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
          IF(ABS(ZMAX - ZU(1)).GT.EPSABS.AND.
     A       ABS(ZMAX - ZU(NPNTU)).GT.EPSABS) THEN
            WRITE(STDOUT,69)
          ENDIF
 
          ZMAX = MIN(ZU(1), ZU(NPNTU))
        ENDIF
 
C       ASSIGN THE VALUES OF N FROM THE SUBSECTIONS TO THE LINE SEGMENT
C       LOCATIONS.
        DO 120 I=1,NPNTU-1
          LSNU(I) = NU(SBU(I))
 120    CONTINUE
 
 
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
            WRITE(TABID,90) TABU
            CALL STRIP_L_BLANKS(
     M                          TABID)
            CALL GET_INTERNAL_TAB_NUMBER
     I                                  (STDOUT, TABID,
     M                                   EFLAG,
     O                                   ITEMP)
c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

      call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

            CALL TABOUT
     I                 (STDOUT, STDTAB, ITEMP, STATT, ZMINU, MESG,
     I                  SAVOPT, OUTOPT, BETOPT,
     i                  zone, hgrid, vdatum, unitsys, basis,
     i                  khflag, alphaflag, betaflag, maflag, mqflag,
     M                  NDEPU, XSTU, FTP)
 200      CONTINUE
        ELSE
C         OUTPUT TO A scratch FILE ATTACHED TO UNIT IDSN
 
          WRITE(IDSN,77)
 
          IF(EXTEND.EQ.0) THEN
            WRITE(IDSN,75) TABU, MONFLG, BETOPT, SAVOPT, OUTOPT
          ELSE
            WRITE(IDSN,76) TABU, MONFLG, BETOPT, SAVOPT, OUTOPT
          ENDIF
 
          TP = STATU/SFAC
          CALL MKFMT
     I              (TP, 10,
     O               F10X)
          STAOUT(13:18) = F10X//')'
          WRITE(IDSN,STAOUT) TP
          IF(NCON.EQ.1) THEN
            WRITE(IDSN,79)
          ELSE
            WRITE(IDSN,78)
          ENDIF
          WRITE(IDSN,80) NSUBU, (NU(J), J=1,NSUBU)
          WRITE(IDSN,81) XSNAME
 
          OLDSUB = 0
          DO 300 J=1,NPNTU
            ISUB = SBU(J)
            IF(ISUB.LT.0) THEN
              WRITE(IDSN,82) XU(J), ZU(J), ISUB
            ELSE
              IF(ISUB.EQ.OLDSUB) THEN
                WRITE(IDSN,82) XU(J), ZU(J), ISUB
               ELSE
 
                IF(NVARU(ISUB).GT.0) THEN
C                 n varies with depth.  Output the details.
                  WRITE(IDSN,82) XU(J), ZU(J), SBU(J), NATYU(1,ISUB),
     A             (YATNU(K,ISUB), NATYU(K,ISUB),K=2,NNYU(ISUB))
                ELSE
                  WRITE(IDSN,82) XU(J), ZU(J), ISUB
                ENDIF
              ENDIF
            ENDIF
            OLDSUB = ISUB
 
 300      CONTINUE
          WRITE(IDSN,50)
        ENDIF
 
        TABU = TABU + TABINC
        GOTO 100
 
      END
