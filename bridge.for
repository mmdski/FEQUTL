C
C
C
      SUBROUTINE   CXSELM
     I                   (ZGIVE, STDOUT, NPNT, NSUB, NAVM, ZMIN, ZMAX,
     I                    X, Z, SB, NFAC, LEFT, RIGHT, LSN, SN,
     M                    EFLAG,
     O                    N, XSV)
 
C     + + + PURPOSE + + +
C     Compute cross section elements at a given elevation with a given
C     subset without changing the cross section data.
 
      IMPLICIT NONE
C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NAVM, NPNT, NSUB, STDOUT
      INTEGER SB(NPNT)
      REAL LEFT, LSN(NPNT), N(*), NFAC, RIGHT, SN(NPNT), X(NPNT),
     A     XSV(PMXELM), Z(NPNT), ZGIVE, ZMAX, ZMIN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ZGIVE  - Water surface elevation for computing elements
C     STDOUT - Fortran unit number for user output and messages
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     NAVM   - Flag for averaging roughness
C     ZMIN   - Minimum elevation
C     ZMAX   - Maximum elevation
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     LEFT   - Defines the left-most offset for a subset to be
C               taken out of a cross section.  If LEFT > RIGHT, then
C               no subset taken.
C     RIGHT  - Right hand limit for subset from a cross section. No
C              subset taken if RIGHT < LEFT.
C     LSN    - Line segment Manning's n value
C     SN     - Sinuousity at a point on a cross section boundary
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     N      - Manning's n values
C     XSV    - Vector of various elements of cross section
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, NPNTS, SNFLGS, WRN557, NSUBS
      INTEGER NNYS(PMXSUB), NVARS(PMXSUB), SBSUB(PMXPNT+2)
      REAL KOLD(PMXSUB), LSNS(PMXPNT+2), NATYS(9,PMXSUB), SNS(PMXPNT+2),
     A     SUBMIN, TSOLD(PMXSUB), XSUB(PMXPNT+2), YATNS(9,PMXSUB),
     B     ZSUB(PMXPNT+2)
      CHARACTER BETOPT*8
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL COMPEL, SUBSET
C***********************************************************************
      BETOPT = 'OLDBETA'
 
C     CHECK FOR SUBSET REQUEST.
 
      SUBMIN = ZMIN
      IF(LEFT.LT.RIGHT) THEN
 
 
C       WE HAVE A SUBSET REQUEST
 
        NPNTS = NPNT + 2
        CALL SUBSET
     I             (STDOUT, NPNT, X, Z, SB, ZMAX, SN, LSN,
     M              LEFT, RIGHT, EFLAG, NPNTS,
     O              XSUB, ZSUB, SBSUB, SNS, LSNS)

C       Added walls are treated as frictionless.  Thus we have added 
C       two subsection numbers. 
        SBSUB(1) = NSUB + 1
        SBSUB(NPNTS-1) = NSUB + 2
        NSUBS = NSUB + 2

C       Clear the values that are used in the cross section
C       computations but that need default values.
 
        DO 100 I=1,NSUBS
          NATYS(1,I) = 0.0
          YATNS(1,I) = 0.0
          NNYS(I) = 0
          NVARS(I) = 0
 100    CONTINUE
 
        SNFLGS = 0
      

C      WRITE(STDOUT,93) LEFT, RIGHT
C93    FORMAT('SUBSET results after subsection adjustment:
C     A      ',/,' LEFT=',F10.2,' RIGHT=',F10.2, ' NSUBS=',I5)
C      DO 500 I=1,NPNTS
C        WRITE(STDOUT,92) XSUB(I),ZSUB(I),SBSUB(I)
C92      FORMAT(2F10.2,I5)
C500   CONTINUE

        SUBMIN = 9999999.0
        DO 140 J=1,NPNTS
          IF(ZSUB(J).LT.SUBMIN) SUBMIN = ZSUB(J)
 140    CONTINUE
        DO 150 I=1,NSUBS
          KOLD(I) = 0.0
          TSOLD(I) = 1.E30
 150    CONTINUE
        WRN557 = 0
        CALL COMPEL
     I             (ZGIVE, NPNTS, NSUBS, NAVM, XSUB, ZSUB, SBSUB, NFAC,
     I              BETOPT, SNFLGS, LSNS, NVARS, NATYS, YATNS, NNYS,
     I              SNS, WRN557,
     M              KOLD, TSOLD,
     O              N, XSV)
        XSV(1) = ZGIVE - SUBMIN
      ELSE
 
        SNFLGS = 0
        DO 160 I=1,NSUB
          KOLD(I) = 0.0
          TSOLD(I) = 1.E30
          NATYS(1,I) = 0.0
          YATNS(1,I) = 0.0
          NNYS(I) = 0
          NVARS(I) = 0
 160    CONTINUE
        WRN557 = 0
        CALL COMPEL
     I             (ZGIVE, NPNT, NSUB, NAVM, X, Z, SB, NFAC, BETOPT,
     I            SNFLGS, LSN, NVARS, NATYS, YATNS, NNYS, SN, WRN557,
     M              KOLD, TSOLD,
     O              N, XSV)
        XSV(1) = ZGIVE - SUBMIN
 
      ENDIF
      RETURN
      END
C
C
C
      REAL FUNCTION   BWFSK
     I                     (SKEW, M, N, SK, SKMAT)
 
C     + + + PURPOSE + + +
C     Compute skew corrections.

      IMPLICIT NONE
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      REAL M, SK(N), SKEW, SKMAT(5,N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     SKEW   - Skew
C     M      - Contraction ratio
C     N      - Number of items in the skew table
C     SK     - Skew
C     SKMAT  - Matrix of skew adjustments
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, JJ
      REAL P, S1, S2
C***********************************************************************
      DO 110 J=2,N
        IF(SKEW.LE.SK(J)) GOTO 120
 110    CONTINUE
 
C     MAY NEED AN ERROR MESSAGE HERE
 
 120  CONTINUE
      P=(SKEW-SK(J-1))/(SK(J)-SK(J-1))
      S1=SKMAT(5,J-1)
      S2=SKMAT(5,J)
      DO 130 I=1,4
        JJ=5-I
        S1=SKMAT(JJ,J-1)+M*S1
        S2=SKMAT(JJ,J)+M*S2
 130    CONTINUE
 
      BWFSK=P*S2+(1.-P)*S1
      RETURN
      END
C
C
C
      SUBROUTINE   STDARG
     I                   (ARGMIN, ARGMAX, N,
     O                    ARG)
 
C     + + + PURPOSE + + +
C     Create a standard argument progression.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      REAL ARG(N), ARGMAX, ARGMIN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ARGMIN - Minimum argument value
C     ARGMAX - Maximum argument value
C     N      - Number of arguments
C     ARG    - Standard argument progression
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL YMAX
 
C     + + + INTRINSICS + + +
      INTRINSIC FLOAT
C***********************************************************************
      YMAX= ARGMAX - ARGMIN
 
      DO 100 I=1,N
        ARG(I) = ARGMIN + YMAX*FLOAT(I-1)/FLOAT(N-1)
 100    CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   INELEV
     I                   (STDIN, STDOUT, ZMIN, ZMAX, NFT,
     O                    ZFT, TAB)
 
C     + + + PURPOSE + + +
C     Input elevations for the flow table.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NFT, STDIN, STDOUT, TAB
      REAL ZFT(NFT), ZMAX, ZMIN
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     ZMIN   - Minimum elevation
C     ZMAX   - Maximum elevation
C     NFT    - Number of values in the bridge loss table
C     ZFT    - Tabulated values of elevation
C     TAB    - Table number
 
C     + + + LOCAL VARIABLES + + +
      INTEGER HEAD(20), EFLAG
      CHARACTER LINE*80, TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL inline, STDARG
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(20A4)
 2    FORMAT(7X,I5)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' ',20A4)
 52   FORMAT(' TABID=',A)
C***********************************************************************
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) HEAD
      WRITE(STDOUT,51) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'TAB',
     O                EFLAG, TABID, TAB)
C      READ(LINE,2,ERR=991) TAB
      WRITE(STDOUT,52) TABID
      CALL STDARG
     I           (ZMIN, ZMAX, NFT,
     O            ZFT)
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   BWCOEF
     I                   (BSCURV, BSKEW, ABTYPE, ADJFAC, PTYPE, M, AP,
     I                    AB,
     O                    BW)
 
C     + + + PURPOSE + + +
C     Compute backwater coef.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ABTYPE, BSCURV, PTYPE
      REAL AB, ADJFAC, AP, BSKEW, BW, M
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     BSCURV - Number for basic backwater coefficient curve(1:3)
C     BSKEW  - Skew of bridge crossing
C     ABTYPE - Code for abutment type
C     ADJFAC - Adjustment factor for backwater coefficient
C     PTYPE  - Code for pier type
C     M      - Contraction ratio
C     AP     - Area of piers
C     AB     - Area of bridge opening excluding piers
C     BW     - Backwater coefficient
 
C     + + + SAVED VALUES + + +
      REAL BSBW(5,3), PRBW1(5,8), PRBW2(3,8), SKA(6), SKB(4),
     A     SKBWA(5,6), SKBWB(5,4)
      SAVE BSBW, PRBW1, PRBW2, SKA, SKB, SKBWA, SKBWB
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL JF, S, S1
 
C     + + + EXTERNAL FUNCTIONS + + +
      REAL BWFSK
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BWFSK
 
C     + + + DATA INITIALIZATIONS + + +
      DATA BSBW/4.139891,-10.52214,14.56526,-14.14747,5.964457,
     1 3.965164,-8.507644,9.610875,-9.256056,4.187662,3.774540,
     2 -5.689646,1.324003,0.289164,0.301939/
      DATA PRBW1/
     1  0.077,0.100,1.69903,5.19906,0.89837,
     2  0.056,0.100,2.29032,9.01086,1.28111,
     3  0.069,0.150,2.57732,5.84647,1.77051,
     4  0.058,0.150,3.12500,9.28954,2.04741,
     5  0.049,0.150,3.78788,14.8297,2.33457,
     6  0.043,0.150,4.62963,26.5409,2.34711,
     7  0.031,0.150,6.57895,56.1367,3.09847,
     8  0.026,0.150,8.06452,88.2802,3.47395/
      DATA PRBW2/
     1 0.233333,1.466666,-0.700000,
     2 0.216111,1.389444,-0.605556,
     3 0.26889,1.57111,-0.70556,
     4 0.052778,1.752778,-0.805556,
     5 -.002222,1.857778,-0.855556,
     6 -.057222,1.962778,-0.905556,
     7 -.098889,1.954444,-0.855556,
     8 -.136111,1.963889,-0.827778/
      DATA SKA/0.,10.,20.,30.,40.,45./
      DATA SKBWA/
     1 0.,0.,0.,0.,0.,
     2 -0.368636,1.49909,-2.550346,2.060606,-0.640692,
     3 -0.664364,2.411550,-3.603463,2.548917,-0.692640,
     4 -0.877909,2.690160,-3.938657,3.286579,-1.160173,
     5 -1.062273,2.230758,-1.222598,-.283550,0.337662,
     6 -1.427182,3.776013,-3.919826,1.822077,-0.251082/
      DATA SKB/0.,15.,30.,45./
      DATA SKBWB/
     1 0.,0.,0.,0.,0.,
     2 0.083000,-2.062048,7.720762,-8.609524,3.380952,
     3 0.737000,-7.943190,22.95381,-24.93810 ,9.190476,
     4 -0.60600,-3.204000,18.04333,-23.56666 ,9.333333/
C***********************************************************************
C     FIND BASE VALUE
 
      S=BSBW(5,BSCURV)
      DO 100 I=1,4
        J=5-I
        S=BSBW(J,BSCURV) +M*S
 100    CONTINUE
 
      IF(S.LT.0.) S=0.
      BW=S
 
C     COMPUTE INCREMENT FOR PIERS
 
 
      IF(PTYPE.EQ.0) GOTO 500
        JF=AP/AB
 
        IF(JF.GT.PRBW1(1,PTYPE)) GOTO 120
          S=JF*(PRBW1(5,PTYPE)+PRBW1(4,PTYPE)*JF)
          GOTO 140
 120    CONTINUE
          S=PRBW1(2,PTYPE)+PRBW1(3,PTYPE)*(JF-PRBW1(1,PTYPE))
 140    CONTINUE
 
        S1=PRBW2(1,PTYPE)+M*(PRBW2(2,PTYPE)+M*PRBW2(3,PTYPE))
        IF(S1.LT.0.) S1=0.
        S=S*S1
        BW=BW+S
 500  CONTINUE
 
C       CORRECT FOR SKEW
 
        IF(BSKEW.LE.0.) GOTO 900
          IF(ABTYPE.EQ.1) GOTO 600
 
C           ABUTMENTS PARALLEL TO FLOW HERE
 
            BW=BW+BWFSK(BSKEW,M,6,SKA,SKBWA)
            GOTO 700
 600      CONTINUE
 
C           ABUTMENTS NOT PARALLEL TO FLOW
 
            BW=BW+BWFSK(BSKEW,M,4,SKB,SKBWB)
 700      CONTINUE
 900    CONTINUE
 
        BW=BW*ADJFAC
      RETURN
      END
C
C
C
      SUBROUTINE   MCOMP
     I                  (STDOUT, CNTR, T, ZW, K, NPNT, NSUB, NAVM, X, Z,
     I                   SB, NFAC, SN, LSN,
     M                   EFLAG,
     O                   N, M)
 
C     + + + PURPOSE + + +
C     Compute the contraction ratio.

      IMPLICIT NONE
      include 'arsize.prm'

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NAVM, NPNT, NSUB, STDOUT
      INTEGER SB(NPNT)
      REAL CNTR, K, LSN(NPNT), M, N(pmxsub), NFAC, SN(NPNT), T, X(NPNT),
     A     Z(NPNT), ZW
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     CNTR   - Offset in the approach section at the main channel
C               point of minimum elevation
C     T      - Bridge opening width
C     ZW     - Water surface elevation for computing contraction ratio
C     K      - conveyance
C     NPNT   - Number of points on boundary of a cross section
C     NSUB   - Number of subsections
C     NAVM   - Flag for averaging roughness
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     SB     - Subsection numbers for the line segments
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     SN     - Sinuousity at a point on a cross section boundary
C     LSN    - Line segment Manning's n value
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     N      - Manning's n values
C     M      - Contraction ratio
 
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, NPNTS, NSUBS, SNFLGS, WRN557
      INTEGER NNYS(PMXSUB), NVARS(PMXSUB), SBS(PMXPNT+2)
      REAL KOLD(PMXSUB), KSUB, LSNS(PMXPNT+2), NATYS(9,PMXSUB),
     A     SNS(PMXPNT+2), TSOLD(PMXSUB), XL, XR, XS(PMXPNT+2), 
     B     XSV(PMXELM), YATNS(9,PMXSUB), ZMAX, ZS(PMXPNT+2)
      CHARACTER BETOPT*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL COMPEL, SUBSET
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' *WRN:509* M > 1 IN MCOMP.K =',1PE12.5,' KSUB =',1PE12.5)
C***********************************************************************
      BETOPT = 'OLDBETA'
C     COMPUTE XL AND XR
 
      XL=CNTR-T/2.
      XR=CNTR+T/2.
 
 
C     DEVELOP SUBSET OF THE BASIC CROSS SECTION
 
      ZMAX = Z(1)

      NPNTS = PMXPNT+2
      CALL SUBSET
     I           (STDOUT, NPNT, X, Z, SB, ZMAX, SN, LSN,
     M            XL, XR, EFLAG, NPNTS, 
     O            XS, ZS, SBS, SNS, LSNS)
 
C      Force added walls to be frictionless
       SBS(1) = NSUB + 1
       SBS(NPNTS-1) = NSUB + 2
       NSUBS = NSUB + 2

C     Set sinuousity and the Manning's n variation values to the
C     default values.
 
      DO 100 I= 1,NSUBS
        NATYS(1,I) = 0.0
        YATNS(1,I) = 0.0
        NNYS(I) = 0
        NVARS(I) = 0
 100  CONTINUE
      SNFLGS = 0
 
C     COMPUTE THE CONVEYANCE IN THE SUBSECTION
 
      DO 110 I=1,NSUBS
        KOLD(I) = 0.0
        TSOLD(I) = 1.E30
 110  CONTINUE
      WRN557 = 0
      CALL COMPEL
     I           (ZW, NPNTS, NSUBS, NAVM, XS, ZS, SBS, NFAC, BETOPT,
     I            SNFLGS, LSNS, NVARS, NATYS, YATNS, NNYS, SNS, WRN557,
     M            KOLD, TSOLD,
     O            N, XSV)
 
      KSUB=XSV(5)**2
      M=KSUB/K
 
 
      IF(M.GT.1.0) THEN
        M=1.0
        IF(ABS(K - KSUB)/K.GT.0.005) THEN
          WRITE(STDOUT,51) K,KSUB
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE   CPTA
     I                 (N, Z, PNUM, PWIDTH,
     O                  TPV, APV)
 
C     + + + PURPOSE + + +
C     Compute top width and area for piers.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER N
      INTEGER PNUM(N)
      REAL APV(N), PWIDTH(N), TPV(N), Z(N)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     N      - Number of items in pier definition
C     Z      - Elevation values for pier widths
C     PNUM   - Number of piers
C     PWIDTH - Pier width
C     TPV    - Vector of values of total width of piers in bridge opening
C     APV    - Vector of values giving area of piers
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
C***********************************************************************
      TPV(1)=PNUM(1)*PWIDTH(1)
      APV(1)=0.
      DO 100 J=2,N
        TPV(J)=PNUM(J)*PWIDTH(J)
        APV(J)=APV(J-1)+(Z(J)-Z(J-1))*(TPV(J-1)+TPV(J))/2.
 100  CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   SKADJ
     I                  (BSKEW, PSKEW, PTYPE, PLEN, NPZ, NPNT,
     M                   PWIDTH, X)
 
C     + + + PURPOSE + + +
C     Adjust opening and piers for skew.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNT, NPZ, PTYPE
      REAL BSKEW, PLEN, PSKEW, PWIDTH(NPZ), X(NPNT)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     BSKEW  - Skew of bridge crossing
C     PSKEW  - Pier skew
C     PTYPE  - Code for pier type
C     PLEN   - Pier length
C     NPZ    - Number of values
C     NPNT   - Number of points on boundary of a cross section
C     PWIDTH - Pier width
C     X      - Offsets defining bridge opening
 
C     + + + SAVED VALUES + + +
      REAL DRAD
      SAVE DRAD
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
      REAL FAC, PW
 
C     + + + INTRINSICS + + +
      INTRINSIC COS, SIN
 
C     + + + DATA INITIALIZATIONS + + +
      DATA DRAD/57.3/
C***********************************************************************
      IF(PSKEW.LE.0.0.OR.PTYPE.EQ.0) GOTO 300
 
C       ADJUST THE WIDTH OF PIER FOR SKEWNESS
 
        FAC=PSKEW/DRAD
        FAC=SIN(FAC)*PLEN
          DO 200 J=1,NPZ
          PW=PWIDTH(J)
          IF(FAC.GT.2.*PW) GOTO 100
            PWIDTH(J)=PW+FAC
            GOTO 110
 100      CONTINUE
            PWIDTH(J)=3.*PW
 110      CONTINUE
 200      CONTINUE
 
 300  CONTINUE
 
      IF(BSKEW.LE.0.0) GOTO 500
 
C       ADJUST BRIDGE OPENING FOR SKEW BY MULTIPLYING
C       ALL OFFSETS BY THE COSINE OF THE SKEW ANGLE.
 
        FAC=BSKEW/DRAD
        FAC=COS(FAC)
        DO 400 J=1,NPNT
         X(J)=X(J)*FAC
 400    CONTINUE
 500  CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   INBRID
     I                   (STDIN, STDOUT,
     M                    TABDIR, EFLAG,
     O                    FTNUM, LABEL)
 
C     + + + PURPOSE + + +
C     Input a bridge description.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTNUM, STDIN, STDOUT
      INTEGER LABEL(20), TABDIR(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTNUM  - function table counter
C     LABEL  - Label for identification
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomb.cmn'
      INCLUDE 'bridge.cmn'
 
C     + + + SAVED VALUES + + +
      CHARACTER RDFLOW*8
      SAVE RDFLOW
 
C     + + + LOCAL VARIABLES + + +
      INTEGER MODE
      INTEGER HEAD(20)
      REAL LEFT, RIGHT, ZMAX, ZOLD
      CHARACTER BETOPT*8, LINE*80, OUTOPT*8, SAVOPT*8, TYPE*8,
     a          zone*8, hgrid*8, vdatum*8, unitsys*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL INFEQX, inline
 
C     + + + DATA INITIALIZATIONS + + +
      DATA RDFLOW/'RDFLOW'/
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(20A4)
 2    FORMAT(12X,I5)
 3    FORMAT(19X,I5)
 4    FORMAT(12X,F10.0)
 5    FORMAT(18X,F10.0)
 6    FORMAT(16X,F10.0)
 7    FORMAT(18X,F10.0)
 8    FORMAT(10X,I5)
 9    FORMAT(12X,F10.0)
 10   FORMAT(10X,F10.0)
 11   FORMAT(F10.0,I5,F10.0)
 12   FORMAT(5X,A8)
 
C     + + + OUTPUT FORMATS + + +
 51   FORMAT(' ',20A4)
 52   FORMAT(' ','BASE CURVE#=',I5)
 53   FORMAT('0*ERR:513 BASE CURVE# > 3 OR < 1.')
 54   FORMAT(' ','ABUTMENT ALIGNMENT=',I5)
 55   FORMAT(' BRIDGE SKEW=',F10.2)
 56   FORMAT('0*WRN:507*BRIDGE SKEW < 0 OR > 45. RESET TO:',F8.1,
     1  ' DEGREES.')
 57   FORMAT(' ADJUSTMENT FACTOR=',F10.2)
 58   FORMAT(' UPSTREAM OFFSET=',F10.2)
 59   FORMAT(' DOWNSTREAM OFFSET=',F10.2)
 60   FORMAT(' PIER TYPE=',I5)
 61   FORMAT(' *ERR:514* PIER TYPE < 0 OR > 8.')
 62   FORMAT(' PIER LENGTH=',F10.2)
 63   FORMAT(' PIER SKEW=',F10.2)
 64   FORMAT(' *WRN:508* PIER SKEW < 0. RESET TO:',F10.2)
 65   FORMAT(' *ERR:515* NUMBER OF ENTRIES IN PIER NUMBER-WIDTH',
     1     ' TABLE >',I5)
 66   FORMAT(' ',F10.2,I5,F10.2)
 67   FORMAT('0*ERR:516* ELEVATION FOR PIER NUMBER-WIDTH TABLE IS',
     1       ' DECREASING AT:',F10.2)
 68   FORMAT(' BRIDGE TYPE=',A8)
 69   FORMAT('0*ERR:517* ONLY TYPE RDFLOW IS VALID FOR A BRIDGE')
C***********************************************************************
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) LABEL
      WRITE(STDOUT,51) LABEL
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,12,ERR=991) TYPE
      WRITE(STDOUT,68) TYPE
      FTYPE=0
      IF(TYPE.EQ.RDFLOW) FTYPE=1
      IF(TYPE.NE.RDFLOW) THEN
        WRITE(STDOUT,69)
        EFLAG = EFLAG + 1
      ENDIF
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,2,ERR=991) BSCURV
      WRITE(STDOUT,52) BSCURV
      IF(BSCURV.GE.1.AND.BSCURV.LE.3) GOTO 100
        BSCURV=1
        WRITE(STDOUT,53)
        EFLAG=EFLAG+1
 100  CONTINUE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,3,ERR=991) ABTYPE
      WRITE(STDOUT,54) ABTYPE
      IF(ABTYPE.GT.0) ABTYPE=1
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) BSKEW
      WRITE(STDOUT,55) BSKEW
      IF(BSKEW.GE.0.AND.BSKEW.LE.45) GOTO 110
        IF(BSKEW.LT.0) BSKEW=ABS(BSKEW)
        IF(BSKEW.GT.45.) BSKEW=45.
        WRITE(STDOUT,56) BSKEW
 110  CONTINUE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,5,ERR=991) ADJFAC
      WRITE(STDOUT,57) ADJFAC
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) CNTRU
      WRITE(STDOUT,58) CNTRU
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,7,ERR=991) CNTRD
      WRITE(STDOUT,59) CNTRD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,8,ERR=991) PTYPE
      WRITE(STDOUT,60) PTYPE
 
      IF(PTYPE.GE.0.AND.PTYPE.LE.8) GOTO 120
        WRITE(STDOUT,61)
        EFLAG=EFLAG+1
        PTYPE=1
 120  CONTINUE
 
      IF(PTYPE.GT.0) GOTO 130
        NPZ=0
        GOTO 180
 130  CONTINUE
 
C     INPUT PIER DESCRIPTION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,9,ERR=991) PLEN
      WRITE(STDOUT,62) PLEN
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,10,ERR=991) PSKEW
      WRITE(STDOUT,63) PSKEW
 
        IF(PSKEW.GE.0.) GOTO 140
          PSKEW=ABS(PSKEW)
          WRITE(STDOUT,64) PSKEW
 140    CONTINUE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) HEAD
      WRITE(STDOUT,51) HEAD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,1,ERR=991) HEAD
      WRITE(STDOUT,51) HEAD
 
C     INPUT THE PIER NUMBER-WIDTH TABLE
 
      NPZ=0
      ZOLD=-1.E30
 
 150  CONTINUE
        NPZ=NPZ+1
        IF(NPZ.LE.MAXNPZ) GOTO 160
          WRITE(STDOUT,65) MAXNPZ
          NPZ=MAXNPZ
          EFLAG=EFLAG+1
 160    CONTINUE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,11,ERR=991) PZ(NPZ),PNUM(NPZ),PWIDTH(NPZ)
      WRITE(STDOUT,66) PZ(NPZ),PNUM(NPZ),PWIDTH(NPZ)
 
      IF(PWIDTH(NPZ).LT.0.) GOTO 180
      IF(PZ(NPZ).GT.ZOLD) GOTO 170
        WRITE(STDOUT,67) ZOLD
        EFLAG=EFLAG+1
 170  CONTINUE
        ZOLD=PZ(NPZ)
        GOTO 150
 
 180  CONTINUE
 
C     ADJUST NPZ FOR THE EXTRA LINE READ
 
      NPZ=NPZ-1
 
C     INPUT A STANDARD CROSS SECTION TABLE-
      MODE = 2
      CALL INFEQX
     I           (STDIN, STDOUT, MXPNTB, MODE,
     M            TABDIR, EFLAG,
     O            TABB, STATB, NPNTB, NSUBB, NAVMB, XB, ZB, SBB, NB,
     O            LEFT, RIGHT, SAVOPT, OUTOPT, BETOPT, ZMAX, zone, 
     O            hgrid, vdatum, unitsys)
 
      FTNUM=TABB
 
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   SPBRID
     I                   (STDIN, STDOUT, STDTAB, NFAC,
     M                    TABDIR, EFLAG)
 
C     + + + PURPOSE + + +
C     Compute special bridge table giving the energy loss coef.
C     for application to the downstream velocity head(FTYPE=0)
C     or to the nominal velocity head in the bridge opening(FTYPE=1).
C     all elevations(in bridge opening, in upstream cross section, and
C     in downstream cross section) are computed from the depth in the
C     bridge opening.  Note that the theory also requires that the
C     elements of the upstream cross section be computed at the
C     downstream depth.  All this implies that the elevation of the
C     minimum point in the bridge opening, in the upstream cross
C     section, and in the downstream cross section must be identical
C     when the generated tables are used in FEQ.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL NFAC
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     NFAC   - Factor in Manning's formula(1.49 or 1.0)
C     TABDIR - Table directory to remember table numbers
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomd.cmn'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'xscomb.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'bridge.cmn'
      INCLUDE 'flotab.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER EFLAGI, FTNUM, FTP, I, J, LTAB, NPZTMP, SNFLGB, WRN557
      INTEGER LABEL(20), NNYB(PMXSUB), NVARB(PMXSUB)
      REAL AB, AD, ALPHAB, ALPHAD, ALPHAU, AP, BWB, BWD, KOLD(PMXSUB),
     A     KU, LEFT, LSNB(PMXPNT), M, NATYB(9,PMXSUB), RIGHT,
     B     SNB(PMXPNT), TB, TBMAX, TP, TSOLD(PMXSUB), XSV(PMXELM),
     C     YATNB(9,PMXSUB), ZDWN, ZUP
      CHARACTER BETOPT*8, OUTOPT*8, SAVOPT*8,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8

 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL BWCOEF, COMPEL, CPTA, CXSTAB, INBRID, INELEV, MCOMP,
     A         SKADJ, SPTOUT, TABOUT, XLKA, XLKALL
 
C     + + + OUTPUT FORMATS + + +
 20     FORMAT('*ERR:511* UPSTREAM CROSS SECTION MISSING FOR BRIDGE')
 30     FORMAT('0*ERR:512* DOWNSTREAM CROSS SECTION MISSING',
     1     ' FOR BRIDGE')
C***********************************************************************
      zone = 'NONE'
      hgrid = 'NONE'
      vdatum = 'NONE'
      unitsys = 'NONE'
      basis = 'NONE'
      
C     ENABLE CONVEYANCE MESSAGES
      NOCM = 0
      SLOT = 1.E30
C     NULLIFY THE SAVING FUNCTION IN TABOUT
      SAVOPT = 'NOSAVE'
      OUTOPT = 'OUT1'
      BETOPT = 'OLDBETA'
      FTP = 0
 
C     ENABLE AUTOMATIC EXTENSION
      EXTEND = 1
 
C     INPUT THE BRIDGE DESCRIPTION
 
      CALL INBRID
     I           (STDIN, STDOUT,
     M            TABDIR, EFLAG,
     O            FTNUM, LABEL)
 
C     CHECK TO MAKE SURE NEEDED UPSTREAM AND DOWNSTREAM CROSS
C     SECTIONS ARE DEFINED.
 
      IF(NPNTU.GT.0) GOTO 50
        WRITE(STDOUT,20)
        STOP 'Abnormal stop. Errors found.'
 50   CONTINUE
        IF(NPNTD.GT.0) GOTO 100
        WRITE(STDOUT,30)
        STOP 'Abnormal stop. Errors found.'
 100  CONTINUE
 
C     SET THE EXTENDED CROSS SECTION DESCRIPTORS TO THE NULL VALUES.
C     THIS ROUTINE IS HOPEFULLY DOOMED SOON!
 
      SNFLGB = 0
      SNFLGU = 0
      DO 90 I=1,NSUBB
        NATYB(1,I) = 0.0
        YATNB(1,I) = 0.0
        NVARB(I) = 0
        NNYB(I) = 0
        NVARU(I) = 0
        NNYU(I) = 0
 90   CONTINUE
      DO 91 I=1,NPNTB-1
        LSNB(I) = NB(SBB(I))
        SNB(I) = 1.0
 91   CONTINUE
      LSNB(NPNTB) = 0
      SNB(NPNTB) = 1.0
      DO 92 I=1,NPNTU-1
        LSNU(I) = NU(SBU(I))
        SNU(I) = 1.0
 92   CONTINUE
      LSNU(NPNTU) = 0
      SNU(NPNTU) = 1.0
 
      EFLAGI = 0
 
C     ADJUST BRIDGE OPENING AND BRIDGE PIERS FOR SKEWNESS
 
C     AVOID PROBLEM W/ ZERO DIMENSION
      NPZTMP=NPZ
      IF(NPZ.LE.0) NPZTMP=1
 
      CALL SKADJ
     I          (BSKEW, PSKEW, PTYPE, PLEN, NPZTMP, NPNTB,
     M           PWIDTH, XB)
 
C     COMPUTE CROSS SECTION TABLE FOR BRIDGE OPENING
 
      RIGHT = 0.0
      LEFT = 0.0
      ZMINB = 9999999.0
      ZMAXB = -9999999.0
      DO 140 J=1,NPNTB
        IF(ZB(J).LT.ZMINB) ZMINB = ZB(J)
        IF(ZB(J).GT.ZMAXB) ZMAXB = ZB(J)
 140  CONTINUE
 
      CALL CXSTAB
     I           (STDOUT, NSUBB, NAVMB, NFAC, MXPNTB, LEFT, RIGHT,
     I            BETOPT, SNFLGB, NVARB, NATYB, YATNB, NNYB,
     M            NPNTB, ZMINB, ZMAXB, XB, ZB, SBB, EFLAGI, LSNB, SNB,
     O            NB, NDEPB, XSTB)
 
      IF(FTYPE.EQ.1) THEN

c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

        call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

        CALL TABOUT
     I             (STDOUT, STDTAB, TABB, STATB, ZMINB, 0,
     I              SAVOPT, OUTOPT, BETOPT,
     i              zone, hgrid, vdatum, unitsys, basis,
     i       khflag, alphaflag, betaflag, maflag, mqflag,
     M              NDEPB, XSTB, FTP)
      ENDIF
 
C     COMPUTE WIDTH AND AREA OF PIERS
 
C      WRITE(STDOUT,*) 'IN SPBRID PTYPE=',PTYPE, ' NPZ=',NPZ
      IF(PTYPE.GT.0) CALL CPTA
     I                        (NPZ, PZ, PNUM, PWIDTH,
     O                         TPV, APV)
 
C      WRITE(STDOUT,*) ' APV(NPZ)=',APV(NPZ)
 
C     INPUT THE ELEVATIONS FOR THE ENERGY LOSS TABLES
 
      NFT=21
 
      CALL INELEV
     I           (STDIN, STDOUT, ZMINB, ZMAXB, NFT,
     O            ZFT, LTAB)
 
C     COMPUTE THE ENTRIES IN THE TABLE. FIRST ELEVATION IS AT
C     BOTTOM OF THE CHANNEL
 
C     USE MAXIMUM BRIDGE OPENING WIDTH ENCOUNTERED FOR COMPUTATION
C     OF M.  TBMAX WILL GIVE THAT VALUE.
 
      TBMAX=0.0
 
      WRN557 = 1
      DO 1000 I=2,NFT
        ZDWN=ZFT(I)
 
C       FIND ELEMENTS AT DOWNSTREAM SECTION
 
        XSV(1)=ZDWN-ZMINB
        CALL XLKALL
     I             (NDEPD, MXPNTD, XSTD,
     M              XSV)
 
        AD=XSV(3)
        ALPHAD=XSV(7)
 
C       FIND ELEMENTS IN BRIDGE OPENING USING DOWNSTREAM DEPTH
 
        XSV(1)=ZDWN-ZMINB
        CALL XLKALL
     I             (NDEPB, MXPNTB, XSTB,
     M              XSV)
        TB=XSV(2)
        AB=XSV(3)
 
C       FIND CURRENT MAXIMUM BRIDGE OPENING WIDTH
 
        IF(TB.GT.TBMAX) TBMAX=TB
 
C       COMPUTE UPSTREAM ELEMENTS AT THE DOWNSTREAM DEPTH
 
        ZUP=ZDWN -ZMINB +ZMINU
        DO 110 J=1,NSUBU
          KOLD(J) = 0.0
          TSOLD(J) = 1.E30
 110    CONTINUE
 
        CALL COMPEL
     I             (ZUP, NPNTU, NSUBU, NAVMU, XU, ZU, SBU, NFAC, BETOPT,
     I          SNFLGU, LSNU, NVARU, NATYU, YATNU, NNYU, SNU, WRN557,
     M              KOLD, TSOLD,
     O              NU, XSV)
        KU=XSV(5)**2
        ALPHAU=XSV(7)
 
C       COMPUTE AREA AND TOP WIDTH OF THE PIERS
 
        ZUP=ZDWN -ZMINB
 
        TP=0.0
        AP=0.0
C        WRITE(STDOUT,*) 'BEFORE CALL TO XLKA PTYPE=',PTYPE, ' NPZ=',NPZ
C        WRITE(STDOUT,*) 'PZ(1)=',PZ(1),' PZ(NPZ)=',PZ(NPZ)
 
        IF(PTYPE.GT.0) CALL XLKA
     I                          (ZUP, NPZ, PZ, TPV, APV,
     O                           TP, AP)
 
C       COMPUTE CONTRACTION RATIO
 
        ZUP=ZDWN-ZMINB+ZMINU
 
        CALL MCOMP
     I            (STDOUT, CNTRU, TBMAX, ZUP, KU, NPNTU, NSUBU, NAVMU,
     I             XU, ZU, SBU, NFAC, SNU, LSNU,
     M             EFLAG,
     O             NU, M)
 
 
C     COMPUTE VELOCITY HEAD COEF IN BRIDGE OPENING
 
        ALPHAB=1.0 +M*(ALPHAU-1.0)
 
C       COMPUTE BACK WATER COEF FOR THE NOMINAL VELOCITY HEAD
C       IN THE BRIDGE OPENING
 
        CALL BWCOEF
     I             (BSCURV, BSKEW, ABTYPE, ADJFAC, PTYPE, M, AP, AB,
     O              BWB)
 
C       COMPUTE THE COEF TO APPLY TO THE DOWNSTREAM VELOCITY HEAD
 
        IF(FTYPE.EQ.1) GOTO 500
        BWD=(ALPHAB*BWB*(AD/AB)**2)/ALPHAD
        GOTO 600
 500    CONTINUE
         BWD=ALPHAB*BWB
 600    CONTINUE
 
C       STORE IN QCRITV
 
        QCRITV(I)=BWD
 1000   CONTINUE
 
C     SET THE ZERO DEPTH ENTRY TO THE SAME VALUE AS THE FIRST
C     NON-ZERO DEPTH
 
      QCRITV(1)=QCRITV(2)
 
C     OUTPUT THE TABLE
 
      CALL SPTOUT
     I           (STDOUT, STDTAB, NFT, LTAB, LABEL, ZFT, QCRITV)
 
      RETURN
      END
