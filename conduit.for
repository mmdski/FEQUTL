C
C
C
      SUBROUTINE   QCLIM
     I                  (STDIN, STDOUT,
     M                   EFLAG)
 
C     + + + PURPOSE + + +
C     Compute an arbitary limit to the critical flow values in
C     a cross section table of type 22 or 25 for a closed conduit.
C     The limit is arbitrary but none the less needed for the
C     computations in EXPCON for closed conduits.
 
C     QCLIM extrapolates from the last two top-widths in the table
C     which are before the start of the slot added to the conduit to
C     maintain a free surface.
 
C     QCLIM does not create any new table.  The table must exist
C     in the storage  system of FEQUTL.  Thus the QCLIM command appears
C     before the EXPCON command to modify the closed conduit tables.
 
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
      INCLUDE 'offcom.cmn'
      INCLUDE 'grvcom.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER ADR, HA, LA, TABLE, TYPE, XOFF
      REAL FACTOR, K, P, QC, QC0, QC1, SC, TNEW, TOLD, Y, Y0, Y1, 
     A     MXSLOT
      CHARACTER LINE*80, NAME*7, TABID*16
 
C     + + + INTRINSICS + + +
      INTRINSIC EXP, LOG
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      EXTERNAL inline, KIL, GET_INTERNAL_TAB_NUMBER, LENSTR,
     A         STRIP_L_BLANKS
 
C     + + + OUTPUT FORMATS + + +
 52   FORMAT(' TABID= ',A)
 56   FORMAT(' *WRN:560* Table may not be of a closed conduit. Slot',
     A  ' width=',F8.4,' > ',F8.4)
 58   FORMAT(' *ERR:644* Table is not for closed conduit. No changes',
     A  ' made.')
 72   FORMAT(' *ERR:645 Cross section table for QCLIMIT must be',
     A  ' TYPE=22 or 25 but TYPE=',I5,' found.')
 81   FORMAT(/,' LIMITING FLOWRATE=',F10.1)
 82   FORMAT(' *ERR:540* TABID= ',A,' NOT FOUND')
 84   FORMAT(/,' CRITICAL SLOPE=',F10.5)
C***********************************************************************
      IF(GRAV.GT.15.0) THEN
        MXSLOT = 0.07
      ELSE
        MXSLOT = 0.02134
      ENDIF
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(7X,A)',ERR=991) TABID
      CALL STRIP_L_BLANKS(
     M                     TABID)
      WRITE(STDOUT,52) TABID(1:LENSTR(TABID))
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, TABID,
     M                             EFLAG,
     O                             TABLE)
      
 
      IF(FTPNT(TABLE).LE.0) THEN
C       MAKE SURE TABLE EXISTS
        WRITE(STDOUT,82)  TABID(1:LENSTR(TABID))
        EFLAG = 1
      ENDIF
 
C     INPUT THE ADJUSTMENT FACTOR FOR THE BUILT IN EXTRAPOLATION
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,'(A7,F10.0)',ERR=991) NAME, FACTOR
      IF(FACTOR.EQ.0.0) FACTOR = 1.0
      WRITE(STDOUT,'(1X,A7,F10.3)') NAME, FACTOR
 
      IF(EFLAG.GT.0) RETURN
 
C     CONVERT THE TABLE NUMBER TO TABLE ADDRESS FOR CONVENIENCE
 
      TABLE = FTPNT(TABLE)
 
C     ENSURE THAT THE TABLES ARE OF THE CORRECT TYPE.
 
      TYPE = ITAB(TABLE+2)
      IF(TYPE.NE.22.AND.TYPE.NE.25) THEN
        WRITE(STDOUT,72) TYPE
        EFLAG = 1
      ENDIF
 
C     SET THE TABLE OFFSET
 
      XOFF = OFFVEC(TYPE)
 
      IF(EFLAG.GT.0) RETURN
 
C     BASIC DATA HAS BEEN INPUT IN VALID FORM.
C     Search for the top width just before the slot is entered.
C     To do so start at the top of the table and search to smaller
C     depths until a top width 1 percent larger than the current
C     top width is found.
 
      LA = TABLE + XTIOFF
 
      HA = ITAB(TABLE)
 
      ADR = HA
 
      TOLD = FTAB(HA+1)
      IF(TOLD.GT.MXSLOT) THEN
        WRITE(STDOUT,56) TOLD, MXSLOT
      ENDIF
 
 100  CONTINUE
 
        ADR = ADR - XOFF
        IF(ADR.EQ.LA) THEN
          WRITE(STDOUT,58)
          EFLAG = 1
          RETURN
        ENDIF
 
        TNEW = FTAB(ADR+1)
        IF((TNEW - TOLD)/TOLD.GT.0.01) THEN
C         ASSUME THAT TNEW IS THE FIRST TOP WIDTH WHICH IS GREATER
C         THAN THE SLOT WIDTH.
 
          QC1 = FTAB(ADR+7)
          Y1 = FTAB(ADR)
          QC0 = FTAB(ADR+7-XOFF)
          Y0 = FTAB(ADR-XOFF)
 
          Y = FTAB(ADR+XOFF)
          P = LOG(Y/Y0)*LOG(QC1/QC0)/LOG(Y1/Y0)
          QC = QC0*EXP(P)*FACTOR
 
          WRITE(STDOUT,81) QC
 
C          WRITE(STDOUT,*) ' Y0=',Y0,' QC0=',QC0,' Y1=',Y1,' QC1=',QC1,
C     A               ' Y=',Y,' QC=',QC
 
          K = FTAB(ADR+3+XOFF)
          SC = QC**2/K**4
          WRITE(STDOUT,84) SC
 
 200      CONTINUE
C           REPLACE VALUES IN THE SLOT
            ADR = ADR + XOFF
            FTAB(ADR+7) = QC
            IF(ADR.LT.HA) GOTO 200
 
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
      SUBROUTINE   MULCON
     I                   (STDIN, STDOUT, STDTAB, NFAC,
     M                    TABDIR, EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute a cross section for multiple conduits.
C     The allowed shapes are:
 
C     CIRC - circular
C     NHE - nominal elliptical with major axis horizontal-ASTM
C     NVE - nominal elliptical with major axis vertical-ASTM
C     RCPA - reinforced concrete pipe arch shape- ASTM
C     THE - true elliptical with major axis horizontal-TYPE=TE
C     TVE - true elliptical with major axis vertical-TYPE=TE
C     BOX - box culvert- true rectangular opening.
C     CMPA- corrugated metal pipe arch from-1980 to present
C     CMPAB - corrugated metal pipe arch before current sizes
C             were introduced. Pre-1980
C     CMPA1- corrugated metal pipe arch with 1 inch corrugations-1980
C            to present
C     CMPA1B- corrugated metal pipe arch with 1 inch corrugations
C            before current sizes were introduced. Pre-1980.
C     SPPA18- structural plate pipe arch with 18 inch corners.
C             6 by 2 inch corrugations
C     SPPA31- structural plate pipe arch with 31 inch corners.
C             6 by 2 inch corrugations.
 
C     Each conduit will have its own slot and the slot width for the
C     aggregate cross section will be the sum of the individual slot
C     widths.
 
C     The number of sides to use for each shape is predetermined and
C     is not under user control.
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
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'xtadd.cmn'
 
C     + + + SAVED VALUES + + +
      INTEGER NURQ
      REAL R1CM(12), R1CM1(15), R1CM1B(15), R1CMB(12), R1RC(17),
     A     R1SP18(34), R1SP31(24), R3CM(12), R3CM1(15), R3CM1B(15),
     B     R3CMB(12), R3RC(17), R3SP18(34), R3SP31(24), RSCM(12),
     C     RSCM1(15), RSCM1B(15), RSCMB(12), RSRC(17), RSSP18(34),
     D     RSSP31(24), SPCM(12), SPCM1(15), SPCM1B(15), SPCMB(12),
     E     SPRC(17), SPSP18(34), SPSP31(24), URQX(11), URQY(11)
      SAVE NURQ, R1CM, R1CM1, R1CM1B, R1CMB, R1RC, R1SP18, R1SP31, R3CM,
     A     R3CM1, R3CM1B, R3CMB, R3RC, R3SP18, R3SP31, RSCM, RSCM1,
     B     RSCM1B, RSCMB, RSRC, RSSP18, RSSP31, SPCM, SPCM1, SPCM1B,
     C     SPCMB, SPRC, SPSP18, SPSP31, URQX, URQY
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, ISTART, IT, ITEND, J, JE, JS, NPIPES, NPNT, NRH,
     A        NSIDES, TAB
      REAL BOTTOM(PMXSUB), CFAC, D, DH, DV, HSLOT, LEFT, MAXK, MAXSOF,
     A     NMUD(PMXSUB), RIGHT, RISE(PMXSUB), ROUGH(PMXSUB),
     B     SPAN(PMXSUB), WSLOT, XOFF, XT(23), YOFF, YT(23), ZEPS,
     C     ZMUD(PMXSUB), ZMUD2, ZOFF
      CHARACTER BETOPT*8, CIN*64, CZEPS*10, LINE*80, MONTON*8, OUTOPT*8,
     A          SAVOPT*8, TYPE(PMXSUB)*10,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8


 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER LENSTR
      CHARACTER GETTOK*10
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CXSTAB, GETTOK, inline, inlineb, MKBOX, MKPIPE, RHARCH, 
     A         RHMAK,  SETMUD, SETOPT, SETSUB, TABCHK, TABOUT, URQMAK,
     B         URQTE, LENSTR, READ_TABID_PLUS
 
C     + + + DATA INITIALIZATIONS + + +
      DATA NURQ/11/
      DATA URQX/0.0000000, 0.0960301, 0.1913097, 0.2850942, 0.3766507,
     A          0.4464486, 0.5086755, 0.5604021, 0.5991936, 0.6232238,
     B          0.6313616/
      DATA URQY/0.4008645, 0.3966157, 0.3839025, 0.3628242, 0.3335457,
     A          0.3060603, 0.2641675, 0.2098394, 0.1456334, 0.0745719,
     B          0.0000000/
      DATA RSRC/
     *  11.0000, 13.5000, 15.5000, 18.0000, 22.5000, 26.6250,
     *  31.3125, 36.0000, 40.0000, 45.0000, 54.0000, 62.0000,
     *  72.0000, 77.2500, 87.1250, 96.8750,106.5000/, SPRC/
     *  18.0000, 22.0000, 26.0000, 28.5000, 36.2500, 43.7500,
     *  51.1250, 58.5000, 65.0000, 73.0000, 88.0000,102.0000,
     * 115.0000,122.0000,138.0000,154.0000,168.7500/, R1RC/
     *  22.8750, 27.5000, 35.5000, 40.6875, 51.0000, 62.0000,
     *  73.0000, 84.0000, 92.5000,105.0000,126.0000,162.5000,
     * 183.0000,218.0000,269.0000,301.3750,329.0000/, R3RC/
     *   4.0313,  5.2500,  5.2500,  4.5938,  6.0313,  6.3750,
     *   7.5625,  8.7500,  9.8125, 11.2188, 12.5625, 13.9688,
     *  19.2813, 20.0625, 22.3750, 24.0000, 26.8750/
      DATA RSCM/
     *  13.0000, 15.0000, 18.0000, 20.0000, 24.0000, 29.0000,
     *  33.0000, 38.0000, 43.0000, 47.0000, 52.0000, 57.0000/, SPCM/
     *  17.0000, 21.0000, 24.0000, 28.0000, 35.0000, 42.0000,
     *  49.0000, 57.0000, 64.0000, 71.0000, 77.0000, 83.0000/,R1CM/
     *  25.6250, 33.1250, 34.6250, 42.2500, 55.1250, 66.1250,
     *  77.2500, 88.2500, 99.2500,110.2500,121.2500,132.2500/,R3CM/
     *   3.5000,  4.1250,  4.8750,  5.5000,  6.8750,  8.2500,
     *   9.6250, 11.0000, 12.3750, 13.7500, 15.1250, 16.5000/
      DATA RSCMB/
     *  11.0000, 13.0000, 16.0000, 18.0000, 22.0000, 27.0000,
     *  31.0000, 36.0000, 40.0000, 44.0000, 49.0000, 54.0000/, SPCMB/
     *  18.0000, 22.0000, 25.0000, 29.0000, 36.0000, 43.0000,
     *  50.0000, 58.0000, 65.0000, 72.0000, 79.0000, 85.0000/, R1CMB/
     *  19.1250, 37.0625, 33.5000, 55.0000, 73.2500, 91.5625,
     *  97.2500,115.6875,129.3125,142.9375,145.5000,154.5000/, R3CMB/
     *   3.5000,  4.0000,  4.0000,  4.5000,  5.0000,  5.5000,
     *   6.0000,  7.0000,  8.0000,  9.0000, 10.0000, 11.0000/
      DATA RSCM1/
     *  32.7500, 38.0000, 43.2500, 48.5000, 54.0000, 58.2500,
     *  62.5000, 67.2500, 71.7500, 76.0000, 80.5000, 84.7500,
     *  89.2500, 93.7500, 98.0000/, SPCM1/
     *  38.5000, 45.0000, 51.7500, 58.5000, 65.0000, 72.5000,
     *  79.0000, 86.5000, 93.5000,101.5000,108.5000,116.5000,
     * 123.5000,131.0000,138.5000/, R1CM1/
     *  31.3750, 36.6250, 43.7500, 51.1250, 56.2500, 63.7500,
     *  82.6250, 92.2500,100.2500,111.6250,120.2500,131.7500,
     * 139.7500,149.5000,162.3750/, R3CM1/
     *  12.5000, 14.5000, 16.6250, 18.7500, 20.7500, 22.8750,
     *  20.8750, 22.6250, 24.3750, 26.1250, 27.7500, 29.5000,
     *  31.2500, 33.0000, 34.7500/
      DATA RSCM1B/
     *  27.0000, 31.0000, 36.0000, 40.0000, 44.0000, 55.0000,
     *  59.0000, 63.0000, 67.0000, 71.0000, 75.0000, 79.0000,
     *  83.0000, 87.0000, 91.0000/, SPCM1B/
     *  43.0000, 50.0000, 58.0000, 65.0000, 72.0000, 73.0000,
     *  81.0000, 87.0000, 95.0000,103.0000,112.0000,117.0000,
     * 128.0000,137.0000,142.0000 /, R1CM1B/
     *  54.7500, 67.0000, 82.0000, 91.2500, 98.5000, 76.2500,
     *  92.7500,100.5000,116.0000,132.2500,151.7500,160.5000,
     * 185.0000,201.0000,210.0000/, R3CM1B/
     *   7.7500,  9.0000, 10.5000, 12.0000, 13.2500, 18.0000,
     *  18.0000, 18.0000, 18.0000, 18.0000, 18.0000, 18.0000,
     *  18.0000, 18.0000, 18.0000/
      DATA RSSP18/
     *  55.0000, 57.0000, 59.0000, 61.0000, 63.0000, 65.0000,
     *  67.0000, 69.0000, 71.0000, 73.0000, 75.0000, 77.0000,
     *  79.0000, 81.0000, 83.0000, 85.0000, 87.0000, 89.0000,
     *  91.0000, 93.0000, 95.0000, 97.0000,100.0000,101.0000,
     * 103.0000,105.0000,107.0000,109.0000,111.0000,113.0000,
     * 115.0000,118.0000,119.0000,121.0000/, SPSP18/
     *  73.0000, 76.0000, 81.0000, 84.0000, 87.0000, 92.0000,
     *  95.0000, 98.0000,103.0000,106.0000,112.0000,114.0000,
     * 117.0000,123.0000,128.0000,131.0000,137.0000,139.0000,
     * 142.0000,148.0000,150.0000,152.0000,154.0000,161.0000,
     * 167.0000,169.0000,171.0000,178.0000,184.0000,186.0000,
     * 188.0000,190.0000,197.0000,199.0000/
      DATA R1SP18/
     *  76.3200, 83.5200, 83.5200,104.1600,136.2000,109.8000,
     * 137.8800,182.8800,141.0000,178.6800,144.6000,177.4800,
     * 227.7600,178.3200,153.2400,180.3600,157.9200,183.2400,
     * 216.3600,186.4800,216.8400,257.4000,314.7600,254.7600,
     * 220.6800,254.1600,297.6000,254.2800,226.8000,255.7200,
     * 291.4800,338.1600,290.8800,332.7600/, R3SP18/
     *  18.0000, 18.0000, 18.0000, 18.0000, 18.0000, 18.0000,
     *  18.0000, 18.0000, 18.0000, 18.0000, 18.0000, 18.0000,
     *  18.0000, 18.0000, 18.0000, 18.0000, 18.0000, 18.0000,
     *  18.0000, 18.0000, 18.0000, 18.0000, 18.0000, 18.0000,
     *  18.0000, 18.0000, 18.0000, 18.0000, 18.0000, 18.0000,
     *  18.0000, 18.0000, 18.0000, 18.0000/
      DATA RSSP31/
     * 112.0000,114.0000,116.0000,118.0000,120.0000,122.0000,
     * 124.0000,126.0000,128.0000,130.0000,132.0000,134.0000,
     * 136.0000,138.0000,140.0000,142.0000,144.0000,146.0000,
     * 148.0000,150.0000,152.0000,154.0000,156.0000,158.0000/, SPSP31/
     * 159.0000,162.0000,168.0000,170.0000,173.0000,179.0000,
     * 184.0000,187.0000,190.0000,195.0000,198.0000,204.0000,
     * 206.0000,209.0000,215.0000,217.0000,223.0000,225.0000,
     * 231.0000,234.0000,236.0000,239.0000,245.0000,247.0000/, R1SP31/
     * 192.6000,219.9600,197.8800,222.6000,256.5600,227.7600,
     * 208.5600,232.0800,260.6400,236.0400,263.1600,240.9600,
     * 266.7600,297.9600,270.6000,299.7600,274.5600,302.2800,
     * 278.6400,305.1600,336.4800,374.2800,338.1600,373.5600/, R3SP31/
     *  31.0000, 31.0000, 31.0000, 31.0000, 31.0000, 31.0000,
     *  31.0000, 31.0000, 31.0000, 31.0000, 31.0000, 31.0000,
     *  31.0000, 31.0000, 31.0000, 31.0000, 31.0000, 31.0000,
     *  31.0000, 31.0000, 31.0000, 31.0000, 31.0000, 31.0000/
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,A)
 4    FORMAT(5X,6A10)
 6    FORMAT(6X,F10.0)
 7    FORMAT(7X,I5,A10)
 8    FORMAT(5X,6F10.0)
 9    FORMAT(5X,6F10.0)
 10   FORMAT(5X,6F10.0)
 12   FORMAT(5X,6F10.0)
 14   FORMAT(5X,6F10.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' TABID=',A,2X,A)
 54   FORMAT(' Type=     ',6(A10),1X,/,(10X,6(A10)))
 56   FORMAT(' Width of slot=',F7.3)
 57   FORMAT(' Number of conduits=',I5,' Mud line EPS=',F7.3)
 58   FORMAT(' Height of slot above datum=',F7.0)
 62   FORMAT(' *ERR:573* Offset space filled in MULCON.')
 64   FORMAT(' *ERR:574* SPAN <=0.0 invalid for TYPE=CIRC.')
 66   FORMAT(' *ERR:528* Number of conduits=',I5,' > ',I5)
 70   FORMAT(' Span=     ',6F10.3)
 71   FORMAT(' Rise=     ',6F10.3)
 72   FORMAT(' Bottom=   ',6F10.3)
 74   FORMAT(' Roughness=',6F10.3)
 76   FORMAT(' *ERR:575* RISE and SPAN are 0.0. one must be >0.0 ')
 78   FORMAT(' *WRN:528* EQIV. DIA. FROM SPAN=',F10.3,
     A   ' DIA. FROM RISE=',F10.3,/,10X,
     B   ' DIFFERENCE INDICATES POSSIBLE ERROR.')
 80   FORMAT(' PROCESSING  TYPE:',1X,A10)
 82   FORMAT(' *ERR:576* TYPE=',A10,' UNKNOWN')
 84   FORMAT(' *ERR:577* RISE AND SPAN MUST BE > 0.0 FOR TYPE TE')
 86   FORMAT('     USING SPAN=',F7.2,' AND RISE=',F7.2)
 92   FORMAT('0*WRN:549* BETA OPTION MUST BE "OLDBETA".  THE OPTION',
     A   ' HAS BEEN RESET TO "OLDBETA".')
 94   FORMAT(' *ERR:616* RISE AND SPAN MUST BE > 0.0 FOR TYPE: BOX')
 96   FORMAT(' MUD LEVL=',6F10.3)
 97   FORMAT(/,' *WRN:562* Slot height=',F8.2,' < ',' maximum soffit',
     A       ' height=',F8.2,/,10X,'Redoing with increased slot',
     B       ' height=',F8.2)
 98   FORMAT(' MUD ROUGH=',6F10.3)
C***********************************************************************
C     ENABLE AUTOMATIC EXTENSION
      EXTEND = 1
 
C     Set the conversion factor
      IF(ABS(NFAC - 1.0).LE.0.0001) THEN
C       Dimensions in meters
        CFAC = 39.37008
      ELSE
        CFAC = 12.0
      ENDIF 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID_PLUS
     I                   (STDOUT, LINE,
     O                    EFLAG, TABID, TAB, CIN)
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), CIN
 
      CALL SETOPT
     I           (STDOUT, CIN,
     O            SAVOPT, OUTOPT, MONTON, BETOPT)
      IF(BETOPT.NE.'OLDBETA') THEN
        WRITE(STDOUT,92)
        BETOPT = 'OLDBETA'
      ENDIF
 
      IF(TAB.LT.0) THEN
        TAB = -TAB
        NOCM = 0
      ELSE
        NOCM = 1
      ENDIF
 
      IF(TAB.GE.0) CALL TABCHK
     I                        (STDOUT, PMXTAB,
     M                         TAB, TABDIR, EFLAG)
 
      TABU = TAB
 
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
      READ(LINE,6,ERR=991) WSLOT
      WRITE(STDOUT,56) WSLOT
 
      SLOT = WSLOT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) HSLOT
      WRITE(STDOUT,58) HSLOT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,7,ERR=991) NPIPES, CZEPS
      IF(CZEPS.EQ.' ') THEN
        IF(NFAC.GT.1.0) THEN
C         English unit set
          ZEPS = 0.02
        ELSE
          ZEPS = 0.006096
        ENDIF
      ELSE
        READ(CZEPS,'(F10.0)') ZEPS
      ENDIF
      WRITE(STDOUT,57) NPIPES, ZEPS
 
      IF(NPIPES.GT.PMXSUB) THEN
        WRITE(STDOUT,66) NPIPES, PMXSUB
        EFLAG = 1
        RETURN
      ENDIF

C     26 October 2006, ddf: Set the gisid to blank.  Uninitialized otherwise!
      gisid = ' '
 
C     Compute the controlling parameter for input of the conduit
C     descriptions.
 
      ITEND = (NPIPES - 1)/6 + 1
      DO 485 IT=1,ITEND
        JS = 1 + 6*(IT - 1)
        JE = JS + 5
        JE = MIN(JE, NPIPES)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,4,ERR=991) (TYPE(J), J=JS,JE)
        WRITE(STDOUT,54) (TYPE(J), J=JS,JE)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,8,ERR=991) (SPAN(J),J=JS,JE)
        WRITE(STDOUT,70) (SPAN(J), J=JS,JE)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,8,ERR=991) (RISE(J),J=JS,JE)
        WRITE(STDOUT,71) (RISE(J), J=JS,JE)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,9,ERR=991) (BOTTOM(J), J=JS,JE)
        WRITE(STDOUT,72) (BOTTOM(J), J=JS,JE)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,10,ERR=991) (ROUGH(J), J=JS,JE)
        WRITE(STDOUT,74) (ROUGH(J), J=JS,JE)
 
C       Look ahead for the mudline input line.  If not there, backspace
C       the input file and set the mud line level to null.
 
        CALL inlineb
     I            (STDIN, STDOUT,
     O             LINE)
        IF(LINE(1:4).NE.'MUDL') THEN
          BACKSPACE STDIN
          DO 480 I=JS,JE
            ZMUD(I) = 0.0
            NMUD(I) = 0.0
 480      CONTINUE
        ELSE
C         READ THE MUDLINE LEVEL
          READ(LINE,12,ERR=991) (ZMUD(J), J=JS,JE)
          WRITE(STDOUT,96) (ZMUD(J), J=JS,JE)
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
          READ(LINE,14,ERR=991) (NMUD(J), J=JS,JE)
          WRITE(STDOUT,98) (NMUD(J), J=JS,JE)
        ENDIF
 
 485  CONTINUE
C     STRIP LEADING BLANKS FROM THE TYPE NAMES
 
      DO 99 I=1,NPIPES
        TYPE(I) = GETTOK(TYPE(I))
 99   CONTINUE
 
C     CONSTRUCT THE POINTS ON THE BOUNDARY OF THE SLOTTED CONDUIT
C     PLACE CONDUITS SIDE BY SIDE EACH WITH A SLOT.  ASSIGN A DIFFERENT
C     SUBSECTION TO EACH.  XOFF WILL BE MADE LARGE ENOUGH
C     SO THAT THE PIPES DO NOT OVERLAP
C     MATCH AREA OF CONDUIT WITH A POLYGON WITH A SUITABLE NUMBER OF
C     SIDES.
C     THE ORIGIN FOR THE COORDINATES WILL BE AT THE DATUM LEVEL. THUS
C     THERE WILL BE SOME NEGATIVE `ELEVATIONS' IN THE DESCRIPTION
C     OF THE PERIMETER.
 
C     CHECK FOR SPACE- ALL SHAPES HAVE THE SAME NUMBER OF SIDES
C     NOT  TRUE FOR BOX CULVERTS  AND ARCH SHAPES.  TAKE MOST
C     DEMANDING VALUE.
 
      NSIDES = 44
      IF((NSIDES+4)*NPIPES.GT.PMXPNT) THEN
        WRITE(STDOUT,62)
        EFLAG = 1
        RETURN
      ENDIF
 
 
C     SINUOSITY CORRECTION NOT NEEDED FOR CLOSED CONDUITS.
      SNFLGU = 0
 
C     Clear values not set elsewhere.  Manning's n is not a variable
C     with depth.
      DO 490 I=1,NPIPES
        NVARU(I) = 0
        NNYU(I) = 0
 490  CONTINUE
 
C     CLEAR THE MAXIMUM SOFFIT VALUE.
 
 
 9000 CONTINUE
        MAXSOF = 0.0
        XOFF = 0.0
        ISTART = 1
        NPNTU = 0
        DO 500 I=1,NPIPES
 
          WRITE(STDOUT,80) TYPE(I)
 
          ZOFF = BOTTOM(I)
 
          IF(TYPE(I).EQ.'CIRC') THEN
 
            NSIDES = 40
 
 
            IF(SPAN(I).LE.0.0) THEN
              WRITE(STDOUT,64)
              EFLAG = 1
              RETURN
            ENDIF
            RISE(I) = SPAN(I)
            WRITE(STDOUT,86) SPAN(I), RISE(I)
 
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL MKPIPE
     I                 (XOFF, ZOFF, NSIDES, SPAN(I), WSLOT, HSLOT,
     O                  NPNT, XU(ISTART), ZU(ISTART))
 
          ELSEIF(TYPE(I).EQ.'NHE') THEN
C           FIND THE EQUIVALENT DIAMETER TO USE
 
            DH = SPAN(I)/1.26
            DV = RISE(I)/0.80
            IF(DH.LE.0.0.AND.DV.LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ELSE IF (DH.LE.0.0) THEN
              D = DV
            ELSE IF(DV.LE.0.0) THEN
              D = DH
            ELSE
C             BOTH > 0.0 HERE- ARE THEY CONSISTENT?
              IF(ABS(DH - DV)/(DH + DV).LE.0.015) THEN
                D = DV
              ELSE
                WRITE(STDOUT,78) DH, DV
                D = DV
              ENDIF
            ENDIF
 
 
            SPAN(I) = 1.26*D
            RISE(I) = 0.8*D
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
C           SCALE THE STANDARD UPPER RIGHT QUADRANT TO THE CURRENT
C           EQUIVALENT DIAMETER
 
            DO 100 J=1,NURQ
              XT(J) = D*URQX(J)
              YT(J) = D*URQY(J)
 100        CONTINUE
 
            YOFF = ZOFF + 0.4*D
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL URQMAK
     I                 (WSLOT, HSLOT, YOFF, XOFF, NURQ, XT, YT,
     O                  XU(ISTART), ZU(ISTART), NPNT)
 
 
          ELSEIF(TYPE(I).EQ.'NVE') THEN
C           FIND THE EQUIVALENT DIAMETER TO USE
 
            DH = SPAN(I)/0.80
            DV = RISE(I)/1.26
            IF(DH.LE.0.0.AND.DV.LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ELSEIF (DH.LE.0.0) THEN
              D = DV
            ELSEIF(DV.LE.0.0) THEN
              D = DH
            ELSE
C             BOTH > 0.0 HERE- ARE THEY CONSISTENT?
              IF(ABS(DH - DV)/(DH + DV).LE.0.015) THEN
                D = DV
              ELSE
                WRITE(STDOUT,78) DH, DV
                D = DV
              ENDIF
            ENDIF
 
            RISE(I) = 1.26*D
            SPAN(I) = 0.8*D
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
C           SCALE THE STANDARD UPPER RIGHT QUADRANT TO THE CURRENT
C           EQUIVALENT DIAMETER AND ROTATE FOR VERTICAL OPTION
 
            DO 102 J=NURQ,1,-1
              XT(NURQ+1-J) = D*URQY(J)
              YT(NURQ+1-J) = D*URQX(J)
 102        CONTINUE
 
            YOFF = ZOFF + 0.63*D
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL URQMAK
     I                 (WSLOT, HSLOT, YOFF, XOFF, NURQ, XT, YT,
     O                  XU(ISTART), ZU(ISTART), NPNT)
 
 
          ELSEIF(TYPE(I).EQ.'RCPA') THEN
C           Reinforced concrete pipe arch
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I                 (STDOUT, 17, RSRC, SPRC, R1RC, R3RC, CFAC,
     M                  RISE(I), SPAN(I),
     O                  EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'CMPA') THEN
C           Corrugated metal pipe arch
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I                 (STDOUT, 12, RSCM, SPCM, R1CM, R3CM, CFAC,
     M                  RISE(I), SPAN(I),
     O                  EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'CMPAB') THEN
C           Corrugated metal pipe arch before current values
C           were introduced.  In use in 1967 and 1971.
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I                 (STDOUT, 12, RSCMB, SPCMB, R1CMB, R3CMB, CFAC,
     M                  RISE(I), SPAN(I),
     O                  EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'CMPA1') THEN
C           Corrugated metal pipe arch with 1 inch corrugations.
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I                 (STDOUT, 15, RSCM1, SPCM1, R1CM1, R3CM1, CFAC,
     M                  RISE(I), SPAN(I),
     O                  EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'CMPA1B') THEN
C           Corrugated metal pipe arch with 1 inch corrugations
C           before current values were introduced.
C           In use in 1967 and 1971.
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I              (STDOUT, 15, RSCM1B, SPCM1B, R1CM1B, R3CM1B, CFAC,
     M               RISE(I), SPAN(I),
     O               EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'SPPA18') THEN
C           Structural plate pipe arch with 18 inch corners.
C           Has 6 inch by 2 inch corrugations.
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I              (STDOUT, 34, RSSP18, SPSP18, R1SP18, R3SP18, CFAC,
     M               RISE(I), SPAN(I),
     O               EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'SPPA31') THEN
C           Structural plate pipe arch with 31 inch corners.
C           Has 6 inch by 2 inch corrugations.
 
            IF(SPAN(I).LE.0.0.AND.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,76)
              EFLAG = 1
              RETURN
            ENDIF
 
C           COMPUTE THE SEMI-PERIMETER FOR THIS SIZE ARCH PIPE
 
            CALL RHARCH
     I               (STDOUT, 24, RSSP31, SPSP31, R1SP31, R3SP31, CFAC,
     M                RISE(I), SPAN(I),
     O                EFLAG, NRH, XT, YT, YOFF)
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            IF(EFLAG.GT.0) RETURN
            YOFF = YOFF + ZOFF
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL RHMAK
     I                (WSLOT, HSLOT, YOFF, XOFF, NRH, XT, YT,
     O                 XU(ISTART), ZU(ISTART), NPNT)
 
          ELSEIF(TYPE(I).EQ.'TE') THEN
 
            IF(SPAN(I).LE.0.0.OR.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,84)
              EFLAG = 1
              RETURN
            ENDIF
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            CALL URQTE
     I                (RISE(I), SPAN(I),
     O                 NURQ, XT, YT)
 
            YOFF = ZOFF + RISE(I)/2.0
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL URQMAK
     I                 (WSLOT, HSLOT, YOFF, XOFF, NURQ, XT, YT,
     O                  XU(ISTART), ZU(ISTART), NPNT)
 
 
          ELSEIF(TYPE(I).EQ.'BOX') THEN
 
            IF(SPAN(I).LE.0.0.OR.RISE(I).LE.0.0) THEN
              WRITE(STDOUT,94) SPAN(I), RISE(I)
              EFLAG = 1
              RETURN
            ENDIF
 
            WRITE(STDOUT,86) SPAN(I), RISE(I)
            MAXSOF = MAX(MAXSOF, RISE(I) + ZOFF)
 
            if(i == 1) then
               xoff = xoff + 0.5*span(i) + 0.25
            else
               xoff = xoff + 0.5*(span(i-1) + span(i)) + 0.25
            endif
            CALL MKBOX
     I                (XOFF, ZOFF, SPAN(I), RISE(I), WSLOT, HSLOT,
     O                 NPNT, XU(ISTART), ZU(ISTART))
 
 
          ELSE
            WRITE(STDOUT,82) TYPE(I)
            EFLAG =1
            RETURN
          ENDIF
 
C         COMPLETE PROCESSING: SET THE SUBSECTION NUMBER, THE LINE
C         SEGMENT ROUGHNESS, AND THE MUDLINE IF ANY.
 
          CALL SETSUB
     I               (I, NPNT, ISTART, ROUGH(I),
     O                SBU, LSNU)
c          write(stdout,91) i, type(i), npnt, istart, rough(i)
c91    format(' i=',i5,' type=',a8,' npnt=',i5,' istart=',i5,
c     a        ' rough(i)=',f10.4)
c          do j=istart,istart+npnt-1
c            write(stdout,89) j, xu(j), zu(j), lsnu(j)
c89    format(i5,f10.4,f10.4,f10.3)
c          enddo
 
          IF(ZMUD(I).GT.0.0) THEN
C           ADJUST THE CROSS SECTION FOR A MUD LINE
            ZMUD2 = ZMUD(I) + ZOFF
            CALL SETMUD
     I                 (STDOUT, ISTART, PMXPNT, NMUD(I), ROUGH(I),
     I                  ZMUD2, ZEPS,
     M                  NPNT, XU, ZU, LSNU, SBU, EFLAG)
            IF(EFLAG.NE.0) RETURN
          ENDIF
 
          ISTART = ISTART + NPNT
          NPNTU = NPNTU + NPNT
c          IF(I.NE.NPIPES) THEN
c            write(stdout,8232) i, span(i), span(i+1)
c8232  format('i=',i5,' span(i)=',f10.4,' span(i+1)=',f10.4)
c            XOFF = XOFF + MAX(SPAN(I), SPAN(I+1)) + 1.0
c          ENDIF
 500    CONTINUE
 
      IF(MAXSOF.GT.HSLOT) THEN
C       PROBLEMS.  INCREASE HSLOT AND DO AGAIN.
        WRITE(STDOUT,97) HSLOT, MAXSOF, MAXSOF + 5.0
        HSLOT = MAXSOF + 5.0
        GOTO 9000
      ENDIF

c     Set sinuosity to 1.0 everywhere
      do i=1,npntu
        snu(i) = 1.0
      end do
 
c      WRITE(STDOUT,*) ' '
c      WRITE(STDOUT,*) ' DUMP OF CROSS-SECTION BOUNDARY'
c      WRITE(STDOUT,*) ' NPNTU=',npntu
c      DO 1000 I=1,NPNTU
c        WRITE(STDOUT,90) XU(I), ZU(I), SBU(I), LSNU(I)
c90    FORMAT(F10.4,F10.3,I5,F10.3)
c1000  CONTINUE
 
      NAVMU = 0
      NSUBU = NPIPES
      DO 600 I=1,NPIPES
        NU(I) = ROUGH(I)
 600  CONTINUE
      LEFT = 1.0
      RIGHT = 0.0
      STATU = 0.0
 
C     COMPUTE ELEMENTS FOR CURRENT CROSS SECTION
 
      IF(NPNTU.GT.1) THEN
C       FIND THE MAXIMUM AND MINIMUM ARGUMENT VALUES
        ZMINU= 9999999.
        ZMAXU = -9999999.
        DO 150 J=1,NPNTU
          IF(ZU(J).GT.ZMAXU) ZMAXU = ZU(J)
          IF(ZU(J).LT.ZMINU) ZMINU = ZU(J)
 150    CONTINUE
 
       CALL CXSTAB
     I            (STDOUT, NSUBU, NAVMU, NFAC, MXPNTU, LEFT, RIGHT,
     I             BETOPT, SNFLGU, NVARU, NATYU, YATNU, NNYU,
     M             NPNTU, ZMINU, ZMAXU, XU, ZU, SBU, EFLAG, LSNU, SNU,
     O             NU, NDEPU, XSTU)
      ENDIF
 
C     FORCE THE CONVEYANCE TO BE MONOTONE INCREASING.
C     limiting the maximum conveyance to the full-flow conveyance.
C     The flow near the crown of the opening is of an uncertain
C     character in most cases.
      IF(NOCM.EQ.1) THEN
C       Find start of the slot.
        DO 168 J=NDEPU-2,2,-1
          IF(XSTU(J-1,2).GT.XSTU(J,2)) THEN
C           J is the entrance to the top slot
            GOTO 169
          ENDIF
 168    CONTINUE
 169    CONTINUE
        MAXK = XSTU(J,5)
        DO 170 J=2,NDEPU
          IF(XSTU(J-1,5).GT.MAXK) THEN
            XSTU(J-1,5) = MAXK
          ENDIF
 170    CONTINUE
      ENDIF
 
C     OUTPUT THE TABLE IF NO ERRORS AND IF OUTPUT IS REQUESTED
 
      IF(EFLAG.NE.0.OR.TABU.EQ.0) GOTO 200
C        IF(ABS(ZMINU).LT.0.001) ZMINU = 0.0
 
C       REMOVE THE TOP VALUE TO ELIMINATE THE
C       ERRATIC TOP WIDTH.
        NDEPU = NDEPU - 1
 
C       FORCE ALPHA AND BETA TO BE 1.00 TO REDUCE PROBLEMS WITH
C       CRITICAL DEPTH.  ALPHA AND BETA ARE CLOSE TO 1.00 IN ANY
C       CASE.
 
        DO 190 J=1,NDEPU
          XSTU(J,6) = 1.0
          XSTU(J,7) = 1.0
 190    CONTINUE
 
c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

        call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

        CALL TABOUT
     I             (STDOUT, STDTAB, TABU, STATU, ZMINU, 0, SAVOPT,
     I              OUTOPT, BETOPT, zone, hgrid, vdatum, unitsys, basis,
     i                    khflag, alphaflag, betaflag, maflag, mqflag,
     M              NDEPU, XSTU, FTP)
 200  CONTINUE
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   PIPES
     I                  (STDIN, STDOUT, STDTAB, NFAC,
     M                   TABDIR, EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute a cross section for several pipes
C     not all of the same diameter nor at the same invert
C     elevation. Each pipe will have its own slot.
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
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'xtadd.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, ISTART, IT, ITEND, J, JE, JS, NPIPES, NPNT, NSIDES,
     A        TAB
      REAL BOTTOM(PMXSUB), DIAM(PMXSUB), HSLOT, LEFT, MAXK, MAXSOF,
     A     NMUD(PMXSUB), RIGHT, ROUGH(PMXSUB), WSLOT, XOFF, ZEPS,
     B     ZMUD(PMXSUB), ZMUD2, ZOFF
      CHARACTER BETOPT*8, CIN*64, CZEPS*10, LINE*80, MONTON*8, OUTOPT*8,
     A          SAVOPT*8,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8


 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, MAX, MIN
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      EXTERNAL CXSTAB, inline, inlineb, MKPIPE, SETMUD, SETOPT, TABCHK,
     A         TABOUT, LENSTR, READ_TABID_PLUS
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,A)
 4    FORMAT(7X,I5)
 6    FORMAT(6X,F10.0)
 7    FORMAT(7X,I5,A10)
 8    FORMAT(5X,6F10.0)
 9    FORMAT(5X,6F10.0)
 10   FORMAT(5X,6F10.0)
 12   FORMAT(5X,6F10.0)
 14   FORMAT(5X,6F10.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' TABID= ',A,2X,A)
 54   FORMAT(' POLYGON APPROX HAS ',I3,' SIDES')
 56   FORMAT(' WIDTH OF SLOT=',F7.3)
 57   FORMAT(' NUMBER OF PIPES=',I5,' MUD LINE EPS=',F7.3)
 58   FORMAT(' HEIGHT OF SLOT ABOVE DATUM=',F7.0)
 62   FORMAT('0*WRN:513* NUMBER OF SIDES=',I3,' TOO SMALL.',
     A       ' NSIDES RESET TO 10')
 64   FORMAT('0*WRN:514* NUMBER OF SIDES=',I4,' TOO LARGE.',
     A       ' NSIDES RESET TO',I5)
 66   FORMAT('0*ERR:528* NUMBER OF PIPES=',I5,' > ',I5)
 70   FORMAT(' DIAMETERS=',6F7.2)
 72   FORMAT(' BOTTOMS=  ',6F7.2)
 74   FORMAT(' ROUGNESS= ',6F7.3)
 76   FORMAT(' MUD LEVEL=',6F7.2)
 78   FORMAT(' MUD ROUGH=',6F7.3)
 92   FORMAT('0*WRN:549* BETA OPTION MUST BE "OLDBETA".  THE OPTION',
     A   ' HAS BEEN RESET TO "OLDBETA".')
C***********************************************************************
C     26 October 2006, ddf: Set the gisid to blank.  Uninitialized otherwise!
      gisid = ' '
C     ALLOW AUTOMATIC EXTENSION
      EXTEND = 1
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID_PLUS
     I                    (STDOUT, LINE,
     O                     EFLAG, TABID, TAB, CIN)
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), CIN
      CALL SETOPT
     I           (STDOUT, CIN,
     O            SAVOPT, OUTOPT, MONTON, BETOPT)
      IF(BETOPT.NE.'OLDBETA') THEN
        WRITE(STDOUT,92)
        BETOPT = 'OLDBETA'
      ENDIF
 
 
      IF(TAB.LT.0) THEN
        TAB = -TAB
        NOCM = 0
      ELSE
        NOCM = 1
      ENDIF
 
      IF(TAB.GE.0) CALL TABCHK
     I                        (STDOUT, PMXTAB,
     M                         TAB, TABDIR, EFLAG)
 
      TABU = TAB
 
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
      READ(LINE,4,ERR=991) NSIDES
      WRITE(STDOUT,54) NSIDES
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) WSLOT
      WRITE(STDOUT,56) WSLOT
 
      SLOT = WSLOT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) HSLOT
      WRITE(STDOUT,58) HSLOT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,7,ERR=991) NPIPES, CZEPS
      IF(CZEPS.EQ.' ') THEN
        ZEPS = EPSDIF
      ELSE
        READ(CZEPS,'(F10.0)') ZEPS
      ENDIF
      WRITE(STDOUT,57) NPIPES, ZEPS
 
      IF(NPIPES.GT.PMXSUB) THEN
        WRITE(STDOUT,66) NPIPES, PMXSUB
        EFLAG = 1
        RETURN
      ENDIF
 
C     Compute the control value for reading the conduit descriptions
 
      ITEND = (NPIPES - 1)/6 + 1
      DO 485 IT=1,ITEND
        JS = 1 + 6*(IT - 1)
        JE = JS + 5
        JE = MIN(JE, NPIPES)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,8,ERR=991) (DIAM(J),J=JS,JE)
        WRITE(STDOUT,70) (DIAM(J), J=JS,JE)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,9,ERR=991) (BOTTOM(J), J=JS,JE)
        WRITE(STDOUT,72) (BOTTOM(J), J=JS,JE)
 
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,10,ERR=991) (ROUGH(J), J=JS,JE)
        WRITE(STDOUT,74) (ROUGH(J), J=JS,JE)
 
C       Look ahead for the mudline input line.  If not there, backspace
C       the input file and set the mud line level to null.
 
        CALL inlineb
     I            (STDIN, STDOUT,
     O             LINE)
        IF(LINE(1:4).NE.'MUDL') THEN
          BACKSPACE STDIN
          DO 480 J=JS,JE
            ZMUD(J) = 0.0
            NMUD(J) = 0.0
 480      CONTINUE
        ELSE
C         READ THE MUDLINE LEVEL
          READ(LINE,12,ERR=991) (ZMUD(J), J=JS,JE)
          WRITE(STDOUT,76) (ZMUD(J), J=JS,JE)
          CALL inline
     I              (STDIN, STDOUT,
     O               LINE)
          READ(LINE,14,ERR=991) (NMUD(J), J=JS,JE)
          WRITE(STDOUT,78) (NMUD(J), J=JS,JE)
        ENDIF
 485  CONTINUE
 
      IF(NSIDES.LT.10) THEN
        WRITE(STDOUT,62) NSIDES
        NSIDES = 10
      ENDIF
      IF(NPIPES*NSIDES.GT.MXPNTU-NPIPES*4) THEN
        WRITE(STDOUT,64) NSIDES*NPIPES, NSIDES*(MXPNTU - 4)
        EFLAG = 1
        RETURN
      ENDIF
 
C     CONSTRUCT THE POINTS ON THE BOUNDARY OF THE SLOTTED PIPE.
C     PLACE PIPES SIDE BY SIDE EACH WITH A SLOT.  ASSIGN A DIFFERENT
C     SUBSECTION TO EACH.  XOFF WILL BE LARGE ENOUGH
C     SO THAT THE PIPES DO NOT OVER LAP
C     MATCH AREA OF PIPE WITH A POLYGON WITH NSIDES SIDES
 
      SNFLGU = 0
      DO 490 I=1,NPIPES
        NNYU(I) = 0
 490  CONTINUE
 
C     MAKE SURE THE SLOT IS HIGH ENOUGH
      MAXSOF = 0.0
      DO 492 I=1,NPIPES
        MAXSOF = MAX(MAXSOF, BOTTOM(I) + DIAM(I))
 492  CONTINUE
      IF(HSLOT.LT.MAXSOF + 5.0) THEN
        HSLOT = MAXSOF + 5.0
      ENDIF
 
      XOFF = 0.0
      ISTART = 1
      NPNTU = 0
      DO 500 I=1,NPIPES
 
        ZOFF = BOTTOM(I)
 
C       CALL MKPIPE WITH THE THE RIGHT OFFSETS USED IN THE
C       ARGUMENTS SO THAT THE VECTORS ARE CREATED INCREMENTALLY
 
        CALL MKPIPE
     I             (XOFF, ZOFF, NSIDES, DIAM(I), WSLOT, HSLOT,
     O              NPNT, XU(ISTART), ZU(ISTART))
 
        DO 400 J=1,NPNT
          SBU(ISTART+J-1) = I
          LSNU(ISTART+J-1) = ROUGH(I)
          SNU(ISTART+J-1) = 1.
 400    CONTINUE
 
        IF(ZMUD(I).GT.0.0) THEN
C         ADJUST THE CROSS SECTION FOR A MUD LINE
 
 
          ZMUD2 = ZMUD(I) + ZOFF
          CALL SETMUD
     I               (STDOUT, ISTART, PMXPNT, NMUD(I), ROUGH(I), ZMUD2,
     I                ZEPS,
     M                NPNT, XU, ZU, LSNU, SBU, EFLAG)
          IF(EFLAG.NE.0) RETURN
 
        ENDIF
 
 
        ISTART = ISTART + NPNT
        NPNTU = NPNTU + NPNT
        IF(I.NE.NPIPES) THEN
          XOFF = XOFF + MAX(DIAM(I), DIAM(I+1)) + 1.0
        ENDIF
 
 500  CONTINUE
 
c      WRITE(STDOUT,*) ' '
c      WRITE(STDOUT,*) ' DUMP OF CROSS-SECTION BOUNDARY'
c      WRITE(STDOUT,*) ' '
c
c      DO 1000 I=1,NPNTU
c        WRITE(STDOUT,90) XU(I), ZU(I), SBU(I), LSNU(I), SNU(I)
c90    FORMAT(F10.4,F10.4,I5, F10.4,f10.2)
c1000  CONTINUE
c      WRITE(STDOUT,*) ' '
 
 
      NAVMU = 0
      NSUBU = NPIPES
      DO 600 I=1,NPIPES
        NU(I) = ROUGH(I)
        NVARU(I) = 0
 600  CONTINUE
      LEFT = 1.0
      RIGHT = 0.0
      STATU = 0.0
 
C     COMPUTE ELEMENTS FOR CURRENT CROSS SECTION
 
      IF(NPNTU.GT.1) THEN
C       FIND THE MAXIMUM AND MINIMUM ARGUMENT VALUES
        ZMINU= 9999999.
        ZMAXU = -9999999.
        DO 150 J=1,NPNTU
          IF(ZU(J).GT.ZMAXU) ZMAXU = ZU(J)
          IF(ZU(J).LT.ZMINU) ZMINU = ZU(J)
 150    CONTINUE
 
 
        CALL CXSTAB
     I             (STDOUT, NSUBU, NAVMU, NFAC, MXPNTU, LEFT, RIGHT,
     I              BETOPT, SNFLGU, NVARU, NATYU, YATNU, NNYU,
     M              NPNTU, ZMINU, ZMAXU, XU, ZU, SBU, EFLAG, LSNU, SNU,
     O              NU, NDEPU, XSTU)
 
      ENDIF
 
 
C     FORCE THE CONVEYANCE TO BE MONOTONE INCREASING.
C     limiting the maximum conveyance to the full-flow conveyance.
C     The flow near the crown of the opening is of an uncertain
C     character in most cases.
      IF(NOCM.EQ.1) THEN
C       Find start of the slot.
        DO 168 J=NDEPU-2,2,-1
          IF(XSTU(J-1,2).GT.XSTU(J,2)) THEN
C           J is the entrance to the top slot
            GOTO 169
          ENDIF
 168    CONTINUE
 169    CONTINUE
        MAXK = XSTU(J,5)
        DO 170 J=2,NDEPU
          IF(XSTU(J-1,5).GT.MAXK) THEN
            XSTU(J-1,5) = MAXK
          ENDIF
 170    CONTINUE
      ENDIF
 
C     OUTPUT THE TABLE IF NO ERRORS AND IF OUTPUT IS REQUESTED
 
      IF(EFLAG.NE.0.OR.TABU.EQ.0) GOTO 200
C        IF(ABS(ZMINU).LT.0.001) ZMINU = 0.0
 
C       REMOVE THE TOP VALUE TO ELIMINATE THE
C       ERRATIC TOP WIDTH.
        NDEPU = NDEPU - 1
 
C       FORCE ALPHA AND BETA TO BE 1.00 TO REDUCE PROBLEMS WITH
C       CRITICAL DEPTH.  ALPHA AND BETA ARE CLOSE TO 1.00 IN ANY
C       CASE.
 
        DO 190 J=1,NDEPU
          XSTU(J,6) = 1.0
          XSTU(J,7) = 1.0
 190    CONTINUE
 
c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

        call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

        CALL TABOUT
     I             (STDOUT, STDTAB, TABU, STATU, ZMINU, 0, SAVOPT,
     I              OUTOPT, BETOPT, zone, hgrid, vdatum, unitsys, basis,
     i                    khflag, alphaflag, betaflag, maflag, mqflag,
     M              NDEPU, XSTU, FTP)
 200  CONTINUE
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   SEWER
     I                  (STDIN, STDOUT, STDTAB, NFAC,
     M                   TABDIR, EFLAG, FTP)
 
C     + + + PURPOSE + + +
C     Compute a sewer cross section complete with upper slot.
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
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     FTP    - next open location in the function table storage
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'xscomu.cmn'
      INCLUDE 'nrdzcm.cmn'
      INCLUDE 'xtadd.cmn'
C     + + + LOCAL VARIABLES + + +
      INTEGER I, ISTART, J, NSIDES, TAB
      REAL D, HSLOT, LEFT, MAXK, N, NMUD, RIGHT, WSLOT, XOFF, ZEPS,
     A     ZMUD, ZOFF
      CHARACTER BETOPT*8, CIN*64, LINE*80, MONTON*8, OUTOPT*8, SAVOPT*8,
     a khflag(pmxpnt)*1, alphaflag(pmxpnt)*1, betaflag(pmxpnt)*1,
     b     maflag(pmxpnt)*1, mqflag(pmxpnt)*1,
     c     zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
     

 
C     + + + INTRINSICS + + +
      INTRINSIC ABS
 
C     + + + EXTERNAL NAMES + + +
      INTEGER LENSTR
      EXTERNAL CXSTAB, inline, inlineb, MKPIPE, SETMUD, SETOPT, TABCHK, 
     A         TABOUT, READ_TABID_PLUS, LENSTR
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,A)
 2    FORMAT(9X,F10.0)
 4    FORMAT(7X,I5)
 6    FORMAT(6X,F10.0)
 8    FORMAT(2X,F10.0)
 10   FORMAT(5X,F10.0)
 12   FORMAT(5X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' TABID= ',A,2X,A)
 52   FORMAT(' DIAMETER=',F8.3)
 54   FORMAT(' POLYGON APPROX HAS',I3,' SIDES')
 56   FORMAT(' WIDTH OF SLOT=',F7.4)
 58   FORMAT(' HEIGHT OF SLOT=',F7.0)
 60   FORMAT(' MANNING N=',F7.3)
 62   FORMAT(/,'*WRN:511* NUMBER OF SIDES=',I3,' TOO SMALL.',
     A       ' NSIDES RESET TO 10')
 64   FORMAT(/,'*WRN:512* NUMBER OF SIDES=',I4,' TOO LARGE.',
     A       ' NSIDES RESET TO',I5)
 66   FORMAT(/,'*WRN:549* BETA OPTION MUST BE "OLDBETA".  THE OPTION',
     A   ' HAS BEEN RESET TO "OLDBETA".')
 68   FORMAT(' MUD LINE LEVEL=',F10.3)
 70   FORMAT(' MUD LINE N=',F10.3)
C***********************************************************************
C     26 October 2006, ddf: Set the gisid to blank.  Uninitialized otherwise!
      gisid = ' '
C     ALLOW AUTOMATIC EXTENSION
 
      EXTEND = 1
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID_PLUS
     I                   (STDOUT, LINE,
     O                    EFLAG, TABID, TAB, CIN)
C      READ(LINE,1,ERR=991) TAB, CIN
      WRITE(STDOUT,50) TABID(1:LENSTR(TABID)), CIN
 
      CALL SETOPT
     I           (STDOUT, CIN,
     O            SAVOPT, OUTOPT, MONTON, BETOPT)
      IF(BETOPT.NE.'OLDBETA') THEN
        WRITE(STDOUT,66)
        BETOPT = 'OLDBETA'
      ENDIF
 
      IF(TAB.LT.0) THEN
        TAB = -TAB
        NOCM = 0
      ELSE
C       SUPPRESS CONVEYANCE MESSAGES AND FORCE MONOTONICITY
        NOCM = 1
      ENDIF
      IF(TAB.GE.0) CALL TABCHK
     I                        (STDOUT, PMXTAB,
     M                         TAB, TABDIR, EFLAG)
      TABU = TAB
 
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
      READ(LINE,2,ERR=991) D
      WRITE(STDOUT,52) D
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,4,ERR=991) NSIDES
      WRITE(STDOUT,54) NSIDES
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) WSLOT
      WRITE(STDOUT,56) WSLOT
 
      SLOT = WSLOT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,6,ERR=991) HSLOT
      WRITE(STDOUT,58) HSLOT
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,8,ERR=991) N
      WRITE(STDOUT,60) N
 
 
C     LOOK AHEAD FOR THE MUDLINE INPUT LINE.  IF NOT THERE, BACKSPACE THE
C     INPUT FILE AND SET THE MUD LINE LEVEL TO NULL.
 
      CALL inlineb
     I          (STDIN, STDOUT,
     O           LINE)
      IF(LINE(1:4).NE.'MUDL') THEN
        BACKSPACE STDIN
        ZMUD = 0.0
        NMUD = 0.0
      ELSE
C       READ THE MUDLINE LEVEL
        READ(LINE,10,ERR=991) ZMUD
        WRITE(STDOUT,68) ZMUD
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        READ(LINE,12,ERR=991) NMUD
        WRITE(STDOUT,70) NMUD
        ZEPS = 0.0
      ENDIF
 
      IF(NSIDES.LT.10) THEN
        WRITE(STDOUT,62) NSIDES
        NSIDES = 10
      ENDIF
      IF(NSIDES.GT.MXPNTU-4) THEN
        WRITE(STDOUT,64) NSIDES, MXPNTU - 4
        NSIDES = MXPNTU - 4
      ENDIF
 
C     CONSTRUCT THE POINTS ON THE BOUNDARY OF THE SLOTTED PIPE.
C     MATCH AREA OF PIPE WITH A POLYGON WITH NSIDES SIDES
 
      IF(HSLOT.LT.3.*D) THEN
        HSLOT = 3.*D 
      ENDIF
      XOFF = 0.0
      ZOFF = 0.0
      CALL MKPIPE
     I           (XOFF, ZOFF, NSIDES, D, WSLOT, HSLOT,
     O            NPNTU, XU, ZU)
 
      NAVMU = 0
      NU(1) = N
      NSUBU = 1
      LEFT = 1.0
      RIGHT = 0.0
      STATU = 0.0
      NVARU(1) = 0
      SNFLGU = 0
      NNYU(1) = 0
 
      DO 100 I=1,NPNTU
        SBU(I) = 1
        LSNU(I) = NU(1)
        SNU(I) = 1.
 100  CONTINUE
 
      IF(ZMUD.GT.0.0) THEN
C       ADJUST THE CROSS SECTION FOR A MUD LINE
 
        ISTART = 1
        CALL SETMUD
     I             (STDOUT, ISTART, PMXPNT, NMUD, N, ZMUD, ZEPS,
     M              NPNTU, XU, ZU, LSNU, SBU, EFLAG)
        IF(EFLAG.NE.0) RETURN
 
      ENDIF
C     COMPUTE ELEMENTS FOR CURRENT CROSS SECTION
 
      IF(NPNTU.GT.1) THEN
C       FIND THE MAXIMUM AND MINIMUM ARGUMENT VALUES
        ZMINU= 9999999.
        ZMAXU = -9999999.
        DO 150 J=1,NPNTU
          IF(ZU(J).GT.ZMAXU) ZMAXU = ZU(J)
          IF(ZU(J).LT.ZMINU) ZMINU = ZU(J)
 150    CONTINUE
 
 
        CALL CXSTAB
     I             (STDOUT, NSUBU, NAVMU, NFAC, MXPNTU, LEFT, RIGHT,
     I              BETOPT, SNFLGU, NVARU, NATYU, YATNU, NNYU,
     M              NPNTU, ZMINU, ZMAXU, XU, ZU, SBU, EFLAG, LSNU, SNU,
     O              NU, NDEPU, XSTU)
      ENDIF
 
 
C     FORCE THE CONVEYANCE TO BE MONOTONE INCREASING.
C     limiting the maximum conveyance to the full-flow conveyance.
C     The flow near the crown of the opening is of an uncertain
C     character in most cases.
      IF(NOCM.EQ.1) THEN
C       Find start of the slot.
        DO 168 J=NDEPU-2,2,-1
          IF(XSTU(J-1,2).GT.XSTU(J,2)) THEN
C           J is the entrance to the top slot
            GOTO 169
          ENDIF
 168    CONTINUE
 169    CONTINUE
        MAXK = XSTU(J,5)
        DO 170 J=2,NDEPU
          IF(XSTU(J-1,5).GT.MAXK) THEN
            XSTU(J-1,5) = MAXK
          ENDIF
 170    CONTINUE
      ENDIF
 
C     OUTPUT THE TABLE IF NO ERRORS AND IF OUTPUT IS REQUESTED
 
      IF(EFLAG.NE.0.OR.TABU.EQ.0) GOTO 200
C        IF(ABS(ZMINU).LT.0.001) ZMINU = 0.0

c     compute derivatives of square root of conveyance, alpha, beta, 
c     da, and dq

          call xsecfit
     I            (STDOUT, 0,
     M             NDEPU, XSTU, 
     o             khflag, alphaflag, betaflag, maflag, mqflag)

        CALL TABOUT
     I             (STDOUT, STDTAB, TABU, STATU, ZMINU, 0, SAVOPT,
     I              OUTOPT, BETOPT, zone, hgrid, vdatum, unitsys, basis,
     i                    khflag, alphaflag, betaflag, maflag, mqflag,
     M              NDEPU, XSTU, FTP)
 200  CONTINUE
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop. Errors found.'
      END
C
C
C
      SUBROUTINE   MKPIPE
     I                   (XOFF, ZOFF, NSIDES, D, WSLOT, HSLOT,
     O                    NPNTS, X, Z)
 
C     + + + PURPOSE + + +
C     Construct polygon with NSIDES with slot added.  The slot
C     height is measured from the same datum as is ZOFF.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNTS, NSIDES
      REAL D, HSLOT, WSLOT, X(*), XOFF, Z(*), ZOFF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     XOFF   - Offset between conduits to prevent overlap
C     ZOFF   - Vertical offset of invert from culvert datum
C     NSIDES - Number of sides in the approximating polygon
C     D      - Maximum vertical extent of a closed conduit
C     WSLOT  - Width of slot in closed conduit to preserve free surface
C     HSLOT  - Height of slot for a closed conduit
C     NPNTS  - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL XL, ZL
      DOUBLE PRECISION A, R, THETA
 
C     + + + INTRINSICS + + +
      INTRINSIC COS, DBLE, SIN, SQRT
C***********************************************************************
C     FIND CENTRAL ANGLE-- 2*PI = 6.283185307179586D0
 
      THETA = 6.283185307179586D0/DBLE(NSIDES)
 
C     FIND RADIUS OF CIRCUMSCRIBED CIRCLE FOR THE POLYGON SO THAT THE
C     POLYGON WILL HAVE THE SAME AREA AS THE SEWER
C     PI/2 = 1.570796...
 
      R = SQRT(1.57079632679489D0*DBLE(D)**2/
     A                          (DBLE(NSIDES)*SIN(THETA)))
 
C     FIND VALUES FOR THE SLOT LOCATION
 
      XL = -WSLOT/2.0
      ZL = R + (1.0 - COS(THETA))*XL/SIN(THETA) + 0.5*D + ZOFF
 
C     ADD AN EXTRA POINT WHICH IS LATER DELETED TO REMOVE NOISE FROM
C     THE TOP OF THE SLOT
 
      X(1) = XL + XOFF
      Z(1) =  HSLOT + 0.01
      X(2) = XL + XOFF
      Z(2) = HSLOT
      X(3) = XL + XOFF
      Z(3) = ZL
 
      A = THETA
      DO 100 I=4,NSIDES+2
        X(I) = -R*SIN(A) + XOFF
        Z(I) = R*COS(A) + 0.5*D + ZOFF
        A = A + THETA
 100  CONTINUE
 
      X(NSIDES+3) = -XL + XOFF
      Z(NSIDES+3) = ZL
      X(NSIDES+4) = -XL + XOFF
      Z(NSIDES+4) =  HSLOT
 
      NPNTS = NSIDES + 4
 
      RETURN
 
      END
C
C
C
      SUBROUTINE   SETMUD
     I                   (STDOUT, ISTART, MAXPNT, NMUD, NCON, ZMUD,
     I                    ZEPS,
     M                    NPNT, X, Z, LSN, SB, EFLAG)
 
C     + + + PURPOSE + + +
C     Set a mud line into a simple cross section shape.  SETMUD
C     assumes that there will be only two intersections: one on the
C     left and one on the right.  It also assumes that the mud line
C     is really present, i. e. it will report a bug if called
C     with a mudline level that matches the minimum point in the
C     cross section boundary.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, ISTART, MAXPNT, NPNT, STDOUT
      INTEGER SB(MAXPNT)
      REAL LSN(MAXPNT), NCON, NMUD, X(MAXPNT), Z(MAXPNT), ZEPS, ZMUD
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     ISTART - Start point for a conduit cross section description
C     MAXPNT - Maximum number of tabulated values in a cross section
C               function table
C     NMUD   - Manning;s n for the mud line surface
C     NCON   - Manning's n value for the conduit
C     ZMUD   - Elevation of mud line
C     ZEPS   - Increment in elevation used to put a slope on the
C               mudline
C     NPNT   - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary
C     LSN    - Line segment Manning's n value
C     SB     - Subsection numbers for the line segments
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
C     + + + LOCAL VARIABLES + + +
      INTEGER EXACTL, EXACTR, I, ILEFT, IRIGHT, NP
      REAL XL, XLEFT, XR, XRIGHT, ZL, ZLEFT, ZMUDL, ZMUDR, ZR, ZRIGHT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:681* Space of',I5,' points in SETMUD too small.')
C***********************************************************************
C     PUT SOME SLOPE ON THE MUD LINE TO PREVENT LATER CODE FROM
C     ADJUSTING ONE END TO PREVENT A HORIZONTAL LINE IN THE CROSS
C     SECTION BOUNDARY THAT IS NOT AT THE MINIMUM ELEVATION IN THE
C     CROSS SECTION.
 
      ZMUDL = ZMUD + ZEPS
      ZMUDR = ZMUD - ZEPS
 
C     SCAN THE LINE SEGMENTS UNTIL A LEFT AND RIGHT INTERSECTION WITH
C     THE MUD LINE LEVEL HAVE BEEN FOUND.  THE LEFT INTERSECTION WILL
C     BE FOUND FIRST.
 
      XL = X(ISTART)
      ZL = Z(ISTART)
      DO 100 I=ISTART+1,ISTART+NPNT-1
        XR = X(I)
        ZR = Z(I)
        IF(ZL.GE.ZMUDL.AND.ZR.LT.ZMUDL) THEN
C         FOUND LEFT INTERSECTION
          IF(ZL.EQ.ZMUDL) THEN
            EXACTL = 1
            ILEFT = I - 1
          ELSE
            EXACTL = 0
            ILEFT = I -1
            ZLEFT = ZMUDL
            XLEFT = XL + (ZLEFT - ZL)*(XR - XL)/(ZR - ZL)
          ENDIF
        ENDIF
        IF(ZR.GE.ZMUDR.AND.ZL.LT.ZMUDR) THEN
C         FOUND INTERSECTION ON RIGHT
          IF(ZR.EQ.ZMUDR) THEN
            EXACTR = 1
            IRIGHT = I
          ELSE
            EXACTR = 0
            IRIGHT = I
            ZRIGHT = ZMUDR
            XRIGHT = XL + (ZRIGHT - ZL)*(XR - XL)/(ZR - ZL)
          ENDIF
C         BOTH INTERSECTIONS HAVE BEEN FOUND. EXIT THE LOOP.
          GOTO 110
        ENDIF
        XL = XR
        ZL = ZR
 100  CONTINUE
 
      WRITE(STDOUT,*) ' *BUG:XXX* Intersection failure in SETMUD.'
      STOP 'Abnormal stop. Errors found.'
 
 110  CONTINUE
 
C      WRITE(STDOUT,*) ' ILEFT=',ILEFT,' EXACTL=',EXACTL
C      IF(EXACTL.EQ.0) THEN
C        WRITE(STDOUT,*) ' XLEFT=',XLEFT,' ZLEFT=',ZLEFT
C      ENDIF
C      WRITE(STDOUT,*) ' IRIGHT=',IRIGHT,' EXACTR=',EXACTR
C      IF(EXACTR.EQ.0) THEN
C        WRITE(STDOUT,*) ' XRIGHT=',XRIGHT,' ZRIGHT=',ZRIGHT
C      ENDIF
 
C     ADJUST THE CROSS SECTION DESCRIPTION TO INCLUDE THE MUD LINE.
C     MAKE ROOM FOR ADDED POINTS.
 
      IF(IRIGHT - ILEFT.EQ.2) THEN
C       SHIFT DATA TO MAKE ROOM FOR POSSIBLE ADDITION OF TWO POINTS
C       WHEN ONLY SPACE FOR ONE EXISTS.
 
        IF(ISTART+NPNT.GT.MAXPNT) THEN
C         NOT ENOUGH ROOM LEFT IN THE VECTORS
          WRITE(STDOUT,50) MAXPNT
          EFLAG = EFLAG + 1
          RETURN
        ENDIF
 
        DO 200 I=ISTART+NPNT-1, IRIGHT,-1
          X(I+1) = X(I)
          Z(I+1) = Z(I)
          SB(I+1) = SB(I)
          LSN(I+1) = LSN(I)
 200    CONTINUE
        IRIGHT = IRIGHT + 1
      ENDIF
 
C     SET INITIAL POINTER TO THE LAST VALID POINT
 
      NP = ILEFT
      IF(EXACTL.EQ.0) THEN
C       TRANSFER THE ADDITIONAL LEFT INTERSECTION POINT
        NP = NP + 1
        X(NP)  = XLEFT
        Z(NP) = ZLEFT
        SB(NP) = SB(NP-1)
        LSN(NP) = NMUD
      ELSE
C       EXACT MATCH-NO ADDITIONAL POINT TO TRANSFER. ADJUST N VALUE.
        LSN(NP) = NMUD
      ENDIF
 
      IF(EXACTR.EQ.0) THEN
C       TRANSFER THE ADDITIONAL RIGHT INTERSECTION POINT
        NP = NP + 1
        X(NP) = XRIGHT
        Z(NP) = ZRIGHT
        SB(NP) = SB(NP-1)
        LSN(NP) = NCON
      ENDIF
 
C     TRANSFER THE REMAINING POINTS
      DO 300 I=IRIGHT,ISTART+NPNT-1
        NP = NP + 1
        X(NP) = X(I)
        Z(NP) = Z(I)
        SB(NP) = SB(I)
        LSN(NP) = LSN(I)
 300  CONTINUE
 
C     CHANGE THE NUMBER OF POINTS
 
      NPNT = NP - ISTART + 1
 
C      WRITE(STDOUT,999)
C999   FORMAT('    OFFSET ELEVATION SUB# ROUGHNESS')
C      DO 400 I=ISTART,ISTART+NPNT-1
C        WRITE(STDOUT,1000) X(I), Z(I), SB(I), LSN(I)
C1000  FORMAT(F10.4,F10.3,I5,F10.3)
 400  CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE   SETSUB
     I                   (ISUB, NPNT, ISTART, N,
     O                    SBU, LSN)
 
C     + + + PURPOSE + + +
C     Set the subsection and line segment roughness
C     for a multiple conduit system.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER ISTART, ISUB, NPNT
      INTEGER SBU(*)
      REAL LSN(*), N
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ISUB   - Subsection number to assign
C     NPNT   - Number of points on boundary of a cross section
C     ISTART - Start point for a conduit cross section description
C     N      - Manning's n value
C     SBU    - Subsection numbers for the line segments-upstream
C               location
C     LSN    - Line segment Manning's n value
 
C     + + + LOCAL VARIABLES + + +
      INTEGER J
C***********************************************************************
      DO 400 J=1,NPNT
        SBU(ISTART+J-1) = ISUB
        LSN(ISTART+J-1) = N
 400  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   MKBOX
     I                  (XOFF, ZOFF, SPAN, RISE, WSLOT, HSLOT,
     O                   NPNTS, X, Z)
 
C     + + + PURPOSE + + +
C     Construct a box culvert of width, SPAN, and height, RISE.
C     RISE is measured from the bottom of the culvert.  ZOFF
C     gives the height of the bottom of the culvert above some
C     datum used to reference the other conduits in the multiple
C     conduit system.  Height of slot, HSLOT, is measured from the
C     same datum as is ZOFF.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNTS
      REAL HSLOT, RISE, SPAN, WSLOT, X(*), XOFF, Z(*), ZOFF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     XOFF   - Offset between conduits to prevent overlap
C     ZOFF   - Vertical offset of invert from culvert datum
C     SPAN   - Maximum horizontal extent of a closed conduit
C     RISE   - Maximum vertical extent of a conduit
C     WSLOT  - Width of slot in closed conduit to preserve free surface
C     HSLOT  - Height of slot for a closed conduit
C     NPNTS  - Number of points on boundary of a cross section
C     X      - Offsets of points on cross section boundary
C     Z      - Elevation at points on cross section boundary

C     + + + LOCAL VARIABLES + + +
      INTEGER I, II, N

      REAL H, DIV, PI

C     + + + DATA INITIALIZATIONS + + +
      DATA N/18/, PI/3.14159265/

C***********************************************************************
C     ADD AN EXTRA POINT WHICH IS LATER DELETED TO REMOVE NOISE FROM
C     THE TOP OF THE SLOT
 
      X(1) = -0.5*WSLOT + XOFF
      Z(1) =  HSLOT + 0.01
 
C     TRAVERSE THE PERIMETER OF THE BOX, INCLUDING THE SLOT WHICH
C     IS CONSIDERED TO RISE FROM THE MIDDLE OF THE BOX.  THE TOP
C     SURFACE OF THE BOX IS CONSIDERED TO BE SLOPED SLIGHTLY.
 
      X(2) = -0.5*WSLOT + XOFF
      Z(2) = HSLOT
 
      X(3) = -0.5*WSLOT + XOFF
      Z(3) = RISE + WSLOT + ZOFF

      H = RISE - WSLOT
      DIV = FLOAT(N-1) 
      DO 100 I=1,N
        ii = n - i + 1
        X(I+3) = -0.5*SPAN 
        Z(I+3) = H*(1.0 - SQRT(FLOAT(I-1)/DIV)) + ZOFF
C                    3        2       1    1      1   12    2       1   123
C        Z(I+3) = 0.5*(h -h*cos(pi*real(2*ii-1)/real(2*n))/
C     a       cos(pi/real(2*n))) + ZOFF
C        WRITE(*,'(I5,F10.3)') I, Z(I+3)
100   CONTINUE

C     Flip to the right-hand side to complete the perimeter
      DO 200 I=1,N
        X(I+N+3) = -X(N+4-I)      
        Z(I+N+3) = Z(N+4-I)
        X(I+N+3) = X(I+N+3) + XOFF
        X(N+4-I) = X(N+4-I) + XOFF
200   CONTINUE
      
C     Add last two points to the slot
      X(N+N+4) = 0.5*WSLOT + XOFF
      Z(N+N+4) = Z(3)
      
      X(N+N+5) = 0.5*WSLOT + XOFF
      Z(N+N+5) = Z(2)
      NPNTS = N + N + 5      
 
      RETURN
 
      END
C
C
C
      SUBROUTINE   RHARCH
     I                   (STDOUT, NP, RSVEC, SPVEC, R1VEC, R3VEC, CFAC,
     M                    RISE, SPAN,
     O                    EFLAG, NRH, RHX, RHY, A)
 
C     + + + PURPOSE + + +
C     Compute the points on the right-hand semi-perimeter of an
C     arch pipe.  The nature of the arch is defined in the
C     four vectors that describe the standard parameters for an
C     arch pipe: Rise, Span, radius of the two corners, and the
C     radius of the bottom.  The other paramters are computed from
C     these to force the parameter set to be consistent with the
C     ideal geometry of these pipes.  Adjust the points so that the
C     area of the resulting ploygon approximating the perimeter matches
C     the area computed from the parameters describing the shape of the
C     arch pipe.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, NP, NRH, STDOUT
      REAL A, CFAC, R1VEC(NP), R3VEC(NP), RHX(23), RHY(23), RISE, 
     A     RSVEC(NP), SPAN, SPVEC(NP)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     NP     - number of points on the boundary
C     RSVEC  - vector of conduit rises
C     SPVEC  - vector of conduit spans
C     R1VEC  - vector of invert radi
C     R3VEC  - vector of corner radi
C     CFAC   - conversion factor from source units(feet or meters)
C              to inches, the unit used in the standard tables
C     RISE   - Maximum vertical extent of a conduit
C     SPAN   - Maximum horizontal extent of a closed conduit
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     NRH    - Number of points on the right-hand semi-perimeter
C     RHX    - Right hand semi-perimeter offset values
C     RHY    - Right hand semi-perimeter ordinate values
C     A      - offset of conduit soffit center from invert
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL ALPHA, B, C, DALPHA, FAC, OMEGA, P, PHI, R1, R2, R3, RS, SP,
     A     TA, THETA
      DOUBLE PRECISION EA
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, ASIN, COS, SIN, SQRT
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(' *ERR:571* RISE=',F8.3,' < MIN RISE=',F8.3,
     A      ' OR  > MAX RISE=',F8.3)
 52   FORMAT(' *WRN:527* GIVEN SPAN=',F8.3,' DIFFERS FROM TABLE SPAN=',
     A        F8.3,' BY > 2%')
 54   FORMAT(' *BUG:505* SPAN AND RISE=0.0 IN SUBROUTINE RHARCH.')
 56   FORMAT(' *ERR:572* SPAN=',F8.3,' < MIN SPAN=',F8.3,
     A     ' OR > MAX SPAN=',F8.3)
C***********************************************************************
C     CONVERT RISE AND SPAN TO INCHES FOR TABLE LOOKUP
      RISE = CFAC*RISE
      SPAN = CFAC*SPAN
 
C     SEARCH FOR THE GIVEN RISE IF RISE IS GREATER THAN 0.0
      IF(RISE.GT.0.0) THEN
        IF(RISE.LT.RSVEC(1).AND.RISE.GT.RSVEC(1) - .1) THEN
          RISE = RSVEC(1)
        ENDIF
        IF(RISE.GT.RSVEC(NP).AND.RISE.LT.RSVEC(NP) + .1) THEN
          RISE = RSVEC(NP)
        ENDIF
        IF(RISE.LT.RSVEC(1).OR.RISE.GT.RSVEC(NP)) THEN
          WRITE(STDOUT,50)  RISE/CFAC, RSVEC(1)/CFAC, RSVEC(NP)/CFAC
          EFLAG = 1
          RETURN
        ENDIF
 
        I = 1
 100    CONTINUE
          I = I + 1
          IF(RSVEC(I).LT.RISE) GOTO 100
 
C       I POINTS TO THE UPPER END OF THE INTERVAL CONTAIN THE GIVEN
C       VALUE OF RISE.  INTERPOLATE FOR OTHER VALUES.
 
        P = (RISE - RSVEC(I-1))/(RSVEC(I) - RSVEC(I-1))
        SP = SPVEC(I-1) + P*(SPVEC(I) - SPVEC(I-1))
        IF(SPAN.GT.0.0) THEN
          IF(ABS(SPAN - SP)/SP.GT.0.02) THEN
            WRITE(STDOUT,52) SPAN/CFAC, SP/CFAC
          ENDIF
        ENDIF
        SPAN = SP
 
        R1 = R1VEC(I-1) + P*(R1VEC(I) - R1VEC(I-1))
        R3 = R3VEC(I-1) + P*(R3VEC(I) - R3VEC(I-1))
      ELSEIF(SPAN.GT.0.0) THEN
C       USE THE SPAN TO LOOK UP VALUES AND DEFINE THE RISE.  IF
C       THE RISE IS GIVEN IT IS TAKEN AS THE DEFINING VALUE.
 
        IF(SPAN.LT.SPVEC(1).OR.SPAN.GT.SPVEC(NP)) THEN
          IF(SPAN.LT.SPVEC(1).AND.SPAN.GT.SPVEC(1) - .1) THEN
            SPAN = SPVEC(1)
          ENDIF
          IF(SPAN.GT.SPVEC(NP).AND.SPAN.LT.SPVEC(NP) + .1) THEN
            SPAN = SPVEC(NP)
          ENDIF
          WRITE(STDOUT,56)  SPAN/CFAC, SPVEC(1)/CFAC, SPVEC(NP)/CFAC
          EFLAG = 1
          RETURN
        ENDIF
 
        I = 1
 105    CONTINUE
          I = I + 1
          IF(SPVEC(I).LT.SPAN) GOTO 105
 
C       I POINTS TO THE UPPER END OF THE INTERVAL CONTAIN THE GIVEN
C       VALUE OF RISE.  INTERPOLATE FOR OTHER VALUES.
 
        P = (SPAN - SPVEC(I-1))/(SPVEC(I) - SPVEC(I-1))
 
        RS = RSVEC(I-1) + P*(RSVEC(I) - RSVEC(I-1))
        RISE = RS
        R1 = R1VEC(I-1) + P*(R1VEC(I) - R1VEC(I-1))
        R3 = R3VEC(I-1) + P*(R3VEC(I) - R3VEC(I-1))
      ELSE
        WRITE(STDOUT,54)
        STOP 'Abnormal stop. Errors found.'
      ENDIF
 
C     CONVERT VALUES TO FEET or meters
 
      RISE = RISE/CFAC
      SPAN = SPAN/CFAC
      R1 = R1/CFAC
      R3 = R3/CFAC
 
C     COMPUTE CONDUIT PARAMETERS
 
      C = SPAN/2.0 - R3
      B = R1 - SQRT((R1 -R3)**2 - C**2)
      R2 = (R3**2 - C**2 - (RISE - B)**2)/(2.*(B -RISE + R3))
      A = RISE - R2
 
 
      THETA = ASIN(C/(R2 - R3))
      IF(A.GT.B) THEN
        FAC = -1.0
        THETA = 3.1415927 - THETA
      ELSE
        FAC = 1.0
      ENDIF
 
      OMEGA = ASIN(C/(R1 -R3))
      PHI = 3.14593 - THETA - OMEGA
C     COMPUTE THE TRUE AREA OF THE SEMI-PERIMETER
 
      TA = 0.5*( THETA*R2**2 - FAC*C*SQRT((R2 -R3)**2 - C**2) +
     A     C*(2.*R1 -R3)*(R1 - B)*R3/(R1 -R3)**2 +
     B     PHI*R3**2 + R1**2*(OMEGA - C*(R1 - B)/(R1 -R3)**2))
 
 
C     NOW COMPUTE THE POINTS ON THE PERIMETER
 
C     FIRST SECTION- 9 SEGMENTS
 
      DO 110 I=1,10
        ALPHA = (I - 1)*THETA/9.0
        RHX(I) = R2*SIN(ALPHA)
        RHY(I) = R2*COS(ALPHA)
 110  CONTINUE
 
C     SECOND SECTION- 8 SEGMENTS
 
      DALPHA = PHI/8.0
      ALPHA = 1.570796 - THETA
      DO 120 I=11,18
        ALPHA = ALPHA - DALPHA
        RHX(I) = C + R3*COS(ALPHA)
        RHY(I) = (B - A) + R3*SIN(ALPHA)
 120  CONTINUE
 
C     THIRD SECTION- 5 SEGMENTS
 
      DO 130 I=19,23
        ALPHA = (23 - I)*OMEGA/5.0
        RHX(I) = R1*SIN(ALPHA)
        RHY(I) = R1*(1.0 - COS(ALPHA)) - A
 130  CONTINUE
 
      NRH = 23
 
C     COMPUTE POLYGONAL AREA AND ADJUST THE POINTS TO MATCH TRUE AREA
C     EXACTLY
 
      EA = 0.D0
      DO 140 I=2,NRH
        EA = EA + 0.5*(RHY(I) + RHY(I-1))*(RHX(I) - RHX(I-1))
 140  CONTINUE
 
      FAC = SQRT(TA/EA)
      DO 150 I=1,NRH
        RHX(I) = FAC*RHX(I)
        RHY(I) = FAC*RHY(I)
 150  CONTINUE
 
      RETURN
      END
C
C
C
      SUBROUTINE   RHMAK
     I                  (WSLOT, HSLOT, YOFF, XOFF, NRH, RHX, RHY,
     O                   X, Y, NPNTS)
 
C     + + + PURPOSE + + +
C     Construct a closed conduit shape given the right half of the
C     shape, i. e. the right hand semi-perimeter.  Use
C     methods similar to those in URQMAK.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNTS, NRH
      REAL HSLOT, RHX(*), RHY(*), WSLOT, X(*), XOFF, Y(*), YOFF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     WSLOT  - Width of slot in closed conduit to preserve free surface
C     HSLOT  - Height of slot for a closed conduit
C     YOFF   - Vertical offset to construct complete shape from
C               the right semi-perimeter
C     XOFF   - Offset for defining the conduit shape from the
C               right hand semi-perimeter
C     NRH    - Number of points on the right-hand semi-perimeter
C     RHX    - Right hand semi-perimeter offset values
C     RHY    - Right hand semi-perimeter ordinate values
C     X      - Offsets of points on cross section boundary
C     Y      - ordinate values on boundary of conduit
C     NPNTS  - Number of points on boundary of a cross section
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'stdun.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL XR, YR
C***********************************************************************
C     RHX(1) must = 0!
 
      IF(RHX(2).EQ.0.0) THEN
        WRITE(STD6,*) ' ZERO DIVIDE 1 IN RHMAK'
        STOP 'Abnormal stop. Errors found.'
      ENDIF
      XR = WSLOT/2.0
      YR = RHY(1) + (RHY(2) - RHY(1))*XR/RHX(2)
 
      X(1) = -XR + XOFF
      Y(1) = HSLOT + 0.01
      X(2) = -XR + XOFF
      Y(2) = HSLOT
      X(3) = -XR + XOFF
      Y(3) = YR + YOFF
 
C     MOVE AROUND THE SHAPE COUNTER CLOCKWISE
 
      J = 4
      DO 100  I=2,NRH
        X(J) = -RHX(I) + XOFF
        Y(J) = RHY(I) + YOFF
        J = J + 1
 100  CONTINUE
 
      DO 110 I=NRH-1,2,-1
        X(J) = RHX(I) + XOFF
        Y(J) = RHY(I) + YOFF
        J = J + 1
 110  CONTINUE
 
      X(J) = XR + XOFF
      Y(J) = YR + YOFF
      X(J+1) = XR + XOFF
      Y(J+1) = HSLOT
 
      NPNTS = J + 1
      RETURN
      END
C
C
C
      SUBROUTINE   URQMAK
     I                   (WSLOT, HSLOT, YOFF, XOFF, NURQ, URQX, URQY,
     O                    X, Y, NPNTS)
 
C     + + + PURPOSE + + +
C     Construct a closed conduit shape given the upper right quadrant,
C     the slot width, etc describing the conduit and its location
C     relative to other conduits and relative to the datum for the
C     collection of conduits.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NPNTS, NURQ
      REAL HSLOT, URQX(*), URQY(*), WSLOT, X(*), XOFF, Y(*), YOFF
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     WSLOT  - Width of slot in closed conduit to preserve free surface
C     HSLOT  - Height of slot for a closed conduit
C     YOFF   - Vertical offset to construct complete shape from
C               the upper right quadrant
C     XOFF   - Offset for defining the conduit shape from the upper
C               right quadrant
C     NURQ   - Number of points in the upper right quadrant
C     URQX   - Offsets for upper right quadrant of a closed conduit
C               shape
C     URQY   - Ordinates for upper right quadrant of a closed conduit
C               shape
C     X      - Offsets of points on cross section boundary
C     Y      - ordinate values on boundary of conduit
C     NPNTS  - Number of points on boundary of a cross section
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J
      REAL XR, YR
C***********************************************************************
C     NOTE: URQX(1) must = 0!
 
      XR = WSLOT/2.0
 
      YR = URQY(1) + (URQY(2) - URQY(1))*XR/URQX(2)
 
      X(1) = -XR + XOFF
      Y(1) = HSLOT + 0.01
      X(2) = -XR + XOFF
      Y(2) = HSLOT
      X(3) = -XR + XOFF
      Y(3) = YR + YOFF
 
C     MOVE AROUND THE SHAPE COUNTER CLOCKWISE
 
      J = 4
      DO 100 I=2,NURQ
        X(J) = -URQX(I) + XOFF
        Y(J) = URQY(I) + YOFF
        J = J + 1
 100  CONTINUE
 
      DO 110 I=NURQ-1,1,-1
        X(J) = -URQX(I) + XOFF
        Y(J) = -URQY(I) + YOFF
        J = J + 1
 110  CONTINUE
 
      DO 120 I=2,NURQ
        X(J) = URQX(I) + XOFF
        Y(J) = -URQY(I) + YOFF
        J = J + 1
 120  CONTINUE
 
      DO 130  I=NURQ-1,2,-1
        X(J) = URQX(I) + XOFF
        Y(J) = URQY(I) + YOFF
        J = J + 1
 130  CONTINUE
 
      X(J) = XR + XOFF
      Y(J) = YR + YOFF
      X(J+1) = XR + XOFF
      Y(J+1) = HSLOT
 
      NPNTS = J + 1
      RETURN
      END
C
C
C
      SUBROUTINE   URQTE
     I                  (RISE, SPAN,
     O                   NURQ, URQX, URQY)
 
C     + + + PURPOSE + + +
C     Construct the upper right quadrant from a true elliptical shape.
C     RISE and SPAN are any postive numbers.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER NURQ
      REAL RISE, SPAN, URQX(*), URQY(*)
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     RISE   - Maximum vertical extent of a conduit
C     SPAN   - Maximum horizontal extent of a closed conduit
C     NURQ   - Number of points in the upper right quadrant
C     URQX   - Offsets for upper right quadrant of a closed conduit
C               shape
C     URQY   - Ordinates for upper right quadrant of a closed conduit
C               shape
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I
      REAL A, B, FAC, R, TA, THETA
      DOUBLE PRECISION EA
 
C     + + + INTRINSICS + + +
      INTRINSIC COS, SIN, SQRT
C***********************************************************************
      A = SPAN/2.0
      B = RISE/2.0
 
      NURQ = 11
      DO 100 I=1,NURQ
        THETA = 1.570796*(NURQ - I)/(NURQ - 1)
        R = SQRT( (A*B)**2/( (A*SIN(THETA))**2 + (B*COS(THETA))**2))
        URQX(I) = R*COS(THETA)
        URQY(I) = R*SIN(THETA)
 100  CONTINUE
 
C     COMPUTE THE TRUE AREA AND THE POLYGONAL AREA AND ADJUST
 
      TA = RISE*SPAN*3.141593/4.0
 
      EA = 0.D0
      DO 110 I=2,NURQ
        EA = EA + 0.5*(URQY(I) + URQY(I-1))*(URQX(I) - URQX(I-1))
 110  CONTINUE
 
      FAC = SQRT(TA/(4.*EA))
 
      DO 120 I=1,NURQ
        URQX(I) = FAC*URQX(I)
        URQY(I) = FAC*URQY(I)
 120  CONTINUE
 
      RETURN
      END
