c     Routines that relate directly to handling zone, hgrid, vdatum, unitsys,
c     easting, and northing, that is, items that define the location of 
c     a function table or point in space.

c     These routines will be called (sooner or later) by the code for the
c     the following commands: MULPIPES, MULCON, SEWER, AXIALPUMP, PUMPLOSS,
c     CULVERT, ORIFICE, EXPCON, CRITQ, EMBANKQ,
c     RISERCLV and maybe more:)
 

C
C
C 
      SUBROUTINE SET_lctn_ITEM_DEFAULTS()

C     Set the default values in the vectors used to establish location
      IMPLICIT NONE
      INCLUDE 'lctnitem.cmn'

C***********************************************************************
C     Default for: ZONE 
      lctnITMCTAB(  1) = 'NONE'           
C     Default for: HGRID
      lctnITMCTAB(  2) = 'NONE'           
C     Default for: VDATUM
      lctnITMCTAB(  3) = 'NONE'
C     Default for: UNITSYS
      lctnITMCTAB(  4) = 'NONE'
C     Default for: BASIS
      lctnITMCTAB(  5) = 'NONE'
c     We set the default for the following items to a value thaat
c     is negative and smaller than any hgrid value that might make
c     sense.  Of course I am assuming that the person defining an hgrid
c     uses either feet or meters.  If they choose some stupid unit,
c     that is too small, then they deserve to get stupid results:)
c     That value is that for 25000/4=6250 miles*5280 ft/mile = 33d6.  Which
c     is slightly larger than one fourth the circumference of the 
c     earth. 
C     Default for: EASTING
      lctnITMDTAB(  1) = -33d6
C     Default for: NORTHING
      lctnITMDTAB(  2) = -33d6

      RETURN
      END
C
C
C
      SUBROUTINE  SET_lctn_ITEMS(
     O             zone, hgrid, vdatum, unitsys, basis, easting,
     O             northing)

C     Set items for location for various commands, many in fact:)
C     All values not set explicitly by user are at their default value.

      IMPLICIT NONE

      real*8 easting, northing
      character*8 zone, hgrid, vdatum, unitsys, basis
    


      INCLUDE 'lctnitem.cmn'
      INCLUDE 'stdun.cmn'
C***********************************************************************
C     Set the value for ZONE
      zone =    lctnITMCTAB( 1)
C     Set the value for HGRID
      hgrid =   lctnITMCTAB( 2)
C     Set the value for VDATUM
      vdatum =  lctnITMCTAB( 3)   
C     Set the value for UNITSYS
      unitsys = lctnITMCTAB( 4)
C     Set the value for BASIS
      basis = lctnITMCTAB( 5)
C     Set the value for EASTING
      easting =  lctnITMDTAB( 1)
C     Set the value for NORTHING
      northing = lctnITMDTAB( 2)

      call chk_vdatum_unitsys(std6, vdatum, unitsys, 
     a          ' during command input')

    
      RETURN
      END

C
C
C
      SUBROUTINE GET_lctn_ITEMS(STDIN, STDOUT,
     M                   EFLAG)

C     Get the location items from those commands that do 
c     not easily support them now.  

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG

      INCLUDE 'lctnitem.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  INTVAL, REAVAL, NONE,
     A         CHRVAL, DPRVAL, EXACT, LOWER,
     B         NUMERIC, CHAR, N_SYMBOL, NXTBLK
      PARAMETER(N_SYMBOL=22, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, NXTBLK=2, NONE=0,
     B          EXACT=0,LOWER=1, NUMERIC=0, CHAR=1)

      INTEGER MAX_LINE
     A        

      EXTERNAL GET_NAMED_ITEMS, SET_lctn_ITEM_DEFAULTS

C     + + + SAVED VALUES + + +
      INTEGER GROUP(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL), GROUP_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE  SYMBOL_TABLE, GROUP, RESPONSE_TYPE, CONVERT_RULE,
     A      GROUP_INDEX

      DATA  SYMBOL_TABLE /
     A 'ZONE','HGRID','VDATUM','UNITSYS',
     b 'EASTING','NORTHING','NSIDES','WSLOT','DIAMETER','QUNIT',
     C 'QUNIT','MODE','FIT_WITH','LABEL','CROSS','APPTAB','PKCWTB',
     d 'PLCWTB','BASIS','SHIFT','WIDTH','HDATUM'/
                                                                      
      DATA GROUP  /
     * 4*CHAR, 2*NUMERIC,12*NXTBLK,CHAR,3*NXTBLK/          

      DATA GROUP_INDEX /
     * 1, 2, 3, 4, 1, 3, 12*0,5,3*0/
                                       
      DATA RESPONSE_TYPE  /
     * 4*CHRVAL,2*DPRVAL,12*NONE,CHRVAL,3*NONE/
    
      DATA CONVERT_RULE /
     * 4*LOWER,2*LOWER,12*EXACT,LOWER,3*EXACT/
   
C***********************************************************************
C     Set Defaults
 
      CALL SET_lctn_ITEM_DEFAULTS()

      MAX_LINE = 6
      CALL GET_NAMED_ITEMS(
     I                     STDIN, STDOUT,  MAX_LINE, N_SYMBOL,
     I  GROUP, RESPONSE_TYPE, CONVERT_RULE, GROUP_INDEX, SYMBOL_TABLE,
     I  MAXR_lctnITM, MAXDP_lctnITM, MAXC_lctnITM, 'Location items',
     O  lctnITMITAB, lctnITMFTAB, lctnITMDTAB, lctnITMCTAB,
     O  EFLAG)
      
      RETURN

      END

