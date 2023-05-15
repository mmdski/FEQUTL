C
C
C
      REAL FUNCTION FIND_Y4_FROM_Y43_RESID(
     I                                       Y4)

C     Find conditions at section 4 given a depth 
C     section 3

      IMPLICIT NONE
      REAL Y4

      INCLUDE 'y43_to_4.cmn'
C     Local

      REAL  M4, A4, T4, DT4, J4, K4, DK4, B4, DB4

c     common

C***********************************************************************
C      WRITE(STDOUTA,*) ' Y4=',Y4,' QSQR=',QSQR
      CALL XLKT21
     I           (TAB4,
     M            Y4,
     O            A4, T4, DT4, J4, K4, DK4, B4, DB4)

      M4 = (B4*QSQR/A4) + G*J4
     
C      WRITE(STDOUTA,*) ' M4=',M4,' MEXIT=',MEXIT
      FIND_Y4_FROM_Y43_RESID = (M4 - MEXIT)/MEXIT
      RETURN
      END
          
C
C
C
      SUBROUTINE FIND_Y4_FROM_Y43(
     I            STDOUT, GRAV, Q, DEPTAB,
     I            Z3B, Z43B, Y3, AFLUX, COSTHETA,
     O            Y4, FLAG )  

C     Find the value of depth at Y3 given the flow and depth
C     in the departure section.  

      IMPLICIT NONE

      INTEGER FLAG, STDOUT, DEPTAB

      REAL GRAV, Q, Z3B, Z43B,  Y3,
     A     COSTHETA, Y4, AFLUX

      INCLUDE 'y43_to_4.cmn'
      INCLUDE 'epscom.cmn'

C     Functions

      REAL FIND_Y4_FROM_Y43_RESID
      EXTERNAL FIND_Y4_FROM_Y43_RESID

C     Local

      INTEGER I
      REAL YL, FL, YR, FR, Y43,
     A     A43, T43, DT43, J43, K43, DK43, B43, DB43


C***********************************************************************
      QSQR = Q*Q
      G = GRAV

      Y43 = Y3 + Z3B - Z43B
      CALL XLKT21
     I           (DEPTAB,
     M            Y43,
     O            A43, T43, DT43, J43, K43, DK43, B43, DB43)

      MEXIT = COSTHETA*QSQR/AFLUX + G*J43

      TAB4 = DEPTAB
      STDOUTA = STDOUT
      
C     Search for a change in sign of the residual.

      I = 50
      YL = Y3 + Z3B - Z43B
      YR = -1.0
 100  CONTINUE
        FL = FIND_Y4_FROM_Y43_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        I = I - 1
        IF(I.EQ.0) THEN
          WRITE(STDOUT,*) ' FIND_Y4_FROM_Y43: No neg. resid.'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
        IF(FL.LE.0.0D0) THEN
          GOTO 110
        ELSE
          YR = YL
          FR = FL
          YL = 0.9*YL
          GOTO 100
        ENDIF
 110  CONTINUE
      IF(YR.LT.0.0D0) THEN
        YR = YL
 120    CONTINUE
          YR = 1.1*YR
          FR = FIND_Y4_FROM_Y43_RESID(YR)
          I = I - 1
          IF(I.EQ.0) THEN
            WRITE(STDOUT,*) ' FIND_Y4_FROM_Y43: No pos. resid.'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
          IF(FR.GE.0.0D0) THEN
            GOTO 130
          ELSE
            FL = FR
            YL = YR
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF

C     Found an interval that contains a root. 
C      WRITE(STDOUT,*) ' Calling REGFLT with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL REGFLT
     I           (EPSARG, EPSF, FIND_Y4_FROM_Y43_RESID,
     M            YL, YR, FL, FR,
     O            Y4, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from REGFLT with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' Y3 =',Y3
      IF(FLAG.EQ.3)  FLAG = 0

      RETURN
      END
C
C
C
      REAL*8 FUNCTION FIND_SEQUENT_DEPTH_RESID(
     I                                         YB)

C     Find the sequent depth for a super-critical flow in
C     conduit.  Ignores friction and gravity. 

      IMPLICIT NONE
      REAL*8 YB

      INCLUDE 'barrel.cmn'

      INCLUDE 'seqd.cmn'


C     Local

      REAL*8 AB, JB, MB

C***********************************************************************
      CALL LKT_JDA(
     I             0, XLOC, YB,
     O             AB, JB)

      MB =  QSQR/AB + G*JB*COS_THETA
      FIND_SEQUENT_DEPTH_RESID = (MB - MVC)/MVC
      RETURN
      END
C
C
C
      SUBROUTINE FIND_SEQUENT_DEPTH(
     I                              STDOUT, GRAV, HG, BG, CC, Q,
     I                              EPSARG, EPSF,
     O                              YBA)

C     Find depth sequent to a super-critical flow

      IMPLICIT NONE
      INTEGER STDOUT

      REAL*8 HG, BG, CC, Q, YBA, GRAV,
     A       EPSARG, EPSF

      INCLUDE 'barrel.cmn'
      INCLUDE 'seqd.cmn'

C     Subroutines and functions
      REAL*8 FIND_SEQUENT_DEPTH_RESID
      EXTERNAL FIND_SEQUENT_DEPTH_RESID

C     Local variables

      INTEGER I, FLAG

      REAL*8 YL, FL, YR, FR, FRDN2, YVC, AVC, JVC, QCVC

C***********************************************************************
      G = GRAV
      QSQR = Q*Q
      XLOC = STATION(1) + HG
      YVC = CC*HG
      CALL LKT_JDA(
     I             1, XLOC, YVC,
     O             AVC, JVC)

      MVC = QSQR/(YVC*BG) + G*JVC*COS_THETA
      CALL LKT_QCD(
     I           0, XLOC, YVC,
     O           QCVC)

      FRDN2 = QSQR/QCVC**2

      YL = 0.5D0*YVC*(SQRT(1.D0 + 8.0D0*FRDN2) - 1.D0)
      YR = -1.D0
      I = 100
 100  CONTINUE
        FL = FIND_SEQUENT_DEPTH_RESID(YL)
        I = I - 1
        IF(I.EQ.0) THEN
          WRITE(STDOUT,*) ' FIND_SEQUENT_DEPTH: No negative residual.'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        IF(FL.LE.0.0D0) THEN
          GOTO 110
        ELSE
          YR = YL
          FR = FL
          YL = 0.95D0*YL
          GOTO 100
        ENDIF
 110  CONTINUE
  
      IF(YR.LT.0.0D0) THEN
        YR = YL
 120    CONTINUE
          YR = 1.1*YR
          FR = FIND_SEQUENT_DEPTH_RESID(YR)
          I = I - 1
          IF(I.EQ.0) THEN
        WRITE(STDOUT,*) ' FIND_SEQUENT_DEPTH: No positive residual.'
            STOP 'Abnormal stop. Errors found.'
          ENDIF
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
          IF(FR.GE.0.0D0) THEN
            GOTO 130
          ELSE
            FL = FR
            YL = YR
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF


C     Found an interval that contains a root. 
C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG, EPSF, FIND_SEQUENT_DEPTH_RESID,
     M             YL, YR, FL, FR,
     O             YBA, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' YBA =',YBA, ' KNT=',INT(FR)
      IF(FLAG.EQ.3)  FLAG = 0

      RETURN
      END
            
C     
C
C
C
      REAL FUNCTION FIND_SO_FLOW_RESID(
     A                                 YVC)

C     Find the SO flow given tailwater at section 4 and the
C     headwater at section 1.  We use the depth at the vena
C     contracta to try to make the solution faster.

      IMPLICIT NONE
      REAL YVC

      INCLUDE 'so_flow.cmn'
      INCLUDE 'barrel.cmn'

C     Functions and subroutines
      REAL*8 E, RHSE
      EXTERNAL E, RHSE, FIND_YCYNYM

C     Local
      INTEGER NS, NE
      REAL*8 YS, YB, AB, JB, ABVC, JBVC, ML, MR,
     A       YCRIT, YNORM, YMAX, EPSARG, EPSF, EPSABS


C***********************************************************************
      EPSARG = 0.01D0*EPSD
      EPSF = 0.00001D0*EPSD
      EPSABS = 0.0001D0*EPSD

C     Compute an estimate of the flow 
C     Compute the flow from the given YVC and the energy equation. 
C     Note that only valid values of YVC should appear here.

C      WRITE(STDOUTA,*) ' FIND_SO_FLOW_RESID: YVC=',YVC
C      WRITE(STDOUTA,*) ' AFAC=',AFAC, ' H1EFF=',H1EFF
C      WRITE(STDOUTA,*) ' DIV=',DIV,' COSTHETA=',COSTHETA
      
      QSQR = AFAC*AFAC*((G+G)*(H1EFF - YVC*COSTHETA)/DIV)

C      WRITE(STDOUTA,*) ' FIND_SO_FLOW_RESID: Q=',SQRT(QSQR)

C     Now compute the conditions at the culvert exit from
C     the tailwater level and the just computed flow.

      CALL FIND_Y43_FROM_Y4(
     I            STDOUTA, G, QSQR, DTAB, ETAB,
     I            Z3B, Z4B, Y4, A4, B4, J4, COSTHETA,
     O            Y3, RFLAG )  


C     Next, compute water surface profile from barrel exit 
C     to the vena contracta in the culvert barrel.

      QD = SQRT(DBLE(QSQR))

C     Find critical depth and normal depth
      YCRIT = Y3
      YNORM = Y3
      CALL FIND_YCYNYM
     I                (STDOUTA, XS, EPSARG, EPSF, EPSABS, 
     M                 YCRIT, YNORM,
     O                 YMAX, RFLAG)

      IF(Y3.LE.YCRIT) THEN
        YS = 1.001*YCRIT
        YC_FLAG = 1
C        WRITE(STDOUTA,*) ' Invalid Y3 in FIND_SO_FLOW_RESID.'
C        WRITE(STDOUTA,*) ' Y3=',Y3,' YCRIT=',YCRIT
C        WRITE(STDOUTA,*) ' QD=',QD
      ELSE
       YC_FLAG = 0
       YS = DBLE(Y3)
      ENDIF
      CALL SFWSP(
     I           STDOUTA, XS, XE, YS, EPSD,
     I           E, RHSE, FIND_YCYNYM, NMAX,
     O           NS, NE, XVEC, YVEC, YCVEC, YNVEC, RFLAG)

C     Check to see if the profile reached the goal
      IF(ABS(XE - XVEC(NS))/ABS(XS - XE).GT.1.D-6) THEN
        RFLAG = -9
        FIND_SO_FLOW_RESID= 0.0
        XE = XVEC(NS)
      ELSE

C       Compute the residual from the momentum balance at the 
C       vena contracta

        YB = YVEC(NS)

        CALL LKT_JDA(
     I               0, XE, YB,
     O               AB, JB)
            
        CALL LKT_JDA(
     I               0, XE, DBLE(YVC),
     O               ABVC, JBVC)

        ML = QSQR/AJET + G*JBVC*COSTHETA
        MR = QSQR/AB + G*JB*COSTHETA

        FIND_SO_FLOW_RESID = (MR - ML)/ML
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE FIND_SO_FLOW(
     I    STDOUT, DEPTAB, EXITTAB, GRAV, H1, HG, BG, CC, CD,
     I    SIN_THETA, COS_THETA, Z3BA, Z4BA, Y4A, A4A, B4A, J4A,
     I    ALPHA1, A1,
     I    NXS, STATION, EPS_D,
     M    YVCOLD,
     O    Q, Y3A, XE_R, FLAG)

C     Compute the SO flows 

      IMPLICIT NONE
      INTEGER STDOUT, DEPTAB, EXITTAB, FLAG, NXS

      REAL GRAV, H1, HG, BG, CC, CD, SIN_THETA, COS_THETA, 
     A     Z3BA, Z4BA, Y4A, A4A, B4A, J4A, Q, ALPHA1, A1,
     B     YVCOLD, Y3A, XE_R

      REAL*8 STATION(NXS), EPS_D

      INCLUDE 'so_flow.cmn'
      INCLUDE 'epscom.cmn'

C     Functions and subroutines
      REAL FIND_SO_FLOW_RESID
      EXTERNAL FIND_SO_FLOW_RESID       

C     Local

      INTEGER I
      REAL YL, FL, YR, FR, YVC_MIN, YVC_MAX, YVC
C***********************************************************************

      FLAG = 0
C      WRITE(STDOUT,*) ' FIND_SO_FLOW: H1=',H1,' Y4=',Y4A
C     Set values in the common block for the residual function
      DTAB = DEPTAB
      ETAB = EXITTAB
      STDOUTA = STDOUT
      RFLAG = 0
      G = GRAV
      
      AJET = DBLE(HG)*DBLE(BG)*DBLE(CC)
      AFAC = HG*BG*CC*CD
      DIV = 1.0 - ALPHA1*(AFAC/A1)**2
      H1EFF = H1 + HG*SIN_THETA
      COSTHETA = COS_THETA

      Z3B = Z3BA
      Z4B = Z4BA
      Y4 = Y4A
      A4 = A4A
      B4 = B4A
      J4 = J4A
      XS = STATION(NXS)
      XE = STATION(1) + HG
      EPSD = EPS_D

      YVC_MIN = YVCOLD
      YVC_MAX = 0.9999*H1EFF/COS_THETA

C      WRITE(STDOUT,*) ' YVC_MIN=',YVC_MIN,' YVC_MAX=',YVC_MAX
      YL = YVC_MIN
      YR = -1.0
      I = 100
 100  CONTINUE
        FL = FIND_SO_FLOW_RESID(YL)
        I = I - 1
        IF(I.EQ.0) THEN
          WRITE(STDOUT,*) ' FIND_SO_FLOW: No negative residual.'
          STOP 'Abnormal stop. Errors found.'
        ENDIF
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        IF(FL.LE.0.0D0) THEN
          GOTO 110
        ELSE
          YR = YL
          FR = FL
          YL = 0.9*YL + 0.1*YVC_MAX
          GOTO 100
        ENDIF
 110  CONTINUE
  
      IF(YR.LT.0.0D0) THEN
        YR = YL
 120    CONTINUE
          YR = 0.9*YR + 0.1*YVC_MAX
          FR = FIND_SO_FLOW_RESID(YR)
          I = I - 1
          IF(I.EQ.0) THEN
C            WRITE(STDOUT,*) ' FIND_SO_FLOW: No positive residual.'
C            WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
C            WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
             FLAG = -10
           RETURN

          ENDIF
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
          IF(FR.GE.0.0D0) THEN
            GOTO 130
          ELSE
            FL = FR
            YL = YR
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF

      IF(RFLAG.EQ.-9) THEN
        FLAG = -11
        XE_R = XE
        RETURN
      ENDIF

C     Found an interval that contains a root. 
C      WRITE(STDOUT,*) ' Calling RGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL RGF
     I        (EPSARG, EPSF, FIND_SO_FLOW_RESID,
     M         YL, YR, FL, FR,
     O         YVC, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from RGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' Y3 =',Y3
      IF(FLAG.EQ.3)  FLAG = 0

      YVCOLD = YVC
      Q = SQRT(QSQR)
      Y3A = Y3
      IF(YC_FLAG.EQ.1) FLAG = -12

      RETURN
      END
        
C
C
C
      REAL FUNCTION FIND_Y43_FROM_Y4_RESID(
     I                                       Y3)

C     Find conditions at culvert exist given a depth in the
C     departure section and with the flow known. 

      IMPLICIT NONE
      REAL Y3

      INCLUDE 'y4_to_43.cmn'
C     Local

      REAL  J43, MEXIT, Y43, A3


C***********************************************************************
C     Find first moment at section 43
      Y43 = Y3 + ZB3 - ZB43
      CALL LKTJ
     I         (TAB43,
     M          Y43,
     O          J43)

C     Find area in barrel exit.
      CALL LKTA
     I          (TABEXIT,
     M           Y3,
     O           A3)

      MEXIT = QSQR/A3*COSTHETA + G*J43
     
      FIND_Y43_FROM_Y4_RESID = (MEXIT - M4)/M4
      RETURN
      END
          
C
C
C
      SUBROUTINE FIND_Y43_FROM_Y4(
     I            STDOUT, GRAV, Q2, DEPTAB, EXITTAB,
     I            Z3B, Z4B, Y4, A4, B4, J4, COS_THETA,
     O            Y3, FLAG )  

C     Find the value of depth at Y3 given the flow and depth
C     in the departure section.  

      IMPLICIT NONE

      INTEGER FLAG, STDOUT, EXITTAB, DEPTAB

      REAL GRAV, Q2, Z3B, Z4B, A4, B4, J4, Y3,
     A     COS_THETA, Y4

      INCLUDE 'y4_to_43.cmn'
      INCLUDE 'epscom.cmn'

C     Functions

      REAL FIND_Y43_FROM_Y4_RESID
      EXTERNAL FIND_Y43_FROM_Y4_RESID

C     Local

      REAL YL, FL, YR, FR


C***********************************************************************
      FLAG = 0 
      QSQR = Q2
      G = GRAV
      M4 = B4*QSQR/A4 + G*J4

      TAB43 = DEPTAB
      TABEXIT = EXITTAB
      ZB43 = Z4B
      ZB3 = Z3B
      COSTHETA = COS_THETA

C     Search for a change in sign of the residual.

      YL = Y4
      YR = -1.D0

 100  CONTINUE
        FL = FIND_Y43_FROM_Y4_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        IF(FL.LE.0.0D0) THEN
          GOTO 110
        ELSE
          YR = YL
          FR = FL
          YL = 0.98*YL
          IF(YL.LT.0.01) THEN
            WRITE(STDOUT,*) ' NO NEG RESIDUAL'
            Y3 = 0.D0
            FLAG = 4
            RETURN
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
      IF(YR.LT.0.0) THEN
        IF(ABS(FL).LE.EPSF) THEN
          Y3 = YL
          RETURN
        ENDIF
        WRITE(STDOUT,*) 'Problem in FIND_Y43_FROM_Y4'
        WRITE(STDOUT,*) 'No positive residual found'
        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        STOP 'Abnormal stop. Errors found.'
      ENDIF

C     Found an interval that contains a root. 
C      WRITE(STDOUT,*) ' Calling REGFLT with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL REGFLT
     I           (EPSARG, EPSF, FIND_Y43_FROM_Y4_RESID,
     M            YL, YR, FL, FR,
     O            Y3, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from REGFLT with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' Y3 =',Y3
      IF(FLAG.EQ.3)  FLAG = 0

      RETURN
      END




C
C
C
      REAL*8 FUNCTION FIND_Y4SW_RESID(
     I                                Y4SW)

C     Residual function for finding the level at section 4 that
C     causes submerged weir flow to contact the gate lip.

      IMPLICIT NONE

      REAL*8 Y4SW

      INCLUDE 'y4sw.cmn'

C     Local

      INTEGER FREE

      REAL Z4, Y2, DQED, DQEU

C***********************************************************************

      Z4 = Y4SW + Z4B

C      WRITE(STD6,*) ' Z4=',Z4,' Z1=',Z1
C      WRITE(STD6,*) ' GETY2=',GETY2
      CALL TDLK13                                                
     I           (STD6, GETY2, 13, 0, 0.0, Z4, Z1,
     I            HDATUM,
     O            Y2, DQED, DQEU, FREE)

      FIND_Y4SW_RESID = (HG - Y2)/HG
      RETURN
      END

C
C
C
      SUBROUTINE FIND_Y4SW(
     I                     STDOUT, ZUP, ZB2, ZB4, TABY2, HG_D,
     I                     EPSARG, EPSF,
     O                     Y4SW)

C     Find the level at section 4 that causes submerged weir flow
C     to contact the gate lip.

      IMPLICIT NONE

      INTEGER STDOUT, TABY2

      REAL ZUP, ZB2, ZB4
      REAL*8 Y4SW, EPSARG, EPSF

      INCLUDE 'y4sw.cmn'

      REAL*8 HG_D

C     Functions
      REAL*8 FIND_Y4SW_RESID
      EXTERNAL FIND_Y4SW_RESID

C     Local

      INTEGER FLAG

      REAL*8 YL, FL, YR, FR

C***********************************************************************
C     Set the common block values for the residual function

      HG = HG_D
      STD6 = STDOUT
      GETY2 = TABY2
      HDATUM = ZB2
      Z4B = ZB4
      Z1 = ZUP

C     The elevation at section 4 will lie between ZUP and
C     and an unknown increment above Z4B.  Start at the upper
C     elevation and  search for a positive residual. 

      YL = ZUP - Z4B
      YR = -1.D0

 100  CONTINUE
        FL = FIND_Y4SW_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
        IF(FL.LE.0.0D0) THEN
          GOTO 110
        ELSE
          YR = YL
          FR = FL
          YL = 0.9D0*YL
          IF(YL.LT.0.01) THEN
            WRITE(STDOUT,*) ' NO NEG RESIDUAL'
            Y4SW = 0.D0
            FLAG = 4
            RETURN
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
C     Negative residual found; search for positive residual.
      
      IF(YR.LT.0.0D0) THEN
        YR = YL
 120    CONTINUE
          YR = 0.9D0*YR
          IF(YR.LT.0.01) THEN
            WRITE(STDOUT,*) ' NO POS RESIDUAL'
            FLAG = 4
            Y4SW = 0.D0
            RETURN
          ENDIF
          FR = FIND_Y4SW_RESID(YR)
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
          IF(FR.GE.0.0D0) THEN
            GOTO 130
          ELSE
            FL = FR
            YL = YR
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF

C     Found an interval that contains a root. 
C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG, EPSF, FIND_Y4SW_RESID,
     M             YL, YR, FL, FR,
     O             Y4SW, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' Y4SW =',Y4SW, ' KNT=',INT(FR)
      IF(FLAG.EQ.3)  FLAG = 0

      RETURN
      END             
C
C
C
      SUBROUTINE PROFILE_UPS(
     I                       YCRIT, YNORM, YVC, CC,
     O                       FLAG)

C     Compute a profile upstream using either full flow
C     or part full flow 

      IMPLICIT NONE

      INTEGER FLAG

      REAL*8 YCRIT, YNORM, YVC, CC

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local variables

      INTEGER NS, NE, NTAB, FULL

      REAL   YOVERD, DERIV

      REAL*8  YS, XS, XE, 
     A      AB, JB, JVC, ML, MR, AJC,
     B     FR, YC_RATIO, ABFB, JBFB, YBFS, ABFS, JBFS,
     C     P

C     Called functions

      REAL*8 E, RHSE

      EXTERNAL E, RHSE, FIND_YCYNYM, GET_YCYNYM, LKT_JDA


C***********************************************************************

      XS = STATION(NXS)
      XE = STATION(1) + HG_D

C     Check on the relative depth to decide if we need to compute
C     a free-surface profile or we can use the full pipe approximation.

      
      YC_RATIO = YCRIT/VERT_D
      IF(YC_RATIO.LT.FB_RATIO) THEN
        FULL = 0
      ELSEIF(YC_RATIO.LT.FS_RATIO) THEN
        FULL = 1
      ELSE
        FULL = 2
      ENDIF

      IF(FULL.GE.1) THEN
C       We compute the flow with the pipe flowing full and with 
C       pieozometric level at the exit a function of the flow.
        
        FR = QD/(AFULL_D*SQRT(GRAV_D*VERT_D))
        CALL LKTAB
     I            (ADRS_YOVERD, SNGL(FR), 0,
     O               YOVERD, NTAB, DERIV)
        YS = YOVERD*VERT_D
C        YBFB =  (ZBEX_D + COS_THETA*YS + (XS - XE)*
C     A                            (QD/KFULL_D)**2 - ZBVC_D)/COS_THETA
        ABFB = AFULL_D
        JBFB = JFULL_D
      ENDIF        
      IF(FULL.LE.1) THEN
C       Compute water-surface profile to vena contracta starting at critical 
C       depth. 
        YS = YCRIT*1.001D0
        CALL SFWSP(
     I             STDOUT_D, XS, XE, YS, EPS_D,
     I             E, RHSE, FIND_YCYNYM, NMAX,
     O             NS, NE, XVEC, YVEC, YCVEC, YNVEC, RFLAG_D)
        NS_RETURN = NS
        NE_RETURN = NE
        XE_RETURN = XE
      
        YBFS = YVEC(NS)
        CALL LKT_JDA(
     I               0, XE, YBFS,
     O               ABFS, JBFS)
      ENDIF

      IF(FULL.EQ.0) THEN
C        YB = YBFS
        AB = ABFS
        JB = JBFS
      ELSEIF(FULL.EQ.1) THEN
        P = (YC_RATIO - FB_RATIO)/(FS_RATIO - FB_RATIO)
C        YB = YBFS + P*(YBFB - YBFS)
        AB = ABFS + P*(ABFB - ABFS)
        JB = JBFS + P*(JBFB - JBFS)
      ELSE
C        YB = YBFB
        AB = ABFB
        JB = JBFB
      ENDIF
      YC_RETURN = YS
      FULL_RETURN = FULL
      CALL LKT_JDA(
     I             0, XE, YVC,
     O             AJC, JVC)

C     Compute residual in the momentum balance.

      ML = QD**2/(HG_D*BG_D*CC) + GRAV_D*JVC
      MR = QD**2/AB + GRAV_D*JB
      IF(ML.GE.MR) THEN
C       Vena contracta is not drowned
        FLAG = 1
      ELSE
C       Vena contracta is drowned.
        FLAG = 0
      ENDIF
      RETURN
      END
C
C
C
      REAL*8 FUNCTION FIND_FB_YCRIT_RESID(
     I                                    YC)

C     Residual function for FIND_FB_YCRIT

      IMPLICIT NONE
      REAL*8 YC

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local

      INTEGER NTAB
      REAL YOVERD, DERIV
      REAL*8 QCD, FR, YS, YB


C***********************************************************************

      CALL LKT_QCD(
     I           1, XS_D, YC,
     O           QCD)

      FR = QCD/(AFULL_D*SQRT(GRAV_D*VERT_D))
      
      IF(FR.GT.55.D0) FR = 55.D0
      CALL LKTAB
     I          (ADRS_YOVERD, SNGL(FR), 0,
     O             YOVERD, NTAB, DERIV)
C      WRITE(STDOUT_D,*) ' FR=',SNGL(FR),' YOVERD=',YOVERD
C      WRITE(STDOUT_D,*) ' VERT_D=',VERT_D
C      WRITE(STDOUT_D,*) ' ZBEX_D=',ZBEX_D,' ZBVC_D=',ZBVC_D
C      WRITE(STDOUT_D,*) ' XS_D=',XS_D,' XE_D=',XE_D
C      WRITE(STDOUT_D,*) ' QCD=',QCD
C      WRITE(STDOUT_D,*) ' KFULL_D=',KFULL_D
      YS = YOVERD*VERT_D
      YB =  (ZBEX_D + COS_THETA*YS + (XS_D - XE_D)*
     A                            (QCD/KFULL_D)**2 - ZBVC_D)/COS_THETA

C      WRITE(STDOUT_D,*) ' YS=',YS,' YB=',YB
      FIND_FB_YCRIT_RESID = (YB - VERT_D)/VERT_D
C      WRITE(STDOUT_D,*) ' RESID=',FIND_FB_YCRIT_RESID
C      WRITE(STDOUT_D,*) ' '
     
      RETURN
      END

C
C
C
      SUBROUTINE FIND_FB_YCRIT
     I                        (STDOUT,
     O                         ROOT) 

C     Find the critical depth  ratio that will cause the barrel to 
C     flow flow with certainty. 

      IMPLICIT NONE
      INTEGER STDOUT
      REAL*8 ROOT

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local

      INTEGER FLAG
      REAL*8 YL, FL, YR, FR, RATIO 


C     Functions
      REAL*8 FIND_FB_YCRIT_RESID

      EXTERNAL FIND_FB_YCRIT_RESID

      

C*********************************************************************** 
      XS_D = STATION(NXS)
      XE_D = STATION(1) + HG_D


C     Find an interval with a sign change. The critical flow at VERT_D
C     must cause a positive residual.  Otherwise there may be no
C     solution

      YR = VERT_D
      FR = FIND_FB_YCRIT_RESID(YR)
      IF(FR.LT.0.D0) THEN
        WRITE(STDOUT,*) 
     A      ' FIND_FB_YCRIT_RESID: No root. Max Qc too small.' 
        STOP 'Abnormal stop. Errors found.'
      ENDIF

      RATIO = 0.8D0
      YR = RATIO*VERT_D
      YL = -1.D0
 100  CONTINUE
        FR = FIND_FB_YCRIT_RESID(YR)
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
        IF(FR.GE.0.0D0) THEN
          GOTO 110
        ELSE
          YL = YR
          FL = FR
          RATIO = 0.5D0*(RATIO + 1.D0)
          YR = RATIO*VERT_D
          GOTO 100
        ENDIF
 110  CONTINUE
C     POSITIVE RESIDUAL FOUND- SEARCH FOR NEGATIVE RESIDUAL
 
      IF(YL.LT.0.0D0) THEN
 120    CONTINUE
          RATIO = RATIO - 0.05D0
          YL = RATIO*VERT_D
          FL = FIND_FB_YCRIT_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
          IF(FL.LE.0.0D0) THEN
            GOTO 130
          ELSE
            FR = FL
            YR = YL
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF

C     Found an interval that contains a root. 
C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG_D, EPSF_D, FIND_FB_YCRIT_RESID,
     M             YL, YR, FL, FR,
     O             ROOT, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' YC =',ROOT, ' KNT=',INT(FR)
      IF(FLAG.EQ.3)  FLAG = 0
      ROOT = ROOT/VERT_D
      RETURN
      END             
C
C
C
      SUBROUTINE BISECT_D
     I                 (STDOUT, EPSX, EPSF, F,
     M                  YL, YR, FL, FR,
     O                  ROOT, FLAG)

C     Try bisection for a root.  

      IMPLICIT NONE
      INTEGER FLAG, STDOUT

      REAL*8 EPSX, EPSF, YL, YR, FL, FR, ROOT

      INCLUDE 'ufgate_d.cmn'

      REAL*8 F
      EXTERNAL F

C     Local
      
      INTEGER I

      REAL*8 YM, FM
C***********************************************************************
      FLAG = 0
      I = 50
100   CONTINUE
        YM = 0.5D0*(YL + YR)
        FM = F(YM)
C        WRITE(STDOUT,*) ' '
C        WRITE(STDOUT,*) ' I=',I,' FULL=',FULL_RETURN
C        WRITE(STDOUT,*) ' YL=',YL,' YR=',YR
C        WRITE(STDOUT,*) ' FL=',FL,' FR=',FR
C        WRITE(STDOUT,*) ' YM=',YM,' FM=',FM
        IF(FM.LT.0.D0) THEN
          IF(FL.LT.0.D0) THEN
            YL = YM
            FL = FM
          ELSE
            YR = YM
            FR = FM
          ENDIF
        ELSE
          IF(FL.LT.0.D0) THEN
            YR = YM
            FR = FM
          ELSE
            YL = YM
            FL = FM
          ENDIF
        ENDIF
        IF(ABS(FM).LT.EPSF) THEN
          ROOT = YM
          FL = FM
          RETURN
        ENDIF
C        IF(ABS(YR - YL)/(ABS(YL) + ABS(YR)).LT.EPSX) THEN
C          ROOT = YM
C          FL = FM
C          FLAG = 3
C          RETURN
C        ENDIF

        I = I -1 
        IF(I.EQ.0) THEN
          FLAG = 1
          ROOT = YM
          FL = FM
          RETURN
        ENDIF
        GOTO 100
      END
C
C
C
      SUBROUTINE CONSTRUCT_YOVERD(
     I                            STDOUT, 
     M                            FTPT,
     O                            ADRS_YOVERD)

C     Construct the table for the piezometric level versus discharge
C     for a closed conduit. Based on work of J. L. French, 1956.
C     Argument has been changed to be dimensionless.  Argument 
C     extends to heads about 40 times the vertical diameter of the
C     conduit.  This should be sufficient for most culverts.  

      IMPLICIT NONE

      INTEGER FTPT, ADRS_YOVERD, STDOUT


C     Called program units

      EXTERNAL PUT1D

C     Local

      INTEGER N
      PARAMETER( N=43)

      REAL ARG(N), F1(N), F2(N)

      DATA  ARG/
     A 0.040797, 0.44876, 0.67314, 0.89752, 1.1219, 1.3463, 1.5707,
     A           1.795  , 2.0194 , 2.2438 , 2.4682, 2.6925, 2.9169,
     B           3.1413 , 3.3657,  3.5901 , 3.8144, 4.0388, 4.2632,
     C           4.4876 , 4.712 ,  4.9363 , 5.1607, 5.3851, 5.6095,
     D           5.8339 , 6.0582,  6.2826 , 6.507 , 6.7314, 6.9557,
     E           7.1801 , 7.4045,  7.6289 , 7.8533, 8.0776, 8.302 ,
     F           8.5264 , 8.7508,  8.9752,  9.1995, 22.438, 56.095/
      DATA  F1/
     A      1.0, 0.9    , 0.845  , 0.773  , 0.684  , 0.641 , 0.615,
     A           0.6    , 0.588  , 0.581  , 0.575  , 0.571 , 0.567,              
     B           0.56238, 0.55899, 0.55599, 0.55331, 0.5509, 0.54872,            
     C           0.54674, 0.54493, 0.54327, 0.54174, 0.54033,0.53902,            
     D           0.5378 , 0.53666, 0.5356 , 0.5346 , 0.53366,0.53278,            
     E           0.53195, 0.53117, 0.53042, 0.52972, 0.52905,0.52841,            
     F           0.5278,  0.52722, 0.52667, 0.52615, 0.51271,0.50605/            
      DATA  F2/  -0.24893,  -0.23751,  -0.26664,  -0.39397,  -0.31008, 
     A           -0.13059,  -0.09012, -0.057113, -0.042426, -0.027218, 
     B          -0.022517, -0.016419, -0.018771, -0.015044,  -0.01239, 
     C          -0.011371, -0.010182,-0.0092525,-0.0084272,-0.0077166, 
     D         -0.0070933,-0.0065453,-0.0060602,-0.0056288,-0.0052432, 
     E         -0.0048971,-0.0045853,-0.0043032,-0.0040472,-0.0038141, 
     F         -0.0036012,-0.0034062,-0.0032271,-0.0030622,  -0.00291, 
     G         -0.0027693,-0.0026389,-0.0025175,-0.0024057,-0.0022973, 
     H         -0.0022114,-0.0003429,-0.00012507/                       
C***********************************************************************
      CALL PUT1D
     I          (STDOUT, -1, 4, N, ARG, F1, F2,                              
     M           FTPT,                                                       
     O           ADRS_YOVERD)                                                
      RETURN                                                                 
      END                                                                    
C
C
C
      REAL*8 FUNCTION FIND_FCQ_RESID(
     I                                YVC)

C     Compute the residual function for subroutine FIND_FCQ.
c     The unknown is the depth at the vena contracta


      IMPLICIT NONE

      REAL*8 YVC

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local variables

      INTEGER NS, NE, NTAB, FULL

      REAL   YOVERD, DERIV

      REAL*8 YCRIT, YNORM, YMAX, YS, XS, XE, 
     A      AB, JB, JVC, ML, MR, AJC, AT,
     B     FR, YC_RATIO,  ABFB,  JBFB, YBFS, ABFS, JBFS,
     C     P

C     Called functions

      REAL*8 E, RHSE

      EXTERNAL E, RHSE, FIND_YCYNYM, GET_YCYNYM, LKT_JDA

C***********************************************************************

      RFLAG_D = 0

C     Compute the flow from the given YVC and the energy equation. 
C     Note that only valid values of YVC should appear here.

      AT = HG_D*BG_D*CC_D*CD_D      
      QD = AT*SQRT(TWOG_D*
     A     (H1_D - YVC*COS_THETA + HG_D*SIN_THETA)/
     B     (1.D0 - ALPHA1_D*(AT/A1_D)**2))


C     Find the normal and critical depth at culvert exit.
      YCRIT = YVC
      YNORM = 0.5D0*(YVC + VERT_D)
      CALL GET_YCYNYM
     I                (STDOUT_D, STATION(NXS), EPSARG_D, EPSF_D, 
     I                 EPSABS_D,
     M                 YCRIT, YNORM,
     O                 YMAX, RFLAG_D)

    
      IF(YNORM.GT.0.D0) THEN
        IF(YCRIT.GT.YNORM) THEN
C         The profile is supercritical. Some problem here.
          RFLAG_D = -1
          GOTO 1000
        ENDIF
      ENDIF

      XS = STATION(NXS)
      XE = STATION(1) + HG_D

C     Check on the relative depth to decide if we need to compute
C     a free-surface profile or we can use the full pipe approximation.

      
      YC_RATIO = YCRIT/VERT_D
      IF(YC_RATIO.LT.FB_RATIO) THEN
        FULL = 0
      ELSEIF(YC_RATIO.LT.FS_RATIO) THEN
        FULL = 1
      ELSE
        FULL = 2
      ENDIF

      IF(FULL.GE.1) THEN
C       We compute the flow with the pipe flowing full and with 
C       pieozometric level at the exit a function of the flow.
        
        FR = QD/(AFULL_D*SQRT(GRAV_D*VERT_D))
        CALL LKTAB
     I            (ADRS_YOVERD, SNGL(FR), 0,
     O               YOVERD, NTAB, DERIV)
        YS = YOVERD*VERT_D
C        YBFB =  (ZBEX_D + COS_THETA*YS + (XS - XE)*
C     A                            (QD/KFULL_D)**2 - ZBVC_D)/COS_THETA
        ABFB = AFULL_D
        JBFB = JFULL_D
        Y_RETURN = YS
      ENDIF        
      IF(FULL.LE.1) THEN
C       Compute water-surface profile to vena contracta starting at critical 
C       depth. 
        YS = YCRIT*1.001D0
        CALL SFWSP(
     I             STDOUT_D, XS, XE, YS, EPS_D,
     I             E, RHSE, FIND_YCYNYM, NMAX,
     O             NS, NE, XVEC, YVEC, YCVEC, YNVEC, RFLAG_D)
        NS_RETURN = NS
        NE_RETURN = NE
        XE_RETURN = XE
      
        YBFS = YVEC(NS)
        CALL LKT_JDA(
     I               0, XE, YBFS,
     O               ABFS, JBFS)
      ENDIF

      IF(FULL.EQ.0) THEN
C        YB = YBFS
        AB = ABFS
        JB = JBFS
      ELSEIF(FULL.EQ.1) THEN
        P = (YC_RATIO - FB_RATIO)/(FS_RATIO - FB_RATIO)
C        YB = YBFS + P*(YBFB - YBFS)
        AB = ABFS + P*(ABFB - ABFS)
        JB = JBFS + P*(JBFB - JBFS)
      ELSE
C        YB = YBFB
        AB = ABFB
        JB = JBFB
      ENDIF
      YC_RETURN = YCRIT
      FULL_RETURN = FULL
      CALL LKT_JDA(
     I             0, XE, YVC,
     O             AJC, JVC)

C     Compute residual in the momentum balance.

      ML = QD**2/(HG_D*BG_D*CC_D) + GRAV_D*JVC*COS_THETA
      MR = QD**2/AB + GRAV_D*JB*COS_THETA


      FIND_FCQ_RESID = (ML - MR)/MR

      RETURN
1000  CONTINUE
C     Force the residual to be such that YC will be reduced.
      FIND_FCQ_RESID =  1.D0
      RETURN

      END

C
C
C
      SUBROUTINE FIND_FCQ(
     I                     STDOUT, CD,  H1,
     O                     FCQ, YVC, RFLAG)

C     Find the flow when FC flow controls.  That is, critical depth
C     at the culvert exit drowns the orifice at the sluice gate.

      IMPLICIT NONE

      INTEGER STDOUT, RFLAG

      REAL*8  CD, H1, FCQ, YVC

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local variables
      
      INTEGER I, FLAG, NEG, POS

      REAL A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1, DALPHA,
     A            QCT, TP
      REAL*8 YCL, FL, YCR, FR, Y1, Y, F, YVC_MIN,
     A       YVC_MAX

C     External functions

      REAL*8 FIND_FCQ_RESID

      EXTERNAL FIND_FCQ_RESID
C     ******************************FORMATS*****************************
50    FORMAT(/,' *BUG:XXX* Positive residual not found for FIND_FCQ')
52    FORMAT(/,' *BUG:XXX* Negative residual not found for FIND_FCQ')
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT FDBLRGF CLAIMS NONE IN',
     A       ' FIND_FIND_FCQ.')
 60   FORMAT(' *BUG:XXX* FDBLRGF ITERATION>100 IN FIND_FIND_FCQ')
C***********************************************************************
C     Set values in common block ufgate_d.cmn
      CD_D = CD
      H1_D = DBLE(H1)
C     Find the conditions at section 1

C      WRITE(STDOUT,*) ' Entering FIND_FCQ: CC=',CC,' CD=',CD
      Y1 = H1 + HDATUM_D - Z1B_D
      TP = SNGL(Y1)
      CALL XLKT22
     I           (APPTAB_D,
     M            TP,
     O            A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1, DALPHA,
     O            QCT)

      A1_D = DBLE(A1)
      ALPHA1_D = DBLE(ALPHA1)

      VERT_D = VERT_DVEC(NXS)

      YVC_MIN = 0.9999D0*HG_D*CC_D

      YVC_MAX = MIN(0.999D0*(H1_D + HG_D*SIN_THETA), VERT_D)

C     Explore the function


C      N = 50
C      DO 90 I= 1,N
C        YCL = YVC_MIN + DBLE(I-1)*(0.8*YVC_MAX - YVC_MIN)/DBLE(N)
C        RESID = FIND_FCQ_RESID(YCL)
C        WRITE(STDOUT,9782) YCL, RESID, FULL_RETURN
C9782  FORMAT(' YCL=',F10.4,' RESID=',1PE10.3,' FULL=',I5)
C90    CONTINUE



      Y = YVC_MIN
      NEG = 0
      POS = 0
      I = 100
100   CONTINUE
      IF(NEG+POS.LT.2) THEN
C       Evaluate the residual
        F = FIND_FCQ_RESID(Y)
        I = I - 1
        IF(I.EQ.0) THEN
          WRITE(STDOUT,*) ' NO SIGN CHANGE FOUND IN FIND_FCQ'
          STOP 'Abnormal stop. Errors found.'
        ENDIF

C        WRITE(STDOUT,9787) I, Y, F, FULL_RETURN
C9787  FORMAT(' I=',I5,' Y=',F10.4,' F=',1PE10.3,' FULL=',I5)
        IF(F.LT.0.D0) THEN
          FL = F
          YCL = Y
          NEG = 1
          Y = 0.7D0*YCL + 0.3D0*YVC_MAX
          GOTO 100
        ELSE
          FR = F
          YCR = Y
          POS = 1
          Y = 1.002D0*YCR
          GOTO 100
        ENDIF
      ENDIF


C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YCL=',YCL,' FL=',FL
C      WRITE(STDOUT,*) ' YCR=',YCR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG_D, EPSF_D, FIND_FCQ_RESID,
     M             YCL, YCR, FL, FR,
     O             YVC, FLAG)
   

C      WRITE(STDOUT,*) ' FL=',FL      
C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' YVC =',YVC, ' KNT=',INT(FR)
      IF(FLAG.EQ.3)  FLAG = 0
      
C     Find the station when the barrel is full.
      FULL_BARREL_STATION = -1.D0
C      WRITE(STDOUT,*) ' XE_RETURN=',XE_RETURN
C
C      DO 800 I=NS_RETURN,NE_RETURN
C        WRITE(STDOUT,5678) XVEC(I), YVEC(I)
C5678  FORMAT(F12.6,F12.6)
C800   CONTINUE

      DO 900 I=NE_RETURN,NS_RETURN,-1
        IF(YVEC(I).GE.VERT_D)  THEN
          FULL_BARREL_STATION = XVEC(I)
          GOTO 901
        ENDIF
900   CONTINUE
901   CONTINUE      
            
      FCQ = QD
    
      RETURN
      END
C
C
C
      REAL*8 FUNCTION FIND_FOLL_RESID(
     I                                H)

C     Compute the residual function for subroutine FIND_FOLL

      IMPLICIT NONE

      REAL*8 H

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local variables

      INTEGER NS, NE, NTAB, FULL

      REAL DERIV, YOVERD, RG

      REAL*8 YCRIT, YNORM, YMAX, YS, XS, XE, AVC, AB,
     A   JVC, JB, YVC, YB, MVC, MB, FR, YC_RATIO, YBFB,
     B   ABFB,JBFB, YBFS, ABFS, JBFS, P

C     Called functions

      REAL*4 FINDCC

      REAL*8 E, RHSE

      EXTERNAL E, RHSE, FIND_YCYNYM, GET_YCYNYM, LKT_JDA, FINDCC

C***********************************************************************
      RG = HG_D/H
      IF(RG.GT.1.0) RG = 1.0
      CC_D = FINDCC(RG, SNGL(CONCC_D), CCTAB_D)      
      RFLAG_D = 0
C     Find FO flow 
      CALL UFG_FNDFOQ
     I            (STDOUT_D, APPTAB_D, H, HDATUM_D, HG_D, BG_D, Z1B_D, 
     I             CC_D, CD_D, TWOG_D,
     O                    QD)

C     Find the critical and normal depths at the culvert exit.
      VERT_D = VERT_DVEC(NXS)
      YCRIT = 0.5D0*VERT_D
      YNORM = 0.5D0*VERT_D
      CALL GET_YCYNYM
     I                (STDOUT_D, STATION(NXS), EPSARG_D, EPSF_D, 
     I                 EPSABS_D,
     M                 YCRIT, YNORM,
     O                 YMAX, RFLAG_D)

      IF(YNORM.GT.0.D0) THEN
        IF(YCRIT.GT.YNORM) THEN
C         The profile is supercritical.  H is probably too large
          RFLAG_D = -1
          GOTO 1000
        ENDIF
      ENDIF
          
      YC_RETURN = YCRIT

      XS = STATION(NXS)
      XE = STATION(1) + HG_D

C     Check on the relative depth to decide if we need to compute
C     a free-surface profile or we can use the full pipe approximation.

      YC_RATIO = YCRIT/VERT_D
      IF(YC_RATIO.LT.FB_RATIO) THEN
        FULL = 0
      ELSEIF(YC_RATIO.LT.FS_RATIO) THEN
        FULL = 1
      ELSE
        FULL = 2
      ENDIF
      IF(FULL.GE.1) THEN
C       We compute the flow with the pipe flowing full and with 
C       pieozometric level at the exit a function of the flow.
        
        FR = QD/(AFULL_D*SQRT(GRAV_D*VERT_D))
        CALL LKTAB
     I            (ADRS_YOVERD, SNGL(FR), 0,
     O               YOVERD, NTAB, DERIV)
        YS = YOVERD*VERT_D
        YBFB =  (ZBEX_D + COS_THETA*YS + (XS - XE)*
     A                            (QD/KFULL_D)**2 - ZBVC_D)/COS_THETA
        ABFB = AFULL_D
        JBFB = JFULL_D
      ENDIF        
      IF(FULL.LE.1) THEN
C       Compute water-surface profile to vena contracta starting at critical 
C       depth. 
        YS = YCRIT*1.001D0
        CALL SFWSP(
     I             STDOUT_D, XS, XE, YS, EPS_D,
     I             E, RHSE, FIND_YCYNYM, NMAX,
     O             NS, NE, XVEC, YVEC, YCVEC, YNVEC, RFLAG_D)
        NS_RETURN = NS
        NE_RETURN = NE
        XE_RETURN = XE
      
        YBFS = YVEC(NS)
        CALL LKT_JDA(
     I               0, XE, YBFS,
     O               ABFS, JBFS)
      ENDIF
      YC_RETURN = YS
      FULL_RETURN = FULL
      IF(FULL.EQ.0) THEN
        YB = YBFS
        AB = ABFS
        JB = JBFS
      ELSEIF(FULL.EQ.1) THEN
        P = (YC_RATIO - FB_RATIO)/(FS_RATIO - FB_RATIO)
        YB = YBFS + P*(YBFB - YBFS)
        AB = ABFS + P*(ABFB - ABFS)
        JB = JBFS + P*(JBFB - JBFS)
      ELSE
        YB = YBFB
        AB = ABFB
        JB = JBFB
      ENDIF

            
C     Compute simple momentum balance at the vena contracta.  Add
C     gravity force later for rectangular channels with free surface 
C     flow.

      Y_RETURN = YB
    
      YVC = CC_D*HG_D
      CALL LKT_JDA(
     I            0, XE, YVC,
     O            AVC, JVC)

      MVC = QD**2/AVC + GRAV_D*JVC*COS_THETA
      MB = QD**2/AB  + GRAV_D*JB*COS_THETA

      FIND_FOLL_RESID = (MVC - MB)/MB
      RETURN
1000  CONTINUE
C     Force the residual to be such that H will be reduced.
      FIND_FOLL_RESID =  1.D50
      RETURN

      END
C
C
C
      SUBROUTINE FIND_FOLL(
     I                     STDOUT, CD, MAX_H1,
     M                     H1,
     O                     RFLAG)

C     Find the head at section 1 for the lower limit of FO flow when
C     FC flow is found to exist. 

      IMPLICIT NONE

      INTEGER STDOUT, RFLAG

      REAL*8 CD, H1, MAX_H1

      INCLUDE 'barrel.cmn'
      INCLUDE 'ufgate_d.cmn'

C     Local variables
      
      INTEGER  FLAG

      REAL*8 HL, FL, HR, FR, RESID_AT_H1FWUL, RESID_AT_MAX_H1

C     External functions

      REAL*8 FIND_FOLL_RESID

      EXTERNAL FIND_FOLL_RESID
C     ******************************FORMATS*****************************
50    FORMAT(/,' *BUG:XXX* Positive residual not found for FIND_FOLL')
52    FORMAT(/,' *BUG:XXX* Negative residual not found for FIND_FOLL')
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT FDBLRGF CLAIMS NONE IN',
     A       ' FIND_FIND_FOLL.')
 60   FORMAT(' *BUG:XXX* FDBLRGF ITERATION>100 IN FIND_FIND_FOLL')
C***********************************************************************
C     Set values in common block ufgate_d.cmn
      CD_D = CD

C     Evaluate the residual at the initial value of H1.  This should 
C     be H1FWUL.  
      RESID_AT_H1FWUL = FIND_FOLL_RESID(H1)
      RESID_AT_MAX_H1 = FIND_FOLL_RESID(MAX_H1)
      IF(RESID_AT_H1FWUL.LT.0.D0) THEN
C       FO with normal CC is submerged at H1FWUL.  H1FOLL exists
        IF(RESID_AT_MAX_H1.GT.0.D0) THEN
C         Maximum head is greater than H1FOLL.  We appear to have
C         a root.  Experiments indicate there is only one in this
C         case.  

          HL = H1
          FL = RESID_AT_H1FWUL
          HR = MAX_H1
          FR = RESID_AT_MAX_H1


C         Try bisection here
C          CALL BISECT_D
C     I                (STDOUT, EPSARG_D, EPSF_D, FIND_FOLL_RESID,
C     M                 HL, HR, FL, FR,
C     O                 H1, FLAG)

C          WRITE(STDOUT,*) ' BISECTION: H1FOLL=',H1,' FLAG=',FLAG
C          WRITE(STDOUT,*) ' FL=',FL

C          WRITE(STDOUT,*) ' Calling FDBLRGF with: HL=',HL,' FL=',FL
C          WRITE(STDOUT,*) ' HR=',HR,' FR=',FR

           CALL FDBLRGF
     I                (EPSARG_D, EPSF_D, FIND_FOLL_RESID,
     M                 HL, HR, FL, FR,
     O                 H1, FLAG)
   

C           WRITE(STDOUT,*) ' FL=',FL,' FULL=',FULL_RETURN
C           WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C           WRITE(STDOUT,*) ' H1FOLL =',H1, ' KNT=',INT(FR)

          
           IF(FLAG.EQ.3)  FLAG = 0

        
        ELSE
C         Maximum head is less than H1FOLL.  Thus the free flow 
C         at this gate opening is always FC. 

          WRITE(STDOUT,*) ' H1FOLL > maximum head at section 1.'
          H1 = MAX_H1
        ENDIF
      ELSE
C       Probable that H1FOLL does not exist.  If the 
C       residual has an initial decline as the head is increased
C       slightly, we may have a two root solution.  Skip 
C       for now unless we can define it more closely

        WRITE(STDOUT,*) ' H1FOLL may not exist. Flow is always FO'
        H1 = 0.D0
      ENDIF

      RETURN

      END

      
C
C
C
      SUBROUTINE   UFG_FNDFOQ
     I            (STDOUT, APPTAB, H1, HDATUM, HG, BG, Z1B, CC, CD,
     I                    TWOG,
     O                    Q)
 
C     + + + PURPOSE + + +
C     Find the free orifice flow. 
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER APPTAB, STDOUT
      REAL*8 H1, HDATUM, Q, HG, BG, Z1B, CC, CD, TWOG
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     H1     - head at section 1
C     HDATUM - Datum for measuring head
C     Q      - Flowrate
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'barrel.cmn'
 
C     + + + LOCAL VARIABLES + + +
      REAL A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1, DALPHA,
     A     QCT, TP
      REAL*8 HVC, Y1, AT
 
C     + + + INTRINSICS + + +
      INTRINSIC SQRT, SNGL
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL XLKT22
C***********************************************************************
      
      Y1 = H1 + HDATUM - Z1B
      TP = SNGL(Y1)
      CALL XLKT22
     I           (APPTAB,
     M            TP,
     O            A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1, DALPHA,
     O            QCT)
      
C      WRITE(STDOUT,*) ' Y1=',Y1,' A1=',A1
C      WRITE(STDOUT,*) ' ALPHA1=',ALPHA1
C      WRITE(STDOUT,*) ' SIN_THETA=',SIN_THETA
      AT = CD*CC*BG*HG
      HVC = CC*HG - SIN_THETA*HG 
C      WRITE(STDOUT,*) ' HVC=',HVC
      Q = AT*SQRT(TWOG*(H1 - HVC)/
     A              (1.D0 - DBLE(ALPHA1)*(AT/DBLE(A1))**2))


      RETURN
      END
C
C
C

C      REAL*8 FUNCTION FW_CC_FO_RESID(
C     I                            CC)
C
CC     Compute values for the free-weir contraction coefficient 
CC     residual function.
C
C      IMPLICIT NONE
C
C      REAL*8 CC
C
C      INCLUDE 'fwcc.cmn'
CC***********************************************************************
C
C      FW_CC_FO_RESID =  M1/CC**2 + CC - E_AT_1
C      RETURN
C      END
CC
CC
CC
C      SUBROUTINE FIND_FW_CC_FO(
C     I                  STDOUT, EPSARG, EPSF, EPSABS, APPTAB, GRAV, 
C     I              H1FWULA, HDATUM, Z1B, STATION_1, HG, BG, Y2HAT,
C     O                  CD, FW_CC, FOFLAG, FLAG)  
C
CC     Find the contraction coefficient that forces the FO equation 
CC     to produce the same flow as the flow at the FW limit.
C
C      IMPLICIT NONE
C
C      INTEGER APPTAB, STDOUT, FLAG, FOFLAG
C
C      REAL*8 EPSARG, EPSF, EPSABS, GRAV, H1FWULA, Z1B, STATION_1, 
C     A       HG, FW_CC, Y2HAT, BG, CD, HDATUM, CC
C
C
C      INCLUDE 'barrel.cmn'
C      INCLUDE 'fwcc.cmn'
C
C
CC     Local variables
C
C      INTEGER I
C
C      REAL*8 TWOG, YL, YR, FL, FR, L12, H1FWUL,
C     A       FR2VC, YSEQ, Q 
C      REAL A1, T1, DT1, J1, K1, DK1, BETA1, DBETA1, ALPHA1,
C     A     DALPHA1, QC1, AB, TB, DTB, JB, KB, DKB, BETAB, DBETAB, 
C     B     ALPHAB, DALPHAB, QCB
C
CC     External program units
C
C      REAL*8 FW_CC_FO_RESID
C      EXTERNAL XLKT22, FW_CC_FO_RESID, FDBLRGF
CC     *****************************FORMATS******************************
C50    FORMAT(/,' *BUG:XXX* Positive residual not found for FW_CC_FO')
C52    FORMAT(/,' *BUG:XXX* Negative residual not found for FW_CC_FO')
C 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT FDBLRGF CLAIMS NONE IN',
C     A       ' FIND_FW_CC_FO.')
C 60   FORMAT(' *BUG:XXX* FDBLRGF ITERATION>100 IN FIND_FW_CC_FO')
CC***********************************************************************
C      FLAG = 0
C      TWOG = 2.D0*GRAV
C      
CC     Increment the section 1 limit 
C
C      H1FWUL = H1FWULA 
CC     Find elements at section 1
C      WRITE(STDOUT,*) ' H1=',H1FWUL
C
C      CALL XLKT22
C     I           (APPTAB,
C     M            SNGL(H1FWUL + HDATUM - Z1B),
C     O            A1, T1, DT1, J1, K1, DK1, BETA1, DBETA1, ALPHA1,
C     O            DALPHA1, QC1)
C
CC     Find elements at section 2      
C      CALL XLKT22
C     I           (XSEC_ADRS(1),
C     M            SNGL(HG),
C     O            AB, TB, DTB, JB, KB, DKB, BETAB, DBETAB, ALPHAB,
C     O            DALPHAB, QCB)
C
CC     Compute loss term for approach reach.
C
C      L12 = STATION(1) - STATION_1
C      LOSS = L12*QD**2/(DBLE(K1)*DBLE(KB))
C
CC     Compute the RHS of the equation.  Put in E_AT_1
C
C      E_AT_1 = DBLE(ALPHA1)*(QD/DBLE(A1))**2/TWOG + H1FWUL 
C     A          + SIN_THETA*HG  - LOSS
C
CC     Compute the velocity head factor at the vena contract
C
C      M1 = (QD/(BG*HG))**2/TWOG
C
C
CC     Rescale by HG
C
C      E_AT_1 = E_AT_1/HG
C
C      M1 = M1/HG
C
C
CC     Now setup for solving for FW_CC_FO
C
CC     Explore the function
C
CC      YR = 0.01D0
CC      DO 90 I= 1,101
CC        WRITE(STDOUT,*) ' YR=',YR,' RESID=',FW_CC_FO_RESID(YR)
CC        YR = YR + 0.01D0
CC90    CONTINUE
CC     SEARCH FOR A positive RESIDUAL
C      YR = 1.1D0
C      YL = -1.D0
C 100  CONTINUE
C        FR = FW_CC_FO_RESID(YR)
CC        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
C        IF(FR.GE.0.0D0) THEN
C          GOTO 110
C        ELSE
C          IF(FR.GT.-1.D5) THEN
C            YL = YR
C            FL = FR
C          ENDIF
C          YR = YR  - 0.01D0
C          IF(YR.LT.0.01D0) THEN
C            WRITE(STDOUT,50) 
C            FLAG = 1
C            GOTO 1000
C          ENDIF
C          GOTO 100
C        ENDIF
C 110  CONTINUE
CC     POSITIVE RESIDUAL FOUND- SEARCH FOR NEGATIVE RESIDUAL
C 
C      IF(YL.LT.0.0D0) THEN
C        YL = YR
C 120    CONTINUE
C          YL = YL - 0.01D0
C          IF(YL.LE.0.01D0) THEN
C            WRITE(STDOUT, 52) 
C            FLAG = 1
C            GOTO 1000
C          ENDIF
C          FL = FW_CC_FO_RESID(YL)
CC        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
C          IF(FL.LE.0.0D0) THEN
C            GOTO 130
C          ELSE
C            FR = FL
C            YR = YL
C            GOTO 120
C          ENDIF
C 130    CONTINUE
C      ENDIF
C 
CC     WE HAVE A SIGN CHANGE OR ONE OR BOTH POINTS HAVE ZERO RESIDUAL
C
CC      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
CC      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
C
C      CALL FDBLRGF
C     I            (EPSARG, EPSF, FW_CC_FO_RESID,
C     M             YL, YR, FL, FR,
C     O             FW_CC, FLAG)
C   
C
CC      WRITE(STDOUT,*) ' FL=',FL      
CC      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
CC      WRITE(STDOUT,*) ' FW_CC_FO =',FW_CC
C      IF(FLAG.EQ.3)  FLAG = 0
C
CC     Compute sequent depth for a jump from the jet implied by 
CC     FW_CC.
C
C      FR2VC = QD**2/(BG**2*(HG*FW_CC)**3*GRAV)
C      YSEQ = FW_CC*HG*0.5D0*(SQRT(1.D0 + 8.D0*FR2VC) - 1.D0)
C
C      IF(YSEQ.GT.HG) THEN
C        WRITE(STDOUT,*) ' FO passed test'
C        FOFLAG = 1
C      ELSE
C        FOFLAG = 0
C      ENDIF
C
CC      WRITE(STDOUT,*) ' YSEQ=',YSEQ,' HG=',HG
CC     Compute a coefficient of discharge applied to the velocity head
CC     at the vena contracta to approximate the losses for all heads for 
CC     this gate opening. 
C
C
C      
C      CD = SQRT(1.D0/(1.D0 + TWOG*L12*(HG*BG*FW_CC)**2/
C     A                      (DBLE(K1)*DBLE(KB))))
C
C
C      WRITE(STDOUT,*) ' CD=',CD
C
CC     Try to compute Q from the basic FO solution
C
C      Q = HG*BG*FW_CC*SQRT((TWOG*(H1FWUL - FW_CC*HG + SIN_THETA*HG))/
C     A   (1. - ALPHA1*(HG*BG*FW_CC/A1)**2 + 
C     B   TWOG*(HG*BG*FW_CC)**2*L12/(K1*KB)))
C
C      WRITE(STDOUT,*) ' Q FO inverse=',Q
C
C
C      IF(FLAG.EQ.1) THEN
C        WRITE(STDOUT, 54)
C      ELSEIF(FLAG.EQ.2) THEN
C        WRITE(STDOUT,60)
C      ENDIF
C 
C1000  CONTINUE 
C      RETURN
C      END

C
C
C
      REAL*8 FUNCTION FW_CC_SO_RESID(
     I                            CC)

C     Compute values for the free-weir contraction coefficient 
C     residual function.

      IMPLICIT NONE

      REAL*8 CC

C     Local variables

      REAL*8 JVCR, JVCL, YVC, AJC

      INCLUDE 'fwcc.cmn'
C***********************************************************************
      YVC = ECON - EFAC/CC**2
      CALL LKT_JDA(
     I            0, STAUP+HGT, YVC,
     O            AJC, JVCL)

      JVCR = M1 - M2/CC
      
      FW_CC_SO_RESID = (JVCR - JVCL)/JVCR
      
      RETURN
      END
C
C
C
      SUBROUTINE FIND_FW_CC_SO(
     I                  STDOUT, EPSARG, EPSF, EPSABS, APPTAB, GRAV, 
     I                  Z1FWULA, Z1B, STATION_1, HG, BG, Y2HAT, 
     O                  CD, FW_CC, SOFLAG, FLAG)  

C     Find the contraction coefficient that forces the SO equation 
C     to produce the same flow as the flow at the FW limit.

      IMPLICIT NONE

      INTEGER APPTAB, STDOUT, FLAG, SOFLAG

      REAL*8 EPSARG, EPSF, EPSABS, GRAV, Z1FWULA, Z1B, STATION_1, 
     A       HG, FW_CC, Y2HAT, BG, CD


      INCLUDE 'barrel.cmn'
      INCLUDE 'fwcc.cmn'


C     Local variables

      REAL*8 TWOG, YL, YR, FL, FR, YVC, L12, CCMIN, Z1FWUL, LOSS,
     A       AB, JB

      REAL A1, T1, DT1, J1, K1, DK1, BETA1, DBETA1, ALPHA1,
     A     DALPHA1, QC1, TB, DTB, DKB, BETAB, DBETAB, 
     B     ALPHAB, DALPHAB, QCB, A, J, KLOSS, TP

C     External program units

      REAL*8 FW_CC_SO_RESID
      EXTERNAL XLKT22, FW_CC_SO_RESID, FDBLRGF
C     *****************************FORMATS******************************
50    FORMAT(/,' *BUG:XXX* Positive residual not found for FW_CC_SO')
52    FORMAT(/,' *BUG:XXX* Negative residual not found for FW_CC_SO')
 54   FORMAT(' *BUG:XXX* SIGN CHNG ON ENTRY BUT FDBLRGF CLAIMS NONE IN',
     A       ' FIND_FW_CC_SO.')
 60   FORMAT(' *BUG:XXX* FDBLRGF ITERATION>100 IN FIND_FW_CC_SO')
C***********************************************************************
      FLAG = 0
      TWOG = 2.D0*GRAV

      STAUP = STATION(1)
      HGT = HG
      Z1FWUL = Z1FWULA 
C     Find elements at section 1

      TP = SNGL(Z1FWUL-Z1B)
      CALL XLKT22
     I           (APPTAB,
     M            TP,
     O            A1, T1, DT1, J1, K1, DK1, BETA1, DBETA1, ALPHA1,
     O            DALPHA1, QC1)

C     Find elements at section 2 at the depth computed from SFWSP 
C     excluding the distributed entrance losses.

      CALL LKT_JDA(
     I            1, STAUP+HGT, Y2HAT,
     O            AB, JB)

C     Find elements at the level used for computing approach friction
C     losses in CULVERT

      TP = SNGL(HG)
      CALL XLKT22
     I           (XSEC_ADRS(1),
     M            TP,
     O            A, TB, DTB, J, KLOSS, DKB, BETAB, DBETAB, ALPHAB,
     O            DALPHAB, QCB)


C     Compute loss term for approach reach.
      L12 = STATION(1) - STATION_1

      LOSS = L12*QD**2/(DBLE(K1)*DBLE(KLOSS))

      ECON = (DBLE(ALPHA1)*((QD/DBLE(A1))**2)/TWOG + Z1FWUL 
     A          + SIN_THETA*HG  - INVERT_Z(1) - LOSS)/COS_THETA

      EFAC = (QD/(HG*BG))**2/(TWOG*COS_THETA)


C     Compute the first momentum based term

      M1 = QD**2/(GRAV*AB*COS_THETA) + JB

C     Compute the second momentum based term

      M2 =  QD**2/(GRAV*HG*BG*COS_THETA)


      CCMIN = MAX(SQRT(EFAC/ECON), M2/M1)

      CCMIN = 1.0001D0*CCMIN



C     Now setup for solving for FW_CC

C     Explore the function

C      YR = 0.01D0
C      DO 90 I= 1,101
C        WRITE(STDOUT,*) ' YR=',YR,' RESID=',FW_CC_SO_RESID(YR)
C        YR = YR + 0.01D0
C90    CONTINUE
C     SEARCH FOR A positive RESIDUAL
      YR = 1.1D0
      YL = -1.D0
 100  CONTINUE
        FR = FW_CC_SO_RESID(YR)
C        WRITE(STDOUT,*) ' YR=',YR,' FR=',FR
        IF(FR.GE.0.0D0) THEN
          GOTO 110
        ELSE
          IF(FR.GT.-1.D5) THEN
            YL = YR
            FL = FR
          ENDIF
          YR = YR  - 0.1D0
          IF(YR.LT.CCMIN) THEN
C            WRITE(STDOUT,50) 
            FLAG = 1
            GOTO 1000
          ENDIF
          GOTO 100
        ENDIF
 110  CONTINUE
C     POSITIVE RESIDUAL FOUND- SEARCH FOR NEGATIVE RESIDUAL
 
      IF(YL.LT.0.0D0) THEN
        YL = YR
 120    CONTINUE
          YL = YL - 0.1D0
          IF(YL.LE.CCMIN) THEN
C            WRITE(STDOUT, 52) 
            FLAG = 1
            GOTO 1000
          ENDIF
          FL = FW_CC_SO_RESID(YL)
C        WRITE(STDOUT,*) ' YL=',YL,' FL=',FL
          IF(FL.LE.0.0D0) THEN
            GOTO 130
          ELSE
            FR = FL
            YR = YL
            GOTO 120
          ENDIF
 130    CONTINUE
      ENDIF
 
C     WE HAVE A SIGN CHANGE OR ONE OR BOTH POINTS HAVE ZERO RESIDUAL

C      WRITE(STDOUT,*) ' Calling FDBLRGF with: YL=',YL,' FL=',FL
C      WRITE(STDOUT,*) ' YR=',YR,' FR=',FR

      CALL FDBLRGF
     I            (EPSARG, EPSF, FW_CC_SO_RESID,
     M             YL, YR, FL, FR,
     O             FW_CC, FLAG)
   
      IF(FLAG.EQ.3) FLAG = 0

      YVC = SQRT(M1 - M2/FW_CC)
C      WRITE(STDOUT,*) ' Return from FDBLRGF with: FLAG=',FLAG
C      WRITE(STDOUT,*) ' FW_CC_SO =',FW_CC,' YVC=',YVC
      SOFLAG = 0
      IF(YVC.GT.HG*FW_CC) THEN
        WRITE(STDOUT,*) ' SO flow passed test.'
        SOFLAG = 1
      ENDIF
 
      IF(FLAG.EQ.1) THEN
        WRITE(STDOUT, 54)
      ELSEIF(FLAG.EQ.2) THEN
        WRITE(STDOUT,60)
      ENDIF
 
C     Compute a coefficient of discharge applied to the velocity head
C     at the vena contracta to approximate the losses for all heads for 
C     this gate opening. 


      
      CD = SQRT(1.D0/(1.D0 + TWOG*L12*(HG*BG*FW_CC)**2/
     A                      (DBLE(K1)*DBLE(KLOSS))))


C      WRITE(STDOUT,*) ' CD=',CD, ' A1=',A1
       RETURN

1000  CONTINUE
C     Compute a loss CD using a contraction coef of 1.0
      CD = SQRT(1.D0/(1.D0 + TWOG*L12*(HG*BG)**2/
     A                      (DBLE(K1)*DBLE(KLOSS))))
      FW_CC = 0.D0
      RETURN
      END


C
C
C
      INTEGER FUNCTION FIND_BARREL_INTERVAL(
     I                                     X)

C     Find the barrel interval containing X.
C     LOC points to the index for the upstream end of the barrel.

      IMPLICIT NONE

      REAL*8  X

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER OLDLOC, LOC
      DATA OLDLOC/1/

C***********************************************************************
100   CONTINUE
        IF(X.GE.STATION(OLDLOC)) THEN
          IF(X.LE.STATION(OLDLOC+1)) THEN
            LOC = OLDLOC
            GOTO 110
          ELSE
            OLDLOC = OLDLOC + 1
          ENDIF
        ELSE
          OLDLOC = OLDLOC - 1
        ENDIF
        GOTO 100
110   CONTINUE
      FIND_BARREL_INTERVAL = LOC
      RETURN
      END
C
C
C
      SUBROUTINE  SET_BARREL_INTERVAL(
     A                                 I)

C     Set the values needed for integration and computations in the
C     barrel cross section interval with upstream-end index I. 

      IMPLICIT NONE
      INTEGER I

      INCLUDE 'barrel.cmn'

C***********************************************************************
      X_L = STATION(I)
      Z_L = INVERT_Z(I)
      ADRS_L = XSEC_ADRS(I)
      YMAX_L = XSEC_YMAX(I)
      X_R = STATION(I+1)
      Z_R = INVERT_Z(I+1)
      ADRS_R = XSEC_ADRS(I+1)
      YMAX_R = XSEC_YMAX(I+1)

C     SLOPE_FLAG =0: slope is constant between cross sections.
C     SLOPE_FLAG >0: slope varies between cross sections.
      SLOPE_FLAG = INVERT_SLOPE_CAT(I) + INVERT_SLOPE_CAT(I+1)
      PRIS_FLAG = CHANNEL_VARIATION(I)
      IF(SLOPE_FLAG.EQ.0) THEN
        SIN_THETA = SINE_THETA(I)
        COS_THETA = COSINE_THETA(I)
      ENDIF
      RETURN
      END      
C
C
C
      SUBROUTINE SFWSP(
     I                 STDOUT, XS, XE, YS, EPS,
     I                 P, F, FIND_YCYNYM, NMAX,
     O                 NS, NE, XVEC, YVEC, YCVEC, YNVEC, EFLAG)

C     Compute a steady-flow water-surface profile.  See IMPTRAP
C     for definition of dummy arguments.

      IMPLICIT NONE

      INTEGER  NS, NE, NMAX, EFLAG, STDOUT
      REAL*8 XS, XE, YS, EPS
      REAL*8 XVEC(NMAX), YVEC(NMAX), YCVEC(NMAX), YNVEC(NMAX)

      REAL*8 F, P
      EXTERNAL F, P, FIND_YCYNYM

      INCLUDE 'barrel.cmn'

C     Local variables

      INTEGER DIR, I, LOC, N, IS, IE, EXTRAP, RFLAG
      REAL*8  EPSARG, EPSF, EPSABS, LOC_XS, LOC_XE, LOC_YS, 
     A        YCRITA, YNORMA, YMAXA, HMINA, HMAXA, YSFAC, 
     B        DHLIM

C     Called program units

      INTEGER FIND_BARREL_INTERVAL
      EXTERNAL FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL


C***********************************************************************
C     Set the local values to control the computation
      EPSARG = 0.01D0*EPS
      EPSF = 0.00001D0*EPS
      EPSABS = 0.0001D0*EPS

      HMINA = 0.01D0*YS
      IF(HMINA.GT.0.001D0) HMINA = 0.0001D0
      HMAXA = 0.2D0*ABS(XE - XS)

      YSFAC = 1.0D0
      DHLIM = 2.0D0
      EXTRAP = 1
      
C     Set the direction computation.  Note that we now require that the
C     stationing increase downstream.  Thus we can infer the direction
C     of computation from the relationship between the starting and
C     ending station.  

      IF(XS.GT.XE) THEN
C       We are computing upstream on a subcritical profile. 
        DIR = -1
        NE = NMAX
        N = NE
      ELSE
C       We are computing downstream on a super critical profile.
        DIR = 1
        NS = 1
        N = NS
      ENDIF

C     Find the interval in the channel description that contains 
C     the starting position. 

       LOC = FIND_BARREL_INTERVAL(
     I                            XS)

C     LOC gives the index that is at or upstream of the starting
C     location.  Start the loop over the channel description.

      IF(DIR.GT.0) THEN
        IS = LOC
        IE = NXS - 1
      ELSE
        IS = LOC 
        IE = 1
      ENDIF

      LOC_XS = XS
      LOC_YS = YS
      DO 500 I=IS, IE, DIR
C       The left point is always set to the ups end of the
C       channel interval being computed.

        CALL  SET_BARREL_INTERVAL(
     A                            I)
         
C       Find the critical and normal depths at the starting 
C       location.  YMAXA is recomputed in this case.

        IF(I.EQ.IS) THEN
C         Provide starting values for the solution.
          YCRITA = YS
          YNORMA = YS
          CALL FIND_YCYNYM
     I                    (STDOUT, XS, EPSARG, EPSF, EPSABS, 
     M                     YCRITA, YNORMA,
     O                     YMAXA, RFLAG)
        ENDIF
        

        IF(DIR.GT.0) THEN
          LOC_XE = X_R
          IF(LOC_XE.GT.XE) LOC_XE = XE
        ELSE
          LOC_XE = X_L
          IF(LOC_XE.LT.XE) LOC_XE = XE
        ENDIF
        CALL IMPTRAP                  
     I              (STDOUT, PRIS_FLAG, LOC_XS, LOC_XE, LOC_YS, EPS,
     I               EPSARG, EPSF, EPSABS, HMINA, HMAXA, 
     I               YCRITA, YNORMA, YMAXA, P, F, FIND_YCYNYM, 
     I               EXTRAP, YSFAC, DHLIM, NMAX,
     M               N,
     O               XVEC, YVEC, YCVEC, YNVEC, RFLAG)

        IF(RFLAG.LT.1) THEN
          EFLAG = 1
          IF(DIR.GT.0) THEN
            NE = N
          ELSE
            NS = N
          ENDIF
          GOTO 1000
        ELSEIF(RFLAG.EQ.2) THEN
C         Stopped at critical depth.
          IF(DIR.GT.0) THEN
            NE = N
          ELSE
            NS = N
          ENDIF
          GOTO 1000
        ENDIF

C       Update values for the next segment. 
        IF(I.NE.IE) THEN
          LOC_XS = LOC_XE
          LOC_YS = YVEC(N)
          YCRITA = YCVEC(N)
          YNORMA = YNVEC(N)
C          N = N + DIR
        ENDIF

500   CONTINUE
      IF(DIR.GT.0) THEN
        NE = N
      ELSE
        NS = N
      ENDIF
1000  CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE INPUT_BARREL(
     I                        STDIN, STDOUT,
     O                        SFAC,  EFLAG)

C     Input the barrel description.  

      IMPLICIT NONE

      INTEGER STDIN, STDOUT, EFLAG

      REAL*8 SFAC


      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'barrel.cmn'

C     Local variables
      
      INTEGER MAXN
      PARAMETER (MAXN=5)

      INTEGER NODE, OLD_NODE, TABN, N
      INTEGER ITEM_START(MAXN), ITEM_END(MAXN)
      
      REAL*8 STAT, ELEVATION, SLOPE

      CHARACTER LINE*80, NODEID*16, JUST*5

C     Called program units
      REAL FMXARG, GETD
      CHARACTER GET_TABID*16
      EXTERNAL inline, CHKTAB, FMXARG, GETD, GET_ITEM_LIMITS,
     A         READ_BARREL_ITEMS, GET_TABID
C     *************************FORMATS**********************************
1     FORMAT(5X,F10.0)
2     FORMAT(I5,1X,A8,1X,I5,2F10.0)
50    FORMAT(/,A)
51    FORMAT(/,' Station factor=',1PE12.5)
52    FORMAT(I5,1X,A8,1X,A16,1X,F10.3,F10.3)
53    FORMAT(' Node   NodeId Table Identifier    Station Elevation')

80    FORMAT(/,' *ERR:XXX* More than ',I5,' cross sections in a ',
     A          'barrel description for UFGCULV.')
C***********************************************************************
C     Get the station factor.

      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)

      READ(LINE,1) SFAC
      IF(SFAC.EQ.0.D0) SFAC = 1.D0
      WRITE(STDOUT,51) SFAC

C     Read the heading line
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      JUST = 'RIGHT'
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     N, ITEM_START, ITEM_END)

C     Output new standard heading
      WRITE(STDOUT,53)

C     Input the cross section references along the barrel.
      NXS = 0
100   CONTINUE

        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
C        READ(LINE,2) NODE, NODEID, TABN, STAT, ELEVATION
        CALL READ_BARREL_ITEMS(
     I                         STDOUT, LINE, N, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, 
     O                         NODE, NODEID, TABN, STAT, 
     O                         ELEVATION)
        IF(NODE.EQ.0) NODE = OLD_NODE + 1

        IF(NODE.GT.0) THEN
          WRITE(STDOUT,52) NODE, NODEID, GET_TABID(TABN), STAT, 
     A                     ELEVATION
          NXS = NXS + 1
          IF(NXS.GT.MXNXS) THEN
            WRITE(STDOUT,80) MXNXS
            EFLAG = 1
            NXS = NXS - 1
          ENDIF
          CALL CHKTAB
     I               (12, STDOUT, FTPNT, MFTNUM,
     M                TABN,
     O                EFLAG)
          XSEC_ADRS(NXS) = TABN
          STATION(NXS) = STAT*SFAC
          INVERT_Z(NXS) = ELEVATION
C         Set the slope category to discontinuous.
          INVERT_SLOPE_CAT(NXS) = 0
          INVERT_DZDX(NXS) = 0.D0
          IF(TABN.GT.0.AND.EFLAG.EQ.0) THEN
C           Get the maximum argument from the table
            XSEC_YMAX(NXS) = DBLE(FMXARG(TABN))
C           Get the vertical diameter
            VERT_DVEC(NXS) = GETD(TABN, STDOUT)
          ENDIF
          IF(NXS.GT.1) THEN
C           Set some of the values that depend on the nature 
C           of the channel between cross sections as infered 
C           from the information given.
            IF(XSEC_ADRS(NXS).EQ.XSEC_ADRS(NXS-1)) THEN
C             This interval is prismatic.  We have not implemented
C             continuous slopes for this command.
              CHANNEL_VARIATION(NXS-1) = 1
            ENDIF
C           Since continuous slope at cross sections has not
C           been implemented for this command, every intra-section
C           interval has a constant slope. Compute the sine and
C           cosine of the angle of inclination treating a
C           decline in elevation from ups to dns as a positive
C           angle.

            SLOPE = (INVERT_Z(NXS-1) - INVERT_Z(NXS))/
     A              (STATION(NXS) - STATION(NXS-1))
            SINE_THETA(NXS-1) = SLOPE/SQRT(1.D0 + SLOPE**2)
            COSINE_THETA(NXS-1) = 1.D0/SQRT(1.D0 + SLOPE**2)
          ENDIF
          OLD_NODE = NODE
          GOTO 100
        ENDIF
        WRITE(STDOUT,*) ' '

      RETURN
      END
        
        

     


C
C
C
      SUBROUTINE   GET_H1_VS_Y2
     I                         (ADRS,
     M                          N, 
     O                          Y2CULV, H1CULV)

C     Extract the values of Y2, depth at culvert entrance, and
C     H1, head at section 1 for the culvert, from the 2-D table
C     of Type 13.  This table has been computed and stored 
C     by a CULVERT command occurring earlier in the input. 

      IMPLICIT NONE

      INTEGER ADRS, N
      REAL*8 Y2CULV(N), H1CULV(N)

C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'

C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     ADRS - address for the table
C     N    - on entry gives the maximum number of elements in the vectors
C     Y2CULV - contains the values of entrance depth
C     H1CULV - contains the values of head at section 1
 
C     Local variables

      INTEGER I, NHU, PHU, PPFD

      REAL Y2OLD, Y2
C***********************************************************************
C     Get the address for zero upstream head, and the address for the first
C     freedop sequence.
C     head.
      PHU = ITAB(ADRS+3)
      PPFD = ITAB(ADRS+5)

C     Compute the number of upstream heads( includes zero head)
      NHU = (PPFD - PHU)/4

C     Now extract the values until two consecutive values for Y2 match
C     or the end of the list has been found.
      Y2OLD = -1.0
      DO 100 I=1,NHU
        Y2 = FTAB(PHU+2)
        IF(ABS((Y2 - Y2OLD)/Y2OLD).LE.1.E-5) THEN
C         Presume we have a match.  The current 
C         value is discarded.
          N = I-1
          GOTO 102
        ENDIF

        Y2CULV(I) = Y2
        H1CULV(I) = FTAB(PHU)
        Y2OLD = Y2
        IF(Y2OLD.EQ.0.0) THEN
          Y2OLD = 0.1
        ENDIF
        PHU = PHU + 4
100   CONTINUE
      N = NHU
    
102   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE MAKE_H1_VS_Y2_FOR_CULV(
     I                            STDOUT, GETY2,
     M                            FTPTEMP,
     O                            H1VSY2, EFLAG)

C     Make a 1-D table of H1 versus Y2 for for culvert flow
C     for use in computations for an underflow gate on a culvert.

      IMPLICIT NONE

      INTEGER EFLAG, FTPTEMP, GETY2, H1VSY2, STDOUT

C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'

C     + + +DUMMY ARGUMENT DEFINITIONS + + +

C     STDOUT - unit for standard output
C     GETY2 - address of table of type 13 containing the Y2 values
C             for the culvert.
C     FTPTEMP- temp value for the function table pointer.
C     H1VSY2- address for the resulting 1-D table.
C     EFLAG- error flag.

C     Local variables.

      INTEGER I, N

      REAL  ARGVEC(PMXNHU), F1(PMXNHU), F2(PMXNHU)

      REAL*8 Y2CULV(PMXNHU), H1CULV(PMXNHU), M(PMXNHU), SLOPE

      character*1 ADJLOC(PMXNHU)

C     + + + EXTERNAL NAMES + + +
      EXTERNAL GET_H1_VS_Y2, VLCHPP, PUT1D

C***********************************************************************
C     Extract the values from the table of Y2
      N = PMXNHU
      CALL GET_H1_VS_Y2
     I                 (GETY2,
     M                  N, 
     O                  Y2CULV, H1CULV)

C     Compute a variation limited cubic spline.  Use linear approx to 
C     derivative for the left end condition and a derivative of zero
C     for the right end condition.

      SLOPE = (H1CULV(2) - H1CULV(1))/(Y2CULV(2) - Y2CULV(1))
      CALL VLCHPP
     I           (STDOUT, N, Y2CULV, H1CULV, 1, SLOPE, 1, 0.D0,
     O            M, ADJLOC)

C     Transfer to single-precision values.
      DO 100 I=1,N
        ARGVEC(I) = Y2CULV(I)
        F1(I) = H1CULV(I)
        F2(I) = M(I)
100   CONTINUE
       
C     Store as table of type 4 using a table number of -1.
      CALL PUT1D
     I          (STDOUT, -1, 4, N, ARGVEC, F1, F2,
     M           FTPTEMP,
     O           H1VSY2)

      RETURN
      END
C
C
C
      SUBROUTINE   UFGCULV
     I                    (GRAV, STDIN, STDOUT, STDTAB, FTP,
     M                     EFLAG, TABDIR)
 
C     + + + PURPOSE + + +
C     Compute the function tables needed to define the flow
C     for an underflow gate on the upstream face of a 
C     culvert.
 
      IMPLICIT NONE

C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, STDIN, STDOUT, STDTAB
      INTEGER TABDIR(*)
      REAL GRAV
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     GRAV   - value of acceleration due to gravity
C     STDIN  - Fortran unit number for user input file
C     STDOUT - Fortran unit number for user output and messages
C     STDTAB - Fortran unit number for output of function tables
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
C     TABDIR - Table directory to remember table numbers
 
C     + + + COMMON BLOCKS + + +
      INCLUDE 'arsize.prm'
      INCLUDE 'ftable.cmn'
      INCLUDE 'epscom.cmn'
      INCLUDE 'ufgate.cmn'
      INCLUDE 'ufgate_d.cmn'
      INCLUDE 'barrel.cmn'
 
C     + + + LOCAL VARIABLES + + +
      INTEGER CCTAB, EF, ELFLAG, I, IHG, IHGERR, IHU, IPFD, J,
     A        JBASE, N, NFRAC, NHG, NHU, NRMS, TAB, TABGT, TABLT, TPTAB,
     B        GETQ, GETY2, NTAB, FTPTEMP, H1VSY2, SOFLAG,
     C        FCFLAG, COLWIDTH, NITEM, IDLEN

      INTEGER MAXN
      PARAMETER (MAXN=4)
      INTEGER TAB2D(PMXNHG), ITEM_START(MAXN), ITEM_END(MAXN)
      REAL ANGLE(PMXNHG), BIGERR, BRKPFD, CCVAL(PMXNHG), CONCC,
     A     DE1TO4, DROP, DT, FDROP, FDVEC(PMXNHU),
     B     FINPOW, FWFOTR, H1, H1FOLL, H1FWUL, H3, H4, H4F,
     C     HDATUM, HERR, HGVEC(PMXNHG),
     D     HSTUFF(PMXNHG,3), HUMAX, HUMIN, HUVEC(PMXNHU), LIMPFD,
     E     LIPREC, LIPVEC(PMXNHG), MAXHU, MINHU, MINPFD, MNHVEC(PMXNHG),
     F     OFFSET, OLDHG, P, PFDTMP(PMXFRC), PFDVEC(PMXFRC), POW,
     G     POWER, PREC, Q, QFREE, QHAT, QMAT(PMXNHU,PMXFRC), RERR,
     H     RG, RGFWUL, RMS, WORK(PMXFRC), XBRK(PMXFRC),  Y2FW,
     I     Y4F, Y4SW, Z1FWUL, ZW4, ZW4F,
     J     DERIV, Z1, zrhufd

      REAL*8 SFAC

      INTEGER NS, NE, RFLAG, FREE
      REAL DQED, DQEU, TP, HDATUM_FOR_H1VSY2, AFULL, KFULL, JFULL,
     A      DK, B, DB, YVCOLD, XE_R, AFLUX, A3, QTT
      REAL*8  XS, XE, YS, YCRIT, YNORM, YMAX, STATION_1,
     A         CD_CULV, H1FWUL_D, YVC, Y2HAT, YB_D,
     B        FW_CC_SO, Y4SW_D, 
     C        H1FOLL_D, easting, northing

      REAL*8 DBLE

      CHARACTER CHAR6*6, FTYPE*2, HGLAB*10, LAB2*50, LABEL*50, LINE*80,
     A          CQ*11, CHAR8*8, TABID*16, KEY*16, APPTABID*16, 
     B          DEPTABID*16, CCTABID*16, JUST*5, IDOUT*32, ID*16,
     c   zone*8, hgrid*8, vdatum*8, unitsys*8, basis*8
 
C     + + + INTRINSICS + + +
      INTRINSIC ABS, FLOAT, LOG, SQRT
 
C     + + + EXTERNAL FUNCTIONS + + +
      INTEGER FIND_BARREL_INTERVAL, LENSTR
      REAL FINDCC,  GETHDD
      REAL*8 E, RHSE, FIND_FOLL_RESID
      CHARACTER PUT10*10, GET_TABID*16
 
C     + + + EXTERNAL NAMES + + +
      EXTERNAL CHKTAB, FDROOT, FINDCC, FNDELV, inline,
     A         INVTSE, LKTJ, LSTOPF, REGFAL, 
     B         TABCHK, TWDOUT, XLKT22, PUT10,
     C         MAKE_H1_VS_Y2_FOR_CULV, LKTAB, E, RHSE, FIND_YCYNYM,
     D         FIND_BARREL_INTERVAL, SET_BARREL_INTERVAL, TDLK13,
     E         FNDSTA, GETHDD, FIND_FW_CC_SO, FIND_FOLL_RESID,
     D         GET_TABID, LENSTR, GET_INTERNAL_TAB_NUMBER,
     E         GET_UFGC_TABIDS, SET_UFGC_TABIDS
 
C     + + + INPUT FORMATS + + +
 1    FORMAT(7X,I5,6X,I5,7X,I5)
 2    FORMAT(6X,A)
 4    FORMAT(7X,I5)
 6    FORMAT(7X,I5)
 8    FORMAT(9X,F10.0)
 10   FORMAT(10X,F10.0)
 12   FORMAT(3X,F10.0)
 14   FORMAT(6X,I5)
 16   FORMAT(7X,F10.0)
 18   FORMAT(7X,F10.0)
 19   FORMAT(9X,F10.0)
 20   FORMAT(8X,F10.0)
 21   FORMAT(8X,F10.0)
 22   FORMAT(10X,F10.0)
 
C     + + + OUTPUT FORMATS + + +
 50   FORMAT(/,' Table id= ',A,' for type 15 table for underflow',
     A     ' gate/culvert.')
 51   FORMAT(/,' *ERR:XXX* Culvert flow TabId= ',A,' given but GETY2',
     A         ' value is missing.')
 52   FORMAT(/,' Label=',A50)
 53   FORMAT(/,' *ERR:XXX* TabId for depth at culvert exit= ',A,
     A       ' given but GETQ value is missing.')
 54   FORMAT(/,' Approach section table id=',A)
 55   FORMAT(/,' *ERR:XXX* GETQ and GETY2 values are missing.')
 56   FORMAT(/,' Departure section table id= ',A)
 57   FORMAT(/,' Table id''s for culvert flow= ',A,
     A       ' and for entrance depth= ',A)
 58   FORMAT(/,' Elevation of gate sill=',F10.3)
59    FORMAT(/,' FC flow cannot match FW flow at the FW limit.'/,
     A       ' May be a discontinuity in flow at the FW limit.')
 60   FORMAT(/,' Total gate opening width=',F10.3)
61    FORMAT(/,' FC flow matches FW flow at FW limit with CC=',F10.3)
 62   FORMAT(/,' User input of discharge coeff. read but ignored',
     A       ' in UFGCULV.')
63    FORMAT(/,' Coefficient of discharge for FC, FO, and SO flows='
     A     ,F10.4)      
 64   FORMAT(/,' Contraction coefficient table id= ',A)
65    FORMAT(/' Critical depth/vert. dia. ratio>=',F10.5,
     A        ' for full-barrel flow.',/,
     B    ' Critical depth/vert. dia. ratio<=',F10.5,
     C    ' upper limit for free-surface flow.')
 66   FORMAT(/,' Minimum partial free drop=',F10.5)
 67   FORMAT(/,' Partial free drop at power breakpoint=',F10.5)
 68   FORMAT(/,' Limiting partial free drop=',F10.3)
 69   FORMAT(/,' Final local power for partial free drops=',F10.3)
 70   FORMAT(/,'Two-D table computations for gate opening=',F8.3,
     A    '  Two-D table number=',I6)
71    FORMAT(5X,'Head at section 1 for FW/FO-FC boundary=', F9.4,/,
     A       5X,'Head-relative gate opening at FW/FO-FC boundary=',F9.4)
 72   FORMAT(/,'Upstream head=',F9.4,
     A  ' Elevation=',F10.4,/,2X,'Depth at section 1=',F8.4,
     B  ' Gate opening=',F8.4)
73    FORMAT(6X,'Processing gate opening=',F8.3)
 74   FORMAT(/,
     A '  Partial  Drop    Head    Head   Flow Cont.  Discharge  Local',
     B '  Energy',/,
     C '   free    sect.   sect.   sect.  type coef.             power',
     D '   loss',/,
     E '   drop    1->4     3       4           Cc  ',17X,'    1->4',/,
     F ' --------  ------  ------  ------  ---  ---- ----------- -----',
     G '  ------')
 75   FORMAT(1X,F8.4,F8.3,A8,F8.3,3X,A2,A6,1X,A11,F6.3,F8.3)
 76   FORMAT(1X,F8.4,F8.3,F8.3,F8.3,3X,A2,A6,1X,A11,6X,F8.3)
 77   FORMAT(/,' *WRN:592* Please review results.  One or more',
     A        ' cases with energy gain found.')
 78   FORMAT(1X,F8.4,F8.3,A8,F8.3,3X,A2,A6,1X,A11,6X,F8.3)
 79   FORMAT(/,' Free weir to free orifice transition fraction=',F8.2)
 80   FORMAT(/,' Maximum upstream head=',F8.2)
 81   FORMAT(/,' Minimum non-zero upstream head=',F8.2)
 82   FORMAT(/,' Linear interpolation precision=',F8.3)
83    FORMAT(/,' Flow at the FW limit=',F10.2)
 86   FORMAT(/,' Maximum relative error=',F6.3,' Gate opening=',
     A     F8.4,/,'   Upstream head=',F9.4,' Partial free drop=',
     B     F8.5)
 87   FORMAT(/,' *ERR:607* TABLE# <= 0')
 88   FORMAT(/' Root-mean-squared error=',F6.3,' N in sample=',I5)
 89   FORMAT('  Processing UFGCULV TabId= ',A)
 90   FORMAT(/,' *ERR:712* Gate sill elevation=',F10.3,' <  elevation',
     A       ' of floor of',/,11X,' departure reach=',F10.3)
 91   FORMAT(/,' *ERR:713* Gate opening width=',F10.3,' <= 0.0')
 92   FORMAT(/,' *ERR:714* Approach loss Cd <=0.0 or > 1.0')
 93   FORMAT(/,' *ERR:715* Contraction coef.=',F10.3,' <= 0.0 or > 1.0')
 94   FORMAT(/,' *ERR:716* In UFGATE: FLAG=',I3,' No solution',A)
 95   FORMAT('TABID=',A,/,'TYPE=  -15',/,'REFL=',7X,'0.0 LABEL=',A50,
     B /,'   OPENING Table Identifier    H1FWULR   H4FWULR H4SWSOMDR')
 96   FORMAT(F10.3,1X,A16,1X,3F10.6)
 97   FORMAT(/,' *BUG:XXX* In UFGATE no root in 16 tries for:',A)
 98   FORMAT(/,' *ERR:717* Gate sill elevation=',F10.3,' <  elevation',
     A       ' of floor of',/,11X,' approach reach=',F10.3)
 99   FORMAT(/,' *ERR:635* Interpolation precision tables missing.',
     A   ' Check for version',/,11X,'of file: TYPE5.TAB.')
C***********************************************************************
      EPS_D = 0.0005D0
      EPSARG_D = 0.01D0*EPS_D
      EPSF_D = 0.00001D0*EPS_D
      EPSABS_D = 0.0001D0*EPS_D
     
      STDOUT_D = STDOUT

C     Establish the function for piezometric level for a full conduit
      FTPTEMP = FTP 
     
      CALL CONSTRUCT_YOVERD(
     I                      STDOUT, 
     M                      FTPTEMP,
     O                      ADRS_YOVERD)

C     Define the linear interpolation precision tables.
      KEY = '10001'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, KEY,
     M                             EFLAG,
     O                             TABLT)
      TABLT = FTPNT(TABLT)
      KEY = '10002'
      CALL GET_INTERNAL_TAB_NUMBER
     I                            (STDOUT, KEY,
     M                             EFLAG,
     O                             TABGT)
      TABGT = FTPNT(TABGT)
      IF(TABLT.LT.1.OR.TABGT.LT.1) THEN
        WRITE(STDOUT,99) 
        EFLAG = 1
        RETURN
      ENDIF
      G = GRAV
      GRAV_D = DBLE(GRAV)
      TWOG = 2.*GRAV
      TWOG_D = 2.D0*GRAV_D
C      CALL inline
C     I          (STDIN, STDOUT,
C     O           LINE)
C      READ(LINE,1,ERR=991) TAB, GETQ, GETY2
      CALL GET_UFGC_TABIDS(STDIN, STDOUT,
     M                     EFLAG)
      CALL SET_UFGC_TABIDS(
     O                     TAB, GETQ, GETY2,
     M                    EFLAG)
      WRITE(STDOUT,50) GET_TABID(TAB)
      WRITE(*,89) GET_TABID(TAB)
 
      CALL TABCHK
     I           (STDOUT, PMXTAB,
     M            TAB, TABDIR, EFLAG)
 
C     Check the GETQ and GETY2 table numbers.  If positive, the tables
C     must already exist in the FTAB/ITAB system.  Also both must be 
C     positive, that is, both tables must be specified.
      IF(GETQ.GT.0) THEN
        IF(GETY2.LE.0) THEN
C         Table for depth at culvert entrance missing.
          WRITE(STDOUT,51) GET_TABID(GETQ)
          EFLAG = 1
        ELSE
C         Both table numbers are positive.  Check to see if the tables
C         exist.
          WRITE(STDOUT,57) GET_TABID(GETQ), GET_TABID(GETY2)
          CALL CHKTAB
     I               (6, STDOUT, FTPNT, MFTNUM,
     M                GETQ,
     O                EFLAG)
          CALL CHKTAB
     I               (6, STDOUT, FTPNT, MFTNUM,
     M                GETY2,
     O                EFLAG)
        ENDIF
      ELSE
        IF(GETY2.LE.0) THEN
C         Both tables missing.
          WRITE(STDOUT,55)
          EFLAG = 1
        ELSE     
          WRITE(STDOUT,53) GET_TABID(GETY2)
C         Table for culvert flow is missing.
          EFLAG = 1
        ENDIF
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
      READ(LINE,2,ERR=991) LAB2
      WRITE(STDOUT,52) LAB2
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'APPTAB',
     O                EFLAG, APPTABID, APPTAB)
      WRITE(STDOUT,54) APPTABID
 

      IF(APPTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        EF = 0
        TPTAB = APPTAB
        CALL CHKTAB
     I             (12, STDOUT, FTPNT, MFTNUM,
     M              APPTAB,
     O              EF)
        IF(EF.EQ.0) THEN
          APPTAB_D = APPTAB
C         Find the invert elevation and station  from the table
          CALL FNDELV
     I               (TPTAB, STDOUT,
     O                EFLAG, Z1B)
          Z1B_D = DBLE(Z1B)

          CALL FNDSTA
     I                (TPTAB, STDOUT,
     O                 EFLAG, TP)
          STATION_1 = DBLE(TP)

        ELSE
          EFLAG = 1
        ENDIF
      ENDIF
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'DEPTAB',
     O                EFLAG, DEPTABID, DEPTAB)
      WRITE(STDOUT,56) DEPTABID
 
      IF(DEPTAB.LE.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSE
        EF = 0
        TPTAB = DEPTAB
        CALL CHKTAB
     I             (1, STDOUT, FTPNT, MFTNUM,
     M              DEPTAB,
     O              EF)
        IF(EF.EQ.0) THEN
C         Find the invert elevation from the table
          CALL FNDELV
     I               (TPTAB, STDOUT,
     O                EFLAG, Z4B)
        ELSE
          EFLAG = 1
        ENDIF
      ENDIF
 
C     For the moment make Z3B the same as Z4B.  Note that
C     the table called the departure table will give the
C     cross section that defines the upstream hydrostatic pressure
C     force on the control volume downstream of the underflow
C     gate.  An optional table, to be added later will give
C     the table for computing the momentum flux exiting from the
C     control volume.
 
      Z3B = Z4B
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,8,ERR=991) Z2B
 
      WRITE(STDOUT,58)  Z2B
 
      IF(Z2B.LT.Z4B) THEN
        WRITE(STDOUT,90) Z2B, Z4B
        EFLAG = 1
      ENDIF
      IF(Z2B.LT.Z1B) THEN
        WRITE(STDOUT,98) Z2B, Z1B
        EFLAG = 1
      ENDIF
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,10,ERR=991)  BG
      WRITE(STDOUT,60) BG
 
      IF(BG.LE.0.0)  THEN
        WRITE(STDOUT,91) BG
        EFLAG = 1
      ENDIF
      BG_D = BG
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,12,ERR=991) CD
      WRITE(STDOUT,62)  
C      IF(CD.LE.0.0.OR.CD.GT.1.0) THEN
C        WRITE(STDOUT,92) CD
C        EFLAG = 1
C      ENDIF
 
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      CALL READ_TABID
     I               (STDOUT, LINE, 'CCTAB',
     O                EFLAG, CCTABID, CCTAB)
      WRITE(STDOUT,64) CCTABID
      IF(CCTAB.LT.0) THEN
        WRITE(STDOUT,87)
        EFLAG = 1
      ELSEIF(CCTAB.GT.0) THEN
        CALL CHKTAB
     I             (2, STDOUT, FTPNT, MFTNUM,
     M              CCTAB,
     O              EFLAG)
      ENDIF
 
      CCTAB_D = CCTAB

C     Get the size of the transition between FW and FO flow.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,19,ERR=991)  FWFOTR
      IF(FWFOTR.EQ.0.0) FWFOTR = 0.1
      WRITE(STDOUT,79) FWFOTR


C     Get the maximum upstream head
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,20,ERR=991)  MAXHU
      WRITE(STDOUT,80) MAXHU
 
C     Get the smallest non-zero upstream head
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,21,ERR=991)  MINHU
      IF(MINHU.EQ.0.0) MINHU = 0.1
      WRITE(STDOUT,81) MINHU
 
C     Get the global linear interpolation precisions
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,22,ERR=991)  LIPREC
      IF(LIPREC.EQ.0.0) LIPREC = 0.02
      WRITE(STDOUT,82) LIPREC
 
C     Process the gate opening table.
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,'(A)') LINE
 
 
      JUST = 'RIGHT'
      CALL GET_ITEM_LIMITS(
     I                     STDOUT, LINE, MAXN, JUST,
     O                     NITEM, ITEM_START, ITEM_END)
C     Set the  column width
      COLWIDTH =  ITEM_END(2) - ITEM_START(2) + 1
 
 
      OLDHG = 0.0
      NHG = 1
      IDLEN = 0
 100  CONTINUE
        CALL inline
     I            (STDIN, STDOUT,
     O             LINE)
        CALL READ_UFGATE_ITEMS(
     I                         STDOUT, LINE, NITEM, ITEM_START,
     I                         ITEM_END,
     M                         EFLAG, IDLEN,  
     O            HGVEC(NHG), TAB2D(NHG), CCVAL(NHG), ANGLE(NHG))
C        READ(LINE,'(F10.0, I10, 4F10.0)') HGVEC(NHG), TAB2D(NHG),
C     A                      CCVAL(NHG), ANGLE(NHG)
C     B                      MNHVEC(NHG), LIPVEC(NHG)
        MNHVEC(NHG) = 0.0
        LIPVEC(NHG) = 0.0
        IF(HGVEC(NHG).LE.OLDHG) THEN
C         Input complete.
          NHG = NHG - 1
        ELSE
          OLDHG = HGVEC(NHG)
          IF(MNHVEC(NHG).EQ.0.0) THEN
            MNHVEC(NHG) = MINHU
          ENDIF
          IF(LIPVEC(NHG).EQ.0.0) THEN
            LIPVEC(NHG) = LIPREC
          ENDIF
          ID = GET_TABID(TAB2D(NHG))
          IDOUT = ' '
          IDOUT(COLWIDTH-IDLEN+1:COLWIDTH) = ID(1:IDLEN)
          WRITE(STDOUT,
     A        '(F10.3,A,F10.3,F10.1,F10.2,F10.3)') HGVEC(NHG),
     B           IDOUT(1:COLWIDTH), CCVAL(NHG), ANGLE(NHG)
C     C           MNHVEC(NHG), LIPVEC(NHG)
          CALL TABCHK
     I               (STDOUT, PMXTAB,
     M                TAB2D(NHG), TABDIR, EFLAG)
          IF(CCTAB.EQ.0.0) THEN
C           If no table is given for contraction coefficient, then
C           a value must be given by the user.
            IF(CCVAL(NHG).LE.0.0.OR.CCVAL(NHG).GT.1.0) THEN
              WRITE(STDOUT,93) CCVAL(NHG)
               EFLAG = 1
            ENDIF
          ENDIF
          NHG = NHG + 1
          GOTO 100
        ENDIF
 
      IF(NHG.LT.2) THEN
        WRITE(STDOUT,98)
        EFLAG = 1
      ENDIF
 

C     Now get the description of the barrel.  The descriptive values
C     are stored in common block in the include file UFGC.COM  
      CALL INPUT_BARREL(
     I                  STDIN, STDOUT,
     O                  SFAC, EFLAG)

C     Adjust the station at the approach section now that SFAC is
C     available.

      STATION_1 =  STATION_1*SFAC



C     INPUT THE FACTORS CONTROLLING THE DISTRIBUTION OF DROPS
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,'(A)') LINE
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) MINPFD
      WRITE(STDOUT,66) MINPFD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) BRKPFD
      WRITE(STDOUT,67) BRKPFD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,18,ERR=991) LIMPFD
      WRITE(STDOUT,68) LIMPFD
 
      CALL inline
     I          (STDIN, STDOUT,
     O           LINE)
      READ(LINE,16,ERR=991) FINPOW
      WRITE(STDOUT,69) FINPOW
 
 
C     COMPUTE THE PROPORTIONS OF FREE DROP
 
      POW = 0.5
      OFFSET = 0.0
      CALL LSTOPF 
     I           (STDOUT, TABLT, TABGT, POW, OFFSET, MINPFD, BRKPFD,
     I            LIPREC, PMXFRC,
     O            N, XBRK, EFLAG)
 
      WORK(1) = 0.0
      DO 195 I=1,N
        WORK(I+1) = XBRK(I)
 195  CONTINUE
      NFRAC = N + 1
C      OFFSET = 0.9*BRKPFD
      OFFSET = 0.0
      CALL LSTOPF 
     I           (STDOUT, TABLT, TABGT, FINPOW, OFFSET, BRKPFD, LIMPFD,
     I            LIPREC, PMXFRC,
     O            N, XBRK, EFLAG)
      DO 196 I=2,N
        WORK(NFRAC+I-1) = XBRK(I)
 196  CONTINUE
      NFRAC = NFRAC + N - 1
 
C      WORK(N+2) = 0.95
C      WORK(N+3) = 0.975
C      WORK(N+4) = 0.98
C      WORK(N+5) = 0.99
      WORK(NFRAC+1) = 1.0
      NFRAC = NFRAC + 1
 
      PFDVEC(1) = WORK(1)
      PFDVEC(2) = WORK(2)
      J = 2
      DO 201 I=3,NFRAC-1
        J = J + 1
        PFDVEC(J) = 0.5*(WORK(I) + WORK(I-1))
        J = J + 1
        PFDVEC(J) = WORK(I)
 201  CONTINUE
      PFDVEC(J+1) = WORK(NFRAC)
      NFRAC = J + 1
C      DO 202 I=1,NFRAC
C        WRITE(STDOUT,*) ' I=',I,' PFDVEC(I)=',PFDVEC(I)
C202   CONTINUE
 
      IF(EFLAG.NE.0) RETURN
 
C     INPUT OF DATA COMPLETE.  BEGIN THE COMPUTATIONS.
      HDATUM = Z2B
      HDATUM_D = DBLE(HDATUM)

C     Prepare a 1-D table giving the head at section 1 as a 
C     function of the depth in the entrance to the culvert
C     for culvert flow. 
      CALL MAKE_H1_VS_Y2_FOR_CULV(
     I                            STDOUT, GETY2,
     M                            FTPTEMP,
     O                            H1VSY2, EFLAG)

C     Set the head datum for the heads in the table, H1VSY2.
C     They may differ from the head here.

      HDATUM_FOR_H1VSY2 = GETHDD(GETY2)

C     Lookup the characteristics of the full barrel.  We assume 
C     a prismatic barrel. 
      VERT_D = VERT_DVEC(NXS)
      TP = SNGL(VERT_D)
      CALL XLKT21
     I           (XSEC_ADRS(NXS),
     M            TP,
     O            AFULL, TT, DTT, JFULL, KFULL, DK, B, DB)
      AFULL_D = DBLE(AFULL)
      JFULL_D = DBLE(JFULL)
      KFULL_D = DBLE(KFULL)

C     Find the elevation of the exit invert
      ZBEX_D = INVERT_Z(NXS)
 
      Z3B = ZBEX_D

      BIGERR = 0.0
      RMS = 0.0
      NRMS = 0
      OLDHG = 0.0
      DO 1000 IHG=1,NHG
        ELFLAG = 0
        HG = HGVEC(IHG)
        HG_D = DBLE(HG)
        WRITE(HGLAB,'(''Hg='',F6.3,'':'')') HG
        LABEL = HGLAB // LAB2
        MINHU = MNHVEC(IHG)
        PREC = LIPVEC(IHG)
        CONCC = CCVAL(IHG)
        CONCC_D = CONCC

        WRITE(*,73) HG
        WRITE(STDOUT,70) HG, TAB2D(IHG)
 
        AG = BG*HG
 
C       Find the free  flow(FW flow ) that exists when critical depth
C       in the gate opening equals the gate opening.
        Y2 = HG
        CALL LKTAB
     I            (H1VSY2, Y2, 0,
     O               H1FWUL, NTAB, DERIV)

        Z1FWUL = H1FWUL + HDATUM_FOR_H1VSY2
        H1FWUL = Z1FWUL - HDATUM
        HSTUFF(IHG,1) = H1FWUL/HG
        H1FWUL_D = DBLE(H1FWUL)
C       Compute the gate-opening ratio for the free-weir flow limit.
        RGFWUL = HG/H1FWUL
 
        WRITE(STDOUT,71) H1FWUL, RGFWUL

        CALL TDLK13
     I             (STDOUT, GETQ, 13, 0, 0.0, 0.0, Z1FWUL,
     I              HDATUM_FOR_H1VSY2,
     O              Q, DQED, DQEU, FREE)
        QD = Q
        WRITE(STDOUT,83)  Q 

C       Find the tailwater elevation at the limit of FW flow 
C       at the upper limit for FW flow at section 1, that is, H1FWUL.

        CALL TD13_FDROP
     I               (STDOUT, GETY2, 13, H1FWUL,
     I                HDATUM_FOR_H1VSY2,
     O                FDROP)

        HSTUFF(IHG,2) = (H1FWUL - FDROP)/HG

C       Find the tailwater level that causes the water to contact
C       the gate lip when the head at section 1 is midway between
C       HG and H1FWUL.  
        Z1 = 0.5*(HG + H1FWUL) + HDATUM

        CALL FIND_Y4SW(
     I                 STDOUT, Z1, Z2B, Z4B, GETY2, HG_D,
     I                 EPSARG_D, EPSF_D,
     O                 Y4SW_D)
                     
        HSTUFF(IHG,3) = (Y4SW_D + Z4B - HDATUM)/HG

C       Find the elevation at the vc for the current gate opening
        CALL LKT_ZB(
     I              1, STATION(1)+HG_D,
     O               ZBVC_D)

 
C       Determine if FC flow exists and find its limits.  FC flow
C       is flow at the barrel exit that drowns flow under the gate.
C       The flow can be critical at the exit or the barrel can be
C       full and barrel resistance can cause the flow under the
C       gate to be drowned. 

C       Find the the contraction coefficient that causes the flow
C       at the FW limit to be matched by the SO equation. 

C       Find the critical depth that will support full barrel 
C       flow
        CALL FIND_FB_YCRIT
     I                    (STDOUT,
     O                     FB_RATIO) 
C       Set an increment above the full-barrel ratio for critical 
C       depth for overlap between full-barrel computations and
C       free-surface computations. 
        
        FS_RATIO = FB_RATIO + 0.125D0*(1.D0 - FB_RATIO)

        WRITE(STDOUT,65) FB_RATIO, FS_RATIO

C       Find the critical and normal depths. Make starting guesses.
        
        YCRIT = DBLE(Y2)
        YNORM = 0.5D0*VERT_D
        CALL GET_YCYNYM
     I                  (STDOUT, STATION(NXS), EPSARG_D, EPSF_D, 
     I                   EPSABS_D,
     M                   YCRIT, YNORM,
     O                   YMAX, RFLAG)
        IF(RFLAG.LT.0) THEN
          WRITE(STDOUT,*) ' Problem with critical depth comp. FC flow'
          STOP 'Abnormal stop. Errors found.'
        ENDIF

C        WRITE(STDOUT,*) ' YCRIT=',YCRIT
C        WRITE(STDOUT,*) ' YNORM=',YNORM

        IF(YNORM.GT.YCRIT.OR.YNORM.LT.0.D0) THEN
          YS = YCRIT*1.0001D0
          XS = STATION(NXS)
          XE = STATION(1) + HG_D
          CALL SFWSP(
     I               STDOUT, XS, XE, YS, EPS_D,
     I               E, RHSE, FIND_YCYNYM, NMAX,
     O               NS, NE, XVEC, YVEC, YCVEC, YNVEC, EFLAG)

          Y2HAT = YVEC(NS)
        ELSE
          Y2HAT = Y2
        ENDIF

C       Current assumption is that FC only occurs on mild or horizontal
C       conduit slopes.  Adverse slopes not supported.  We make also
C       try for steep slopes to define CD_CULV.
C        WRITE(STDOUT,*) ' Seeking FC limit: YCRIT=',SNGL(YCRIT)
C        WRITE(STDOUT,*) ' YNORM=',SNGL(YNORM)

        CALL FIND_FW_CC_SO(
     I                  STDOUT, EPSARG_D, EPSF_D, EPSABS_D, APPTAB, 
     I        GRAV_D, DBLE(Z1FWUL), Z1B_D, STATION_1, 
     I        HG_D, BG_D, Y2HAT, 
     O        CD_CULV,  FW_CC_SO, SOFLAG, RFLAG)
        IF(RFLAG.NE.0) THEN
C          WRITE(STDOUT,*) ' Problem finding FW_CC_SO'
        
        ENDIF 
        IF(FW_CC_SO.EQ.0.0) THEN
          WRITE(STDOUT,59) 
        ELSE
          WRITE(STDOUT,61) FW_CC_SO
        ENDIF

        WRITE(STDOUT,63) CD_CULV
        CD_D = CD_CULV

        IF(YNORM.GT.YCRIT.OR.YNORM.LT.0.D0) THEN
C         Seek the lower limit of FO flow assuming the standard 
C         contraction coefficient

          H1FOLL_D = H1FWUL_D
          CALL FIND_FOLL(
     I                   STDOUT, CD_CULV, DBLE(MAXHU),
     M                   H1FOLL_D,
     O                   RFLAG)
        ELSE
C         Barrel slope is steep.  FO assumed to start immediately
          H1FOLL_D = 0.D0

        ENDIF



        H1FOLL = SNGL(H1FOLL_D)
        IF(H1FOLL_D.EQ.0.D0) THEN
C         Do transition over a user-given fraction of the gate
C         opening
          H1FOLL_D = H1FWUL_D + FWFOTR*HG_D
          H1FOLL = H1FOLL_D
C         We cannot compute flows between H1FWUL and H1FOLL 
          FCFLAG = 0
        ELSE
          FCFLAG = 1
        ENDIF
 
C       Assign the upstream heads for this gate opening.  The
C       ranges are: MINHU to HG, HG to H1FWUL, H1FWUL to H1FOLL,
C       and H1FOLL to MAXHU. 
        POWER = 1.5
        OFFSET = 0.0
        IF(MINHU.GE.HG) THEN
          MINHU = 0.5*HG
        ENDIF
        HUMIN = MINHU
        CALL LSTOPF 
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, HUMIN, HG,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        DO 210 I=1,N
          HUVEC(I) = XBRK(I)
 210    CONTINUE
        NHU = N
 
        CALL LSTOPF 
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, HG, H1FWUL,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        IF(N.EQ.2) THEN
          XBRK(3) = XBRK(2)
          XBRK(2) = 0.5*(XBRK(1) + XBRK(3))
          N = 3
        ENDIF
        DO 212 I=2,N
          NHU = NHU + 1
          HUVEC(NHU) = XBRK(I)
 212    CONTINUE
        OFFSET = -CC*HG
C       Flow varies as the .5 power of the head - offset.  However,
C       the drop to free flow varies more nearly like the 1.5 power.
C       The spacing for the 1.5 power is smaller than for .5 power;
C       so we use 1.5 power.
        POWER = 1.5
        IF(FCFLAG.EQ.1) THEN
          CALL LSTOPF 
     I               (STDOUT, TABLT, TABGT, POWER, OFFSET, H1FWUL, 
     I                H1FOLL, PREC, PMXNHU,
     O                N, XBRK, EFLAG)
          IF(N.EQ.2) THEN
            XBRK(3) = XBRK(2)
            XBRK(2) = 0.5*(XBRK(1) + XBRK(3))
            N = 3
          ENDIF
        ELSE
C         Can only compute at the end points here.
          N = 2
          XBRK(1) = H1FWUL
          XBRK(2) = H1FOLL
        ENDIF
        DO 213 I=2,N
          NHU = NHU + 1
          HUVEC(NHU) = XBRK(I)
 213    CONTINUE
        IF(OLDHG.GT.0.0) THEN
          HUMAX = MAXHU*HG/OLDHG
        ELSE
          HUMAX = MAXHU
        ENDIF
        IF(HUMAX.LE.H1FOLL) THEN
C         The maximum requested head is less than the head at section
C         1 at the lower limit of FO flow.
          HUMAX = H1FOLL*1.10
        ENDIF
        OLDHG = HG
        CALL LSTOPF 
     I             (STDOUT, TABLT, TABGT, POWER, OFFSET, H1FOLL, HUMAX,
     I              PREC, PMXNHU,
     O              N, XBRK, EFLAG)
        DO 214 I=2,N
          NHU = NHU + 1
          HUVEC(NHU) = XBRK(I)
 214    CONTINUE
 
        
C       Now compute the free flow for each of the upstream heads
        DO 900 IHU=1,NHU
          Z1 = HUVEC(IHU) + HDATUM
          H1 = Z1 - Z2B
          H1_D = DBLE(H1)
          Y1 = H1 + Z2B - Z1B
          RG = HG/H1
          IF(RG.GT.1.0) RG = 1.0
          CC = FINDCC(RG, CONCC, CCTAB)
          
          WRITE(STDOUT,72) HUVEC(IHU), H1 + Z2B, Y1, HG
          WRITE(STDOUT,74)
C         Set the flow for zero partial free drop to 0.0
          QMAT(IHU,1) = 0.0
C         Find conditions at section 1
          CALL XLKT22
     I               (APPTAB,
     M                Y1,
     O                A1, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA1,
     O                DALPHA, QCT)
          A1_D = DBLE(A1)
          ALPHA1_D = DBLE(ALPHA1)

          IF(H1.LE.H1FWUL) THEN
C           Flow is free weir flow. 
            FTYPE = 'FW'

            CALL TDLK13
     I                 (STDOUT, GETQ, 13, 0, 0.0, 0.0, Z1,
     I                  HDATUM_FOR_H1VSY2,
     O                  QFREE, DQED, DQEU, FREE)
            CALL TDLK13
     I                 (STDOUT, GETY2, 13, 0, 0.0, 0.0, Z1,
     I                  HDATUM_FOR_H1VSY2,
     O                  Y2FW, DQED, DQEU, FREE)

C            WRITE(STDOUT,*) ' H1=',H1,' FW Q=', QFREE,' Y2FW=',Y2FW
 
C           Now find the level at section 4 that defines
C           the limit of free weir flow, Y4F. Define the constant
C           values.

            CALL TD13_FDROP
     I                   (STDOUT, GETY2, 13, H1,
     I                    HDATUM_FOR_H1VSY2,
     O                    FDROP)

C           Submerged weir flow can transition to submerged orifice
C           flow if the water level at section 1 is above the
C           elevation of the gate lip.
            IF((H1 - HG)/HG.GT.5.E-5) THEN
C            IF(Y1 + Z1B.GT.HG+Z2B) THEN
C              WRITE(STDOUT,*) ' H1=',H1,' H1FWUL=',H1FWUL
              IF((H1 - H1FWUL)/H1FWUL.LE.-EPSF) THEN
C               Find the tailwater level that causes the water to
C               contact the lip. 
                CALL FIND_Y4SW(
     I                       STDOUT, Z1, Z2B, Z4B, GETY2, HG_D,
     I                       EPSARG_D, EPSF_D,
     O                       Y4SW_D)
                Y4SW = Y4SW_D
              ELSE
C               Water is in contact already.
                Y4SW = Z1 - FDROP - Z4B
              ENDIF
            ELSE
C             Set level to the zero flow level, that is, to the level
C             at section 1. 
              Y4SW = Y1 + Z1B - Z4B
            ENDIF
C            WRITE(STDOUT,*) ' Y4SW=',Y4SW,' FDROP=',FDROP
           

            H4F = Z1 - FDROP - HDATUM
            Y4F = H4F + HDATUM - Z4B
            CALL XLKT22
     I                 (DEPTAB,
     M                  Y4F,
     O                  A4, TT, DT, JT, KT, DKT, BETA, DBETA, ALPHA4,
     O                  DALPHA, QCT)

            DE1TO4 = H1 + ALPHA1*(QFREE/A1)**2/TWOG -
     A         (H4F + ALPHA4*(QFREE/A4)**2/TWOG)
            IF(DE1TO4.LT.0.0) ELFLAG = 1
            
            CHAR8 = '   ---  '
            QMAT(IHU,NFRAC) = QFREE
            FDVEC(IHU) = FDROP
            CHAR6 = '  --- '
            CQ(1:10) = PUT10(QFREE)
            CQ(11:11) = ' '
            WRITE(STDOUT,78) 1.00, FDROP,
     A          CHAR8, H4F, FTYPE, CHAR6, CQ, DE1TO4
C            OLDH3 = Y3 + Z3B - HDATUM

C           Compute the submerged flows.  They may be of two types:
C           SO or SW.
            YVCOLD = CC*HG
            DO 400 J=NFRAC-1,2,-1
              DROP = FDROP*PFDVEC(J)
              ZW4 = Z1B + Y1 - DROP
              Y4 = ZW4 - Z4B
              H4 = ZW4 - HDATUM
              CALL XLKT22
     I                   (DEPTAB,
     M                    Y4,
     O                    A4, TT, DT, J4, KT, DKT, BETA4, DBETA, ALPHA4,
     O                    DALPHA, QCT)
              IF((Y4-Y4SW)/Y4SW.LE.1.E-5) THEN
C               Flow can only be submerged weir flow.  Look it up
C               in the table transfered from CULVERT.

                FTYPE = 'SW'
                CHAR6 = '  --- '
                CHAR8 = '   ---  '
                CALL TDLK13                                                
     I                     (STDOUT, GETQ, 13, 0, 0.0, ZW4, Z1,
     I                      HDATUM,
     O                      Q, DQED, DQEU, FREE)
                RFLAG = 0
              ELSE
C               The flow is SO (we think).  We might have problems 
C               because the CULVERT routines thinks it is but the
C               computations here may not.  In this case the tailwater
C               is at section 4
C              
                FTYPE = 'SO'
                WRITE(CHAR6,'(F6.3)') CC

                CALL FIND_SO_FLOW(
     I    STDOUT, DEPTAB, XSEC_ADRS(NXS), GRAV, H1, HG, BG, CC, 
     I    SNGL(CD_CULV),
     I    SNGL(SIN_THETA), SNGL(COS_THETA), SNGL(INVERT_Z(NXS)), 
     I    Z4B, Y4, A4, BETA4, J4,
     I    ALPHA1, A1,
     I    NXS, STATION, EPS_D,
     M    YVCOLD,
     O    Q, Y3, XE_R, RFLAG)
                IF(RFLAG.EQ.-11) THEN
                  WRITE(STDOUT,*) ' SO FLOW is not SO!: XE_R=',XE_R
                  Q = QFREE
                  YVCOLD = HG*CC
                ELSEIF(RFLAG.EQ.-10) THEN
                  WRITE(STDOUT,*) ' SO FLOW is not SO!'
                  Q = QFREE
                  YVCOLD = HG*CC
                ENDIF
                H3 = Y3 + Z3B - HDATUM
                WRITE(CHAR8,'(F8.3)') H3

              ENDIF
              QMAT(IHU,J) = Q

C             Prevent problems wiht local power, P, if one
C             of the flows should be zero.
              IF(Q.LE.0.0) THEN
                QTT = 0.1*QMAT(IHU,J+1)
              ELSE
                QTT = Q
              ENDIF
              P = LOG(QTT/QMAT(IHU,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
              DE1TO4 = H1 + ALPHA1*(Q/A1)**2/TWOG -
     A           (H4 + ALPHA4*(Q/A4)**2/TWOG)
              IF(DE1TO4.LT.0.0) ELFLAG = 1
 
              CQ(1:10) = PUT10(Q)
              IF(RFLAG.EQ.-12) THEN
                CQ(11:11) = '*'
              ELSE
                CQ(11:11) = ' '
              ENDIF
              WRITE(STDOUT,75) PFDVEC(J), DROP,
     A          CHAR8, H4, FTYPE, CHAR6, CQ, P, DE1TO4
 400        CONTINUE
          ELSEIF(H1.LT.H1FOLL.AND.FCFLAG.EQ.1) THEN
            IF(FCFLAG.EQ.1) THEN
C             FC flow 
              FTYPE = 'FC'

C             Define the contraction coefficient to use. 
C             If FW_CC_SO > 0 then use linear interpolation
C             between it and CC over  H1FWUL TO H1FOLL
  

              IF(FW_CC_SO.GT.0) THEN
                CC_D = FW_CC_SO + (H1_D - H1FWUL_D)
     A                 *(DBLE(CC) - FW_CC_SO)/(H1FOLL_D - H1FWUL_D) 
              ELSE
                CC_D = DBLE(CC)
              ENDIF

              CALL FIND_FCQ(
     I                    STDOUT, CD_CULV, H1_D,
     O                    QD, YVC, RFLAG)
              
              QFREE = QD

C             Compute the tailwater for the current flow.  The
C             exit contion can be either of 3 cases: FULL=0,
C             FULL=1, or FULL=2.  When FULL=1 we have a weighted
C             average of full-barrel flow and part-full-barrel flow. 
              IF(FULL_RETURN.EQ.0) THEN
C               We have critical depth at the exit and a free
C               surface.
                Y3 = YC_RETURN
                CALL LKTA
     I                  (XSEC_ADRS(NXS),
     M                   Y3,
     O                   AFLUX)
              ELSE
C               We have  a full pipe with piezometric surface less
C               than vertical diameter from invert.  Ignore the 
C               the small region of averaging the free surface
C               and the full flow.
                Y3 = Y_RETURN
                AFLUX = AFULL
              ENDIF
               
C             Next compute the water level at section 4 from the
C             water level at section 3.
              CALL FIND_Y4_FROM_Y43(
     I               STDOUT, GRAV, QFREE, DEPTAB, 
     I               Z3B, Z4B, Y3, AFLUX, SNGL(COS_THETA),
     O               Y4F, EFLAG )

              ZW4F = Y4F + Z4B
              H4F = ZW4F - HDATUM
              CALL XLKT22
     I                   (DEPTAB,
     M                    Y4F,
     O                    A4, TT, DT, J4, KT, DKT, BETA4, DBETA, 
     O                    ALPHA4, DALPHA, QCT)
              DE1TO4 = H1 + ALPHA1*(QFREE/A1)**2/TWOG -
     A           (H4F + ALPHA4*(QFREE/A4)**2/TWOG)
              IF(DE1TO4.LT.0.0) ELFLAG = 1
              QMAT(IHU,NFRAC) = QFREE
              FDROP = Y1 + Z1B - ZW4F
              FDVEC(IHU) = FDROP
              CQ(1:10) = PUT10(QFREE)
              CQ(11:11) = ' '
              WRITE(CHAR6,'(F6.3)') CC_D
              WRITE(STDOUT,76) 1.00, FDROP,
     A            Y3 + Z3B - HDATUM, H4F, FTYPE, CHAR6, CQ, DE1TO4
              FTYPE = 'SO'
              YVCOLD = CC*HG
              DO 450 J=NFRAC-1,2,-1
                DROP = FDROP*PFDVEC(J)
                ZW4 = Z1B + Y1 - DROP
                Y4 = ZW4 - Z4B
                CALL XLKT22
     I                     (DEPTAB,
     M                      Y4,
     O                      A4, TT, DT, J4, KT, DKT, BETA4, DBETA, 
     O                      ALPHA4, DALPHA, QCT)
                H4 = ZW4 - HDATUM
                CALL FIND_SO_FLOW(
     I    STDOUT, DEPTAB, XSEC_ADRS(NXS), GRAV, H1, HG, BG, CC, 
     I    SNGL(CD_CULV),
     I    SNGL(SIN_THETA), SNGL(COS_THETA), SNGL(INVERT_Z(NXS)), 
     I    Z4B, Y4, A4, BETA4, J4,
     I    ALPHA1, A1,
     I    NXS, STATION, EPS_D,
     M    YVCOLD,
     O    Q, Y3, XE_R, RFLAG)
                IF(RFLAG.EQ.-10) THEN
                  WRITE(STDOUT,*) ' SO FLOW is not SO!: XE_R=',XE_R
                  Q = QFREE
                  YVCOLD = HG*CC
                ENDIF
                QMAT(IHU,J) = Q
                H3 = Y3 + Z3B - HDATUM
                WRITE(CHAR8,'(F8.3)') H3
                P = LOG(Q/QMAT(IHU,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
                DE1TO4 = H1 + ALPHA1*(Q/A1)**2/TWOG -
     A             (H4 + ALPHA4*(Q/A4)**2/TWOG)
                IF(DE1TO4.LT.0.0) ELFLAG = 1
 
                CQ(1:10) = PUT10(Q)
                IF(RFLAG.EQ.-12) THEN
                  CQ(11:11) = '*'
                ELSE
                  CQ(11:11) = ' '
                ENDIF
                WRITE(STDOUT,75) PFDVEC(J), DROP,
     A            CHAR8, H4, FTYPE, CHAR6, CQ, P, DE1TO4
 450        CONTINUE

            ELSE
C             Should not get here. 
              WRITE(STDOUT,*) ' no FC flow wrong branch.'
              STOP 'Abnormal stop. Errors found.'
            ENDIF

          ELSE
C           Free flow is FO here
 
            FTYPE = 'FO'
            RG = HG/H1

             
            WRITE(CHAR6,'(F6.3)') CC

            CALL UFG_FNDFOQ
     I                   (STDOUT, APPTAB, H1_D, DBLE(HDATUM), 
     I                    DBLE(HG), DBLE(BG), DBLE(Z1B), DBLE(CC), 
     I                    CD_CULV, TWOG_D,
     O                    QD)

            QFREE = SNGL(QD)
            Y2 = CC*HG

C           Find the water depth at the culvert exit and also check
C           to make sure the vena contracta is not drowned. 

            YCRIT = 0.5D0*VERT_D
            YNORM = YCRIT
            CALL GET_YCYNYM
     I                   (STDOUT, STATION(NXS), EPSARG_D, EPSF_D, 
     I                    EPSABS_D,
     M                    YCRIT, YNORM,
     O                    YMAX, RFLAG)

            YVC = HG_D*DBLE(CC)
            IF(YNORM.LT.0.D0.OR.YCRIT.LT.YNORM) THEN
C             The slope is mild or flat.  To check for possible 
C             drowning of the vena contract, compute a profile
C             upstream using the same methods as in FIND_FOLL and
C             FIND_FCQ.  

              Y3 = YCRIT

              CALL PROFILE_UPS(
     I                         YCRIT, YNORM, YVC, DBLE(CC),
     O                         RFLAG)

              IF(RFLAG.EQ.0) THEN
C                WRITE(STDOUT,*) ' Vena contracta is drowned.'
              ELSE
C                WRITE(STDOUT,*) ' Vena contracta is free.'
              ENDIF
            ELSE
C             The slope is steep.  Compute a free surface
C             profile from the vena contracta towards the exit. 
              YS = YVC
              XE = STATION(NXS)
              XS = STATION(1) + HG_D
              CALL SFWSP(
     I               STDOUT, XS, XE, YS, EPS_D,
     I               E, RHSE, FIND_YCYNYM, NMAX,
     O               NS, NE, XVEC, YVEC, YCVEC, YNVEC, EFLAG)
              IF(ABS((XVEC(NE)- XE)
     A                /XE).LT.0.001) THEN
C                WRITE(STDOUT,*) ' Super critical flow to end. Y=',
C     A                     SNGL(YVEC(NE))
                Y3 = YVEC(NE)
              ELSE
                WRITE(STDOUT,*) ' Jump expected when none should appear'
                WRITE(STDOUT,*) ' YCRIT=',YCRIT,' YNORM=',YNORM
                WRITE(STDOUT,*) ' XVEC(NE)=',XVEC(NE),' YE=',
     A                           YVEC(NE)
                WRITE(STDOUT,*) ' QD=',QD
                STOP 'Abnormal stop. Errors found.'
              ENDIF
              
            ENDIF

C           Now find the level at section 4 that defines
C           the limit of free orifice flow, Y4F. Define the constant
C           values.  First, find the sequent depth at the vena
C           contracta

            CALL FIND_SEQUENT_DEPTH(
     I                 STDOUT, GRAV_D, HG_D, BG_D, DBLE(CC), QD,
     I                 EPSARG_D, EPSF_D,
     O                 YB_D)
 
C            WRITE(STDOUT,*) ' SEQUENT DEPTH=',YB_D
 
C           Now try to compute a water-surface profile
C           downstream from the sequent depth to the
C           culvert exit. 
            YS = YB_D
            XE = STATION(NXS)
            XS = STATION(1) + HG_D
            CALL SFWSP(
     I             STDOUT, XS, XE, YS, EPS_D,
     I             E, RHSE, FIND_YCYNYM, NMAX,
     O             NS, NE, XVEC, YVEC, YCVEC, YNVEC, EFLAG)
C            WRITE(STDOUT,*) ' Y3F=',YVEC(NE)
            Y3 = YVEC(NE)
            CALL LKTA
     I               (XSEC_ADRS(NXS),
     M                Y3,
     O                A3)

C           Next compute the water level at section 4 from the
C           water level at section 3.
            CALL FIND_Y4_FROM_Y43(
     I             STDOUT, GRAV, QFREE, DEPTAB,
     I             Z3B, Z4B, Y3, A3, SNGL(COS_THETA),
     O             Y4F, EFLAG )

C            WRITE(STDOUT,*) ' Y4F=',Y4F  



            ZW4F = Y4F + Z4B
            H4F = ZW4F - HDATUM
            CALL XLKT22
     I                 (DEPTAB,
     M                  Y4F,
     O                  A4, TT, DT, J4, KT, DKT, BETA4, DBETA, 
     O                  ALPHA4, DALPHA, QCT)
            
            DE1TO4 = H1 + ALPHA1*(QFREE/A1)**2/TWOG -
     A           (H4F + ALPHA4*(QFREE/A4)**2/TWOG)
            IF(DE1TO4.LT.0.0) ELFLAG = 1
            QMAT(IHU,NFRAC) = QFREE
            FDROP = Y1 + Z1B - ZW4F
            FDVEC(IHU) = FDROP
            CQ = PUT10(QFREE)
            CQ(11:11) = ' '
            WRITE(STDOUT,76) 1.00, FDROP,
     A          Y3 + Z3B - HDATUM, H4F, FTYPE, CHAR6, CQ, DE1TO4
C            OLDH3 = Y3 + Z3B - HDATUM
            FTYPE = 'SO'
C            AT = CD*CC*AG
            YVCOLD = CC*HG
            DO 500 J=NFRAC-1,2,-1
              DROP = FDROP*PFDVEC(J)
              ZW4 = Z1B + Y1 - DROP
              Y4 = ZW4 - Z4B
              CALL XLKT22
     I                   (DEPTAB,
     M                    Y4,
     O                    A4, TT, DT, J4, KT, DKT, BETA4, DBETA, ALPHA4,
     O                    DALPHA, QCT)
              H4 = ZW4 - HDATUM
              CALL FIND_SO_FLOW(
     I    STDOUT, DEPTAB, XSEC_ADRS(NXS), GRAV, H1, HG, BG, CC, 
     I    SNGL(CD_CULV),
     I    SNGL(SIN_THETA), SNGL(COS_THETA), SNGL(INVERT_Z(NXS)), 
     I    Z4B, Y4, A4, BETA4, J4,
     I    ALPHA1, A1,
     I    NXS, STATION, EPS_D,
     M    YVCOLD,
     O    Q, Y3, XE_R, RFLAG)
              IF(RFLAG.EQ.-10) THEN
                WRITE(STDOUT,*) ' SO FLOW is not SO!: XE_R=',XE_R
                Q = QFREE
                YVCOLD = HG*CC
              ENDIF
              QMAT(IHU,J) = Q
              H3 = Y3 + Z3B - HDATUM
              WRITE(CHAR8,'(F8.3)') H3
              P = LOG(Q/QMAT(IHU,J+1))/LOG(PFDVEC(J)/PFDVEC(J+1))
              DE1TO4 = H1 + ALPHA1*(Q/A1)**2/TWOG -
     A           (H4 + ALPHA4*(Q/A4)**2/TWOG)
              IF(DE1TO4.LT.0.0) ELFLAG = 1
 
              CQ = PUT10(Q)
              IF(RFLAG.EQ.-12) THEN
                CQ(11:11) = '*'
              ELSE
                CQ(11:11) = ' '
              ENDIF
              WRITE(STDOUT,75) PFDVEC(J), DROP,
     A          CHAR8, H4, FTYPE, CHAR6, CQ, P, DE1TO4
 500        CONTINUE
 
          ENDIF

C         Compute approximate maximum error and report
          DO 600 J=3,NFRAC-2,2
            QHAT = 0.5*(QMAT(IHU,J-1) + QMAT(IHU,J+1))
            RERR = ABS(QHAT - QMAT(IHU,J))/QMAT(IHU,J)
            RMS = RMS + RERR*RERR
            NRMS = NRMS + 1
            IF(RERR.GT.BIGERR) THEN
              BIGERR = RERR
              HERR = H1
              IPFD = J
              IHGERR = IHG
            ENDIF
 600      CONTINUE
 
CC         Eliminate the checking values from QMAT
C          JBASE = 3
C          DO 700 J=4,NFRAC-1,2
C            QMAT(IHU,JBASE) = QMAT(IHU,J)
C            JBASE = JBASE + 1
C 700      CONTINUE
C          QMAT(IHU,JBASE) = QMAT(IHU,NFRAC)

900    CONTINUE
C        PFDTMP(1) = PFDVEC(1)
C        PFDTMP(2) = PFDVEC(2)
C        JBASE = 3
C        DO 910 J=4,NFRAC-1,2
C          PFDTMP(JBASE) = PFDVEC(J)
C          JBASE = JBASE + 1
C 910    CONTINUE
C        PFDTMP(JBASE) = PFDVEC(NFRAC)
 
        zrhufd = 0.0
        CALL TWDOUT
     I             (STDOUT, STDTAB, TAB2D(IHG), LABEL, NHU, NFRAC,
     I              HUVEC, FDVEC, PFDVEC, QMAT, HDATUM,
     I              13,' UFGCULV', zrhufd,
     i              zone, hgrid, vdatum, unitsys, basis,
     i              easting, northing,
     O              EFLAG)
 
        IF(ELFLAG.GT.0) THEN
          WRITE(STDOUT,77)
        ENDIF
 1000 CONTINUE



      WRITE(STDOUT,86) BIGERR, HGVEC(IHGERR), HERR,
     A                 PFDVEC(IPFD)
 
      RMS = SQRT(RMS/FLOAT(NRMS))
      WRITE(STDOUT,88) RMS, NRMS
 
C     Output the table of type 15 to the table file
      TABID = GET_TABID(TAB)
      WRITE(STDTAB,95) TABID(1:LENSTR(TABID)), LAB2
      DO 1005 I=1,NHG
        WRITE(STDTAB,96) HGVEC(I), GET_TABID(TAB2D(I)), 
     A              (HSTUFF(I,J), J=1,3)
 1005 CONTINUE
      WRITE(STDTAB,96) -1.0
      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) LINE
        STOP 'Abnormal stop.  Errors found.'
      END
C
C
C
      SUBROUTINE SET_UFGC_TABIDS_DEFAULTS()

C     Set the default values in the vectors used to process a
C     deparature reach table spec.
      IMPLICIT NONE
      INCLUDE 'ufgcitm.cmn'

C***********************************************************************
C     Default for: TABID 
      UFGCITMCTAB(  1) = '    '           
C     Default for: TABLE - note # is ignored in the standard scanner 
      UFGCITMCTAB(  2) = '    '           
C     Default for: GETQ
      UFGCITMCTAB(  3) = '    '
C     Default for: GETY2
      UFGCITMCTAB(  4) = '    '
      RETURN
      END
C
C
C
      SUBROUTINE  SET_UFGC_TABIDS(
     O                       TAB, GETQ, GETY2,
     M                       EFLAG)

C     Set items UFGCULV
C     All values not set explicitly by user are at their default value.

      IMPLICIT NONE

      INTEGER TAB, GETQ, GETY2, EFLAG
    

C     Local
      CHARACTER*16 KEY, KEY1, KEY2

      INCLUDE 'ufgcitm.cmn'
      INCLUDE 'stdun.cmn'
C***********************************************************************
      KEY1 = UFGCITMCTAB(  1)
      KEY2 = UFGCITMCTAB(  2)
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
      KEY = UFGCITMCTAB(  3)
      IF(KEY.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY,
     M                                EFLAG,
     O                                GETQ)
      ELSE
        GETQ = 0
      ENDIF
      KEY = UFGCITMCTAB(  4)
      IF(KEY.NE.' ') THEN
        CALL GET_INTERNAL_TAB_NUMBER
     I                               (STD6, KEY,
     M                                EFLAG,
     O                                GETY2)
      ELSE
        GETY2 = 0
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE GET_UFGC_TABIDS(STDIN, STDOUT,
     M                   EFLAG)

C     Get the table ids for the a culvert with an underflow gate

      IMPLICIT NONE

      INCLUDE 'arsize.prm'

      INTEGER STDIN, STDOUT, EFLAG

      INCLUDE 'ufgcitm.cmn'

C     Local

C     + + + LOCAL PARAMETERS + + +
      INTEGER  INTVAL, REAVAL, 
     A         CHRVAL, DPRVAL, EXACT, LOWER,
     B         NUMERIC, CHAR, N_SYMBOL
      PARAMETER(N_SYMBOL=4, INTVAL=1, REAVAL=2,
     A          DPRVAL=3, CHRVAL=4, 
     B          EXACT=0,LOWER=1, NUMERIC=0, CHAR=1)

      INTEGER MAX_LINE
     A        

      EXTERNAL GET_NAMED_ITEMS, SET_UFGC_TABIDS_DEFAULTS

C     + + + SAVED VALUES + + +
      INTEGER GROUP(N_SYMBOL), RESPONSE_TYPE(N_SYMBOL),
     A        CONVERT_RULE(N_SYMBOL), GROUP_INDEX(N_SYMBOL)
      CHARACTER SYMBOL_TABLE(N_SYMBOL)*16

      SAVE  SYMBOL_TABLE, GROUP, RESPONSE_TYPE, CONVERT_RULE,
     A      GROUP_INDEX

      DATA  SYMBOL_TABLE /
     *'TABID','TABLE','GETQ','GETY2'/
                                                                      
      DATA GROUP  /
     *CHAR, CHAR, CHAR, CHAR/          

      DATA GROUP_INDEX /
     *1, 2, 3, 4/
                                       
      DATA RESPONSE_TYPE  /
     *CHRVAL,CHRVAL,CHRVAL, CHRVAL/
    
      DATA CONVERT_RULE /
     *LOWER,LOWER,LOWER,LOWER/
   
C***********************************************************************
C     Set Defaults
 
      CALL SET_UFGC_TABIDS_DEFAULTS()

      MAX_LINE = 1
      CALL GET_NAMED_ITEMS(
     I                     STDIN, STDOUT,  MAX_LINE, N_SYMBOL,
     I  GROUP, RESPONSE_TYPE, CONVERT_RULE, GROUP_INDEX, SYMBOL_TABLE,
     I  MAXR_UFGCITM, MAXDP_UFGCITM, MAXC_UFGCITM, 'UFGCULV items',
     O  UFGCITMITAB, UFGCITMFTAB, UFGCITMDTAB, UFGCITMCTAB,
     O  EFLAG)
      
      RETURN

      END
C
C
C
      SUBROUTINE READ_BARREL_ITEMS(
     I                            STDOUT, LINE, NITEM, ITEM_START,
     I                            ITEM_END,
     M                            EFLAG,
     O                            NODE, NAME, TABN, STAT, 
     O                            ELEV)

C     Get the items of data from a barrel input line in UFGCULV

      IMPLICIT NONE
      INTEGER STDOUT, NITEM, ITEM_START(NITEM), ITEM_END(NITEM),
     A        NODE,  EFLAG, TABN

      REAL*8 STAT, ELEV

      CHARACTER LINE*80, NAME*8

C     Local

      INTEGER IE, IS, ITAB, LKEY, N
      CHARACTER TPC*20, KEY*16

C     Called program units
      EXTERNAL STRIP_L_BLANKS, GET_INTERNAL_TAB_NUMBER
C     ***********************FORMATS************************************
50    FORMAT(/,' *ERR:XXX* Only ',I3,' items given in ',
     A   'ufgculv-barrel description line.  Need at least five items.')
52    FORMAT(/,' *ERR:XXX* No interpolation allowed for a culvert',
     A          ' barrel for UFGCULV.')
54    FORMAT(/,' *ERR:XXX* Cross-section id missing for a culvert',
     A         ' barrel for UFGCULV.')
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

      IF(NODE.LT.0) RETURN

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
C     Convert from the table id to an internal number.
      IF(KEY.NE.' ') THEN
        IF(KEY(1:1).NE.'-') THEN
C         We have an id given.
          CALL GET_INTERNAL_TAB_NUMBER
     I                                (STDOUT, KEY,
     M                                 EFLAG,
     O                                 TABN)
          
        ELSE
          WRITE(STDOUT,52)
          STOP 'Abnormal stop.  Errors found.' 
        ENDIF
      ELSE
        WRITE(STDOUT,54)
        STOP 'Abnormal stop.  Errors found.' 
      ENDIF
      
C     Process the station 
      N = 4
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F12.0)',ERR=991) STAT

C     Process the invert elevation
      N = 5
      IS = ITEM_START(N)
      IE = ITEM_END(N)
      TPC = LINE(IS:IE)
      CALL STRIP_L_BLANKS(
     M                    TPC)
      READ(TPC,'(F12.0)',ERR=991) ELEV

      RETURN
 991  CONTINUE
        WRITE(STDOUT,*) ' *ERR:500* Conversion error in line:'
        WRITE(STDOUT,*) TPC
        STOP 'Abnormal stop.  Errors found.' 

      END
