c     28 Feb 2003:  Start to develop ability to fit existing tables
c     with functions that will provide at least a continuous first derivative
c     and often a continuous second derivative.  Do this in the hope that 
c     the Newton numerical solution will be more predictable and controllable. 




      

C
C
C
      subroutine   twodfit
     I                   (STDOUT, TABLE, NHU, NFRAC,
     I                    HUVEC, FDROP, PFDVEC, Q, HDATUM,
     I                    TYPE, SOURCE, zrhufd, verbose,
     m                    ftp,
     O                    EFLAG, ftpup)
 
C     + + + PURPOSE + + +
C     Experiment with fitting 2-d tables of type 13.
 
      IMPLICIT NONE

C     + + + PARAMETERS + + +
      INCLUDE 'arsize.prm'
 
C     + + + DUMMY ARGUMENTS + + +
      INTEGER EFLAG, FTP, ftpup,  NFRAC, NHU, STDOUT,  TABLE, 
     a        TYPE, verbose
      REAL FDROP(PMXNHU), HDATUM, HUVEC(PMXNHU), PFDVEC(PMXFRC),
     A     Q(PMXNHU,PMXFRC), zrhufd
      CHARACTER  SOURCE*8
 
C     + + +DUMMY ARGUMENT DEFINITIONS + + +
C     STDOUT - Fortran unit number for user output and messages
C     TABLE  - Table number
C     NHU    - Number of upstream heads
C     NFRAC  - Number of fractions for defining partial free drop
C     HUVEC  - Vector of upstream heads
C     FDROP  - Free drop values
C     PFDVEC - Partial free drop vector
C     Q      - Flowrate
C     HDATUM - Datum for measuring head
C     TYPE   - Table type.
C     SOURCE - source for the table.
C     EFLAG  - Error flag: EFLAG=0-no errors; else one or more errors
 
   

C     + + + SAVED VALUES + + +
      INTEGER IOFF
      SAVE IOFF
 
C     + + + LOCAL VARIABLES + + +
      INTEGER I, J, K, L, KNT, LIM, ihu, ipfd, free
      real ed, eu, qt, dqted, dqteu, dfl, dfm, dfr, dqpfd,
     a     qtl, qtm, qtr, diff, spfd
      CHARACTER DUMMY*7, LINE(10)*7, TABID*16
 
      integer n
      real*8 
     b       huall(pmxnhu+1),
     e       freedrop(pmxnhu+1),
     f       qmat(pmxnhu+1,pmxfrc), pfd(pmxfrc)
C***********************************************************************
      IF(TYPE.EQ.13) THEN

c       Need to add zero upstream head to the mix.  Flow is zero 
c       but some partial derivatives are not!
        huall(1) = 0.d0
        freedrop(1) = dble(zrhufd)
        n = nhu + 1
        do i=1,nhu
          huall(i+1) = huvec(i)
          freedrop(i+1) = fdrop(i)
        end do
      
        do ipfd=1,nfrac
          pfd(ipfd) = dble(pfdvec(ipfd))
        end do

        do ihu=1,n
          if(ihu.eq.1) then
c           special case.
            do ipfd=1,nfrac
              qmat(ihu,ipfd) = 0.d0
            end do
          else
            do ipfd=1,nfrac
              qmat(ihu,ipfd) = dble(q(ihu-1,ipfd))
            end do
          endif
        end do

        call   twodtabfit
     I                   (stdout, table, n, nfrac, pmxnhu+1, hdatum,
     i                    huall, freedrop, pfd, qmat,
     I                    source, verbose,
     m                    ftp,
     O                    EFLAG, ftpup)


      endif
      RETURN
      END


     
c
c
c
      subroutine locate_extreme_point(stdout, yl, yr, ml, mr, h,
     o                            nex, t_at_ex)

c     Given the function values, yl, yr, and the derivative values, 
c     ml, mr, at each end of an interval, h, compute the location of
c     extreme point for the cubic that matches the four data points. 
c     Return the extreme point in terms of the relative distance from
c     the left end point.  

      implicit none

      integer nex, stdout

      real*8 yl, yr, ml, mr, h, t_at_ex(2)

c     Local

      real*8 a, b, c, t1, t2

c***********************************************************************

      c = h*ml

      b = 2.*( (-2.*ml - mr)*h + 3.*(yr - yl))

      a = 3.*(h*(ml + mr) - 2.*(yr - yl))

      if(b*b - 4.*a*c.lt.0.d0) then
c        write(stdout,*) ' Imaginary solution in locate_extreme_point'
        nex = 0
        
        return
      endif

      if(a.eq.0.d0) then
        if(b.ne.0.d0) then
          t1 = -c/b
          t2 = -1.0
        else
c          write(stdout,*) 
c     a 'both a and b are zero in locate_extreme_point'
c          write(stdout,50) yl, yr, ml, mr, h
c50    format(' ql=',1pe12.5,' qr=',1pe12.5,' ml=',1pe12.5,
c     a       ' mr=',1pe12.5,' h=',f10.4)
c          stop 'Abnormal stop. Bug found.'
          nex = 0
          return
        endif
      else
c       a ne to zero here!
        if(c.eq.0.d0) then
          t1 = -b/a
          t2 = 0.d0
        else
c         neither a nor c are zero here. 
          if(b.lt.0.d0) then
            t1 = (-b + sqrt(b*b - 4.*a*c))/(2.*a)
            t2 =  2.*c/(-b + sqrt(b*b - 4.*a*c))
          else
            t1 = 2.*c/(-b - sqrt(b*b -4.*a*c))
            t2 = (-b - sqrt(b*b - 4.*a*c))/(2.*a)
          endif
        endif
      endif
      if(t1.ge.0.d0.and.t1.le.1.0d0) then
        nex = 1
        t_at_ex(nex) = t1
      else
        nex = 0
      endif
      if(t2.ge.0.d0.and.t2.le.1.0d0) then
        nex = nex + 1
        t_at_ex(nex) = t2
      endif

      
      return
      end
