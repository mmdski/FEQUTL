c     Program units related to tracking of files.  Will 
c     expand over time as we extend coverage to a greater
c     depth in the files involved in a model.

c
c 
c 
      subroutine update_message(line,
     m                              m, message)

c     Add the non-blank contents of line to the end of message .
c     On entry m points to the last character added to message

      implicit none

      integer m
      character*(*) line
      character*(*) message

c     Local
      integer n, i
      character*256 string
c***********************************************************************
      string = adjustl(line)
      
      n = len_trim(string)
      do i=1,n
        if(string(i:i) /= ' ') then
          m = m + 1
          message(m:m) = string(i:i)
        endif
      end do
      return
      end


c
c 
c 
      subroutine copy_master_input_file(stdin, stdscr)

c     Copy the master-input to a scratch file so that we can process
c     it for computing md5 digests of various blocks of input after we
c     have processed the header information.    

      implicit none

      integer stdin, stdscr

c     Local
      integer it

      character*256 line
c***********************************************************************
c     open the scratch file
      open(unit=stdscr, file='scratch', status='unknown')  

      do
        read(stdin,'(a)', end=999) line
        it = len_trim(line)
        write(stdscr,'(a)') line(1:it)
        if(line(1:6) == 'FINISH') then
          exit
        endif
      end do 
c     Do not close the scratch file because if we do it will be deleted.
      rewind stdscr
      rewind stdin
      return  

999   continue
c     We should reach here only if the user has failed to properly terminate the 
c     master-input file for FEQUTL.  Supply a FINISH statement for the scratch 
c     file to prevent problems later.  The end of file on the master-input file 
c     is detected elsewhere during processing.
      line = 'FINISH'
      it = 6
      write(stdscr,'(a)') line(1:it)
      rewind stdscr
      return

      end
c
c 
c 
      subroutine find_md5_for_fequtl_input(stdscr, stdout, stdsys, 
     i                           stdin, fname, ncmd, cmdtab, cmdval,
     i                           conf_file, ghome )

c     Compute the standard md5 digests for various parts of the master-input 
c     file for fequtl and place the results in a special file.  The master-input
c     file is open at this point but nothing has been read from it.

      implicit none

      integer stdscr, stdout, stdsys, stdin, ncmd

      character*(*) fname, conf_file, ghome
      integer cmdval(ncmd)
      character*8 cmdtab(ncmd)

c     Called program units
      integer get_unit
      external get_unit

c     Local

      integer MAXM
      parameter (MAXM=256000)

      integer i, ip, it, j, is, ie,  m, stdmd5, stdtmp,
     a            stdmd5_header, stdmd5_mif

      character line*192, message*(MAXM), cmd_in_process*8, 
     a           cmd*8, md5_digest*32, md5_file_name*256,
     b          tabid*16, class*8, md5_header_file_name*256,
     c          name*9, command(44)*8, standard_name*8,
     d          cwd*256, mif_with_path*256, conf_file_with_path*256

c     *****************************formats******************************      
50    format(a16,1x,a8,1x,i5,i5,1x,a8,1x,a32, 1x,a)
c***********************************************************************
c     Setup the standard command names.  Needed because a user can change
c     the name of a command!
      command = 'null'
      command( 1) = 'FEQX    '
      command( 2) = 'FLOODWAY'    
      command( 3) = 'BRIDGE  '    
      command( 4) = 'CULVERT '    
      command( 5) = 'FINISH  '    
      command( 8) = 'FEQXLST '    
      command(10) = 'SEWER   '    
      command(11) = 'MULPIPES'    
      command(12) = 'FTABIN  '    
      command(13) = 'EMBANKQ '    
      command(15) = 'CRITQ   '    
      command(16) = 'GRITTER '    
      command(18) = 'MULCON  '    
      command(19) = 'CHANRAT '    
      command(20) = 'EXPCON  '    
      command(21) = 'HEC2X   '    
      command(22) = 'QCLIMIT '    
      command(23) = 'XSINTERP'    
      command(25) = 'FEQXEXT '    
      command(26) = 'CHANNEL '    
      command(27) = 'WSPROX  '    
      command(28) = 'WSPROQZ '    
      command(29) = 'WSPROT14'    
      command(30) = 'UFGATE  '    
      command(31) = 'RISERCLV'    
      command(32) = 'ORIFICE '    
      command(33) = 'AXIALPMP'    
      command(34) = 'PUMPLOSS'    
      command(35) = 'SETSLOT '    
      command(36) = 'CLRSLOT '    
      command(39) = 'MKEMBANK'    
      command(40) = 'MKWSPRO '    
      command(41) = 'wsprot14'    
      command(42) = 'LPRFIT  '    
      command(44) = 'SETSLOTE'    

c     Create the file name for the md5 digest and related info.  It is
c     formed from the master-input file name by appending '.md5' to 
c     that name.  The file will be placed in the current working directory
c     from which fequtl was invoked.  The standard usage for fequtl is to execute it in the
c     directory that contains the files on its command line.  
c     Find the current working directory.
      call pwd(
     o         cwd)
      write(stdout,*) ' Current working directory=',cwd

c     Add the path to the master-input file (mif) name. 
      it = len_trim(cwd)
      
      mif_with_path = cwd(1:it)//'/'//fname
      call os_file_style(
     m                   mif_with_path)
      it = len_trim(mif_with_path)
      
      md5_file_name = mif_with_path(1:it)//'.md5'
      stdmd5_mif = get_unit(0)
      open(unit=stdmd5_mif, file=md5_file_name,status='unknown')


      it = len_trim(ghome)
      mif_with_path = mif_with_path(it+1:)

      write(stdout,*) 
     a ' MIF name relative to global home name=',
     b   mif_with_path
      
c     General outline of the approach:  For each command block of the 
c     input:  
c       record first line number and last line number
c       compute the message as a concatenation of all characters in 
c       block with all spaces removed.  We also delete blank lines,
c       and comment lines.  However, trailing comments on a line are 
c       retained.  However, these lines are retained when establishing
c       the starting and ending line numbers.  
c     
c       When the block has a TABID it will be used to reference the 
c       md5 digest and related info.  When it does not have a TABID,
c       it may stand along, like FTABIN, or it may be a modifier of
c       a command to come later.  This applies to SETSLOTx, and CHANNEL,
c       for example.  We need to remember these later two and then 
c       add them to the description for the block or blocks they affect
c       and that do have a TABID.  There are some cases in which a 
c       block has a TABID but yet is also a modifier to a subsequent
c       command.  An example of this would be the FEQxxx family of commands
c       that compute cross section tables.  These are often used to compute
c       the cross-section function tables required for a subsequent 
c       command.  There may be no table output from such a command and
c       then it is a pure modifier.  But it is also possible that the 
c       FEQxxx command is both a modifier and also a creater of a 
c       function table that is placed in the *.tab file.  Thus the 
c       prescence of a SAVE option creates a modifer state and the 
c       prescence of NOOUT makes it a pure modifier.  As is already 
c       clear, we will have to make some rather fine distinctions to 
c       make sense of the tracking system. 
c 
c       Many details yet to be defined: 13 Dec 2007
c       14 Dec 2007:  After a few test runs it is apparent that
c       we must implement command specific process here.  Each command
c       has its nuances.  Thus we will start with implementing the 
c       simpler commands.  This also requires that we allow for a 
c       there being no tracking for some function tables.  
c 
c     Process the header block. If stdin = stdsys, then the header block
c     is in the file pointed to by stdscr, otherwise, it is in the 
c     file pointed to by stdsys and has ALREADY been read to the end. 
c     Thus we must rewind the file before we start to process it!  We are 
c     interested in the following lines for the reasons given:

c     UNITS=  ENGLISH  NOMINAL      45.0       0.0
c      The above line affects the value of g used.  Thus affects
c      anything that includes the value of the gravitational acceleration. 
c 
c     DZLIM= 
c       The above line affects cross-section function tables
c 
c     NRZERO=
c       The above line affects cross-section function tables
c 
c     EPSARG
c       Affects any commands that use an iteration to find a solution.
c 
c     EPSF=      EPSABS=
c       Affects any commands that us an iteration to find a solution.
c 
c     EXTEND=
c       Affects cross-section function tables
c 
c     MINQ=
c      Affects flow tables from EMBANKQ and CHANRAT
c 
c     We will compute the md5sum for each of these lines, if present,
c     and store each using the first name on the line.  These will all
c     be modifiers,  that is, they are values whose change affect all 
c     tables of a given class.

      if(stdsys == stdin) then
c       The header information is in the file pointed to by stdscr
c       and has not been read yet!  The header block is in the master-input
c       file to fequtl
        stdtmp = stdscr
        stdmd5 = stdmd5_mif
        conf_file_with_path = mif_with_path
      else
c       The header information is in the file pointed to by stdsys and
c       has been read already!   THe header block is in a distinct configuration
c       file. 
        rewind(stdsys)
        stdtmp = stdsys

c       Get the path name for the conf file.  The file was given 
c       with a full path name (perhaps missing a drive letter under
c       MSW) or it was in the current working directory.  We first
c       try to determine if the full path name was given, in which
c       case we have what we want.  If that fails, then the file 
c       was in the current working directory, so we prepend the path
c       we found for the MIF. 
        conf_file_with_path = conf_file
        call os_file_style(
     m                     conf_file_with_path)
        it = index(conf_file_with_path, ':') !If drive letter is presenst so is a colon.
        if(it > 0 .or. 
     a     conf_file_with_path(1:1) == '~' .or. 
     b     conf_file_with_path(1:1) == '/' .or. 
     c     conf_file_with_path(1:1) == '\') then
c         A full path name was given when the file was referenced
          continue
        else
c         Prepend the cwd found above
          it = len_trim(cwd)
          conf_file_with_path = cwd(1:it)//'/'//conf_file
          call os_file_style(
     m                       conf_file_with_path)
        endif

c       create the md5 name for the header
        it = len_trim(conf_file_with_path)
        md5_header_file_name = conf_file_with_path(1:it)//'.md5'
        stdmd5_header = get_unit(0)
        open(unit=stdmd5_header, file=md5_header_file_name,
     a           status='unknown')
        stdmd5 = stdmd5_header
c       Get the file name for the standard header relative to the global
c       home name.
        it = len_trim(ghome)
        conf_file_with_path = conf_file_with_path(it+1:)
      endif
      write(stdout,*) 
     a 'conf file name relative to global home name=',
     b  conf_file_with_path

      i = 0   !Clear the line counter
      do 
        read(stdtmp,'(a)')  line
        i = i + 1

        if(line(1:6) == 'UNITS=') then
          name = 'UNITS='
        elseif(line(1:6) == 'DZLIM=') then
          name = 'DZLIM='
        elseif(line(1:7) == 'NRZERO=') then
          name = 'NRZERO='
        elseif(line(1:7) == 'EPSARG=') then
          name = 'EPSARG='
        elseif(line(1:5) == 'EPSF=') then
          name = 'EPSF='
        elseif(line(1:7) == 'EXTEND=') then
          name = 'EXTEND='
        elseif(line(1:5) == 'MINQ=') then
          name = 'MINQ='
        elseif(line == ' ') then
          name = 'done'
        else
          name = ' '
        endif

        if(name == 'done') then
c         We have reached the end of the header block. Tidy up afterward!
          if( stdtmp == stdsys) then
c           The header block was in its own file. Close stdmd5, reset
c           to the master-input file value, clear the line counter.
            close(stdmd5)
            stdmd5 = stdmd5_mif
            i = 0
          endif
          exit
        endif

        if(name /= ' ') then
          m = 0
          message = ''
          is = i
          ie = i
          tabid = name
          class = 'modifier'
          cmd_in_process = ' '
          call update_message(line,
     m                        m, message)
          call md5(message(1:m),
     o              md5_digest)
          it = len_trim(conf_file_with_path)
          write(stdmd5,50) tabid, class, is, ie, cmd_in_process, 
     a             md5_digest, conf_file_with_path(1:it) 
        endif
      end do

      
      cmd_in_process = ' ' !Clear the command we are processing
      tabid = ' '
      class = ' '
      do
        read(stdscr,'(a)') line
       
        i = i + 1
        write(stdout,*) 'From scratch i=',i,' line=',line


c       Is this line a comment, or is it blank?  Note: old style comments are not
c       detected and become part of the command block that precedes them--DO NOT use
c       old style comments!

        if(line(1:1) == '*' .or. line(1:1) == '+' .or. line(1:1) == ';' 
     a      .or. line == ' ') then
          write(stdout,*) ' i=',i,' Skipping line:',line
          cycle
        else
c         Look for a command.  It should be the only value on the line.
          cmd = line(1:8)
          CALL BINSER
     I           (cmd, NCMD, CMDTAB,
     O            ip)
          if(ip /= 0) then
c           We have found a command.  Is it FINISH?
            standard_name = command(cmdval(ip))
            if(standard_name /= 'FINISH') then
c             No, it is not FINISH!
c             Are we already processing a cmd?
              if(cmd_in_process /= ' ') then
c               Yes we are.  Wrap it up and output the result.   The message
c               is complete and ready to process. 
                call md5(message(1:m),
     o                    md5_digest)
         
                ie = i - 1
         
                it = len_trim(mif_with_path)
                write(stdmd5,50) tabid, class, is, ie, cmd_in_process, 
     a                   md5_digest, mif_with_path(1:it) 
         
c               Clear the message buffer
                m = 0
                message = ''
                is = i
                cmd_in_process = cmd
                tabid = ' '
                class = ' '
                call update_message(line,
     m                              m, message)
              else
c               No, this is the first command
                m = 0
                message = ''
                is = i
                cmd_in_process = cmd
                call update_message(line,
     m                              m, message)
                write(stdout,*) ' Trking: First command is:', cmd
                write(stdout,*) 'is=',is, ' m=',m
              endif
            else
c             We have found the end of the input.  Wrap up the md5 file, and 
c             close the scratch file. 
              call free_unit(0, stdscr)
              call free_unit(0, stdmd5)
              return
            endif


          else
c           We do not have a command.
            if(cmd_in_process /= ' ') then
c             We are processing a command. 
              call update_message(line,
     m                            m, message)
              it = index(line,'TABID')
              if(it > 0) then
c               We have a tabid.  Extract it.
                line = line(it+5:)
                it = index(line,'=')
                if(it >0) then
                  line = adjustl(line(it+1:))
                  it = index(line,' ')
                  tabid = line(1:it-1)
                  class = 'user'
                  if(tabid == '-1') then
                    tabid = ' '
                    class = ' '
                  endif
                  if(tabid(1:1) == '-') then
                    tabid = tabid(2:)
                  endif
                endif
              endif
            else
              write(stdout,*) ' Skipping line;',line
            endif
          endif
        endif
      end do

      end

              
          

