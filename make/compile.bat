ECHO OFF
copy ..\*.for
copy ..\*.cmn
copy ..\*.prm
copy ..\..\share\*.for
copy ..\..\share\*.f90
copy ..\..\share\*.cmn
copy ..\..\share\arsize.prm
copy ..\..\share\morgux.prm  morg.prm
copy ..\..\share\lf95_ms_s_npf.who lf95_ms_s_npf.for
copy ..\..\svn_report\svn_reportdmy.for svn_report.for
copy pwd_msw_lf95.for pwd.for
del fqshrarg.for
del locsubux.for
del locsubs.for
del timer90.for
del pwd_msw_lf95.for
del pwd_msw_g95.for
del pwd_lx_lf95.for
del pwd_lx_g95.for
del getsvn_lx_g95.for
del getsvn_lx_lf95.for
del getsvn_msw_g95.for
call am > am.out
del *.for
del *.f90
del *.cmn
del *.prm
del *.map
del *.obj

