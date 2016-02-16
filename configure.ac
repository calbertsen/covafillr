AC_INIT([covafillr], [0])

INSTALLED_WITH_JAGS="FALSE"
JAGS_ROOT=""
JAGS_INCL="" 

AC_PATH_PROG([PKG_CONFIG], [pkg-config], [], [$PATH:/usr/local/bin])

if test [ -n "${PKG_CONFIG}" ] ; then
   if `${PKG_CONFIG} --exists jags`; then
         JAGS_LIBS=`${PKG_CONFIG} --libs jags`
         JAGS_CXX=`${PKG_CONFIG} --cflags jags`
         if test [ -z "${JAGS_LIBS}" -a -z "${JAGS_CXX}" ] ; then
            AC_MSG_WARN([Cannot get information about jags from pkg-config])
         else
            AC_MSG_NOTICE([Setting compile and link flags according to pkg-config])
            AC_MSG_NOTICE([Compile flags are $JAGS_CFLAGS])
            AC_MSG_NOTICE([Link flags are $JAGS_LIBS])
            INSTALLED_WITH_JAGS="TRUE"
	    JAGS_INCL="#include <covafill/JAGS>"
         fi
   else
      JAGS_LIBS=``
      JAGS_CXX=``
      AC_MSG_WARN([pkg-config file for jags unavailable])
   fi
else
	JAGS_LIBS=``
	JAGS_CXX=``
fi

AC_SUBST(JAGS_LIBS)
AC_SUBST(JAGS_CXX)
AC_SUBST(INSTALLED_WITH_JAGS)
AC_SUBST(JAGS_INCL)
AC_SUBST(JAGS_ROOT)

AC_CONFIG_FILES(src/Makevars)
AC_CONFIG_FILES(src/jags.cpp)
AC_CONFIG_FILES(R/zzz.R)
AC_OUTPUT
