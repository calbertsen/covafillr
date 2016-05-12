AC_INIT([covafillr], [0])

AC_SUBST(INSTALLED_WITH_JAGS,"FALSE")
AC_SUBST(JAGS_ROOT,"")
AC_SUBST(JAGS_INCL,"")
AC_SUBST(JAGS_LIBS,"")
AC_SUBST(JAGS_CXX,"")

AC_PATH_PROG([PKG_CONFIG], [pkg-config], [], [$PATH:/usr/local/bin])

AC_ARG_WITH([jags],
AS_HELP_STRING([--with-jags],[Will install covafillr with links to JAGS. This will enable the covafill JAGS module to be used from R.]),
[AS_IF([test [ -n "${PKG_CONFIG}" ]],
	    [AS_IF([${PKG_CONFIG} --exists jags],
	      [AC_SUBST(JAGS_LIBS,[`${PKG_CONFIG} --libs jags`])
	      AC_SUBST(JAGS_CXX,[`${PKG_CONFIG} --cflags jags`])
	      AC_MSG_NOTICE([Link flags are $JAGS_LIBS])
	      AS_IF([test [ -z "${JAGS_LIBS}" -a -z "${JAGS_CXX}" ]],
	        [AC_MSG_WARN([Cannot get information about jags from pkg-config])],
		[AC_MSG_NOTICE([Setting compile and link flags according to pkg-config])
                 AC_MSG_NOTICE([Compile flags are $JAGS_CXX])
                 AC_MSG_NOTICE([Link flags are $JAGS_LIBS])
		 AC_SUBST(INSTALLED_WITH_JAGS,["TRUE"])
 		 AC_SUBST(JAGS_INCL,["#include <covafill/JAGS>"])
		 ]
		 )
	      ],
	      [
	      AC_MSG_WARN([pkg-config file for jags unavailable])
	    ]
	    )],
	    [	    
	    ]
	    )
	    ],
	    [AC_MSG_NOTICE([Installing without JAGS])]
	    )
	    

AC_CONFIG_FILES(src/Makevars)
AC_CONFIG_FILES(src/jags.cpp)
AC_CONFIG_FILES(R/zzz.R)
AC_OUTPUT
