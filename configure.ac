AC_PREREQ(2.59)
AC_INIT([buRden], 0.1.0)
# AC_CONFIG_MACRO_DIR([m4])

CC=`R CMD config CC`
CPP=`R CMD config CPP`
CFLAGS=`R CMD config CFLAGS`
CPPFLAGS="`R CMD config CPPFLAGS` $CPPFLAGS"
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_LANG(C++)

AC_CHECK_HEADER(boost/accumulators/accumulators.hpp,, [AC_MSG_ERROR([boost/accumulators/accumulators.hpp not found])])
AC_CHECK_HEADER(boost/accumulators/statistics/stats.hpp,, [AC_MSG_ERROR([boost/accumulators/statistics/stats.hpp not found])])
AC_CHECK_HEADER(boost/accumulators/statistics/tail.hpp,, [AC_MSG_ERROR([boost/accumulators/statistics/tail.hpp not found])])

AC_CONFIG_FILES([src/Makevars])
XTRA_CPPFLAGS=$CPPFLAGS
XTRA_LDFLAGS=$LDFLAGS
AC_SUBST(XTRA_CPPFLAGS)
AC_SUBST(XTRA_LDFLAGS)
AC_OUTPUT
