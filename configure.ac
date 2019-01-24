# Perlude.
AC_INIT([vcfstatsalive], [1.0.0], [atrusqiao@gmail.com])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#Fixing file creation time
find ${srcdir} -name aclocal.m4 -exec touch {} \;
find ${srcdir} -name \*.in -exec touch {} \;
find ${srcdir} -name configure -exec touch {} \;

# Checks for programs
AC_PROG_CXX
AM_PROG_AR
LT_PREREQ([2.2])
LT_INIT()
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_14

AC_CONFIG_HEADERS([config.h])

CFLAGS=
CPPFLAGS=
CXXFLAGS=

## START: DEBUG / RELEASE ##
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging, default: no])],
              [case "${enableval}" in
                yes) debug=true ;;
                no)  debug=false ;;
                *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
               esac],
              [debug=false])

AM_CONDITIONAL([DEBUG], [test x"$debug" = x"true"])
AS_IF([test "x$debug" = xtrue], [AC_DEFINE([LOGLVL], [LOGLV_DEBUG], [Define log level to be debug])])
## END: DEBUG / RELEASE ##

## START: PROFILING ##
AC_ARG_ENABLE([profiling],
              [AS_HELP_STRING([--enable-profiling],
                              [enable profiling, default: no])],
              [case "${enableval}" in
                yes) profiling=true ;;
                no)  profiling=false ;;
                *)   AC_MSG_ERROR([bad value ${enableval} for --enable-profiling]) ;;
               esac],
              [profiling=false])

AM_CONDITIONAL([PROFILING], [test x"$profiling" = x"true"])
## END: DEBUG / RELEASE ##


## START: Checking htslib availability ##
AC_ARG_WITH([htslib],
            [AS_HELP_STRING([--with-htslib=<path-to-htslib>],
                            [Specify the path to which htslib this program will link against])])
AC_ARG_WITH([htslib-static],
            [AS_HELP_STRING([--with-htslib-static=<path-to-libhts.a>],
                            [Specify the path to which libhts.a this program will link against])])
AC_ARG_WITH([htslib-header],
            [AS_HELP_STRING([--with-htslib-header=<path-to-htslib-headers>],
                            [Specify the path in which htslib headers are searched])])

AS_IF([test "x$with_htslib" != x -a "x$with_htslib_static" != x], [AC_MSG_ERROR([--with-htslib and --with-htslib-static cannot be specified together])])

AS_IF([test "x$with_htslib" != x], [CPPFLAGS="$CPPFLAGS -I$with_htslib/include"])
AS_IF([test "x$with_htslib" != x], [LDFLAGS="$LDFLAGS -L$with_htslib/lib"])


AS_IF([test "x$with_htslib_static" != x], [LIBS="$with_htslib_static $LIBS"])
AS_IF([test "x$with_htslib_static" != x], [AS_IF([test "x$with_htslib_header" != x],
      [CPPFLAGS="-I $with_htslib_header $CPPFLAGS"],
      [AC_MSG_ERROR([--with-htslib-static and --with-htslib-header need to be specified together])])
])
## END: Checking htslib availability ##

## START: Checking libjansson availibility ##
AC_ARG_WITH([jansson],
            [AS_HELP_STRING([--with-jansson=<path-to-jansson>],
                            [Specify the path to which lib-jansson this program will link against])])
AC_ARG_WITH([jansson-static],
            [AS_HELP_STRING([--with-jansson-static=<path-to-libjansson.a>],
                            [Specify the path to which libjansson.a this program will link against])])
AC_ARG_WITH([jansson-header],
            [AS_HELP_STRING([--with-jansson-header=<path-to-libjansson-headers>],
                            [Specify the path in which jansson headers are searched])])

AS_IF([test "x$with_jansson" != x -a "x$with_jansson_static" != x],
      [AC_MSG_ERROR([--with-jansson and --with-jansson-static cannot be specified together])])

AS_IF([test "x$with_jansson" != x], [CPPFLAGS="$CPPFLAGS -I$with_jansson/include"])
AS_IF([test "x$with_jansson" != x], [LDFLAGS="$LDFLAGS -L$with_jansson/lib"])

AS_IF([test "x$with_jansson_static" != x], [LIBS="$with_jansson_static $LIBS"])
AS_IF([test "x$with_jansson_static" != x], [AS_IF([test "x$with_jansson_header" != x],
      [CPPFLAGS="-I $with_jansson_header $CPPFLAGS"],
      [AC_MSG_ERROR([--with-jansson-static and --with-jansson-header need to be specified together])])
])
## END: Checking libjansson availability ##

echo "$CPPFLAGS"

## Checks for headers
AC_CHECK_HEADERS([htslib/hts.h htslib/vcf.h htslib/kstring.h htslib/kseq.h], [], [AC_MSG_ERROR([Could not find htslib header files])])
AC_CHECK_HEADER([jansson.h], [], [AC_MSG_ERROR([Could not find jansson header files])])

## Checks for libraries
AC_CHECK_LIB([m], [floor])
AS_IF([test "x$with_htslib_static" == x], [
       AC_CHECK_LIB([hts], [bcf_read], [], [AC_MSG_ERROR([Could not find htslib shared library. Please specify htslib installation path using --with-htslib, or link statically using --with-htslib-static])])
       ])
AS_IF([test "x$with_jansson_static" == x], [
       AC_CHECK_LIB([jansson], [json_dumps], [], [AC_MSG_ERROR([Could not find jansson shared library. Please specify jansson installation path using --with-jansson, or link statically using --with-jansson-static])])
       ])

## Config files
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

echo \
"----------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'
Compiler: '${CXX}'
CFLAGS: '${CFLAGS}'
CPPFLAGS: '${CPPFLAGS}'
CXXFLAGS: '${CXXFLAGS}'

Extra variables:
LDFLAGS='${LDFLAGS}'
LDADD='${LDADD}'
LIBS='${LIBS}'

PTHREAD_LIBS='${PTHREAD_LIBS}'
PTHREAD_CFLAGS='${PTHREAD_CFLAGS}'
PTHREAD_CC='${PTHREAD_CC}'

DATADIR=${DATADIR}
debug=${debug}
"
