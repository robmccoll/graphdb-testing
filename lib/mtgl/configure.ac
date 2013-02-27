#
# Script to be processed by autoconf to produce a configure script
# and config.h header file.
#

###############################################################
# Init Autoconf / Automake
#

AC_PREREQ(2.59)
AC_INIT([mtgl], [1.1.1], [gemacke@sandia.gov], [],
        [https://software.sandia.gov/trac/mtgl/wiki])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.9.6 foreign dist-zip no-define tar-pax nostdinc])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Add a header guard to the config file.  We need this because we include
# the config file in header files by necessity since MTGL is comprised of only
# header files.  This way we don't redefine things when the config file needs
# to be included in multiple header files.
AH_TOP([#ifndef MTGL_CONFIG_H
#define MTGL_CONFIG_H])
AH_BOTTOM([#endif])

# Information on this package.
AC_CONFIG_HEADER([mtgl/mtgl_config.h])

# Host information.
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

# Configure command line options.
MTGL_COMPILER_OPTIONS
MTGL_COMPILER_OPTION

# --with-traceview flag (default=off)
AC_ARG_WITH(traceview,
            AC_HELP_STRING([--with-traceview],[compile with traceivew enabled]),
            [case "${withval}" in
               yes)
                 echo "--with-traceview = yes"
                 traceview=true
                 ;;
               no)
                 echo "--with-traceview = no"
                 traceview=false
                 ;;
               *)
                 AC_MSG_ERROR(bad value in ${withval} for --with-traceview)
                 ;;
            esac],
            [traceview=false]
)
AM_CONDITIONAL([WITH_TRACEVIEW],[ test "$traceview" = true ])

# --with_debug_arraybounds (default=off)
#   Compiles with -DCHECK_BOUNDS to control bounds checking in
#   dynamicarrays.
AC_ARG_WITH(debug_arraybounds,
            AC_HELP_STRING([--with-debug_arraybounds],
                           [Enable array bounds check in DynamicArray]),
            [case "${withval}" in
               yes)
                 echo "--with-debug_arraybounds = yes"
                 debug_arraybounds=true
                 ;;
               no)
                 echo "--with-debug_arraybounds = no"
                 debug_arraybounds=false
                 ;;
               *)
                 AC_MSG_ERROR(bad value in ${withval} for --with-debug_arraybounds)
                 ;;
            esac],
            [debug_arraybounds=false]
)
AM_CONDITIONAL([WITH_DEBUG_ARRAYBOUNDS],[ test "$debug_arraybounds" = true ])

###############################################################
# Checks for programs.
#

# Get rid of any default compiler flags, so we can start with nothing and
# only add what we want.
CFLAGS=''
CXXFLAGS=''

# List cc/c++ first so that we default to building for the xmt
# compute nodes rather than the head nodes.
AC_PROG_CC([cc gcc])
AC_PROG_CC_C_O
AC_PROG_CXX([c++ CC g++])
AC_PROG_CPP

AC_PROG_LN_S
#AM_PROG_AS
#AC_PROG_RANLIB

# Check if we are compiling on an XMT.
AC_CACHE_CHECK([whether targeting MTA/XMT],
               [mtgl_cv_on_mta],
               [AC_TRY_COMPILE([], [ return __MTA__; ],
                               [mtgl_cv_on_mta=yes], [mtgl_cv_on_mta=no])])
AM_CONDITIONAL(ON_MTA, [test "$mtgl_cv_on_mta" = "yes"])

# Check if we are compiling on a mac.
AC_CACHE_CHECK([whether on OS X],
               [mtgl_cv_on_mac],
               [AC_TRY_COMPILE([], [ return __APPLE__; ],
                               [mtgl_cv_on_mac=yes], [mtgl_cv_on_mac=no])])
AM_CONDITIONAL(ON_MAC, [test "$mtgl_cv_on_mac" = "yes"])

# If we are compiling on a unix machine that is not a mac and not an xmt,
# we need to include the rt library.
AS_IF([test "${mtgl_cv_on_mta}" = yes],[],
      [test "${mtgl_cv_on_mac}" = yes],[],
      [LIBS="${LIBS} -lrt"])

# Add the m and prand libraries if we are on the XMT.  Otherwise set the
# optimization level to O2, and search for a couple of libraries.
AS_IF([test "${mtgl_cv_on_mta}" = yes],[LIBS="${LIBS} -lm -lprand" ],
      [CFLAGS='-O2']
      [CXXFLAGS='-O2']
      [AC_SEARCH_LIBS([sin], [m m_ppc])])

###############################################################
# Checks for libraries.
#

MTGL_CHECK_OMP_ATOMICS
MTGL_CHECK_QTHREADS
MTGL_CHECK_OPENMP

###############################################################
# Checks for headers.
#

#AC_HEADER_STDC

###############################################################
# Checks for types.
#

# Checking for the size of types doesn't work on the XMT, so
# we manually set them for that platform.
AS_IF([test "${mtgl_cv_on_mta}" = yes],
      [
        AC_DEFINE([SIZEOF_INT], [8],
                  [The size of 'int', as computed by sizeof.])
        AC_DEFINE([SIZEOF_LONG], [8],
                  [The size of 'long', as computed by sizeof.])
        AC_DEFINE([SIZEOF_LONG_LONG], [8],
                  [The size of 'long long', as computed by sizeof.])
        AC_DEFINE([SIZEOF_UNSIGNED], [8],
                  [The size of 'unsigned', as computed by sizeof.])
        AC_DEFINE([SIZEOF_UNSIGNED_LONG], [8],
                  [The size of 'unsigned long', as computed by sizeof.])
        AC_DEFINE([SIZEOF_UNSIGNED_LONG_LONG], [8],
                  [The size of 'unsigned long long', as computed by sizeof.])
      ],
      [
        AC_CHECK_SIZEOF([int])
        AC_CHECK_SIZEOF([unsigned])
        AC_CHECK_SIZEOF([long])
        AC_CHECK_SIZEOF([unsigned long])
        AC_CHECK_SIZEOF([long long])
        AC_CHECK_SIZEOF([unsigned long long])
      ])

###############################################################
# Checks for structures.
#


###############################################################
# Checks for compiler characteristics.


###############################################################
# Checks for library functions.


###############################################################
# Output of autoconf/automake.
#

# Last minute fix-up so that #include <mtgl/*.hpp> always works.
CPPFLAGS='-I$(top_srcdir)'" $CPPFLAGS"

AC_CONFIG_FILES([Makefile mtgl/Makefile test/Makefile tutorial/Makefile])

AC_MSG_NOTICE(CC=${CC})
AC_MSG_NOTICE(CXX=${CXX})
AC_MSG_NOTICE(CPPFLAGS=${CPPFLAGS})
AC_MSG_NOTICE(CFLAGS=${CFLAGS})
AC_MSG_NOTICE(CXXFLAGS=${CXXFLAGS})
AC_MSG_NOTICE(LDFLAGS=${LDFLAGS})
AC_MSG_NOTICE(LIBS=${LIBS})

AC_CONFIG_COMMANDS_POST(AS_IF([test "${mtgl_cv_on_mta}" = yes],
                              [rm -f a.out.pl],
                              []))

AC_OUTPUT
