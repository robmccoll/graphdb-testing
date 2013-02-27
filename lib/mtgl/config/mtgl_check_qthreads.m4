dnl -*- Autoconf -*-
dnl

AC_DEFUN([MTGL_CHECK_QTHREADS],
[
  AC_ARG_WITH([qthreads],
              [AC_HELP_STRING([--with-qthreads=@<:@DIR@:>@],
                              [Enable qthreads threading support, optionally looking for qthreads in DIR])],
              [], [with_qthreads=no])

  AS_IF([test -n "$with_qthreads" -a "x$with_qthreads" != "xno"],
        [AS_IF([test "x$with_qthreads" != "xyes"],
               [CPPFLAGS="$CPPFLAGS -I$with_qthreads/include"
                LDFLAGS="$LDFLAGS -L$with_qthreads/lib"])
         [CPPFLAGS="$CPPFLAGS -DUSING_QTHREADS"]
         AC_CHECK_LIB([qthread], [qthread_init], [LIBS="$LIBS -lqthread"],
                      [AC_MSG_ERROR([Qthreads could not be found])])])
])

AC_DEFUN([MTGL_CHECK_OPENMP],
[
  AC_ARG_WITH([openmp],
              [AC_HELP_STRING([--with-openmp],
                              [Enable OpenMP threading support])],
              [], [with_openmp=no])

  m4_ifdef([AC_OPENMP],
           [AS_IF([test -n "$with_openmp" -a "x$with_openmp" != "xno"],
                   [AC_OPENMP
                    [CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"]
                    [CFLAGS="$CFLAGS $OPENMP_CFLAGS"]])])
])
