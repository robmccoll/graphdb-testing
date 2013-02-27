dnl -*- Autoconf -*-
dnl

AC_DEFUN([MTGL_CHECK_OMP_ATOMICS],
[
  AC_CACHE_CHECK([whether compiler supports builtin atomic incr],
                 [mtgl_cv_atomic_incr],
                 [AS_IF([test "${mtgl_cv_on_mta}" = yes],
                        [mtgl_cv_atomic_incr="no"],
                        [AC_LINK_IFELSE([AC_LANG_SOURCE([[
#ifdef HAVE_IA64INTRIN_H
# include <ia64intrin.h>
#elif HAVE_IA32INTRIN_H
# include <ia32intrin.h>
#endif
#include <stdlib.h>
#include <stdint.h> /* for uint64_t */

int main()
{
uint64_t bar=1;
uint64_t foo = __sync_fetch_and_add(&bar, 1);
return foo;
}]])],
                                        [mtgl_cv_atomic_incr="yes"],
                                        [mtgl_cv_atomic_incr="no"])])])

AS_IF([test "${mtgl_cv_atomic_incr}" = yes],
      [AC_DEFINE([MTGL_ATOMIC_INCR], [1],
                 [Define if the compiler supports __sync_fetch_and_add().])])

])
