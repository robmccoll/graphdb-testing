#ifndef  _OMP_H
#define  _OMP_H

#if defined(_OPENMP)
#define OMP(x) _Pragma(x)
#else
#define OMP(x)
#endif

#endif  /*_OMP_H*/
