#ifndef  STATIC_PAGERANK_H
#define  STATIC_PAGERANK_H

#include <stdint.h>

#include "stinger.h"

int64_t
pagerank(stinger_t * S, int64_t NV, double * pr, double epsilon, double dampingfactor, int64_t maxiter);

#endif  /*STATIC_PAGERANK_H*/
