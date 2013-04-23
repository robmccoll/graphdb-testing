#include "static_pagerank.h"
#include "xmalloc.h"

int64_t
pagerank(stinger_t * S, int64_t NV, double * pr, double epsilon, double dampingfactor, int64_t maxiter)
{
  double * tmp_pr = xcalloc(NV, sizeof(double));

  OMP("omp parallel for")
  for(uint64_t v = 0; v < NV; v++) {
    pr[v] = 1 / ((double)NV);
  }
  
  int64_t iter = maxiter;
  double delta = 1;

  while(delta > epsilon && iter > 0) {
    OMP("omp parallel for")
    for(uint64_t v = 0; v < NV; v++) {
      tmp_pr[v] = 0;

      STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
	tmp_pr[v] += (((double)pr[STINGER_EDGE_DEST]) / 
	  ((double) stinger_outdegree(S, STINGER_EDGE_DEST)));
      } STINGER_FORALL_EDGES_OF_VTX_END();
    }

    OMP("omp parallel for")
    for(uint64_t v = 0; v < NV; v++) {
      tmp_pr[v] = tmp_pr[v] * dampingfactor + (((double)(1-dampingfactor)) / ((double)NV));
    }

    delta = 0;
    OMP("omp parallel for reduction(+:delta)")
    for(uint64_t v = 0; v < NV; v++) {
      double mydelta = tmp_pr[v] - pr[v];
      if(mydelta < 0)
	mydelta = -mydelta;
      delta += mydelta;
    }

    OMP("omp parallel for")
    for(uint64_t v = 0; v < NV; v++) {
      pr[v] = tmp_pr[v];
    }
  }
}
