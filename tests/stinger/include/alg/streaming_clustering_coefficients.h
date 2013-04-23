#ifndef  STREAMING_CLUSTERING_COEFFICIENTS_H
#define  STREAMING_CLUSTERING_COEFFICIENTS_H

#include "stinger.h"
#include "stinger-return.h"
#include "stinger-workflow.h"

#include <stdint.h>

typedef struct streaming_clustering_coefficients_workpace {
  char * ntri_name;
  char * local_cc_name;
  stinger_named_result_t * ntri_nr;
  stinger_named_result_t * local_cc_nr;

  int64_t * ntri;
  double  * local_cc;
  int64_t * affected;
  int64_t global_ntri;

  uint64_t nv;
  uint8_t print_global;

} streaming_clustering_coefficients_workpace_t;

stinger_return_t
streaming_clustering_coefficients_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
streaming_clustering_coefficients_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
streaming_clustering_coefficients_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, 
  int64_t batch, edge_action_t * actions, int * result, int64_t count);

#endif  /*STREAMING_CLUSTERING_COEFFICIENTS_H*/
