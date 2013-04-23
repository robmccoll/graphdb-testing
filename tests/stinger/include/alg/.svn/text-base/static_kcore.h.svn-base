#ifndef  STATIC_KCORE_H
#define  STATIC_KCORE_H

#include "stinger-workflow.h"
#include "xmalloc.h"

typedef struct static_kcore_workpace {
  uint8_t print_k;
  uint8_t histogram_k_to_file;
  uint8_t histogram_k_sizes_to_file;
  char *  path;
  char *  filename;
  int64_t nv;
  int64_t * labels;
  int64_t * counts;
  int64_t k;
} static_kcore_workpace_t;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * IMPLEMENTATION FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void
kcore_find(stinger_t *S, int64_t * labels, int64_t * counts, int64_t nv, int64_t * k_out);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

stinger_return_t
static_kcore_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
static_kcore_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
static_kcore_before_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count);

stinger_return_t
static_kcore_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count);

#endif  /*STATIC_KCORE_H*/
