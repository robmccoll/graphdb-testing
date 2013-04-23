#ifndef  STATIC_CLUSTERING_H
#define  STATIC_CLUSTERING_H

#include "stinger.h"
#include "stinger-return.h"
#include "stinger-workflow.h"

#include <stdint.h>

typedef struct static_clustering_workpace {
  uint8_t print_count;
  uint8_t histogram_sizes_to_file;
  char *  path;
  char *  filename;
  double * scores;
  uint64_t labels[0];
} static_clustering_workpace_t;


int64_t
parallel_shiloach_vishkin_clustering (struct stinger * S, int64_t nv,
                                      int64_t * component_map);

stinger_return_t
static_clustering_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
static_clustering_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
static_clustering_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count);

#endif  /*STATIC_CLUSTERING_H*/
