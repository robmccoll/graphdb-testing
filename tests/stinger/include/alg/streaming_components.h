#ifndef  STREAMING_COMPONENTS_H
#define  STREAMING_COMPONENTS_H

#include "stinger.h"
#include "stinger-return.h"
#include "stinger-workflow.h"

#include <stdint.h>

typedef struct streaming_components_workpace {
  uint64_t nv;
  uint8_t print_count;
  uint8_t histogram_sizes_to_file;
  char *  path;
  char *  filename;

  int64_t * queue;
  int64_t * level;
  int64_t * found;
  int64_t * same_level_queue;

  int64_t   parentsPerVertex;
  int64_t * parentArray;
  int64_t * parentCounter;
  int64_t * bfs_components;
  int64_t * bfs_component_sizes;
} streaming_components_workpace_t;

stinger_return_t
streaming_components_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
streaming_components_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
streaming_components_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count);


#endif  /*STREAMING_COMPONENTS_H*/
