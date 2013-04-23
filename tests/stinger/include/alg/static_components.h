#ifndef  STATIC_COMPONENTS_H
#define  STATIC_COMPONENTS_H

#include "stinger.h"
#include "stinger-return.h"
#include "stinger-workflow.h"

#include <stdint.h>

typedef struct static_components_workpace {
  uint8_t print_count;
  uint8_t histogram_sizes_to_file;
  char *  path;
  char *  filename;
  int64_t components[0];
} static_components_workpace_t;


stinger_return_t
static_components_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
static_components_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
static_components_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count);

int64_t
parallel_shiloach_vishkin_components (struct stinger * S, int64_t nv,
                                      int64_t * component_map);

#endif  /*STATIC_COMPONENTS_H*/
