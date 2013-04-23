#ifndef  RANDOM_STREAM_H
#define  RANDOM_STREAM_H

#include "stinger-workflow.h"

#include <stdint.h>

typedef struct random_stream {
  uint64_t edges_init;
  uint64_t edges_per_batch;
  uint64_t nv;
  uint64_t num_batches;
} random_stream_t;

random_stream_t *
random_stream_new(uint64_t edges_init, uint64_t edges_per_batch, uint64_t nv, uint64_t num_batches);

stinger_return_t
random_stream_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
random_stream_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t ** actions, int64_t * count);

stinger_return_t
random_stream_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

#endif  /*RANDOM_STREAM_H*/
