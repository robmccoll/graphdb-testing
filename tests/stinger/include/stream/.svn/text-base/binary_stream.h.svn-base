#ifndef  BINARY_STREAM_H
#define  BINARY_STREAM_H

#include "stinger-workflow.h"

#include <stdint.h>

typedef struct binary_stream {
  int64_t batch_size;
  int64_t nbatch;
  char * init_file;
  char * stream_file;

  int64_t nv, ne, naction;
  int64_t *off;
  int64_t *from;
  int64_t *ind;
  int64_t *weight;
  int64_t *action;
  int64_t actno;

  int64_t *graphmem;
  int64_t *actionmem;
} binary_stream_t;

binary_stream_t *
binary_stream_new(char * init_file, char * stream_file, uint64_t edges_per_batch, uint64_t num_batches);

stinger_return_t
binary_stream_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
binary_stream_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t ** actions, int64_t * count);

stinger_return_t
binary_stream_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

#endif  /*binary_STREAM_H*/
