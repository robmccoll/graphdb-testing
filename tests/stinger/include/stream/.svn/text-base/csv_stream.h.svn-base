#ifndef  CSV_STREAM_H
#define  CSV_STREAM_H

#include "stinger-workflow.h"

#include <stdint.h>

typedef struct csv_emap {
  int64_t src_field;
  int64_t src_type;
  int64_t dest_field;
  int64_t dest_type;
  int64_t weight_field;
  int64_t time_field;
  int64_t edge_type;
  struct csv_emap * next;
} csv_emap_t;

typedef struct csv_stream {
  uint64_t edges_init;
  uint64_t edges_per_batch;
  uint64_t num_batches;
  char * file;
  char sep;
  FILE * fp;
  csv_emap_t * list;

  char      *   buf;
  uint64_t      bufSize;
  char      **  fields;
  uint64_t  *   lengths;
  uint64_t      fieldsSize;
  uint64_t      count;
} csv_stream_t;

csv_stream_t *
csv_stream_new(char * file, uint64_t edges_init, uint64_t edges_per_batch, uint64_t num_batches);

stinger_return_t
csv_stream_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
csv_stream_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t ** actions, int64_t * count);

stinger_return_t
csv_stream_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

#endif  /*CSV_STREAM_H*/
