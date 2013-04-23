#ifndef  RESULT_WRITER_H
#define  RESULT_WRITER_H

typedef struct result_writer_workspace {
  char * path;
} result_writer_workspace_t;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
stinger_return_t
result_writer_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
result_writer_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
result_writer_before_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count);

stinger_return_t
result_writer_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count);

#endif  /*result_writer_H*/
