#ifndef  ALGNAME_H
#define  ALGNAME_H

typedef struct ALGNAME_workspace {
} ALGNAME_workspace_t;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
stinger_return_t
ALGNAME_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
ALGNAME_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
ALGNAME_before_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count);

stinger_return_t
ALGNAME_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count);

#endif  /*ALGNAME_H*/
