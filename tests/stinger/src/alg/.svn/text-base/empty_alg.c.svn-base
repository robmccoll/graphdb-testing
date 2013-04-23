#include "stinger-workflow.h"
#include "xmalloc.h"
#include "ALGNAME.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ * 
 * IMPLEMENTATION FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* TODO Put any implementation support functions here.  You can expose them in
 * the header if you want to. */

ALGNAME_workspace_t *
ALGNAME_workspace_from_void(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace)
{
  ALGNAME_workspace_t * ws =
    *((ALGNAME_workspace_t**)workspace);

  if(!ws) {
    ws = xmalloc(sizeof(ALGNAME_workspace_t));
    if(!ws) return NULL;

    *workspace = ws;
    ws->nv = stinger_vertices_max_vertices_get(stinger_vertices_get(S));
  }

  return ws;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
stinger_return_t
ALGNAME_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace)
{
  ALGNAME_workspace_t * ws =
    ALGNAME_workspace_from_void(S, wkflow, workspace);
  /* TODO */
}

void
ALGNAME_help(void ** workspace)
{
  printf(
    /* TODO */
    "Algorithm: \n"
    "======================================\n\n"
    "The scores are recalculated from scratch each iteration.\n\n"
    "Options:\n"
    "\tpath:         path to store results files\n"
    "\tfilename:     prefix for the file\n");
  exit(-1);
}

stinger_return_t
ALGNAME_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings)
{
  char ** keys = NULL, ** values = NULL;
  int num = 0;
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);
  ALGNAME_workspace_t * ws =
    ALGNAME_workspace_from_void(S, wkflow, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);

    /*
      INPUT TO STRINGSET TO GENERATE CODE BELOW:
      TODO

      CODE_TOO=1 ./tools/stringset/stringset \
      help "ALGNAME_help(workspace);"\
      ? "ALGNAME_help(workspace);"\
      -? "ALGNAME_help(workspace);"\
      --? "ALGNAME_help(workspace);"\
      -h "ALGNAME_help(workspace);"\
      --help "ALGNAME_help(workspace);"\
      path "ws->path = values[i];"\
      filename "ws->filename = values[i];"\

    */
    /* GENERATED WITH stringset.c */
    /* END GENERATED CODE */
  }

  free(keys);
  free(values);
  return STINGER_SUCCESS;
}

stinger_return_t
ALGNAME_before_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count)
{
  batch += 1;
  ALGNAME_workspace_t * ws =
    ALGNAME_workspace_from_void(S, wkflow, workspace);

  /* TODO */
}

stinger_return_t
ALGNAME_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count)
{
  batch += 1;
  ALGNAME_workspace_t * ws =
    ALGNAME_workspace_from_void(S, wkflow, workspace);

  /* TODO */
}
