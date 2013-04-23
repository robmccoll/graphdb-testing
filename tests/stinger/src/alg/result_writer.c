#include "stinger-workflow.h"
#include "xmalloc.h"
#include "result_writer.h"

result_writer_workspace_t *
result_writer_workspace_from_void(stinger_t * S, void ** workspace) {
  result_writer_workspace_t * ws = *((result_writer_workspace_t **)workspace);
  if(!ws) {
    ws = xmalloc(sizeof(result_writer_workspace_t) + sizeof(int64_t) *  
      stinger_vertices_max_vertices_get(stinger_vertices_get(S)));
    if(!ws) return NULL;
    *workspace = ws;
    ws->path = "./";
  }
  return ws;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
stinger_return_t
result_writer_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace)
{
  result_writer_workspace_t * ws =
    result_writer_workspace_from_void(S, workspace);

  stinger_workflow_write_named_results(wkflow, ws->path, 0);
}

void
result_writer_help(void ** workspace)
{
  printf(
    "Algorithm: Result Writer\n"
    "========================\n\n"
    "Writes all of the named results in the workflow to files every iteration.\n\n"
    "Options:\n"
    "\tpath:         path to store results files\n");
  exit(-1);
}

stinger_return_t
result_writer_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings)
{
  char ** keys = NULL, ** values = NULL;
  int num = 0;
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);
  result_writer_workspace_t * ws =
    result_writer_workspace_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);

    /*
      INPUT TO STRINGSET TO GENERATE CODE BELOW:

      CODE_TOO=1 ./tools/stringset/stringset \
      help "result_writer_help(workspace);"\
      ? "result_writer_help(workspace);"\
      -? "result_writer_help(workspace);"\
      --? "result_writer_help(workspace);"\
      -h "result_writer_help(workspace);"\
      --help "result_writer_help(workspace);"\
      path "ws->path = values[i];"\

    */

/* GENERATED WITH stringset.c */
    if(len) switch(*str) {
      case '-':
	{
	  str++; len--;
	  if(len) switch(*str) {
	    case '-':
	      {
		str++; len--;
		if(len) switch(*str) {
		  case '?':
		    {
		      str++; len--;
		      if(len == 0) {
			/* --? */
			result_writer_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  result_writer_help(workspace);
			}
		      }
		    } break;
		}
	      } break;
	    case '?':
	      {
		str++; len--;
		if(len == 0) {
		  /* -? */
		  result_writer_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  result_writer_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    result_writer_help(workspace);
	  }
	} break;
      case 'h':
	{
	  str++; len--;
	  if(!strncmp(str, "elp", len)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* help */
	      result_writer_help(workspace);
	    }
	  }
	} break;
      case 'p':
	{
	  str++; len--;
	  if(!strncmp(str, "ath", len)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* path */
	      ws->path = values[i];
	    }
	  }
	} break;
    }
    /* END GENERATED CODE */
  }

  free(keys);
  free(values);
  return STINGER_SUCCESS;
}

stinger_return_t
result_writer_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count)
{
  batch += 1;
  result_writer_workspace_t * ws =
    result_writer_workspace_from_void(S, workspace);

  stinger_workflow_write_named_results(wkflow, ws->path, batch);
}
