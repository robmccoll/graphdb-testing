#include "binary_stream.h"
#include "stinger-utils.h"
#include "xmalloc.h"

binary_stream_t *
binary_stream_from_void(stinger_t * S, void ** workspace) {
  binary_stream_t * stream = *((binary_stream_t **)workspace);
  if(!stream) {
    stream = binary_stream_new("initial-graph.bin", "action-stream.bin", 10000, 10);
    *workspace = stream;
  }
  return stream;
}

binary_stream_t *
binary_stream_new(char * init_file, char * stream_file, uint64_t edges_per_batch, uint64_t num_batches) {
  binary_stream_t * stream = xmalloc(sizeof(binary_stream_t));
  if(stream) {
    stream->batch_size = edges_per_batch;
    stream->nbatch = num_batches;
    stream->init_file = init_file;
    stream->stream_file = stream_file;
    stream->actno = 0;
  }
  return stream;
}

stinger_return_t
binary_stream_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  binary_stream_t * stream = binary_stream_from_void(S, workspace);

  snarf_graph (stream->init_file, &stream->nv, &stream->ne, &stream->off, &stream->ind, &stream->weight, &stream->graphmem);

  snarf_actions (stream->stream_file, &stream->naction, &stream->action, &stream->actionmem);

  stinger_set_initial_edges (S, stream->nv, 0, stream->off, stream->ind, stream->weight, NULL, NULL, 0);

  if(stream->naction < stream->nbatch * stream->batch_size) {			
    fprintf (stderr, "WARNING: not enough actions\n");	
    stream->nbatch = (stream->naction + stream->batch_size - 1) / stream->batch_size;		
  }								

  OMP("omp parallel for")
  for(uint64_t v = 0; v < stream->nv; v++) {
    stinger_vtype_set(S, v, 1);
  }

  free(stream->graphmem);

  return STINGER_SUCCESS;
}

stinger_return_t
binary_stream_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t ** actions, int64_t * count) {
  binary_stream_t * stream = binary_stream_from_void(S, workspace);

  if(stream->actno < stream->nbatch * stream->batch_size) {

    const int64_t endact = (stream->actno + stream->batch_size > stream->naction ? stream->naction : stream->actno + stream->batch_size);
    int64_t *acts = &stream->action[2*stream->actno];
    int64_t numActions = endact - stream->actno;

    if(*count < numActions) {
      *actions = xrealloc(*actions, sizeof(edge_action_t) * numActions);
      if(!(*actions))
	return STINGER_ALLOC_FAILED;
    }
    *count = numActions;

    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    OMP("omp parallel for")
    for(uint64_t k = 0; k < numActions; k++) {
      const int64_t i = acts[2 * k];
      const int64_t j = acts[2 * k + 1];
      
      (*actions)[k].type	  = 0;
      (*actions)[k].source	  = i;
      (*actions)[k].dest	  = j;
      (*actions)[k].weight	  = 1;
      (*actions)[k].time	  = batch;
    }

    stream->actno += stream->batch_size;
    if(stream->actno >= stream->nbatch * stream->batch_size || stream->actno >= stream->naction)
      return STINGER_REMOVE;
    else
      return STINGER_SUCCESS;
  } else {
    return STINGER_REMOVE;
  }
}

void
binary_stream_help(void ** workspace) {
  printf(
"Stream: Binary stream\n"
"=====================\n\n"
"Loads the graph from binary format graph and actions files created by \n"
"gen-streams.  The user may specify the number of batches to load and the size\n"
"of each batch.\n\n"
"Options:\n"
"\tbatch_size:	   Number of edges per batch\n"
"\tnbatch:	   Number of batches\n"
"\tinit_file:	   Name of the initial graph file\n"
"\tstream_file:    Name of the stream actions file\n");
exit(-1);
}

stinger_return_t
binary_stream_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  binary_stream_t * ws = binary_stream_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  CODE_TOO=1 ./tools/stringset/stringset help "binary_stream_help(workspace);" \
  ? "binary_stream_help(workspace);"                                           \
  -? "binary_stream_help(workspace);"                                          \
  --? "binary_stream_help(workspace);"                                         \
  -h "binary_stream_help(workspace);"                                          \
  --help "binary_stream_help(workspace);"                                      \
  batch_size "ws->batch_size = atol(values[i]);" \
  nbatch "ws->nbatch = atol(values[i]);" \
  init_file "ws->init_file = values[i];" \
  stream_file "ws->stream_file = values[i];" \

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
			binary_stream_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  binary_stream_help(workspace);
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
		  binary_stream_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  binary_stream_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    binary_stream_help(workspace);
	  }
	} break;
      case 'b':
	{
	  str++; len--;
	  if(!strncmp(str, "atch_size", len)) {
	    str += 9; len -= 9;
	    if(len == 0) {
	      /* batch_size */
	      ws->batch_size = atol(values[i]);
	    }
	  }
	} break;
      case 'h':
	{
	  str++; len--;
	  if(!strncmp(str, "elp", len)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* help */
	      binary_stream_help(workspace);
	    }
	  }
	} break;
      case 'i':
	{
	  str++; len--;
	  if(!strncmp(str, "nit_file", len)) {
	    str += 8; len -= 8;
	    if(len == 0) {
	      /* init_file */
	      ws->init_file = values[i];
	    }
	  }
	} break;
      case 'n':
	{
	  str++; len--;
	  if(!strncmp(str, "batch", len)) {
	    str += 5; len -= 5;
	    if(len == 0) {
	      /* nbatch */
	      ws->nbatch = atol(values[i]);
	    }
	  }
	} break;
      case 's':
	{
	  str++; len--;
	  if(!strncmp(str, "tream_file", len)) {
	    str += 10; len -= 10;
	    if(len == 0) {
	      /* stream_file */
	      ws->stream_file = values[i];
	    }
	  }
	} break;
    }

/* END GENERATED CODE */
  }

  free(keys); free(values);
  return STINGER_SUCCESS;
}
