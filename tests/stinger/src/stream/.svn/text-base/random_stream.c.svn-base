#include "random_stream.h"
#include "xmalloc.h"

random_stream_t *
random_stream_from_void(stinger_t * S, void ** workspace) {
  random_stream_t * stream = *((random_stream_t **)workspace);
  if(!stream) {
    stream = random_stream_new(1000, 1000, 1000, 10);
  }
  return stream;
}

random_stream_t *
random_stream_new(uint64_t edges_init, uint64_t edges_per_batch, uint64_t nv, uint64_t num_batches) {
  random_stream_t * stream = xmalloc(sizeof(random_stream_t));
  if(stream) {
    stream->edges_init	  = edges_init;
    stream->edges_per_batch  = edges_per_batch;
    stream->nv		  = nv;
    stream->num_batches	  = num_batches;
  }
  return stream;
}

stinger_return_t
random_stream_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  random_stream_t * stream = *((random_stream_t **)workspace);
  OMP("omp parallel for")
  for(uint64_t v = 0; v < stream->nv; v++) {
    stinger_vtype_set(S, v, 1);
  }

  OMP("omp parallel for")
  for(uint64_t i = 0; i < stream->edges_init; i++) {
    stinger_incr_edge_pair(S, 0, rand() % stream->nv, rand() % stream->nv, 1, 0);
  }
  return STINGER_SUCCESS;
}

stinger_return_t
random_stream_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t ** actions, int64_t * count) {
  random_stream_t * stream = *((random_stream_t **)workspace);

  if(*count < stream->edges_per_batch) {
    *actions = xrealloc(*actions, sizeof(edge_action_t) * stream->edges_per_batch);
    if(!(*actions))
      return STINGER_ALLOC_FAILED;
  }
  *count = stream->edges_per_batch;

  OMP("omp parallel for")
  for(uint64_t i = 0; i < stream->edges_per_batch; i++) {
    (*actions)[i].type	  = 0;
    (*actions)[i].source  = rand() % stream->nv;
    (*actions)[i].dest	  = rand() % stream->nv;
    (*actions)[i].weight  = 1;
    (*actions)[i].time	  = batch + 1;
  }

  if(batch == (stream->num_batches - 1))
    return STINGER_REMOVE;
  else
    return STINGER_SUCCESS;
}

void
random_stream_help(void ** workspace) {
  printf(
"Stream: Random stream\n"
"=====================n\n"
"Generates random edges for testing purposes (uses rand() - NOT A TRUE \n"
"UNIFORMLY RANDOM GENERATOR).\n\n"
"Options:\n"
"\tedges_init	   Number of initial edges to generate\n"
"\tedges_per_batch Number of edges to generate per batch\n"
"\tnv		   Number of vertices (edge end points will be 0..nv-1)\n"
"\tnum_batches	   Number of batches to generate before stopping\n");
  exit(-1);
}

stinger_return_t
random_stream_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  random_stream_t * ws = random_stream_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  help "random_stream_help(workspace);"
  ? "random_stream_help(workspace);"
  -? "random_stream_help(workspace);"
  --? "random_stream_help(workspace);"
  -h "random_stream_help(workspace);"
  --help "random_stream_help(workspace);"
  edges_init	  "ws->edges_init = atol(values[i]);"
  edges_per_batch "ws->edges_per_batch = atol(values[i]);"
  nv		  "ws->nv = atol(values[i]);"
  num_batches	  "ws->num_batches = atol(values[i]);"
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
			random_stream_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  random_stream_help(workspace);
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
		  random_stream_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  random_stream_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    random_stream_help(workspace);
	  }
	} break;
      case 'e':
	{
	  str++; len--;
	  if(!strncmp(str, "dges_", 5)) {
	    str += 5; len -= 5;
	    if(len) switch(*str) {
	      case 'i':
		{
		  str++; len--;
		  if(!strncmp(str, "nit", len)) {
		    str += 3; len -= 3;
		    if(len == 0) {
		      /* edges_init */
		      printf("we are here: %ld\n", atol(values[i]));
		      ws->edges_init = atol(values[i]);
		    }
		  }
		} break;
	      case 'p':
		{
		  str++; len--;
		  if(!strncmp(str, "er_batch", len)) {
		    str += 8; len -= 8;
		    if(len == 0) {
		      /* edges_per_batch */
		      ws->edges_per_batch = atol(values[i]);
		    }
		  }
		} break;
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
	      random_stream_help(workspace);
	    }
	  }
	} break;
      case 'n':
	{
	  str++; len--;
	  if(len) switch(*str) {
	    case 'u':
	      {
		str++; len--;
		if(!strncmp(str, "m_batches", len)) {
		  str += 9; len -= 9;
		  if(len == 0) {
		    /* num_batches */
		    ws->num_batches = atol(values[i]);
		  }
		}
	      } break;
	    case 'v':
	      {
		str++; len--;
		if(len == 0) {
		  /* nv */
		  ws->nv = atol(values[i]);
		}
	      } break;
	  }
	} break;
    }

/* END GENERATED CODE */
  }

  free(keys); free(values);
  return STINGER_SUCCESS;
}
