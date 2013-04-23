#include "static_kcore.h"
#include "histogram.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ * 
 * IMPLEMENTATION FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* TODO Put any implementation support functions here.  You can expose them in
 * the header if you want to. */

void
kcore_find(stinger_t *S, int64_t * labels, int64_t * counts, int64_t nv, int64_t * k_out) {
  int64_t k = 0;
  
  for(int64_t v; v < nv; v++) {
    if(stinger_outdegree_get(S,v)) {
      labels[v] = 1;
    } else {
      labels[v] = 0;
    }
  }

  int changed = 1;
  while(changed) {
    changed = 0;
    k++;

    OMP("omp parallel for")
    for(int64_t v = 0; v < nv; v++) {
      if(labels[v] == k) {
	int64_t count = 0;
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
	  count += (labels[STINGER_EDGE_DEST] >= k);
	} STINGER_FORALL_EDGES_OF_VTX_END();
	if(count > k)
	  counts[v] = count;
      }
    }

    OMP("omp parallel for")
    for(int64_t v = 0; v < nv; v++) {
      if(labels[v] == k) {
	int64_t count = 0;
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
	  count += (labels[STINGER_EDGE_DEST] >= k && 
	    counts[STINGER_EDGE_DEST] > k);
	} STINGER_FORALL_EDGES_OF_VTX_END();
	if(count > k) {
	  labels[v] = k+1;
	  changed = 1;
	}
      }
    }
  }

  *k_out = k;
}

static_kcore_workpace_t *
static_kcore_workspace_from_void(stinger_t * S, void ** workspace) {
  static_kcore_workpace_t * ws = *((static_kcore_workpace_t **)workspace);
  if(!ws) {
    ws = xmalloc(sizeof(static_kcore_workpace_t));
    if(!ws) return NULL;
    *workspace = ws;
    ws->nv = stinger_vertices_max_vertices_get(stinger_vertices_get(S));
    ws->path = "./";
    ws->filename = "kcore";
    ws->print_k = 1;
    ws->histogram_k_to_file = 1;
    ws->histogram_k_sizes_to_file = 1;
    ws->labels = xmalloc(sizeof(uint64_t) * ws->nv);
    ws->counts = xmalloc(sizeof(uint64_t) * ws->nv);
  }
  return ws;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
stinger_return_t
static_kcore_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace)
{
  static_kcore_workpace_t * ws =
    static_kcore_workspace_from_void(S, workspace);

  kcore_find(S, ws->labels, ws->counts, ws->nv, &(ws->k));

  if(ws->print_k) {
    printf("Static k-core: k = %ld after %d\n", ws->k, 0);
  }
  if(ws->histogram_k_to_file) {
    char filename[1024];
    sprintf(filename, "%s.%s", ws->filename, "klabels");
    histogram_int64(S,ws->labels,ws->nv,ws->path,filename,0);
  }
  if(ws->histogram_k_sizes_to_file) {
    char filename[1024];
    sprintf(filename, "%s.%s", ws->filename, "sizes");
    histogram_int64(S,ws->counts,ws->nv,ws->path,filename,0);
  }
  return STINGER_SUCCESS;
}

void
static_kcore_help(void ** workspace)
{
  printf(
    /* TODO */
    "Algorithm: K-Core\n\n"
    "=================\n\n"
    "Finds the largest k-core(s) of the graph by iterating over the graph. The result\n"
    "is an array of labels and counts where the label of each vertex will be the k of \n"
    "the largest k-core to which that vertex belongs, and the count will be the number\n"
    "of neighbors also in that k-core.\n\n"
    "Options:\n"
    "\tprint_k:      1 or 0, print out the largest k found\n"
    "\thistogram_k_to_file: histogram the max k of the vertices\n"
    "\tpath:         path to store results files\n"
    "\tfilename:     prefix for the file\n");
  exit(-1);
}

stinger_return_t
static_kcore_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings)
{
  char ** keys = NULL, ** values = NULL;
  int num = 0;
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);
  static_kcore_workpace_t * ws =
    static_kcore_workspace_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);

    /*
      INPUT TO STRINGSET TO GENERATE CODE BELOW:
      TODO

      CODE_TOO=1 ./tools/stringset/stringset \
      help "static_kcore_help(workspace);"\
      ? "static_kcore_help(workspace);"\
      -? "static_kcore_help(workspace);"\
      --? "static_kcore_help(workspace);"\
      -h "static_kcore_help(workspace);"\
      --help "static_kcore_help(workspace);"\
      path "ws->path = values[i];"\
      filename "ws->filename = values[i];"\
      print_k "ws->print_k = atol(values[i]);" \
      histogram_k_to_file "ws->histogram_k_to_file = atol(values[i]);" \
      histogram_k_sizes_to_file "ws->histogram_k_sizes_to_file = atol(values[i]);" \

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
			static_kcore_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  static_kcore_help(workspace);
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
		  static_kcore_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  static_kcore_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    static_kcore_help(workspace);
	  }
	} break;
      case 'f':
	{
	  str++; len--;
	  if(!strncmp(str, "ilename", len)) {
	    str += 7; len -= 7;
	    if(len == 0) {
	      /* filename */
	      ws->filename = values[i];
	    }
	  }
	} break;
      case 'h':
	{
	  str++; len--;
	  if(len) switch(*str) {
	    case 'e':
	      {
		str++; len--;
		if(!strncmp(str, "lp", len)) {
		  str += 2; len -= 2;
		  if(len == 0) {
		    /* help */
		    static_kcore_help(workspace);
		  }
		}
	      } break;
	    case 'i':
	      {
		str++; len--;
		if(!strncmp(str, "stogram_k_", 10)) {
		  str += 10; len -= 10;
		  if(len) switch(*str) {
		    case 's':
		      {
			str++; len--;
			if(!strncmp(str, "izes_to_file", len)) {
			  str += 12; len -= 12;
			  if(len == 0) {
			    /* histogram_k_sizes_to_file */
			    ws->histogram_k_sizes_to_file = atol(values[i]);
			  }
			}
		      } break;
		    case 't':
		      {
			str++; len--;
			if(!strncmp(str, "o_file", len)) {
			  str += 6; len -= 6;
			  if(len == 0) {
			    /* histogram_k_to_file */
			    ws->histogram_k_to_file = atol(values[i]);
			  }
			}
		      } break;
		  }
		}
	      } break;
	  }
	} break;
      case 'p':
	{
	  str++; len--;
	  if(len) switch(*str) {
	    case 'a':
	      {
		str++; len--;
		if(!strncmp(str, "th", len)) {
		  str += 2; len -= 2;
		  if(len == 0) {
		    /* path */
		    ws->path = values[i];
		  }
		}
	      } break;
	    case 'r':
	      {
		str++; len--;
		if(!strncmp(str, "int_k", len)) {
		  str += 5; len -= 5;
		  if(len == 0) {
		    /* print_k */
		    ws->print_k = atol(values[i]);
		  }
		}
	      } break;
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
static_kcore_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count)
{
  batch += 1;
  static_kcore_workpace_t * ws =
    static_kcore_workspace_from_void(S, workspace);

  kcore_find(S, ws->labels, ws->counts, ws->nv, &(ws->k));

  if(ws->print_k) {
    printf("Static k-core: k = %ld after %ld\n", ws->k, batch);
  }
  if(ws->histogram_k_to_file) {
    char filename[1024];
    sprintf(filename, "%s.%s", ws->filename, "klabels");
    histogram_int64(S,ws->labels,ws->nv,ws->path,filename,batch);
  }
  if(ws->histogram_k_sizes_to_file) {
    char filename[1024];
    sprintf(filename, "%s.%s", ws->filename, "sizes");
    histogram_int64(S,ws->counts,ws->nv,ws->path,filename,batch);
  }
  return STINGER_SUCCESS;
}
