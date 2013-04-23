#include "static_components.h"
#include "histogram.h"
#include "xmalloc.h"

int64_t
parallel_shiloach_vishkin_components (struct stinger * S, int64_t nv,
                                      int64_t * component_map)
{
  /* Initialize each vertex with its own component label in parallel */
  OMP ("omp parallel for")
    for (uint64_t i = 0; i < nv; i++) {
      component_map[i] = i;
    }

  /* Iterate until no changes occur */
  while (1) {
    int changed = 0;

    /* For all edges in the STINGER graph of type 0 in parallel, attempt to assign
       lesser component IDs to neighbors with greater component IDs */
    STINGER_PARALLEL_FORALL_EDGES_BEGIN (S, 0) {
      if (component_map[STINGER_EDGE_DEST] <
          component_map[STINGER_EDGE_SOURCE]) {
        component_map[STINGER_EDGE_SOURCE] = component_map[STINGER_EDGE_DEST];
        changed++;
      }
    }
    STINGER_PARALLEL_FORALL_EDGES_END ();

    /* if nothing changed */
    if (!changed)
      break;

    /* Tree climbing with OpenMP parallel for */
    OMP ("omp parallel for")
      MTA ("mta assert nodep")
      for (uint64_t i = 0; i < nv; i++) {
        while (component_map[i] != component_map[component_map[i]])
          component_map[i] = component_map[component_map[i]];
      }
  }

  /* Count components */
  uint64_t components = 1;
  MTA ("mta assert nodep")
    OMP ("omp parallel for reduction(+:components)")
    for (uint64_t i = 1; i < nv; i++) {
      if (component_map[i] == i) {
	components += (stinger_vtype_get(S,i) != 0);
      }
    }

  return components;
}

static_components_workpace_t *
static_components_workspace_from_void(stinger_t * S, void ** workspace) {
  static_components_workpace_t * ws = *((static_components_workpace_t **)workspace);
  if(!ws) {
    ws = xmalloc(sizeof(static_components_workpace_t) + sizeof(int64_t) *  
      stinger_vertices_max_vertices_get(stinger_vertices_get(S)));
    if(!ws) return NULL;
    *workspace = ws;
    ws->print_count = 1;
    ws->histogram_sizes_to_file = 1;
    ws->path = "./";
    ws->filename = "components";
  }
  return ws;
}

stinger_return_t
static_components_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  static_components_workpace_t * ws = static_components_workspace_from_void(S, workspace);

  int64_t comp_count = parallel_shiloach_vishkin_components(S, stinger_vertices_max_vertices_get(stinger_vertices_get(S)), ws->components);

  if(ws->print_count) {
    printf("Static Components: %ld after %d\n", comp_count, 0);
  }
  if(ws->histogram_sizes_to_file) {
    histogram_label_counts(S,ws->components,stinger_vertices_max_vertices_get(stinger_vertices_get(S)),ws->path,ws->filename,0);
  }
  return STINGER_SUCCESS;
}

void
static_components_help(void ** workspace) {
  printf(
"Algorithm: Static Connected Components\n"
"======================================\n\n"
"Calculates the connected components in the graph using a parallel static \n"
"implementation of the Shiloach-Vishkin connected components algorithms.\n"
"The components are recalculated from scratch each iteration.\n\n"
"Options:\n"
"\tprint_count: 1 or 0 will turn on or off printing the number of \n"
"\t             components found in each iteration of the algorithm\n"
"\thistogram_sizes_to_file: 1 or 0 will turn on or off histogramming\n"
"\t             the component sizes and writing out to file\n"
"\tpath:         path to store results files\n"
"\tfilename:     prefix for the file\n");
  exit(-1);
}

stinger_return_t
static_components_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  static_components_workpace_t * ws = static_components_workspace_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  CODE_TOO=1 ./tools/stringset/stringset \
  help "static_components_help(workspace);"\
  ? "static_components_help(workspace);"\
  -? "static_components_help(workspace);"\
  --? "static_components_help(workspace);"\
  -h "static_components_help(workspace);"\
  --help "static_components_help(workspace);"\
  print_count "ws->print_count = atoi(values[i]);"\
  histogram_sizes_to_file "ws->histogram_sizes_to_file = atoi(values[i]);"\
  path "ws->path = values[i];"\
  filename "ws->filename = values[i];"\

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
			static_components_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  static_components_help(workspace);
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
		  static_components_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  static_components_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    static_components_help(workspace);
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
		    static_components_help(workspace);
		  }
		}
	      } break;
	    case 'i':
	      {
		str++; len--;
		if(!strncmp(str, "stogram_sizes_to_file", len)) {
		  str += 21; len -= 21;
		  if(len == 0) {
		    /* histogram_sizes_to_file */
		    ws->histogram_sizes_to_file = atoi(values[i]);
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
		if(!strncmp(str, "int_count", len)) {
		  str += 9; len -= 9;
		  if(len == 0) {
		    /* print_count */
		    ws->print_count = atoi(values[i]);
		  }
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

stinger_return_t
static_components_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count) {
  batch += 1;

  static_components_workpace_t * ws = *((static_components_workpace_t **)workspace);
  int64_t comp_count = parallel_shiloach_vishkin_components(S, stinger_vertices_max_vertices_get(stinger_vertices_get(S)), ws->components);

  if(ws->print_count) {
    printf("Static Components: %ld after %ld\n", comp_count, batch);
  }
  if(ws->histogram_sizes_to_file) {
    histogram_label_counts(S,ws->components,stinger_vertices_max_vertices_get(stinger_vertices_get(S)),ws->path,ws->filename,batch);
  }
  return STINGER_SUCCESS;
}

