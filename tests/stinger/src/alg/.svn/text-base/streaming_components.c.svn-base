#include "streaming_components.h"
#include "histogram.h"
#include "xmalloc.h"
#include "stinger-atomics.h"
#include "x86-full-empty.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * DEFINITIONS FOR INTERNAL USE
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#define LEVEL_IS_NEG(k) ((level[(k)] < 0) && (level[(k)] != INFINITY_MY))
#define LEVEL_IS_POS(k) ((level[(k)] >= 0) && (level[(k)] != INFINITY_MY))
#define LEVEL_IS_INF(k) (level[(k)] == INFINITY_MY)
#define LEVEL_EQUALS(k,y) ((level[(k)] == (y)) && (level[(k)] != INFINITY_MY))
#define SWAP_UINT64(x,y) {uint64_t tmp = (x); (x) = (y); (y) = tmp;}

#define INFINITY_MY 1073741824
#define EMPTY_NEIGHBOR -1073741824

inline void 
update_tree_for_delete_directed (
    int64_t * parentArray, 
    int64_t * parentCounter, 
    int64_t * level,
    int64_t parentsPerVertex, 
    int64_t i, 
    int64_t j);

inline int64_t
update_tree_for_insert_directed (
    int64_t * parentArray, 
    int64_t * parentCounter, 
    int64_t * level,
    int64_t * component,
    int64_t parentsPerVertex, 
    int64_t i, 
    int64_t j);

inline int64_t
is_delete_unsafe (
    int64_t * parentArray, 
    int64_t * parentCounter, 
    int64_t * level,
    int64_t parentsPerVertex, 
    int64_t i); 

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * STREAMING BFS FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

uint64_t bfs_build_component (
    struct stinger* S, 
    uint64_t currRoot, 
    uint64_t* queue,
    uint64_t* level, 
    uint64_t* parentArray, 
    uint64_t parentsPerVertex, 
    uint64_t* parentCounter, 
    int64_t * component)
{
  component[currRoot] = currRoot;
  level[currRoot] = 0;

  queue[0] = currRoot;
  int64_t qStart  = 0; 
  int64_t qEnd	  = 1;

  /* while queue is not empty */
  while(qStart != qEnd) {
    uint64_t old_qEnd = qEnd;


    OMP("omp parallel for")
    for(int64_t i = qStart; i < old_qEnd; i++) {
      uint64_t currElement = queue[i];
      uint64_t myLevel = level[currElement];
      uint64_t nextLevel = myLevel+1;

      STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, currElement) {
	uint64_t k = STINGER_EDGE_DEST;

	/* if k hasn't been found */
	if(LEVEL_IS_INF(k)) {
	  /* add k to the frontier */
	  if(INFINITY_MY == stinger_int64_cas(level + k, INFINITY_MY, nextLevel)) {
	    uint64_t which = stinger_int64_fetch_add(&qEnd, 1);
	    queue[which] = k;
	    component[k] = currRoot;
	  }
	}

	/* if k has space */
	if(parentCounter[k] < parentsPerVertex) {
	  if(LEVEL_EQUALS(k,nextLevel)) {
	    uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	    if(which < parentsPerVertex) {
	      /* add me to k's parents */
	      parentArray[k*parentsPerVertex+which] = currElement;
	    } else {
	      parentCounter[k] = parentsPerVertex;
	    }
	  } else if(LEVEL_EQUALS(k, myLevel)) {
	    uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	    if(which < parentsPerVertex) {
	    /* add me to k as a neighbor (bitwise negate for vtx 0) */
	      parentArray[k*parentsPerVertex+which] = ~currElement;
	    } else {
	      parentCounter[k] = parentsPerVertex;
	    }
	  }
	}
      } STINGER_FORALL_EDGES_OF_VTX_END();
    }
    qStart = old_qEnd;
  }

  return qEnd;
}

uint64_t bfs_build_new_component (
    struct stinger* S, 
    uint64_t currRoot, 
    uint64_t component_id, 
    int64_t start_level, 
    uint64_t* queue,
    uint64_t* level, 
    uint64_t* parentArray, 
    uint64_t parentsPerVertex, 
    uint64_t* parentCounter, 
    int64_t * component)
{
  component[currRoot] = component_id;
  level[currRoot] = start_level;

  queue[0] = currRoot;
  int64_t qStart  = 0; 
  int64_t qEnd	  = 1;


  /* while queue is not empty */
  while(qStart != qEnd) {
    uint64_t old_qEnd = qEnd;


    OMP("omp parallel for")
    for(int64_t i = qStart; i < old_qEnd; i++) {
      uint64_t currElement = queue[i];
      int64_t myLevel = level[currElement];
      int64_t nextLevel = myLevel+1;

      STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, currElement) {
	uint64_t k = STINGER_EDGE_DEST;

	/* if k hasn't been found */
	if(component[k] != component_id) {
	  /* local level */
	  int64_t k_level = readfe(level + k);
	  if(component[k] != component_id) {
	    uint64_t which = stinger_int64_fetch_add(&qEnd, 1);
	    queue[which] = k;
	    component[k] = component_id;
	    parentCounter[k] = 0;
	    k_level = nextLevel;
	  }
	  writeef(level + k, k_level);
	}

	/* if k has space */
	if(parentCounter[k] < parentsPerVertex) {
	  if(readff(level + k) == nextLevel) {
	    uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	    if(which < parentsPerVertex) {
	      /* add me to k's parents */
	      parentArray[k*parentsPerVertex+which] = currElement;
	    } else {
	      parentCounter[k] = parentsPerVertex;
	    }
	  } else if(readff(level+k) == myLevel) {
	    uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	    if(which < parentsPerVertex) {
	    /* add me to k as a neighbor (bitwise negate for vtx 0) */
	      parentArray[k*parentsPerVertex+which] = ~currElement;
	    } else {
	      parentCounter[k] = parentsPerVertex;
	    }
	  }
	}
      } STINGER_FORALL_EDGES_OF_VTX_END();
    }
    qStart = old_qEnd;
  }

  return qEnd;
}

#define LDB(...) //printf(__VA_ARGS__);

/** @brief Fix a component after an unsafe deletion. 
 * Attempts to find a path back up the tree (to a vertex that can reach a 
 * higher level). If one is found, the tree is rebuilt downward from there.
 * Otherwise, this function also builds a new component along the way.
 **/
uint64_t bfs_rebuild_component (
    struct stinger* S, 
    uint64_t currRoot, 
    uint64_t component_id, 
    int64_t start_level, 
    uint64_t* queue,
    uint64_t* level, 
    uint64_t* parentArray, 
    uint64_t parentsPerVertex, 
    uint64_t* parentCounter, 
    int64_t * component,
    int64_t * start_level_queue)
{
  int64_t old_component = component[currRoot];
  component[currRoot] = component_id;
  level[currRoot] = 0;
  parentCounter[currRoot] = 0;

  queue[0] = currRoot;
  int64_t qStart    = 0; 
  int64_t qEnd	    = 1;

  int64_t slq_start = 0;
  int64_t slq_end   = 0;

  int path_up_found = 0;

  /* while queue is not empty */
  while(qStart != qEnd) {
    uint64_t old_qEnd = qEnd;

    OMP("omp parallel for")
    for(int64_t i = qStart; i < old_qEnd; i++) {
      uint64_t currElement = queue[i];
      int64_t myLevel = level[currElement];
      int64_t nextLevel = myLevel+1;

      STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, currElement) {
	uint64_t k = STINGER_EDGE_DEST;

	/* if k hasn't been found */
	if(component[k] != component_id) {
	  /* local level */
	  int64_t k_level = readfe(level + k);
	  if(component[k] != component_id) {
	    if(k_level == start_level) {
	      LDB("\n\tFOUND ROOT: %ld", k);
	      component[k] = component_id;
	      path_up_found = 1;
	      uint64_t which = stinger_int64_fetch_add(&slq_end, 1);
	      start_level_queue[which] = k;
	    } else {
	      if((k_level < start_level) && (k_level > ~start_level)) {
		LDB("\n\tFOUND ROOT: %ld", k);
		component[k] = component_id;
		path_up_found = 1;
		uint64_t which = stinger_int64_fetch_add(&slq_end, 1);
		start_level_queue[which] = k;
	      } else {
		LDB("\n\tFOUND OTHER: %ld", k);
		parentCounter[k] = 0;
		uint64_t which = stinger_int64_fetch_add(&qEnd, 1);
		queue[which] = k;
		component[k] = component_id;
		k_level = nextLevel;
	      }
	    }
	  }
	  writeef(level + k, k_level);
	}

	/* if k has space */
	if(parentCounter[k] < parentsPerVertex) {
	  if(readff(level + k) == nextLevel) {
	    uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	    if(which < parentsPerVertex) {
	      /* add me to k's parents */
	      parentArray[k*parentsPerVertex+which] = currElement;
	    } else {
	      parentCounter[k] = parentsPerVertex;
	    }
	  } else if(readff(level+k) == myLevel) {
	    uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	    if(which < parentsPerVertex) {
	    /* add me to k as a neighbor (bitwise negate for vtx 0) */
	      parentArray[k*parentsPerVertex+which] = ~currElement;
	    } else {
	      parentCounter[k] = parentsPerVertex;
	    }
	  }
	}
      } STINGER_FORALL_EDGES_OF_VTX_END();
    }
    qStart = old_qEnd;
  }

  if(!path_up_found) {
    return qEnd;
  } else {

    OMP("omp parallel for")
    for(int64_t i = slq_start; i < slq_end; i++) {
      int64_t currElement = component[start_level_queue[i]] = old_component;
    }

    while(slq_end != slq_start) {
      /* rebuilt from the current level down */
      uint64_t old_slq_end= slq_end;

      OMP("omp parallel for")
      for(int64_t i = slq_start; i < old_slq_end; i++) {
	int64_t currElement = start_level_queue[i];
	int64_t myLevel = level[currElement];

	if(myLevel < 0)
	  myLevel = ~myLevel;

	int64_t nextLevel = myLevel+1;

	LDB("\n\tSTARTING FROM: %ld", currElement);

	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, currElement) {
	  uint64_t k = STINGER_EDGE_DEST;

	  /* if k hasn't been found */
	  if(component[k] != old_component) {
	    /* local level */
	    int64_t k_level = readfe(level + k);
	    if(component[k] != old_component) {
	      parentCounter[k] = 0;
	      uint64_t which = stinger_int64_fetch_add(&slq_end, 1);
	      start_level_queue[which] = k;
	      component[k] = old_component;
	      k_level = nextLevel;
	    }
	    writeef(level + k, k_level);
	  }

	  /* if k has space */
	  if(parentCounter[k] < parentsPerVertex) {
	    if(readff(level + k) == nextLevel) {
	      uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	      if(which < parentsPerVertex) {
		/* add me to k's parents */
		parentArray[k*parentsPerVertex+which] = currElement;
	      } else {
		parentCounter[k] = parentsPerVertex;
	      }
	    } else if(readff(level+k) == myLevel) {
	      uint64_t which = stinger_int64_fetch_add(parentCounter + k, 1);
	      if(which < parentsPerVertex) {
	      /* add me to k as a neighbor (bitwise negate for vtx 0) */
		parentArray[k*parentsPerVertex+which] = ~currElement;
	      } else {
		parentCounter[k] = parentsPerVertex;
	      }
	    }
	  }
	} STINGER_FORALL_EDGES_OF_VTX_END();
      }
      slq_start = old_slq_end;
    }
  }
}

void 
update_tree_for_delete_directed (
    int64_t * parentArray, 
    int64_t * parentCounter, 
    int64_t * level,
    int64_t parentsPerVertex, 
    int64_t i, 
    int64_t j) 
{
  int64_t local_parent_counter = readfe((uint64_t *)(parentCounter + i));
    int i_parents = 0;
    for(int64_t p = 0; p < local_parent_counter; p++) {
      /* if j is a neighbor or parent of i */
      if(parentArray[(i)*parentsPerVertex+p]==(j) ||
	  parentArray[(i)*parentsPerVertex+p]==(~j)) {
	/* replace with last element */
	local_parent_counter--;
	parentArray[(i)*parentsPerVertex + p] = parentArray[(i)*parentsPerVertex + local_parent_counter];
	parentArray[(i)*parentsPerVertex + local_parent_counter] = EMPTY_NEIGHBOR;
      }
      /* also search to see if we still have parents */
      if(parentArray[(i)*parentsPerVertex+p] >= 0) {
	i_parents = 1;
      }
    }
    if(!i_parents && LEVEL_IS_POS(i)) {
      level[i] = ~level[i];
    }
  writeef((uint64_t *)(parentCounter + i), local_parent_counter);
}

int64_t 
update_tree_for_insert_directed (
    int64_t * parentArray, 
    int64_t * parentCounter, 
    int64_t * level,
    int64_t * component,
    int64_t parentsPerVertex, 
    int64_t i, 
    int64_t j) 
{
  if(component[i] == component[j]) {
    int64_t local_parent_counter = readfe((uint64_t *)(parentCounter + j));
      if(LEVEL_IS_NEG(j)) {
	if(LEVEL_IS_POS(i) && level[i] < ~level[j]) {
	  if(local_parent_counter < parentsPerVertex) {
	    parentArray[j*parentsPerVertex + local_parent_counter] = i;
	    local_parent_counter++;
	  } else {
	    parentArray[j*parentsPerVertex] = i;
	  }
	  level[j] = ~level[j];
	} 
      } else if(LEVEL_IS_POS(i)) {
	if(local_parent_counter < parentsPerVertex) {
	  if(level[i] < level[j]) {
	    parentArray[j*parentsPerVertex + local_parent_counter] = i;
	    local_parent_counter++;
	  } else if(level[i] == level[j]) {
	    parentArray[j*parentsPerVertex + local_parent_counter] = ~i;
	    local_parent_counter++;
	  }
	} else if(level[i] < level[j]) {
	  /* search backward - more likely to have neighbors at end */
	  for(int64_t p = parentsPerVertex - 1; p >= 0; p--) {
	    if(parentArray[j*parentsPerVertex + p] < 0) {
	      parentArray[j*parentsPerVertex + p] = i;
	      break;
	    }
	  }
	}
      }
    writeef((uint64_t *)(parentCounter + j), local_parent_counter);
    return 0;
  } else {
    return 1;
  }
}

int64_t
is_delete_unsafe (
    int64_t * parentArray, 
    int64_t * parentCounter, 
    int64_t * level,
    int64_t parentsPerVertex, 
    int64_t i)
{
  int i_parents = 0;
  int i_neighbors = 0;
  for(int p = 0; p < parentCounter[i]; p++) {
    if(parentArray[(i)*parentsPerVertex+p] >= 0) {
      i_parents = 1;
    } else if(level[~parentArray[(i)*parentsPerVertex+p]] >= 0) {
      i_neighbors= 1;
    }
  }
  if(!i_parents) {
    /* if no parents, I isn't safe */
    if(level[i] >= 0) {
      fprintf(stderr,"\n\t%s %d You shouldn't see this run.", __func__, __LINE__);
      level[i] = ~level[i];
    }
    /* if no safe neighbors this delete is dangerous */
    if(!i_neighbors) {
      return 1;
    }
  }
  return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

streaming_components_workpace_t *
streaming_components_workspace_from_void(stinger_t * S, void ** workspace) {
  streaming_components_workpace_t * ws = *((streaming_components_workpace_t **)workspace);
  if(!ws) {
    ws = xmalloc(sizeof(streaming_components_workpace_t));

    if(!ws) return NULL;

    *workspace = ws;
    ws->nv = stinger_vertices_max_vertices_get(stinger_vertices_get(S));
    ws->print_count = 1;
    ws->histogram_sizes_to_file = 1;
    ws->path = "./";
    ws->filename = "streaming_components";

    ws->queue = xmalloc(4*ws->nv*sizeof(uint64_t));
    ws->level = &(ws->queue[ws->nv]);
    ws->found = &(ws->queue[2*ws->nv]);
    ws->same_level_queue = &(ws->queue[3*ws->nv]);

    ws->parentsPerVertex = 4;
    ws->parentArray	= xmalloc((ws->parentsPerVertex + 3) * ws->nv * sizeof(uint64_t));
    ws->parentCounter	= &(ws->parentArray[(ws->parentsPerVertex) * ws->nv]);
    ws->bfs_components	= &(ws->parentArray[(ws->parentsPerVertex + 1) * ws->nv]);
    ws->bfs_component_sizes = &(ws->parentArray[(ws->parentsPerVertex + 2) * ws->nv]);

    OMP("omp parallel for")
    for(int i = 0; i < ws->nv; i++) {
      ws->bfs_components[i] = i;
      ws->bfs_component_sizes[i] = 0;
      ws->level[i] = INFINITY_MY;
      ws->parentCounter[i] = 0;
    }
  }
  return ws;
}

stinger_return_t
streaming_components_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  streaming_components_workpace_t * ws = streaming_components_workspace_from_void(S, workspace);

  OMP("omp parallel for")
  for(int i = 0; i < ws->nv; i++) {
    ws->bfs_components[i] = i;
    ws->bfs_component_sizes[i] = 0;
    ws->level[i] = INFINITY_MY;
    ws->parentCounter[i] = 0;
  }

  uint64_t bfs_num_components = 0;
  for(uint64_t i = 0; i < ws->nv; i++) {
    if(ws->level[i] == INFINITY_MY) {
      if(stinger_outdegree(S, i)) {
	ws->bfs_component_sizes[i] = 
	    bfs_build_component (
		S, i, ws->queue,ws->level, ws->parentArray, 
		ws->parentsPerVertex, ws->parentCounter, ws->bfs_components);
	bfs_num_components += stinger_vtype_get(S,i) != 0;
      } 
    }
  }

  OMP("omp parallel for reduction(+:bfs_num_components)")
  for(uint64_t i = 0; i < ws->nv; i++) {
    if(ws->level[i] == INFINITY_MY) {
      ws->bfs_components[i] = i;
      ws->bfs_component_sizes[i] = 1;
      ws->level[i] = 0;
      bfs_num_components += stinger_vtype_get(S,i) != 0;
    }
  }

  if(ws->print_count) {
    printf("Components: %ld after %d\n", bfs_num_components, 0);
  }
  if(ws->histogram_sizes_to_file) {
    histogram_label_counts(S,ws->bfs_components,stinger_vertices_max_vertices_get(stinger_vertices_get(S)),ws->path,ws->filename,0);
  }
  return STINGER_SUCCESS;
}

void
streaming_components_help(void ** workspace) {
  printf(
"Algorithm: Streaming Connected Components\n"
"======================================\n\n"
"Calculates the connected components in the graph using a parallel dynamic\n"
"implementation of the Parent-Neighbor subgraph components tracking \n"
"algorithm that avoids full recompute by using a subgraph and limited\n"
"breadth-first traversals of the local graph. Components are maintained "
"as the graph changes.\n\n"
"Options:\n"
"\tprint_count: 1 or 0 will turn on or off printing the number of \n"
"\t             components found in each iteration of the algorithm\n"
"\thistogram_sizes_to_file 1 or 0 will turn on or off histogramming\n"
"\t             the component sizes and writing out to file\n"
"\tpath         path to store results files\n"
"\tfilename     prefix for the file\n");
  exit(-1);
}

stinger_return_t
streaming_components_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  streaming_components_workpace_t * ws = streaming_components_workspace_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  CODE_TOO=1 ./tools/stringset/stringset \
  help "streaming_components_help(workspace);"\
  ? "streaming_components_help(workspace);"\
  -? "streaming_components_help(workspace);"\
  --? "streaming_components_help(workspace);"\
  -h "streaming_components_help(workspace);"\
  --help "streaming_components_help(workspace);"\
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
			streaming_components_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  streaming_components_help(workspace);
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
		  streaming_components_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  streaming_components_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    streaming_components_help(workspace);
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
		    streaming_components_help(workspace);
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
streaming_components_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count) {
  batch += 1;
  streaming_components_workpace_t * ws = *((streaming_components_workpace_t **)workspace);

  int64_t * action_stack = malloc(sizeof(int64_t) * count * 2 * 2);
  int64_t * action_stack_components = malloc(sizeof(int64_t) * count * 2 * 2);
  int64_t delete_stack_top = 0;
  int64_t insert_stack_top = count * 2 * 2 - 2;

  for(uint64_t a = 0; a < count; a++) {
    if(actions[a].source >= 0) {
      /* if a new edge in stinger, update parents */
      if(result[a] & 0x1 == 0x1) {
	if(update_tree_for_insert_directed(ws->parentArray, ws->parentCounter, ws->level, ws->bfs_components, 
	  ws->parentsPerVertex, actions[a].source, actions[a].dest)) {
	  /* if component ids are not the same       *
	   * tree update not handled, push onto queue*/
	  int64_t which = stinger_int64_fetch_add(&insert_stack_top, -2);
	  action_stack[which] = actions[a].source;
	  action_stack[which+1] = actions[a].dest;
	}
      }
      if(result[a] & 0x2 == 0x2) {
	if(update_tree_for_insert_directed(ws->parentArray, ws->parentCounter, ws->level, ws->bfs_components, 
	  ws->parentsPerVertex, actions[a].dest, actions[a].source)) {
	  /* if component ids are not the same       *
	   * tree update not handled, push onto queue*/
	  int64_t which = stinger_int64_fetch_add(&insert_stack_top, -2);
	  action_stack[which] = actions[a].dest;
	  action_stack[which+1] = actions[a].source;
	}
      }
    }
  }

  for(int64_t k = count * 2 * 2 - 2; k > insert_stack_top; k -= 2) {
    int64_t i = action_stack[k]; 
    int64_t j = action_stack[k+1];

    int64_t Ci = ws->bfs_components[i];
    int64_t Cj = ws->bfs_components[j];

    if(Ci == Cj)
      continue;


    int64_t Ci_size = ws->bfs_component_sizes[Ci];
    int64_t Cj_size = ws->bfs_component_sizes[Cj];

    if(Ci_size > Cj_size) {
      SWAP_UINT64(i,j)
      SWAP_UINT64(Ci, Cj)
      SWAP_UINT64(Ci_size, Cj_size)
    }

    ws->parentArray[i*ws->parentsPerVertex] = j;
    ws->parentCounter[i] = 1;
    /* handle singleton */
    if(Ci_size == 1) {
      ws->bfs_component_sizes[Ci] = 0;
      ws->bfs_component_sizes[Cj]++;
      ws->bfs_components[i] = Cj;
    } else {
      ws->bfs_component_sizes[Cj] += bfs_build_new_component(S, i, Cj, (ws->level[j] >= 0 ? ws->level[j] : ~ws->level[j])+1, ws->queue, ws->level, ws->parentArray, ws->parentsPerVertex, ws->parentCounter, ws->bfs_components);
      ws->bfs_component_sizes[Ci] = 0;
    }
  }

  for(uint64_t a = 0; a < count; a++) {
    if(actions[a].source < 0) {
      /* if a new edge in stinger, update parents */
      if(result[a] & 0x1 == 0x1) {
	uint64_t which = stinger_int64_fetch_add(&delete_stack_top, 2);
	action_stack[which] = -actions[a].source;
	action_stack[which+1] = -actions[a].dest;
	action_stack_components[which] = ws->bfs_components[-actions[a].source];
	action_stack_components[which+1] = ws->bfs_components[-actions[a].dest];
	update_tree_for_delete_directed(ws->parentArray, ws->parentCounter, ws->level, ws->parentsPerVertex, -actions[a].source, -actions[a].dest);
      }
      if(result[a] & 0x2 == 0x2) {
	uint64_t which = stinger_int64_fetch_add(&delete_stack_top, 2);
	action_stack[which] = -actions[a].dest;
	action_stack[which+1] = -actions[a].source;
	action_stack_components[which] = ws->bfs_components[-actions[a].dest];
	action_stack_components[which+1] = ws->bfs_components[-actions[a].source];
	update_tree_for_delete_directed(ws->parentArray, ws->parentCounter, ws->level, ws->parentsPerVertex, -actions[a].dest, -actions[a].source);
      }
    }
  }

  /* explicitly not parallel for all unsafe deletes */
  for(uint64_t k = 0; k < delete_stack_top; k += 2) {
    int64_t i = action_stack[k];
    int64_t j = action_stack[k+1];
    int64_t Ci_prev = action_stack_components[k];
    int64_t Cj_prev = action_stack_components[k+1];
    if(i != -1)  {
      int64_t Ci = ws->bfs_components[i];
      int64_t Cj = ws->bfs_components[j];

      if(Ci == Cj && Ci == Ci_prev) {
	if(Ci != i && stinger_outdegree(S, i) == 0) {
	  ws->bfs_component_sizes[Ci]--;
	  ws->bfs_components[i] = i;
	  ws->bfs_component_sizes[i] = 1;
	} else if(Cj != j && stinger_outdegree(S, j) == 0) {
	  ws->bfs_component_sizes[Cj]--;
	  ws->bfs_components[j] = j;
	  ws->bfs_component_sizes[j] = 1;
	} else {
	  int64_t level_i = ws->level[i];
	  int64_t level_j = ws->level[j];

	  if(level_i < 0)
	    level_i = ~level_i;
	  if(level_j < 0)
	    level_j = ~level_j;

	  if(level_i > level_j) {
	    if(is_delete_unsafe(ws->parentArray, ws->parentCounter, ws->level, ws->parentsPerVertex, i)) {
	      ws->parentCounter[i] = 0;
	      ws->bfs_component_sizes[i] = bfs_rebuild_component(S, i, i, (ws->level[i] >= 0 ? ws->level[i] : ~ws->level[i]), ws->queue, ws->level, ws->parentArray, ws->parentsPerVertex, ws->parentCounter, ws->bfs_components, ws->same_level_queue);
	      if(ws->bfs_components[i] != i) {
		ws->bfs_component_sizes[i] = 0;
	      } else {
		ws->bfs_component_sizes[Cj] -= ws->bfs_component_sizes[i];
	      }
	    }
	  } else {
	    if(is_delete_unsafe(ws->parentArray, ws->parentCounter, ws->level, ws->parentsPerVertex, j)) {
	      ws->parentCounter[j] = 0;
	      ws->bfs_component_sizes[j] = bfs_rebuild_component(S, j, j, (ws->level[j] >= 0 ? ws->level[j] : ~ws->level[j]), ws->queue, ws->level, ws->parentArray, ws->parentsPerVertex, ws->parentCounter, ws->bfs_components, ws->same_level_queue);
	      if(ws->bfs_components[j] != j) {
		ws->bfs_component_sizes[j] = 0;
	      } else {
		ws->bfs_component_sizes[Ci] -= ws->bfs_component_sizes[j];
	      }
	    }
	  }
	}
      }

      Ci = ws->bfs_components[i];
      Cj = ws->bfs_components[j];
    }
  }

  OMP("omp parallel for")
  for(uint64_t k = 0; k < delete_stack_top; k += 2) {
    int64_t i = action_stack[k];
    int64_t j = action_stack[k+1];
    if(!(is_delete_unsafe(ws->parentArray, ws->parentCounter, ws->level, ws->parentsPerVertex, i))) {
      action_stack[k] = -1;
      action_stack[k+1] = -1;
    }
  }

  if(ws->print_count) {
    uint8_t * exist = xcalloc(sizeof(uint8_t), ws->nv);
    uint64_t comp_count = 0;
    for(uint64_t v = 0; v < ws->nv; v++) {
      if(stinger_vtype_get(S,v) != 0 && !exist[ws->bfs_components[v]]) {
	exist[ws->bfs_components[v]] = 1;
	comp_count++;
      }
    }
    free(exist);
    printf("Components: %ld after %ld\n", comp_count, batch);
  }
  if(ws->histogram_sizes_to_file) {
    histogram_label_counts(S,ws->bfs_components,stinger_vertices_max_vertices_get(stinger_vertices_get(S)),ws->path,ws->filename,batch);
  }
  return STINGER_SUCCESS;
}

