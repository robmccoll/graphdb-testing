#include "static_betweenness_centrality.h"
#include "xmalloc.h"
#include "histogram.h"

#include <omp.h>

#define INFINITY_MY INT64_MAX

/**
* @file static_betweenness_centrality.c
* @brief Computes Approximate Betweenness Centrality
* @author Rob Mccoll
* @date 2013-02-22
*
* Selecting a number of roots (num_roots), this algorithm will compute the Betweenness
* Centrality contributions of the shortest paths from those num_roots roots to all other 
* vertices in the graph (keeping in mind that if a root is not picked in a given component, 
* no approximate centralities will be calculated for that component resulting in zeros for 
* all vertices in that component).  Also, keep in mind that normally BC will range from 
* [0, (n-1)(n-2)] where n is the number of vertices.  In this approximation, values will 
* range from [0, k(n-2)].  There is no normalization.
*
* The value of this approximation is intended to be found in the relative ranking.  The
* accuracy of the ranking should be very high for the top scoring vertices and low
* for low scoring vertices.  See " Approximating Betweenness Centrality" by Bader et al.
*/


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * IMPLEMENTATION FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void REDUCTION(double_t* finalResultArray, double_t** parallelArray, int64_t threadCount)
{
  OMP("omp parallel for")
  for(int64_t v = 0; v < STINGER_MAX_LVERTICES; v++) {
    for(int64_t t = 0; t < threadCount; t++)
      finalResultArray[v] += parallelArray[t][v];
  }
}

/**
* @brief Computes betweenness centrality based on Brandes's algorithm.
*
* @param tree The tree that needs to be unallocated.
* @param someGraph The tree that needs to be unallocated.
* @param currRoot The vertex in the graph that will be used as the root of the BFS traversal.
* @param totalBC The BC score for each vertex (array) .
* @param someQueue The queue that will be used by the algorithm. It is possible to use a regular queue or a multi-level queue.
* @param parentParam The data structure that will be used to maintain the parents that are in the BFS traversal.
* @param someDataStructure
*
* @return None
*/

void bfsBrandesPerTreeSTINGERSINGLEPA(bcTree* tree, stinger_t * sStinger, uint32_t currRoot, double* totalBC,
                                      uint64_t * Queue, uint64_t ** parentArray, int64_t * parentCounter)
{
  for(uint64_t j = 0; j < tree->NV; j++) {
    parentCounter[j] = 0;
    tree->level[j] = INFINITY_MY;
    tree->pathsToRoot[j] = INFINITY_MY;
    tree->delta[j] = 0;
  }

  tree->level[currRoot] = 0;
  tree->pathsToRoot[currRoot] = 1;
  Queue[0] = currRoot;
  int64_t qStart = 0, qEnd = 1;

  // While queue is not empty
  while(qStart != qEnd) {
    uint64_t currElement;
    currElement = Queue[qStart];
    qStart++;
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
      uint64_t k = STINGER_EDGE_DEST;

      // If this is a neighbor and has not been found
      if(tree->level[k] > tree->level[currElement]) {
        // Checking if "k" has been found.
        if(tree->level[k] == INFINITY_MY) {
          tree->level[k] = tree->level[currElement] + 1;
          tree->delta[k] = 0;
          Queue[qEnd++] = k;
        }

        if(tree->pathsToRoot[k] == INFINITY_MY) {
          // k has not been found and therefore its paths to the roots are through its parent.
          tree->pathsToRoot[k] = tree->pathsToRoot[currElement];
          parentArray[k][parentCounter[k]++] = currElement;
        } else {
          // k has been found and has multiple paths to the root as it has multiple parents.
          tree->pathsToRoot[k] += tree->pathsToRoot[currElement];
          parentArray[k][parentCounter[k]++] = currElement;
        }
      }
    }
    STINGER_FORALL_EDGES_OF_VTX_END();
  }

  // Using Brandes algorithm to compute BC for a specific tree.
  // Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
  // up the tree, all the way to the root.
  qEnd = qStart - 1;

  while(qEnd >= 0) {
    uint64_t currElement = Queue[qEnd];

    for(int j = 0; j < parentCounter[currElement]; j++) {
      uint64_t k = parentArray[currElement][j];
      tree->delta[k] +=
        ((bc_t)tree->pathsToRoot[k] / (bc_t)tree->pathsToRoot[currElement]) *
        (bc_t)(tree->delta[currElement] + 1);
    }

    if(currElement != currRoot) {
      totalBC[currElement] += tree->delta[currElement];
    }

    qEnd--;
  }

  return;
}


void bfsBrandesSTINGERSINGLEPA(bcTree** treeArray, stinger_t * S, double** totalBCArray,
                               uint64_t ** Queues, uint64_t*** parentArrayParam,
                               double* timePerThread, int64_t* selectedRoots, int64_t rootsPerThread)
{
  OMP("omp parallel")
  {
#if defined(_OPENMP)
    int thread = omp_get_thread_num();
#else
    int thread = 0;
#endif

    int start, stop;
    start = (thread) * rootsPerThread;
    stop     = (thread + 1) * rootsPerThread;

    int64_t * parentCounter = (int64_t*)malloc(sizeof(uint64_t) * (*treeArray)->NV);

    double* totalBC = totalBCArray[thread];
    bcTree* tree = treeArray[thread];
    uint64_t** parentArray = NULL;

    if(parentArrayParam != NULL)
      parentArray = parentArrayParam[thread];


    for(int i = 0; i < tree->NV; i++) {
      totalBC[i] = 0.0;
    }

    for(uint64_t i = start; i < stop; i++) {
      bfsBrandesPerTreeSTINGERSINGLEPA(tree, S, selectedRoots[i], totalBC, Queues[thread], parentArray, parentCounter);
    }


    free(parentCounter);
  }
}

// Destroys the parent array
void destroyParentArray(uint64_t** parentArray,uint64_t NV)
{

  for(uint64_t v=0; v<NV;v++) {
    free(parentArray[v]);
  }

  free(parentArray);
}

// Destroys the parent array of each thread/core.
void destroyParallelParentArray(uint64_t*** parallelParentArray,uint64_t NV,uint64_t threadCount)
{

  for(uint64_t t=0; t<threadCount;t++) {
    destroyParentArray(parallelParentArray[t],NV);
  }

  free(parallelParentArray);
}

// Creates the parent array needed for the computation of BC.
// The output of this function is an array of arrays. Each vertex is allocated an array based on it's adjaency.
uint64_t** createParentArrayStinger(struct stinger* S, uint64_t NV)
{
  uint64_t** parentArray = (uint64_t**)malloc(NV * sizeof(uint64_t*));

  for(uint64_t v = 0; v < NV; v++) {
    parentArray[v] = (uint64_t*)malloc(stinger_outdegree(S, v) * sizeof(uint64_t));
  }

  return parentArray;
}

uint64_t*** createParallelParentArrayStinger(struct stinger* S, uint64_t NV, uint64_t threadCount)
{
  uint64_t*** parallelParentArray = (uint64_t***)malloc(threadCount * sizeof(uint64_t**));

  for(uint64_t t = 0; t < threadCount; t++) {
    parallelParentArray[t] = createParentArrayStinger(S, NV);

    if(parallelParentArray[t] == NULL)
      printf("Failed to allocated memory for parallel parent array\n");
  }

  return parallelParentArray;
}

double** createParallelBetweennessArray(int64_t threadCount, int64_t NV)
{
  double** totalBC = (double**)malloc((threadCount) * sizeof(double*));

  for(int64_t i = 0; i < threadCount; i++) {
    totalBC[i] = malloc(sizeof(double) * NV);
  }

  return totalBC;
}

/**
 * @brief Creates the data structures needed for computing BC.
 *       These include the level of each vertex in the BFS tree.
*       The number of shortest paths each vertex has to the root.
*       The delta value computed in the dependency accumulation stage/
*
* @param numVertices The number of vertices in the graph
*
* @return Returns the created data structure.
*/
bcTree* CreateTree(int64_t numVertices)
{
  bcTree* newTree;
  newTree = (bcTree*)xmalloc(sizeof(bcTree));
  newTree->NV = numVertices;
  newTree->level = (int64_t*)xcalloc(numVertices, sizeof(int64_t));
  newTree->pathsToRoot = (int64_t*)xcalloc(numVertices, sizeof(int64_t));
  newTree->delta = (bc_t*)xcalloc(numVertices, sizeof(bc_t));
  return newTree;
}

bcTree** createParallelForest(int64_t threadCount, int64_t NV)
{
  bcTree** parallelForest = (bcTree**)xcalloc(threadCount, sizeof(bcTree*));

  for(int64_t i = 0; i < threadCount; i++) {
    parallelForest[i] = CreateTree(NV);
  }

  return parallelForest;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * UTIL FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static_betweenness_centrality_workpace_t *
static_betweenness_centrality_workspace_from_void(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace)
{
  static_betweenness_centrality_workpace_t * ws =
    *((static_betweenness_centrality_workpace_t**)workspace);

  if(!ws) {
    ws = xmalloc(sizeof(static_betweenness_centrality_workpace_t));

    if(!ws) return NULL;

    *workspace = ws;
    ws->nv = stinger_vertices_max_vertices_get(stinger_vertices_get(S));
    ws->print = 1;
    ws->histogram_scores_to_file = 1;
    ws->reselect_roots = 0;
    ws->path = "./";
    ws->filename = "static_bc";
    ws->BC_ROOTS = 256;
    ws->threadCount = 1;

#if defined(_OPENMP)
    OMP("omp parallel") {
      OMP("omp master")
      ws->threadCount = omp_get_num_threads();
    }
#endif

    ws->selectedRoots = NULL;
    ws->parallelForest = NULL;
    ws->totalBCSS = NULL;
    ws->finalBC = NULL;
    ws->parallelSingleQueue = NULL;
    ws->ppArray = NULL;

    ws->parallelForest  = createParallelForest(ws->threadCount, ws->nv);
    ws->totalBCSS	= createParallelBetweennessArray(ws->threadCount, ws->nv);
    ws->finalBC_name	= "betweenness_centrality";
    ws->finalBC_nr	= stinger_workflow_new_named_result(wkflow, ws->finalBC_name, NR_DBL, ws->nv);

    ws->parallelSingleQueue = (uint64_t**)malloc(sizeof(uint64_t*) * ws->threadCount);
    for(int i = 0; i < ws->threadCount; i++)
      ws->parallelSingleQueue[i] = (uint64_t*)malloc(sizeof(uint64_t) * ws->nv);

  }

  return ws;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
stinger_return_t
static_betweenness_centrality_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace)
{
  static_betweenness_centrality_workpace_t * ws =
    static_betweenness_centrality_workspace_from_void(S, wkflow, workspace);

  ws->finalBC = stinger_named_result_write_data(ws->finalBC_nr);

  ws->rootsPerThread = ws->BC_ROOTS / ws->threadCount;
  ws->ppArray = createParallelParentArrayStinger(S , ws->nv, ws->threadCount);

  if(!ws->selectedRoots) {
    ws->selectedRoots = (int64_t*)malloc(sizeof(int64_t) * ws->BC_ROOTS);
  }

  int64_t r1 = 0;

  for(int64_t v = 0; v < STINGER_MAX_LVERTICES && r1 < ws->BC_ROOTS; v++) {
    if(stinger_outdegree(S, v)) {
      ws->selectedRoots[r1++] = v;
    }
  }

  for(int v = 0; v < STINGER_MAX_LVERTICES; v++)
    ws->finalBC[v] = 0;

  bfsBrandesSTINGERSINGLEPA(ws->parallelForest, S, ws->totalBCSS, ws->parallelSingleQueue,
                            ws->ppArray, NULL, ws->selectedRoots, ws->rootsPerThread);
  REDUCTION(ws->finalBC, ws->totalBCSS, ws->threadCount);

  if(ws->histogram_scores_to_file) {
    histogram_double(S, ws->finalBC, stinger_vertices_max_vertices_get(stinger_vertices_get(S)), ws->path,  ws->filename, 0);
  }

  if(ws->print) {
    printf("Betweenness Centrality Init\n");
  }

  destroyParallelParentArray(ws->ppArray, ws->nv, ws->threadCount);
  stinger_named_result_commit_data(ws->finalBC_nr);
}

void
static_betweenness_centrality_help(void ** workspace)
{
  printf(
    "Algorithm: Static Approximate Betweenness Centrality\n"
    "======================================\n\n"
    "Calculates the approximate betweenness centrality of all nodes in the \n"
    "graph using a parallel static implementation of the Brandes approach.\n"
    "The scores are recalculated from scratch each iteration.\n\n"
"Selecting a number of roots (num_roots), this algorithm will compute the Betweenness     \n"
"Centrality contributions of the shortest paths from those num_roots roots to all other   \n"
"vertices in the graph (keeping in mind that if a root is not picked in a given component,\n"
"no approximate centralities will be calculated for that component resulting in zeros for \n"
"all vertices in that component).  Also, keep in mind that normally BC will range from    \n"
"[0, (n-1)(n-2)] where n is the number of vertices.  In this approximation, values will   \n"
"range from [0, k(n-2)].  There is no normalization.                                      \n"
"                                                                                         \n"
"The value of this approximation is intended to be found in the relative ranking.  The    \n"
"accuracy of the ranking should be very high for the top scoring vertices and low         \n"
"for low scoring vertices.  See \"Approximating Betweenness Centrality\" by Bader et al.   \n"
    "Options:\n"
    "\tprint: 1 or 0 will turn on or off printing when the computation completes\n"
    "\thistogram_scores_to_file: 1 or 0 will turn on or off histogramming\n"
    "\t             the BC scores and writing out to file. Scores are truncated\n"
    "\t             to the nearest integer.\n"
    "\treselect_roots: 1 or 0 will force selecting a new set of roots at recompute\n"
    "\tnum_roots: the number of roots to be selected by the algorithm\n"
    "\tpath:         path to store results files\n"
    "\tfilename:     prefix for the file\n");
  exit(-1);
}

stinger_return_t
static_betweenness_centrality_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings)
{
  char ** keys = NULL, ** values = NULL;
  int num = 0;
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);
  static_betweenness_centrality_workpace_t * ws =
    static_betweenness_centrality_workspace_from_void(S, wkflow, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);

    /*
      INPUT TO STRINGSET TO GENERATE CODE BELOW:

      CODE_TOO=1 ./tools/stringset/stringset \
      help "static_betweenness_centrality_help(workspace);"\
      ? "static_betweenness_centrality_help(workspace);"\
      -? "static_betweenness_centrality_help(workspace);"\
      --? "static_betweenness_centrality_help(workspace);"\
      -h "static_betweenness_centrality_help(workspace);"\
      --help "static_betweenness_centrality_help(workspace);"\
      print "ws->print = atoi(values[i]);"\
      histogram_scores_to_file "ws->histogram_scores_to_file = atoi(values[i]);"\
      path "ws->path = values[i];"\
      filename "ws->filename = values[i];"\
      reselect_roots "ws->reselect_roots = atol(values[i]);"\
      num_roots "ws->BC_ROOTS = atol(values[i]);"\
      finalBC_name "stinger_workflow_delete_named_result(wkflow, ws->finalBC_name); \
	ws->finalBC_name = values[i]; \
	ws->finalBC_name = stinger_workflow_new_named_result(wkflow, ws->finalBC_name, NR_DBL, ws->nv); "

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
			static_betweenness_centrality_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  static_betweenness_centrality_help(workspace);
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
		  static_betweenness_centrality_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  static_betweenness_centrality_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    static_betweenness_centrality_help(workspace);
	  }
	} break;
      case 'f':
	{
	  str++; len--;
	  if(!strncmp(str, "i", 1)) {
	    str += 1; len -= 1;
	    if(len) switch(*str) {
	      case 'l':
		{
		  str++; len--;
		  if(!strncmp(str, "ename", len)) {
		    str += 5; len -= 5;
		    if(len == 0) {
		      /* filename */
		      ws->filename = values[i];
		    }
		  }
		} break;
	      case 'n':
		{
		  str++; len--;
		  if(!strncmp(str, "alBC_name", len)) {
		    str += 9; len -= 9;
		    if(len == 0) {
		      /* finalBC_name */
		      stinger_workflow_delete_named_result(wkflow, ws->finalBC_name); 	ws->finalBC_name = values[i]; 	ws->finalBC_name = stinger_workflow_new_named_result(wkflow, ws->finalBC_name, NR_DBL, ws->nv); 
		    }
		  }
		} break;
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
		    static_betweenness_centrality_help(workspace);
		  }
		}
	      } break;
	    case 'i':
	      {
		str++; len--;
		if(!strncmp(str, "stogram_scores_to_file", len)) {
		  str += 22; len -= 22;
		  if(len == 0) {
		    /* histogram_scores_to_file */
		    ws->histogram_scores_to_file = atoi(values[i]);
		  }
		}
	      } break;
	  }
	} break;
      case 'n':
	{
	  str++; len--;
	  if(!strncmp(str, "um_roots", len)) {
	    str += 8; len -= 8;
	    if(len == 0) {
	      /* num_roots */
	      ws->BC_ROOTS = atol(values[i]);
	    }
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
		if(!strncmp(str, "int", len)) {
		  str += 3; len -= 3;
		  if(len == 0) {
		    /* print */
		    ws->print = atoi(values[i]);
		  }
		}
	      } break;
	  }
	} break;
      case 'r':
	{
	  str++; len--;
	  if(!strncmp(str, "eselect_roots", len)) {
	    str += 13; len -= 13;
	    if(len == 0) {
	      /* reselect_roots */
	      ws->reselect_roots = atol(values[i]);
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
static_betweenness_centrality_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch,
    edge_action_t * actions, int * result, int64_t count)
{
  batch += 1;
  static_betweenness_centrality_workpace_t * ws =
    static_betweenness_centrality_workspace_from_void(S, wkflow, workspace);

  ws->finalBC = stinger_named_result_write_data(ws->finalBC_nr);

  ws->rootsPerThread = ws->BC_ROOTS / ws->threadCount;
  ws->ppArray = createParallelParentArrayStinger(S , ws->nv, ws->threadCount);

  if(!ws->selectedRoots) {
    ws->selectedRoots = (int64_t*)malloc(sizeof(int64_t) * ws->BC_ROOTS);
  }

  if(ws->reselect_roots) {
    int64_t r1 = 0;

    for(int64_t v = 0; v < STINGER_MAX_LVERTICES && r1 < ws->BC_ROOTS; v++) {
      if(stinger_outdegree(S, v)) {
        ws->selectedRoots[r1++] = v;
      }
    }
  }

  for(int v = 0; v < STINGER_MAX_LVERTICES; v++)
    ws->finalBC[v] = 0;

  bfsBrandesSTINGERSINGLEPA(ws->parallelForest, S, ws->totalBCSS, ws->parallelSingleQueue,
                            ws->ppArray, NULL, ws->selectedRoots, ws->rootsPerThread);
  REDUCTION(ws->finalBC, ws->totalBCSS, ws->threadCount);
  histogram_double(S, ws->finalBC, stinger_vertices_max_vertices_get(stinger_vertices_get(S)), ws->path,  ws->filename, batch);

  if(ws->print) {
    printf("Betweenness Centrality Batch %ld\n", batch);
  }

  destroyParallelParentArray(ws->ppArray, ws->nv, ws->threadCount);
  stinger_named_result_commit_data(ws->finalBC_nr);
}
