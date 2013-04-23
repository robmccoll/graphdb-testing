#ifndef  STATIC_BETWEENNESS_CENTRALITY_H
#define  STATIC_BETWEENNESS_CENTRALITY_H

#include "stinger.h"
#include "stinger-return.h"
#include "stinger-workflow.h"

#include <stdint.h>

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * IMPLEMENTATION DEFS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

typedef double bc_t;

typedef struct {
    int64_t NV;
    int64_t* level;
    int64_t* pathsToRoot;
    bc_t* delta;
} bcTree;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * INTERFACE DEFS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

typedef struct static_betweenness_centrality_workpace {
  uint64_t nv;
  uint8_t print;
  uint8_t histogram_scores_to_file;
  uint8_t reselect_roots;
  char *  path;
  char *  filename;
  double * scores;

  bcTree**		    parallelForest;
  double**    		    totalBCSS;
  double*     		    finalBC;
  stinger_named_result_t *  finalBC_nr;
  char *		    finalBC_name;
  uint64_t**  		    parallelSingleQueue;
  int64_t     		    rootsPerThread;
  uint64_t*** 		    ppArray;
  uint64_t    		    threadCount;
  uint64_t *  		    selectedRoots;
  uint64_t    		    BC_ROOTS;
} static_betweenness_centrality_workpace_t;

stinger_return_t
static_betweenness_centrality_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace);

stinger_return_t
static_betweenness_centrality_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings);

stinger_return_t
static_betweenness_centrality_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, 
  edge_action_t * actions, int * result, int64_t count);



#endif  /*STATIC_BETWEENNESS_CENTRALITY_H*/
