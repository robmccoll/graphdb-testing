#include "static_multicontract_clustering.h"
#include "xmalloc.h"
#include "histogram.h"
#include "stinger-return.h"
#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "xmalloc.h"
#include "timer.h"

#include <stdint.h>

double
sum_all_edgeweights(
    struct stinger * S, 
    int64_t type)
{
  double sum = 0;
  struct stinger_eb * ebpool_priv = S->ebpool->ebpool;

  /* for each edge block */
  OMP("omp parallel for reduction(+:sum)")						
  MTA("mta assert parallel")						
  for(uint64_t eb_index = 0; eb_index < S->ETA[(type)].high; eb_index++) {	
    struct stinger_eb *  cur_eb = ebpool_priv + S->ETA[(type)].blocks[eb_index]; 
    uint64_t stop = stinger_eb_high(cur_eb);
    for(uint64_t e = 0; e < stop; e++) { 
      if(!stinger_eb_is_blank(cur_eb, e)) {                   
	sum += cur_eb->edges[e].weight;
      }								
    } 								
  }									
  return sum;
}

void
sequence(uint64_t * arr, uint64_t count) {
  OMP("omp parallel for")
  for(uint64_t i = 0; i < count; i++) {
    arr[i] = i;
  }
}

void
zero(double * arr, uint64_t count) {
  OMP("omp parallel for")
  for(uint64_t i = 0; i < count; i++) {
    arr[i] = 0;
  }
}

void
score_mean_variance_first_match(
    struct stinger * S, 
    uint64_t nv, 
    double volume, 
    double * scores, 
    uint64_t * matches, 
    double * sum, 
    double * sum_squares) 
{
  double local_sum = 0;
  double local_sum_squares = 0;

  /* precompute as much of the scores as possible */
  double half_volume = volume / 2;
  double two_div_vol_sqd = 2 / ((volume) * (volume));

  /* for each vertex */
  OMP("omp parallel for reduction(+:local_sum) reduction(+:local_sum_squares)")
  for(uint64_t u = 0; u < nv; u++) {
    double best_found = 0;
    uint64_t best_match = matches[u];

    /* precompute more of the score (vtx specific part) */
    double wt_u_x2_div_vol_sqd = stinger_vweight_get(S,u) * two_div_vol_sqd;

    /* for all edge blocks of that vertex */
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, u) {
      double wt_v = stinger_vweight_get(S, STINGER_EDGE_DEST);
      double score = (((double)STINGER_EDGE_WEIGHT) / half_volume) - (wt_v * wt_u_x2_div_vol_sqd);
      /* check for the best score for this vertex*/
      if(score > best_found) {
	best_found = score;
	best_match = STINGER_EDGE_DEST;
      }
      /* sum the score and its square */
      local_sum += score;
      local_sum_squares += (score * score);
    } STINGER_FORALL_EDGES_OF_VTX_END();

    /* writeback the best score */
    scores[u] = best_found;
    matches[u] = best_match;
  }

  *sum = local_sum;
  *sum_squares = local_sum_squares;
}

void
filter_scores(
    uint64_t nv, 
    double * scores, 
    uint64_t * matches, 
    double sum, 
    double sum_squares) 
{
  /* compute the threshold */
  double mean = (sum) / nv;
  double mean_of_sq = (sum_squares) / nv;
  double variance = mean_of_sq - (mean * mean);
  double st_dev = sqrt(variance);
  double threshold = mean - 1.5 * st_dev;

  /* for all vertices */
  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    double score = scores[v];
    /* if score is below the threshold, don't match */
    if(score != 0 && score < threshold) {
      matches[v] = v;
    } else {
      /* if my match is trying to match on a lower-scoring edge, remove its match */
      uint64_t match = matches[v];
      if(scores[match] <= score) {
	matches[match] = match;
      }
    }
  }
}

void
tree_climb(
    uint64_t nv, 
    uint64_t * matches) 
{
  /* for all vertices */
  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    uint64_t older_match, old_match, match = v;
    old_match = v;
    /* climb the tree of matchings until we reach a root or cycle */
    do {
      older_match = old_match;
      old_match = match;
      match = matches[match];
      /* found a cycle - pick the lesser ID */
      if(match == older_match) {
	match = match > old_match ? old_match : match;
	break;
      }
    } while(old_match != match);
    matches[v] = match;
  }
}

int
multi_contract_root(
    struct stinger * S, 
    uint64_t nv,
    uint64_t * matches,
    int64_t timestamp)
{
  int work_remaining = 0;

  struct stinger_eb * ebpool_priv = S->ebpool->ebpool;

  /* for all vertices */
  OMP("omp parallel for reduction(+:work_remaining)")
  for(uint64_t u = 0; u < nv; u++) {
    uint64_t match_u = matches[u];
    /* if it's a root */
    if(match_u == u) {
      /* for all edge blocks of the vertex */
      struct stinger_eb * currentBlock = ebpool_priv + stinger_adjacency_get(S,u);
      while(currentBlock != ebpool_priv) {
	struct stinger_edge * edge = currentBlock->edges;
	struct stinger_edge * last_edge = edge + currentBlock->high;
	/* for each edge in a block */
	for(; edge < last_edge; ++edge) {
	  int64_t v = edge->neighbor;
	  /* if the edge exists */
	  if(v >= 0) {
	    uint64_t match_v = matches[v];
	    /* and the edge is to be contracted, remove and increment vtx weight */
	    if(match_u == match_v) {
	      work_remaining = 1;
	      stinger_vweight_increment_atomic(S, match_u, edge->weight);
	      edge->neighbor = ~v;
	      currentBlock->numEdges--;
	      stinger_indegree_increment_atomic(S, v, -1);
	      stinger_outdegree_increment_atomic(S, u, -1);
	    /* otherwise remove, remap, reinsert */
	    } else if(match_v != v) {
	      work_remaining = 1;
	      edge->neighbor = ~v;
	      currentBlock->numEdges--;
	      stinger_indegree_increment_atomic(S, v, -1);
	      stinger_outdegree_increment_atomic(S, u, -1);
	      stinger_incr_edge(S, currentBlock->etype, u, match_v, edge->weight, timestamp);
	    }
	  }
	}
	currentBlock = currentBlock->next + ebpool_priv;
      }
    }
  }

  return work_remaining;
}

int
multi_contract_tree(
    struct stinger * S, 
    uint64_t nv,
    uint64_t * matches,
    int64_t timestamp)
{
  int work_remaining = 0;

  struct stinger_eb * ebpool_priv = S->ebpool->ebpool;

  /* for all vertices */
  OMP("omp parallel for reduction(+:work_remaining)")
  for(uint64_t u = 0; u < nv; u++) {
    uint64_t match_u = matches[u];
    struct stinger_eb * currentBlock = ebpool_priv + stinger_adjacency_get(S,u);
    if(match_u != u) {
      while(currentBlock != ebpool_priv) {
	struct stinger_edge * edge = currentBlock->edges;
	struct stinger_edge * last_edge = edge + currentBlock->high;
	/* for each edge in a block */
	for(; edge < last_edge; ++edge) {
	  int64_t v = edge->neighbor;
	  if(v >= 0) {
	    work_remaining = 1;
	    uint64_t match_v = matches[v];
	    if(match_u == match_v) {
	      stinger_vweight_increment_atomic(S, match_u, edge->weight);
	    } else {
	      stinger_incr_edge(S, currentBlock->etype, match_u, match_v, edge->weight, timestamp);
	    }
	    edge->neighbor = ~v;
	    currentBlock->numEdges--;
	    stinger_indegree_increment_atomic(S, v, -1);
	    stinger_outdegree_increment_atomic(S, u, -1);
	  }
	}
	currentBlock = currentBlock->next + ebpool_priv;
      }
      stinger_vweight_increment_atomic(S, match_u, stinger_vweight_get(S,u));
      stinger_vweight_set(S, u, 0);
    }
  }

  return work_remaining;
}

void
static_multi_contract_clustering (
    uint64_t * matches,
    double * scores,
    uint64_t nv,
    struct stinger * S,
    struct stinger * S_orig)
{

  double volume = 0;

  for(uint64_t t= 0; t < STINGER_NUMETYPES; t++) {
    volume += sum_all_edgeweights(S, t);
  }

  int work_remaining = 1;

  uint64_t iteration = 0;
  double sum, sum_squares;

  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    stinger_vweight_set(S,v,0);
  }

  sequence(matches, nv);

  int64_t max_iterations = 1000;
  if(getenv("MAX_ITERATIONS") && atoi(getenv("MAX_ITERATIONS")) > 0)
     max_iterations = atoi(getenv("MAX_ITERATIONS")) - 1;

  while(work_remaining && iteration < max_iterations) {
    zero(scores, nv);
    work_remaining = 0;
    sum = 0;
    sum_squares = 0;

    score_mean_variance_first_match(S, nv, volume, scores, matches, &sum, &sum_squares);

    filter_scores(nv, scores, matches, sum, sum_squares);

    tree_climb(nv, matches); 

    work_remaining += multi_contract_root(S, nv, matches, 0);

    work_remaining += multi_contract_tree(S, nv, matches, 0);

    iteration++;
  }
}

static_clustering_workpace_t *
static_clustering_workspace_from_void(stinger_t * S, void ** workspace) {
  static_clustering_workpace_t * ws = *((static_clustering_workpace_t **)workspace);
  if(!ws) {
    ws = xmalloc(sizeof(static_clustering_workpace_t) + sizeof(int64_t) *  
      stinger_vertices_max_vertices_get(stinger_vertices_get(S)));
    if(!ws) return NULL;
    *workspace = ws;
    ws->print_count = 1;
    ws->histogram_sizes_to_file = 1;
    ws->path = "./";
    ws->filename = "clusters";
    ws->scores = xmalloc(sizeof(double) * stinger_vertices_max_vertices_get(stinger_vertices_get(S)));
  }
  return ws;
}

stinger_return_t
static_clustering_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  static_clustering_workpace_t * ws = static_clustering_workspace_from_void(S, workspace);

  stinger_t * S_cluster = stinger_new();
  for(uint64_t t = 0; t < STINGER_NUMETYPES; t++) {
    STINGER_PARALLEL_FORALL_EDGES_BEGIN(S, t) {
      stinger_insert_edge(S_cluster, STINGER_EDGE_TYPE, STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_WEIGHT, STINGER_EDGE_TIME_RECENT);
    } STINGER_PARALLEL_FORALL_EDGES_END();
  }
  static_multi_contract_clustering(ws->labels, ws->scores, stinger_vertices_max_vertices_get(stinger_vertices_get(S)), S_cluster, S);
  stinger_free_all(S_cluster);

  if(ws->print_count) {
    uint64_t comm_count = 0;
    for(uint64_t v = 0; v < stinger_vertices_max_vertices_get(stinger_vertices_get(S)); v++) {
      if(stinger_vtype_get(S, v) != 0 && ws->labels[v] == v)
	comm_count++;
    }
    printf("Clusters: %ld after %d\n", comm_count, 0);
  }
  if(ws->histogram_sizes_to_file) {
    histogram_label_counts(S,ws->labels,stinger_vertices_max_vertices_get(stinger_vertices_get(S)),ws->path,ws->filename,0);
  }
  return STINGER_SUCCESS;
}

void
static_clustering_help(void ** workspace) {
  printf(
"Algorithm: Static Multicontract Clustering\n"
"==========================================\n\n"
"Finds and labels independent communities in the graph using a static  \n"
"tree-based multi-contraction agglomerative clustering method that optimizes\n"
"modularity.  The clusters are recalculated each iteration.\n\n"
"Options:\n"
"\tprint_count: 1 or 0 will turn on or off printing the number of \n"
"\t             clusters found in each iteration of the algorithm\n"
"\thistogram_sizes_to_file 1 or 0 will turn on or off histogramming\n"
"\t             the cluster sizes and writing out to file\n"
"\tpath         path to store results files\n"
"\tfilename     prefix for the file\n");
  exit(-1);
}

stinger_return_t
static_clustering_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  static_clustering_workpace_t * ws = static_clustering_workspace_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  CODE_TOO=1 ./tools/stringset/stringset \
  help "static_clustering_help(workspace);"\
  ? "static_clustering_help(workspace);"\
  -? "static_clustering_help(workspace);"\
  --? "static_clustering_help(workspace);"\
  -h "static_clustering_help(workspace);"\
  --help "static_clustering_help(workspace);"\
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
			static_clustering_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  static_clustering_help(workspace);
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
		  static_clustering_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  static_clustering_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    static_clustering_help(workspace);
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
		    static_clustering_help(workspace);
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
static_clustering_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count) {
  batch += 1;

  static_clustering_workpace_t * ws = static_clustering_workspace_from_void(S, workspace);

  stinger_t * S_cluster = stinger_new();
  for(uint64_t t = 0; t < STINGER_NUMETYPES; t++) {
    STINGER_PARALLEL_FORALL_EDGES_BEGIN(S, t) {
      stinger_insert_edge(S_cluster, STINGER_EDGE_TYPE, STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_WEIGHT, STINGER_EDGE_TIME_RECENT);
    } STINGER_PARALLEL_FORALL_EDGES_END();
  }
  static_multi_contract_clustering(ws->labels, ws->scores, stinger_vertices_max_vertices_get(stinger_vertices_get(S)), S_cluster, S);
  stinger_free_all(S_cluster);

  if(ws->print_count) {
    uint64_t comm_count = 0;
    for(uint64_t v = 0; v < stinger_vertices_max_vertices_get(stinger_vertices_get(S)); v++) {
      if(stinger_vtype_get(S, v) != 0 && ws->labels[v] == v)
	comm_count++;
    }
    printf("Clusters: %ld after %ld\n", comm_count, batch);
  }
  if(ws->histogram_sizes_to_file) {
    histogram_label_counts(S,ws->labels,stinger_vertices_max_vertices_get(stinger_vertices_get(S)),ws->path,ws->filename,batch);
  }
  return STINGER_SUCCESS;
}

