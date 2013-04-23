#include "int-ht-seq.h"
#include "streaming_clustering_coefficients.h"
#include "xmalloc.h"
#include "stinger-atomics.h"
#include "stinger-utils.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * IMPLEMENTATION FUNCTIONS
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int csr_cmp(const void * a, const void * b);

static uint64_t
count_triangles(int_ht_seq_t * ht, stinger_t * S, int64_t v) {
  uint64_t count = 0;
  uint64_t count_check = 0;


  STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
    uint64_t neigh = STINGER_EDGE_DEST;

    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, neigh) {
      if(STINGER_EDGE_DEST != v) {
	count += int_ht_seq_exists(ht, STINGER_EDGE_DEST);
      }
    } STINGER_FORALL_EDGES_OF_VTX_END();

    /*
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, neigh) {
      if(STINGER_EDGE_DEST != v) {
	uint64_t neigh_neigh = STINGER_EDGE_DEST;
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
	  count_check += (neigh_neigh == STINGER_EDGE_DEST);
	} STINGER_FORALL_EDGES_OF_VTX_END();
      }
    } STINGER_FORALL_EDGES_OF_VTX_END();
    */

    int_ht_seq_insert(ht, STINGER_EDGE_DEST);
  } STINGER_FORALL_EDGES_OF_VTX_END();

  /*
  if(count_check / 2 != count) {
    printf("%s %d - hash table not working? %ld and %ld deg == %ld\n", __func__, __LINE__, count, count_check / 2, stinger_outdegree_get(S,v));
  }
  */

  int_ht_seq_empty(ht);
  return count * 2;
}

static void
count_all_triangles(stinger_t * S, int64_t nv, uint64_t * triangles) {
  OMP("omp parallel") 
  {
    int_ht_seq_t ht;
    int_ht_seq_init(&ht, 64);

    OMP("omp for nowait") 
    for(uint64_t v = 0; v < nv; v++) {
      if(stinger_outdegree_get(S,v)) {
	triangles[v] = count_triangles(&ht, S, v);
      } else {
	triangles[v] = 0;
      }
    }

    int_ht_seq_free_internal(&ht);
  }
}


streaming_clustering_coefficients_workpace_t *
streaming_clustering_coefficients_workspace_from_void(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  streaming_clustering_coefficients_workpace_t * ws = *((streaming_clustering_coefficients_workpace_t **)workspace);
  if(!ws) {
    ws = xmalloc(sizeof(streaming_clustering_coefficients_workpace_t) + sizeof(int64_t) *  
      stinger_vertices_max_vertices_get(stinger_vertices_get(S)));
    if(!ws) return NULL;
    *workspace = ws;
    ws->nv = stinger_vertices_max_vertices_get(stinger_vertices_get(S));
    ws->ntri_name = "triangles";
    ws->ntri_nr = stinger_workflow_new_named_result(wkflow, ws->ntri_name, NR_I64, ws->nv);
    ws->local_cc_name = "local_clustering_coefficients";
    ws->local_cc_nr = stinger_workflow_new_named_result(wkflow, ws->local_cc_name, NR_DBL, ws->nv);
    ws->affected = xcalloc (ws->nv, sizeof (int64_t));
    ws->global_ntri = 0;
    ws->print_global = 1;
  }
  return ws;
}

void
streaming_clustering_coefficients_help(void ** workspace) {
  printf(
"Algorithm: Streaming Clustering Coefficients\n"
"============================================\n\n"
"Track the clustering coefficients of the graph.  Performs dynamic updates\n"
"to triangle counts as edges are inserted and removed.\n"
"Options:\n"
"\tprint_global: 1 or 0 print the global clustering coefficient after \n"
"\t           each batch\n"
"\tntri_name: The name of the result that will store the number \n"
"\t           of triangles\n"
"\tlocal_cc_name: The name of the result that will store the local\n"
"\t           clustering coefficients\n");
  exit(-1);
}

stinger_return_t
streaming_clustering_coefficients_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  streaming_clustering_coefficients_workpace_t * ws = streaming_clustering_coefficients_workspace_from_void(S, wkflow, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  CODE_TOO=1 ./tools/stringset/stringset \
  help "streaming_clustering_coefficients_help(workspace);"\
  ? "streaming_clustering_coefficients_help(workspace);"\
  -? "streaming_clustering_coefficients_help(workspace);"\
  --? "streaming_clustering_coefficients_help(workspace);"\
  -h "streaming_clustering_coefficients_help(workspace);"\
  --help "streaming_clustering_coefficients_help(workspace);"\
  print_global "ws->print_global = atol(values[i]);" \
  ntri_name "stinger_workflow_delete_named_result(wkflow, ws->ntri_name); \
    ws->ntri_name = values[i]; \
    ws->ntri_nr = stinger_workflow_new_named_result(wkflow, ws->ntri_name, NR_I64, ws->nv); "\
  local_cc_name "stinger_workflow_delete_named_result(wkflow, ws->local_cc_name); \
    ws->local_cc_name = values[i]; \
    ws->local_cc_nr = stinger_workflow_new_named_result(wkflow, ws->local_cc_name, NR_I64, ws->nv); "
    
*/
/* GENERATED WITH stringset.c */
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
			streaming_clustering_coefficients_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  streaming_clustering_coefficients_help(workspace);
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
		  streaming_clustering_coefficients_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  streaming_clustering_coefficients_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    streaming_clustering_coefficients_help(workspace);
	  }
	} break;
      case 'h':
	{
	  str++; len--;
	  if(!strncmp(str, "elp", len)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* help */
	      streaming_clustering_coefficients_help(workspace);
	    }
	  }
	} break;
      case 'l':
	{
	  str++; len--;
	  if(!strncmp(str, "ocal_cc_name", len)) {
	    str += 12; len -= 12;
	    if(len == 0) {
	      /* local_cc_name */
	      stinger_workflow_delete_named_result(wkflow, ws->local_cc_name);     ws->local_cc_name = values[i];     ws->local_cc_nr = stinger_workflow_new_named_result(wkflow, ws->local_cc_name, NR_I64, ws->nv); 
	    }
	  }
	} break;
      case 'n':
	{
	  str++; len--;
	  if(!strncmp(str, "tri_name", len)) {
	    str += 8; len -= 8;
	    if(len == 0) {
	      /* ntri_name */
	      stinger_workflow_delete_named_result(wkflow, ws->ntri_name);     ws->ntri_name = values[i];     ws->ntri_nr = stinger_workflow_new_named_result(wkflow, ws->ntri_name, NR_I64, ws->nv); 
	    }
	  }
	} break;
      case 'p':
	{
	  str++; len--;
	  if(!strncmp(str, "rint_global", len)) {
	    str += 11; len -= 11;
	    if(len == 0) {
	      /* print_global */
	      ws->print_global = atol(values[i]);
	    }
	  }
	} break;
    }

  }

  free(keys); free(values);

  return STINGER_SUCCESS;
}

stinger_return_t
streaming_clustering_coefficients_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  streaming_clustering_coefficients_workpace_t * ws = streaming_clustering_coefficients_workspace_from_void(S, wkflow, workspace);

  ws->ntri = stinger_named_result_write_data(ws->ntri_nr);
  ws->local_cc = stinger_named_result_write_data(ws->local_cc_nr);

  OMP("omp parallel for")
  for (int64_t i = 0; i < ws->nv; i++) {
    ws->ntri[i] = 0;
    ws->local_cc[i] = 0;
  }

  count_all_triangles(S, ws->nv, ws->ntri);

  int64_t global_degsum = 0;
  OMP("omp parallel for")
  MTA("mta assert nodep")
  for (int64_t i = 0; i < ws->nv; i++) {
    int64_t deg = stinger_outdegree_get(S,i);
    if(deg) {
      int64_t d = deg * (deg-1);
      ws->local_cc[i] = (d ? ((double)ws->ntri[i]) / ((double)d) : 0.0);
      global_degsum += d;
      ws->global_ntri += ws->ntri[i];
    }
  }

  stinger_named_result_commit_data(ws->ntri_nr);
  stinger_named_result_commit_data(ws->local_cc_nr);

  if(ws->print_global) {
    double global_cc = (global_degsum ? ws->global_ntri / (double) global_degsum : 0.0);
    printf("Global Clustering Coeff. is %lf after %d\n", global_cc, 0);
  }
}

stinger_return_t
streaming_clustering_coefficients_after_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t * actions, int * result, int64_t count) {
  streaming_clustering_coefficients_workpace_t * ws = streaming_clustering_coefficients_workspace_from_void(S, wkflow, workspace);

  ws->ntri = stinger_named_result_write_data(ws->ntri_nr);
  ws->local_cc = stinger_named_result_write_data(ws->local_cc_nr);

  int64_t min_ts = INT64_MAX;
  OMP("omp parallel")
  {
    int64_t my_min_ts = min_ts;

    OMP("omp for nowait")
    MTA("mta assert parallel")
    for (int64_t k = 0; k < count; k++)
    {
      if(actions[k].time < my_min_ts)
	my_min_ts = actions[k].time;
    }

    if(my_min_ts < min_ts) {
      OMP("omp critical")
      if(my_min_ts < min_ts) {
	min_ts = my_min_ts;
      }
    }
  }

  OMP("omp parallel for")
  MTA("mta assert parallel")
  for (int64_t k = 0; k < ws->nv; k++)
  {
    ws->affected[k] = 0;
  }

  for (uint64_t k = 0; k < count; k++) {
    const int64_t i = actions[k].source;
    const int64_t j = actions[k].dest;

    if (result[k] == 1 && i > j) continue;

    if (i != j && i >= 0 && result[k]) {
      int64_t i_deg, j_deg;
      i_deg = stinger_outdegree_get(S, i);

      int64_t * i_neighbors = xmalloc (i_deg * sizeof(*i_neighbors));
      int64_t * i_neighbors_ts = xmalloc (i_deg * sizeof(*i_neighbors_ts));
      size_t i_degree;
      stinger_gather_successors (S, i, &i_degree, i_neighbors, NULL, i_neighbors_ts, NULL, NULL, i_deg);
      assert (i_degree == i_deg);

      int64_t ** pointers = xmalloc(i_degree * sizeof(int64_t *));
      for(int64_t j = 0; j < i_degree; j++) {
	pointers[j] = i_neighbors+j; 
      }

      qsort (pointers, i_degree, sizeof(int64_t *), csr_cmp);

      int64_t * tmp_buffer = xmalloc(i_degree * sizeof(int64_t));
      for(int64_t j = 0; j < i_degree; j++) {
	tmp_buffer[j] = *pointers[j];
      }
      for(int64_t j = 0; j < i_degree; j++) {
	i_neighbors[j] = tmp_buffer[j];
      }
      for(int64_t i = 0; i < i_degree; i++) {
	tmp_buffer[i] = i_neighbors_ts[pointers[i] - i_neighbors];
      }
      for(int64_t i = 0; i < i_degree; i++) {
	i_neighbors_ts[i] = tmp_buffer[i];
      }

      int64_t incr = 0;

      j_deg = stinger_outdegree_get(S, j);
      size_t j_degree;

      if (!j_deg)
	continue;

      int64_t * j_neighbors = xmalloc (j_deg * sizeof(*j_neighbors));
      int64_t * j_neighbors_ts = xmalloc (j_deg * sizeof(*j_neighbors_ts));
      stinger_gather_successors (S, j, &j_degree, j_neighbors, NULL, j_neighbors_ts, NULL, NULL, j_deg);
      assert (j_degree == j_deg);

      for (int64_t k2 = 0; k2 < j_degree; k2++)
      {
	int64_t u = j_neighbors[k2];
	int64_t where = find_in_sorted (u, i_degree, i_neighbors);
	if (where >= 0) {
	  int64_t neighbor_i_ts = i_neighbors_ts[where];
	  int64_t neighbor_j_ts = j_neighbors_ts[k2];
	  if (neighbor_i_ts < min_ts && neighbor_j_ts < min_ts) {
	    incr += 1;
	    stinger_int64_fetch_add (&ws->ntri[u], 2);
	    stinger_int64_fetch_add (&ws->affected[u], 1);
	  } else if (i < u && neighbor_i_ts < min_ts) {
	    incr += 1;
	    stinger_int64_fetch_add (&ws->ntri[u], 2);
	    stinger_int64_fetch_add (&ws->affected[u], 1);
	  } else if (j < u && neighbor_j_ts < min_ts) {
	    incr += 1;
	    stinger_int64_fetch_add (&ws->ntri[u], 2);
	    stinger_int64_fetch_add (&ws->affected[u], 1);
	  } else if (neighbor_i_ts >= min_ts && neighbor_j_ts >= min_ts) {
	    if (i < u && j < u) {
	      incr += 1;
	      stinger_int64_fetch_add (&ws->ntri[u], 2);
	      stinger_int64_fetch_add (&ws->affected[u], 1);
	    }
	  }
	}
      }

      free (j_neighbors_ts);
      free (j_neighbors);
      free (pointers);
      free (i_neighbors_ts);
      free (i_neighbors);

      stinger_int64_fetch_add (&(ws->ntri[i]), 2*incr);
      stinger_int64_fetch_add (&(ws->ntri[j]), 2*incr);
      ws->global_ntri += (6*incr);
    }
  }

  uint64_t global_degsum = 0;
  OMP("omp parallel for")
  MTA("mta assert nodep")
  for (int64_t i = 0; i < ws->nv; i++) {
    int64_t deg = stinger_outdegree_get(S, i);
    int64_t d = deg * (deg-1);
    global_degsum += d;
  }

  /* Update the local clustering coefficient for affected vertices */
  OMP("omp parallel for")
  MTA("mta assert parallel")
  for (int64_t k = 0; k < ws->nv; k++) {
    if (ws->affected[k]) {
      double d = (double) stinger_outdegree_get(S, k);
      double new_local_cc = (d > 1 ? ws->ntri[k] / (d * (d-1)) : 0.0);
      ws->local_cc[k] = new_local_cc;
    }
  }

  stinger_named_result_commit_data(ws->ntri_nr);
  stinger_named_result_commit_data(ws->local_cc_nr);

  if(ws->print_global) {
    double global_cc = (global_degsum ? ws->global_ntri / (double) global_degsum : 0.0);
    printf("global cc is %lf\n", global_cc);
    printf("Global Clustering Coeff. is %lf after %ld\n", global_cc, batch);
  }
}
