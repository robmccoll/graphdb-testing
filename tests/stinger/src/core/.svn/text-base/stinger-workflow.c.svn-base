#include "stinger.h"
#include "stinger-workflow.h"
#include "xmalloc.h"
#include "csv.h"
#include "timer.h"

/**
* @file stinger-workflow.c
* @brief Implementation of the stinger_workflow and stinger_named_result
* @author Rob McColl
* @date 2013-02-13
*/

/**
* @brief Create a new STINGER workflow.
*
* @param S The STINGER that will store the graph in this workflow.
*
* @return A pointer to the newly created workflow.
*/
stinger_workflow_t * 
stinger_workflow_new(stinger_t * S) {

  stinger_workflow_t * workflow = xmalloc(sizeof(stinger_workflow_t));

  workflow->integrated		= 0;
  workflow->S			= S;
  workflow->alg_count		= 0;
  workflow->alg_size		= 10;
  workflow->algorithms		= xmalloc(workflow->alg_size * sizeof(stinger_alg_t));
  workflow->stream_count	= 0;
  workflow->stream_size		= 10;
  workflow->streams		= xmalloc(workflow->stream_size * sizeof(stinger_stream_t));
  workflow->named_result_count	= 0;
  workflow->named_result_size	= 10;
  workflow->named_results	= xmalloc(workflow->named_result_size * sizeof(stinger_named_result_t *));

  return workflow;
}

/**
* @brief Free a STINGER workflow and its internal memory. Sets the pointer to NULL.
*
* @param workflow The workflow.
*
* @return Returns a NYLL pointer.
*/
stinger_workflow_t * 
stinger_workflow_free(stinger_workflow_t ** workflow) {

  if(*workflow) {
    if((*workflow)->algorithms) {
      free((*workflow)->algorithms);
    }
    if((*workflow)->streams) {
      free((*workflow)->streams);
    }
    free(*workflow);
  }

  return *workflow;
}

/**
* @brief Register an algorithm to the workflow.
*
* @param workflow The workflow to which you want to register the algorithm.
* @param name A unique name to refer to the algorithm.
* @param workspace A pointer to the workspace or NULL.
* @param init A function pointer to the init function or NULL.
* @param settings A fuction pointer to the settings function or NULL.
* @param before_batch A function pointer that will be called before each batch is applied or NULL.
* @param after_batch A function pointer that will be called after each batch is applied or NULL.
* @param cleanup A function pointer that will be called after one of the batch or initial functions returns STINGER_REMOVE 
*		or when all streams have ended or NULL.
*
* @return The same workflow as the input for chaining.
*/
stinger_workflow_t * 
stinger_workflow_register_alg(stinger_workflow_t * workflow, char * name, void * workspace,
  stinger_return_t (*init)(stinger_t *, stinger_workflow_t *, void **),
  stinger_return_t (*settings)	(stinger_t *, stinger_workflow_t *, void **, char *),
  stinger_return_t (*before_batch)(stinger_t *, stinger_workflow_t *, void **, int64_t, edge_action_t *, int64_t),
  stinger_return_t (*after_batch)(stinger_t *, stinger_workflow_t *, void **, int64_t, edge_action_t *, int *, int64_t),
  stinger_return_t (*cleanup)(stinger_t *, stinger_workflow_t *, void **)) {
    
    if(workflow->alg_count >= workflow->alg_size) {
      workflow->alg_size *= 2;
      workflow->algorithms = xrealloc(workflow->algorithms, workflow->alg_size * sizeof(stinger_alg_t));
    }

    workflow->algorithms[workflow->alg_count].name		= name;
    workflow->algorithms[workflow->alg_count].workspace		= workspace;
    workflow->algorithms[workflow->alg_count].init		= init;
    workflow->algorithms[workflow->alg_count].settings		= settings;
    workflow->algorithms[workflow->alg_count].before_batch	= before_batch;
    workflow->algorithms[workflow->alg_count].after_batch	= after_batch;
    workflow->algorithms[workflow->alg_count].cleanup		= cleanup;

    workflow->alg_count++;

    return workflow;
}

/**
* @brief Remove an algorithm from the workflow by its index.
*
* @param workflow The workflow from which to remove the algorithm.
* @param index Index of the algorithm within the workflow.
*
* @return The same workflow as the input for chaining.
*/
stinger_workflow_t *
stinger_workflow_unregister_alg_by_index(stinger_workflow_t * workflow, int64_t index) {
  for(uint64_t i = index; i < workflow->alg_count-1; i++) {
    workflow->algorithms[i].name	  = workflow->algorithms[i+1].name;
    workflow->algorithms[i].workspace	  = workflow->algorithms[i+1].workspace;
    workflow->algorithms[i].init	  = workflow->algorithms[i+1].init;
    workflow->algorithms[i].settings	  = workflow->algorithms[i+1].settings;
    workflow->algorithms[i].before_batch  = workflow->algorithms[i+1].before_batch;
    workflow->algorithms[i].after_batch	  = workflow->algorithms[i+1].after_batch;
    workflow->algorithms[i].cleanup	  = workflow->algorithms[i+1].cleanup;
  }
  workflow->alg_count--;
  return workflow;
}

/**
* @brief Remove an algorithm from the workflow by its name.
*
* Names that aren't found are silently ignored.  This function <strong>will not</strong>
* call the provided cleanup function.
*
* @param workflow The workflow from which to remove the algorithm.
* @param name The name of the workflow to be removed (must match).
*
* @return The same workflow as the input for chaining.
*/
stinger_workflow_t *
stinger_workflow_unregister_alg(stinger_workflow_t * workflow, char * name) {
  for(uint64_t i = 0; i < workflow->alg_count; i++) {
    if(0 == strcmp(workflow->algorithms[i].name, name)) {
      return stinger_workflow_unregister_alg_by_index(workflow, i);
    }
  }
  return workflow;
}

/**
* @brief Register stream to the workflow.
*
* @param workflow The workflow to which you want to register the stream.
* @param name A unique name to refer to the stream.
* @param workspace A pointer to the workspace or NULL.
* @param init A function pointer to the init function or NULL.
* @param settings A fuction pointer to the settings function or NULL.
* @param stream_batch A function pointer to a function that will generate batches.
* @param cleanup A function pointer that will be called after one of the batch or initial functions returns STINGER_REMOVE 
*		or when all streams have ended or NULL.
*
* @return The same workflow as the input for chaining.
*/
stinger_workflow_t * 
stinger_workflow_register_stream(stinger_workflow_t * workflow, char * name, void * workspace,
  stinger_return_t	(*init)		(stinger_t *, stinger_workflow_t *, void **),
  stinger_return_t	(*settings)	(stinger_t *, stinger_workflow_t *, void **, char *),
  stinger_return_t	(*stream_batch)	(stinger_t *, stinger_workflow_t *, void **, int64_t, edge_action_t **, int64_t*),
  stinger_return_t	(*cleanup)	(stinger_t *, stinger_workflow_t *, void **)) {
    
    if(workflow->stream_count >= workflow->stream_size) {
      workflow->stream_size *= 2;
      workflow->streams = xrealloc(workflow->streams, workflow->stream_size * sizeof(stinger_stream_t));
    }

    workflow->streams[workflow->stream_count].name	    = name;
    workflow->streams[workflow->stream_count].workspace	    = workspace;
    workflow->streams[workflow->stream_count].init	    = init;
    workflow->streams[workflow->stream_count].settings	    = settings;
    workflow->streams[workflow->stream_count].stream_batch  = stream_batch;
    workflow->streams[workflow->stream_count].cleanup	    = cleanup;

    workflow->stream_count++;

    return workflow;
}

/**
* @brief Remove a stream from a workflow by its index in the workflow.
*
* @param workflow The workflow from which to remove the stream.
* @param index Which stream to remove in the internal order.
*
* @return The same workflow as the input for chaining.
*/
stinger_workflow_t *
stinger_workflow_unregister_stream_by_index(stinger_workflow_t * workflow, int64_t index) {
  for(uint64_t i = index; i < workflow->stream_count-1; i++) {
    workflow->streams[i].name	      = workflow->streams[i+1].name;
    workflow->streams[i].workspace    = workflow->streams[i+1].workspace;
    workflow->streams[i].init         = workflow->streams[i+1].init;
    workflow->streams[i].settings     = workflow->streams[i+1].settings;
    workflow->streams[i].stream_batch = workflow->streams[i+1].stream_batch;
    workflow->streams[i].cleanup      = workflow->streams[i+1].cleanup;
  }
  workflow->stream_count--;
  return workflow;
}

/**
* @brief Remove a stream from a workflow by its name.
*
* @param workflow The workflow from which to remove the stream.
* @param name Name of the algorithm to remove.  <string>Must match</strong>.
*
* @return The same workflow as the input for chaining.
*/
stinger_workflow_t *
stinger_workflow_unregister_stream(stinger_workflow_t * workflow, char * name) {
  for(uint64_t i = 0; i < workflow->stream_count; i++) {
    if(0 == strcmp(workflow->streams[i].name, name)) {
      return stinger_workflow_unregister_stream_by_index(workflow, i);
    }
  }
  return workflow;
}

void
stinger_workflow_help(stinger_workflow_t * workflow) {
  printf(
"STINGER Workflow Help\n"
"=====================\n\n"
"This program is using a STINGER workflow which will execute a combination \n"
"of various streams and algorithms.  streams create the initial static \n"
"graph and can add and remove edges and vertices to and from the graph as \n"
"time goes on.  Algorithms can run on the initial graph and can update their \n"
"metrics as the graph changes.  By default, all algorithms and streams will\n"
"be run.  You can select a subset of algorithms, streams, or both to run \n"
"using the following options:\n"
"\talgs algorithm1,algorithm2\n\tstreams stream1,stream2\n\n"
"Available streams:\n");
  for(uint64_t g = 0; g < workflow->stream_count; g++) {
    printf("\t%s\n", workflow->streams[g].name);
  }
  printf("\nAvailable Algorithms:\n");
  for(uint64_t a = 0; a < workflow->alg_count; a++) {
    printf("\t%s\n", workflow->algorithms[a].name);
  }
  printf("\nAdditionally, you can pass configuration data to algorithms and streams\n"
  "that support it via this syntax:\n"
  "\t@alg: <algname> key1:value1,key2,value2\n\t@stream: <streamname> key1:value1,key2:value2\n");
  printf("\nSome algorithms and streams may provide their own help through:\n\t@alg: <algname> help\n");
  exit(-1);
}

void
stinger_workflow_to_keyvalue(char * str, int len, char *** keys, char *** values, int * num) {
  char ** ks;
  char ** vs;
  int n;

  char ** fields = NULL; uint64_t * lengths = NULL; uint64_t fieldsSize = 0; uint64_t count = 0;
  char ** s_fields = NULL; uint64_t * s_lengths = NULL; uint64_t s_fieldsSize = 0; uint64_t s_count = 0;
  splitLineCSVDynamicInPlace(',', str, len, &fields, &lengths, &fieldsSize, &count);

  n = count;
  ks = xmalloc(sizeof(char *) * count);
  vs = xmalloc(sizeof(char *) * count);

  for(uint64_t i = 0; i < count; i++) {
    splitLineCSVDynamicInPlace(':', fields[i], lengths[i], &s_fields, &s_lengths, &s_fieldsSize, &s_count);
    ks[i] = s_fields[0];
    if(s_count >= 2) {
      vs[i] = s_fields[1];
    } 
  }

  (*keys) = ks; (*values) = vs; (*num) = n;
  free(fields); free(lengths);
  free(s_fields); free(s_lengths);
}

void
stinger_workflow_parse_args(int argc, char ** argv, stinger_workflow_t * workflow) {

  for(uint64_t i = 1; i < argc; i += 2) {
    char * str = argv[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:
  CODE_TOO=1 ./tools/stringset/stringset  algs " \
  { if(i+1 < argc) { \
    char * buf = NULL; int64_t bufSize = 0; char ** fields = NULL; uint64_t * lengths = NULL; uint64_t fieldsSize = 0; uint64_t count = 0; \
    splitLineCSVDynamic(',', argv[i+1], strlen(argv[i+1]), &buf, &bufSize, &fields, &lengths, &fieldsSize, &count); \
    for(uint64_t a = 0; a < workflow->alg_count; a++) { \
      int found = 0; \
      for(uint64_t f = 0; f < count; f++) { \
	if(0 == strcmp(workflow->algorithms[a].name, fields[f])) { \
	  found = 1; \
	} \
      } \
      if(0 == found) { \
	stinger_workflow_unregister_alg_by_index(workflow, a); \
	a--;
      } \
    } \
    free(buf); free(fields); free(lengths); \
  } }" \
  streams " \
  { if(i+1 < argc) { \
    char * buf = NULL; int64_t bufSize = 0; char ** fields = NULL; uint64_t * lengths = NULL; uint64_t fieldsSize = 0; uint64_t count = 0; \
    splitLineCSVDynamic(',', argv[i+1], strlen(argv[i+1]), &buf, &bufSize, &fields, &lengths, &fieldsSize, &count); \
    for(uint64_t g = 0; g < workflow->stream_count; g++) { \
      int found = 0; \
      for(uint64_t f = 0; f < count; f++) { \
	if(0 == strcmp(workflow->streams[g].name, fields[f])) { \
	  found = 1; \
	} \
      } \
      if(0 == found) { \
	stinger_workflow_unregister_stream_by_index(workflow, g); \
	g--;
      } \
    } \
    free(buf); free(fields); free(lengths); \
  } }" \
  @alg: " \
  { if(i+2 < argc) { \
    for(uint64_t a = 0; a < workflow->alg_count; a++) { \
      if(0 == strcmp(workflow->algorithms[a].name, argv[i+1])) { \
	if(workflow->algorithms[a].settings) { \
	  workflow->algorithms[a].settings(workflow->S, workflow, &(workflow->algorithms[a].workspace), argv[i+2]); \
	} \
	i++; \
	break; \
      } \
    } \
  } }" \
  @stream: " \
  { if(i+2 < argc) { \
    for(uint64_t g = 0; g < workflow->stream_count; g++) { \
      if(0 == strcmp(workflow->streams[g].name, argv[i+1])) { \
	if(workflow->streams[g].settings) { \
	  workflow->streams[g].settings(workflow->S, workflow, &(workflow->streams[g].workspace), argv[i+2]); \
	} \
	i++; \
	break; \
      } \
    } \
  } }" \
  help "stinger_workflow_help(workflow);" \
  ? "stinger_workflow_help(workflow);" \
  -? "stinger_workflow_help(workflow);" \
  --? "stinger_workflow_help(workflow);" \
  -h "stinger_workflow_help(workflow);" \
  --help "stinger_workflow_help(workflow);"
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
			stinger_workflow_help(workflow);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  stinger_workflow_help(workflow);
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
		  stinger_workflow_help(workflow);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  stinger_workflow_help(workflow);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    stinger_workflow_help(workflow);
	  }
	} break;
      case '@':
	{
	  str++; len--;
	  if(len) switch(*str) {
	    case 'a':
	      {
		str++; len--;
		if(!strncmp(str, "lg:", len)) {
		  str += 3; len -= 3;
		  if(len == 0) {
		    /* @alg: */
		    { if(i+2 < argc) {     for(uint64_t a = 0; a < workflow->alg_count; a++) {       if(0 == strcmp(workflow->algorithms[a].name, argv[i+1])) { 	if(workflow->algorithms[a].settings) { 	  workflow->algorithms[a].settings(workflow->S, workflow, &(workflow->algorithms[a].workspace), argv[i+2]); 	} 	i++; 	break;       }     }   } }
		  }
		}
	      } break;
	    case 's':
	      {
		str++; len--;
		if(!strncmp(str, "tream:", len)) {
		  str += 6; len -= 6;
		  if(len == 0) {
		    /* @stream: */
		    { if(i+2 < argc) {     for(uint64_t g = 0; g < workflow->stream_count; g++) {       if(0 == strcmp(workflow->streams[g].name, argv[i+1])) { 	if(workflow->streams[g].settings) { 	  workflow->streams[g].settings(workflow->S, workflow, &(workflow->streams[g].workspace), argv[i+2]); 	} 	i++; 	break;       }     }   } }
		  }
		}
	      } break;
	  }
	} break;
      case 'a':
	{
	  str++; len--;
	  if(!strncmp(str, "lgs", len)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* algs */
	      { if(i+1 < argc) {     char * buf = NULL; int64_t bufSize = 0; char ** fields = NULL; uint64_t * lengths = NULL; uint64_t fieldsSize = 0; uint64_t count = 0;     splitLineCSVDynamic(',', argv[i+1], strlen(argv[i+1]), &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);     for(uint64_t a = 0; a < workflow->alg_count; a++) {       int found = 0;       for(uint64_t f = 0; f < count; f++) { 	if(0 == strcmp(workflow->algorithms[a].name, fields[f])) { 	  found = 1; 	}       }       if(0 == found) { 	stinger_workflow_unregister_alg_by_index(workflow, a); 	a--;
																																																																	       }     }     free(buf); free(fields); free(lengths);   } }
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
	      stinger_workflow_help(workflow);
	    }
	  }
	} break;
      case 's':
	{
	  str++; len--;
	  if(!strncmp(str, "treams", len)) {
	    str += 6; len -= 6;
	    if(len == 0) {
	      /* streams */
	      { if(i+1 < argc) {     char * buf = NULL; int64_t bufSize = 0; char ** fields = NULL; uint64_t * lengths = NULL; uint64_t fieldsSize = 0; uint64_t count = 0;     splitLineCSVDynamic(',', argv[i+1], strlen(argv[i+1]), &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);     for(uint64_t g = 0; g < workflow->stream_count; g++) {       int found = 0;       for(uint64_t f = 0; f < count; f++) { 	if(0 == strcmp(workflow->streams[g].name, fields[f])) { 	  found = 1; 	}       }       if(0 == found) { 	stinger_workflow_unregister_stream_by_index(workflow, g); 	g--;
																																																																	       }     }     free(buf); free(fields); free(lengths);   } }
	    }
	  }
	} break;
    }
    /* END GENERATED CODE */
  }
}

void
stinger_workflow_run(int argc, char ** argv, stinger_workflow_t * workflow, int print_time) {
  stinger_workflow_parse_args(argc, argv, workflow);

  if((!workflow->stream_count) || (!workflow->alg_count))
    return;

  /* init all streams */
  for(uint64_t i = 0; i < workflow->stream_count; i++) {
    if(workflow->streams[i].init) {
      workflow->streams[i].init(workflow->S, workflow, &(workflow->streams[i].workspace));
    }
  }

  /* init all algorithms */
  for(uint64_t i = 0; i < workflow->alg_count; i++) {
    if(workflow->algorithms[i].init) {
      workflow->algorithms[i].init(workflow->S, workflow, &(workflow->algorithms[i].workspace));
    }
  }

  /* main run loops */
  if(!workflow->integrated) {
    edge_action_t * actions = NULL;
    int64_t actions_size    = 0;
    int64_t actions_len	    = 0;

    int * result	= NULL;
    int64_t result_size = 0;

    int64_t batch_num = 0;

    /* while streams and algs are still running...
     * for each stream { stream, pre-batch each alg, apply batch, post-batch each alg } */
    while(workflow->stream_count && workflow->alg_count) {
      for(uint64_t i = 0; i < workflow->stream_count; i++) {
	if(workflow->streams[i].stream_batch) {
	  printf("Starting batch: %s\n", workflow->streams[i].name);
	  tic();
	  /* gather actions */
	  if(actions_size < actions_len) {
	    actions_size = actions_len;
	  }
	  if(STINGER_REMOVE == workflow->streams[i].stream_batch(workflow->S, workflow,
	    &(workflow->streams[i].workspace), batch_num, &actions, &actions_len)) {
	      if(workflow->streams[i].cleanup) {
		workflow->streams[i].cleanup(workflow->S, workflow, &(workflow->streams[i].workspace));
	      }
	      stinger_workflow_unregister_stream_by_index(workflow, i);
	  }

	  /* batch all algorithms */
	  for(uint64_t i = 0; i < workflow->alg_count; i++) {
	    if(workflow->algorithms[i].before_batch) {
	      if(STINGER_REMOVE == workflow->algorithms[i].before_batch(workflow->S, workflow,
		&(workflow->algorithms[i].workspace), batch_num, actions, actions_len)) {
		  if(workflow->algorithms[i].cleanup) {
		    workflow->algorithms[i].cleanup(workflow->S, workflow, &(workflow->algorithms[i].workspace));
		  }
		  stinger_workflow_unregister_alg_by_index(workflow, i);
	      }
	    }
	  }

	  /* apply to stinger */
	  if(actions_len > result_size) {
	    result_size = actions_len;
	    result = xrealloc(result, result_size * sizeof(int));
	  }

	  OMP("omp parallel for")
	  for(uint64_t a = 0; a < actions_len; a++) {
	    if(actions[a].source >= 0) {
	      result[a] = stinger_incr_edge_pair(workflow->S, actions[a].type, actions[a].source, 
		actions[a].dest, actions[a].weight, actions[a].time);
	    }
	  }
	  OMP("omp parallel for")
	  for(uint64_t a = 0; a < actions_len; a++) {
	    if(actions[a].source < 0) {
	      result[a] = stinger_remove_edge_pair(workflow->S, actions[a].type, 
	        ~actions[a].source, ~actions[a].dest);
	    }
	  }

	  /* post processes algorithms */
	  for(uint64_t i = 0; i < workflow->alg_count; i++) {
	    if(workflow->algorithms[i].after_batch) {
	      if(STINGER_REMOVE == workflow->algorithms[i].after_batch(workflow->S, workflow,
		&(workflow->algorithms[i].workspace), batch_num, actions, result, actions_len)) {
		  if(workflow->algorithms[i].cleanup) {
		    workflow->algorithms[i].cleanup(workflow->S, workflow, &(workflow->algorithms[i].workspace));
		  }
		  stinger_workflow_unregister_alg_by_index(workflow, i);
	      }
	    }
	  }
	  double bt = toc();
	  printf("Finished: %s in %20.15e seconds.  %20.15e edges per second.\n", 
	    workflow->streams[i].name, bt, ((double)actions_len) / bt);

	  batch_num++;
	}
      }
    }
  } else {
    /* TODO */
    fprintf(stderr, "NOT IMPLEMENTED (YET)... sorry!\n");
  }

  /* cleanup all streams */
  for(uint64_t i = 0; i < workflow->stream_count; i++) {
    if(workflow->streams[i].cleanup) {
      workflow->streams[i].cleanup(workflow->S, workflow, &(workflow->streams[i].workspace));
    }
  }

  /* cleanup all algorithms */
  for(uint64_t i = 0; i < workflow->alg_count; i++) {
    if(workflow->algorithms[i].cleanup) {
      workflow->algorithms[i].cleanup(workflow->S, workflow, &(workflow->algorithms[i].workspace));
    }
  }
}

stinger_named_result_t *
stinger_workflow_new_named_result(stinger_workflow_t * workflow, char * name, stinger_named_result_type_t type, uint64_t count) {
    if(workflow->named_result_count >= workflow->named_result_size) {
      workflow->named_result_size *= 2;
      workflow->named_results = xrealloc(workflow->named_results, workflow->named_result_size * sizeof(stinger_named_result_t));
    }

    uint64_t size = sizeof(stinger_named_result_t);
    switch (type) {
      default :
      case NR_I64:
      size += count * sizeof(int64_t);
      break;
      case NR_DBL:
      size += count * sizeof(double);
      break;
      case NR_U8:
      size += count * sizeof(uint8_t);
      break;
      case NR_I64PAIRS:
      size += count * sizeof(int64_t)*2;
      break;
    }

    uint64_t cur = workflow->named_result_count;
    workflow->named_result_count++;

    workflow->named_results[cur] = xmalloc(size);
    workflow->named_results[cur]->type = type;
    workflow->named_results[cur]->elements= count;

    strncpy(workflow->named_results[cur]->name, name, 1023);

    return workflow->named_results[cur];
}

void
stinger_workflow_write_named_results(stinger_workflow_t * workflow, char * path, uint64_t batch) {
  for(uint64_t n = 0; n < workflow->named_result_count; n++) {
    char filename[1024];
    sprintf(filename, "%s/%s.%ld.csv", path, workflow->named_results[n]->name, batch);
    FILE * fp = fopen(filename, "w");
    if(fp) {
      switch (workflow->named_results[n]->type) {
	default :
	case NR_I64:
	  csvIfIDExistsint64(fp, ',', workflow->S, NULL, stinger_vertices_max_vertices_get(stinger_vertices_get(workflow->S)), (int64_t *)workflow->named_results[n]->data);
	break;

	case NR_DBL:
	  csvIfIDExistsdouble(fp, ',', workflow->S, NULL, stinger_vertices_max_vertices_get(stinger_vertices_get(workflow->S)), (double *)workflow->named_results[n]->data);
	break;

	case NR_U8:
	  csvIfIDExistsint8(fp, ',', workflow->S, NULL, stinger_vertices_max_vertices_get(stinger_vertices_get(workflow->S)), workflow->named_results[n]->data);
	break;

	case NR_I64PAIRS:
	{
	  int64_t * data = (int64_t *)workflow->named_results[n]->data;
	  for(uint64_t v = 0; v < workflow->named_results[n]->elements; v++) {
	    if(data[v*2] || data[v*2+1]) {
	      fprintf(fp, "%ld, %ld\n", data[v*2], data[v*2+1]);
	    }
	  }
	}
	break;
      }
      fclose(fp);
    }
  }
}


stinger_return_t
stinger_workflow_delete_named_result(stinger_workflow_t * workflow, char * name) {
  for(uint64_t a = 0; a < workflow->named_result_count; a++) {
    if(0 == strcmp(workflow->named_results[a]->name, name)) {
      free(workflow->named_results[a]);
      for(uint64_t i = a; i < workflow->named_result_count-1; i++) {
	workflow->named_results[i] = workflow->named_results[i+1];
      }
      workflow->named_result_count--;
      return STINGER_SUCCESS;
    }
  }
  return STINGER_NOT_FOUND;
}

stinger_named_result_t *
stinger_workflow_get_named_result(stinger_workflow_t * workflow, char * name) {
  for(uint64_t i = 0; i < workflow->named_result_count; i++) {
    if(0 == strcmp(workflow->named_results[i]->name, name)) {
      return workflow->named_results[i];
    }
  }
  return NULL;
}


stinger_named_result_type_t
stinger_named_result_type_get(stinger_named_result_t * nr) {
  return nr->type;
}

char *
stinger_named_result_name_get(stinger_named_result_t * nr) {
  return nr->name;
}

uint64_t
stinger_named_result_count_get(stinger_named_result_t * nr) {
  return nr->elements;
}

void *
stinger_named_result_read_data(stinger_named_result_t * nr) {
  return nr->data;
}

void *
stinger_named_result_write_data(stinger_named_result_t * nr) {
  return nr->data;
}

void
stinger_named_result_commit_data(stinger_named_result_t * nr) {
  /* Nothing for now */
}
