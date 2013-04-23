#include "csv_stream.h"
#include "xmalloc.h"
#include "stinger-physmap.h"
#include "csv.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * TODO:
 * - change input config to k-v pairs of vfield:vtype (or similar)
 * - make time, weight, etc, optional
 * - go to string mapping for types and such (where applicable)
 * - enable using functions where applicable
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

csv_stream_t *
csv_stream_from_void(stinger_t * S, void ** workspace) {
  csv_stream_t * stream = *((csv_stream_t **)workspace);
  if(!stream) {
    stream = csv_stream_new("edges.csv", 1000, 1000, 1000);
    *workspace = stream;
  }
  return stream;
}

csv_stream_t *
csv_stream_new(char * file, uint64_t edges_init, uint64_t edges_per_batch, uint64_t num_batches) {
  csv_stream_t * stream = xmalloc(sizeof(csv_stream_t));
  if(stream) {
    stream->edges_init = edges_init;
    stream->edges_per_batch = edges_per_batch;
    stream->num_batches = num_batches;
    stream->file = file;
    stream->fp = fopen(file, "r");
    stream->sep = ',';
    stream->list = NULL;

    stream->buf = NULL;
    stream->bufSize = 0;
    stream->fields = NULL;
    stream->lengths = NULL;
    stream->fieldsSize = 0;
    stream->count = 0;
  }
  return stream;
}

stinger_return_t
csv_stream_init(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace) {
  csv_stream_t * stream = csv_stream_from_void(S, workspace);

  uint64_t src_vtx, dest_vtx;

  if(stream->fp) {
    for(uint64_t e = 0; e < stream->edges_init && !feof(stream->fp); e++) {
      readCSVLineDynamic(stream->sep, stream->fp, &stream->buf, &stream->bufSize, &stream->fields, &stream->lengths, &stream->fieldsSize, &stream->count);

      csv_emap_t * cur = stream->list;
      while(cur) {
	if(stinger_mapping_create(S, stream->fields[cur->src_field], stream->lengths[cur->src_field], &src_vtx) == 1) {
	  stinger_vtype_set(S, src_vtx, cur->src_type);
	}

	if(stinger_mapping_create(S, stream->fields[cur->dest_field], stream->lengths[cur->dest_field], &dest_vtx) == 1) {
	  stinger_vtype_set(S, dest_vtx, cur->dest_type);
	}
	if(src_vtx != -1 && dest_vtx != -1) {
	  stinger_incr_edge_pair(S, cur->edge_type, src_vtx, dest_vtx, cur->weight_field != -1 ? atol(stream->fields[cur->weight_field]) : 1, 
	    cur->time_field != -1 ? atol(stream->fields[cur->time_field]) : 0);
	}
	cur = cur->next;
      }
    }
    return STINGER_SUCCESS;
  } else {
    return STINGER_GENERIC_ERR;
  }
}

stinger_return_t
csv_stream_batch(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, int64_t batch, edge_action_t ** actions, int64_t * count) {
  csv_stream_t * stream = csv_stream_from_void(S, workspace);

  if(stream->fp) {
    if(*count < stream->edges_per_batch) {
      *actions = xrealloc(*actions, sizeof(edge_action_t) * stream->edges_per_batch);
      if(!(*actions))
	return STINGER_ALLOC_FAILED;
    }

    *count = stream->edges_per_batch;

    uint64_t src_vtx, dest_vtx;

    for(uint64_t e = 0; e < stream->edges_per_batch && !feof(stream->fp); e++) {
      readCSVLineDynamic(stream->sep, stream->fp, &stream->buf, &stream->bufSize, &stream->fields, &stream->lengths, &stream->fieldsSize, &stream->count);

      csv_emap_t * cur = stream->list;
      while(cur) {
	if(stinger_mapping_create(S, stream->fields[cur->src_field], stream->lengths[cur->src_field], &src_vtx) == 1) {
	  stinger_vtype_set(S, src_vtx, cur->src_type);
	}

	if(stinger_mapping_create(S, stream->fields[cur->dest_field], stream->lengths[cur->dest_field], &dest_vtx) == 1) {
	  stinger_vtype_set(S, dest_vtx, cur->dest_type);
	}
	if(src_vtx != -1 && dest_vtx != -1) {
	  (*actions)[e].type	  = cur->edge_type;
	  (*actions)[e].source	  = src_vtx;
	  (*actions)[e].dest	  = dest_vtx;
	  (*actions)[e].weight	  = cur->weight_field != -1 ? atol(stream->fields[cur->weight_field]) : 1;
	  (*actions)[e].time	  = cur->time_field != -1 ? atol(stream->fields[cur->time_field]) : 1;
	}
	cur = cur->next;
      }
    }

    if(batch == (stream->num_batches - 1))
      return STINGER_REMOVE;
    else
      return STINGER_SUCCESS;
  } else {
    return STINGER_GENERIC_ERR;
  }
}

void
csv_stream_help(void ** workspace) {
  printf(
"Stream: CSV stream\n"
"=====================\n\n"
"Generates the graph from an input CSV file.\n\n"
"Options:\n"
"\tedges_init:	   Number of initial edges to generate\n"
"\tedges_per_batch:Number of edges to generate per batch\n"
"\tnum_batches:	   Number of batches to generate before stopping\n"
"\tfile:           Name of the file including path\n"
"\tseparator:      The separator character\n"
"\tedge:	   An edge, currently defined as a set of integers \n"
"\t                separated by |'s that represent the field of text\n"
"\t                that maps to the source vertex, the integer type \n"
"\t                of the source the field to map the destination, \n"
"\t                the type of the destination, the field to parse the \n"
"\t                integer weight or -1, the field to parse the integer\n"
"\t                timestamp or -1, and the integer edge type, -1's indicate\n"
"\t                that the field does not exist\n\n"
"src_field|src_type|dest_field|dest_type|weight_field|time_field|edge_type\n\n"
"Each edge will be mapped from each row of the CSV input file.  You may map as\n"
"many edges from each file as you see fit.\n");
exit(-1);
}

stinger_return_t
csv_stream_settings(stinger_t * S, stinger_workflow_t * wkflow, void ** workspace, char * settings) {
  char ** keys = NULL, ** values = NULL; int num = 0; 
  stinger_workflow_to_keyvalue(settings, strlen(settings), &keys, &values, &num);

  csv_stream_t * ws = csv_stream_from_void(S, workspace);

  for(uint64_t i = 0; i < num; i++) {
    char * str = keys[i];
    int len = strlen(str);
/*
  INPUT TO STRINGSET TO GENERATE CODE BELOW:

  CODE_TOO=1 ./tools/stringset/stringset help "csv_stream_help(workspace);" \
  ? "csv_stream_help(workspace);"                                           \
  -? "csv_stream_help(workspace);"                                          \
  --? "csv_stream_help(workspace);"                                         \
  -h "csv_stream_help(workspace);"                                          \
  --help "csv_stream_help(workspace);"                                      \
  edges_init	  "ws->edges_init = atol(values[i]);"                    \
  edges_per_batch "ws->edges_per_batch = atol(values[i]);"               \
  num_batches	  "ws->num_batches = atol(values[i]);"                   \
  file		  "ws->file = values[i]; if(ws->fp) fclose(ws->fp); ws->fp = fopen(ws->file, \"r\");" \
  separator	  "ws->sep = values[i][0];" \
  edge	" \
    splitLineCSVDynamicInPlace('|', values[i], strlen(values[i]), &ws->fields, &ws->lengths, &ws->fieldsSize, &ws->count); \
    csv_emap_t * emap = xmalloc(sizeof(csv_emap_t));                                                                       \
    emap->src_field = atol(ws->fields[0]);                                                                                 \
    emap->src_type = atol(ws->fields[1]);                                                                                  \
    emap->dest_field = atol(ws->fields[2]);                                                                                \
    emap->dest_type = atol(ws->fields[3]);                                                                                 \
    emap->weight_field = atol(ws->fields[4]);                                                                              \
    emap->time_field = atol(ws->fields[5]);                                                                                \
    emap->edge_type = atol(ws->fields[6]);                                                                                 \
    emap->next = ws->list;                                                                                                 \
    ws->list = emap;"                                                                                                      \


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
			csv_stream_help(workspace);
		      }
		    } break;
		  case 'h':
		    {
		      str++; len--;
		      if(!strncmp(str, "elp", len)) {
			str += 3; len -= 3;
			if(len == 0) {
			  /* --help */
			  csv_stream_help(workspace);
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
		  csv_stream_help(workspace);
		}
	      } break;
	    case 'h':
	      {
		str++; len--;
		if(len == 0) {
		  /* -h */
		  csv_stream_help(workspace);
		}
	      } break;
	  }
	} break;
      case '?':
	{
	  str++; len--;
	  if(len == 0) {
	    /* ? */
	    csv_stream_help(workspace);
	  }
	} break;
      case 'e':
	{
	  str++; len--;
	  if(!strncmp(str, "dge", 3)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* edge */
	      splitLineCSVDynamicInPlace('|', values[i], strlen(values[i]), &ws->fields, &ws->lengths, &ws->fieldsSize, &ws->count);     csv_emap_t * emap = xmalloc(sizeof(csv_emap_t));                                                                           emap->src_field = atol(ws->fields[0]);                                                                                     emap->src_type = atol(ws->fields[1]);                                                                                      emap->dest_field = atol(ws->fields[2]);                                                                                    emap->dest_type = atol(ws->fields[3]);                                                                                     emap->weight_field = atol(ws->fields[4]);                                                                                  emap->time_field = atol(ws->fields[5]);                                                                                    emap->edge_type = atol(ws->fields[6]);                                                                                     emap->next = ws->list;                                                                                                     ws->list = emap;
	    }
	    if(!strncmp(str, "s_", 2)) {
	      str += 2; len -= 2;
	      if(len) switch(*str) {
		case 'i':
		  {
		    str++; len--;
		    if(!strncmp(str, "nit", len)) {
		      str += 3; len -= 3;
		      if(len == 0) {
			/* edges_init */
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
	  }
	} break;
      case 'f':
	{
	  str++; len--;
	  if(!strncmp(str, "ile", len)) {
	    str += 3; len -= 3;
	    if(len == 0) {
	      /* file */
	      ws->file = values[i]; if(ws->fp) fclose(ws->fp); ws->fp = fopen(ws->file, "r");
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
	      csv_stream_help(workspace);
	    }
	  }
	} break;
      case 'n':
	{
	  str++; len--;
	  if(!strncmp(str, "um_batches", len)) {
	    str += 10; len -= 10;
	    if(len == 0) {
	      /* num_batches */
	      ws->num_batches = atol(values[i]);
	    }
	  }
	} break;
      case 's':
	{
	  str++; len--;
	  if(!strncmp(str, "eparator", len)) {
	    str += 8; len -= 8;
	    if(len == 0) {
	      /* separator */
	      ws->sep = values[i][0];
	    }
	  }
	} break;
    }

/* END GENERATED CODE */
  }

  free(keys); free(values);
  return STINGER_SUCCESS;
}
