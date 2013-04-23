#include "jsmn_extension.h"

void 
parse_array(jsmntok_t ** tok, string_t * s, string_t * loc, void (*process)(void *, string_t *, char *, int), void * workspace) {
  int remaining = (*tok)->size;
  (*tok)++;
  int len = string_length(loc);
  string_append_cstr_len(loc, "@.", 2);

  /* all values*/
  while(remaining) {
    switch((*tok)->type) {
      case JSMN_PRIMITIVE:
	process(workspace, loc, s->str + (*tok)->start, (*tok)->end - (*tok)->start);
	(*tok)++;
      break;
      case JSMN_OBJECT:
	parse_object(tok, s, loc, process, workspace);
      break;
      case JSMN_ARRAY:
	parse_array(tok, s, loc, process, workspace);
      break;
      case JSMN_STRING:
	process(workspace, loc, s->str + (*tok)->start, (*tok)->end - (*tok)->start);
	(*tok)++;
      break;
    }
    remaining--;
  }
  string_truncate(loc, len);
}

void
parse_object(jsmntok_t ** tok, string_t * s, string_t * loc, void (*process)(void *, string_t *, char *, int), void * workspace) {
  int remaining = (*tok)->size;
  (*tok)++;

  int len = string_length(loc);
  string_append_cstr_len(loc, "$.", 2);

  /* parse keys and values */
  while(remaining) {
    int inner_len = string_length(loc);
    /* key */
    if((*tok)->type == JSMN_STRING) {
      string_append_cstr_len(loc, s->str + (*tok)->start, (*tok)->end - (*tok)->start);
      string_append_cstr_len(loc, ".", 1);
    } else {
      printf("%s %d ERROR!!!\n", __func__, __LINE__);
      printf("token is: type %d size %d start %d end %d string %.*s\n", (*tok)->type, (*tok)->size, (*tok)->start, (*tok)->end,
	(*tok)->end - (*tok)->start, s->str + (*tok)->start);
      exit(-1);
    }
    (*tok)++;
    remaining--;

    /* value */
    switch((*tok)->type) {
      case JSMN_PRIMITIVE:
	process(workspace, loc, s->str + (*tok)->start, (*tok)->end - (*tok)->start);
	(*tok)++;
      break;
      case JSMN_OBJECT:
	parse_object(tok, s, loc, process, workspace);
      break;
      case JSMN_ARRAY:
	parse_array(tok, s, loc, process, workspace);
      break;
      case JSMN_STRING:
	process(workspace, loc, s->str + (*tok)->start, (*tok)->end - (*tok)->start);
	(*tok)++;
      break;
    }
    string_truncate(loc, inner_len);
    remaining--;
  }
  string_truncate(loc, len);
}

void
process_json(string_t * s, string_t * loc, long int * token_count, jsmntok_t ** tokens, void (*process)(void *, string_t *, char *, int), void * workspace) {
  string_append_char(s, '\0');
  string_truncate(loc, 0);

  jsmn_parser p;
  jsmn_init(&p);

  int result = jsmn_parse(&p, s->str, *tokens, *token_count);
  while(result == JSMN_ERROR_NOMEM) {
    (*token_count) *= 2;
    (*tokens) = realloc((*tokens), (*token_count) * sizeof(jsmntok_t));
    result = jsmn_parse(&p, s->str, (*tokens), (*token_count));
  }
  if(JSMN_SUCCESS == result) {
    if((*tokens)[0].type == JSMN_OBJECT && (*tokens)[0].size > 0) {
      jsmntok_t * tok = *tokens;
      parse_object(&tok, s, loc, process, workspace);
    }
  }

  string_truncate(s, 0);
}

void
process_json_stream(FILE * file, int count, void (*process)(void *, string_t *, char *, int), void (*post_process)(void*), void * workspace) {
  string_t * s = string_new();
  string_t * loc = string_new();

  long int token_count = 4096;
  jsmntok_t * tokens = malloc(token_count * sizeof(jsmntok_t));

  int c;
  int in_quotes = 0;
  int bracket_count = 0;
  int object_count = 0;
  if(count > 0) {
    if(post_process) {
      do {
	c = getc(file);

	string_append_char(s, c);

	/* tracking when we hit the end of the outermost object */
	/* we must ignore anything inside of quotes */
	if(c == '"') {
	  in_quotes = !in_quotes;
	} else if(!in_quotes) {
	  if(c == '{') {
	    bracket_count++;
	  } else if (c == '}') {
	    bracket_count--;
	    if(bracket_count == 0 && s->len > 0) {
	      process_json(s, loc, &token_count, &tokens, process, workspace);
	      object_count++;
	      post_process(workspace);
	    }
	  }
	}
      } while (c != EOF && object_count < count);
    } else {
      void ** w = (void **)workspace;
      do {
	c = getc(file);

	string_append_char(s, c);

	/* tracking when we hit the end of the outermost object */
	/* we must ignore anything inside of quotes */
	if(c == '"') {
	  in_quotes = !in_quotes;
	} else if(!in_quotes) {
	  if(c == '{') {
	    bracket_count++;
	  } else if (c == '}') {
	    bracket_count--;
	    if(bracket_count == 0 && s->len > 0) {
	      process_json(s, loc, &token_count, &tokens, process, w[object_count]);
	      object_count++;
	    }
	  }
	}
      } while (c != EOF && object_count < count);
    }
  } else if(count == -1) {
    if(post_process) {
      do {
	c = getc(file);

	string_append_char(s, c);

	/* tracking when we hit the end of the outermost object */
	/* we must ignore anything inside of quotes */
	if(c == '"') {
	  in_quotes = !in_quotes;
	} else if(!in_quotes) {
	  if(c == '{') {
	    bracket_count++;
	  } else if (c == '}') {
	    bracket_count--;
	    if(bracket_count == 0 && s->len > 0) {
	      process_json(s, loc, &token_count, &tokens, process, workspace);
	      post_process(workspace);
	    }
	  }
	}
      } while (c != EOF);
    } else {
      do {
	c = getc(file);

	string_append_char(s, c);

	/* tracking when we hit the end of the outermost object */
	/* we must ignore anything inside of quotes */
	if(c == '"') {
	  in_quotes = !in_quotes;
	} else if(!in_quotes) {
	  if(c == '{') {
	    bracket_count++;
	  } else if (c == '}') {
	    bracket_count--;
	    if(bracket_count == 0 && s->len > 0) {
	      process_json(s, loc, &token_count, &tokens, process, workspace);
	    }
	  }
	}
      } while (c != EOF);
    }
  }
}
