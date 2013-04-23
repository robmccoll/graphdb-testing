#include "jsmn.h"
#include "astring.h"
#include <stdio.h>
#include <stdlib.h>

void
parse_object(jsmntok_t ** tok, string_t * s, string_t * loc, void (*process)(void *, string_t *, char *, int), void * workspace);

void
parse_array(jsmntok_t ** tok, string_t * s, string_t * loc, void (*process)(void *, string_t *, char *, int), void * workspace);

void
process_json(string_t * s, string_t * loc, long int * token_count, jsmntok_t ** tokens, void (*process)(void *, string_t *, char *, int), void * workspace);

void
process_json_stream(FILE * file, int count, void (*process)(void *, string_t *, char *, int), void (*post_process)(void*), void * workspace);
