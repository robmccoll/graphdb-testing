#ifndef STINGER_ITERATOR_H 
#define STINGER_ITERATOR_H 

#include    "stinger.h"

struct stinger_iterator_internal;

typedef struct stinger_iterator {
  struct stinger_iterator_internal i;
  /* Accessible data for current edge */
  int64_t source;
  int64_t dest;
  int64_t weight;
  int64_t type;
  int64_t timefirst;
  int64_t timerecent;
} stinger_iterator_t;

stinger_iterator_t *
stinger_iterator_new(struct stinger * s);

stinger_iterator_t *
stinger_iterator_renew(stinger_iterator_t * iter, struct stinger * s);

void
stinger_iterator_free(stinger_iterator_t * iter);

int64_t
stinger_iterator_consistency_check(stinger_iterator_t * iter, uint64_t nv);

int
stinger_iterator_is_active(stinger_iterator_t * iter);

void
stinger_iterator_deactivate(stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_vertex_filter(int64_t * vertices, int64_t count, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_vertex_type_filter(int64_t * types, int64_t count, int both, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_edge_type_filter(int64_t * types, int64_t count, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_vertex_filter_no_copy(int64_t * vertices, int64_t count, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_vertex_type_filter_no_copy(int64_t * types, int64_t count, int both, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_edge_type_filter_no_copy(int64_t * types, int64_t count, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_created_before(int64_t time, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_created_after(int64_t time, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_modified_before(int64_t time, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_modified_after(int64_t time, stinger_iterator_t * iter);

stinger_iterator_t *
stinger_iterator_custom_filter(int (*predicate)(stinger_iterator_t *), stinger_iterator_t * iter);

int
stinger_iterator_next(stinger_iterator_t * iter);

/*
 * IDEA These functions will enable some level of parallelism via this iterator
 * using for loops like:
 * int64_t end = stinger_iterator_loop_end(iter);
 * parallel_for(i = 0; i < end; i++) {
 *   stinger_iterator_loop_jump(iter, i);
 *   while(stinger_iterator_loop_next(iter, i) {
 *     // do something
 *   }
 * }
 */
int64_t
stinger_iterator_loop_end(stinger_iterator_t * iter);

int
stinger_iterator_loop_jump(stinger_iterator_t * iter, int64_t i);

int
stinger_iterator_loop_next(stinger_iterator_t * iter);

#endif  /*STINGER-ITERATOR_H*/
