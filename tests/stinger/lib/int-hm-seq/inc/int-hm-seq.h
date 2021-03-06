#ifndef  INT_HM_SEQ_H
#define  INT_HM_SEQ_H

#include <stdint.h>

typedef struct int_hm_seq {
  int64_t size;
  int64_t mask;
  int64_t elements;
  int64_t removed;
  int64_t *keys;
  int64_t *vals;
} int_hm_seq_t;

#define INT_HT_SEQ_EMPTY INT64_MIN
#define INT_HT_SEQ_REMOVED INT64_MAX


int_hm_seq_t *
int_hm_seq_new(int64_t size);

int_hm_seq_t *
int_hm_seq_free(int_hm_seq_t * ht);

int_hm_seq_t * 
int_hm_seq_expand(int_hm_seq_t * ht, int64_t new_size);

int_hm_seq_t * 
int_hm_seq_expand_versioned(int_hm_seq_t * ht, int64_t new_size, int64_t v);

int_hm_seq_t *
int_hm_seq_insert_versioned(int_hm_seq_t * ht, int64_t k, int64_t v);

int_hm_seq_t *
int_hm_seq_insert(int_hm_seq_t * ht, int64_t k, int64_t v);

int64_t *
int_hm_seq_get_location(int_hm_seq_t * ht, int64_t k);

int_hm_seq_t *
int_hm_seq_remove(int_hm_seq_t * ht, int64_t k);

int64_t
int_hm_seq_get(int_hm_seq_t * ht, int64_t k);

#endif  /*INT-HM-SEQ_H*/
