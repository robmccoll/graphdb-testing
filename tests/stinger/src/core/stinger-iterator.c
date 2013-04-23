#include    "stinger-iterator.h"
#include    "stinger-internal.h"
#include    "xmalloc.h"
#include    "stinger.h"

/** @brief Allocate and initialize a new stinger_iterator_t
 *
 * Allocates and initializes a new iterator.  Apply one or more filters
 * to the iterator, then use stinger_iterator_next() in a while loop to
 * advance through each edge in the graph that passes the filter.  You 
 * are responsible for freeing this iterator via stinger_iterator_free().
 *
 * @param s A valid STINGER structure to traverse
 * @return A pointer to a new stinger_iterator_t
 */
stinger_iterator_t *
stinger_iterator_new(struct stinger * s) {
  stinger_iterator_t * iter = xmalloc(sizeof(stinger_iterator_t));
  iter->i.vtx_filter_copy	= 1;
  iter->i.edge_type_filter_copy	= 1;
  iter->i.vtx_type_filter_copy	= 1;
  return stinger_iterator_renew(iter, s);
}

/** @brief Clears filters from and resets an existing iterator.
 *
 * @param iter The iterator to reset
 * @param s A valid STINGER structure to traverse
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_renew(stinger_iterator_t * iter, struct stinger * s) {
  if(!iter->i.vtx_filter_copy && iter->i.vtx_filter)
    free(iter->i.vtx_filter);
  if(!iter->i.vtx_type_filter_copy && iter->i.vtx_type_filter)
    free(iter->i.vtx_type_filter);
  if(!iter->i.vtx_filter_copy && iter->i.vtx_filter)
    free(iter->i.vtx_filter);
  iter->i.s			= s;
  iter->i.flags			= 0;
  iter->i.modified_before	= INT64_MAX;
  iter->i.modified_after	= INT64_MIN; 
  iter->i.created_before	= INT64_MAX; 
  iter->i.created_after		= INT64_MIN; 
  iter->i.vtx_filter_copy	= 1;
  iter->i.edge_type_filter_copy	= 1;
  iter->i.vtx_type_filter_copy	= 1;
  iter->i.active		= 0;
  return iter;
}

/** @brief Tests to see if an iterator is in a valid state.
 *
 * Should return 0 for a valid iterator.
 *
 * @param iter The iterator
 * @param nv Maximum number of vertices
 * @return A bitfield of flags that each indicate an error
 */
int64_t
stinger_iterator_consistency_check(stinger_iterator_t * iter, uint64_t nv) {
  uint64_t  rtn = 0;

  if(!iter->i.s)
    return 0x1;

  /* times */
  if(iter->i.flags & 0x8) {
    if(iter->i.modified_before < iter->i.modified_after)
      rtn &= 0x2;
    if(iter->i.created_before < iter->i.created_after)
      rtn &= 0x4;
  }

  /* vertices */
  if(iter->i.flags & 0x1) {
    if(!iter->i.vtx_filter)
      rtn &= 0x8;
    else if(iter->i.vtx_index >= iter->i.vtx_filter_count)
      rtn &= 0x10;
    else for(uint64_t i = iter->i.vtx_index; i < iter->i.vtx_filter_count; i++)
      if(iter->i.vtx_filter[i] >= nv) 
	rtn &= 0x20;
  }

  /* edge types */
  if(iter->i.flags & 0x2) {
    if(!iter->i.edge_type_filter)
      rtn &= 0x40;
    else if(iter->i.edge_type_index >= iter->i.edge_type_filter_count)
      rtn &= 0x80;
    else for(uint64_t i = iter->i.edge_type_index; i < iter->i.edge_type_filter_count; i++)
      if(iter->i.edge_type_filter[i] >= STINGER_NUMETYPES) 
	rtn &= 0x100;
  }

  /* vtx types */
  if(iter->i.flags & 0x4) 
    if(!iter->i.vtx_type_filter)
      rtn &= 0x200;

  /* predicate */
  if(iter->i.flags & 0x10) 
    if(!iter->i.predicate)
      rtn &= 0x400;

  if(!iter->i.cur_eb)
    rtn &= 0x800;
  else {
    rtn &= 0x1000;
    for(uint64_t i = 0; i < iter->i.s->ETA[iter->i.cur_eb->etype].high; i++)
      if(iter->i.s->ebpool->ebpool + iter->i.s->ETA[iter->i.cur_eb->etype].blocks[i] == iter->i.cur_eb) {
	rtn &= ~0x1000;
        if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  rtn &= 0x2000;
      }
  }

  return rtn;
}

/** @brief Frees an iterator and its internal data.
 * 
 * @param iter Pointer to the iterator
 */
void
stinger_iterator_free(stinger_iterator_t * iter) {
  if(!iter->i.vtx_filter_copy && iter->i.vtx_filter)
    free(iter->i.vtx_filter);
  if(!iter->i.vtx_type_filter_copy && iter->i.vtx_type_filter)
    free(iter->i.vtx_type_filter);
  if(!iter->i.vtx_filter_copy && iter->i.vtx_filter)
    free(iter->i.vtx_filter);
  free(iter);
}

/** @brief Is the iterator in the middle of an iteration.
 *
 * @param iter The iterator in quesiton
 * @return A boolean int (active is 1)
 */
int
stinger_iterator_is_active(stinger_iterator_t * iter) {
  return iter->i.active;
}

/** @brief Deactivate an iterator
 *
 * Effectively resets an iterator to the beginning of its iteration without 
 * clearing any existing filters on the iterator.  Useful if you want to 
 * restart the iteration, break from a current iteration, or change a filter.
 *
 * NOTE: May not work for parallel iterators.
 *
 * @param iter The iterator to be deactivated
 */
void
stinger_iterator_deactivate(stinger_iterator_t * iter) {
  iter->i.active = 0;
}

/** @brief Filter an iterator's traversal to only include a set of vertices and 
 * their neighbors.
 *
 * Should only be applied to an inactive iterator.  Copies the array of 
 * vertexIDs internally.
 * 
 * @param vertices An array of vertices
 * @param count The number of vertices in the array
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_vertex_filter(int64_t * vertices, int64_t count, stinger_iterator_t * iter) {
  if(!iter->i.vtx_filter_copy && iter->i.vtx_filter)
    free(iter->i.vtx_filter);
  iter->i.vtx_filter = xmalloc(count * sizeof(int64_t));
  xelemcpy(iter->i.vtx_filter, vertices, count);
  iter->i.vtx_filter_count = count;
  iter->i.vtx_filter_copy = 0;
  iter->i.flags |= 0x1;
  return iter;
}

/** @brief Filter an iterator's traversal to only include vertices of a certain type 
 * and optionally their neighbors.
 *
 * Should only be applied to an inactive iterator.  Copies the array of 
 * vertex types internally.
 *
 * @param types An array of vertex types
 * @param count The number of vertex types in the array
 * @param both Boolean indicating if both vertices must match the filter (1) or at least one (0).
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_vertex_type_filter(int64_t * types, int64_t count, int both, stinger_iterator_t * iter) {
  if(!iter->i.vtx_type_filter_copy && iter->i.vtx_type_filter)
    free(iter->i.vtx_type_filter);
  iter->i.vtx_type_filter = xmalloc(count * sizeof(int64_t));
  xelemcpy(iter->i.vtx_type_filter, types, count);
  iter->i.vtx_type_filter_count = count;
  iter->i.vtx_type_filter_both = both;
  iter->i.vtx_type_filter_copy = 0;
  iter->i.flags |= 0x4;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges of a certain type 
 *
 * Should only be applied to an inactive iterator.  Copies the array of 
 * edge types internally.
 *
 * @param types An array of edge types
 * @param count The number of edge types in the array
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_edge_type_filter(int64_t * types, int64_t count, stinger_iterator_t * iter) {
  if(!iter->i.edge_type_filter_copy && iter->i.edge_type_filter)
    free(iter->i.edge_type_filter);
  iter->i.edge_type_filter = xmalloc(count * sizeof(int64_t));
  xelemcpy(iter->i.edge_type_filter, types, count);
  iter->i.edge_type_filter_count = count;
  iter->i.edge_type_filter_copy = 0;
  iter->i.flags |= 0x2;
  return iter;
}

/** @brief Filter an iterator's traversal to only include a set of vertices and 
 * their neighbors.
 *
 * Should only be applied to an inactive iterator.  Does not copy the array of 
 * vertexIDs internally.  The input array must remain constant while
 * this iterator is in use.
 * 
 * @param vertices An array of vertices
 * @param count The number of vertices in the array
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_vertex_filter_no_copy(int64_t * vertices, int64_t count, stinger_iterator_t * iter) {
  if(!iter->i.vtx_filter_copy && iter->i.vtx_filter)
    free(iter->i.vtx_filter);
  iter->i.vtx_filter = vertices;
  iter->i.vtx_filter_count = count;
  iter->i.vtx_filter_copy = 1;
  iter->i.flags |= 0x1;
  return iter;
}

/** @brief Filter an iterator's traversal to only include vertices of a certain type 
 * and optionally their neighbors.
 *
 * Should only be applied to an inactive iterator.  Does not copy the array of 
 * vertex types internally.  The input array must remain constant while
 * this iterator is in use.
 *
 * @param types An array of vertex types
 * @param count The number of vertex types in the array
 * @param both Boolean indicating if both vertices must match the filter (1) or at least one (0).
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_vertex_type_filter_no_copy(int64_t * types, int64_t count, int both, stinger_iterator_t * iter) {
  if(!iter->i.vtx_type_filter_copy && iter->i.vtx_type_filter)
    free(iter->i.vtx_type_filter);
  iter->i.vtx_type_filter = types;
  iter->i.vtx_type_filter_count = count;
  iter->i.vtx_type_filter_both = both;
  iter->i.vtx_type_filter_copy = 1;
  iter->i.flags |= 0x4;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges of a certain type 
 *
 * Should only be applied to an inactive iterator.  Does not copy the array of 
 * edge types internally.  The input array must remain constant while
 * this iterator is in use.
 *
 * @param types An array of edge types
 * @param count The number of edge types in the array
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_edge_type_filter_no_copy(int64_t * types, int64_t count, stinger_iterator_t * iter) {
  if(!iter->i.edge_type_filter_copy && iter->i.edge_type_filter)
    free(iter->i.edge_type_filter);
  iter->i.edge_type_filter = types;
  iter->i.edge_type_filter_count = count;
  iter->i.edge_type_filter_copy = 1;
  iter->i.flags |= 0x2;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges with a first time
 * before the specified time.
 *
 * @param time The timestamp
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_created_before(int64_t time, stinger_iterator_t * iter) {
  iter->i.created_before = time;
  iter->i.flags |= 0x8;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges with a first time
 * after the specified time.
 *
 * @param time The timestamp
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_created_after(int64_t time, stinger_iterator_t * iter) {
  iter->i.created_after = time;
  iter->i.flags |= 0x8;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges with a recent time
 * before the specified time.
 *
 * @param time The timestamp
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_modified_before(int64_t time, stinger_iterator_t * iter) {
  iter->i.modified_before = time;
  iter->i.flags |= 0x8;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges with a recent time
 * after the specified time.
 *
 * @param time The timestamp
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_modified_after(int64_t time, stinger_iterator_t * iter) {
  iter->i.modified_after = time;
  iter->i.flags |= 0x8;
  return iter;
}

/** @brief Filter an iterator's traversal to only include edges that pass the 
 * user-specified predicate function.
 *
 * The function pointer will be passed an iterator with an edge assigned to it 
 * to be tested.
 *
 * @param predicate A function pointer that returns non-zero if the edge should be visited
 * @param iter The iterator to be filtered
 * @return A pointer to the same iterator to allow function chaining
 */
stinger_iterator_t *
stinger_iterator_custom_filter(int (*predicate)(stinger_iterator_t *), stinger_iterator_t * iter) {
  iter->i.predicate = predicate;
  iter->i.flags |= 0x10;
  return iter;
}

/********************************************** 
 * helper functions for stinger_iterator_next * 
 **********************************************/
inline int
stinger_iterator_check_time(stinger_iterator_t * iter);
inline void
stinger_iterator_get_edge(stinger_iterator_t * iter);
inline void
stinger_iterator_get_metadata(stinger_iterator_t * iter);
inline int
stinger_iterator_next_block_by_type(stinger_iterator_t * iter);
inline int
stinger_iterator_next_block_by_vtx(stinger_iterator_t * iter);
inline int
stinger_iterator_next_block_all_types(stinger_iterator_t * iter);
inline int
stinger_iterator_check_vtype(stinger_iterator_t * iter);
inline int
stinger_iterator_check_etype(stinger_iterator_t * iter);
inline int
stinger_iterator_check_predicate(stinger_iterator_t * iter);

int
stinger_iterator_check_time(stinger_iterator_t * iter) {
  return (iter->i.cur_eb->edges[iter->i.cur_edge].timeFirst > iter->i.created_after &&
	  iter->i.cur_eb->edges[iter->i.cur_edge].timeFirst < iter->i.created_before &&
	  iter->i.cur_eb->edges[iter->i.cur_edge].timeRecent > iter->i.modified_after &&
	  iter->i.cur_eb->edges[iter->i.cur_edge].timeRecent < iter->i.modified_before);
}

void
stinger_iterator_get_edge(stinger_iterator_t * iter) {
  iter->dest = iter->i.cur_eb->edges[iter->i.cur_edge].neighbor;
  iter->weight = iter->i.cur_eb->edges[iter->i.cur_edge].weight;
  iter->timerecent= iter->i.cur_eb->edges[iter->i.cur_edge].timeRecent;
  iter->timefirst = iter->i.cur_eb->edges[iter->i.cur_edge].timeFirst;
}

void
stinger_iterator_get_metadata(stinger_iterator_t * iter) {
  iter->source = iter->i.cur_eb->vertexID;
  iter->type = iter->i.cur_eb->etype;
}

int
stinger_iterator_next_block_by_type(stinger_iterator_t * iter) {
  int64_t curtype = iter->i.edge_type_filter[iter->i.edge_type_index];
  iter->i.cur_edge = 0;
  while(1) {
    iter->i.edge_block_index++;
    if(iter->i.edge_block_index < iter->i.s->ETA[curtype].high) {
      iter->i.cur_eb = iter->i.s->ebpool->ebpool + iter->i.s->ETA[curtype].blocks[iter->i.edge_block_index];
      stinger_iterator_get_metadata(iter);
      return 0;
    } else {
      iter->i.edge_type_index++;
      if(iter->i.edge_type_index < iter->i.edge_type_filter_count) {
	iter->i.edge_block_index = 0;
	curtype = iter->i.edge_type_filter[iter->i.edge_type_index];
      } else {
	iter->i.active = 0;
	return -1;
      }
    }
  }
}

int
stinger_iterator_next_block_by_vtx(stinger_iterator_t * iter) {
  iter->i.cur_edge = 0;
  if(iter->i.cur_eb && iter->i.cur_eb->next) {
    iter->i.cur_eb = iter->i.s->ebpool->ebpool + iter->i.cur_eb->next;
  } else {
    while(1) {
      iter->i.vtx_index++;
      if(iter->i.vtx_index < iter->i.vtx_filter_count) {
	iter->i.cur_eb = iter->i.s->ebpool->ebpool + stinger_vertex_edges_get(iter->i.s->vertices,iter->i.vtx_filter[iter->i.vtx_index]);
	if(iter->i.cur_eb) {
	  stinger_iterator_get_metadata(iter);
	  return 0;
	}
      } else {
	iter->i.active = 0;
	return -1;
      }
    }
  }
  return 0;
}

int
stinger_iterator_next_block_all_types(stinger_iterator_t * iter) {
  int64_t curtype = iter->i.edge_type_index;
  iter->i.cur_edge = 0;
  while(1) {
    iter->i.edge_block_index++;
    if(iter->i.edge_block_index < iter->i.s->ETA[curtype].high) {
      iter->i.cur_eb = iter->i.s->ebpool->ebpool + iter->i.s->ETA[curtype].blocks[iter->i.edge_block_index];
      stinger_iterator_get_metadata(iter);
      return 0;
    } else {
      iter->i.edge_type_index++;
      curtype = iter->i.edge_type_index;
      if(iter->i.edge_type_index < STINGER_NUMETYPES) {
	iter->i.edge_block_index = 0;
      } else {
	iter->i.active = 0;
	return -1;
      }
    }
  }
}

int
stinger_iterator_check_vtype(stinger_iterator_t * iter) {
  int64_t src_type = stinger_vtype_get(iter->i.s, iter->i.cur_eb->vertexID);
  int64_t dst_type = stinger_vtype_get(iter->i.s, iter->i.cur_eb->edges[iter->i.cur_edge].neighbor);
  if(iter->i.vtx_type_filter_both) {
    int found_src = 0;
    int found_dst = 0;
    for(uint64_t i = 0; i < iter->i.vtx_type_filter_count; i++) {
      if(src_type == iter->i.vtx_type_filter[i])
	found_src = 1;
      if(dst_type == iter->i.vtx_type_filter[i])
	found_dst = 1;
      if(found_src && found_dst)
	return 1;
    }
    return 0;
  } else {
    for(uint64_t i = 0; i < iter->i.vtx_type_filter_count; i++) {
      if(src_type == iter->i.vtx_type_filter[i] || dst_type == iter->i.vtx_type_filter[i])
	return 1;
    }
    return 0;
  }
}

int
stinger_iterator_check_etype(stinger_iterator_t * iter) {
  int64_t type = iter->i.cur_eb->etype;
  for(uint64_t i = 0; i < iter->i.edge_type_filter_count; i++) {
    if(type == iter->i.edge_type_filter[i])
      return 1;
  }
  return 0;
}

int
stinger_iterator_check_predicate(stinger_iterator_t * iter) {
  stinger_iterator_get_edge(iter);
  return iter->i.predicate(iter);
}

/** @brief Advance an iterator to the next edge that matches the internal filter.
 *
 * This function will advance the iterator through the STINGER data in the most
 * efficient manner possible based on the filters.  It is intended to be used in
 * a serial context.  No particular ordering is guaranteed, but the iterator 
 * will visit all edges that pass the filter exactly one.
 *
 * \code
 * while(stinger_iterator_nex(iter)) {
 *     // do something with the current edge
 * }
 * \endcode
 *
 * @param iter The iterator to be advanced
 * @return A boolean int indicating more edges to be visited.
 */
int
stinger_iterator_next(stinger_iterator_t * iter) {
  switch(iter->i.flags) { 
    default: { 
      fprintf(stderr, "WARNING: No iterator filter provided\n"); 
      return 0;
    } break; 

    /* VTX */  
    case 1: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* ETYPE */ 
    case 2: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* VTX + etype */ 
    case 3: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_etype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* VTYPE */ 
    case 4: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* VTX + vtype */ 
    case 5: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* ETYPE + vtype */ 
    case 6: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* vtx + etype + vtype */ 
    case 7: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_etype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time */
    case 8: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time + vtx */
    case 9: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time + ETYPE */ 
    case 10: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time + vtx + etype */ 
    case 11: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_etype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 
    
    /* time + vtype */ 
    case 12: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time + vtx + vtype */ 
    case 13: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time + ETYPE + vtype */ 
    case 14: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* time+ vtx + etype + vtype */ 
    case 15: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_etype(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate */ 
    case 16: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + vtx */  
    case 17: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + ETYPE */ 
    case 18: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + VTX + etype */ 
    case 19: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_etype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + VTYPE */ 
    case 20: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + VTX + vtype */ 
    case 21: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + ETYPE + vtype */ 
    case 22: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + vtx + etype + vtype */ 
    case 23: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_etype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time */
    case 24: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time + vtx */
    case 25: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time + ETYPE */ 
    case 26: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time + vtx + etype */ 
    case 27: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_etype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 
    
    /* predicate + time + vtype */ 
    case 28: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_all_types(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_all_types(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time + vtx + vtype */ 
    case 29: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time + ETYPE + vtype */ 
    case 30: { 
      if(!iter->i.active) {
	iter->i.edge_type_index = 0;
	iter->i.edge_block_index = -1;
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	if(stinger_iterator_next_block_by_type(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_type(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 

    /* predicate + time+ vtx + etype + vtype */ 
    case 31: { 
      if(!iter->i.active) {
	iter->i.cur_edge = 0;
	iter->i.active = 1;
	iter->i.vtx_index = -1;
	iter->i.cur_eb = NULL;
	if(stinger_iterator_next_block_by_vtx(iter))
	  return 0;
      } else {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high)
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
      }
      while(stinger_eb_is_blank(iter->i.cur_eb, iter->i.cur_edge) || !stinger_iterator_check_time(iter) || !stinger_iterator_check_vtype(iter) || !stinger_iterator_check_etype(iter) || !stinger_iterator_check_predicate(iter)) {
	iter->i.cur_edge++;
	if(iter->i.cur_edge >= iter->i.cur_eb->high) {
	  if(stinger_iterator_next_block_by_vtx(iter))
	    return 0;
	}
      }
      stinger_iterator_get_edge(iter);
    } break; 
  } 
  return iter->i.active;
}

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
stinger_iterator_loop_end(stinger_iterator_t * iter) {
  //TODO This will enable parallel for loops over blocks through this iterator
  fprintf(stderr, "ERROR: %s not implemented yet\n", __func__);
  return 0;
}

int
stinger_iterator_loop_jump(stinger_iterator_t * iter, int64_t i) {
  //TODO This will enable parallel for loops over blocks through this iterator
  fprintf(stderr, "ERROR: %s not implemented yet\n", __func__);
  return 0;
}

int
stinger_iterator_loop_next(stinger_iterator_t * iter) {
  //TODO This will enable parallel for loops over blocks through this iterator
  fprintf(stderr, "ERROR: %s not implemented yet\n", __func__);
  return 0;
}
