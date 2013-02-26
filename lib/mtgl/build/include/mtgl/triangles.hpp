/*  _________________________________________________________________________
 *
 *  MTGL: The MultiThreaded Graph Library
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top MTGL directory.
 *  _________________________________________________________________________
 */

/****************************************************************************/
/*! \file triangles.hpp

    \author Jon Berry (jberry@sandia.gov)

    \description: Cohen's triangle enumeration algorithm for simple graphs.
                  This doesn't work on multigraphs.

    \date 11/2006
*/
/****************************************************************************/

#ifndef MTGL_TRIANGLES_HPP
#define MTGL_TRIANGLES_HPP

#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cassert>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/partitioning.hpp>

#define MY_BLOCK_SIZE 1024
#define LOAD_BALANCED

namespace mtgl {

template <typename Graph>
class default_triangles_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;

  default_triangles_visitor() {}
  void operator()(size_type v1, size_type v2, size_type v3) {}

  void accumulate() {}
};

namespace detail {

template <typename size_type>
class tri_edge_struct {
public:
  tri_edge_struct(uint64_t k = 0, size_type et = 0) : key(k), eid(et) {}
  uint64_t key;    // s * order + t   : a hash/search key
  size_type eid;   // mtgl edge id: good for property map lookups
};

template <typename size_type>
static inline int e_cmp(const void* u, const void* v)
{
  uint64_t ukey = ((tri_edge_struct<size_type>*)u)->key;
  uint64_t vkey = ((tri_edge_struct<size_type>*)v)->key;
  return static_cast<int>(ukey - vkey);
}

template <typename size_type>
static inline
size_type low_endpoint_winner(size_type* deg, size_type src, size_type dest)
{
  size_type winner = 0;

  if (deg[src] < deg[dest])
  {
    winner = src;
  }
  else if (deg[src] > deg[dest])
  {
    winner = dest;
  }
  else if (src < dest)
  {
    winner = src;
  }
  else
  {
    winner = dest;
  }

  return winner;
}

#ifdef USING_QTHREADS
template <typename Graph, typename Visitor>
class triangles_loop_functor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef size_t qt_size_t;

  triangles_loop_functor(Graph& gg, Visitor& myv,
                         size_type* dg,
                         size_type* accum_wrk,
                         size_type* strt,
                         size_type num_blcks,
                         size_type ordr,
                         size_type* my_edgs,
                         tri_edge_struct<size_type>* ekys
#ifdef DEBUG
                         , size_type& tr
#endif
                         ) :
    g(gg), tri_visitor(myv), deg(dg), accum_work(accum_wrk),
    start(strt), num_blocks(num_blcks), order(ordr), my_edges(my_edgs),
    ekeys(ekys)
#ifdef DEBUG
    , tri(tr)
#endif
  {
    vid_map = get(_vertex_id_map, g);
    edgs = edges(g);
  }

private:
  Graph& g;
  Visitor& tri_visitor;
  size_type* deg;
  size_type* accum_work;
  size_type* start;
  size_type num_blocks;
  size_type order;
  size_type* my_edges;
  tri_edge_struct<size_type>* ekeys;
#ifdef DEBUG
  size_type& tri;
#endif
  vertex_id_map<Graph> vid_map;
  edge_iterator edgs;

public:
  inline void operator()(qt_size_t start_at, qt_size_t stop_at)
  {
    for (qt_size_t stream_id = start_at; stream_id < stop_at; ++stream_id)
    {
      Visitor my_visitor = tri_visitor;

#ifdef DEBUG
      size_type my_tri = 0;
#endif

      size_type start_pos = begin_block_range(accum_work[order],
                                              stream_id, num_blocks);
      size_type end_pos = end_block_range(accum_work[order],
                                          stream_id, num_blocks);

      // *********************************************************************
      // *********************************************************************
      // start_pos is an index into total work.  We need to find the index of
      // the vertex that owns this work and the index of the vertex just after
      // the end of the work.
      // *********************************************************************
      // *********************************************************************
      size_type start_outer =
        begin_manhattan_outer_range(accum_work, order, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_work, order, end_pos);

      // *********************************************************************
      // *********************************************************************
      // Now, we know which vertex or vertices own the work assocted with this
      // stream.  Iterate over them. Note that this stream might be responsible
      // for only a fraction of the work associated with one vertex.  For the
      // first vertex, find the offset into its work for the portion of work
      // that this stream will do.
      // *********************************************************************
      // *********************************************************************
      for (size_type i = start_outer; i != end_outer; ++i)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_work, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_work, end_pos, end_outer, i);

        // *********************************************************************
        // *********************************************************************
        // Now we are actually going to do some work.  start_inner is the index
        // into the work associated with vertex i.  We decode this to find the
        // pair of edges out of the bucket to test for triangle-ness.  There is
        // w = k(k-1)/2 work for a vertex with bucket size k.
        //
        // We unrank position w via the following code, and the result is the
        // pair (index1, index2), which we use to select two edges out of
        // vertex i's bucket to test for triangle-ness.
        // *********************************************************************
        // *********************************************************************
        size_type k = my_edges[i];
        for (size_type j = start_inner; j != end_inner; ++j)
        {
          size_type x = 0;
          size_type p = 1;

          while (x <= j)
          {
            x += (k - p);
            p++;
          }

          size_type index1 = p - 2;
          size_type oldx = x - (k - index1);
          size_type index2 = index1 + (j - oldx);

          // *********** (index1, index2) specifies our pair *******************

          // *******************************************************************
          // *******************************************************************
          // Retrieve the edge corresponding to index1.
          // *******************************************************************
          // *******************************************************************
          size_type eid1 = ekeys[start[i] + index1].eid;
          size_type src1 = get(vid_map, source(edgs[eid1], g));
          size_type dest1 = get(vid_map, target(edgs[eid1], g));

          // *******************************************************************
          // *******************************************************************
          // Retrieve the edge corresponding to index2.
          // *******************************************************************
          // *******************************************************************
          size_type eid2 = ekeys[start[i] + index2].eid;
          size_type src2 = get(vid_map, source(edgs[eid2], g));
          size_type dest2 = get(vid_map, target(edgs[eid2], g));

          // *******************************************************************
          // *******************************************************************
          // Test the triangle.
          // *******************************************************************
          // *******************************************************************
          size_type s = (src1 == i) ? dest1 : src1;
          size_type t = (src2 == i) ? dest2 : src2;
          order_pair(s, t);
          tri_edge_struct<size_type> ekey(s * order + t, 0);
          size_type v = low_endpoint_winner(deg, s, t);

          tri_edge_struct<size_type>* eresult = (tri_edge_struct<size_type>*)
            bsearch((void*) &ekey, (void*) &ekeys[start[v]], my_edges[v],
                    sizeof(tri_edge_struct<size_type>), e_cmp<size_type>);

          if (eresult != 0)
          {
            my_visitor(eid1, eid2, eresult->eid);

#ifdef DEBUG
            ++my_tri;
#endif
          }
        } // Done with our portion of work associated with vertex i.
      }

      my_visitor.accumulate();

#ifdef DEBUG
      mt_incr(tri, my_tri);
#endif
    } // Done with for loop over streams
  } // operator()
};
#endif

}

/*! \brief Counts the number of triangles in a graph.

    \author Jon Berry (jberry@sandia.gov)

    \date 11/2006

    The algorithm is the one given by Cohen in "Graph Twiddling in a
    MapReduce World".
*/
template <typename Graph, typename Visitor>
void find_triangles(Graph& g, Visitor& tri_visitor)
{
  #pragma mta noalias g
  #pragma mta noalias tri_visitor

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  if (!is_undirected(g) && !is_bidirectional(g))
  {
    fprintf(stderr, "error: find_triangles() runs on undirected or "
            "bidirectional graphs.  Please use a transpose_adapter to apply "
            "this algorithm to directed graphs\n");
    exit(290);
  }

  size_type order = num_vertices(g);
  size_type size = num_edges(g) * 2;

  //size_type threshold = order;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> eid_map = get(_edge_id_map, g);

#ifdef DEBUG
  if (order < 100) print(g);
#endif

  size_type* deg = (size_type*) malloc(sizeof(size_type) * order);
  size_type* next = (size_type*) malloc(sizeof(size_type) * order);
  size_type* my_edges = (size_type*) malloc(sizeof(size_type) * order);

  vertex_iterator verts = vertices(g);

  for (size_type i = 0; i < order; i++)
  {
    deg[i] = out_degree(verts[i], g);
#ifdef DEBUG
  printf("out_degree[%lu]: %lu\n", verts[i], deg[i]);
#endif
    next[i] = 0;
    my_edges[i] = 0;
  }

  size_type* prefix_counts = 0;
  size_type* started_nodes = 0;

  size_type m = num_edges(g);
  edge_iterator edgs = edges(g);
#ifdef DEBUG
  printf("m: %lu\n", m);
#endif
  for (size_type i = 0; i < m; i++)
  {
    edge_descriptor e = edgs[i];
    size_type src_id = get(vid_map, source(e, g));
    size_type dest_id = get(vid_map, target(e, g));
//    size_type dsrc = deg[src_id];
//    size_type ddst = deg[dest_id];

    // We assume we'll see each undir. edge as two dir. edg.
//    if ((dsrc > threshold) || (ddst > threshold)) continue;

    size_type smaller_deg_v = detail::low_endpoint_winner(deg, src_id, dest_id);

    mt_incr(my_edges[smaller_deg_v], 1);
  }

  size_type* start = (size_type*) malloc(sizeof(size_type*) * (order + 1));
  start[0] = 0;

  for (size_type i = 1; i <= order; i++)
  {
    start[i] = start[i - 1] + my_edges[i - 1];
  }

#ifdef DEBUG
  size_type max_bucket_size = 0;
  size_type max_degree = 0;

  #pragma mta assert nodep
  for (size_type i = 0; i < order; i++)
  {
    if (my_edges[i] > max_bucket_size) max_bucket_size = my_edges[i];
  }

  #pragma mta assert nodep
  for (size_type i = 0; i < order; i++)
  {
    if (deg[i] > max_degree) max_degree = deg[i];
  }

  printf("max b: %lu\n", max_bucket_size);
  printf("max d: %lu\n", max_degree);
  fflush(stdout);
#endif

  detail::tri_edge_struct<size_type>* ekeys =
    (detail::tri_edge_struct<size_type>*)
      malloc(sizeof(detail::tri_edge_struct<size_type>) * size);

  #pragma mta assert noalias *ekeys, *start
  #pragma mta assert nodep
  for (size_type i = 0; i < m; i++)
  {
    edge_descriptor e = edgs[i];
    size_type eid = get(eid_map, e);
    size_type src_id = get(vid_map, source(e, g));
    size_type dest_id = get(vid_map, target(e, g));
//    size_type dsrc = deg[src_id];
//    size_type ddst = deg[dest_id];

#ifdef DEBUG
    printf("bucketing %lu (%lu, %lu) deg: (%lu, %lu)\n", eid, src_id,dest_id,
           deg[src_id], deg[dest_id]);
#endif

//    if ((dsrc > threshold) || (ddst > threshold)) continue;

    size_type smaller_deg_v = detail::low_endpoint_winner(deg, src_id, dest_id);
    order_pair(src_id, dest_id);
    uint64_t ekey = src_id * order + dest_id;
    size_type pos = mt_incr(next[smaller_deg_v], 1);
    ekeys[start[smaller_deg_v] + pos] =
      detail::tri_edge_struct<size_type>(ekey, eid);
#ifdef DEBUG
    printf("bucket[%lu] add %lu (%lu, %lu) (key: %lu)\n", smaller_deg_v,
           eid, src_id, dest_id, ekey);
#endif
  }

  free(prefix_counts);
  free(started_nodes);

  #pragma mta assert parallel
  for (size_type i = 0; i < order; i++)
  {
    qsort((void*) &ekeys[start[i]], my_edges[i],
          sizeof(detail::tri_edge_struct<size_type>), detail::e_cmp<size_type>);
  }

#ifdef LOAD_BALANCED
  // *************************************************************************
  // *************************************************************************
  // At this point, the buckets of Cohen's algorithm are all constructed.
  // my_edges[i] is the size of vertex i's bucket.  To explore all potential
  // triangles incident on vertex i, we need to look at all pairs of edges
  // in its bucket.  Up until mtgl v3283, this was done with a parallel loop
  // over all vertices.  However, even with the load balanced bucket strategy,
  // there can be XMT load balancing issues.
  //
  // We'll use a load balancing strategy that eventually will be available
  // via "thread_iterator"s.  Since these were not yet checked into the
  // trunk at the time of this writing, the strategy will be implemented
  // manually.
  //
  // First, compute the total work, which is the sum of k(k-1)/2, where
  // the bucket size of the next vertex is k.
  // *************************************************************************
  // *************************************************************************
  size_type* work = (size_type*) malloc(order * sizeof(size_type));
  size_type* accum_work = (size_type*) malloc((order + 1) * sizeof(size_type));

  for (size_type i = 0; i < order; ++i)
  {
    work[i] = my_edges[i] * (my_edges[i] - 1) / 2;
  }

  accum_work[0] = 0;
  for (size_type i = 0; i < order; ++i)
  {
    accum_work[i + 1] = accum_work[i] + work[i];
  }

  free(work);

  // *************************************************************************
  // *************************************************************************
  // Now figure out how many blocks (work units) to use.
  // *************************************************************************
  // *************************************************************************
#if defined(USING_QTHREADS)
  size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
  size_type num_blocks = (accum_work[order] + MY_BLOCK_SIZE - 1) /
                         MY_BLOCK_SIZE;
#else
  size_type num_blocks = 1;
#endif

  // *************************************************************************
  // *************************************************************************
  // Now begins the load balanced loop.
  // *************************************************************************
  // *************************************************************************

#ifdef DEBUG
  size_type tri = 0;
#endif

#ifndef USING_QTHREADS

  #pragma mta assert parallel
  #pragma mta assert noalias *ekeys, *my_edges, *start
  #pragma mta no scalar expansion
  for (size_type stream_id = 0; stream_id < num_blocks; ++stream_id)
  {
#ifdef DEBUG
    size_type my_tri = 0;
#endif

    Visitor my_visitor = tri_visitor;
    size_type start_pos = begin_block_range(accum_work[order],
                                            stream_id, num_blocks);
    size_type end_pos = end_block_range(accum_work[order],
                                        stream_id, num_blocks);
#ifdef DEBUG
    printf("block %lu: start: %lu, end: %lu\n", stream_id, start_pos,end_pos);
#endif

    // ***********************************************************************
    // ***********************************************************************
    // start_pos is an index into total work.  We need to find the index of
    // the vertex that owns this work and the index of the vertex just after
    // the end of the work.
    // ***********************************************************************
    // ***********************************************************************
    size_type start_outer =
      begin_manhattan_outer_range(accum_work, order, start_pos);
    size_type end_outer = end_manhattan_outer_range(accum_work, order, end_pos);
#ifdef DEBUG
    printf("block %lu: start_outer: %lu, end_outer: %lu\n", stream_id,
           start_outer,end_outer);
#endif

    // ***********************************************************************
    // ***********************************************************************
    // Now, we know which vertex or vertices own the work associated with this
    // stream.  Iterate over them. Note that this stream might be responsible
    // for only a fraction of the work associated with one vertex.  For the
    // first vertex, find the offset into its work for the portion of work
    // that this stream will do.
    // ***********************************************************************
    // ***********************************************************************
    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_work, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_work, end_pos, end_outer, i);
#ifdef DEBUG
    printf("block %lu: outer_iter: %lu, start_inner: %lu, end_inner: %lu\n",
           stream_id, i, start_inner,end_inner);
#endif

      // *********************************************************************
      // *********************************************************************
      // Now we are actually going to do some work.  start_inner is the index
      // into the work associated with vertex i.  We decode this to find the
      // pair of edges out of the bucket to test for triangle-ness.  There is
      // w = k(k-1)/2 work for a vertex with bucket size k.
      //
      // We unrank position w via the following code, and the result is the
      // pair (index1, index2), which we use to select two edges out of
      // vertex i's bucket to test for triangle-ness.
      // *********************************************************************
      // *********************************************************************
      size_type k = my_edges[i];

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type x = 0;
        size_type p = 1;

        while (x <= j)
        {
          x += (k - p);
          p++;
        }

        size_type index1 = p - 2;
        size_type oldx = x - (k - index1);
        size_type index2 = index1 + (j - oldx);

        // *********** (index1, index2) specifies our pair *******************

        // *******************************************************************
        // *******************************************************************
        // Retrieve the edge corresponding to index1.
        // *******************************************************************
        // *******************************************************************
        size_type eid1 = ekeys[start[i] + index1].eid;
        size_type src1 = get(vid_map, source(edgs[eid1], g));
        size_type dest1 = get(vid_map, target(edgs[eid1], g));
#ifdef DEBUG
    printf("block %lu: outer_iter: %lu, inner_iter: %lu, edge at index1: %lu (%lu, %lu)\n", stream_id, i, j, eid1, src1, dest1);
#endif

        // *******************************************************************
        // *******************************************************************
        // Retrieve the edge corresponding to index2.
        // *******************************************************************
        // *******************************************************************
        size_type eid2 = ekeys[start[i] + index2].eid;
        size_type src2 = get(vid_map, source(edgs[eid2], g));
        size_type dest2 = get(vid_map, target(edgs[eid2], g));
#ifdef DEBUG
    printf("block %lu: outer_iter: %lu, inner_iter: %lu, edge at index2: %lu (%lu, %lu)\n", stream_id, i, j, eid2, src2, dest2);
#endif

        // *******************************************************************
        // *******************************************************************
        // Test the triangle.
        // *******************************************************************
        // *******************************************************************
        size_type s = (src1 == i) ? dest1 : src1;
        size_type t = (src2 == i) ? dest2 : src2;
        order_pair(s, t);

        detail::tri_edge_struct<size_type> ekey(s * order + t, 0);
        size_type v = detail::low_endpoint_winner(deg, s, t);
#ifdef DEBUG
        printf("searching for (%lu, %lu) in %lu's bucket\n", s, t, v);
#endif

        detail::tri_edge_struct<size_type>* eresult =
          (detail::tri_edge_struct<size_type>*)
          bsearch((void*) &ekey, (void*) &ekeys[start[v]], my_edges[v],
                  sizeof(detail::tri_edge_struct<size_type>),
                  detail::e_cmp<size_type>);

        if (eresult != 0)
        {
          my_visitor(eid1, eid2, eresult->eid);

#ifdef DEBUG
          ++my_tri;
#endif
        }
      }  // Done with our portion of work associated with vertex i.
    }

    my_visitor.accumulate();

#ifdef DEBUG
    mt_incr(tri, my_tri);
#endif
  }
#else
  detail::triangles_loop_functor<Graph, Visitor>
    tlf(g, tri_visitor, deg, accum_work, start, num_blocks,
        order, my_edges, ekeys
#ifdef DEBUG
        , tri
#endif
        );
  qt_loop_balance(0, num_blocks, tlf);
#endif

  free(accum_work);

  // *************************************************************************
  // *************************************************************************
  // This is the original (non-load-balanced) way it was.
  // *************************************************************************
  // *************************************************************************
#else
  size_type this_stream = 0;
  size_type num_streams = 1;

  #pragma mta for all streams this_stream of num_streams
  // for (size_type i = 0; i < order; i++)
  {
    size_type beg = begin_block_range(order, this_stream, num_streams);
    size_type end = end_block_range(order, this_stream, num_streams);

    for (size_type i = beg; i < end; i++)
    {
      Visitor my_visitor = tri_visitor;

      for (size_type j = 0; j < my_edges[i]; j++)
      {
        size_type eid1 = ekeys[start[i] + j].eid;
        size_type src1 = get(vid_map, source(edgs[eid1], g));
        size_type dest1 = get(vid_map, target(edgs[eid1], g));

        for (size_type k = j + 1; k < my_edges[i]; k++)
        {
          size_type eid2 = ekeys[start[i] + k].eid;
          size_type src2 = get(vid_map, source(edgs[eid2], g));
          size_type dest2 = get(vid_map, target(edgs[eid2], g));

          size_type s = (src1 == i) ? dest1 : src1;
          size_type t = (src2 == i) ? dest2 : src2;
          order_pair(s, t);

          detail::tri_edge_struct<size_type> ekey(s * order + t, 0);
          size_type v = detail::low_endpoint_winner(deg, s, t);

          detail::tri_edge_struct<size_type>* eresult =
            (detail::tri_edge_struct<size_type>*)
            bsearch((void*) &ekey, (void*) &ekeys[start[v]], my_edges[v],
                    sizeof(detail::tri_edge_struct<size_type>),
                    detail::e_cmp<size_type>);

          if (eresult != 0)
          {
            my_visitor(eid1, eid2, eresult->eid);

#ifdef DEBUG
            mt_incr(tri, 1);
#endif
          }
        }
      }
    }
  }
#endif // LOAD_BALANCED

#ifdef DEBUG
  #pragma mta trace "done"
  printf("num tri: %lu\n", tri);
  fflush(stdout);
#endif

  free(deg);
  free(start);
  free(next);
  free(my_edges);
  free(ekeys);
}

}

#undef MY_BLOCK_SIZE
#undef LOAD_BALANCED

#endif
