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
/*! \file badrank.hpp

    \author Greg Mackey (gemacke@sandia.gov)

    \date 11/22/2011
*/
/****************************************************************************/

#ifndef MTGL_BADRANK_HPP
#define MTGL_BADRANK_HPP

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/partitioning.hpp>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

#define BR_CHUNK 256
#define PUB_ALG

namespace mtgl {

namespace detail {

#ifdef USING_QTHREADS
template <typename Graph, typename RankMap, typename InitMap, typename AccMap>
class init_propmaps {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  init_propmaps(Graph& gg, RankMap& r, InitMap& i, AccMap& a) :
    g(gg), rank(r), init(i), acc(a) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;
      rank[v] = 0;
      init[v] = 0.0;
      acc[v] = 0.0;
    }
  }

private:
  Graph& g;
  RankMap& rank;
  InitMap& init;
  AccMap& acc;
};

template <typename Graph, typename RankMap, typename AccMap, typename DegreeMap>
class compute_acc_qt {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  compute_acc_qt(Graph& gg, RankMap& r, AccMap& a,
                 DegreeMap& id, size_type* ad, size_type nb) :
    g(gg), rank(r), acc(a), in_degrees(id), accum_deg(ad), num_blocks(nb) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type order = num_vertices(g);

    for (size_t block_id = start; block_id < stop; ++block_id)
    {
      size_type start_pos =
        begin_block_range(accum_deg[order], block_id, num_blocks);
      size_type end_pos =
        end_block_range(accum_deg[order], block_id, num_blocks);

      size_type start_outer =
        begin_manhattan_outer_range(accum_deg, order, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_deg, order, end_pos);

      thread_vertex_iterator tverts = thread_vertices(start_outer, g);

      for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

        vertex_descriptor u = *tverts;
        thread_adjacency_iterator tadj_verts =
          thread_adjacent_vertices(u, start_inner, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor v = *tadj_verts;
          mt_incr(acc[u], rank[v] / in_degrees[v]);
        }
      }
    }
  }

private:
  Graph& g;
  RankMap& rank;
  AccMap& acc;
  DegreeMap& in_degrees;
  size_type* accum_deg;
  size_type num_blocks;
};

template <typename Graph, typename RankMap, typename InitMap, typename AccMap>
class badrank_body_qt {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  badrank_body_qt(Graph& gg, RankMap& r, InitMap& i, AccMap& a, double d,
                  double* m) :
    g(gg), rank(r), init(i), acc(a), dampen(d), maxdiff_vec(m) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;

      double oldval = rank[v];
      double newval;

#ifdef PUB_ALG
      rank[v] = init[v] + dampen * acc[v];
      newval = rank[v];
#else
      if (out_degree(v, g) == 0)
      {
        newval = init[v];
      }
      else
      {
        rank[v] = init[v] + dampen * acc[v] / out_degree(v, g);
        newval = rank[v];
      }
#endif

      maxdiff_vec[start_pos] =
        (oldval > newval) ? oldval - newval : newval - oldval;

      acc[v] = 0.0;
    }
  }

private:
  Graph& g;
  RankMap& rank;
  InitMap& init;
  AccMap& acc;
  double dampen;
  double* maxdiff_vec;
};

template <typename Graph, typename DegreeMap>
class init_in_degrees_qt {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  init_in_degrees_qt(Graph& gg, DegreeMap& id) : g(gg), in_degrees(id) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      in_degrees[*tverts] = 0;
    }
  }

private:
  Graph& g;
  DegreeMap& in_degrees;
};

template <typename Graph, typename DegreeMap>
class compute_in_degrees_qt {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  compute_in_degrees_qt(Graph& gg, DegreeMap& id, size_type* ad, size_type nb) :
    g(gg), in_degrees(id), accum_deg(ad), num_blocks(nb) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type order = num_vertices(g);

    for (size_t block_id = start; block_id < stop; ++block_id)
    {
      size_type start_pos =
        begin_block_range(accum_deg[order], block_id, num_blocks);
      size_type end_pos =
        end_block_range(accum_deg[order], block_id, num_blocks);

      size_type start_outer =
        begin_manhattan_outer_range(accum_deg, order, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_deg, order, end_pos);

      thread_vertex_iterator tverts = thread_vertices(start_outer, g);

      for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

        thread_adjacency_iterator tadj_verts =
          thread_adjacent_vertices(*tverts, start_inner, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          mt_incr(in_degrees[*tadj_verts], 1);
        }
      }
    }
  }

private:
  Graph& g;
  DegreeMap& in_degrees;
  size_type* accum_deg;
  size_type num_blocks;
};
#endif

template <typename Graph, typename RankMap, typename AccMap, typename DegreeMap>
void compute_acc(Graph& g, RankMap& rank, AccMap& acc, DegreeMap& in_degrees,
                 typename graph_traits<Graph>::size_type* accum_deg,
                 typename graph_traits<Graph>::size_type num_blocks)
{
  #pragma mta noalias g
  #pragma mta noalias rank
  #pragma mta noalias acc
  #pragma mta noalias in_degrees
  #pragma mta noalias *accum_deg

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  size_type order = num_vertices(g);

  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos =
      begin_block_range(accum_deg[order], block_id, num_blocks);
    size_type end_pos =
      end_block_range(accum_deg[order], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start_pos);
    size_type end_outer =
      end_manhattan_outer_range(accum_deg, order, end_pos);

    thread_vertex_iterator tverts = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

      vertex_descriptor u = *tverts;
      thread_adjacency_iterator tadj_verts =
        thread_adjacent_vertices(u, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
      {
        vertex_descriptor v = *tadj_verts;
        mt_incr(acc[u], rank[v] / in_degrees[v]);
      }
    }
  }
}

template <typename Graph, typename DegreeMap>
void compute_in_degrees(Graph& g, DegreeMap& in_degrees,
                        typename graph_traits<Graph>::size_type* accum_deg,
                        typename graph_traits<Graph>::size_type num_blocks)
{
  #pragma mta noalias g
  #pragma mta noalias in_degrees
  #pragma mta noalias *accum_deg

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  size_type order = num_vertices(g);

#ifdef USING_QTHREADS
  detail::init_in_degrees_qt<Graph, DegreeMap> iid(g, in_degrees);
  qt_loop_balance(0, order, iid);
  detail::compute_in_degrees_qt<Graph, DegreeMap>
    cidq(g, in_degrees, accum_deg, num_blocks);
  qt_loop_balance(0, num_blocks, cidq);
#else
  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      in_degrees[*tverts] = 0;
    }
  }

  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos =
      begin_block_range(accum_deg[order], block_id, num_blocks);
    size_type end_pos =
      end_block_range(accum_deg[order], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start_pos);
    size_type end_outer =
      end_manhattan_outer_range(accum_deg, order, end_pos);

    thread_vertex_iterator tverts = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

      vertex_descriptor u = *tverts;
      thread_adjacency_iterator tadj_verts =
        thread_adjacent_vertices(u, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
      {
        mt_incr(in_degrees[*tadj_verts], 1);
      }
    }
  }
#endif
}

template <typename Graph, typename DegreeMap>
void compute_in_degrees(Graph& g, DegreeMap& in_degrees)
{
  #pragma mta noalias g
  #pragma mta noalias in_degrees

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  // Get the accumulation array of the vertex out degrees.
  size_type* accum_deg = (size_type*) malloc(sizeof(size_type) * (order + 1));
  accumulate_out_degree(accum_deg, g);

#ifdef USING_QTHREADS
  size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
  size_type num_blocks = (accum_deg[order] + BR_CHUNK - 1) / BR_CHUNK;
#else
  size_type num_blocks = 1;
#endif

  compute_in_degrees(g, in_degrees, accum_deg, num_blocks);

  free(accum_deg);
}

}

template <typename Graph, typename RankMap, typename VertexArray>
void badrank(Graph& g, RankMap& rank, VertexArray& bad_verts,
             typename graph_traits<Graph>::size_type num_bad,
             double delta = .00001, double dampen = .8)
{
  #pragma mta noalias g
  #pragma mta noalias rank
  #pragma mta noalias bad_verts

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  typedef vertex_property_map<Graph, double> InitMap;
  InitMap init(g);

  typedef vertex_property_map<Graph, size_type> DegreeMap;
  DegreeMap in_degrees(g);

  typedef vertex_property_map<Graph, double> AccMap;
  AccMap acc(g);

  // Initialize all vertices as good vertices.
#ifdef USING_QTHREADS
  detail::init_propmaps<Graph, RankMap, InitMap, AccMap> ip(g, rank, init, acc);
  qt_loop_balance(0, order, ip);
#else
  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos =
      begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;
      rank[v] = 0;
      init[v] = 0.0;
      acc[v] = 0.0;
    }
  }
#endif

  // Reinitialize the bad vertices.
  for (size_type i = 0; i < num_bad; ++i)
  {
    vertex_descriptor v = bad_verts[i];
    rank[v] = 1;
    init[v] = 1.0;
  }

  // Get the accumulation array of the vertex out degrees.
  size_type* accum_deg = (size_type*) malloc(sizeof(size_type) * (order + 1));
  accumulate_out_degree(accum_deg, g);

#ifdef USING_QTHREADS
    size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
    size_type num_blocks = (accum_deg[order] + BR_CHUNK - 1) / BR_CHUNK;
#else
    size_type num_blocks = 1;
#endif

  detail::compute_in_degrees(g, in_degrees, accum_deg, num_blocks);

  double maxdiff = 0;

#ifdef USING_QTHREADS
  // Used to store the indidual values for the maxdiff reduction.
  double* maxdiff_vec = (double*) malloc(sizeof(double) * (order));
#endif

  do
  {
#ifdef USING_QTHREADS
    detail::compute_acc_qt<Graph, RankMap, AccMap, DegreeMap>
      caq(g, rank, acc, in_degrees, accum_deg, num_blocks);
    qt_loop_balance(0, num_blocks, caq);
#else
    detail::compute_acc(g, rank, acc, in_degrees, accum_deg, num_blocks);
#endif

#ifdef USING_QTHREADS
    detail::badrank_body_qt<Graph, RankMap, InitMap, AccMap>
      bbq(g, rank, init, acc, dampen, maxdiff_vec);
    qt_loop_balance(0, order, bbq);
    maxdiff = qt_double_max(maxdiff_vec, order, 0);
#else
    maxdiff = 0;

    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;

        double oldval = rank[v];
        double newval;

#ifdef PUB_ALG
        rank[v] = init[v] + dampen * acc[v];
        newval = rank[v];
#else
        if (out_degree(v, g) == 0)
        {
          newval = init[v];
        }
        else
        {
          rank[v] = init[v] + dampen * acc[v] / out_degree(v, g);
          newval = rank[v];
        }
#endif

        double absdiff = (oldval > newval) ? oldval - newval : newval - oldval;

        if (absdiff > maxdiff) maxdiff = absdiff;

        acc[v] = 0.0;
      }
    }
#endif
  } while (maxdiff > delta);

#ifdef USING_QTHREADS
  free(maxdiff_vec);
#endif
  free(accum_deg);
}

}

#undef BR_CHUNK

#endif
