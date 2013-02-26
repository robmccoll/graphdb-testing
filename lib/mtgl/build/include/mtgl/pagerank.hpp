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
/*! \file pagerank.hpp

    \author Jon Berry (jberry@sandia.gov)
    \author Karen Devine (kddevin@sandia.gov)

    \date 8/8/2008
*/
/****************************************************************************/

#ifndef MTGL_PAGERANK_HPP
#define MTGL_PAGERANK_HPP

#include <iostream>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/partitioning.hpp>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

#define PR_CHUNK 256

namespace mtgl {

namespace detail {

#ifdef USING_QTHREADS
template <typename Graph, typename RankMap>
class initialize_rank {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  initialize_rank(Graph& gg, RankMap& rm) : g(gg), rank(rm) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type order = num_vertices(g);

    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      rank[*tverts] = 1.0 / order;
    }
  }

private:
  Graph& g;
  RankMap& rank;
};

template <typename Graph, typename RankMap, typename AccMap>
class prepare_compute_acc {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  prepare_compute_acc(Graph& gg, RankMap& r, AccMap& a, double& s) :
    g(gg), rank(r), acc(a), sum(s) {}

  void operator()(const size_t start, const size_t stop)
  {
    double my_sum = 0.0;

    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;
      acc[v] = 0.0;
      my_sum += rank[v];
    }

    mt_incr(sum, my_sum);
  }

private:
  Graph& g;
  RankMap& rank;
  AccMap& acc;
  double& sum;
};

template <typename Graph, typename RankMap, typename AccMap,
          typename Direction = typename graph_traits<Graph>::directed_category>
class compute_acc_outer {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  compute_acc_outer(Graph& gg, RankMap& r, AccMap& a,
                    size_type* ad, size_type nb) :
    g(gg), rank(r), acc(a), accum_deg(ad), num_blocks(nb) {}

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
        size_type deg = out_degree(u, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor v = *tadj_verts;
          mt_incr(acc[v], rank[u] / deg);
        }
      }
    }
  }

private:
  Graph& g;
  RankMap& rank;
  AccMap& acc;
  size_type* accum_deg;
  size_type num_blocks;
};

template <typename Graph, typename RankMap, typename AccMap>
class compute_acc_outer<Graph, RankMap, AccMap, bidirectionalS> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_in_adjacency_iterator
          thread_in_adjacency_iterator;

  compute_acc_outer(Graph& gg, RankMap& r, AccMap& a,
                    size_type* ad, size_type nb) :
    g(gg), rank(r), acc(a), accum_deg(ad), num_blocks(nb) {}

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

        vertex_descriptor v = *tverts;
        thread_in_adjacency_iterator tadj_verts =
          thread_in_adjacent_vertices(v, start_inner, g);

        double total = 0.0;

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor u = *tadj_verts;

          total += rank[u] / out_degree(u, g);
        }

        mt_incr(acc[v], total);
      }
    }
  }

private:
  Graph& g;
  RankMap& rank;
  AccMap& acc;
  size_type* accum_deg;
  size_type num_blocks;
};

template <typename Graph, typename RankMap>
class adjust_zero_outdeg {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  adjust_zero_outdeg(Graph& gg, RankMap& r, double& s) :
    g(gg), rank(r), sum(s) {}

  void operator()(const size_t start, const size_t stop)
  {
    double my_sum = 0.0;

    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;
      my_sum += (out_degree(v, g) == 0) * rank[v];
    }

    mt_incr(sum, my_sum);
  }

private:
  Graph& g;
  RankMap& rank;
  double& sum;
};

template <typename Graph, typename AccMap>
class compute_norm_vec {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  compute_norm_vec(Graph& gg, AccMap& am, double* v, double a, double d) :
    g(gg), acc(am), vec(v), adjustment(a), dampen(d) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;
      acc[v] = adjustment + (dampen * acc[v]);
      vec[start_pos] = acc[v] >= 0 ? acc[v] : -1 * acc[v];
    }
  }

private:
  Graph& g;
  AccMap& acc;
  double* vec;
  double adjustment;
  double dampen;
};

template <typename Graph, typename RankMap, typename AccMap>
class compute_diff_vec {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  compute_diff_vec(Graph& gg, RankMap& r, AccMap& a, double* vec, double n) :
    g(gg), rank(r), acc(a), vec(vec), norm(n) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;

      double oldval = rank[v];
      double newval = acc[v] / norm;

      rank[v] = newval;
      vec[start_pos] = oldval > newval ?  oldval - newval : newval - oldval;
    }
  }

private:
  Graph& g;
  RankMap& rank;
  AccMap& acc;
  double* vec;
  double norm;
};
#endif

template <typename Graph,
          typename Direction = typename graph_traits<Graph>::directed_category>
class accumulate_dir_degree {
public:
  template <typename T>
  void
  operator()(T* accum_deg, Graph& g) { accumulate_out_degree(accum_deg, g); }
};

template <typename Graph>
class accumulate_dir_degree<Graph, bidirectionalS> {
public:
  template <typename T>
  void
  operator()(T* accum_deg, Graph& g) { accumulate_in_degree(accum_deg, g); }
};

template <typename Graph, typename RankMap, typename AccMap,
          typename Direction = typename graph_traits<Graph>::directed_category>
class compute_acc {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  void operator()(Graph& g, RankMap& rank, AccMap& acc,
                  size_type* accum_deg, size_type num_blocks)
  {
    #pragma mta noalias g
    #pragma mta noalias rank
    #pragma mta noalias acc
    #pragma mta noalias *accum_deg

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
        size_type deg = out_degree(u, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor v = *tadj_verts;
          mt_incr(acc[v], rank[u] / deg);
        }
      }
    }
  }
};

// There is a hotspot on accumulating the acc of vertices with high indegree.
// This version uses the in edges present in a bidirectional graph to avoid
// the hotspot.
template <typename Graph, typename RankMap, typename AccMap>
class compute_acc<Graph, RankMap, AccMap, bidirectionalS> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_in_adjacency_iterator
          thread_in_adjacency_iterator;

  void operator()(Graph& g, RankMap& rank, AccMap& acc,
                  size_type* accum_deg, size_type num_blocks)
  {
    #pragma mta noalias g
    #pragma mta noalias rank
    #pragma mta noalias acc
    #pragma mta noalias *accum_deg

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

        vertex_descriptor v = *tverts;
        thread_in_adjacency_iterator tadj_verts =
          thread_in_adjacent_vertices(v, start_inner, g);

        double total = 0.0;

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor u = *tadj_verts;
          total += rank[u] / out_degree(u, g);
        }

        mt_incr(acc[v], total);
      }
    }
  }
};

}

template <typename Graph, typename RankMap>
void pagerank(Graph& g, RankMap& rank,
              double delta = .00001, double dampen = .8)
{
  #pragma mta noalias g
  #pragma mta noalias rank

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  typedef vertex_property_map<Graph, double> AccMap;
  AccMap acc(g);

#ifdef USING_QTHREADS
  // Used to store the indidual values for the norm and maxdiff reductions.
  double* rank_update = (double*) malloc(sizeof(double) * (order));
#endif

#ifdef DEBUG
  double time_prep = 0.0;
  double time_acc = 0.0;
  double time_zdc = 0.0;
  double time_norm = 0.0;
  double time_maxdiff = 0.0;
  mt_timer timer;
#endif

#ifdef USING_QTHREADS
  detail::initialize_rank<Graph, RankMap> cz(g, rank);
  qt_loop_balance(0, order, cz);
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
      rank[*tverts] = 1.0 / order;
    }
  }
#endif

  // Get the accumulation array of the vertex degrees.
  //   undirectedS, directedS - out_degree()
  //   bidirectionalS - in_degree()
  size_type* accum_deg = (size_type*) malloc(sizeof(size_type) * (order + 1));
  detail::accumulate_dir_degree<Graph> add;
  add(accum_deg, g);

#ifdef USING_QTHREADS
    size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
    size_type num_blocks = (accum_deg[order] + PR_CHUNK - 1) / PR_CHUNK;
#else
    size_type num_blocks = 1;
#endif

  int iter_cnt = 0;
  double maxdiff = 0.0;

  do
  {
    ++iter_cnt;

    double sum = 0.0;

#ifdef DEBUG
    timer.start();
#endif

#ifdef USING_QTHREADS
    detail::prepare_compute_acc<Graph, RankMap, AccMap> pca(g, rank, acc, sum);
    qt_loop_balance(0, order, pca);
#else
    stream_id = 0;
    num_streams = 1;

    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;

        acc[v] = 0.0;
        sum += rank[v];
      }
    }
#endif

#ifdef DEBUG
    timer.stop();
    time_prep += timer.getElapsedSeconds();
    timer.start();
#endif

#ifdef USING_QTHREADS
    detail::compute_acc_outer<Graph, RankMap, AccMap>
      cacc(g, rank, acc, accum_deg, num_blocks);
    qt_loop_balance(0, num_blocks, cacc);
#else
    detail::compute_acc<Graph, RankMap, AccMap> cacc;
    cacc(g, rank, acc, accum_deg, num_blocks);
#endif

    double adjustment = (1 - dampen) / order * sum;

#ifdef DEBUG
    timer.stop();
    time_acc += timer.getElapsedSeconds();
    timer.start();
#endif

    // Adjustment for zero-outdegree vertices.
    sum = 0;

#ifdef USING_QTHREADS
    detail::adjust_zero_outdeg<Graph, RankMap> azo(g, rank, sum);
    qt_loop_balance(0, order, azo);
#else
    stream_id = 0;
    num_streams = 1;

    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;
        sum += (out_degree(v, g) == 0) * rank[v];
      }
    }
#endif

    adjustment += dampen * sum / order;

#ifdef DEBUG
    timer.stop();
    time_zdc += timer.getElapsedSeconds();
    timer.start();
#endif

    // Compute new solution vector and scaling factor.
    double norm = 0.0;

#ifdef USING_QTHREADS
    detail::compute_norm_vec<Graph, AccMap> cnv(g, acc, rank_update,
                                                adjustment, dampen);
    qt_loop_balance(0, order, cnv);
    norm = qt_double_max(rank_update, order, 0);
#else
    stream_id = 0;
    num_streams = 1;

    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;

        acc[v] = adjustment + (dampen * acc[v]);

        double tmp = acc[v] >= 0 ? acc[v] : -1 * acc[v];
        if (tmp > norm) norm = tmp;
      }
    }
#endif

#ifdef DEBUG
    timer.stop();
    time_norm += timer.getElapsedSeconds();
    timer.start();
#endif

    maxdiff = 0;

#ifdef USING_QTHREADS
    detail::compute_diff_vec<Graph, RankMap, AccMap> cdv(g, rank, acc,
                                                         rank_update, norm);
    qt_loop_balance(0, order, cdv);
    maxdiff = qt_double_max(rank_update, order, 0);
#else
    stream_id = 0;
    num_streams = 1;

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
        double newval = acc[v] / norm;
        rank[v] = newval;

        double absdiff = oldval > newval ? oldval - newval : newval - oldval;
        if (absdiff > maxdiff) maxdiff = absdiff;
      }
    }
#endif

#ifdef __MTA__
    // The next line is only here to get this to compile on the XMT.
    #pragma mta fence
#endif

#ifdef DEBUG
    timer.stop();
    time_maxdiff += timer.getElapsedSeconds();
#endif
  } while (maxdiff > delta);

#ifdef DEBUG
  std::cout << "        Iterations: " << iter_cnt << std::endl
            << "         Prep time: " << time_prep << std::endl
            << " Accumulation time: " << time_acc << std::endl
            << "    Zero comp time: " << time_zdc << std::endl
            << "Normalization time: " << time_norm << std::endl
            << "      Maxdiff time: " << time_maxdiff << std::endl;
#endif

#ifdef USING_QTHREADS
  free(rank_update);
#endif
  free(accum_deg);
}

}

//#undef DEBUG
#undef PR_CHUNK

#endif
