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
/*! \file test_graph_parallel.cpp

    \brief Tests how the graph iterators parallelize.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 10/3/2011
*/
/****************************************************************************/

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/partitioning.hpp>
#include <mtgl/mtgl_test.hpp>

#include <iostream>
#include <iomanip>

#define CHUNK_SIZE 128

using namespace mtgl;

template <typename Graph>
void test_vector_iterators_func(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << std::endl
            << "Testing vertex_iterator:" << std::endl;
  vertex_iterator vIter = vertices(g);
  size_type sum = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = vIter[i];
    size_type uid = get(vipm, u);
    mt_incr(sum, uid);
  }
  std::cout << "Vertex id sum: " << sum << std::endl;

  std::cout << std::endl << "Testing edge_iterator:" << std::endl;
  edge_iterator eIter = edges(g);
  sum = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = eIter[i];
    size_type uid = get(vipm, source(e, g));
    size_type vid = get(vipm, target(e, g));
    mt_incr(sum, uid + vid);
  }
  std::cout << "Endpoint id sum: " << sum << std::endl;

  std::cout << std::endl << "Testing adjacency_iterator:" << std::endl;
  sum = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = vIter[i];

    size_type out_deg = out_degree(u, g);
    adjacency_iterator adjIter = adjacent_vertices(u, g);
    #pragma mta assert nodep
    for (size_type j = 0; j < out_deg; ++j)
    {
      vertex_descriptor v = adjIter[j];
      size_type vid = get(vipm, v);
      mt_incr(sum, vid);
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;

  std::cout << std::endl << "Testing out_edge_iterator:" << std::endl;
  sum = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = vIter[i];

    size_type out_deg = out_degree(u, g);
    out_edge_iterator oeIter = out_edges(u, g);
    #pragma mta assert nodep
    for (size_type j = 0; j < out_deg; ++j)
    {
      edge_descriptor e = oeIter[j];
      size_type vid = get(vipm, target(e, g));
      mt_incr(sum, vid);
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;
}

template <typename Graph>
void test_bid_vector_iterators_func(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::in_adjacency_iterator
          in_adjacency_iterator;
  typedef typename graph_traits<Graph>::in_edge_iterator in_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << std::endl << "Testing in_adjacency_iterator:" << std::endl;
  vertex_iterator vIter = vertices(g);
  size_type sum = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = vIter[i];

    size_type in_deg = in_degree(v, g);
    in_adjacency_iterator adjIter = in_adjacent_vertices(v, g);
    #pragma mta assert nodep
    for (size_type j = 0; j < in_deg; ++j)
    {
      vertex_descriptor u = adjIter[j];
      size_type uid = get(vipm, u);
      mt_incr(sum, uid);
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;

  std::cout << std::endl << "Testing in_edge_iterator:" << std::endl;
  sum = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = vIter[i];

    size_type in_deg = in_degree(v, g);
    in_edge_iterator ieIter = in_edges(v, g);
    #pragma mta assert nodep
    for (size_type j = 0; j < in_deg; ++j)
    {
      edge_descriptor e = ieIter[j];
      size_type uid = get(vipm, target(e, g));
      mt_incr(sum, uid);
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;
}

template <typename Graph>
void test_thread_iterators_func(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;
  typedef typename graph_traits<Graph>::thread_edge_iterator
          thread_edge_iterator;
  typedef typename graph_traits<Graph>::thread_out_edge_iterator
          thread_out_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  std::cout << std::endl << "Testing thread_vertex_iterator:" << std::endl;
  size_type sum = 0;
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start = begin_block_range(order, stream_id, num_streams);
    size_type end = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tvIter = thread_vertices(start, g);

    for (size_type i = start; i != end; ++i, ++tvIter)
    {
      size_type uid = get(vipm, *tvIter);
      mt_incr(sum, uid);
    }
  }
  std::cout << "Vertex id sum: " << sum << std::endl;

  std::cout << std::endl << "Testing thread_edge_iterator:" << std::endl;
  sum = 0;
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start = begin_block_range(size, stream_id, num_streams);
    size_type end = end_block_range(size, stream_id, num_streams);

    thread_edge_iterator teIter = thread_edges(start, g);

    for (size_type i = start; i != end; ++i, ++teIter)
    {
      edge_descriptor e = *teIter;
      size_type uid = get(vipm, source(e, g));
      size_type vid = get(vipm, target(e, g));
      mt_incr(sum, uid + vid);
    }
  }
  std::cout << "Endpoint id sum: " << sum << std::endl;

  // Get a vector of the degrees of all the vertices.  This is required to
  // make the compiler perform the linear reduction for the accumulation.
  size_type* accum_deg = (size_type*) malloc(sizeof(size_type) * (order + 1));
  accum_deg[0] = 0;
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start = begin_block_range(order, stream_id, num_streams);
    size_type end = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tvIter = thread_vertices(start, g);

    for (size_type i = start; i != end; ++i, ++tvIter)
    {
      accum_deg[i + 1] = out_degree(*tvIter, g);
    }
  }

  // Calculate the accumulated degree which is the number of adjacencies that
  // occur before a given vertex.
  for (size_type i = 1; i <= order; ++i) accum_deg[i] += accum_deg[i - 1];

#ifdef __MTA__
    size_type num_blocks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
#else
    size_type num_blocks = 1;
#endif

  std::cout << std::endl << "Testing thread_adjacency_iterator:" << std::endl;
  sum = 0;
  #pragma mta assert parallel
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start = begin_block_range(size, block_id, num_blocks);
    size_type end = end_block_range(size, block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start);
    size_type end_outer = end_manhattan_outer_range(accum_deg, order, end);

    thread_vertex_iterator tvIter = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tvIter)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end, end_outer, i);

      thread_adjacency_iterator tadjIter =
        thread_adjacent_vertices(*tvIter, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tadjIter)
      {
        size_type vid = get(vipm, *tadjIter);
        mt_incr(sum, vid);
      }
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;

  std::cout << std::endl << "Testing thread_out_edge_iterator:" << std::endl;
  sum = 0;
  #pragma mta assert parallel
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start = begin_block_range(size, block_id, num_blocks);
    size_type end = end_block_range(size, block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start);
    size_type end_outer = end_manhattan_outer_range(accum_deg, order, end);

    thread_vertex_iterator tvIter = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tvIter)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end, end_outer, i);

      thread_out_edge_iterator toeIter =
        thread_out_edges(*tvIter, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++toeIter)
      {
        edge_descriptor e = *toeIter;
        size_type vid = get(vipm, target(e, g));
        mt_incr(sum, vid);
      }
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;

  free(accum_deg);
}

template <typename Graph>
void test_bid_thread_iterators_func(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_in_adjacency_iterator
          thread_in_adjacency_iterator;
  typedef typename graph_traits<Graph>::thread_in_edge_iterator
          thread_in_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Get a vector of the degrees of all the vertices.  This is required to
  // make the compiler perform the linear reduction for the accumulation.
  size_type* accum_deg = (size_type*) malloc(sizeof(size_type) * (order + 1));
  accum_deg[0] = 0;
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start = begin_block_range(order, stream_id, num_streams);
    size_type end = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tvIter = thread_vertices(start, g);

    for (size_type i = start; i != end; ++i, ++tvIter)
    {
      accum_deg[i + 1] = in_degree(*tvIter, g);
    }
  }

  // Calculate the accumulated degree which is the number of adjacencies that
  // occur before a given vertex.
  for (size_type i = 1; i <= order; ++i) accum_deg[i] += accum_deg[i - 1];

#ifdef __MTA__
    size_type num_blocks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
#else
    size_type num_blocks = 1;
#endif

  std::cout << std::endl
            << "Testing thread_in_adjacency_iterator:" << std::endl;
  size_type sum = 0;
  #pragma mta assert parallel
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start = begin_block_range(size, block_id, num_blocks);
    size_type end = end_block_range(size, block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start);
    size_type end_outer = end_manhattan_outer_range(accum_deg, order, end);

    thread_vertex_iterator tvIter = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tvIter)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end, end_outer, i);

      thread_in_adjacency_iterator tadjIter =
        thread_in_adjacent_vertices(*tvIter, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tadjIter)
      {
        size_type uid = get(vipm, *tadjIter);
        mt_incr(sum, uid);
      }
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;

  std::cout << std::endl << "Testing thread_in_edge_iterator:" << std::endl;
  sum = 0;
  #pragma mta assert parallel
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start = begin_block_range(size, block_id, num_blocks);
    size_type end = end_block_range(size, block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start);
    size_type end_outer = end_manhattan_outer_range(accum_deg, order, end);

    thread_vertex_iterator tvIter = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tvIter)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end, end_outer, i);

      thread_in_edge_iterator tieIter =
        thread_in_edges(*tvIter, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tieIter)
      {
        edge_descriptor e = *tieIter;
        size_type uid = get(vipm, target(e, g));
        mt_incr(sum, uid);
      }
    }
  }
  std::cout << "Adjacency sum: " << sum << std::endl;

  free(accum_deg);
}

template <typename Graph,
          typename Direction = typename Graph::directed_category>
struct test_vector_iterators {
  void operator()(Graph& g)
  {
    test_vector_iterators_func(g);
  }
};

template <typename Graph>
struct test_vector_iterators<Graph, bidirectionalS> {
  void operator()(Graph& g)
  {
    test_vector_iterators_func(g);
    test_bid_vector_iterators_func(g);
  }
};

template <typename Graph,
          typename Direction = typename Graph::directed_category>
struct test_thread_iterators {
  void operator()(Graph& g)
  {
    test_thread_iterators_func(g);
  }
};

template <typename Graph>
struct test_thread_iterators<Graph, bidirectionalS> {
  void operator()(Graph& g)
  {
    test_thread_iterators_func(g);
    test_bid_thread_iterators_func(g);
  }
};

template <typename Graph, typename Iterator = typename Graph::iterator_category>
struct test_iterators {
  void operator()(Graph& g)
  {
    test_vector_iterators<Graph> tvi;
    tvi(g);
    test_thread_iterators<Graph> tti;
    tti(g);
  }
};

template <typename Graph>
struct test_iterators<Graph, vector_iterators> {
  void operator()(Graph& g)
  {
    test_vector_iterators<Graph> tvi;
    tvi(g);
  }
};

template <typename Graph>
struct test_iterators<Graph, thread_iterators> {
  void operator()(Graph& g)
  {
    test_thread_iterators<Graph> tti;
    tti(g);
  }
};

//#define TEST_STINGER

int main(int argc, char *argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef directedS DIRECTION;
//  typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;
#endif

  typedef graph_traits<Graph>::size_type size_type;

  init_test(argc, argv);

#ifdef TEST_STINGER
  SGraph sg(1);
  Graph g(sg);
#else
  Graph g;
#endif

  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "order: " << std::setw(2) << order << std::endl
            << " size: " << std::setw(2) << size << std::endl;

  if (order < 10)
  {
    std::cout << "The graph:" << std::endl;
    print(g);
  }

  test_iterators<Graph> ti;
  ti(g);

  return 0;
}
