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
/*! \file test_graph.cpp

    \brief Tests various graph functions and iterators.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 10/3/2011
*/
/****************************************************************************/

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/partitioning.hpp>

#include <iostream>
#include <iomanip>

using namespace mtgl;

template <typename Graph,
          typename Direction = typename Graph::directed_category>
struct print_degree {
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  void operator()(Graph& g, vertex_descriptor& u)
  {
    size_type deg = degree(u, g);
    std::cout << std::setw(2) << deg << ", ";

    size_type out_deg = out_degree(u, g);
    std::cout << std::setw(2) << out_deg << ")" << std::endl;
  }
};

template <typename Graph>
struct print_degree<Graph, bidirectionalS> {
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  void operator()(Graph& g, vertex_descriptor& u)
  {
    size_type deg = degree(u, g);
    std::cout << std::setw(2) << deg << ", ";

    size_type in_deg = in_degree(u, g);
    std::cout << std::setw(2) << in_deg << ", ";

    size_type out_deg = out_degree(u, g);
    std::cout << std::setw(2) << out_deg << ")" << std::endl;
  }
};

template <typename Graph>
void print_vector_iterators(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);
  edge_id_map<Graph> eipm = get(_edge_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << std::endl
            << "Testing vertex_iterator and degree functions:" << std::endl;
  vertex_iterator vIter = vertices(g);
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = vIter[i];
    size_type uid = get(vipm, u);
    std::cout << std::setw(2) << uid << ": (";

    print_degree<Graph> pd;
    pd(g, u);
  }

  std::cout << std::endl << "Testing edge_iterator:" << std::endl;
  edge_iterator eIter = edges(g);
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = eIter[i];
    size_type eid = get(eipm, e);

    vertex_descriptor src = source(e, g);
    vertex_descriptor trg = target(e, g);

    size_type sid = get(vipm, src);
    size_type tid = get(vipm, trg);

    std::cout << std::setw(2) << eid << ": (" << std::setw(2) << sid << " -> "
              << std::setw(2) << tid << ")" << std::endl;
  }

  std::cout << std::endl << "Testing adjacency_iterator:" << std::endl;
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = vIter[i];

    size_type uid = get(vipm, u);
    std::cout << std::setw(2) << uid << ":";

    size_type out_deg = out_degree(u, g);
    adjacency_iterator adjIter = adjacent_vertices(u, g);
    for (size_type j = 0; j < out_deg; ++j)
    {
      vertex_descriptor v = adjIter[j];
      size_type vid = get(vipm, v);

      std::cout << "  " << std::setw(2) << vid;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << "Testing out_edge_iterator:" << std::endl;
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = vIter[i];

    size_type uid = get(vipm, u);
    std::cout << std::setw(2) << uid << ":";

    size_type out_deg = out_degree(u, g);
    out_edge_iterator oeIter = out_edges(u, g);
    for (size_type j = 0; j < out_deg; ++j)
    {
      edge_descriptor e = oeIter[j];
      size_type eid = get(eipm, e);

      vertex_descriptor src = source(e, g);
      vertex_descriptor trg = target(e, g);

      size_type sid = get(vipm, src);
      size_type tid = get(vipm, trg);

      std::cout << "  {" << std::setw(2) << eid << "}: (" << std::setw(2)
                << sid << " -> " << std::setw(2) << tid << ")";
    }
    std::cout << std::endl;
  }
}

template <typename Graph>
void print_bid_vector_iterators(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::in_adjacency_iterator
          in_adjacency_iterator;
  typedef typename graph_traits<Graph>::in_edge_iterator in_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);
  edge_id_map<Graph> eipm = get(_edge_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << std::endl << "Testing in_adjacency_iterator:" << std::endl;
  vertex_iterator vIter = vertices(g);
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = vIter[i];

    size_type vid = get(vipm, v);
    std::cout << std::setw(2) << vid << ":";

    size_type in_deg = in_degree(v, g);
    in_adjacency_iterator adjIter = in_adjacent_vertices(v, g);
    for (size_type j = 0; j < in_deg; ++j)
    {
      vertex_descriptor u = adjIter[j];
      size_type uid = get(vipm, u);

      std::cout << "  " << std::setw(2) << uid;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << "Testing in_edge_iterator:" << std::endl;
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = vIter[i];

    size_type vid = get(vipm, v);
    std::cout << std::setw(2) << vid << ":";

    size_type in_deg = in_degree(v, g);
    in_edge_iterator ieIter = in_edges(v, g);
    for (size_type j = 0; j < in_deg; ++j)
    {
      edge_descriptor e = ieIter[j];
      size_type eid = get(eipm, e);

      vertex_descriptor src = source(e, g);
      vertex_descriptor trg = target(e, g);

      size_type sid = get(vipm, src);
      size_type tid = get(vipm, trg);

      std::cout << "  {" << std::setw(2) << eid << "}: (" << std::setw(2)
                << sid << " -> " << std::setw(2) << tid << ")";
    }
    std::cout << std::endl;
  }
}

template <typename Graph>
void print_thread_iterators(Graph& g)
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
  edge_id_map<Graph> eipm = get(_edge_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  size_type num_threads = 2;

  std::cout << std::endl
            << "Testing thread_vertex_iterator and degree functions:"
            << std::endl;
  for (size_type i = 0; i < num_threads; ++i)
  {
    size_type start = begin_block_range(order, i, num_threads);
    size_type end = end_block_range(order, i, num_threads);

    thread_vertex_iterator tvIter = thread_vertices(start, g);

    for (size_type j = start; j != end; ++j, ++tvIter)
    {
      vertex_descriptor u = *tvIter;
      size_type uid = get(vipm, u);
      std::cout << std::setw(2) << uid << ": (";

      print_degree<Graph> pd;
      pd(g, u);
    }
  }

  std::cout << std::endl << "Testing thread_edge_iterator:" << std::endl;
  for (size_type i = 0; i < num_threads; ++i)
  {
    size_type start = begin_block_range(size, i, num_threads);
    size_type end = end_block_range(size, i, num_threads);

    thread_edge_iterator teIter = thread_edges(start, g);

    for (size_type j = start; j != end; ++j, ++teIter)
    {
      edge_descriptor e = *teIter;
      size_type eid = get(eipm, e);

      vertex_descriptor src = source(e, g);
      vertex_descriptor trg = target(e, g);

      size_type sid = get(vipm, src);
      size_type tid = get(vipm, trg);

      std::cout << std::setw(2) << eid << ": (" << std::setw(2) << sid << " -> "
                << std::setw(2) << tid << ")" << std::endl;
    }
  }

  std::cout << std::endl << "Testing thread_adjacency_iterator:" << std::endl;
  thread_vertex_iterator tvIter = thread_vertices(0, g);
  for (size_type i = 0; i < order; ++i, ++tvIter)
  {
    vertex_descriptor u = *tvIter;

    size_type uid = get(vipm, u);
    std::cout << std::setw(2) << uid << ":";

    size_type out_deg = out_degree(u, g);
    for (size_type j = 0; j < num_threads; ++j)
    {
      size_type start = begin_block_range(out_deg, j, num_threads);
      size_type end = end_block_range(out_deg, j, num_threads);

      thread_adjacency_iterator tadjIter =
        thread_adjacent_vertices(u, start, g);

      for (size_type k = start; k != end; ++k, ++tadjIter)
      {
        vertex_descriptor v = *tadjIter;
        size_type vid = get(vipm, v);

        std::cout << "  " << std::setw(2) << vid;
      }
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << "Testing thread_out_edge_iterator:" << std::endl;
  tvIter = thread_vertices(0, g);
  for (size_type i = 0; i < order; ++i, ++tvIter)
  {
    vertex_descriptor u = *tvIter;

    size_type uid = get(vipm, u);
    std::cout << std::setw(2) << uid << ":";

    size_type out_deg = out_degree(u, g);
    for (size_type j = 0; j < num_threads; ++j)
    {
      size_type start = begin_block_range(out_deg, j, num_threads);
      size_type end = end_block_range(out_deg, j, num_threads);

      thread_out_edge_iterator toeIter = thread_out_edges(u, start, g);

      for (size_type k = start; k != end; ++k, ++toeIter)
      {
        edge_descriptor e = *toeIter;
        size_type eid = get(eipm, e);

        vertex_descriptor src = source(e, g);
        vertex_descriptor trg = target(e, g);

        size_type sid = get(vipm, src);
        size_type tid = get(vipm, trg);

        std::cout << "  {" << std::setw(2) << eid << "}: (" << std::setw(2)
                  << sid << " -> " << std::setw(2) << tid << ")";
      }
    }
    std::cout << std::endl;
  }
}

template <typename Graph>
void print_bid_thread_iterators(Graph& g)
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
  edge_id_map<Graph> eipm = get(_edge_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  size_type num_threads = 2;

  std::cout << std::endl
            << "Testing thread_in_adjacency_iterator:" << std::endl;
  thread_vertex_iterator tvIter = thread_vertices(0, g);
  for (size_type i = 0; i < order; ++i, ++tvIter)
  {
    vertex_descriptor v = *tvIter;

    size_type vid = get(vipm, v);
    std::cout << std::setw(2) << vid << ":";

    size_type in_deg = in_degree(v, g);
    for (size_type j = 0; j < num_threads; ++j)
    {
      size_type start = begin_block_range(in_deg, j, num_threads);
      size_type end = end_block_range(in_deg, j, num_threads);

      thread_in_adjacency_iterator tadjIter =
        thread_in_adjacent_vertices(v, start, g);

      for (size_type k = start; k != end; ++k, ++tadjIter)
      {
        vertex_descriptor u = *tadjIter;
        size_type uid = get(vipm, u);

        std::cout << "  " << std::setw(2) << uid;
      }
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << "Testing thread_in_edge_iterator:" << std::endl;
  tvIter = thread_vertices(0, g);
  for (size_type i = 0; i < order; ++i, ++tvIter)
  {
    vertex_descriptor v = *tvIter;

    size_type vid = get(vipm, v);
    std::cout << std::setw(2) << vid << ":";

    size_type in_deg = in_degree(v, g);
    for (size_type j = 0; j < num_threads; ++j)
    {
      size_type start = begin_block_range(in_deg, j, num_threads);
      size_type end = end_block_range(in_deg, j, num_threads);

      thread_in_edge_iterator tieIter = thread_in_edges(v, start, g);

      for (size_type k = start; k != end; ++k, ++tieIter)
      {
        edge_descriptor e = *tieIter;
        size_type eid = get(eipm, e);

        vertex_descriptor src = source(e, g);
        vertex_descriptor trg = target(e, g);

        size_type sid = get(vipm, src);
        size_type tid = get(vipm, trg);

        std::cout << "  {" << std::setw(2) << eid << "}: (" << std::setw(2)
                  << sid << " -> " << std::setw(2) << tid << ")";
      }
    }
    std::cout << std::endl;
  }
}

template <typename Graph,
          typename Direction = typename Graph::directed_category>
struct test_vector_iterators {
  void operator()(Graph& g)
  {
    print_vector_iterators(g);
  }
};

template <typename Graph>
struct test_vector_iterators<Graph, bidirectionalS> {
  void operator()(Graph& g)
  {
    print_vector_iterators(g);
    print_bid_vector_iterators(g);
  }
};

template <typename Graph,
          typename Direction = typename Graph::directed_category>
struct test_thread_iterators {
  void operator()(Graph& g)
  {
    print_thread_iterators(g);
  }
};

template <typename Graph>
struct test_thread_iterators<Graph, bidirectionalS> {
  void operator()(Graph& g)
  {
    print_thread_iterators(g);
    print_bid_thread_iterators(g);
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

int main(int argc, char *argv[])
{
  typedef directedS DIRECTION;
//  typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;

  typedef graph_traits<Graph>::size_type size_type;

  Graph g;

  // Initialize graph.

#if 0
  // Test case 1.
  const size_type numVerts = 8;
  const size_type numEdges = 9;

  size_type sources[numEdges] = { 0, 1, 2, 2, 3, 4, 5, 6, 6 };
  size_type targets[numEdges] = { 7, 7, 0, 2, 0, 3, 4, 0, 2 };

//  size_type sources[numEdges] = { 1, 2, 4, 6, 0, 2, 3, 5, 6 };
//  size_type targets[numEdges] = { 7, 2, 3, 0, 7, 0, 0, 4, 2 };
#endif

#if 0
  // Test case 2.
  const size_type numVerts = 12;
  const size_type numEdges = 11;

  size_type sources[numEdges] = { 0, 1, 2, 3, 5, 6, 6, 7, 7,  9, 10 };
  size_type targets[numEdges] = { 1, 2, 3, 4, 4, 3, 9, 6, 8, 10, 11 };
#endif

  // Test case 3.
  const size_type numVerts = 6;
  const size_type numEdges = 8;

  size_type sources[numEdges] = { 0, 0, 1, 1, 2, 3, 3, 4 };
  size_type targets[numEdges] = { 1, 2, 2, 3, 4, 4, 5, 5 };

  init(numVerts, numEdges, sources, targets, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "order: " << std::setw(2) << order << std::endl
            << " size: " << std::setw(2) << size << std::endl << std::endl
            << "The graph:" << std::endl;

  print(g);

  std::cout << std::endl << "Testing copy constructor:" << std::endl;

  Graph g2(g);
  print(g2);

  std::cout << std::endl << "Testing assignment:" << std::endl;

  Graph g3(g);
  g3 = g;
  print(g3);

  test_iterators<Graph> ti;
  ti(g);

  std::cout << std::endl
            << "Testing is_directed(): " << is_directed(g) << std::endl
            << "Testing is_undirected(): " << is_undirected(g) << std::endl
            << "Testing is_bidirectional(): " << is_bidirectional(g)
            << std::endl;

  return 0;
}
