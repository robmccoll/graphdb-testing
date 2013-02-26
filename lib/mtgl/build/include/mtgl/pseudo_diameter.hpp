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
/*! \file pseudo_diameter.hpp

    \brief Finds an approximate diameter of a graph.

    \author Vitus Leung (vjleung@sandia.gov)

    \date 6/9/2008
*/
/****************************************************************************/

#ifndef MTGL_PSEUDO_DIAMETER_HPP
#define MTGL_PSEUDO_DIAMETER_HPP

#include <cmath>

#include <mtgl/breadth_first_search.hpp>
#include <mtgl/random.hpp>

namespace mtgl {

namespace detail {

template <typename Graph, typename NodeLevelMap>
class pd_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  pd_visitor(NodeLevelMap& nl) : node_level(nl) {}

  void tree_edge(edge_descriptor& e, Graph& g) const
  {
    node_level[target(e, g)] = node_level[source(e, g)] + 1;
  }

private:
  NodeLevelMap& node_level;
};

}

/// \brief Finds an approximate diameter of a graph.
///
/// The algorithm picks u, the starting vertex, randomly and uses a breadth
/// first search to find v, the vertex farthest from u.  Next, it uses v as
/// the starting vertex and finds the vertex farthest from it.  The algorithm
/// keeps finding the farthest vertex from the previous farthest vertex as
/// long as the distance increases.  The pseudo-diameter is the final
/// distance.
template <typename Graph>
typename graph_traits<Graph>::size_type
pseudo_diameter(Graph& g)
{
  #pragma mta noalias g

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  typedef vertex_property_map<Graph, unsigned long> NodeLevelMap;
  NodeLevelMap node_level(g);
  vertex_property_map<Graph, unsigned long> visited(g);

  vertex_descriptor farthest_vertex = *thread_vertices(
    static_cast<size_type>(rint(mt_drand48() * (order - 1))), g);

  node_level[farthest_vertex] = 0;

  detail::pd_visitor<Graph, NodeLevelMap> pdv(node_level);
  breadth_first_search(g, farthest_vertex, pdv, visited);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Make sure all the vertices in the graph were visited.
  unsigned long total_visited = 0;
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      total_visited += visited[*verts] > 0;
    }
  }

  if (total_visited != order)
  {
    if (is_undirected(g))
    {
      printf("UNDIRECTED GRAPH NOT CONNECTED!\n");
    }
    else
    {
      printf("DIRECTED GRAPH NOT STRONGLY CONNECTED!\n");
    }
    exit(1);
  }

  farthest_vertex = *thread_vertices(0, g);

  // Find the vertex farthest from the starting vertex.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      if (node_level[v] > node_level[farthest_vertex]) farthest_vertex = v;
    }
  }

  unsigned long prev_diameter = 0;

  while (prev_diameter < node_level[farthest_vertex])
  {
    prev_diameter = node_level[farthest_vertex];
    node_level[farthest_vertex] = 0;

    breadth_first_search(g, farthest_vertex, pdv, visited);

    // Make sure all the vertices in the graph were visited.
    total_visited = 0;
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++verts)
      {
        total_visited += visited[*verts] > 0;
      }
    }

    if (total_visited != order)
    {
      if (is_undirected(g))
      {
        printf("UNDIRECTED GRAPH NOT CONNECTED!\n");
      }
      else
      {
        printf("DIRECTED GRAPH NOT STRONGLY CONNECTED!\n");
      }
      exit(1);
    }

    farthest_vertex = *thread_vertices(0, g);

    // Find the vertex farthest from the starting vertex.
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++verts)
      {
        vertex_descriptor v = *verts;
        if (node_level[v] > node_level[farthest_vertex]) farthest_vertex = v;
      }
    }

#ifdef __MTA__
    // The next line is only here to get this to compile on the XMT.
    int xmt_compiler_is_retarded = 0;
#endif
  }

  return prev_diameter;
}

}

#endif
