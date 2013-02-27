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
/*! \file vertex_betweenness.hpp

    \brief This will be the class that will hold the methods and structures
           to do vertex betweenness.

    \author Vitus Leung (vjleung@sandia.gov)

    \date 7/9/2008
*/
/****************************************************************************/

#ifndef MTGL_VERTEX_BETWEENNESS_HPP
#define MTGL_VERTEX_BETWEENNESS_HPP

#include <mtgl/breadth_first_search.hpp>

#ifdef __MTA__
#include <sys/mta_task.h>
#include <machine/runtime.h>
#endif

namespace mtgl {

namespace detail {

template <typename Graph, typename DistanceMap, typename SigmaMap,
          typename SucMap, typename SucCountMap>
class vb_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  vb_visitor(DistanceMap& d, SigmaMap& s, SucMap& su, SucCountMap& sc) :
    distance(d), sigma(s), successors(su), successor_count(sc) {}

  void tree_edge(edge_descriptor& e, Graph& g)
  {
    vertex_descriptor u = source(e, g);
    vertex_descriptor v = target(e, g);

    mt_write(distance[v], distance[u] + 1);

    // Add w to v's successor list.
    mt_incr(sigma[v], sigma[u]);
    successors[u][mt_incr(successor_count[u], 1)] = v;
  }

  void non_tree_edge(edge_descriptor& e, Graph& g)
  {
    vertex_descriptor u = source(e, g);
    vertex_descriptor v = target(e, g);

    size_type dist = mt_readff(distance[v]);

    // Add w to v's successor list.
    if (dist == distance[u] + 1)
    {
      mt_incr(sigma[v], sigma[u]);
      successors[u][mt_incr(successor_count[u], 1)] = v;
    }
  }

  DistanceMap& distance;
  SigmaMap& sigma;
  SucMap& successors;
  SucCountMap& successor_count;
};

}

/// \brief An implementation of vertex betweenness centrality for unweighted
///        graphs.
///
/// \param g A graph.
/// \param vbc A vertex property map of doubles that will hold the vertex
///            betweenness centrality values on exit.
template <typename Graph, typename BetweennessMap>
void vertex_betweenness(Graph& g, BetweennessMap& vbc)
{
  #pragma mta noalias g
  #pragma mta noalias vbc

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  // Memory allocation and initialization.

  size_type order = num_vertices(g);

  vertex_descriptor* Q =
    (vertex_descriptor*) malloc(order * sizeof(vertex_descriptor));

  // Distance from source vertex.
  typedef vertex_property_map<Graph, unsigned long> DistanceMap;
  DistanceMap d(g);

  // Number of shortest paths from source to a given vertex.
  typedef vertex_property_map<Graph, size_type> SigmaMap;
  SigmaMap sigma(g);

  // List of successors for a vertex. The successor list size is bounded
  // by the out-degree (degree in the case of undirected graphs).
  typedef vertex_property_map<Graph, vertex_descriptor*> SucMap;
  SucMap successors(g);

  // Needed to index the successor lists.
  typedef vertex_property_map<Graph, size_type> SucCountMap;
  SucCountMap successor_count(g);

  // Partial dependencies.
  vertex_property_map<Graph, double> delta(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      successors[v] = (vertex_descriptor*) malloc(out_degree(v, g) *
                                                  sizeof(vertex_descriptor));
    }
  }

  stream_id = 0;
  num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts) vbc[*verts] = 0;
  }

  detail::vb_visitor<Graph, DistanceMap, SigmaMap, SucMap, SucCountMap>
    vbv(d, sigma, successors, successor_count);
  vertex_property_map<Graph, size_type> visited(g);
  vertex_descriptor* buffer = 0;
  size_type buf_size = 0;
  size_type* accum_deg = 0;
  size_type accum_deg_size = 0;

  // The betweenness centrality algorithm begins here.

  // Loop: order iterations.
  thread_vertex_iterator verts = thread_vertices(0, g);
  for (size_type i = 0; i < order; ++i, ++verts)
  {
    // s is the source vertex.
    vertex_descriptor s = *verts;

    // BC of a degree-0 and degree-1 vertex is zero.
    // Don't run BFS from degree-0 and 1 vertices.
//        if (out_degree(s, g) < 2) continue; // This perhaps for XMT, big data.
    if (out_degree(s, g) < 1) continue;

    // Initialize all internal variables for iteration.
    stream_id = 0;
    num_streams = 1;

    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator my_verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++my_verts)
      {
        vertex_descriptor v = *my_verts;
        mt_purge(d[v]);
        sigma[v] = 0;
        successor_count[v] = 0;
        delta[v] = 0.0;
      }
    }

    // Run BFS from source vertex s.
    mt_write(d[s], 0);
    sigma[s] = 1;

    size_type count =
      breadth_first_search(g, s, vbv, visited, Q,
                           buffer, buf_size, accum_deg, accum_deg_size);

    if (count <= 1) continue;

    // Backward pass: accumulate dependencies.
    for (unsigned long j = d[Q[count - 1]] - 1; j > 0; --j)
    {
      stream_id = 0;
      num_streams = 1;

      #pragma mta assert nodep
      #pragma mta for all streams stream_id of num_streams
      {
        size_type start_pos = begin_block_range(order, stream_id, num_streams);
        size_type end_pos = end_block_range(order, stream_id, num_streams);

        thread_vertex_iterator my_verts = thread_vertices(start_pos, g);
        for ( ; start_pos != end_pos; ++start_pos, ++my_verts)
        {
          vertex_descriptor v = *my_verts;

          if (d[v] == j)
          {
            double sigma_v = (double) sigma[v];
            for (size_type k = 0; k < successor_count[v]; ++k)
            {
              vertex_descriptor w = successors[v][k];
              delta[v] += sigma_v / sigma[w] * (1 + delta[w]);
            }

            // s not in Q.
            vbc[v] += delta[v];
          }
        }
      }
    }
  }

  free(Q);
  free(buffer);
  free(accum_deg);

  stream_id = 0;
  num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      free(successors[*verts]);
    }
  }
}

}

#endif
