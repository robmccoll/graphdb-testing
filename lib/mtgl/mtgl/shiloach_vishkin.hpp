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
/*! \file shiloach_vishkin.hpp

    \author Jon Berry (jberry@sandia.gov)

    \date 4/22/2005
*/
/****************************************************************************/

#ifndef MTGL_SHILOACH_VISHKIN_HPP
#define MTGL_SHILOACH_VISHKIN_HPP

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/partitioning.hpp>

#define SV_BLOCK_SIZE 128

namespace mtgl {

namespace detail {

template <typename Graph>
class default_sv_filter {
public:
  default_sv_filter() {}

  bool operator[](typename graph_traits<Graph>::edge_descriptor e)
  { return true; }
};

}

/*! \brief A simple implementation of the Shiloach-Vishkin connected
           components algorithm.  This code scales through about 10 MTA
           processors, then scaling ends due to hot-spotting.  Before 10
           processors, it is faster then the Kahan and bully algorithms;
           beyond that, it is slower.

    The algorithm itself is described well in JaJa's "Introduction to
    Parallel Algorithms."
*/
template <typename Graph, typename ComponentMap, typename EdgeFilter>
void shiloach_vishkin_no_init(Graph& g, ComponentMap& result, EdgeFilter& epm)
{
  #pragma mta noalias g
  #pragma mta noalias result
  #pragma mta noalias epm

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_edge_iterator
          thread_edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  vertex_property_map<Graph, vertex_descriptor> D(g);

  size_type size = num_edges(g);
  size_type order = num_vertices(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      D[v] = *(thread_vertices(result[v], g));
    }
  }

  int graft = 1;

  while (graft != 0)
  {
    graft = 0;

    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(size, stream_id, num_streams);
      size_type end_pos = end_block_range(size, stream_id, num_streams);

      thread_edge_iterator edgs = thread_edges(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++edgs)
      {
        edge_descriptor e = *edgs;

        if (epm[e])
        {
          vertex_descriptor u = source(e, g);
          vertex_descriptor v = target(e, g);

          if (get(vid_map, D[u]) < get(vid_map, D[v]) && D[v] == D[D[v]])
          {
            D[D[v]] = D[u];
            graft = 1;
          }

          if (get(vid_map, D[v]) < get(vid_map, D[u]) && D[u] == D[D[u]])
          {
            D[D[u]] = D[v];
            graft = 1;
          }
        }
      }
    }

    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++verts)
      {
        vertex_descriptor v = *verts;

        while (D[v] != D[D[v]])
        {
          D[v] = D[D[v]];
        }
      }
    }
  }

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      result[v] = get(vid_map, D[v]);
    }
  }
}

template <typename Graph, typename ComponentMap>
void shiloach_vishkin_no_init(Graph& g, ComponentMap& result)
{
  #pragma mta noalias g
  #pragma mta noalias result

  detail::default_sv_filter<Graph> df;
  shiloach_vishkin_no_init(g, result, df);
}

template <typename Graph, typename ComponentMap, typename EdgeFilter>
void shiloach_vishkin(Graph& g, ComponentMap& result, EdgeFilter& epm)
{
  #pragma mta noalias g
  #pragma mta noalias result
  #pragma mta noalias epm

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Initialize each vertex to have itself as its leader.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      put(result, *verts, start_pos);
    }
  }

  shiloach_vishkin_no_init(g, result, epm);
}

template <typename Graph, typename ComponentMap>
void shiloach_vishkin(Graph& g, ComponentMap& result)
{
  #pragma mta noalias g
  #pragma mta noalias result

  detail::default_sv_filter<Graph> df;

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Initialize each vertex to have itself as its leader.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      put(result, *verts, start_pos);
    }
  }

  shiloach_vishkin_no_init(g, result, df);
}

}

#undef SV_BLOCK_SIZE

#endif
