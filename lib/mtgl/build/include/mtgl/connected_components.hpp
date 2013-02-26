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
/*! \file connected_components.hpp

    \brief This file contains the main connected components algorithm.  The
           main algorithm uses a BFS on the vertex with the highest degree
           to find the largest component and then uses the Shiloach Vishkin
           algorithm to assign the remainder of the components.

    \author Jon Berry (jberry@sandia.gov)

    \date 9/2005
*/
/****************************************************************************/

#ifndef MTGL_CONNECTED_COMPONENTS_HPP
#define MTGL_CONNECTED_COMPONENTS_HPP

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/breadth_first_search.hpp>
#include <mtgl/shiloach_vishkin.hpp>
#include <mtgl/edge_array_adapter.hpp>
#include <mtgl/partitioning.hpp>
#include <mtgl/algorithm.hpp>
#include <mtgl/dynamic_array.hpp>

namespace mtgl {

namespace detail {

template <typename Graph, typename ComponentMap>
class news_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  news_visitor(ComponentMap& res, size_type to_writ) :
    result(res), to_write(to_writ) {}

  void tree_edge(edge_descriptor& e, Graph &g)
  {
    vertex_descriptor v = target(e, g);
    result[v] = to_write;
  }

protected:
  ComponentMap& result;
  size_type to_write;
};

}

/// Search the giant component, then Shiloach-Vishkin on the result.
template <typename Graph, typename ComponentMap>
void connected_components(Graph& g, ComponentMap& result)
{
  #pragma mta noalias g
  #pragma mta noalias result

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_edge_iterator
          thread_edge_iterator;

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

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

  // Find the highest degree vertex.
  size_type* degrees = new size_type[order];

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      degrees[start_pos] = out_degree(*verts, g);
    }
  }

  size_type max_deg_vert_id = 0;
  size_type maxdeg = degrees[0];

  #pragma mta assert nodep
  for (size_type i = 1; i < order; i++)
  {
    if (degrees[i] > maxdeg)
    {
      maxdeg = degrees[i];
      max_deg_vert_id = i;
    }
  }

  delete [] degrees;

  vertex_descriptor max_deg_vert = *(thread_vertices(max_deg_vert_id, g));

  // Search from the highest degree vertex to label the giant component (GCC).

  detail::news_visitor<Graph, ComponentMap> nvis(result, max_deg_vert);
  breadth_first_search(g, max_deg_vert, nvis);

  size_type orphan_edges = 0;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(size, stream_id, num_streams);
    size_type end_pos = end_block_range(size, stream_id, num_streams);

    size_type my_orphan_edges = 0;

    thread_edge_iterator tedgs = thread_edges(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tedgs)
    {
      edge_descriptor e = *tedgs;
      vertex_descriptor u = source(e, g);
      vertex_descriptor v = target(e, g);

      if (result[u] != max_deg_vert_id || result[v] != max_deg_vert_id)
      {
        ++my_orphan_edges;
      }
    }

    mt_incr(orphan_edges, my_orphan_edges);
  }

  size_type next_index = 0;

  size_type* srcs = new size_type[orphan_edges];
  size_type* dests = new size_type[orphan_edges];

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(size, stream_id, num_streams);
    size_type end_pos = end_block_range(size, stream_id, num_streams);

    thread_edge_iterator tedgs = thread_edges(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tedgs)
    {
      edge_descriptor e = *tedgs;
      vertex_descriptor u = source(e, g);
      vertex_descriptor v = target(e, g);

      if (result[u] != max_deg_vert_id || result[v] != max_deg_vert_id)
      {
        size_type my_ind = mt_incr(next_index, 1);
        srcs[my_ind] = get(vid_map, u);
        dests[my_ind] = get(vid_map, v);
      }
    }
  }

  edge_array_adapter<size_type> eaa(srcs, dests, order, orphan_edges);
  vertex_property_map<edge_array_adapter<size_type>, size_type>
    eaa_components(eaa);

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      eaa_components[v] = result[v];
    }
  }

  shiloach_vishkin_no_init(eaa, eaa_components);

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      result[v] = eaa_components[v];
    }
  }

  delete [] srcs;
  delete [] dests;
}

template <typename Graph, typename ComponentMap>
typename graph_traits<Graph>::vertex_descriptor
largest_connected_component(Graph& g, const ComponentMap& C)
{
  #pragma mta noalias g
  #pragma mta noalias C

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type* rcount = (size_type*) malloc(order * sizeof(size_type));
  for (size_type i = 0; i < order; ++i) rcount[i] = 0;

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Count the number of vertices in each component.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      mt_incr(rcount[C[*verts]], 1);
    }
  }

  // Find the id of the largest component.
  size_type largest_comp_id = 0;
  size_type largest_size = 0;

  for (size_type i = 0; i < order; ++i)
  {
    if (rcount[i] > largest_size)
    {
      largest_comp_id = i;
      largest_size = rcount[i];
    }
  }

  free(rcount);

  return largest_comp_id;
}

template <typename Graph, typename ComponentMap>
typename graph_traits<Graph>::size_type
count_connected_components(Graph& g, const ComponentMap& C)
{
  #pragma mta noalias g
  #pragma mta noalias C

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type* rcount = (size_type*) malloc(order * sizeof(size_type));
  for (size_type i = 0; i < order; ++i) rcount[i] = 0;

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Count the number of vertices in each component.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      mt_incr(rcount[C[*verts]], 1);
    }
  }

  size_type num_components = 0;

  // Count the number of components.
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    if (rcount[i] > 0) mt_incr(num_components, 1);
  }

  free(rcount);

  return num_components;
}

template <typename Graph, typename ComponentMap>
void
component_leaders(
    Graph& g, const ComponentMap& C,
    dynamic_array<typename graph_traits<Graph>::vertex_descriptor>& leaders)
{
  #pragma mta noalias g
  #pragma mta noalias C
  #pragma mta noalias leaders

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type* rcount = (size_type*) malloc(order * sizeof(size_type));
  for (size_type i = 0; i < order; ++i) rcount[i] = 0;

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Count the number of vertices in each component.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      mt_incr(rcount[C[*verts]], 1);
    }
  }

  size_type num_components = 0;

  // Count the number of components.
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i)
  {
    if (rcount[i] > 0) mt_incr(num_components, 1);
  }

  // Resize the leaders array.
  leaders.resize(num_components);

  num_components = 0;

  // Put the leaders for each component in the leaders array.
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      if (rcount[start_pos] > 0)
      {
        size_type pos = mt_incr(num_components, 1);
        leaders[pos] = *verts;
      }
    }
  }

  free(rcount);
}

/// The vertex_partition_iterator trades space for O(n log n) time.  Useful
/// if the number of partitions is O(n).
template <typename Graph, typename ComponentMap>
class vertex_partition_iterator {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  vertex_partition_iterator(Graph& gg, ComponentMap& cmps) :
    g(gg), order(num_vertices(g)), C(cmps),
    comp_array((size_type*) malloc(sizeof(size_type) * num_vertices(gg))),
    v_array((vertex_descriptor*) malloc(sizeof(vertex_descriptor) *
                                        num_vertices(gg)))
  {
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
        v_array[start_pos] = v;
        comp_array[start_pos] = C[v];
      }
    }

    sort(comp_array, order, v_array);
    num_parts = 1;

    for (size_type i = 1; i < order; ++i)
    {
      if (comp_array[i] != comp_array[i-1]) mt_incr(num_parts, 1);
    }

    index = (size_type*) malloc(sizeof(size_type) * num_parts);
    index[0] = 0;
    size_type next = 0;

    for (size_type i = 1; i < order; ++i)
    {
      if (comp_array[i] != comp_array[i-1])
      {
        // TODO: These next two lines are NOT threadsafe!
        mt_incr(next, 1);
        index[next] = i;
      }
    }
  }

  ~vertex_partition_iterator()
  {
    free(v_array);
    free(comp_array);
    free(index);
  }

  size_type num_partitions() const { return num_parts; }

  pair<size_type, vertex_descriptor*> operator[](size_type i)
  {
    assert(i < num_parts);

    size_type comp_size = i < num_parts - 1 ? index[i+1] - index[i] :
                                              order - index[i];

    return pair<size_type, vertex_descriptor*>(comp_size, &v_array[index[i]]);
  }

private:
  Graph& g;
  size_type order;
  ComponentMap& C;
  size_type* comp_array;
  vertex_descriptor* v_array;
  size_type* index;
  size_type num_parts;
};

}

#endif
