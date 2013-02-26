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
/*! \file subgraph_adapter.hpp

    \brief This adapter provides an interface to create a subgraph of
           an existing graph adapter.  It stores a reference to the original
           graph and an entire copy of the subgraph.  The subgraph could
           possibly be stored more efficiently, but storing a full copy of the
           subgraph allows this to work with any graph adapter that follows
           the standard interface.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/17/2008

    The associations between the vertices and edges in the original graph and
    the vertices and edges in the subgraph are stored explicitly.  Subgraph to
    original graph associations are stored using dynamic arrays.  Original
    graph to subgraph associations are stored using hash tables to reduce the
    memory usage while keeping almost constant access times.

    This adapter produces deterministic subgraphs in structure, but the ids
    of the local edges may change between parallel runnings.

    The base_adapter_type is expected to correctly implement a deep copy
    for both the copy constructor and the assignment operator.
*/
/****************************************************************************/

#ifndef MTGL_SUBGRAPH_ADAPTER_HPP
#define MTGL_SUBGRAPH_ADAPTER_HPP

#include <iostream>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>

namespace mtgl {

template <typename Graph>
class subgraph_adapter {
public:
  typedef graph_traits<Graph> base_traits;
  typedef typename base_traits::size_type base_size_type;
  typedef typename base_traits::vertex_descriptor base_vertex_descriptor;
  typedef typename base_traits::edge_descriptor base_edge_descriptor;
  typedef typename base_traits::vertex_iterator base_vertex_iterator;
  typedef typename base_traits::adjacency_iterator base_adjacency_iterator;
  typedef typename base_traits::edge_iterator base_edge_iterator;
  typedef typename base_traits::out_edge_iterator base_out_edge_iterator;
  typedef typename base_traits::directed_category base_directed_category;

  typedef compressed_sparse_row_graph<base_directed_category> wrapper_adapter;
  typedef graph_traits<wrapper_adapter> traits;
  typedef typename traits::size_type size_type;
  typedef typename traits::vertex_descriptor vertex_descriptor;
  typedef typename traits::edge_descriptor edge_descriptor;
  typedef typename traits::vertex_iterator vertex_iterator;
  typedef typename traits::adjacency_iterator adjacency_iterator;
  typedef typename traits::in_adjacency_iterator in_adjacency_iterator;
  typedef typename traits::edge_iterator edge_iterator;
  typedef typename traits::out_edge_iterator out_edge_iterator;
  typedef typename traits::in_edge_iterator in_edge_iterator;
  typedef typename traits::thread_vertex_iterator thread_vertex_iterator;
  typedef typename traits::thread_adjacency_iterator thread_adjacency_iterator;
  typedef typename traits::thread_in_adjacency_iterator
          thread_in_adjacency_iterator;
  typedef typename traits::thread_edge_iterator thread_edge_iterator;
  typedef typename traits::thread_out_edge_iterator thread_out_edge_iterator;
  typedef typename traits::thread_in_edge_iterator thread_in_edge_iterator;
  typedef typename traits::directed_category directed_category;
  typedef typename traits::iterator_category iterator_category;

  subgraph_adapter(Graph& g) : original_graph(&g) {}
  subgraph_adapter(const subgraph_adapter& sg) { deep_copy(sg); }

  subgraph_adapter& operator=(const subgraph_adapter& rhs)
  {
    deep_copy(rhs);

    return *this;
  }

  const wrapper_adapter& get_adapter() const { return subgraph; }

  vertex_descriptor
  global_to_local(const base_vertex_descriptor& u_global) const
  {
    // We can make the interchange of size_type and vertex_descriptor because
    // they are the same type for compressed_sparse_row_graph.
    vertex_descriptor u_local = null_vertex(subgraph);
    m_local_vertex.lookup(u_global, u_local);
    return u_local;
  }

  base_vertex_descriptor local_to_global(const vertex_descriptor& u_local) const
  {
    return m_global_vertex[get(get(_vertex_id_map, subgraph), u_local)];
  }

  edge_descriptor global_to_local(const base_edge_descriptor& e_global) const
  {
    // We can make the interchange of size_type and vertex_descriptor because
    // they are the same type for compressed_sparse_row_graph.
    edge_descriptor e_local = null_edge(subgraph);
    m_local_edge.lookup(e_global, e_local);
    return e_local;
  }

  base_edge_descriptor local_to_global(const edge_descriptor& e_local) const
  {
    return m_global_edge[get(get(_edge_id_map, subgraph), e_local)];
  }

  void print()
  {
    mtgl::print(subgraph);

    vertex_id_map<Graph> vid_map = get(_vertex_id_map, *original_graph);
    edge_id_map<Graph> eid_map = get(_edge_id_map, *original_graph);

    size_type mgv_size = m_global_vertex.size();
    std::cout << std::endl
              << "m_global_vertex: " << mgv_size << std::endl;
    for (size_type i = 0; i < mgv_size; ++i)
    {
      std::cout << "  " << i << " -> " << get(vid_map, m_global_vertex[i])
                << std::endl;
    }

    size_type mlv_size = m_local_vertex.size();
    std::cout << std::endl
              << "m_local_vertex: " << mlv_size << std::endl;
    for (size_type i = 0; i < mlv_size; ++i)
    {
      size_type local = 0;
      m_local_vertex.lookup(m_global_vertex[i], local);

      std::cout << "  " << get(vid_map, m_global_vertex[i]) << " -> "
                << local << std::endl;
    }

    size_type mge_size = m_global_edge.size();
    std::cout << std::endl
              << "m_global_edge: " << mge_size << std::endl;
    for (size_type i = 0; i < mge_size; ++i)
    {
      std::cout << "  " << i << " -> " << get(eid_map, m_global_edge[i])
                << std::endl;
    }

    size_type mle_size = m_local_edge.size();
    std::cout << std::endl
              << "m_local_edge: " << mle_size << std::endl;
    for (size_type i = 0; i < mle_size; ++i)
    {
      size_type local = 0;
      m_local_edge.lookup(m_global_edge[i], local);

      std::cout << "  " << get(eid_map, m_global_edge[i]) << " -> "
                << local << std::endl;
    }
  }

  template <typename HT>
  class uniqueid_visitor {
  public:
    uniqueid_visitor(dynamic_array<typename HT::key_type>& g) :
      count(0), global(g) {}

    void operator()(typename HT::value_type* i)
    {
      i->second = mt_incr(count, 1);
      global[i->second] = i->first;
    }

  private:
    size_type count;
    dynamic_array<typename HT::key_type>& global;
  };

  // Creates an edge induced subgraph.  We assume that emask is the size of
  // original_graph.
  template <typename EdgeMaskMap>
  void internal_init_edges(EdgeMaskMap& emask)
  {
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, *original_graph);

    base_edge_iterator edgs = edges(*original_graph);
    base_size_type num_global_edges = num_edges(*original_graph);

    // Count the number of edges.
    size_type num_local_edges = 0;
    for (size_type i = 0; i < num_global_edges; ++i)
    {
      num_local_edges += emask[edgs[i]];
    }

    m_local_vertex.resize(2 * num_local_edges);
    m_local_edge.resize(static_cast<size_type>(1.7 * num_local_edges));

    // Create the unique sets of vertices and edges and set them as the
    // local vertices and edges.
    #pragma mta assert parallel
    for (size_type i = 0; i < num_global_edges; ++i)
    {
      if (emask[edgs[i]])
      {
        base_vertex_descriptor src = source(edgs[i], *original_graph);
        base_vertex_descriptor dest = target(edgs[i], *original_graph);

        // Add the source and destination vertex ids to m_local_vertex.
        if (get(vid_map, src) < get(vid_map, dest)) {
          m_local_vertex.insert(src, 0);
          m_local_vertex.insert(dest, 0);
        } else {
          m_local_vertex.insert(dest, 0);
          m_local_vertex.insert(src, 0);
        }

        // Add the edge ids to m_local_edge.
        m_local_edge.insert(edgs[i], 0);
      }
    }

    size_type num_unique_verts = m_local_vertex.size();
    size_type num_unique_edges = m_local_edge.size();

    // Create unique ids for the local vertices, and initialize
    // m_global_vertex.
    m_global_vertex.resize(num_unique_verts);
    uniqueid_visitor<xmt_hash_table<base_vertex_descriptor, size_type> >
        v_uidvis(m_global_vertex);
    m_local_vertex.visit(v_uidvis);

    // Create unique ids for the local edges, and initialize m_global_edge.
    m_global_edge.resize(num_unique_edges);
    uniqueid_visitor<xmt_hash_table<base_edge_descriptor, size_type> >
        e_uidvis(m_global_edge);
    m_local_edge.visit(e_uidvis);

    // Create the arrays to hold the local edge sources and dests.  These will
    // be passed to subgraph's init().
    size_type* sources = new size_type[num_unique_edges];
    size_type* dests = new size_type[num_unique_edges];

    #pragma mta assert parallel
    for (size_type i = 0; i < num_unique_edges; ++i)
    {
      base_vertex_descriptor src = source(m_global_edge[i], *original_graph);
      base_vertex_descriptor dest = target(m_global_edge[i], *original_graph);

      m_local_vertex.lookup(src, sources[i]);
      m_local_vertex.lookup(dest, dests[i]);
    }

    init(num_unique_verts, num_unique_edges, sources, dests, subgraph);

    delete [] sources;
    delete [] dests;
  }

  // Creates an edge induced subgraph.
  void internal_init_edges(base_size_type numEdges,
                           base_edge_descriptor* e_globals)
  {
    m_local_vertex.resize(2 * numEdges);
    m_local_edge.resize(static_cast<size_type>(1.7 * numEdges));

    // Create the unique sets of vertices and edges and set them as the
    // local vertices and edges.
    #pragma mta assert parallel
    for (base_size_type i = 0; i < numEdges; ++i)
    {
      base_vertex_descriptor src = source(e_globals[i], *original_graph);
      base_vertex_descriptor dest = target(e_globals[i], *original_graph);

      // Add the source and destination vertex ids to m_local_vertex.
      m_local_vertex.insert(src, 0);
      m_local_vertex.insert(dest, 0);

      // Add the edge ids to m_local_edge.
      m_local_edge.insert(e_globals[i], 0);
    }

    size_type num_unique_verts = m_local_vertex.size();
    size_type num_unique_edges = m_local_edge.size();

    // Create unique ids for the local vertices, and initialize
    // m_global_vertex.
    m_global_vertex.resize(num_unique_verts);
    uniqueid_visitor<xmt_hash_table<base_vertex_descriptor, size_type> >
        v_uidvis(m_global_vertex);
    m_local_vertex.visit(v_uidvis);

    // Create unique ids for the local edges, and initialize m_global_edge.
    m_global_edge.resize(num_unique_edges);
    uniqueid_visitor<xmt_hash_table<base_edge_descriptor, size_type> >
        e_uidvis(m_global_edge);
    m_local_edge.visit(e_uidvis);

    // Create the arrays to hold the local edge sources and dests.  These will
    // be passed to subgraph's init().
    size_type* sources = new size_type[num_unique_edges];
    size_type* dests = new size_type[num_unique_edges];

    #pragma mta assert parallel
    for (size_type i = 0; i < num_unique_edges; ++i)
    {
      base_vertex_descriptor src = source(m_global_edge[i], *original_graph);
      base_vertex_descriptor dest = target(m_global_edge[i], *original_graph);

      m_local_vertex.lookup(src, sources[i]);
      m_local_vertex.lookup(dest, dests[i]);
    }

    init(num_unique_verts, num_unique_edges, sources, dests, subgraph);

    delete [] sources;
    delete [] dests;
  }

  // Creates a vertex induced subgraph.
  template <typename VertexMaskMap>
  void internal_init_vertices(VertexMaskMap& vmask)
  {
    base_vertex_iterator verts = vertices(*original_graph);
    base_size_type num_global_vertices = num_vertices(*original_graph);

    // Count the number of vertices.
    size_type num_local_vertices = 0;
    for (base_size_type i = 0; i < num_global_vertices; ++i)
    {
      num_local_vertices += vmask[verts[i]];
    }

    m_local_vertex.resize(static_cast<size_type>(1.7 * num_local_vertices));

    // Create the unique set of vertices and set them as the local vertices.
    #pragma mta assert parallel
    for (base_size_type i = 0; i < num_global_vertices; ++i)
    {
      // Add this vertex if it is in the vmask.
      base_vertex_descriptor v = verts[i];
      if (vmask[v]) m_local_vertex.insert(v, 0);
    }

    // Count the number of unique edges.  Only count edges that have both
    // endpoints in v_globals.
    size_type num_local_edges = 0;
    for (base_size_type i = 0; i < num_global_vertices; ++i)
    {
      base_vertex_descriptor u = verts[i];

      if (vmask[u])
      {
        base_adjacency_iterator adj_iter =
          adjacent_vertices(u, *original_graph);
        base_size_type out_deg = out_degree(u, *original_graph);

        for (base_size_type j = 0; j < out_deg; ++j)
        {
          num_local_edges += vmask[adj_iter[j]];
        }
      }
    }

    m_local_edge.resize(static_cast<size_type>(1.7 * num_local_edges));

    // Create the unique set of edges and set them as the local edges.  Only
    // add edges that have both endpoints in v_globals.
    #pragma mta assert parallel
    for (base_size_type i = 0; i < num_global_vertices; ++i)
    {
      base_vertex_descriptor u = verts[i];

      if (vmask[u])
      {
        base_out_edge_iterator oe_iter = out_edges(u, *original_graph);
        base_size_type out_deg = out_degree(u, *original_graph);

        #pragma mta assert parallel
        for (size_type j = 0; j < out_deg; ++j)
        {
          base_edge_descriptor e = oe_iter[j];

          // Only add this edge if its target is in vmask.  We already know
          // that the source is in vmask.
          if (vmask[target(e, *original_graph)]) m_local_edge.insert(e, 0);
        }
      }
    }

    size_type num_unique_verts = m_local_vertex.size();
    size_type num_unique_edges = m_local_edge.size();

    // Create unique ids for the local vertices, and initialize
    // m_global_vertex.
    m_global_vertex.resize(num_unique_verts);
    uniqueid_visitor<xmt_hash_table<base_vertex_descriptor, size_type> >
        v_uidvis(m_global_vertex);
    m_local_vertex.visit(v_uidvis);

    // Create unique ids for the local edges, and initialize m_global_edge.
    m_global_edge.resize(num_unique_edges);
    uniqueid_visitor<xmt_hash_table<base_edge_descriptor, size_type> >
        e_uidvis(m_global_edge);
    m_local_edge.visit(e_uidvis);

    // Create the arrays to hold the local edge sources and dests.  These will
    // be passed to subgraph's init().
    size_type* sources = new size_type[num_unique_edges];
    size_type* dests = new size_type[num_unique_edges];

    #pragma mta assert parallel
    for (size_type i = 0; i < num_unique_edges; ++i)
    {
      base_vertex_descriptor src = source(m_global_edge[i], *original_graph);
      base_vertex_descriptor dest = target(m_global_edge[i], *original_graph);

      m_local_vertex.lookup(src, sources[i]);
      m_local_vertex.lookup(dest, dests[i]);
    }

    init(num_unique_verts, num_unique_edges, sources, dests, subgraph);

    delete [] sources;
    delete [] dests;
  }

  // Creates a vertex induced subgraph.
  void internal_init_vertices(base_size_type numVertices,
                              base_vertex_descriptor* v_globals)
  {
    m_local_vertex.resize(static_cast<size_type>(1.7 * numVertices));

    // Create the unique set of vertices and set them as the local vertices.
    #pragma mta assert parallel
    for (base_size_type i = 0; i < numVertices; ++i)
    {
      // Add the vertex id to m_local_vertex.
      m_local_vertex.insert(v_globals[i], 0);
    }

    // Count the number of unique edges.  Only count edges that have both
    // endpoints in v_globals.
    size_type num_local_edges = 0;
    #pragma mta assert parallel
    for (size_type i = 0; i < numVertices; ++i)
    {
      base_adjacency_iterator adj_iter = adjacent_vertices(v_globals[i],
                                                           *original_graph);
      base_size_type out_deg = out_degree(v_globals[i], *original_graph);

      #pragma mta assert parallel
      for (base_size_type j = 0; j < out_deg; ++j)
      {
        num_local_edges += m_local_vertex.member(adj_iter[j]);
      }
    }

    m_local_edge.resize(static_cast<size_type>(1.7 * num_local_edges));

    // Create the unique set of edges and set them as the local edges.  Only
    // add edges that have both endpoints in v_globals.
    #pragma mta assert parallel
    for (base_size_type i = 0; i < numVertices; ++i)
    {
      base_out_edge_iterator oe_iter = out_edges(v_globals[i], *original_graph);
      base_size_type out_deg = out_degree(v_globals[i], *original_graph);

      #pragma mta assert parallel
      for (base_size_type j = 0; j < out_deg; ++j)
      {
        base_edge_descriptor e = oe_iter[j];

        // Only add this edge if its target is in v_globals.  We already know
        // that the source is in v_globals.
        if (m_local_vertex.member(target(e, *original_graph)))
        {
          m_local_edge.insert(e, 0);
        }
      }
    }

    size_type num_unique_verts = m_local_vertex.size();
    size_type num_unique_edges = m_local_edge.size();

    // Create unique ids for the local vertices, and initialize
    // m_global_vertex.
    m_global_vertex.resize(num_unique_verts);
    uniqueid_visitor<xmt_hash_table<base_vertex_descriptor, size_type> >
        v_uidvis(m_global_vertex);
    m_local_vertex.visit(v_uidvis);

    // Create unique ids for the local edges, and initialize m_global_edge.
    m_global_edge.resize(num_unique_edges);
    uniqueid_visitor<xmt_hash_table<base_edge_descriptor, size_type> >
        e_uidvis(m_global_edge);
    m_local_edge.visit(e_uidvis);

    // Create the arrays to hold the local edge sources and dests.  These will
    // be passed to subgraph's init().
    size_type* sources = new size_type[num_unique_edges];
    size_type* dests = new size_type[num_unique_edges];

    #pragma mta assert parallel
    for (size_type i = 0; i < num_unique_edges; ++i)
    {
      base_vertex_descriptor src = source(m_global_edge[i], *original_graph);
      base_vertex_descriptor dest = target(m_global_edge[i], *original_graph);

      m_local_vertex.lookup(src, sources[i]);
      m_local_vertex.lookup(dest, dests[i]);
    }

    init(num_unique_verts, num_unique_edges, sources, dests, subgraph);

    delete [] sources;
    delete [] dests;
  }

private:
  void deep_copy(const subgraph_adapter& rhs)
  {
    original_graph = rhs.original_graph;
    subgraph = rhs.subgraph;

    // Since the local <-> global links use ids for the subgraph vertices and
    // edges, simply copying the data structures will work.  The original
    // graph vertex and edge descriptors are fine because the copied graph
    // will point to the same original graph.
    m_global_vertex = rhs.m_global_vertex;
    m_global_edge = rhs.m_global_edge;
    m_local_vertex = rhs.m_local_vertex;
    m_local_edge = rhs.m_local_edge;
  }

private:
  Graph* original_graph;
  wrapper_adapter subgraph;

  // local -> global
  dynamic_array<base_vertex_descriptor> m_global_vertex;
  dynamic_array<base_edge_descriptor> m_global_edge;

  // global -> local
  xmt_hash_table<base_vertex_descriptor, size_type> m_local_vertex;
  xmt_hash_table<base_edge_descriptor, size_type> m_local_edge;
};

/***/

template <typename G>
inline
typename subgraph_adapter<G>::size_type
num_vertices(const subgraph_adapter<G>& sg)
{
  return num_vertices(sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::size_type
num_edges(const subgraph_adapter<G>& sg)
{
  return num_edges(sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::vertex_descriptor
source(const typename subgraph_adapter<G>::edge_descriptor& e,
       const subgraph_adapter<G>& sg)
{
  return source(e, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::vertex_descriptor
target(const typename subgraph_adapter<G>::edge_descriptor& e,
       const subgraph_adapter<G>& sg)
{
  return target(e, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::size_type
degree(const typename subgraph_adapter<G>::vertex_descriptor& u_local,
       const subgraph_adapter<G>& sg)
{
  return degree(u_local, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::size_type
out_degree(const typename subgraph_adapter<G>::vertex_descriptor& u_local,
           const subgraph_adapter<G>& sg)
{
  return out_degree(u_local, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::size_type
in_degree(const typename subgraph_adapter<G>::vertex_descriptor& u_local,
          const subgraph_adapter<G>& sg)
{
  return in_degree(u_local, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::vertex_iterator
vertices(const subgraph_adapter<G>& sg)
{
  return vertices(sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::edge_iterator
edges(const subgraph_adapter<G>& sg)
{
  return edges(sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::adjacency_iterator
adjacent_vertices(const typename subgraph_adapter<G>::vertex_descriptor& v,
                  const subgraph_adapter<G>& sg)
{
  return adjacent_vertices(v, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::in_adjacency_iterator
in_adjacent_vertices(const typename subgraph_adapter<G>::vertex_descriptor& v,
                     const subgraph_adapter<G>& sg)
{
  return in_adjacent_vertices(v, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::out_edge_iterator
out_edges(const typename subgraph_adapter<G>::vertex_descriptor& v,
          const subgraph_adapter<G>& sg)
{
  return out_edges(v, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::in_edge_iterator
in_edges(const typename subgraph_adapter<G>::vertex_descriptor& v,
         const subgraph_adapter<G>& sg)
{
  return in_edges(v, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::thread_vertex_iterator
thread_vertices(typename subgraph_adapter<G>::size_type pos,
                const subgraph_adapter<G>& sg)
{
  return thread_vertices(pos, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::thread_adjacency_iterator
thread_adjacent_vertices(
    const typename subgraph_adapter<G>::vertex_descriptor& v,
    typename subgraph_adapter<G>::size_type pos,
    const subgraph_adapter<G>& sg)
{
  return thread_adjacent_vertices(v, pos, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::thread_in_adjacency_iterator
thread_in_adjacent_vertices(
    const typename subgraph_adapter<G>::vertex_descriptor& v,
    typename subgraph_adapter<G>::size_type pos,
    const subgraph_adapter<G>& sg)
{
  return thread_in_adjacent_vertices(v, pos, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::thread_edge_iterator
thread_edges(typename subgraph_adapter<G>::size_type pos,
             const subgraph_adapter<G>& sg)
{
  return thread_edges(pos, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::thread_out_edge_iterator
thread_out_edges(const typename subgraph_adapter<G>::vertex_descriptor& v,
                 typename subgraph_adapter<G>::size_type pos,
                 const subgraph_adapter<G>& sg)
{
  return thread_out_edges(v, pos, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::thread_in_edge_iterator
thread_in_edges(const typename subgraph_adapter<G>::vertex_descriptor& v,
                typename subgraph_adapter<G>::size_type pos,
                const subgraph_adapter<G>& sg)
{
  return thread_in_edges(v, pos, sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::vertex_descriptor
null_vertex(const subgraph_adapter<G>& sg)
{
  return null_vertex(sg.get_adapter());
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::edge_descriptor
null_edge(const subgraph_adapter<G>& sg)
{
  return null_edge(sg.get_adapter());
}

/***/

template <typename ITERATOR, typename G>
inline
bool
is_valid(ITERATOR& iter, typename subgraph_adapter<G>::size_type p,
         const subgraph_adapter<G>& sg)
{
  return is_valid(iter, p, sg.get_adapter());
}

/***/

template <typename G>
inline
bool
is_directed(const subgraph_adapter<G>& sg)
{
  return is_directed(sg.get_adapter());
}

/***/

template <typename G>
inline
bool
is_undirected(const subgraph_adapter<G>& sg)
{
  return is_undirected(sg.get_adapter());
}

/***/

template <typename G>
inline
bool
is_bidirectional(const subgraph_adapter<G>& sg)
{
  return is_bidirectional(sg.get_adapter());
}

/***/

template <typename G>
class vertex_id_map<subgraph_adapter<G> > :
  public vertex_id_map<typename subgraph_adapter<G>::wrapper_adapter> {
public:
  typedef typename subgraph_adapter<G>::wrapper_adapter SG;

  vertex_id_map() : vertex_id_map<SG>() {}
  vertex_id_map(const vertex_id_map& vm) : vertex_id_map<SG>(vm) {}
};

/***/

template <typename G>
class edge_id_map<subgraph_adapter<G> > :
  public edge_id_map<typename subgraph_adapter<G>::wrapper_adapter> {
public:
  typedef typename subgraph_adapter<G>::wrapper_adapter SG;

  edge_id_map() : edge_id_map<SG>() {}
  edge_id_map(const edge_id_map& em) : edge_id_map<SG>(em) {}
};

/***/

template <typename G, typename EdgeMaskMap>
inline
void
init_edges(EdgeMaskMap& emask, subgraph_adapter<G>& sg)
{
  sg.internal_init_edges(emask);
}

/***/

template <typename G>
inline
void
init_edges(typename subgraph_adapter<G>::base_edge_descriptor* e_globals,
           typename subgraph_adapter<G>::size_type numEdges,
           subgraph_adapter<G>& sg)
{
  sg.internal_init_edges(numEdges, e_globals);
}

/***/

template <typename G, typename VertexMaskMap>
inline
void
init_vertices(VertexMaskMap& vmask, subgraph_adapter<G>& sg)
{
  sg.internal_init_vertices(vmask);
}

/***/

template <typename G>
inline
void
init_vertices(typename subgraph_adapter<G>::base_vertex_descriptor* v_globals,
              typename subgraph_adapter<G>::base_size_type numVertices,
              subgraph_adapter<G>& sg)
{
  sg.internal_init_vertices(numVertices, v_globals);
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::vertex_descriptor
global_to_local(const typename subgraph_adapter<G>::base_vertex_descriptor& v,
                const subgraph_adapter<G>& sg)
{
  return sg.global_to_local(v);
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::base_vertex_descriptor
local_to_global(const typename subgraph_adapter<G>::vertex_descriptor& v,
                const subgraph_adapter<G>& sg)
{
  return sg.local_to_global(v);
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::edge_descriptor
global_to_local(const typename subgraph_adapter<G>::base_edge_descriptor& e,
                const subgraph_adapter<G>& sg)
{
  return sg.global_to_local(e);
}

/***/

template <typename G>
inline
typename subgraph_adapter<G>::base_edge_descriptor
local_to_global(const typename subgraph_adapter<G>::edge_descriptor& e,
                const subgraph_adapter<G>& sg)
{
  return sg.local_to_global(e);
}

}

#endif
