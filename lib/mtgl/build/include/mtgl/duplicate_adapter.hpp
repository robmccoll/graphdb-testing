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
/*! \file duplicate_adapter.hpp

    \brief This adapter provides an interface to take an existing graph
           adapter and create a new one that has two copies of every edge in
           the original graph.  It stores a reference to the original graph
           and an entire copy of the duplicate graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/2008

    For undirected graphs each edge is just duplicated.  For directed graphs
    each edge is duplicated, but the duplicate's direction is the reverse of
    the original edge's direction.

    The associations between the vertices and edges in the original graph and
    the vertices and edges in the subgraph are stored explicitly.  Subgraph to
    original graph associations are stored using dynamic arrays.  Original
    graph to subgraph associations are stored using hash tables to reduce the
    memory usage while keeping almost constant access times.

    The mapping from global edges to local edges is 1 -> 2.  So, calling
    global_to_local() on an edge returns a pair of edge descriptors.

    This adapter produces deterministic duplicate graphs in structure, but the
    ids of the local edges may change between parallel runnings.

    The base_adapter_type is expected to correctly implement a deep copy
    for both the copy constructor and the assignment operator.
*/
/****************************************************************************/

#ifndef MTGL_DUPLICATE_ADAPTER_HPP
#define MTGL_DUPLICATE_ADAPTER_HPP

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>

namespace mtgl {

template <typename Graph>
class duplicate_adapter {
public:
  typedef graph_traits<Graph> base_traits;
  typedef typename base_traits::size_type base_size_type;
  typedef typename base_traits::vertex_descriptor base_vertex_descriptor;
  typedef typename base_traits::edge_descriptor base_edge_descriptor;
  typedef typename base_traits::vertex_iterator base_vertex_iterator;
  typedef typename base_traits::edge_iterator base_edge_iterator;
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

  duplicate_adapter(Graph& g) :
    original_graph(&g),
    m_local_vertex(static_cast<size_type>(1.7 * num_vertices(g)))
  {
    base_size_type order = num_vertices(*original_graph);
    base_size_type size = num_edges(*original_graph);
    base_size_type dsize = 2 * size;

    base_vertex_iterator verts = vertices(*original_graph);
    base_edge_iterator edgs = edges(*original_graph);

    m_global_vertex.resize(order);
    m_global_edge.resize(dsize);

    size_type* sources = new size_type[size * 2];
    size_type* dests = new size_type[size * 2];

    #pragma mta assert parallel 
    for (size_type i = 0; i < order; ++i)
    {
      base_vertex_descriptor v = verts[i]; 
      m_local_vertex.insert(v, i); 
      m_global_vertex[i] = v;
    }

    // Get edge source and target ids.
    #pragma mta assert parallel
    for (base_size_type i = 0; i < size; ++i)
    {
      base_edge_descriptor e = edgs[i];

      // Store the first copy of the original edge.
      m_local_vertex.lookup(source(e, *original_graph), sources[i]);
      m_local_vertex.lookup(target(e, *original_graph), dests[i]);

      // Store the duplicate copy of the original edge.
      if (is_undirected(g))
      {
        sources[i + size] = sources[i];
        dests[i + size] = dests[i];
      }
      else
      {
        sources[i + size] = dests[i];
        dests[i + size] = sources[i];
      }
    }

    // Initialize the duplicate graph.  This is assuming that the init()
    // method of the underlying graph implements parallelization correctly
    // and efficiently.
    init(order, dsize, sources, dests, duplicate_graph);

    #pragma mta assert parallel 
    for (size_type i = 0; i < dsize; ++i)
    {
      if (i < size)
      {
        base_edge_descriptor e = edgs[i];
        m_local_edge.insert(e, i); 
        m_global_edge[i] = e;
      }
      else
      {
        m_global_edge[i] = edgs[i - size];
      }
    }

    delete [] sources;
    delete [] dests;
  }

  duplicate_adapter(const duplicate_adapter& dg) { deep_copy(dg); }

  duplicate_adapter& operator=(const duplicate_adapter& rhs)
  {
    deep_copy(rhs);

    return *this;
  }

  const wrapper_adapter& get_adapter() const { return duplicate_graph; }

  void print()
  {
    mtgl::print(duplicate_graph);

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

      std::cout << "  " << get(eid_map, m_global_edge[i]) << " -> ("
                << local << ", " << local + mle_size << ")" << std::endl;
    }
  }

  vertex_descriptor
  global_to_local(const base_vertex_descriptor& u_global) const
  {
    // We can make the interchange of size_type and vertex_descriptor because
    // they are the same type for compressed_sparse_row_graph.
    vertex_descriptor u_local = null_vertex(duplicate_graph);
    m_local_vertex.lookup(u_global, u_local);
    return u_local;
  }

  base_vertex_descriptor local_to_global(const vertex_descriptor& u_local) const
  {
    return m_global_vertex[get(get(_vertex_id_map, duplicate_graph), u_local)];
  }

  pair<edge_descriptor, edge_descriptor>
  global_to_local(const base_edge_descriptor& e_global) const
  {
    size_type e_local_id;
    if (m_local_edge.lookup(e_global, e_local_id))
    {
      edge_iterator edgs = edges(duplicate_graph);
      return pair<edge_descriptor, edge_descriptor>(
               edgs[e_local_id],
               edgs[e_local_id + num_edges(duplicate_graph) / 2]);
    }
    else
    {
      return pair<edge_descriptor, edge_descriptor>(
               null_edge(duplicate_graph), null_edge(duplicate_graph));
    }
  }

  base_edge_descriptor local_to_global(const edge_descriptor& e_local) const
  {
    return m_global_edge[get(get(_edge_id_map, duplicate_graph), e_local)];
  }

private:
  void deep_copy(const duplicate_adapter& rhs)
  {
    original_graph = rhs.original_graph;
    duplicate_graph = rhs.duplicate_graph;

    // Since the local <-> global links use ids for the duplicate graph
    // vertices and edges, simply copying the data structures will work.  The
    // original graph vertex and edge descriptors are fine because the copied
    // graph will point to the same original graph.
    m_global_vertex = rhs.m_global_vertex;
    m_global_edge = rhs.m_global_edge;
    m_local_vertex = rhs.m_local_vertex;
    m_local_edge = rhs.m_local_edge;
  }

private:
  Graph* original_graph;
  wrapper_adapter duplicate_graph; 

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
typename duplicate_adapter<G>::size_type
num_vertices(const duplicate_adapter<G>& dg)
{
  return num_vertices(dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::size_type
num_edges(const duplicate_adapter<G>& dg)
{
  return num_edges(dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::vertex_descriptor
source(const typename duplicate_adapter<G>::edge_descriptor& e,
       const duplicate_adapter<G>& dg)
{
  return source(e, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::vertex_descriptor
target(const typename duplicate_adapter<G>::edge_descriptor& e,
       const duplicate_adapter<G>& dg)
{
  return target(e, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::size_type
degree(const typename duplicate_adapter<G>::vertex_descriptor& u_local,
       const duplicate_adapter<G>& dg)
{
  return degree(u_local, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::size_type
out_degree(const typename duplicate_adapter<G>::vertex_descriptor& u_local,
           const duplicate_adapter<G>& dg)
{
  return out_degree(u_local, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::size_type
in_degree(const typename duplicate_adapter<G>::vertex_descriptor& u_local,
          const duplicate_adapter<G>& dg)
{
  return in_degree(u_local, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::vertex_iterator
vertices(const duplicate_adapter<G>& dg)
{
  return vertices(dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::edge_iterator
edges(const duplicate_adapter<G>& dg)
{
  return edges(dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::adjacency_iterator
adjacent_vertices(const typename duplicate_adapter<G>::vertex_descriptor& v,
                  const duplicate_adapter<G>& dg)
{
  return adjacent_vertices(v, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::in_adjacency_iterator
in_adjacent_vertices(const typename duplicate_adapter<G>::vertex_descriptor& v,
                     const duplicate_adapter<G>& dg)
{
  return in_adjacent_vertices(v, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::out_edge_iterator
out_edges(const typename duplicate_adapter<G>::vertex_descriptor& v,
          const duplicate_adapter<G>& dg)
{
  return out_edges(v, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::in_edge_iterator
in_edges(const typename duplicate_adapter<G>::vertex_descriptor& v,
         const duplicate_adapter<G>& dg)
{
  return in_edges(v, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::thread_vertex_iterator
thread_vertices(typename duplicate_adapter<G>::size_type pos,
                const duplicate_adapter<G>& dg)
{
  return thread_vertices(pos, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::thread_adjacency_iterator
thread_adjacent_vertices(
    const typename duplicate_adapter<G>::vertex_descriptor& v,
    typename duplicate_adapter<G>::size_type pos,
    const duplicate_adapter<G>& dg)
{
  return thread_adjacent_vertices(v, pos, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::thread_in_adjacency_iterator
thread_in_adjacent_vertices(
    const typename duplicate_adapter<G>::vertex_descriptor& v,
    typename duplicate_adapter<G>::size_type pos,
    const duplicate_adapter<G>& dg)
{
  return thread_in_adjacent_vertices(v, pos, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::thread_edge_iterator
thread_edges(typename duplicate_adapter<G>::size_type pos,
             const duplicate_adapter<G>& dg)
{
  return thread_edges(pos, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::thread_out_edge_iterator
thread_out_edges(const typename duplicate_adapter<G>::vertex_descriptor& v,
                 typename duplicate_adapter<G>::size_type pos,
                 const duplicate_adapter<G>& dg)
{
  return thread_out_edges(v, pos, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::thread_in_edge_iterator
thread_in_edges(const typename duplicate_adapter<G>::vertex_descriptor& v,
                typename duplicate_adapter<G>::size_type pos,
                const duplicate_adapter<G>& dg)
{
  return thread_in_edges(v, pos, dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::vertex_descriptor
null_vertex(const duplicate_adapter<G>& dg)
{
  return null_vertex(dg.get_adapter());
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::edge_descriptor
null_edge(const duplicate_adapter<G>& dg)
{
  return null_edge(dg.get_adapter());
}

/***/

template <typename ITERATOR, typename G>
inline
bool
is_valid(ITERATOR& iter, typename duplicate_adapter<G>::size_type p,
         const duplicate_adapter<G>& dg)
{
  return is_valid(iter, p, dg.get_adapter());
}

/***/

template <typename G>
inline
bool
is_directed(const duplicate_adapter<G>& dg)
{
  return is_directed(dg.get_adapter());
}

/***/

template <typename G>
inline
bool
is_undirected(const duplicate_adapter<G>& dg)
{
  return is_undirected(dg.get_adapter());
}

/***/

template <typename G>
inline
bool
is_bidirectional(const duplicate_adapter<G>& dg)
{
  return is_bidirectional(dg.get_adapter());
}

/***/

template <typename G>
class vertex_id_map<duplicate_adapter<G> > :
  public vertex_id_map<typename duplicate_adapter<G>::wrapper_adapter> {
public:
  typedef typename duplicate_adapter<G>::wrapper_adapter DG;

  vertex_id_map() : vertex_id_map<DG>() {}
  vertex_id_map(const vertex_id_map& vm) : vertex_id_map<DG>(vm) {}
};

/***/

template <typename G>
class edge_id_map<duplicate_adapter<G> > :
  public edge_id_map<typename duplicate_adapter<G>::wrapper_adapter> {
public:
  typedef typename duplicate_adapter<G>::wrapper_adapter DG;

  edge_id_map() : edge_id_map<DG>() {}
  edge_id_map(const edge_id_map& em) : edge_id_map<DG>(em) {}
};

/***/

template <typename G>
inline
typename duplicate_adapter<G>::vertex_descriptor
global_to_local(const typename duplicate_adapter<G>::base_vertex_descriptor& v,
                const duplicate_adapter<G>& dg)
{
  return dg.global_to_local(v);
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::base_vertex_descriptor
local_to_global(const typename duplicate_adapter<G>::vertex_descriptor& v,
                const duplicate_adapter<G>& dg)
{
  return dg.local_to_global(v);
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::edge_descriptor
global_to_local(const typename duplicate_adapter<G>::base_edge_descriptor& e,
                const duplicate_adapter<G>& dg)
{
  return dg.global_to_local(e);
}

/***/

template <typename G>
inline
typename duplicate_adapter<G>::base_edge_descriptor
local_to_global(const typename duplicate_adapter<G>::edge_descriptor& e,
                const duplicate_adapter<G>& dg)
{
  return dg.local_to_global(e);
}

}

#endif
