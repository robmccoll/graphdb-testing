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
/*! \file stinger_graph_adapter.hpp

    \brief The graph adapter wrapper for stinger_graph.hpp.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 5/18/2011
*/
/****************************************************************************/

#ifndef MTGL_STINGER_GRAPH_ADAPTER_HPP
#define MTGL_STINGER_GRAPH_ADAPTER_HPP

#include <mtgl/stinger_graph.hpp>
#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/util.hpp>

namespace mtgl {

namespace detail {

template <typename stinger_graph>
class stinger_edge_adapter {
private:
  typedef typename stinger_graph::edge_block_t edge_block_t;
  typedef typename stinger_graph::edge_t edge_t;

public:
  stinger_edge_adapter() : eblock(0), e(0) {}

  stinger_edge_adapter(edge_block_t* eb, edge_t* _e) : eblock(eb), e(_e) {}

  edge_block_t* eblock;
  edge_t* e;
};

template <typename stinger_graph>
class stinger_thread_vertex_iterator {
private:
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::l_vertex_t l_vertex_t;

public:
  stinger_thread_vertex_iterator() :
    id((std::numeric_limits<size_type>::max)()) {}

  stinger_thread_vertex_iterator(l_vertex_t v) : id(v) {}

  stinger_thread_vertex_iterator(const stinger_thread_vertex_iterator& rhs) :
    id(rhs.id) {}

  stinger_thread_vertex_iterator&
  operator=(const stinger_thread_vertex_iterator& rhs)
  {
    id = rhs.id;
    return *this;
  }

  stinger_thread_vertex_iterator& operator++()
  {
    ++id;
    return *this;
  }

  stinger_thread_vertex_iterator& operator++(int)
  {
    stinger_thread_vertex_iterator temp(*this);

    ++id;
    return temp;
  }

  l_vertex_t operator*() const { return id; }

  bool operator==(const stinger_thread_vertex_iterator& rhs) const
  { return id == rhs.id; }
  bool operator!=(const stinger_thread_vertex_iterator& rhs) const
  { return id != rhs.id; }

  l_vertex_t id;
};

template <typename stinger_graph>
class stinger_thread_edge_iterator {
private:
  typedef typename stinger_graph::l_vertex_t l_vertex_t;
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::vertex_t vertex_t;
  typedef typename stinger_graph::edge_block_t edge_block_t;

public:
  stinger_thread_edge_iterator(stinger_graph& _g, size_type pos) :
    g(_g), vertex_entry(0), edge_block(0), edge_entry(0)
  {
    // Find the vertex containing the pos'th edge.
    size_type n_edges = 0;
    while (vertex_entry < g.num_vertices() &&
           n_edges + g.get_vertex(vertex_entry).out_deg <= pos)
    {
      n_edges += g.get_vertex(vertex_entry).out_deg;
      ++vertex_entry;
    }

    if (vertex_entry < g.num_vertices())
    {
      // Find the edge block in the vertex's adjacency list that contains
      // the pos'th edge.
      edge_block = g.get_vertex(vertex_entry).edge_blocks;
      while (edge_block && (n_edges + edge_block->num_edges) < pos)
      {
        n_edges += edge_block->num_edges;
        edge_block = edge_block->next;
      }

      // Set edge_entry to point to the pos'th edge.
      edge_entry = pos - n_edges;
    }

#ifdef DEBUG
    std::cout << "thread_edge_iterator(" << pos << "): " << vertex_entry
              << "    " << edge_block << "    " << edge_entry << std::endl;
#endif
  }

  stinger_thread_edge_iterator(const stinger_thread_edge_iterator& rhs) :
      g(rhs.g), vertex_entry(rhs.vertex_entry), edge_block(rhs.edge_block),
      edge_entry(rhs.edge_entry) {}

  stinger_thread_edge_iterator&
  operator=(const stinger_thread_edge_iterator& rhs)
  {
    g = rhs.g;
    vertex_entry = rhs.vertex_entry;
    edge_block = rhs.edge_block;
    edge_entry = rhs.edge_entry;
    return *this;
  }

  stinger_thread_edge_iterator& operator++()
  {
    if (edge_entry < edge_block->num_edges - 1)
    {
      // The current edge block still has edges avaiable.  Just go to its
      // next edge.
      ++edge_entry;
    }
    else if (edge_block->next)
    {
      // The current vertex still has edge blocks avaiable.  Just go to its
      // next edge block.
      edge_block = edge_block->next;
      edge_entry = 0;
    }
    else
    {
      // Move to the next vertex that has out edges.
      ++vertex_entry;
      while (vertex_entry < g.num_vertices() &&
             g.get_vertex(vertex_entry).out_deg == 0)
      {
        ++vertex_entry;
      }

      // If we hit the end of the vertex list, set ending iterator values,
      // otherwise point to the first edge block of the current vertex.
      edge_entry = 0;
      edge_block = vertex_entry == g.num_vertices() ? 0 :
                   edge_block = g.get_vertex(vertex_entry).edge_blocks;
    }

#ifdef DEBUG
    std::cout << "          operator++(): " << vertex_entry << "    "
              << edge_block << "    " << edge_entry << std::endl;
#endif

    return *this;
  }

  stinger_thread_edge_iterator& operator++(int)
  {
    stinger_thread_edge_iterator temp(*this);

    if (edge_entry < edge_block->num_edges - 1)
    {
      // The current edge block still has edges avaiable.  Just go to its
      // next edge.
      ++edge_entry;
    }
    else if (edge_block->next)
    {
      // The current vertex still has edge blocks avaiable.  Just go to its
      // next edge block.
      edge_block = edge_block->next;
      edge_entry = 0;
    }
    else
    {
      // Move to the next vertex that has out edges.
      ++vertex_entry;
      while (vertex_entry < g.num_vertices() &&
             g.get_vertex(vertex_entry).out_deg == 0)
      {
        ++vertex_entry;
      }

      // If we hit the end of the vertex list, set ending iterator values,
      // otherwise point to the first edge block of the current vertex.
      edge_entry = 0;
      edge_block = vertex_entry == g.num_vertices() ? 0 :
                   edge_block = g.get_vertex(vertex_entry).edge_blocks;
    }

#ifdef DEBUG
    std::cout << "       operator++(int): " << vertex_entry << "    "
              << edge_block << "    " << edge_entry << std::endl;
#endif

    return temp;
  }

  stinger_edge_adapter<stinger_graph> operator*() const
  {
#ifdef DEBUG
    std::cout << "           operator*(): " << vertex_entry << "    "
              << edge_block << "    " << edge_entry << std::endl;
#endif

    return stinger_edge_adapter<stinger_graph>(edge_block,
                                               edge_block->edges + edge_entry);
  }

  bool operator==(const stinger_thread_edge_iterator& rhs)
  {
    return vertex_entry == rhs.vertex_entry && edge_block == rhs.edge_block &&
           edge_entry == rhs.edge_entry;
  }

  bool operator!=(const stinger_thread_edge_iterator& rhs)
  {
    return vertex_entry != rhs.vertex_entry || edge_block != rhs.edge_block ||
           edge_entry != rhs.edge_entry;
  }

private:
  stinger_graph& g;
  size_type vertex_entry;
  edge_block_t* edge_block;
  size_type edge_entry;
};

template <typename stinger_graph>
class stinger_thread_adjacency_iterator {
private:
  typedef typename stinger_graph::l_vertex_t l_vertex_t;
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::vertex_t vertex_t;
  typedef typename stinger_graph::edge_block_t edge_block_t;

public:
  stinger_thread_adjacency_iterator() : edge_block(0), edge_entry(0) {}

  stinger_thread_adjacency_iterator(stinger_graph& g,
                                    l_vertex_t src, size_type pos)
  {
    vertex_t& v = g.get_vertex(src);

    if (pos == v.out_deg)
    {
      edge_block = 0;
      edge_entry = 0;
    }
    else
    {
      assert(pos < v.out_deg);

      // Find the edge block and edge list entry for the pos'th adjacency.
      edge_block_t* cur_edge_block = v.edge_blocks;

      size_type num_edges_so_far = cur_edge_block->num_edges;
      while (num_edges_so_far <= pos)
      {
        cur_edge_block = cur_edge_block->next;
        num_edges_so_far += cur_edge_block->num_edges;
      }

      edge_block = cur_edge_block;
      edge_entry = pos - (num_edges_so_far - cur_edge_block->num_edges);
    }
#ifdef DEBUG
    std::cout << "thread_adjacency_iterator(" << pos << "): " << edge_block
              << "    " << edge_entry << std::endl;
#endif
  }

  stinger_thread_adjacency_iterator(
    const stinger_thread_adjacency_iterator& rhs) :
      edge_block(rhs.edge_block), edge_entry(rhs.edge_entry) {}

  stinger_thread_adjacency_iterator&
  operator=(const stinger_thread_adjacency_iterator& rhs)
  {
    edge_block = rhs.edge_block;
    edge_entry = rhs.edge_entry;
    return *this;
  }

  stinger_thread_adjacency_iterator& operator++()
  {
    if (edge_entry < edge_block->num_edges - 1)
    {
      ++edge_entry;
    }
    else
    {
      edge_block = edge_block->next;
      edge_entry = 0;
    }
#ifdef DEBUG
    std::cout << "               operator++(): " << edge_block << "    "
              << edge_entry << std::endl;
#endif

    return *this;
  }

  stinger_thread_adjacency_iterator& operator++(int)
  {
    stinger_thread_adjacency_iterator temp(*this);

    if (edge_entry < edge_block->num_edges - 1)
    {
      ++edge_entry;
    }
    else
    {
      edge_block = edge_block->next;
      edge_entry = 0;
    }
#ifdef DEBUG
    std::cout << "            operator++(int): " << edge_block << "    "
              << edge_entry << std::endl;
#endif

    return temp;
  }

  l_vertex_t operator*() const
  {
#ifdef DEBUG
    std::cout << "                operator*(): " << edge_block << "    "
              << edge_entry << std::endl;
#endif

    return edge_block->edges[edge_entry].to_id;
  }

  bool operator==(const stinger_thread_adjacency_iterator& rhs)
  { return edge_block == rhs.edge_block && edge_entry == rhs.edge_entry; }
  bool operator!=(const stinger_thread_adjacency_iterator& rhs)
  { return edge_block != rhs.edge_block || edge_entry != rhs.edge_entry; }

private:
  edge_block_t* edge_block;
  size_type edge_entry;
};

template <typename stinger_graph>
class stinger_thread_out_edge_iterator {
private:
  typedef typename stinger_graph::l_vertex_t l_vertex_t;
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::vertex_t vertex_t;
  typedef typename stinger_graph::edge_block_t edge_block_t;

public:
  stinger_thread_out_edge_iterator() : edge_block(0), edge_entry(0) {}

  stinger_thread_out_edge_iterator(stinger_graph& g,
                                   l_vertex_t src, size_type pos)
  {
    vertex_t& v = g.get_vertex(src);

    if (pos == v.out_deg)
    {
      edge_block = 0;
      edge_entry = 0;
    }
    else
    {
      assert(pos < v.out_deg);

      // Find the edge block and edge list entry for the pos'th adjacency.
      edge_block_t* cur_edge_block = v.edge_blocks;
      size_type num_edges_so_far = cur_edge_block->num_edges;
      while (num_edges_so_far <= pos)
      {
        cur_edge_block = cur_edge_block->next;
        num_edges_so_far += cur_edge_block->num_edges;
      }

      edge_block = cur_edge_block;
      edge_entry = pos - (num_edges_so_far - cur_edge_block->num_edges);
    }
#ifdef DEBUG
    std::cout << "thread_out_edge_iterator(" << pos << "): " << edge_block
              << "    " << edge_entry << std::endl;
#endif
  }

  stinger_thread_out_edge_iterator(
    const stinger_thread_out_edge_iterator& rhs) :
      edge_block(rhs.edge_block), edge_entry(rhs.edge_entry) {}

  stinger_thread_out_edge_iterator&
  operator=(const stinger_thread_out_edge_iterator& rhs)
  {
    edge_block = rhs.edge_block;
    edge_entry = rhs.edge_entry;
    return *this;
  }

  stinger_thread_out_edge_iterator& operator++()
  {
    if (edge_entry < edge_block->num_edges - 1)
    {
      ++edge_entry;
    }
    else
    {
      edge_block = edge_block->next;
      edge_entry = 0;
    }
#ifdef DEBUG
    std::cout << "               operator++(): " << edge_block << "    "
              << edge_entry << std::endl;
#endif

    return *this;
  }

  stinger_thread_out_edge_iterator& operator++(int)
  {
    stinger_thread_out_edge_iterator temp(*this);

    if (edge_entry < edge_block->num_edges - 1)
    {
      ++edge_entry;
    }
    else
    {
      edge_block = edge_block->next;
      edge_entry = 0;
    }
#ifdef DEBUG
    std::cout << "            operator++(int): " << edge_block << "    "
              << edge_entry << std::endl;
#endif

    return temp;
  }

  stinger_edge_adapter<stinger_graph> operator*() const
  {
#ifdef DEBUG
    std::cout << "                operator*(): " << edge_block << "    "
              << edge_entry << std::endl;
#endif

    return stinger_edge_adapter<stinger_graph>(edge_block,
                                               edge_block->edges + edge_entry);
  }

  bool operator==(const stinger_thread_out_edge_iterator& rhs)
  { return edge_block == rhs.edge_block && edge_entry == rhs.edge_entry; }
  bool operator!=(const stinger_thread_out_edge_iterator& rhs)
  { return edge_block != rhs.edge_block || edge_entry != rhs.edge_entry; }

private:
  edge_block_t* edge_block;
  size_type edge_entry;
};

}

template <typename stinger_graph>
class stinger_graph_adapter {
public:
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::l_vertex_t vertex_descriptor;
  typedef detail::stinger_edge_adapter<stinger_graph> edge_descriptor;
  typedef void vertex_iterator;
  typedef void adjacency_iterator;
  typedef void in_adjacency_iterator;
  typedef void edge_iterator;
  typedef void out_edge_iterator;
  typedef void in_edge_iterator;
  typedef detail::stinger_thread_vertex_iterator<stinger_graph>
          thread_vertex_iterator;
  typedef detail::stinger_thread_adjacency_iterator<stinger_graph>
          thread_adjacency_iterator;
  typedef void thread_in_adjacency_iterator;
  typedef detail::stinger_thread_edge_iterator<stinger_graph>
          thread_edge_iterator;
  typedef detail::stinger_thread_out_edge_iterator<stinger_graph>
          thread_out_edge_iterator;
  typedef void thread_in_edge_iterator;
  typedef directedS directed_category;
  typedef thread_iterators iterator_category;

  typedef typename stinger_graph::vertex_t vertex_t;
  typedef typename stinger_graph::edge_t edge_t;
  typedef typename stinger_graph::l_vertex_t l_vertex_t;
  typedef typename stinger_graph::edge_block_t edge_block_t;
  typedef typename stinger_graph::vtype_t vtype_t;
  typedef typename stinger_graph::vweight_t vweight_t;
  typedef typename stinger_graph::etype_t etype_t;
  typedef typename stinger_graph::eweight_t eweight_t;

  stinger_graph_adapter(stinger_graph& g) : graph(g) {}
  ~stinger_graph_adapter() {}

  stinger_graph_adapter& operator=(const stinger_graph_adapter& rhs)
  {
    clear();
    graph = rhs.graph;
    return *this;
  }

  void clear() { graph.clear(); }

  void init(size_type order, size_type size, size_type* srcs, size_type* dests)
  {
    clear();

    // Declare temporary memory for the types and weights, and initialize
    // them all to be 0.
    vtype_t* vtypes = (vtype_t*) calloc(order, sizeof(vtype_t));
    vweight_t* vweights = (vweight_t*) calloc(order, sizeof(vweight_t));
    etype_t* etypes = (etype_t*) calloc(size, sizeof(etype_t));
    eweight_t* eweights = (eweight_t*) calloc(size, sizeof(eweight_t));

    // Initialize the stinger graph.
    graph.init_vertices(order, vtypes, vweights);
    graph.init_edges(size, srcs, dests, etypes, eweights);

    free(vtypes);
    free(vweights);
    free(etypes);
    free(eweights);
  }

  size_type get_order() const { return graph.num_vertices(); }
  size_type get_size() const { return graph.num_edges(); }

  size_type get_degree(const vertex_descriptor& i) const
  {
    return get_out_degree(i);
  }

  size_type get_out_degree(const vertex_descriptor& i) const
  {
    return graph.out_degree(i);
  }

  thread_vertex_iterator
  thread_vertices(size_type pos) const
  {
    return thread_vertex_iterator(static_cast<l_vertex_t>(pos));
  }

  thread_adjacency_iterator
  thread_adjacent_vertices(const vertex_descriptor& v, size_type pos) const
  {
    return thread_adjacency_iterator(graph, v, pos);
  }

  thread_edge_iterator
  thread_edges(size_type pos) const
  {
    return thread_edge_iterator(graph, pos);
  }

  thread_out_edge_iterator
  thread_out_edges(const vertex_descriptor& v, size_type pos) const
  {
    return thread_out_edge_iterator(graph, v, pos);
  }

  void print() const
  {
    size_type num_verts = get_order();
    size_type num_edges = get_size();

    std::cout << "num_verts: " << num_verts << std::endl;
    std::cout << "num_edges: " << num_edges << std::endl;

    for (l_vertex_t u = 0; u < num_verts; ++u)
    {
      vertex_t& u_data = graph.get_vertex(u);

      std::cout << u << " : { ";

      for (edge_block_t* b = u_data.edge_blocks; b != 0; b = b->next)
      {
        for (size_type i = 0; i < b->num_edges; ++i)
        {
          l_vertex_t v = b->edges[i].to_id;
          std::cout << "(" << u << ", " << v << ") ";
        }
      }

      std::cout << "}" << std::endl;
    }
  }

  unsigned long get_mmap_size()
  {
    return graph.get_mmap_size();
  }

  void write_mmap(void* mapped_mem)
  {
    return graph.write_mmap(mapped_mem,
                            mmap_traits<stinger_graph_adapter>::type);
  }

  void read_mmap(void* mapped_mem)
  {
    return graph.read_mmap(mapped_mem);
  }

private:
  stinger_graph& graph;
};

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::size_type
num_vertices(const stinger_graph_adapter<stinger_graph>& g)
{
  return g.get_order();
}


template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::size_type
num_edges(const stinger_graph_adapter<stinger_graph>& g)
{
  return g.get_size();
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::vertex_descriptor
source(
    const typename stinger_graph_adapter<stinger_graph>::edge_descriptor& e,
    const stinger_graph_adapter<stinger_graph>& g)
{
  return e.eblock->from_id;
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::vertex_descriptor
target(
    const typename stinger_graph_adapter<stinger_graph>::edge_descriptor& e,
    const stinger_graph_adapter<stinger_graph>& g)
{
  return e.e->to_id;
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::size_type
degree(
  const typename stinger_graph_adapter<stinger_graph>::vertex_descriptor& v,
  const stinger_graph_adapter<stinger_graph>& g)
{
  return g.get_degree(v);
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::size_type
out_degree(
  const typename stinger_graph_adapter<stinger_graph>::vertex_descriptor& v,
  const stinger_graph_adapter<stinger_graph>& g)
{
  return g.get_out_degree(v);
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::thread_vertex_iterator
thread_vertices(
  const typename stinger_graph_adapter<stinger_graph>::size_type pos,
  const stinger_graph_adapter<stinger_graph>& g)
{
  return g.thread_vertices(pos);
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::thread_adjacency_iterator
thread_adjacent_vertices(
  const typename stinger_graph_adapter<stinger_graph>::vertex_descriptor& v,
  typename stinger_graph_adapter<stinger_graph>::size_type pos,
  const stinger_graph_adapter<stinger_graph>& g)
{
  return g.thread_adjacent_vertices(v, pos);
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::thread_edge_iterator
thread_edges(
  const typename stinger_graph_adapter<stinger_graph>::size_type pos,
  const stinger_graph_adapter<stinger_graph>& g)
{
  return g.thread_edges(pos);
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::thread_out_edge_iterator
thread_out_edges(
  const typename stinger_graph_adapter<stinger_graph>::vertex_descriptor& v,
  typename stinger_graph_adapter<stinger_graph>::size_type pos,
  const stinger_graph_adapter<stinger_graph>& g)
{
  return g.thread_out_edges(v, pos);
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::vertex_descriptor
null_vertex(const stinger_graph_adapter<stinger_graph>& g)
{
  return (std::numeric_limits<typename stinger_graph_adapter<stinger_graph>::
                                       vertex_descriptor>::max)();
}

template <typename stinger_graph>
inline
typename stinger_graph_adapter<stinger_graph>::edge_descriptor
null_edge(const stinger_graph_adapter<stinger_graph>& g)
{
  return typename stinger_graph_adapter<stinger_graph>::edge_descriptor();
}

template <typename ITERATOR, typename stinger_graph>
inline
bool
is_valid(ITERATOR& iter,
         typename stinger_graph_adapter<stinger_graph>::size_type p,
         const stinger_graph_adapter<stinger_graph>& tg)
{
  return true;
}

template <typename stinger_graph>
inline
bool
is_directed(const stinger_graph_adapter<stinger_graph>& g)
{
  return true;
}

template <typename stinger_graph>
inline
bool
is_undirected(const stinger_graph_adapter<stinger_graph>& g)
{
  return false;
}

template <typename stinger_graph>
inline
bool
is_bidirectional(const stinger_graph_adapter<stinger_graph>& g)
{
  return false;
}

template <typename stinger_graph>
inline
void
print(stinger_graph_adapter<stinger_graph>& g)
{
  typedef stinger_graph_adapter<stinger_graph> Graph;
  typedef typename Graph::size_type size_type;
  typedef typename Graph::vertex_descriptor vertex_descriptor;
  typedef typename Graph::thread_vertex_iterator thread_vertex_iterator;
  typedef typename Graph::thread_adjacency_iterator thread_adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  thread_vertex_iterator verts = thread_vertices(0, g);
  thread_vertex_iterator vertsEnd = thread_vertices(num_vertices(g), g);

  for ( ; verts != vertsEnd; ++verts)
  {
    vertex_descriptor u = *verts;

    size_type uid = get(vid_map, u);

    std::cout << uid << " : { ";
    std::cout.flush();

    thread_adjacency_iterator adjs = thread_adjacent_vertices(u, 0, g);
    thread_adjacency_iterator adjsEnd =
      thread_adjacent_vertices(u, out_degree(u, g), g);

    for ( ; adjs != adjsEnd; ++adjs)
    {
      vertex_descriptor v = *adjs;
      size_type vid = get(vid_map, v);

      std::cout << "(" << uid << ", " << vid << ") ";
    }

    std::cout << "}" << std::endl;
  }
}

template <typename stinger_graph>
class vertex_id_map<stinger_graph_adapter<stinger_graph> > :
  public put_get_helper<
             typename stinger_graph_adapter<stinger_graph>::size_type,
             vertex_id_map<stinger_graph_adapter<stinger_graph> > > {
public:
  typedef typename stinger_graph_adapter<stinger_graph>::vertex_descriptor
          key_type;
  typedef typename stinger_graph_adapter<stinger_graph>::size_type value_type;

  vertex_id_map() {}
  value_type operator[] (const key_type& k) const { return k; }
};

template <typename stinger_graph>
class vertex_id_map<const stinger_graph_adapter<stinger_graph> > :
  public put_get_helper<
             typename stinger_graph_adapter<stinger_graph>::size_type,
             vertex_id_map<const stinger_graph_adapter<stinger_graph> > > {
public:
  typedef typename stinger_graph_adapter<stinger_graph>::vertex_descriptor
          key_type;
  typedef typename stinger_graph_adapter<stinger_graph>::size_type value_type;

  vertex_id_map() {}
  value_type operator[] (const key_type& k) const { return k; }
};

template <typename TYPE, typename WEIGHT, typename TS>
class mmap_traits<stinger_graph_adapter<stinger_graph<TYPE, WEIGHT, TS> > > {
public:
  // If TYPE, WEIGHT, or TS is an undefined type, then the graph adapter's
  // type should also be undefined.
  static const unsigned long type =
    (mmap_traits<TYPE>::type != MMAP_TYPE_NOT_DEFINED) *
    (mmap_traits<WEIGHT>::type != MMAP_TYPE_NOT_DEFINED) *
    (mmap_traits<TS>::type != MMAP_TYPE_NOT_DEFINED) *
    (10000 + mmap_traits<TYPE>::type * 256 + mmap_traits<WEIGHT>::type * 16 +
     mmap_traits<TS>::type);
};

template <typename stinger_graph>
inline
void init(typename stinger_graph_adapter<stinger_graph>::size_type n,
          typename stinger_graph_adapter<stinger_graph>::size_type m,
          typename stinger_graph_adapter<stinger_graph>::size_type* srcs,
          typename stinger_graph_adapter<stinger_graph>::size_type* dests,
          stinger_graph_adapter<stinger_graph>& g)
{
  g.init(n, m, srcs, dests);
}

template <typename stinger_graph>
inline
void clear(stinger_graph_adapter<stinger_graph>& g)
{
  return g.clear();
}

}

#endif
