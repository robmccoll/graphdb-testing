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
/*! \file edge_array_adapter.hpp

    \author Jon Berry (jberry@sandia.gov)

    \date 12/11/2008
*/
/****************************************************************************/

#ifndef MTGL_EDGE_ARRAY_ADAPTER_HPP
#define MTGL_EDGE_ARRAY_ADAPTER_HPP

#include <cstdlib>
#include <limits>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/util.hpp>

namespace mtgl {

namespace detail {

template <typename uinttype>
class edge_array_edge_adapter {
public:
  edge_array_edge_adapter() : first((std::numeric_limits<uinttype>::max)()),
                              second((std::numeric_limits<uinttype>::max)()),
                              id((std::numeric_limits<uinttype>::max)()) {}

  template <typename T>
  edge_array_edge_adapter(uinttype v1, uinttype v2, T eid) :
    first(v1), second(v2), id((T)eid) {}

public:
  uinttype first;
  uinttype second;
  uinttype id;
};

template <typename uinttype>
class edge_array_vertex_iterator {
public:
  // We assume that vertices are packed 0..n-1.
  edge_array_vertex_iterator() {}

  uinttype operator[](uinttype p) const { return p; }
};

template <typename uinttype>
class edge_array_edge_iterator {
public:
  edge_array_edge_iterator() : srcs(0), dests(0) {}  // Doesn't work!
  edge_array_edge_iterator(uinttype* scs, uinttype* dsts) :
    srcs(scs), dests(dsts) {}

  edge_array_edge_adapter<uinttype> operator[](uinttype p) const
  {
    return edge_array_edge_adapter<uinttype>(srcs[p], dests[p], p);
  }

private:
  uinttype* srcs;
  uinttype* dests;
};

template <typename uinttype>
class edge_array_thread_vertex_iterator {
public:
  edge_array_thread_vertex_iterator() :
    id((std::numeric_limits<uinttype>::max)()) {}

  edge_array_thread_vertex_iterator(uinttype v) : id(v) {}

  edge_array_thread_vertex_iterator(
      const edge_array_thread_vertex_iterator& rhs) : id(rhs.id) {}

  edge_array_thread_vertex_iterator&
  operator=(const edge_array_thread_vertex_iterator& rhs)
  {
    id = rhs.id;
    return *this;
  }

  edge_array_thread_vertex_iterator& operator++()
  {
    ++id;
    return *this;
  }

  edge_array_thread_vertex_iterator& operator++(int)
  {
    edge_array_thread_vertex_iterator temp(*this);

    ++id;
    return temp;
  }

  uinttype operator*() const { return id; }

  bool operator==(const edge_array_thread_vertex_iterator& rhs) const
  { return id == rhs.id; }
  bool operator!=(const edge_array_thread_vertex_iterator& rhs) const
  { return id != rhs.id; }

private:
  uinttype id;
};

template <typename uinttype>
class edge_array_thread_edge_iterator {
public:
  edge_array_thread_edge_iterator() :
    srcs(0), dests(0), id((std::numeric_limits<uinttype>::max)()) {}

  edge_array_thread_edge_iterator(uinttype* scs, uinttype* dsts, uinttype v) :
    id(v), srcs(scs), dests(dsts) {}

  edge_array_thread_edge_iterator(const edge_array_thread_edge_iterator& rhs) :
    id(rhs.id), srcs(rhs.srcs), dests(rhs.dests) {}

  edge_array_thread_edge_iterator&
  operator=(const edge_array_thread_edge_iterator& rhs)
  {
    id = rhs.id;
    srcs = rhs.srcs;
    dests = rhs.dests;
    return *this;
  }

  edge_array_thread_edge_iterator& operator++()
  {
    ++id;
    return *this;
  }

  edge_array_thread_edge_iterator& operator++(int)
  {
    edge_array_thread_edge_iterator temp(*this);

    ++id;
    return temp;
  }

  edge_array_edge_adapter<uinttype> operator*() const
  {
    return edge_array_edge_adapter<uinttype>(srcs[id], dests[id], id);
  }

  bool operator==(const edge_array_thread_edge_iterator& rhs) const
  { return id == rhs.id; }
  bool operator!=(const edge_array_thread_edge_iterator& rhs) const
  { return id != rhs.id; }

private:
  uinttype id;
  uinttype* srcs;
  uinttype* dests;
};

}

template<typename uinttype>
class edge_array_adapter {
public:
  typedef uinttype size_type;
  typedef uinttype vertex_descriptor;
  typedef detail::edge_array_edge_adapter<uinttype> edge_descriptor;
  typedef detail::edge_array_vertex_iterator<uinttype> vertex_iterator;
  typedef void adjacency_iterator;            // none!
  typedef void in_adjacency_iterator;         // none!
  typedef detail::edge_array_edge_iterator<uinttype> edge_iterator;
  typedef void out_edge_iterator;             // none!
  typedef void in_edge_iterator;              // none!
  typedef detail::edge_array_thread_vertex_iterator<uinttype>
          thread_vertex_iterator;
  typedef void thread_adjacency_iterator;     // none!
  typedef void thread_in_adjacency_iterator;  // none!
  typedef detail::edge_array_thread_edge_iterator<uinttype>
          thread_edge_iterator;
  typedef void thread_out_edge_iterator;      // none!
  typedef void thread_in_edge_iterator;       // none!
  typedef directedS directed_category;
  typedef vector_iterators iterator_category;

  // We assume that vertices are packed 0..n-1.
  edge_array_adapter(size_type* scs, size_type* dsts, size_type nn,
                     size_type mm) : srcs(scs), dests(dsts),  n(nn), m(mm) {}

  size_type get_order() const { return n; }
  size_type get_size() const { return m; }

  vertex_iterator vertices() const { return vertex_iterator(); }
  edge_iterator edges() const { return edge_iterator(srcs, dests); }

  thread_vertex_iterator thread_vertices(size_type pos) const
  {
    return thread_vertex_iterator(pos);
  }

  thread_edge_iterator thread_edges(size_type pos) const
  {
    return thread_edge_iterator(srcs, dests, pos);
  }

private:
  size_type* srcs;
  size_type* dests;
  size_type n;
  size_type m;
};

template <typename uinttype>
inline
uinttype
num_vertices(const edge_array_adapter<uinttype>& g)
{
  return g.get_order();
}

template <typename uinttype>
inline
uinttype
num_edges(const edge_array_adapter<uinttype>& g)
{
  return g.get_size();
}

template <typename uinttype>
inline
uinttype
source(const typename edge_array_adapter<uinttype>::edge_descriptor& e,
       const edge_array_adapter<uinttype>& g)
{
  return e.first;
}

template <typename uinttype>
inline
uinttype
target(const typename edge_array_adapter<uinttype>::edge_descriptor& e,
       const edge_array_adapter<uinttype>& g)
{
  return e.second;
}

template <typename uinttype>
inline
typename edge_array_adapter<uinttype>::vertex_iterator
vertices(const edge_array_adapter<uinttype>& g)
{
  return g.vertices();
}

template <typename uinttype>
inline
typename edge_array_adapter<uinttype>::edge_iterator
edges(const edge_array_adapter<uinttype>& g)
{
  return g.edges();
}

template <typename uinttype>
inline
typename edge_array_adapter<uinttype>::thread_vertex_iterator
thread_vertices(uinttype pos, const edge_array_adapter<uinttype>& g)
{
  return g.thread_vertices(pos);
}

template <typename uinttype>
inline
typename edge_array_adapter<uinttype>::thread_edge_iterator
thread_edges(uinttype pos, const edge_array_adapter<uinttype>& g)
{
  return g.thread_edges(pos);
}

template <typename uinttype>
inline
typename edge_array_adapter<uinttype>::vertex_descriptor
null_vertex(const edge_array_adapter<uinttype>& g)
{
  return (std::numeric_limits<typename edge_array_adapter<uinttype>::
                             vertex_descriptor>::max)();
}

template <typename uinttype>
inline
typename edge_array_adapter<uinttype>::edge_descriptor
null_edge(const edge_array_adapter<uinttype>& g)
{
  return typename edge_array_adapter<uinttype>::edge_descriptor();
}

template <typename ITERATOR, typename uinttype>
inline
bool
is_valid(ITERATOR& iter, typename edge_array_adapter<uinttype>::size_type p,
         const edge_array_adapter<uinttype>& g)
{
  return true;
}

template <typename uinttype>
inline
bool
is_directed(const edge_array_adapter<uinttype>& g)
{
  return true;
}

template <typename uinttype>
inline
bool
is_undirected(const edge_array_adapter<uinttype>& g)
{
  return false;
}

template <typename uinttype>
inline
bool
is_bidirectional(const edge_array_adapter<uinttype>& g)
{
  return false;
}

template <typename uinttype>
class vertex_id_map<edge_array_adapter<uinttype> > :
  public put_get_helper<uint64_t, vertex_id_map<edge_array_adapter
                                                <uinttype> > > {
public:
  typedef typename
    graph_traits<edge_array_adapter<uinttype> >::vertex_descriptor key_type;
  typedef uinttype value_type;

  vertex_id_map() {}
  value_type operator[] (const key_type& k) const { return k; }
};

template <typename uinttype>
class edge_id_map<edge_array_adapter<uinttype> > :
  public put_get_helper<uinttype, edge_id_map<edge_array_adapter
                                              <uinttype> > > {
public:
  typedef typename
    graph_traits<edge_array_adapter<uinttype> >::edge_descriptor key_type;
  typedef uinttype value_type;

  edge_id_map() {}
  value_type operator[] (const key_type& k) const { return k.id; }
};

}

#endif
