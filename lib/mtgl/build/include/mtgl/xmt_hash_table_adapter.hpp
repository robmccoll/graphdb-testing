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
/*! \file xmt_hash_table_adapter.hpp

    \author Jon Berry (jberry@sandia.gov)

    \date 4/5/2008
*/
/****************************************************************************/

#ifndef MTGL_XMT_HASH_TABLE_ADAPTER_HPP
#define MTGL_XMT_HASH_TABLE_ADAPTER_HPP

#include <cstdlib>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/util.hpp>
#include <mtgl/xmt_hash_table.hpp>

namespace mtgl {

namespace detail {

class xmt_hash_table_edge_adapter : public pair<uint64_t, uint64_t> {
public:
  xmt_hash_table_edge_adapter() :
    pair<uint64_t, uint64_t>((std::numeric_limits<uint64_t>::max)(),
                             (std::numeric_limits<uint64_t>::max)()),
    id((std::numeric_limits<uint64_t>::max)()) {}

  xmt_hash_table_edge_adapter(uint64_t v1, uint64_t v2, uint64_t eid) :
    pair<uint64_t, uint64_t>(v1, v2), id(eid) {}
  xmt_hash_table_edge_adapter(pair<uint64_t, uint64_t> p, uint64_t eid) :
    pair<uint64_t, uint64_t>(p), id(eid) {}

  uint64_t get_id() const { return id; }

private:
  uint64_t id;
};

class xmt_hash_table_vertex_iterator {
public:
  xmt_hash_table_vertex_iterator(uint64_t* verts = 0) : vertices(verts) {}

  uint64_t operator[](uint64_t p) const { return vertices[p]; }

private:
  uint64_t* vertices;
};

class xmt_hash_table_edge_iterator {
public:
  xmt_hash_table_edge_iterator() : the_edges(0) {}
  xmt_hash_table_edge_iterator(pair<uint64_t, uint64_t>* edgs) :
    the_edges(edgs) {}

  xmt_hash_table_edge_adapter operator[](uint64_t p) const
  {
    return xmt_hash_table_edge_adapter(the_edges[p], p);
  }

private:
  pair<uint64_t, uint64_t>* the_edges;
};

template <typename HT>
class xmt_hash_table_adapter_init_visitor {
public:
  xmt_hash_table_adapter_init_visitor(pair<uint64_t, uint64_t>* te) : index(0), the_edges(te) {}

  void operator()(typename HT::value_type* i)
  {
    uint64_t one = 1;
    uint64_t ind = mt_incr(index, one);
    the_edges[ind] = i->second;
  }

  uint64_t index;
  pair<uint64_t, uint64_t>* the_edges;
};

}

class xmt_hash_table_adapter {
public:
  typedef uint64_t size_type;
  typedef uint64_t vertex_descriptor;
  typedef detail::xmt_hash_table_edge_adapter edge_descriptor;
  typedef detail::xmt_hash_table_vertex_iterator vertex_iterator;
  typedef void adjacency_iterator;            // none!
  typedef void in_adjacency_iterator;         // none!
  typedef detail::xmt_hash_table_edge_iterator edge_iterator;
  typedef void out_edge_iterator;             // none!
  typedef void in_edge_iterator;              // none!
  typedef void thread_vertex_iterator;        // none!
  typedef void thread_adjacency_iterator;     // none!
  typedef void thread_in_adjacency_iterator;  // none!
  typedef void thread_edge_iterator;          // none!
  typedef void thread_out_edge_iterator;      // none!
  typedef void thread_in_edge_iterator;       // none!
  typedef directedS directed_category;
  typedef vector_iterators iterator_category;

  xmt_hash_table_adapter(xmt_hash_table<uint64_t,
                                        pair<uint64_t, uint64_t> >* edgetable,
                         uint64_t* vertices, uint64_t nn) :
    the_table(edgetable), n(nn), m(edgetable->size()), the_vertices(vertices)
  {
    the_edges = new pair<uint64_t, uint64_t>[m];

    detail::xmt_hash_table_adapter_init_visitor<
      xmt_hash_table<uint64_t, pair<uint64_t, uint64_t> > >
        ivis(the_edges);

    the_table->visit(ivis);
  }

  ~xmt_hash_table_adapter() { delete [] the_edges; }

  size_type get_order() const { return n; }
  size_type get_size() const { return m; }

  vertex_iterator vertices() const { return vertex_iterator(the_vertices); }
  edge_iterator edges() const { return edge_iterator(the_edges); }

private:
  xmt_hash_table<uint64_t, pair<uint64_t, uint64_t> >* the_table;
  size_type n;
  size_type m;
  uint64_t* the_vertices;
  pair<uint64_t, uint64_t>* the_edges;
};

inline
uint64_t
num_vertices(const xmt_hash_table_adapter& g)
{
  return g.get_order();
}

inline
uint64_t
num_edges(const xmt_hash_table_adapter& g)
{
  return g.get_size();
}

inline
xmt_hash_table_adapter::vertex_iterator
vertices(const xmt_hash_table_adapter& g)
{
  return g.vertices();
}

inline
xmt_hash_table_adapter::edge_iterator
edges(const xmt_hash_table_adapter& g)
{
  return g.edges();
}

inline
uint64_t
source(const xmt_hash_table_adapter::edge_descriptor& e,
       const xmt_hash_table_adapter& g)
{
  return e.first;
}

inline
uint64_t
target(const xmt_hash_table_adapter::edge_descriptor& e,
       const xmt_hash_table_adapter& g)
{
  return e.second;
}

inline
xmt_hash_table_adapter::vertex_descriptor
null_vertex(const xmt_hash_table_adapter& g)
{
  return (std::numeric_limits<xmt_hash_table_adapter::vertex_descriptor>::max)();
}

inline
xmt_hash_table_adapter::edge_descriptor
null_edge(const xmt_hash_table_adapter& g)
{
  return xmt_hash_table_adapter::edge_descriptor();
}

template <typename ITERATOR>
inline
bool
is_valid(ITERATOR& iter, xmt_hash_table_adapter::size_type p,
         const xmt_hash_table_adapter& g)
{
  return true;
}

inline
bool
is_directed(const xmt_hash_table_adapter& g)
{
  return true;
}

inline
bool
is_undirected(const xmt_hash_table_adapter& g)
{
  return false;
}

inline
bool
is_bidirectional(const xmt_hash_table_adapter& g)
{
  return false;
}

template <>
class vertex_id_map<xmt_hash_table_adapter> :
  public put_get_helper<uint64_t, vertex_id_map<xmt_hash_table_adapter> > {
public:
  typedef graph_traits<xmt_hash_table_adapter>::vertex_descriptor key_type;
  typedef uint64_t value_type;

  vertex_id_map() {}
  value_type operator[] (const key_type& k) const { return k; }
};

template <>
class edge_id_map<xmt_hash_table_adapter> :
  public put_get_helper<uint64_t, edge_id_map<xmt_hash_table_adapter> > {
public:
  typedef graph_traits<xmt_hash_table_adapter>::edge_descriptor key_type;
  typedef uint64_t value_type;

  edge_id_map() {}
  value_type operator[] (const key_type& k) const { return k.get_id(); }
};

}

#endif
