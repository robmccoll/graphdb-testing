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
/*! \file filter_adapter.hpp

    \brief This adapter provides an interface to access a filtered version
           of a graph.  The graph is filtered based on user-supplied vertex
           and edge filters.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 8/17/2010

    See test/test_filtering.cpp for an example of using the filter adapter.
*/
/****************************************************************************/

#ifndef MTGL_FILTER_ADAPTER_HPP
#define MTGL_FILTER_ADAPTER_HPP

#include <mtgl/mtgl_adapter.hpp>

namespace mtgl {

namespace detail {

template <typename vertex_filter, typename iterator, typename graph_adapter>
class filter_vertex_iterator {
public:
  typedef typename graph_traits<graph_adapter>::size_type size_type;
  typedef typename graph_traits<graph_adapter>::vertex_descriptor
          vertex_descriptor;

  filter_vertex_iterator() {}
  filter_vertex_iterator(const vertex_filter& vf, const iterator& i,
                         const graph_adapter& _g) :
    vfilter(vf), base_iter(i), g(_g) {}

  bool is_valid(const vertex_descriptor& v) { return vfilter(v); }

  vertex_descriptor operator[](size_type p)
  {
    vertex_descriptor v = base_iter[p];
    return is_valid(v) ? v : null_vertex(g);
  }

public:
  vertex_filter vfilter;
  iterator base_iter;
  const graph_adapter& g;
};

template <typename vertex_filter, typename edge_filter, typename graph_adapter>
class filter_edge_iterator {
public:
  typedef typename graph_traits<graph_adapter>::size_type size_type;
  typedef typename graph_traits<graph_adapter>::vertex_descriptor
          vertex_descriptor;
  typedef typename graph_traits<graph_adapter>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<graph_adapter>::edge_iterator edge_iterator;

  filter_edge_iterator() {}
  filter_edge_iterator(const vertex_filter& vf, const edge_filter& ef,
                       const edge_iterator& i, const graph_adapter& _g) :
    vfilter(vf), efilter(ef), base_iter(i), g(_g) {}

  bool is_valid(const edge_descriptor& e)
  {
    return vfilter(source(e, g)) && vfilter(target(e, g)) && efilter(e);
  }

  edge_descriptor operator[](size_type p)
  {
    edge_descriptor e = base_iter[p];
    return is_valid(e) ? e : null_edge(g);
  }

public:
  vertex_filter vfilter;
  edge_filter efilter;
  edge_iterator base_iter;
  const graph_adapter& g;
};

template <typename vertex_filter, typename edge_filter, typename graph_adapter>
class filter_in_edge_iterator {
public:
  typedef typename graph_traits<graph_adapter>::size_type size_type;
  typedef typename graph_traits<graph_adapter>::vertex_descriptor
          vertex_descriptor;
  typedef typename graph_traits<graph_adapter>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<graph_adapter>::in_edge_iterator
          in_edge_iterator;

  filter_in_edge_iterator() {}
  filter_in_edge_iterator(const vertex_filter& vf, const edge_filter& ef,
                          const in_edge_iterator& i, const graph_adapter& _g) :
    vfilter(vf), efilter(ef), base_iter(i), g(_g) {}

  bool is_valid(const edge_descriptor& e)
  {
    return vfilter(source(e, g)) && efilter(e);
  }

  edge_descriptor operator[](size_type p)
  {
    edge_descriptor e = base_iter[p];
    return is_valid(e) ? e : null_edge(g);
  }

public:
  vertex_filter vfilter;
  edge_filter efilter;
  in_edge_iterator base_iter;
  const graph_adapter& g;
};

template <typename vertex_filter, typename edge_filter, typename graph_adapter>
class filter_out_edge_iterator {
public:
  typedef typename graph_traits<graph_adapter>::size_type size_type;
  typedef typename graph_traits<graph_adapter>::vertex_descriptor
          vertex_descriptor;
  typedef typename graph_traits<graph_adapter>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<graph_adapter>::out_edge_iterator
          out_edge_iterator;

  filter_out_edge_iterator() {}
  filter_out_edge_iterator(const vertex_filter& vf, const edge_filter& ef,
                           const out_edge_iterator& i,
                           const graph_adapter& _g) :
    vfilter(vf), efilter(ef), base_iter(i), g(_g) {}

  bool is_valid(const edge_descriptor& e)
  {
    return efilter(e) && vfilter(target(e, g));
  }

  edge_descriptor operator[](size_type p)
  {
    edge_descriptor e = base_iter[p];
    return is_valid(e) ? e : null_edge(g);
  }

public:
  vertex_filter vfilter;
  edge_filter efilter;
  out_edge_iterator base_iter;
  const graph_adapter& g;
};

/// This struct is used to define the bidirectional iterators as void for
/// the undirectedS and directedS case.
template <typename direction, typename graph_adapter, typename edge_filter,
          typename vertex_filter>
struct get_iterator_type {
  typedef void in_adjacency_iterator;
  typedef void in_edge_iterator;
};

/// This partial specialization is used to define the bidirectional iterators
/// for the bidirectionalS case.
template <typename graph_adapter, typename edge_filter, typename vertex_filter>
struct get_iterator_type<bidirectionalS, graph_adapter,
                         edge_filter, vertex_filter> {
  typedef filter_vertex_iterator<
              vertex_filter,
              typename graph_traits<graph_adapter>::in_adjacency_iterator,
              graph_adapter>
          in_adjacency_iterator;

  typedef filter_in_edge_iterator<vertex_filter, edge_filter, graph_adapter>
          in_edge_iterator;
};

}

template <typename graph_adapter, typename edge_filter,
          typename vertex_filter = keep_all>
class filter_adapter {
public:
  typedef graph_traits<graph_adapter> traits;
  typedef typename traits::size_type size_type;
  typedef typename traits::vertex_descriptor vertex_descriptor;
  typedef typename traits::edge_descriptor edge_descriptor;

  typedef detail::filter_vertex_iterator<vertex_filter,
                                         typename traits::vertex_iterator,
                                         graph_adapter>
          vertex_iterator;

  typedef detail::filter_vertex_iterator<vertex_filter,
                                         typename traits::adjacency_iterator,
                                         graph_adapter>
          adjacency_iterator;

  typedef typename detail::get_iterator_type<
                       typename traits::directed_category,
                       graph_adapter,
                       edge_filter,
                       vertex_filter>::in_adjacency_iterator
          in_adjacency_iterator;

  typedef detail::filter_edge_iterator<vertex_filter, edge_filter,
                                       graph_adapter>
          edge_iterator;

  typedef detail::filter_out_edge_iterator<vertex_filter, edge_filter,
                                           graph_adapter>
          out_edge_iterator;

  typedef typename detail::get_iterator_type<
                       typename traits::directed_category,
                       graph_adapter,
                       edge_filter,
                       vertex_filter>::in_edge_iterator
          in_edge_iterator;

  typedef void thread_vertex_iterator;
  typedef void thread_adjacency_iterator;
  typedef void thread_in_adjacency_iterator;
  typedef void thread_edge_iterator;
  typedef void thread_out_edge_iterator;
  typedef void thread_in_edge_iterator;

  typedef typename traits::directed_category directed_category;
  typedef typename traits::iterator_category iterator_category;

  filter_adapter(graph_adapter& g, edge_filter ef) :
    original_graph(g), efilt(ef) {}

  filter_adapter(graph_adapter& g, edge_filter ef, vertex_filter vf) :
    original_graph(g), efilt(ef), vfilt(vf) {}

  // The compiler-synthesized copy control works fine for this class, so we
  // don't implement it ourselves.  We assume that the graph adapter passed
  // as a template parameter has a correctly implemented deep copy which makes
  // the transpose adapter perform a deep copy.

  const graph_adapter& get_original_graph() const { return original_graph; }

public:
  const graph_adapter& original_graph;
  edge_filter efilt;
  vertex_filter vfilt;
};

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::size_type
num_vertices(const filter_adapter<G, EF, VF>& fg)
{
  return num_vertices(fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::size_type
num_edges(const filter_adapter<G, EF, VF>& fg)
{
  return num_edges(fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::vertex_descriptor
source(const typename filter_adapter<G, EF, VF>::edge_descriptor& e,
       const filter_adapter<G, EF, VF>& fg)
{
  return source(e, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::vertex_descriptor
target(const typename filter_adapter<G, EF, VF>::edge_descriptor& e,
       const filter_adapter<G, EF, VF>& fg)
{
  return target(e, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::size_type
degree(const typename filter_adapter<G, EF, VF>::vertex_descriptor& u_local,
       const filter_adapter<G, EF, VF>& fg)
{
  return degree(u_local, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::size_type
out_degree(
    const typename filter_adapter<G, EF, VF>::vertex_descriptor& u_local,
    const filter_adapter<G, EF, VF>& fg)
{
  return out_degree(u_local, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::size_type
in_degree(const typename filter_adapter<G, EF, VF>::vertex_descriptor& u_local,
          const filter_adapter<G, EF, VF>& fg)
{
  return in_degree(u_local, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::vertex_iterator
vertices(const filter_adapter<G, EF, VF>& fg)
{
  typedef typename filter_adapter<G, EF, VF>::vertex_iterator iterator;
  typename graph_traits<G>::vertex_iterator viter =
    vertices(fg.get_original_graph());
  return iterator(fg.vfilt, viter, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::edge_iterator
edges(const filter_adapter<G, EF, VF>& fg)
{
  typedef typename filter_adapter<G, EF, VF>::edge_iterator iterator;
  typename graph_traits<G>::edge_iterator eiter =
    edges(fg.get_original_graph());
  return iterator(fg.vfilt, fg.efilt, eiter, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::adjacency_iterator
adjacent_vertices(
    const typename filter_adapter<G, EF, VF>::vertex_descriptor& v,
    const filter_adapter<G, EF, VF>& fg)
{
  typedef typename filter_adapter<G, EF, VF>::adjacency_iterator iterator;
  typename graph_traits<G>::adjacency_iterator viter =
    adjacent_vertices(v, fg.get_original_graph());
  return iterator(fg.vfilt, viter, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::in_adjacency_iterator
in_adjacent_vertices(
    const typename filter_adapter<G, EF, VF>::vertex_descriptor& v,
    const filter_adapter<G, EF, VF>& fg)
{
  typedef typename filter_adapter<G, EF, VF>::in_adjacency_iterator iterator;
  typename graph_traits<G>::in_adjacency_iterator viter =
    in_adjacent_vertices(v, fg.get_original_graph());
  return iterator(fg.vfilt, viter, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::out_edge_iterator
out_edges(const typename filter_adapter<G, EF, VF>::vertex_descriptor& v,
          const filter_adapter<G, EF, VF>& fg)
{
  typedef typename filter_adapter<G, EF, VF>::out_edge_iterator iterator;
  typename graph_traits<G>::out_edge_iterator eiter =
    out_edges(v, fg.get_original_graph());
  return iterator(fg.vfilt, fg.efilt, eiter, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::in_edge_iterator
in_edges(const typename filter_adapter<G, EF, VF>::vertex_descriptor& v,
         const filter_adapter<G, EF, VF>& fg)
{
  typedef typename filter_adapter<G, EF, VF>::in_edge_iterator iterator;
  typename graph_traits<G>::in_edge_iterator eiter =
    in_edges(v, fg.get_original_graph());
  return iterator(fg.vfilt, fg.efilt, eiter, fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::vertex_descriptor
null_vertex(const filter_adapter<G, EF, VF>& fg)
{
  return null_vertex(fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
typename filter_adapter<G, EF, VF>::edge_descriptor
null_edge(const filter_adapter<G, EF, VF>& fg)
{
  return null_edge(fg.get_original_graph());
}

/***/

template <typename ITERATOR, typename G, typename EF, typename VF>
inline
bool
is_valid(ITERATOR& iter, typename filter_adapter<G, EF, VF>::size_type p,
         const filter_adapter<G, EF, VF>& fg)
{
  return iter.is_valid(iter.base_iter[p]);
}

/***/

template <typename G, typename EF, typename VF>
inline
bool
is_directed(const filter_adapter<G, EF, VF>& fg)
{
  return is_directed(fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
bool
is_undirected(const filter_adapter<G, EF, VF>& fg)
{
  return is_undirected(fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
inline
bool
is_bidirectional(const filter_adapter<G, EF, VF>& fg)
{
  return is_bidirectional(fg.get_original_graph());
}

/***/

template <typename G, typename EF, typename VF>
class vertex_id_map<filter_adapter<G, EF, VF> > : public vertex_id_map<G> {
public:
  vertex_id_map() : vertex_id_map<G>() {}
  vertex_id_map(const vertex_id_map& vm) : vertex_id_map<G>(vm) {}
};

/***/

template <typename G, typename EF, typename VF>
class edge_id_map<filter_adapter<G, EF, VF> > : public edge_id_map<G> {
public:
  edge_id_map() : edge_id_map<G>() {}
  edge_id_map(const edge_id_map& em) : edge_id_map<G>(em) {}
};

}

#endif
