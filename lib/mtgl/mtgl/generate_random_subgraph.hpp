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
/*! \file generate_random_subgraph.hpp

    \brief Randomly generates a subgraph of a large graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 9/11/2008

    This algorithm generates an edge induced subgraph that has no more than m
    edges.  The algorithm attempts to limit the number of vertices to n, but
    a slightly higher number of vertices could be included in the subgraph.

    The algorithm uses a psearch to create a list of edges and vertices
    included in the subgraph.  The number of vertices returned by the psearch
    will be n or less, and the number of edges returned by the psearch will be
    m or less.  An edge induced subgraph is created from the returned edges.
    The psearch allows for an edge to be included with one endpoint not
    included in the vertex list.  In this case, the extra vertex is added to
    the subgraph increasing the number of vertices.

    This algorithm visits edges based on DIRECTED travel.  From a vertex
    in a directed or bidirectional graph, all out edges are visited.  The
    result for a directed and a bidirectional graph should be the same when
    the algorithm is run in serial.  From a vertex in an undirected graph,
    all adjacent edges are visited which means that each edge is visited twice,
    once from each endpoint.
*/
/****************************************************************************/

#ifndef MTGL_GENERATE_RANDOM_SUBGRAPH_HPP
#define MTGL_GENERATE_RANDOM_SUBGRAPH_HPP

#include <mtgl/psearch.hpp>
#include <mtgl/subgraph_adapter.hpp>

namespace mtgl {

namespace detail {

template <typename Graph>
class generate_random_subgraph_visitor : public default_psearch_visitor<Graph>
{
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex;
  typedef typename graph_traits<Graph>::edge_descriptor edge;
  typedef typename graph_traits<Graph>::size_type size_type;

  generate_random_subgraph_visitor(size_type n, size_type m,
                                   vertex_property_map<Graph, bool>& vmsk,
                                   edge_property_map<Graph, bool>& emsk,
                                   Graph& graph, size_type& cv,
                                   size_type& ce) :
    desired_verts(n), desired_edges(m), vmask(vmsk), emask(emsk),
    created_verts(cv), created_edges(ce), g(graph) {}

  void pre_visit(vertex& v)
  {
    size_type next = mt_incr(created_verts, 1);

    if (next < desired_verts)  vmask[v] = true;
  }

  int visit_test(edge& e, vertex& v)
  {
    size_type vert_count = created_verts + (vmask[target(e, g)] ? 0 : 1);

    return (vmask[v] && vert_count <= desired_verts &&
            created_edges < desired_edges);
  }

  void tree_edge(edge& e, vertex& v)
  {
    if (emask[e])  return;

    size_type next = mt_incr(created_edges, 1);

    if (next < desired_edges)  emask[e] = true;
  }

  void back_edge(edge& e, vertex& v)
  {
    if (emask[e])  return;

    size_type next = mt_incr(created_edges, 1);

    if (next < desired_edges)  emask[e] = true;
  }

private:
  size_type desired_verts;
  size_type desired_edges;
  vertex_property_map<Graph, bool>& vmask;
  edge_property_map<Graph, bool>& emask;
  size_type& created_verts;
  size_type& created_edges;
  Graph& g;
};

}

/***/

template <typename Graph>
void
generate_random_subgraph(typename graph_traits<Graph>::size_type desired_verts,
                         typename graph_traits<Graph>::size_type desired_edges,
                         Graph& g, subgraph_adapter<Graph>& subg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  vertex_iterator verts = vertices(g);
  edge_iterator edgs = edges(g);

  // The bitmap indicating which vertices are in the subgraph.
  vertex_property_map<Graph, bool> subgraph_verts(g);
  for (size_type i = 0; i < size; ++i) subgraph_verts[verts[i]] = false;

  // The bitmap indicating which edges are in the subgraph.
  edge_property_map<Graph, bool> subgraph_edges(g);
  for (size_type i = 0; i < size; ++i) subgraph_edges[edgs[i]] = false;

  // These are counters that are shared by all the copies of the visitors.
  size_type created_verts = 0;
  size_type created_edges = 0;

  detail::generate_random_subgraph_visitor<Graph>
      gen_sg_vis(desired_verts, desired_edges, subgraph_verts,
                 subgraph_edges, g, created_verts, created_edges);

  psearch_high_low<Graph, detail::generate_random_subgraph_visitor<Graph>,
                   AND_FILTER, DIRECTED> psrch(g, gen_sg_vis, 1);
  psrch.run();

  init_edges(subgraph_edges, subg);

  delete [] subgraph_verts;
}

}

#endif
