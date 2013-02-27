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
/*! \file edmonds_karp_max_flow.hpp

    \brief Performs a parallel implementation of the Edmonds-Karp algorithm.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 8/2/2011

    This is a parallel implementation of the Edmonds-Karp algorithm based
    on the implementation used in the Boost Graph Library.
*/
/****************************************************************************/

#ifndef MTGL_EDMONDS_KARP_MAX_FLOW_HPP
#define MTGL_EDMONDS_KARP_MAX_FLOW_HPP

#include <cstdio>
#include <limits>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/util.hpp>
#include <mtgl/breadth_first_search.hpp>
#include <mtgl/duplicate_adapter.hpp>

#define BFS_FORWARD 0
#define BFS_BACKWARD 1

namespace mtgl {

#ifdef DEBUG
unsigned long vertices_examined;
unsigned long edges_tested;
unsigned long edges_visited;
#endif

namespace detail {

template <typename Graph, typename ResCapMap, typename PredEdgeMap,
          typename ColorMap>
class bfs_s_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::size_type size_type;

  bfs_s_visitor(ResCapMap& rc, PredEdgeMap& p, ColorMap& c, ColorMap& oc,
                size_type& pf, vertex_descriptor& sv, size_type os) :
    residual_capacity(rc), parents(p), color(c), other_color(oc),
    path_found(pf), stop_vertex(sv), orig_size(os) {}

  void examine_vertex(vertex_descriptor& v, Graph& g)
  {
#ifdef DEBUG
    mt_incr(vertices_examined, 1);
#endif
  }

  bool visit_test(bool isWhite, edge_descriptor& e, Graph &g)
  {
#ifdef DEBUG
    mt_incr(edges_tested, 1);
#endif

    if (get(residual_capacity, e) > 0)
    {
      size_type my_color = mt_incr(color[target(e, g)], 1);
      return my_color == 0;
    }

    return false;
  }

  void tree_edge(edge_descriptor& e, Graph& g)
  {
#ifdef DEBUG
    mt_incr(edges_visited, 1);
#endif

    if (get(other_color, target(e, g)) > 0)
    {
      size_type val = mt_incr(path_found, 1);
      if (val == 0) stop_vertex = target(e, g);
    }

    put(parents, target(e, g), get(get(_edge_id_map, g), e));
  }

  ResCapMap& residual_capacity;
  PredEdgeMap& parents;
  ColorMap& color;
  ColorMap& other_color;
  size_type& path_found;
  vertex_descriptor& stop_vertex;
  size_type orig_size;
};

template <typename Graph, typename ResCapMap, typename PredEdgeMap,
          typename ColorMap>
class bfs_t_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::size_type size_type;

  bfs_t_visitor(ResCapMap& rc, PredEdgeMap& p, ColorMap& c, ColorMap& oc,
                size_type& pf, vertex_descriptor& sv, size_type os) :
    residual_capacity(rc), parents(p), color(c), other_color(oc),
    path_found(pf), stop_vertex(sv), orig_size(os) {}

  void examine_vertex(vertex_descriptor& v, Graph& g)
  {
#ifdef DEBUG
    mt_incr(vertices_examined, 1);
#endif
  }

  bool visit_test(bool isWhite, edge_descriptor& e, Graph &g)
  {
#ifdef DEBUG
    mt_incr(edges_tested, 1);
#endif

    size_type eid = get(get(_edge_id_map, g), e);
    size_type rev_eid = eid < orig_size ? eid + orig_size : eid - orig_size;

    if (get(residual_capacity, edges(g)[rev_eid]) > 0)
    {
      size_type my_color = mt_incr(color[target(e, g)], 1);
      return my_color == 0;
    }

    return false;
  }

  void tree_edge(edge_descriptor& e, Graph& g)
  {
#ifdef DEBUG
    mt_incr(edges_visited, 1);
#endif

    if (get(other_color, target(e, g)) > 0)
    {
      size_type val = mt_incr(path_found, 1);
      if (val == 0) stop_vertex = target(e, g);
    }

    size_type eid = get(get(_edge_id_map, g), e);
    size_type rev_eid = eid < orig_size ? eid + orig_size : eid - orig_size;
    put(parents, target(e, g), rev_eid);
  }

  ResCapMap& residual_capacity;
  PredEdgeMap& parents;
  ColorMap& color;
  ColorMap& other_color;
  size_type& path_found;
  vertex_descriptor& stop_vertex;
  size_type orig_size;
};

}

template <typename Graph, typename FlowType>
FlowType
edmonds_karp_max_flow(
    Graph& orig_g,
    typename graph_traits<Graph>::vertex_descriptor s,
    typename graph_traits<Graph>::vertex_descriptor t,
    FlowType* capacity)
{
  #pragma mta noalias orig_g
  #pragma mta noalias *capacity

  typedef duplicate_adapter<Graph> DGraph;
  typedef typename graph_traits<DGraph>::size_type size_type;
  typedef typename graph_traits<DGraph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<DGraph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<DGraph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<DGraph>::edge_iterator edge_iterator;
  typedef typename graph_traits<DGraph>::out_edge_iterator out_edge_iterator;

  size_type src_id = get(get(_vertex_id_map, orig_g), s);
  size_type sink_id = get(get(_vertex_id_map, orig_g), t);

  if (src_id == sink_id) return static_cast<FlowType>(0);

  DGraph g(orig_g);

#ifdef DEBUG
  printf("DGraph: (%lu, %lu)\n\n", num_vertices(g), num_edges(g));
  unsigned long total_vertices_examined = 0;
  unsigned long total_vertices_in_Q = 0;
  unsigned long total_edges_tested = 0;
  unsigned long total_edges_visited = 0;
#endif

  size_type order = num_vertices(g);
  size_type orig_size = num_edges(g) / 2;

  vertex_iterator verts = vertices(g);
  edge_iterator edgs = edges(g);
  vertex_descriptor src = verts[src_id];
  vertex_descriptor sink = verts[sink_id];

  edge_property_map<DGraph, FlowType> residual_capacity(g);
  #pragma mta assert nodep
  for (size_type i = 0; i < orig_size; ++i)
  {
    put(residual_capacity, edgs[i], capacity[i]);
    put(residual_capacity, edgs[i + orig_size], 0);
  }

  dynamic_array<size_type> path;

  size_type path_found = 1;
  vertex_descriptor stop_vertex = null_vertex(g);

  // Declare the variables used by the forward search.
  vertex_property_map<DGraph, size_type> s_parents(g);

  vertex_property_map<DGraph, size_type> s_color(g);
  put(s_color, sink, 1);

  vertex_property_map<DGraph, size_type> s_bfs_color(g);

  vertex_descriptor* s_bfs_Q =
    (vertex_descriptor*) malloc(num_vertices(g) * sizeof(vertex_descriptor));

  vertex_descriptor* s_bfs_buffer = 0;
  size_type s_bfs_buf_size = 0;
  size_type* s_bfs_accum_deg = 0;
  size_type s_bfs_accum_deg_size = 0;

  // Declare the variables used by the backward search.
  vertex_property_map<DGraph, size_type> t_parents(g);

  vertex_property_map<DGraph, size_type> t_color(g);
  put(t_color, sink, 1);

  vertex_property_map<DGraph, size_type> t_bfs_color(g);

  vertex_descriptor* t_bfs_Q =
    (vertex_descriptor*) malloc(num_vertices(g) * sizeof(vertex_descriptor));

  vertex_descriptor* t_bfs_buffer = 0;
  size_type t_bfs_buf_size = 0;
  size_type* t_bfs_accum_deg = 0;
  size_type t_bfs_accum_deg_size = 0;

  // Declare the visitors used by both search directions.
  detail::bfs_s_visitor<DGraph, edge_property_map<DGraph, FlowType>,
                        vertex_property_map<DGraph, size_type>,
                        vertex_property_map<DGraph, size_type> >
    s_bpv(residual_capacity, s_parents, s_color, t_color, path_found,
          stop_vertex, orig_size);

  detail::bfs_t_visitor<DGraph, edge_property_map<DGraph, FlowType>,
                        vertex_property_map<DGraph, size_type>,
                        vertex_property_map<DGraph, size_type> >
    t_bpv(residual_capacity, t_parents, t_color, s_color, path_found,
          stop_vertex, orig_size);

  while (path_found > 0)
  {
    // Do an st search using two BFS's, one starting at the source and the
    // other starting at the sink, to find the next path to augment flow.

    // Do a BFS to find the path from the source to the sink.  We do the
    // code for the BFS (calling expand_one_edge()) directly because we
    // want to stop the BFS once the sink has been reached.

#ifdef DEBUG
    size_type maxDist = 0;
    vertices_examined = 0;
    edges_tested = 0;
    edges_visited = 0;
#endif

    path_found = 0;
    int bfs_direction = BFS_FORWARD;

    // Initialize the variables used by the forward search.
    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i) put(s_color, verts[i], 0);
    put(s_color, src, 1);

    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i) put(s_bfs_color, verts[i], 0);

    size_type s_head = 0;
    size_type s_tail = 1;

    // For the forward search, put the root in the queue and set its color.
    s_bpv.discover_vertex(src, g);
    s_bfs_Q[s_head] = src;
    put(s_bfs_color, src, 1);

    // Initialize the variables used by the backward search.
    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i) put(t_color, verts[i], 0);
    put(t_color, sink, 1);

    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i) put(t_bfs_color, verts[i], 0);

    size_type t_head = 0;
    size_type t_tail = 1;

    // For the backward search, put the root in the queue and set its color.
    t_bpv.discover_vertex(sink, g);
    t_bfs_Q[t_head] = sink;
    put(t_bfs_color, sink, 1);

    while (path_found == 0 && (s_tail > s_head || t_tail > t_head))
    {
      if (bfs_direction == BFS_FORWARD)
      {
        expand_one_edge(g, s_bfs_Q, s_head, s_tail, s_bpv, s_bfs_color,
                        s_bfs_buffer, s_bfs_buf_size, s_bfs_accum_deg,
                        s_bfs_accum_deg_size);

        if (t_tail > t_head) bfs_direction = BFS_BACKWARD;
      }
      else
      {
        expand_one_edge(g, t_bfs_Q, t_head, t_tail, t_bpv, t_bfs_color,
                        t_bfs_buffer, t_bfs_buf_size, t_bfs_accum_deg,
                        t_bfs_accum_deg_size);

        if (s_tail > s_head) bfs_direction = BFS_FORWARD;
      }

#ifdef DEBUG
      ++maxDist;
#endif
    }

#ifdef DEBUG
    total_vertices_examined += vertices_examined;
    total_vertices_in_Q += s_tail + t_tail;
    total_edges_tested += edges_tested;
    total_edges_visited += edges_visited;

    printf("Number of levels: %lu\n", maxDist);
    printf("    Vertices in Q: %lu\n", s_tail + t_tail);
    printf("Examined vertices: %lu\n", vertices_examined);
    printf("    Tested edges: %lu\n", edges_tested);
    printf("   Visited edges: %lu\n", edges_visited);
    printf("\n");
#endif

    if (path_found > 0)
    {
      path.clear();

      // Get an array that contains the edges in the path found by the BFS
      // from source to sink.  This allows the next two loops to be performed
      // in parallel.

      // First, get the edges from the forward search.
      vertex_descriptor u = stop_vertex;
      while (u != src)
      {
        path.push_back(get(s_parents, u));
        u = source(edgs[get(s_parents, u)], g);
      }

      // Now, get the edges from the reverse search.
      u = stop_vertex;
      while (u != sink)
      {
        path.push_back(get(t_parents, u));
        u = target(edgs[get(t_parents, u)], g);
      }

      // Find the minimum residual capacity in the path.
      FlowType min_residual_capacity = (std::numeric_limits<FlowType>::max)();
      size_type path_size = path.size();
      for (size_type i = 0; i < path_size; ++i)
      {
        FlowType rescap_e = get(residual_capacity, edgs[path[i]]);
        if (rescap_e < min_residual_capacity) min_residual_capacity = rescap_e;
      }

      // Push the minimum residual capacity along the path.
      #pragma mta assert nodep
      for (size_type i = 0; i < path_size; ++i)
      {
        size_type rev_eid = path[i] < orig_size ? path[i] + orig_size :
                                                  path[i] - orig_size;

        residual_capacity[edgs[path[i]]] -= min_residual_capacity;
        residual_capacity[edgs[rev_eid]] += min_residual_capacity;
      }
    }
  }

#ifdef DEBUG
    printf("    Total Vertices in Q: %lu\n", total_vertices_in_Q);
    printf("Total Vertices Examined: %lu\n", total_vertices_examined);
    printf("    Total tested edges: %lu\n", total_edges_tested);
    printf("   Total visited edges: %lu\n", total_edges_visited);
    printf("\n");
#endif

  out_edge_iterator oedgs = out_edges(src, g);
  size_type odeg = out_degree(src, g);

  FlowType total_flow = 0;
  for (size_type i = 0; i < odeg; ++i)
  {
    edge_descriptor e = oedgs[i];
    size_type eid = get(get(_edge_id_map, g), e);

    if (eid < orig_size)
    {
      total_flow += capacity[eid] - residual_capacity[e];
    }
    else
    {
      total_flow -= residual_capacity[e];
    }
  }

  free(s_bfs_buffer);
  free(s_bfs_accum_deg);
  free(s_bfs_Q);
  free(t_bfs_buffer);
  free(t_bfs_accum_deg);
  free(t_bfs_Q);

  return total_flow;
}

}

#undef BFS_FORWARD
#undef BFS_BACKWARD

#endif
