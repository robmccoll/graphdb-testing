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
/*! \file random_walk.hpp

    \brief Algorithms for finding random and weighted random walks though
           a graph.

    \author Jon Berry (jberry@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/4/2007
*/
/****************************************************************************/

#ifndef MTGL_RANDOM_WALK_HPP
#define MTGL_RANDOM_WALK_HPP

#include <mtgl/util.hpp>
#include <mtgl/random.hpp>

namespace mtgl {

/*! \brief Returns a random walk of specified length through the graph.

    \author Jon Berry (jberry@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/4/2007

    \tparam Graph The type of the graph.
    \param g The graph.
    \param src The starting vertex for the walk.
    \param num_path_verts The desired path length which is the sum of the
                          number of vertices and edges in the walk.
    \param path_verts The vertices in the resulting path.
    \param path_edges The edges in the resulting path.
    \return The number of vertices in the resulting path.  If this is less
            than num_path_verts, the algorithm encountered a sink.

    The algorithm starts at the given vertex and randomly chooses one of its
    out edges as the next edge in the walk.  The algorithm moves to the chosen
    edge's other vertex and repeats the process until the walk has reached the
    desired length.  If a sink is reached before the desired length is
    reached, the algorithm returns the path found up to that point.

    For the algorithm to be able to visit all vertices in an undirected graph,
    the graph must be connected.  If the graph is not connected, the algorithm
    will still complete, but it will only visit the vertices in the component
    in which it starts.  For the algorithm to be able to visit all vertices in
    a directed or bidirectional graph, the graph must be strongly connected.

    A weakly connected directed graph can be made strongly connected by
    creating a new graph that contains each edge from the original graph and
    its reverse.  In fact, the resulting graph is Eulerian, as well.  The
    duplicate_adapter is an easy way to do this.
*/
template <typename Graph>
typename graph_traits<Graph>::size_type
random_walk(
  Graph& g,
  typename graph_traits<Graph>::vertex_descriptor src,
  typename graph_traits<Graph>::size_type num_path_verts,
  dynamic_array<typename graph_traits<Graph>::vertex_descriptor>& path_verts,
  dynamic_array<typename graph_traits<Graph>::edge_descriptor>& path_edges)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  path_verts.clear();
  path_edges.clear();
  path_verts.reserve(num_path_verts);
  path_edges.reserve(num_path_verts - 1);

  // Put the source vertex in the path.
  path_verts.unsafe_push_back(src);

  size_type visited_verts = 1;
  vertex_descriptor u = src;
  size_type u_deg = out_degree(u, g);

  while (visited_verts < num_path_verts && u_deg > 0)
  {
    out_edge_iterator o_edges = out_edges(u, g);

    // Get a random number in the range [0, u_deg - 1].  The chosen edge is
    // the (randnum)th adjacent edge of u.
    long randnum = mt_lrand48() % u_deg;

    // Set the next edge and vertex in the walk to the chosen edge and its
    // other vertex.
    edge_descriptor e = o_edges[randnum];
    vertex_descriptor v = target(e, g);
    path_edges.unsafe_push_back(e);
    path_verts.unsafe_push_back(v);

    ++visited_verts;
    u = v;
    u_deg = out_degree(u, g);
  }

  return visited_verts;
}

template <typename Graph, typename WeightMap>
void compute_weights_for_random_walk(Graph& g, WeightMap& eweights,
                                     WeightMap& cumulative_eweights)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_iterator verts = vertices(g);
  size_type order = num_vertices(g);
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> eid_map = get(_edge_id_map, g);

  #pragma mta assert parallel
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = verts[i];
    size_type out_deg = out_degree(v, g);
    size_type total = 0;
    out_edge_iterator edgs = out_edges(v, g);
#ifdef DEBUG
    printf("compute_cum: processing vertex %lu\n", get(vid_map, v));
#endif

    // Perhaps do the following with a linear recurrence on the XMT?
    for (size_type j = 0; j < out_deg; ++j)
    {
      edge_descriptor e = edgs[j];
      total += eweights[e];
      cumulative_eweights[e] = total;
#ifdef DEBUG
      printf("compute_cum: set cw[%lu (%lu,%lu)] to %lu\n",
             get(eid_map, e), get(vid_map, source(e, g)), 
             get(vid_map, target(e, g)), cumulative_eweights[e]);
#endif
    }
  }
}


// linear search used for debugging binary search
template <typename Graph, typename WeightMap>
typename graph_traits<Graph>::size_type
select_weighted_out_edge2(Graph& g,
                          typename graph_traits<Graph>::vertex_descriptor u,
                          WeightMap& cumulative_eweights, long rval)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  size_type u_deg = out_degree(u, g);
  size_type total_weight = 0;

  out_edge_iterator edgs = out_edges(u, g);
  edge_descriptor last = edgs[u_deg-1];

  size_type my_val = rval % cumulative_eweights[last];
  size_type winning_adj_id = u_deg;

  for (size_type i = 0; i < u_deg; ++i)
  {
    edge_descriptor e = edgs[i];
    size_type my_sum = cumulative_eweights[e];

    if (my_sum > my_val)
    {
      winning_adj_id = i;
      break;
    }
  }

  edge_descriptor e = edgs[winning_adj_id];
  size_type value = cumulative_eweights[e];

  return winning_adj_id;
}

template <typename Graph, typename WeightMap>
typename graph_traits<Graph>::size_type
select_weighted_out_edge(Graph& g,
                         typename graph_traits<Graph>::vertex_descriptor v,
                         WeightMap& cumulative_eweights, long rval)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  size_type out_deg = out_degree(v, g);

  out_edge_iterator edgs = out_edges(v, g);
  edge_descriptor last = edgs[out_deg - 1];

  size_type my_val = (rval % cumulative_eweights[last]) + 1;

  int64_t lower = 0;
  int64_t upper = out_deg;
  int64_t mid = (lower + upper) / 2;

#ifdef DOUBLE_DEBUG
  printf("wrw2: rval: %ld, cum: %lu, my_rval: %lu\n", rval,
         cumulative_eweights[last], my_val);
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  printf("wrw2: vertex: %lu:  out_deg: %lu\n", get(vid_map, v), out_deg);
  fflush(stdout);
#endif

  while (lower <= upper)
  {
    mid = (lower + upper) / 2;
    edge_descriptor e = edgs[static_cast<size_type>(mid)];
    size_type value = cumulative_eweights[e];

#ifdef DOUBLE_DEBUG
    printf("wrw2: binsrch: %d %d %d\n", lower, mid, upper);
    fflush(stdout);
    printf("wrw2: binsrch: my_val %lu val[mid]: %lu\n", my_val, value);
    fflush(stdout);
#endif

    // non-standard binary search since each cell represents a range.
    // The appropriate range might lie in cell[mid], so if we're cutting
    // off the top range (adjusting upper), we can't lose mid.
    if (my_val < value && lower < upper)
    {
      upper = mid;
    }
    else if (my_val > value)
    {
      lower = mid + 1;
    }
    else
    {
      break;
    }
  }

#ifdef DOUBLE_DEBUG
  edge_descriptor e = edgs[0];
  size_type value = cumulative_eweights[e];
  printf("wrw: (%lu", value); 
  for (size_type i = 1; i<out_deg; i++) {
     e = edgs[i];
     value = cumulative_eweights[e];
     printf(",%lu", value);
  }
  printf("), rval: %lu, ind: %lu\n", my_val, mid);
#endif

  return static_cast<size_type>(mid);
}

/*! \brief Returns a weighted random walk of specified length through the
           graph.

    \author Jon Berry (jberry@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/4/2007

    \tparam Graph The type of the graph.
    \tparam WeightMap The type of property map containing the weights.
    \param g The graph.
    \param src The starting vertex for the walk.
    \param num_path_verts The desired path length which is the sum of the
                          number of vertices and edges in the walk.
    \param path_verts The vertices in the resulting path.
    \param path_edges The edges in the resulting path.
    \param eweights A property map containing the weights for each edge.
    \return The number of vertices in the resulting path.  If this is less
            than num_path_verts, the algorithm encountered a sink.

    The algorithm starts at the given vertex and randomly chooses one of its
    out edges as the next edge in the walk.  The random choice is weighted
    based on the edge weights.  The algorithm moves to the chosen edge's other
    vertex and repeats the process until the walk has reached the desired
    length.  If a sink is reached before the desired length is reached, the
    algorithm returns the path found up to that point.

    For the algorithm to be able to visit all vertices in an undirected graph,
    the graph must be connected.  If the graph is not connected, the algorithm
    will still complete, but it will only visit the vertices in the component
    in which it starts.  For the algorithm to be able to visit all vertices in
    a directed or bidirectional graph, the graph must be strongly connected.

    A weakly connected directed graph can be made strongly connected by
    creating a new graph that contains each edge from the original graph and
    its reverse.  In fact, the resulting graph is Eulerian, as well.  The
    duplicate_adapter is an easy way to do this.
*/
template <typename Graph, typename WeightMap>
typename graph_traits<Graph>::size_type
weighted_random_walk(
  Graph& g, typename graph_traits<Graph>::vertex_descriptor src,
  typename graph_traits<Graph>::size_type num_path_verts,
  dynamic_array<typename graph_traits<Graph>::vertex_descriptor>& path_verts,
  dynamic_array<typename graph_traits<Graph>::edge_descriptor>& path_edges,
  WeightMap& cumulative_eweights, lrand48_generator& rvals,
  typename graph_traits<Graph>::size_type& rval_start)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  path_verts.clear();
  path_edges.clear();
  path_verts.reserve(num_path_verts);
  path_edges.reserve(num_path_verts - 1);

  // Put the source vertex in the path.
  path_verts.unsafe_push_back(src);

  size_type visited_verts = 1;
  vertex_descriptor u = src;
  size_type u_deg = out_degree(u, g);

  while (visited_verts < num_path_verts && u_deg > 0)
  {
    out_edge_iterator o_edges = out_edges(u, g);

    long next_rand = rvals[rval_start];
    size_type winning_adj_id =
      select_weighted_out_edge(g, u, cumulative_eweights, next_rand);

    // Set the next edge and vertex in the walk to the winning adjacent edge
    // and its other vertex.
    edge_descriptor e = o_edges[winning_adj_id];
    vertex_descriptor v = target(e, g);
    path_edges.unsafe_push_back(e);
    path_verts.unsafe_push_back(v);

    ++rval_start;
    ++visited_verts;
    u = v;
    u_deg = out_degree(u, g);
  }

  return visited_verts;
}

}

#endif
