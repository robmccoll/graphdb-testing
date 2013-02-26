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
/*! \file euler_tour.hpp

    \brief Creates an Euler tour through the graph using Hierholzer's
           algorithm.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 4/2008 (last revised: 10/2011)

    This algorithm works for both directed and undirected graphs.  The
    algorithm assumes that the graph is Eulerian (i.e. an Euler tour exists).
    For undirected graphs this requires that the graph is connected and every
    vertex has even degree.  For directed graphs, this requires that the graph
    is strongly connected and the indegree of every vertex is equal to its
    outdegree.  The algorithm returns a boolean indicating if it was
    successful or not (i.e. if the graph is eulerian or not).  If the graph
    isn't eulerian, the returned path is invalid.

    The path returned by the algorithm is two arrays, one of vertex_descriptors
    and one of edge_descriptors.  They represent the vertices and edges
    encountered in the path.  A path with n vertices is represented by

      (vertices[0], edges[0], vertices[1], ... , edges[n-2] , vertices[n-1])

    The algorithm has a parameter num_tours that specifies how many euler
    tours to take of the graph.

    Any connected undirected graph can be made Eulerian by creating a new
    graph with two copies of each edge in the original graph.  Any weakly
    connected directed graph can be made Eulerian by creating a new graph
    that contains each edge from the original graph and its reverse.  The
    duplicate_adapter is an easy way to do this for both undirected and
    directed graphs.
*/
/****************************************************************************/

#ifndef MTGL_EULER_TOUR_HPP
#define MTGL_EULER_TOUR_HPP

#include <limits>
#include <mtgl/random.hpp>
#include <mtgl/duplicate_adapter.hpp>

namespace mtgl {

template <typename Graph>
bool euler_tour(
  Graph& g,
  typename graph_traits<Graph>::vertex_descriptor start_vertex,
  dynamic_array<typename graph_traits<Graph>::vertex_descriptor>& path_verts,
  dynamic_array<typename graph_traits<Graph>::edge_descriptor>& path_edges, 
  typename graph_traits<Graph>::size_type num_tours = 1)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex;
  typedef typename graph_traits<Graph>::edge_descriptor edge;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator
          out_edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> eid_map = get(_edge_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  size_type tour_length = 2 * size;
  size_type last_tour_length = tour_length;

  // Keeps track of the number of edges adjacent to a vertex that have not
  // been used yet in the partial Euler tour.  This lets us easily find a
  // vertex in the partial tour to start the next cycle from.
  //size_type* num_unused_adj_edges = new size_type[order];
  vertex_property_map<Graph, size_type> num_unused_adj_edges(g);

  // Indexed by edge, each entry in the vector contains the next edge in the
  // euler tour.  An identity value indicates that the edge is not yet on the
  // tour.
  edge_property_map<Graph, edge> visited_edges(g);
  //edge* visited_edges = new size_type[size];

  // Unused edges is temporary storage to hold the unused edges adjacent to a
  // vertex.  The number of adjacent edges to a vertex is bounded by the
  // number of edges in the graph.
  edge* unused_edges = new edge[size];
  size_type num_unused_edges;

  edge_iterator edgs = edges(g);
  vertex_iterator verts = vertices(g);

  for (size_type tour_id = 0; tour_id < num_tours; ++tour_id)
  {
    // Reset num_unused_adj_edges and visited_edges.
    for (size_type i = 0; i < order; ++i)
    {
      num_unused_adj_edges[verts[i]] = out_degree(verts[i], g);
    }

    for (size_type i = 0; i < size; ++i) visited_edges[edgs[i]] = edgs[i];

    // Add the first edge to the tour.
    vertex cycle_start_vertex = start_vertex;
    --num_unused_adj_edges[start_vertex];

    out_edge_iterator eit = out_edges(start_vertex, g);
    size_type deg = out_degree(start_vertex, g);

    if (deg == 0) return false;

    // Get a random entry from the out edges of the starting vertex.
    size_type index = static_cast<size_type>(mt_lrand48()) % deg;
    edge cur_edge = eit[index];

    edge start_edge = cur_edge;
    edge next_cycle_edge = start_edge;
    size_type num_added_edges = 1;

    edge prev_edge = start_edge;
    vertex cur_vert = start_vertex == source(cur_edge, g) ?
                      target(cur_edge, g) : source(cur_edge, g);

    if (is_undirected(g)) --num_unused_adj_edges[cur_vert];

    // Now, add the remaining edges to the tour.
    while (num_added_edges < size)
    {
      // Create a cycle.  (The first edge of the cycle has already been
      // added.)
      while (cur_vert != cycle_start_vertex)
      {
        --num_unused_adj_edges[cur_vert];

        // Randomly pick an edge from the unvisited edges of the vertex.
        num_unused_edges = 0;
        size_type deg = out_degree(cur_vert, g);
        out_edge_iterator eit2 = out_edges(cur_vert, g);

        // Put the unvisited edges of the vertex in an array.
        #pragma mta assert nodep
        for (size_type i = 0; i < deg; ++i)
        {
          edge e = eit2[i];
          size_type eid = get(eid_map, e);
          size_type prev_eid = get(eid_map, prev_edge);
          size_type vis_eid = get(eid_map, visited_edges[e]);

          if (vis_eid == eid && eid != prev_eid)
          {
            size_type idx = mt_incr(num_unused_edges, 1);
            unused_edges[idx] = e;
          }
        }

        if (num_unused_edges == 0) return false;

        // Get a random entry into the unvisited edges array.
        index = static_cast<size_type>(mt_lrand48()) % num_unused_edges;
        cur_edge = unused_edges[index];

        visited_edges[prev_edge] = cur_edge;
        prev_edge = cur_edge;
        cur_vert = cur_vert == source(cur_edge, g) ? target(cur_edge, g) :
                                                     source(cur_edge, g);

        if (is_undirected(g)) --num_unused_adj_edges[cur_vert];

        ++num_added_edges;
      }

      // Add this edge to close the loop of the cycle.
      visited_edges[cur_edge] = next_cycle_edge;

      // If you haven't traveled all the edges yet, find a vertex in the
      // cycle that has an adjacent edge that isn't in the cycle.
      if (num_added_edges < size)
      {
        size_type edge_count = 0;

        while (num_unused_adj_edges[cur_vert] == 0)
        {
          // Move the current edge.
          cur_edge = visited_edges[prev_edge];

          // TODO: Is this statement necessary?
          size_type cur_eid = get(eid_map, cur_edge);
          size_type prev_eid = get(eid_map, prev_edge);
          if (cur_eid == prev_eid)
          {
            cur_edge = start_edge;
          }

          // Move to the next vertex and update the previous edge.
          cur_vert = cur_vert == source(cur_edge, g) ? target(cur_edge, g) :
                                                       source(cur_edge, g);

          prev_edge = cur_edge;

          ++edge_count;

          if (edge_count > num_added_edges + 1) return false;
        }

        next_cycle_edge = visited_edges[cur_edge];

        // Add the first edge of the new cycle.
        cycle_start_vertex = cur_vert;
        --num_unused_adj_edges[cycle_start_vertex];

        // Randomly pick an edge from the unvisited edges of the vertex.
        num_unused_edges = 0;
        size_type deg = out_degree(cur_vert, g);
        out_edge_iterator eit2 = out_edges(cur_vert, g);

        // Put the unvisited edges of the vertex in an array.
        #pragma mta assert nodep
        for (size_type i = 0; i < deg; ++i)
        {
          edge e = eit2[i];
          size_type eid = get(eid_map, e);
          size_type prev_eid = get(eid_map, prev_edge);
          size_type vis_eid = get(eid_map, visited_edges[e]);

          if (vis_eid == eid && eid != prev_eid)
          {
            size_type idx = mt_incr(num_unused_edges, 1);
            unused_edges[idx] = e;
          }
        }

        if (num_unused_edges == 0) return false;

        // Get a random entry into the unvisited edges array.
        index = static_cast<size_type>(mt_lrand48()) % num_unused_edges;
        cur_edge = unused_edges[index];

        visited_edges[prev_edge] = cur_edge;
        prev_edge = cur_edge;
        ++num_added_edges;
        cur_vert = cur_vert == source(cur_edge, g) ? target(cur_edge, g) :
                                                     source(cur_edge, g);

        if (is_undirected(g)) --num_unused_adj_edges[cur_vert];
      }
    }

    // Travel the tour to fill the arrays that are returned with the result.

    path_verts.push_back(start_vertex);

    cur_edge = start_edge;
    cur_vert = start_vertex == source(cur_edge, g) ? target(cur_edge, g) :
                                                     source(cur_edge, g);

    size_type this_tour_length = tour_id + 1 == num_tours ? last_tour_length :
                                 tour_length;

    for (size_type i = 2; i <= this_tour_length; i += 2)
    {
      path_edges.push_back(cur_edge);
      path_verts.push_back(cur_vert);

      // Move edge and vertex forward.
      cur_edge = visited_edges[cur_edge];

      cur_vert = (cur_vert == source(cur_edge, g)) ? target(cur_edge, g) :
                                                     source(cur_edge, g);
    }
  }

  delete [] unused_edges;

  return true;
}

template <typename Graph>
bool has_eulerian_circuit(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_iterator verts = vertices(g);
  edge_iterator edgs = edges(g);
  size_type order = num_vertices(g);
  size_type size  = num_edges(g);

  bool eulerian = true;

  // This function doesn't check that an undirected graph is connected or a
  // directed graph is strongly connected.  These are also requirements for
  // a graph to be eulerian.

  if (is_undirected(g))
  {
    // An undirected graph is not eulerian unless the degree of every vertex
    // is even.
    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i)
    {
      if (out_degree(verts[i], g) % 2 == 1) eulerian = false;
    }
  }
  else if (is_directed(g))
  {
    vertex_property_map<Graph, size_type> indeg(g);

    // Initialize the indegree to 0.
    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i) put(indeg, verts[i], 0);

    // Count the indegree.
    #pragma mta assert nodep
    for (size_type i = 0; i < size; ++i) mt_incr(indeg[target(edgs[i], g)], 1);

    // A directed graph is not eulerian unless the outdegree and indegree are
    // equal for all vertices.
    #pragma mta assert nodep
    for (size_type i = 0; i < order; ++i)
    {
      if (get(indeg, verts[i]) != out_degree(verts[i], g)) eulerian = false;
    }
  }

  return eulerian;
}

template <typename Graph>
bool duplicate_euler_tour(
  Graph& g,
  typename graph_traits<Graph>::vertex_descriptor start_vertex,
  dynamic_array<typename graph_traits<Graph>::vertex_descriptor>& path_verts,
  dynamic_array<typename graph_traits<Graph>::edge_descriptor>& path_edges, 
  typename graph_traits<Graph>::size_type num_tours = 1)
{
  typedef duplicate_adapter<Graph> DGraph;
  typedef typename graph_traits<DGraph>::size_type size_type_d; 
  typedef typename graph_traits<DGraph>::vertex_descriptor vertex_descriptor_d; 
  typedef typename graph_traits<DGraph>::edge_descriptor edge_descriptor_d; 

  if (has_eulerian_circuit(g))
  { 
    bool success = euler_tour(g, start_vertex, path_verts, path_edges); 
    assert(success);
    return success;
  }

  // The duplicate_adapter always has dense edge ids.
  DGraph dg(g);

#ifdef DEBUG
  printf("\nThe duplicate graph:\n");
  print(dg);
  printf("Size of duplicate graph: %d\n", num_edges(dg));
#endif

  dynamic_array<vertex_descriptor_d> path_verts_d;
  dynamic_array<edge_descriptor_d> path_edges_d;

  vertex_descriptor_d start_vertex_d = global_to_local(start_vertex, dg);

  bool success = euler_tour(dg, start_vertex_d, path_verts_d, path_edges_d);
  assert(success);

  // Convert back to original graph descriptors.
  path_verts.clear();
  path_edges.clear();

  size_type_d num_path_verts = path_verts_d.size();
  for (size_type_d i = 0; i < num_path_verts; i += 1)
  {
    path_verts.push_back(local_to_global(path_verts_d[i], dg));
  }

  size_type_d num_path_edges = path_edges_d.size();
  for (size_type_d i = 0; i < num_path_edges; i += 1)
  {
    path_edges.push_back(local_to_global(path_edges_d[i], dg));
  }

  return success;
}

}

#endif
