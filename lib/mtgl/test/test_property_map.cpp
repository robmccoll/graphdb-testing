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
/*! \file test_property_map.cpp

    \brief Tests property maps in MTGL.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 6/17/2010
*/
/****************************************************************************/

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/transpose_adapter.hpp>

#define WRAPPER

using namespace mtgl;

typedef directedS DIRECTION;

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<DIRECTION> BaseGraph;
//  typedef adjacency_list<DIRECTION> BaseGraph;
#ifdef WRAPPER
  typedef transpose_adapter<BaseGraph> Graph;
#else
  typedef BaseGraph Graph;
#endif
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;
  typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;
  typedef graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  const size_type numVerts = 6;
  const size_type numEdges = 8;

  size_type sources[numEdges] = { 0, 0, 1, 1, 2, 3, 3, 4 };
  size_type targets[numEdges] = { 1, 2, 2, 3, 4, 4, 5, 5 };

  // Initialize the graph.
  BaseGraph bg;
  init(numVerts, numEdges, sources, targets, bg);
#ifdef WRAPPER
  Graph g(bg);
#else
  Graph& g = bg;
#endif

  // Get the id maps.
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> eid_map = get(_edge_id_map, g);

  // Declare the C-style array array property maps.
  int* vprop1 = new int[numVerts];
  int* eprop1 = new int[numEdges];
  array_property_map<int, vertex_id_map<Graph> > vPropmap1(vprop1, vid_map);
  array_property_map<int, edge_id_map<Graph> > ePropmap1(eprop1, eid_map);

  // Declare the dynamic_array array property maps.
  dynamic_array<int> vprop2(numVerts);
  dynamic_array<int> eprop2(numEdges);
  array_property_map<int, vertex_id_map<Graph> > vPropmap2(vprop2, vid_map);
  array_property_map<int, edge_id_map<Graph> > ePropmap2(eprop2, eid_map);

  // Declare the xmt_hash_table map property maps.
  xmt_hash_table<size_type, int> vprop3(2 * numVerts);
  xmt_hash_table<size_type, int> eprop3(2 * numEdges);
  map_property_map<int, vertex_id_map<Graph> > vPropmap3(vprop3, vid_map);
  map_property_map<int, edge_id_map<Graph> > ePropmap3(eprop3, eid_map);

  // Declare the vertex_property_map and edge_property_map.
  vertex_property_map<Graph, int> vPropmap4(g);
  edge_property_map<Graph, int> ePropmap4(g);

  // Initialize the property maps.
  vertex_iterator verts = vertices(g);
  size_type order = num_vertices(g);
  #pragma mta assert parallel
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];

    vPropmap1[u] = i * 2;
    vPropmap2[u] = i * 3;
    vPropmap3[u] = i * 4;
    vPropmap4[u] = i * 5;
//    put(vPropmap1, u, i * 2);
//    put(vPropmap2, u, i * 3);
//    put(vPropmap3, u, i * 4);
//    put(vPropmap4, u, i * 5);
  }

  // Test mt_incr(), mt_readfe(), and mt_write().
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];

//    int val = mt_readfe(vPropmap2[u]);
//    printf("    %d\n", val);
//    mt_write(vPropmap2[u], val + 2);

    mt_incr(vPropmap3[u], 2);
  }

  edge_iterator edgs = edges(g);
  size_type size = num_edges(g);
  #pragma mta assert parallel
  for (size_type i = 0; i < size; i++)
  {
    edge_descriptor e = edgs[i];

    ePropmap1[e] = i * 2;
    ePropmap2[e] = i * 3;
    ePropmap3[e] = i * 4;
    ePropmap4[e] = i * 5;
//    put(ePropmap1, e, i * 2);
//    put(ePropmap2, e, i * 3);
//    put(ePropmap3, e, i * 4);
//    put(ePropmap4, e, i * 5);
  }

  size_type id_sum = 0;

  // Test vertex id map.
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      vertex_descriptor v = adjs[j];
      size_type vid = vid_map[v];
//      size_type vid = get(vid_map, v);
      id_sum += vid;
    }
  }

  printf("  Sum of adjacency ids: %lu\n", id_sum);

  int prop_sum = 0;

  // Test vertex array property map.
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      vertex_descriptor v = adjs[j];
      int vtype = vPropmap1[v];
//      int vtype = get(vPropmap1, v);
      prop_sum += vtype;
    }
  }

  printf("Sum of adjacency prop1: %d\n", prop_sum);

  prop_sum = 0;

  // Test vertex dynamic_array property map.
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      vertex_descriptor v = adjs[j];
      int vtype = vPropmap2[v];
//      int vtype = get(vPropmap2, v);
      prop_sum += vtype;
    }
  }

  printf("Sum of adjacency prop2: %d\n", prop_sum);

  prop_sum = 0;

  // Test vertex xmt_hash_table property map.
  #pragma mta assert parallel
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    #pragma mta assert parallel
    for (size_type j = 0; j < end; j++)
    {
      vertex_descriptor v = adjs[j];
      int vtype = vPropmap3[v];
//      int vtype = get(vPropmap3, v);
      prop_sum += vtype;
    }
  }

  printf("Sum of adjacency prop3: %d\n", prop_sum);

  prop_sum = 0;

  // Test vertex_property_map.
  #pragma mta assert parallel
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    #pragma mta assert parallel
    for (size_type j = 0; j < end; j++)
    {
      vertex_descriptor v = adjs[j];
      int vtype = vPropmap4[v];
//      int vtype = get(vPropmap4, v);
      prop_sum += vtype;
    }
  }

  printf("Sum of adjacency prop4: %d\n", prop_sum);

  id_sum = 0;

  // Test edge id map.
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    out_edge_iterator oEdges = out_edges(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      edge_descriptor e = oEdges[j];
      size_type eid = eid_map[e];
//      size_type eid = get(eid_map, e);
      id_sum += eid;
    }
  }

  printf("       Sum of edge ids: %lu\n", id_sum);

  prop_sum = 0;

  // Test edge array property map.
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    out_edge_iterator oEdges = out_edges(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      edge_descriptor e = oEdges[j];
      int etype = ePropmap1[e];
//      int etype = get(ePropmap1, e);
      prop_sum += etype;
    }
  }

  printf("     Sum of edge prop1: %d\n", prop_sum);

  prop_sum = 0;

  // Test edge dynamic_array property map.
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    out_edge_iterator oEdges = out_edges(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      edge_descriptor e = oEdges[j];
      int etype = ePropmap2[e];
//      int etype = get(ePropmap2, e);
      prop_sum += etype;
    }
  }

  printf("     Sum of edge prop2: %d\n", prop_sum);

  prop_sum = 0;

  // Test edge xmt_hash_table property map.
  #pragma mta assert parallel
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    out_edge_iterator oEdges = out_edges(u, g);
    size_type end = out_degree(u, g);
    #pragma mta assert parallel
    for (size_type j = 0; j < end; j++)
    {
      edge_descriptor e = oEdges[j];
      int etype = ePropmap3[e];
//      int etype = get(ePropmap3, e);
      prop_sum += etype;
    }
  }

  printf("     Sum of edge prop3: %d\n", prop_sum);

  prop_sum = 0;

  // Test edge_property_map.
  #pragma mta assert parallel
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    out_edge_iterator oEdges = out_edges(u, g);
    size_type end = out_degree(u, g);
    #pragma mta assert parallel
    for (size_type j = 0; j < end; j++)
    {
      edge_descriptor e = oEdges[j];
      int etype = ePropmap4[e];
//      int etype = get(ePropmap4, e);
      prop_sum += etype;
    }
  }

  printf("     Sum of edge prop4: %d\n", prop_sum);

  delete [] vprop1;
  delete [] eprop1;

  return 0;
}
