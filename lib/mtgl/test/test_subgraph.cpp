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
/*! \file test_subgraph.cpp

    \brief Tests the subgraph adapter.

    \author Greg Mackey (gemacke@sandia.gov)
*/
/****************************************************************************/

#include <iostream>

//#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>

#include <mtgl/subgraph_adapter.hpp>

using namespace mtgl;

typedef directedS DIRECTION;

int main()
{
//  typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;

  typedef subgraph_adapter<Graph> SGraph;
  typedef graph_traits<SGraph>::edge_iterator edge_iterator_sg;

  // Initialize graph.
  Graph g;

  const size_type numVertices = 4;
  const size_type numEdges = 7;

  size_type sources[numEdges] = { 0, 1, 2, 3, 3, 0, 0 };
  size_type targets[numEdges] = { 1, 2, 3, 0, 0, 2, 3 };

  init(numVertices, numEdges, sources, targets, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "Graph: (" << order << ", " << size << ")" << std::endl;
  print(g);

  subgraph_adapter<Graph> sub(g);

  size_type num_subgraph_edges = 3;
  edge_descriptor subEdges[num_subgraph_edges];
  subEdges[0] = edges(g)[0];
  subEdges[1] = edges(g)[3];
  subEdges[2] = edges(g)[4];
//  subEdges[3] = edges(g)[6];

  init_edges(subEdges, num_subgraph_edges, sub);

/*
  size_type num_subgraph_verts = 3;
  vertex_descriptor subVerts[num_subgraph_verts];
  subVerts[0] = vertices(g)[0];
  subVerts[1] = vertices(g)[2];
  subVerts[2] = vertices(g)[3];

  init_vertices(subVerts, num_subgraph_verts, sub);
*/

  std::cout << std::endl << "Subgraph:" << std::endl;
  sub.print();

  std::cout << std::endl << "Edge iterator test:" << std::endl;
  edge_id_map<subgraph_adapter<Graph> > eipm = get(_edge_id_map, sub);
  size_type nse = num_edges(sub);
  edge_iterator_sg sg_eIter = edges(sub);
  for (size_type i = 0; i < nse; ++i)
  {
    edge_descriptor e = sg_eIter[i];
    size_type eid = get(eipm, e);
    std::cout << eid << std::endl;
  }

  subgraph_adapter<Graph> sub_new(g);

  vertex_iterator vIter = vertices(g);
  edge_iterator eIter = edges(g);

/*
  vertex_property_map<Graph, bool> vmask(g);
  vmask[vIter[0]] = true;
  vmask[vIter[1]] = false;
  vmask[vIter[2]] = true;
  vmask[vIter[3]] = true;
*/

  edge_property_map<Graph, bool> emask(g);
  emask[eIter[0]] = true;
  emask[eIter[1]] = false;
  emask[eIter[2]] = false;
  emask[eIter[3]] = true;
  emask[eIter[4]] = true;
  emask[eIter[5]] = false;
  emask[eIter[6]] = false;

  init_edges(emask, sub_new);
//  init_vertices(vmask, sub_new);

  std::cout << std::endl << "Subgraph New:" << std::endl;
  sub_new.print();

  std::cout << std::endl << "Testing copy constructor." << std::endl;
  subgraph_adapter<Graph> sub_copy(sub);
  sub_copy.print();

  std::cout << std::endl << "Testing assignment." << std::endl;
  subgraph_adapter<Graph> sub_assign(g);
  sub_assign = sub;
  sub_assign.print();

  return 0;
}
