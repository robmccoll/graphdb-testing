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
/*! \file test_duplicate.cpp

    \brief Tests the duplicate adapter.

    \author Greg Mackey (gemacke@sandia.gov)
*/
/****************************************************************************/

#include <iostream>

//#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>

#include <mtgl/subgraph_adapter.hpp>
#include <mtgl/duplicate_adapter.hpp>

using namespace mtgl;

typedef directedS DIRECTION;

int main(int argc, char *argv[])
{
//  typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;

  typedef graph_traits<Graph>::edge_descriptor edge;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<duplicate_adapter<Graph> >::edge_descriptor edge2;

  Graph g;

  // Initialize graph.

#if 0
  const size_type numVerts = 7;
  const size_type numEdges = 5;

  size_type sources[numEdges] = { 0, 1, 2, 4, 5 };
  size_type targets[numEdges] = { 1, 2, 3, 5, 6 };
#endif

#if 0
  const size_type numVerts = 3;
  const size_type numEdges = 3;

  size_type sources[numEdges] = { 0, 1, 1 };
  size_type targets[numEdges] = { 1, 1, 2 };
#endif

#if 0
  const size_type numVerts = 4;
  const size_type numEdges = 3;

  size_type sources[numEdges] = { 0, 1, 2 };
  size_type targets[numEdges] = { 3, 3, 3 };
#endif

#if 0
  const size_type numVerts = 2;
  const size_type numEdges = 1;

  size_type sources[numEdges] = { 0 };
  size_type targets[numEdges] = { 1 };
#endif

#if 0
  const size_type numVerts = 4;
  const size_type numEdges = 3;

  size_type sources[numEdges] = { 0, 1, 2 };
  size_type targets[numEdges] = { 1, 2, 3 };
#endif

#if 0
  const size_type numVerts = 3;
  const size_type numEdges = 2;

  size_type sources[numEdges] = { 0, 1 };
  size_type targets[numEdges] = { 1, 2 };
#endif

  // Test case 1.
  const size_type numVerts = 8;
  const size_type numEdges = 9;

  size_type sources[numEdges] = { 0, 1, 2, 2, 3, 4, 5, 6, 6 };
  size_type targets[numEdges] = { 7, 7, 0, 2, 0, 3, 4, 0, 2 };

#if 0
  // Test case 2.
  const size_type numVerts = 12;
  const size_type numEdges = 11;

  size_type sources[numEdges] = { 0, 1, 2, 3, 5, 6, 6, 7, 7,  9, 10 };
  size_type targets[numEdges] = { 1, 2, 3, 4, 4, 3, 9, 6, 8, 10, 11 };
#endif

#if 0
  // Test case 3.
  const size_type numVerts = 13;
  const size_type numEdges = 13;

  size_type sources[numEdges] = { 0, 1, 2, 3, 5, 6, 6, 7,  7, 8, 8, 10, 11 };
  size_type targets[numEdges] = { 1, 2, 3, 4, 4, 3, 6, 6, 10, 7, 9, 11, 12 };
#endif

#if 0
  // Test case 4.
  const size_type numVerts = 15;
  const size_type numEdges = 15;

  size_type sources[numEdges] = { 0, 1, 2, 3, 5, 6, 6, 6, 7, 8, 10, 11, 12, 13, 14 };
  size_type targets[numEdges] = { 1, 2, 3, 4, 4, 3, 6, 7, 8, 9,  9,  7, 11, 12, 13 };
#endif

#if 0
  // Test case 5.
  const size_type numVerts = 8;
  const size_type numEdges = 15;

  size_type sources[numEdges] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 4, 5, 5, 6, 6 };
  size_type targets[numEdges] = { 2, 7, 0, 1, 7, 0, 2, 2, 0, 4, 0, 3, 4, 2, 0 };
#endif

#if 0
  const size_type numVerts = 8;
  const size_type numEdges = 15;

  size_type sources[numEdges] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6 };
  size_type targets[numEdges] = { 2, 7, 0, 1, 7, 0, 2, 2, 0, 3, 4, 0, 4, 2, 0 };
#endif

#if 0
  const size_type numVerts = 8;
  const size_type numEdges = 14;

  size_type sources[numEdges] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 4, 5, 6, 6 };
  size_type targets[numEdges] = { 2, 7, 0, 1, 7, 0, 2, 2, 0, 0, 3, 4, 2, 0 };
#endif

#if 0
  const size_type numVerts = 8;
  const size_type numEdges = 7;

  size_type sources[numEdges] = { 0, 1, 2, 3, 4, 5, 6 };
  size_type targets[numEdges] = { 1, 2, 3, 4, 5, 6, 7 };
#endif

#if 0
  const size_type numVerts = 8;
  const size_type numEdges = 8;

  size_type sources[numEdges] = { 0, 1, 2, 3, 4, 4, 5, 6 };
  size_type targets[numEdges] = { 1, 2, 3, 4, 3, 5, 6, 7 };
#endif

  init(numVerts, numEdges, sources, targets, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "Graph: (" << order << ", " << size << ")" << std::endl;
  print(g);

  // Create duplicate graph.
  duplicate_adapter<Graph> dg(g);

  std::cout << std::endl << "Duplicate Graph:" << std::endl;
  dg.print();

  std::cout << std::endl
            << "degree(0,  g): " << out_degree(vertices(g)[0], g) << std::endl
            << "degree(0,  dg): " << out_degree(vertices(dg)[0], dg)
            << std::endl;

  std::cout << std::endl << "Testing copy constructor." << std::endl;
  duplicate_adapter<Graph> dg2(dg);
  dg2.print();

  std::cout << std::endl << "Testing assignment." << std::endl;
  duplicate_adapter<Graph> dg3(g);
  dg3 = dg;
  dg3.print();

  return 0;
}
