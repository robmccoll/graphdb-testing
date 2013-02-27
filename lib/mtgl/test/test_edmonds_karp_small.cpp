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
/*! \file test_edmonds_karp_small.cpp

    \brief This is a driver to test the Edmonds-Karp maxflow code using
           simple, small graphs.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 8/3/2011
*/
/****************************************************************************/

#define DEBUG

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/edmonds_karp_max_flow.hpp>

using namespace mtgl;

#define TEST_CASE 5

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;

#if TEST_CASE == 1
  // Case 1.
  const size_type n_verts = 6;
  const size_type n_edges = 10;

  size_type srcs[n_edges] =  { 0, 0, 1, 2, 1, 3, 2, 4, 3, 4 };
  size_type dests[n_edges] = { 1, 2, 2, 1, 3, 2, 4, 3, 5, 5 };
  int caps_init[n_edges] = { 16, 13, 10, 4, 12, 9, 14, 7, 20, 4 };

  size_type source_v = 0;
  size_type target_v = 5;
#endif

#if TEST_CASE == 2
  // Case 2.
  const size_type n_verts = 10;
  const size_type n_edges = 11;

  size_type srcs[n_edges] =  { 1, 2, 2, 4, 3, 5, 4, 5, 5, 6, 6 };
  size_type dests[n_edges] = { 3, 3, 4, 3, 5, 4, 6, 6, 7, 8, 9 };
  int caps_init[n_edges] = { 12, 8, 10, 5, 19, 4, 20, 3, 11, 11, 15 };

  size_type source_v = 2;
  size_type target_v = 9;
#endif

#if TEST_CASE == 3
  // Case 3.
  const size_type n_verts = 11;
  const size_type n_edges = 16;

  size_type srcs[n_edges] =  { 0, 0, 1, 2, 2, 4, 3, 5, 4, 5, 5, 6, 6, 7, 8,
                                 9 };
  size_type dests[n_edges] = { 1, 2, 3, 3, 4, 3, 5, 4, 6, 6, 7, 8, 9, 10,
                                 10, 10 };
  int caps_init[n_edges] = { 10, 18, 12, 8, 10, 5, 19, 4, 20, 3, 11, 11, 15,
                                8, 6, 14 };

  size_type source_v = 0;
  size_type target_v = 10;
#endif

#if TEST_CASE == 4
  // Case 4.
  const size_type n_verts = 7;
  const size_type n_edges = 8;

  size_type srcs[n_edges] =  { 0, 0, 0, 1, 2, 3, 4, 5 };
  size_type dests[n_edges] = { 1, 2, 3, 4, 4, 5, 6, 6 };
  int caps_init[n_edges] = { 2, 3, 4, 1, 2, 3, 3, 4 };

  size_type source_v = 0;
  size_type target_v = 6;
#endif

#if TEST_CASE == 5
  // Case 5.
  const size_type n_verts = 7;
  const size_type n_edges = 8;

  size_type srcs[n_edges] =  { 0, 0, 1, 1, 2, 3, 4, 5 };
  size_type dests[n_edges] = { 1, 2, 3, 4, 5, 6, 6, 6 };
  int caps_init[n_edges] = { 3, 4, 1, 2, 3, 2, 3, 4 };

  size_type source_v = 0;
  size_type target_v = 6;
#endif

  // Copy the capacities into a dynamically allocated array because disjoint
  // paths max flow requires it.
  int* caps = (int*) malloc(n_edges * sizeof(int));
  for (size_type i = 0; i < n_edges; ++i) caps[i] = caps_init[i];

  Graph ga;
  init(n_verts, n_edges, srcs, dests, ga);

#ifdef DEBUG
  printf(" ");
#endif
  printf("Graph: (%lu, %lu)\n", n_verts, n_edges);

  vertex_descriptor s = vertices(ga)[source_v];
  vertex_descriptor t = vertices(ga)[target_v];

  int flow = edmonds_karp_max_flow(ga, s, t, caps);

  printf("maxflow: %d\n", flow);

  free(caps);
}
