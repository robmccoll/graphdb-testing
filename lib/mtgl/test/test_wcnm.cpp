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
/*! \file test_wcnm.cpp

    \author Jon Berry (jberry@sandia.gov)

    \date 8/2010
*/
/****************************************************************************/

#include <cstdlib>

#define RECTANGLES
#define MB
//#define WT

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/wcnm.hpp>
#include <mtgl/random.hpp>

// Written weights are int(w*10^7), retrieve w:
//#define READ_SCALED_DIMACS_WEIGHTS 1

using namespace mtgl;

int main(int argc, char* argv[])
{
  // Use the following when the input file lists all edges (both directions).
  //typedef compressed_sparse_row_graph<directedS> Graph;

  // Use the following when the input file lists only directed edges
  // (one direction only).
  typedef compressed_sparse_row_graph<undirectedS> Graph;

  typedef graph_traits<Graph>::size_type size_type;

  mt_srand48(0);

  const int num_ts_args = 2;
  const char* ts_arg_names[num_ts_args] = { "num_communities",
                                            "num_iterations" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "The number of resulting communities.  A value of 0 lets the algorithm "
      "decide the stoppping point.",
    "The number of weighted passes in wcnm."
  };
  char** ts_argv;
  int ts_argc;

  init_test(argc, argv, num_ts_args, ts_arg_names,
            ts_arg_descs, ts_argc, ts_argv);

  Graph g;
  dynamic_array<double> weights;
  create_test_graph(g, weights, argc, argv);

  int num_communities = atoi(ts_argv[0]);
  int num_iterations = atoi(ts_argv[1]);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  std::cout << "Order: " << order << ", Size: " << size << std::endl;

  size_type* leader = (size_type*) malloc(sizeof(size_type) * order);

  double* dweights = 0;

  if (weights.size() == size)
  {
    dweights = (double*) malloc(sizeof(double) * size);

    #pragma mta assert parallel
    for (size_type i = 0; i < size; ++i)
    {
#ifdef READ_SCALED_DIMACS_WEIGHTS
      dweights[i] = weights[i] / 1e7;
#else
//      dweights[i] = weights[i];
      dweights[i] = 1.0;
#endif
    }
  }

  mt_timer timer;
  int issues, memrefs, concur, streams;
  init_mta_counters(timer, issues, memrefs, concur, streams);

  std::cout << std::endl << "desired communities: " << num_communities
            << std::endl;

  wcnm(g, leader, num_communities, num_iterations, dweights);

  std::cout << "  found communities: " << num_communities
            << std::endl;

  if (order < 100)
  {
    std::cout << std::endl;
    for (size_type i = 0; i < order; ++i)
    {
      std::cout << "leader[ " << i << "]: " << leader[i] << std::endl;
    }
  }

  sample_mta_counters(timer, issues, memrefs, concur, streams);
  std::cout << std::endl << "wcnm performance stats:" << std::endl
            << "---------------------------------------------" << std::endl;
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams);

  free(leader);
  if (weights.size() == size) free(dweights);

  return 0;
}
