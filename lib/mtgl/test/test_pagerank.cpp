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
/*! \file test_pagerank.cpp

    \author Jon Berry (jberry@sandia.gov)

    \date 12/20/2007
*/
/****************************************************************************/

//#define DEBUG
//#define TEST_STINGER

#include <iostream>
#include <iomanip>

#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/pagerank.hpp>
#include <mtgl/algorithm.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/random.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef compressed_sparse_row_graph<bidirectionalS> Graph;
//  typedef adjacency_list<directedS> Graph;
#endif
  typedef graph_traits<Graph>::size_type size_type;

  mt_srand48(0);

  const int num_ts_args = 1;
  const char* ts_arg_names[num_ts_args] = { "delta" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "The conversion tolerance of the power iteration."
  };
  char** ts_argv;
  int ts_argc;

  init_test(argc, argv, num_ts_args, ts_arg_names,
            ts_arg_descs, ts_argc, ts_argv);

#ifdef TEST_STINGER
  SGraph sg(1);
  Graph g(sg);
#else
  Graph g;
#endif

  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "Order: " << order << "    Size: " << size << std::endl
            << std::endl;

  double* ranks = new double[order];

  pagerank(g, ranks, atof(ts_argv[0]));

  mt_timer timer;
  int issues, memrefs, concur, streams;
  init_mta_counters(timer, issues, memrefs, concur, streams);

  pagerank(g, ranks, atof(ts_argv[0]));

  sample_mta_counters(timer, issues, memrefs, concur, streams);
  std::cout << "---------------------------------------------" << std::endl
            << "Page rank performance stats:" << std::endl
            << "---------------------------------------------" << std::endl;
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams);

  // Print at most 20 top page ranks.
  size_type maxprint = (order > 20 ? 20 : order);

  int* ranksidx = new int[order];
  for (size_type i = 0; i < order; ++i) ranksidx[i] = i;

  sort(ranks, order, ranksidx, std::greater<double>());

  std::cout << std::endl
            << "---------------------------------------------" << std::endl
            << "Top " << maxprint << " page ranks:" << std::endl
            << "---------------------------------------------" << std::endl;

  for (size_type i = 0; i < maxprint; ++i)
  {
    std::cout << std::setw(8) << ranksidx[i] << "    " << std::fixed
              << std::setprecision(6) << ranks[i] << std::endl;
  }

  std::cout << std::endl << "Page Rank Stats:" << std::endl
            << "      Max Value: " << ranks[order - 1] << std::endl
            << "      Min Value: " << ranks[0] << std::endl;

  delete [] ranksidx;
  delete [] ranks;

  return 0;
}
