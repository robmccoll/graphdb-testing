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
/*! \file test_badrank.cpp

    \author Greg Mackey (gemacke@sandia.gov)

    \date 11/22/2011
*/
/****************************************************************************/

#define PUB_ALG
//#define TEST_STINGER

#include <iostream>
#include <iomanip>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/badrank.hpp>
#include <mtgl/algorithm.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef compressed_sparse_row_graph<directedS> Graph;
#endif
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::thread_vertex_iterator thread_vertex_iterator;

  const int num_ts_args = 2;
  const char* ts_arg_names[num_ts_args] = { "num_bad", "delta" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "The number of bad.",
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

  double delta = static_cast<double>(atof(ts_argv[1]));
  size_type num_bad = atoi(ts_argv[0]);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "Order: " << order << "    Size: " << size << std::endl
            << std::endl;

  if (num_bad > order)
  {
    std::cerr << "Badness greater than order.  Exiting." << std::endl;
    exit(-1);
  }

  if (order < 20) print(g);

  double* ranks = new double[order];
  size_type num_not_bad = order - num_bad;

  vertex_descriptor* bad_verts = new vertex_descriptor[num_bad];

  thread_vertex_iterator tverts = thread_vertices(num_not_bad, g);
  for (size_type i = 0; i < num_bad; ++i, ++tverts)
  {
    bad_verts[i] = *tverts;
  }

  mt_timer timer;
  int issues, memrefs, concur, streams, traps;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  badrank(g, ranks, bad_verts, num_bad, delta);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);
  std::cout << "---------------------------------------------" << std::endl
            << "Bad rank performance stats:" << std::endl
            << "---------------------------------------------" << std::endl;
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

  std::cout << std::endl
            << "---------------------------------------------" << std::endl
            << "Bad Guy Indegree:" << std::endl
            << "---------------------------------------------" << std::endl;

  vertex_property_map<Graph, size_type> in_degrees(g);
  detail::compute_in_degrees(g, in_degrees);

  tverts = thread_vertices(num_not_bad, g);
  for (size_type i = num_not_bad; i < order; ++i, ++tverts)
  {
    std::cout << in_degrees[*tverts] << std::endl;
  }

  // Print at most 20 top bad ranks.
  size_type maxprint = (order > 20 ? 20 : order);

  sort(ranks, num_not_bad, std::greater<double>());

  std::cout << std::endl
            << "---------------------------------------------" << std::endl
            << "Top " << maxprint << " bad ranks:" << std::endl
            << "---------------------------------------------" << std::endl;

  tverts = thread_vertices(0, g);
  for (size_type i = 0; i < maxprint; ++i, ++tverts)
  {
    std::cout << std::fixed << std::setw(8) << std::setprecision(6)
              << ranks[*tverts] << std::endl;
  }

  delete [] ranks;
  delete [] bad_verts;

  return 0;
}
