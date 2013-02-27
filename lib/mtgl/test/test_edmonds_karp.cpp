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
/*! \file test_dpmf.cpp

    \brief Tests the Edmonds-Karp maxflow code.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 8/3/2011

    When running this test in parallel, you must use the snapshot format (or
    an mmap created from the snapshot format) as input to get the same flow
    results between runs.  The ordering of the edges in the sources and
    dests arrays can be different when reading from text files in parallel.
    This test generates edge capacities based on the order of the edges.
    Combine these two facts, and you get different results between runs.
*/
/****************************************************************************/

//#define DEBUG

#include <cstdlib>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/edmonds_karp_max_flow.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/util.hpp>

//#define DOUBLE

using namespace mtgl;

int main(int argc, char *argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  const int num_ts_args = 2;
  const char* ts_arg_names[num_ts_args] = { "src", "dest" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "Source vertex for maxflow.",
    "Target vertex for maxflow."
  };
  char** ts_argv;
  int ts_argc;

  init_test(argc, argv, num_ts_args, ts_arg_names,
            ts_arg_descs, ts_argc, ts_argv);

  Graph ga;
  create_test_graph(ga, argc, argv);

  size_type order = num_vertices(ga);
  size_type size = num_edges(ga);

  size_type sid = atoi(ts_argv[0]);
  size_type tid = atoi(ts_argv[1]);

  if (sid == tid)
  {
    fprintf(stderr, "Error: Source can't equal sink.\n");
    exit(1);
  }

  printf("sid: %lu, tid: %lu\n", sid, tid);
#ifdef DEBUG
  printf(" ");
#endif
  printf("Graph: (%lu, %lu)\n", order, size);

  vertex_descriptor s = vertices(ga)[sid];
  vertex_descriptor t = vertices(ga)[tid];

#ifdef DOUBLE
  typedef double flow_t;
#else
  typedef int flow_t;
#endif

  flow_t* caps = (flow_t*) malloc(size * sizeof(flow_t));
  for (size_type i = 0; i < size; i++) caps[i] = 5 + i % 10;

  mt_timer timer;
  timer.start();

  flow_t flow = edmonds_karp_max_flow(ga, s, t, caps);

  timer.stop();

  std::cout << "Flow: " << flow << std::endl;
  std::cout << "Secs: " << timer.getElapsedSeconds() << std::endl;

  free(caps);

  return 0;
}
