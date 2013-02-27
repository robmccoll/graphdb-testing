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
/*! \file test_vertex_betweenness.cpp

    \author Vitus Leung (vjleung@sandia.gov)

    \date 7/10/2008
*/
/****************************************************************************/

//#define TEST_STINGER

#include <mtgl/vertex_betweenness.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/util.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

int main(int argc, char *argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef compressed_sparse_row_graph<directedS> Graph;
#endif
  typedef graph_traits<Graph>::size_type size_type;

  init_test(argc, argv);

#ifdef TEST_STINGER
  SGraph sg(1);
  Graph g(sg);
#else
  Graph g;
#endif

  create_test_graph(g, argc, argv);

  std::cout << "Graph: " << num_vertices(g) << ", " << num_edges(g)
            << std::endl << std::endl;

  size_type order = num_vertices(g);
  double* res = new double[order];

  mt_timer timer;
  timer.start();

  vertex_betweenness(g, res);

  timer.stop();

  std::cout << "Vertex betweenness centrality:" << std::endl;

  size_type print_num = order > 40 ? 40 : order;
  for (size_type i = 0; i < print_num; ++i)
  {
    std::cout << i << "\t" << std::fixed << std::setprecision(6) << res[i]
              << std::endl;
  }

  std::cout << std::endl << "Time: " << timer.getElapsedSeconds() << std::endl;

  delete [] res;
}
