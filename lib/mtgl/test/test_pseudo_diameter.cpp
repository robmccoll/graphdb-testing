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
/*! \file test_pseudo_diameter.cpp

    \author Vitus Leung (vjleung@sandia.gov)

    \date 6/9/2008
*/
/****************************************************************************/

//#define TEST_STINGER

#include <mtgl/util.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/pseudo_diameter.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

int main(int argc, char *argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef compressed_sparse_row_graph<undirectedS> Graph;
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

  mt_timer timer;
  timer.start();
  size_type pd = pseudo_diameter(g);
  timer.stop();

  std::cout << "Pseudo-diameter: " << pd << std::endl
            << "           Time: " << timer.getElapsedSeconds() << std::endl;
}
