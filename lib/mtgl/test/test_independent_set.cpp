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
/*! \file test_independent_set.cpp

    \author Jon Berry (jberry@sandia.gov)

    \date 4/26/2009
*/
/****************************************************************************/

#include <iostream>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/independent_set.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/random.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  mt_srand48(0);

  init_test(argc, argv);

  Graph g;
  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "ORDER: " << order << ", SIZE: " << size << std::endl;
//  print(g);

  bool* active = new bool[order];

  size_type ind_set_size;

  mt_timer timer;

  std::cout << std::endl
            << "---------------------------------" << std::endl
            << "MIS" << std::endl
            << "---------------------------------" << std::endl;

  timer.start();

  ind_set_size = maximal_independent_set(g, active);

  timer.stop();
  std::cout << "Time: " << timer.getElapsedSeconds() << std::endl
            << "Size: " << ind_set_size << std::endl;

//  print_independent_set(active, order);
  validate_maximal_independent_set(g, active);

  delete [] active;

  return 0;
}
