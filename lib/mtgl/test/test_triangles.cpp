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
/*! \file test_triangles.cpp

    \author Jon Berry (jberry@sandia.gov)

    \date 1/19/2009
*/
/****************************************************************************/

#include <iostream>

#include <mtgl/util.hpp>
#include <mtgl/triangles.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/random.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

template <typename Graph>
class count_triangles_visitor : public default_triangles_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;

  count_triangles_visitor(int& c) : my_count(0), count(c) {}

  void operator()(size_type e1, size_type e2, size_type e3) { ++my_count; }

  void accumulate() { mt_incr(count, my_count); }

private:
  int my_count;
  int& count;
};

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<undirectedS> Graph;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::size_type size_type;

  mt_srand48(0);

  init_test(argc, argv);

  Graph ga;
  create_test_graph(ga, argc, argv);

//  print(ga);
  size_type order = num_vertices(ga);
  size_type size = num_edges(ga);

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, ga);

  std::cout << "Starting test on " << order << " vertices, " << size
            << " edges" << std::endl;

  int count = 0;
  count_triangles_visitor<Graph> ctv(count);
  default_triangles_visitor<Graph> dtv;

  mt_timer tri_time;
  tri_time.start();
  find_triangles(ga, ctv);
  tri_time.stop();

  std::cout << "RESULT: find_triangles " << count << std::endl;
  std::cout << "tri time: " << tri_time.getElapsedSeconds() << std::endl;

  return 0;
}
