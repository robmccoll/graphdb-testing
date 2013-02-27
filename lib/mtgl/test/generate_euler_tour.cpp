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
/*! \file generate_euler_tour.cpp

    \brief Generates an euler tour through a graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/15/2011

    The function reads in the graph as undirected.  This works only for
    graphs that are connected.
*/
/****************************************************************************/

//#define DEBUG

#include <iostream>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/euler_tour.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;

  const int num_ts_args = 1;
  const char* ts_arg_names[num_ts_args] = { "src" };
  const char* ts_arg_descs[num_ts_args] = { "Starting vertex of euler tour." };
  char** ts_argv;
  int ts_argc;

  init_test(argc, argv, num_ts_args, ts_arg_names,
            ts_arg_descs, ts_argc, ts_argv);

  Graph ga;
  create_test_graph(ga, argc, argv);

  size_type order = num_vertices(ga);

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, ga);
  edge_id_map<Graph> eid_map = get(_edge_id_map, ga);

  if (order == 0)
  {
    std::cerr << "Error reading or creating file.  Exiting." << std::endl;
    exit(1);
  }

  dynamic_array<vertex_descriptor> path_verts;
  dynamic_array<edge_descriptor> path_edges;

  bool success =
    duplicate_euler_tour(ga, atoi(ts_argv[0]), path_verts, path_edges);

  assert(success);

  for (size_type i = 0; i < path_edges.size(); ++i)
  {
    vertex_descriptor v = path_verts[i];
    edge_descriptor e = path_edges[i];

    std::cout << "v: " << get(vid_map, v) << ",  e: {" << get(eid_map, e)
              << "}(" <<  get(vid_map, source(e, ga)) << ", "
              << get(vid_map, target(e, ga)) << ")" << std::endl;
  }

  std::cout << "v: " << get(vid_map, path_verts[path_verts.size()-1])
            << std::endl;

  return 0;
}
