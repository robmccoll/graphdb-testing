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
/*! \file test_connected_components.cpp

    \brief Tests the connected component algorithm.

    \author Jon Berry (jberry@sandia.gov)

    \date 3/3/2009
*/
/****************************************************************************/

#include <iostream>
#include <list>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/connected_components.hpp>
#include <mtgl/util.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/random.hpp>
#include <mtgl/dynamic_array.hpp>

//#define DEBUG

using namespace mtgl;

int main(int argc, char* argv[])
{
  // The connected components algorithm used here is INVALID over
  // non-undirected graphs.
  //typedef compressed_sparse_row_graph<undirectedS> Graph;
  typedef adjacency_list<undirectedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;

  mt_srand48(0);

  init_test(argc, argv);

  Graph ga;
  create_test_graph(ga, argc, argv);

  size_type order = num_vertices(ga);
  size_type size = num_edges(ga);
  if (order < 100)
  {
    print(ga);
  }
  else
  {
    std::cout << "Graph constructed with " << order << " vertices and "
              << size << " edges." << std::endl;
  }

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, ga);
  vertex_iterator verts = vertices(ga);

  vertex_property_map<Graph, size_type> components(ga);

  mt_timer cc_time;

  cc_time.start();
  connected_components(ga, components);
  cc_time.stop();

  std::cout << "gcc_sv time: " << cc_time.getElapsedSeconds() << std::endl;
  std::cout << "There are " << count_connected_components(ga, components)
            << " connected components." << std::endl;

  cc_time.start();
  shiloach_vishkin(ga, components);
  cc_time.stop();

  std::cout << "sv time (going through edges): " << cc_time.getElapsedSeconds()
            << std::endl;
  std::cout << "There are " << count_connected_components(ga, components)
            << " connected components." << std::endl;

#ifdef DEBUG
  vertex_descriptor largest_leader =
    largest_connected_component(ga, components);

  std::cout << "The largest has leader " << get(vid_map, largest_leader)
            << std::endl;

#ifndef __MTA__
  std::list<vertex_descriptor> the_component;

  for (int i = 0; i < order; i++)
  {
    if (components[verts[i]] == largest_leader)
    {
      the_component.push_back(verts[i]);
    }
  }

  the_component.sort();
  the_component.unique();
  std::cout << "comp.order: " << the_component.size() << std::endl;

  for (std::list<vertex_descriptor>::iterator it = the_component.begin();
       it != the_component.end(); ++it)
  {
    std::cout << get(vid_map, *it) << std::endl;
  }
#endif

  dynamic_array<vertex_descriptor> leaders;
  component_leaders(ga, components, leaders);
  std::cout << "There are " << leaders.size() << " leaders." << std::endl;
#if 0
  for (size_type i = 0; i < leaders.size(); ++i)
  {
    std::cout << "\t" << get(vid_map, leaders[i]) << std::endl;
  }
#endif
#endif

  return 0;
}
