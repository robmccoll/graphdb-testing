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
/*! \file test_subgraph_isomorphism.cpp

    \brief This is a driver to test the subgraph isomorphism code.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/17/2008
*/
/****************************************************************************/

#include <iostream>
#include <iomanip>

#include <mtgl/subgraph_isomorphism.hpp>
#include <mtgl/duplicate_adapter.hpp>
#include <mtgl/random_walk.hpp>
#include <mtgl/filter_graph.hpp>

//#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>

using namespace mtgl;

#define TEST1

template <typename Graph>
class match_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename dynamic_array<size_type>::size_type d_size_type;

  match_visitor(Graph& gg) : g(gg) {}

  void operator()(dynamic_array<vertex_descriptor>& path_verts,
                  dynamic_array<edge_descriptor>& path_edges)
  {
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g); 
    edge_id_map<Graph> eid_map = get(_edge_id_map, g); 

    d_size_type num_path_edges = path_edges.size();

    std::cout << "Path:" << std::endl;
    for (d_size_type i = 0; i < num_path_edges; ++i)
    {
      std::cout << std::fixed << std::setw(6) << get(vid_map, path_verts[i])
                << std::fixed << std::setw(6) << get(eid_map, path_edges[i]);
      if (i % 6 == 5 && i + 1 != num_path_edges) std::cout << std::endl;
    }
    std::cout << std::fixed << std::setw(6)
              << get(vid_map, path_verts[num_path_edges]) << std::endl;
  }

private:
  Graph& g;
};

/***/

int main(int argc, char *argv[])
{
  typedef directedS DIRECTION;
//  typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;

  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge;
  typedef graph_traits<Graph>::size_type size_type;
  typedef array_property_map<int, vertex_id_map<Graph> > vertex_property_map;
  typedef array_property_map<int, edge_id_map<Graph> > edge_property_map;

  // Initialize the big graph.
  Graph g;

#ifdef TEST1
  // Test 1.
  const size_type numVerts = 8;
  const size_type numEdges = 10;

  int vTypes[numVerts] = { 0, 0, 0, 0, 0, 0, 0, 0 };

  size_type sources[numEdges] = { 0, 0, 1, 2, 2, 3, 4, 4, 5, 6 };
  size_type targets[numEdges] = { 1, 2, 3, 3, 5, 5, 5, 6, 7, 7 };
  int eTypes[numEdges] = { 0, 1, 1, 0, 1, 2, 0, 0, 1, 1 };
#else
  // Test 2.
  const size_type numVerts = 6;
  const size_type numEdges = 7;

  int vTypes[numVerts] = { 1, 2, 3, 2, 3, 1 };

  size_type sources[numEdges] = { 1, 2, 1, 2, 3, 3, 4 };
  size_type targets[numEdges] = { 0, 0, 2, 3, 4, 5, 5 };
  int eTypes[numEdges] = { 1, 2, 3, 3, 3, 1, 2 };
#endif

  init(numVerts, numEdges, sources, targets, g);

  // Create the big graph vertex and edge type property maps.
  vertex_id_map<Graph> g_vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> g_eid_map = get(_edge_id_map, g);
  vertex_property_map vTypemap(vTypes, g_vid_map);
  edge_property_map eTypemap(eTypes, g_eid_map);

  std::cout << std::endl << "The graph (" << num_vertices(g) << ", "
            << num_edges(g) << ")" << std::endl;
  print(g);

  std::cout << std::endl << "Vertex types:" << std::endl;
  for (size_type i = 0; i < numVerts; ++i)
  {
    std::cout << i << ": (" << vTypes[i] << ")" << std::endl;
  }

  std::cout << std::endl << "Edge types:" << std::endl;
  for (size_type i = 0; i < numEdges; ++i)
  {
    std::cout << i << ": (" << eTypes[i] << ")" << std::endl;
  }

  // Initialize the target graph.
  Graph target;

#ifdef TEST1
  // Test 1.
  const size_type numVertsTarget = 4;
  const size_type numEdgesTarget = 4;

  int vTypesTarget[numVerts] = { 0, 0, 0, 0 };

  size_type sourcesTarget[numEdgesTarget] = { 0, 0, 1, 2 };
  size_type targetsTarget[numEdgesTarget] = { 1, 2, 3, 3 };
  int eTypesTarget[numEdgesTarget]  = { 0, 1, 1, 0 };
#else
  // Test 2.
  const size_type numVertsTarget = 3;
  const size_type numEdgesTarget = 3;

  int vTypesTarget[numVerts] = { 1, 2, 3 };

  size_type sourcesTarget[numEdgesTarget] = { 1, 2, 1 };
  size_type targetsTarget[numEdgesTarget] = { 0, 0, 2 };
  int eTypesTarget[numEdgesTarget]  = { 1, 2, 3 };
#endif

  init(numVertsTarget, numEdgesTarget, sourcesTarget, targetsTarget, target);

  // Create the target graph vertex and edge type property maps.
  vertex_id_map<Graph> trg_vid_map = get(_vertex_id_map, target);
  edge_id_map<Graph> trg_eid_map = get(_edge_id_map, target);
  vertex_property_map vTypemapTarget(vTypesTarget, trg_vid_map);
  edge_property_map eTypemapTarget(eTypesTarget, trg_eid_map);

  std::cout << std::endl << "The target (" << num_vertices(target) << ", "
            << num_edges(target) << ")" << std::endl;
  print(target);

  std::cout << std::endl << "Vertex types:" << std::endl;
  for (size_type i = 0; i < numVertsTarget; ++i)
  {
    std::cout << i << ": (" << vTypesTarget[i] << ")" << std::endl;
  }

  std::cout << std::endl << "Edge types:" << std::endl;
  for (size_type i = 0; i < numEdgesTarget; ++i)
  {
    std::cout << i << ": (" << eTypesTarget[i] << ")" << std::endl;
  }
  std::cout << std::endl;

  typedef si_default_comparator<Graph, vertex_property_map, edge_property_map>
          si_comparator;

  si_comparator si_comp(vTypemap, eTypemap, vTypemapTarget, eTypemapTarget);

  match_visitor<Graph> mvisit(g);

  subgraph_isomorphism(g, target, si_comp, mvisit);

  return 0;
}
