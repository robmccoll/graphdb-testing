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

//#define HASH_TABLE_COLOR
//#define DEBUG_COUNT_VERTS

#include <mtgl/util.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/mtgl_test.hpp>

#define DEBUG

#include <mtgl/breadth_first_search.hpp>

#ifdef HASH_TABLE_COLOR
#include <mtgl/xmt_hash_table.hpp>
#else
#include <mtgl/dynamic_array.hpp>
#endif

using namespace mtgl;

template <typename Graph>
class counting_bfs_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  counting_bfs_visitor(size_type& nv) : num_visited(nv) {}

  void tree_edge(edge_descriptor& e, Graph &g) { mt_incr(num_visited, 1); }

private:
  size_type& num_visited;
};

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
//  typedef adjacency_list<undirectedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;

  init_test(argc, argv);

  Graph ga;
  create_test_graph(ga, argc, argv);

#ifdef DEBUG
  printf("\n");
#endif

  size_type order = num_vertices(ga);
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, ga);
  vertex_iterator verts = vertices(ga);

#ifdef DEBUG_COUNT_VERTS
  size_type num_visited = 0;
  counting_bfs_visitor<Graph> bfsv(num_visited);
#else
  default_bfs_visitor<Graph> bfsv;
#endif

  const size_type level_size = 2;

  // This is both the set of vertices to visit and the Q space for
  // their neighbors.
  vertex_descriptor* to_visit =
    (vertex_descriptor*) malloc(sizeof(vertex_descriptor) * order);

#ifdef HASH_TABLE_COLOR
  xmt_hash_table<size_type, size_type> htcolors(2 * order);
  typedef map_property_map<size_type, vertex_id_map<Graph> > colormap;
  colormap color(htcolors, vid_map);
#else
  size_type* array_colors = (size_type*) malloc(sizeof(size_type) * order);
  typedef array_property_map<size_type, vertex_id_map<Graph> > colormap;
  colormap color(array_colors, vid_map);
#endif

  #pragma mta assert parallel
  #pragma mta assert nodep
  for (size_type i = 0; i < order; ++i) color[verts[i]] = 0;

  #pragma mta assert parallel
  #pragma mta assert nodep
  for (size_t i = 0; i < level_size; ++i)
  {
    vertex_descriptor v = verts[i];
    to_visit[i] = v;
    put(color, v, 1);
  }

  size_type head = 0;
  size_type tail = level_size;

  expand_one_edge(ga, to_visit, head, tail, bfsv, color);

  // Head and tail now delimit the next level.

  printf("next level = (%lu, %lu)\n", head, tail);
#ifdef DEBUG_COUNT_VERTS
  printf("num_visited: %lu\n", num_visited);
#endif
  printf("\n");

  free(to_visit);
#ifndef HASH_TABLE_COLOR
  free(array_colors);
#endif

  return 0;
}
