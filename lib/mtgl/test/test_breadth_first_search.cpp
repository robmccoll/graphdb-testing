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
/*! \file test_breadth_first_search.cpp

    \brief Tests the breadth_first_search code.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 9/15/2010
*/
/****************************************************************************/

//#define DEBUG
//#define PHASE_DEBUG
//#define TEST_STINGER

#include <cstdlib>
#include <iostream>

#include <mtgl/util.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/breadth_first_search.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef compressed_sparse_row_graph<directedS> Graph;
//  typedef adjacency_list<directedS> Graph;
#endif
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::thread_vertex_iterator thread_vertex_iterator;

  const int num_ts_args = 1;
  const char* ts_arg_names[num_ts_args] = { "src" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "Index of the bfs source.  A value of 'm' causes the vertex of "
      "maximum degree to be the source."
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

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  std::cout << "Graph: (" << order << ", " << size << ")" << std::endl;
//  print(g);

  size_type source_id;

  if (ts_argv[0][0] == 'm')
  {
    size_type* degrees = (size_type*) malloc(order * sizeof(size_type));

    size_type stream_id = 0;
    size_type num_streams = 1;

    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++verts)
      {
        degrees[start_pos] = out_degree(*verts, g);
      }
    }

    source_id = 0;
    size_type maxdeg = degrees[0];

    #pragma mta assert nodep
    for (size_type i = 1; i < order; ++i)
    {
      if (degrees[i] > maxdeg)
      {
        maxdeg = degrees[i];
        source_id = i;
      }
    }

    std::cout << "index " << source_id << " max deg " << maxdeg << std::endl;

    free(degrees);
  }
  else
  {
    source_id = atoi(ts_argv[0]);
  }
  std::cout << std::endl;

  vertex_descriptor null_vert = null_vertex(g);
  vertex_descriptor source = *thread_vertices(source_id, g);

  typedef vertex_property_map<Graph, vertex_descriptor> VertexPredMap;
  VertexPredMap bfs_parents(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      bfs_parents[*verts] = null_vert;
    }
  }

  bfs_parents[source] = source;

  mt_timer timer;
  int issues, memrefs, concur, streams, traps;

  parents_bfs_visitor<Graph, VertexPredMap> pbv(bfs_parents);
//  default_bfs_visitor<Graph> pbv;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  breadth_first_search(g, source, pbv);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  // On the XMT you get slightly better results the second time when running
  // an algorithm twice consecutively.  So, we do that here.

  std::cout << std::endl;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      bfs_parents[*verts] = null_vert;
    }
  }

  bfs_parents[source] = source;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  breadth_first_search(g, source, pbv);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

  std::cout << std::endl;
#endif

  tree_check(g, source, bfs_parents);

  return 0;
}
