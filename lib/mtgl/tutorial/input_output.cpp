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

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_io.hpp>

using namespace mtgl;

template <typename Graph>
void print_my_graph(Graph& g)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  vertex_iterator verts = vertices(g);
  size_type order = num_vertices(g);
  for (size_type i = 0; i < order; i++)
  {
    vertex_descriptor u = verts[i];
    size_type uid = get(vid_map, u);
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    for (size_type j = 0; j < end; j++)
    {
      vertex_descriptor v = adjs[j];
      size_type vid = get(vid_map, v);
      printf("(%lu, %lu)\n", uid, vid);
    }
  }
}

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  // Usage error message.
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <fileroot>" << std::endl << std::endl
              << "    fileroot: filename root for output file" << std::endl;
    exit(1);
  }

  size_type n;
  size_type m;

  // Read in the number of vertices and edges.
  std::cin >> n;
  std::cin >> m;

  size_type* srcs = new size_type[m];
  size_type* dests = new size_type[m];

  // Read in the ids of each edge's vertices.
  for (size_type i = 0; i < m; ++i)
  {
    std::cin >> srcs[i] >> dests[i];
  }

  // Initialize the graph.
  Graph g;
  init(n, m, srcs, dests, g);

  delete [] srcs;
  delete [] dests;

  // Get the output filenames from the file root supplied as a program
  // argument.
  char srcs_fname[256];
  char dests_fname[256];

  strcpy(srcs_fname, argv[1]);
  strcpy(dests_fname, argv[1]);
  strcat(srcs_fname, ".srcs");
  strcat(dests_fname, ".dests");

  // Write the graph to disk.
  write_binary(g, srcs_fname, dests_fname);

  // Restore the graph from disk to a different graph.
  Graph dg;
  read_binary(dg, srcs_fname, dests_fname);

  // Print the graphs.
  printf("Graph (g):\n");
  print_my_graph(g);
  printf("\n");

  printf("Graph (dg):\n");
  print_my_graph(dg);
  printf("\n");

  return 0;
}
