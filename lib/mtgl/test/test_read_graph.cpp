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
/*! \file test_read_graph.cpp

    \brief Tests the graph read functions.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 5/20/2010
*/
/****************************************************************************/

#define DEBUG
//#define MM_REAL
#define MM_INT

#include <iostream>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_io.hpp>
#include <mtgl/dynamic_array.hpp>

using namespace mtgl;

int main(int argc, char *argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <filename>\n"
              << "       where filename is in dimacs, matrix market, or "
              << "srcs-dests format.\n\n"
              << "DIMACS files must end with suffix '.dimacs.'\n"
              << "MatrixMarket files must end with suffix '.mtx.'\n"
              << "Srcs-dests files must end with suffixes '.srcs' and '.dests'."
              << std::endl;
    exit(1);
  }

  Graph g;

  int llen = strlen(argv[1]);

#ifdef MM_REAL
  dynamic_array<double> values;
#elif defined MM_INT
  dynamic_array<int> values;
#endif

  if (strcmp(&argv[1][llen - 6], "dimacs") == 0)
  {
    // DIMACS input.
#ifdef MM_INT
    read_dimacs(g, argv[1], values);
#else
    read_dimacs(g, argv[1]);
#endif
  }
  else if (strcmp(&argv[1][llen - 3], "mtx") == 0)
  {
    // Matrix-market input.
#if defined MM_REAL || defined MM_INT
    read_matrix_market(g, argv[1], values);
#else
    read_matrix_market(g, argv[1]);
#endif
  }
  else if (strcmp(&argv[1][llen - 4], "srcs") == 0)
  {
    // Snapshot input: <file>.srcs, <file>.dests
    char* srcs_fname = argv[1];
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    read_binary(g, srcs_fname, dests_fname);
  }
  else
  {
    std::cerr << "Invalid graph file type.\nAvailable types:\n"
              << "    dimacs (.dimacs)\n"
              << "    matrix market (.mtx)\n"
              << "    srcs-dests (.srcs, .dests)" << std::endl;
    exit(1);
  }

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "Vertices: " << order << "    Edges: " << size << std::endl
            << std::endl;

  if (size < 50)
  {
    edge_iterator edgs = edges(g);

    for (size_type i = 0; i < size; ++i)
    {
      edge_descriptor e = edgs[i];
      vertex_descriptor u = source(e, g);
      vertex_descriptor v = target(e, g);

      size_type eid = get(get(_edge_id_map, g), e);
      size_type uid = get(get(_vertex_id_map, g), u);
      size_type vid = get(get(_vertex_id_map, g), v);

      std::cout << std::setw(6) << eid << ": (" << std::setw(6) << uid << ", "
                << std::setw(6) << vid << ")";

#if defined MM_REAL || defined MM_INT
      if (values.size() == size)
      {
        std::cout << "    " << std::setw(6) << std::setprecision(3)
                  << std::fixed << values[i];
      }
      else if (values.size() == 2 * size)
      {
        std::cout << "    " << std::setw(6) << std::setprecision(3)
                  << std::fixed << values[i] << std::setw(6)
                  << std::setprecision(3) << std::fixed << values[i + size];
      }
#endif

      std::cout << std::endl;
    }
  }

  return 0;
}
