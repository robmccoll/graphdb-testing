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
/*! \file test_write_graph.cpp

    \brief Tests the graph write functions.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 1/5/2011
*/
/****************************************************************************/

#define DEBUG
//#define MM_COMPLEX
//#define MM_REAL
//#define MM_INT

#include <iostream>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_io.hpp>
#include <mtgl/dynamic_array.hpp>

using namespace mtgl;

int main(int argc, char *argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <in filename> <out filename>\n"
              << "       where the files are in dimacs, matrix market, or "
              << "srcs-dests format.\n\n"
              << "DIMACS files must end with suffix '.dimacs.'\n"
              << "MatrixMarket files must end with suffix '.mtx.'\n"
              << "Srcs-dests files must end with suffixes '.srcs' and '.dests'."
              << std::endl;
    return 1;
  }

  int llen = strlen(argv[1]);

  if (!(llen > 3 && strcmp(&argv[1][llen - 3], "mtx") == 0) &&
      !(llen > 6 && strcmp(&argv[1][llen - 6], "dimacs") == 0) &&
      !(llen > 4 && strcmp(&argv[1][llen - 4], "srcs") == 0))
  {
    std::cerr << "Invalid graph input file type.  Available types:\n"
              << "    dimacs (.dimacs)\n"
              << "    matrix market (.mtx)\n"
              << "    srcs-dests (.srcs, .dests)" << std::endl;
    return 1;
  }

  llen = strlen(argv[2]);

  if (!(llen > 3 && strcmp(&argv[2][llen - 3], "mtx") == 0) &&
      !(llen > 6 && strcmp(&argv[2][llen - 6], "dimacs") == 0) &&
      !(llen > 4 && strcmp(&argv[2][llen - 4], "srcs") == 0))
  {
    std::cerr << "Invalid graph output file type.  Available types:\n"
              << "    dimacs (.dimacs)\n"
              << "    matrix market (.mtx)\n"
              << "    srcs-dests (.srcs, .dests)" << std::endl;
    return 1;
  }

  Graph g;

  llen = strlen(argv[1]);

#if defined MM_COMPLEX || defined MM_REAL
  dynamic_array<double> values;
#elif defined MM_INT
  dynamic_array<int> values;
#endif

  std::cout << "Reading graph." << std::endl;

  if (llen > 3 && strcmp(&argv[1][llen - 3], "mtx") == 0)
  {
    // Matrix-market input.
#if defined MM_COMPLEX || defined MM_REAL || defined MM_INT
    read_matrix_market(g, argv[1], values);
#else
    read_matrix_market(g, argv[1]);
#endif
  }
  else if (llen > 6 && strcmp(&argv[1][llen - 6], "dimacs") == 0)
  {
    // DIMACS input.
#ifdef MM_INT
    read_dimacs(g, argv[1], values);
#else
    read_dimacs(g, argv[1]);
#endif
  }
  else if (llen > 4 && strcmp(&argv[1][llen - 4], "srcs") == 0)
  {
    // Snapshot input: <file>.srcs, <file>.dests
    char* srcs_fname = argv[1];
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    read_binary(g, srcs_fname, dests_fname);
  }

  std::cout << std::endl << "Vertices: " << num_vertices(g) << "    Edges: "
            << num_edges(g) << std::endl << std::endl
            << "Writing graph:" << std::endl;

  llen = strlen(argv[2]);

  if (llen > 3 && strcmp(&argv[2][llen - 3], "mtx") == 0)
  {
    // Matrix-market input.
#if defined MM_COMPLEX || defined MM_REAL || defined MM_INT
    write_matrix_market(g, argv[2], values);
#else
    write_matrix_market(g, argv[2]);
#endif
  }
  else if (llen > 6 && strcmp(&argv[2][llen - 6], "dimacs") == 0)
  {
    // DIMACS input.
#ifdef MM_INT
    write_dimacs(g, argv[2], values);
#else
    write_dimacs(g, argv[2]);
#endif
  }
  else if (llen > 4 && strcmp(&argv[2][llen - 4], "srcs") == 0)
  {
    // Snapshot input: <file>.srcs, <file>.dests
    char* srcs_fname = argv[2];
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    write_binary(g, srcs_fname, dests_fname);
  }

  return 0;
}
