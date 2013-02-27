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
/*! \file graph2file.cpp

    \brief Generates graphs in a particular format or translates from one
           format to another.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 10/24/2011
*/
/****************************************************************************/

#include <iostream>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

int main(int argc, char *argv[])
{
  typedef directedS DIRECTION;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  const int num_ts_args = 1;
  const char* ts_arg_names[num_ts_args] = { "outfile" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "Name of the output file.  Suffix determines type:\n"
    "  .mtx      Matrix Market\n"
    "  .dimacs   Dimacs\n"
    "  .srcs     Snapshot\n"
  };
  char** ts_argv;
  int ts_argc;

  init_test(argc, argv, num_ts_args, ts_arg_names,
            ts_arg_descs, ts_argc, ts_argv);

  int llen = strlen(ts_argv[0]);

  if (!(llen > 3 && strcmp(&ts_argv[0][llen - 3], "mtx") == 0) &&
      !(llen > 6 && strcmp(&ts_argv[0][llen - 6], "dimacs") == 0) &&
      !(llen > 4 && strcmp(&ts_argv[0][llen - 4], "srcs") == 0))
  {
    std::cerr << "Invalid graph output file type.  Available types:\n"
              << "    dimacs (.dimacs)\n"
              << "    matrix market (.mtx)\n"
              << "    srcs-dests (.srcs, .dests)" << std::endl;
    return 1;
  }

  Graph g;

  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  std::cout << "Created graph: (" << order << ", " << size << ")" << std::endl
            << "Writing graph to: " << ts_argv[0] << std::endl;

  if (llen > 3 && strcmp(&ts_argv[0][llen - 3], "mtx") == 0)
  {
    write_matrix_market(g, ts_argv[0]);
  }
  else if (llen > 6 && strcmp(&ts_argv[0][llen - 6], "dimacs") == 0)
  {
    write_dimacs(g, ts_argv[0]);
  }
  else if (llen > 4 && strcmp(&ts_argv[0][llen - 4], "srcs") == 0)
  {
    char* srcs_fname = ts_argv[0];
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    write_binary(g, srcs_fname, dests_fname);
  }

  return 0;
}
