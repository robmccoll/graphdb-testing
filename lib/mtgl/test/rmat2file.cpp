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
/*! \file rmat2file.cpp

    \author Jon Berry (jberry@sandia.gov)

    \date 12/20/2008
*/
/****************************************************************************/

#include <iostream>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/generate_rmat_graph.hpp>
#include <mtgl/mtgl_io.hpp>
#include <mtgl/random.hpp>

#ifdef USING_QTHREADS
#include <qthread.h>
#endif

using namespace mtgl;

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;

  if (argc < 8)
  {
    std::cerr << "Usage: " << argv[0] << " <outfile> <scale> <deg> <a> <b> "
              << "<c> <d>" << std::endl;
    exit(1);
  }

  char fname[256];
  strcpy(fname, argv[1]);

  int llen = strlen(fname);

  if (!(llen > 3 && strcmp(&fname[llen-3], "mtx") == 0) &&
      !(llen > 6 && strcmp(&fname[llen-6], "dimacs") == 0) &&
      !(llen > 4 && strcmp(&fname[llen-4], "srcs") == 0))
  {
    std::cerr << "Invalid graph output file type.  Available types:\n"
              << "    dimacs (.dimacs)\n"
              << "    matrix market (.mtx)\n"
              << "    srcs-dests (.srcs, .dests)" << std::endl;
    return 1;
  }

  int scale = atoi(argv[2]);            // 2^scale vertices
  int deg = atoi(argv[3]);              // ave degree <deg>
  double a = atof(argv[4]);
  double b = atof(argv[5]);
  double c = atof(argv[6]);
  double d = atof(argv[7]);

  mt_srand48(0);

#ifdef USING_QTHREADS
  qthread_initialize();
#endif

  Graph g;
  generate_rmat_graph(g, scale, deg, a, b, c, d);

  std::cout << "Generated graph: (" << num_vertices(g) << ", "
            << num_edges(g) << ")" << std::endl
            << "Writing graph to: " << fname << std::endl;

  if (llen > 3 && strcmp(&fname[llen-3], "mtx") == 0)
  {
    write_matrix_market(g, fname);
  }
  else if (llen > 6 && strcmp(&fname[llen-6], "dimacs") == 0)
  {
    write_dimacs(g, fname);
  }
  else if (llen > 4 && strcmp(&fname[llen-4], "srcs") == 0)
  {
    char* srcs_fname = fname;
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    write_binary(g, srcs_fname, dests_fname);
  }

  return 0;
}
