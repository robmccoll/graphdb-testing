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
/*! \file custom2file.cpp

    \brief Writes out a custom graph file format as described below.

    \author Jon Berry (jberry@sandia.gov)

    \date 3/19/2009

    Let the user define a custom graph in a file of 

    src dest
    src dest
    ...

    then snapshot that out to be read in as the other files are.
*/
/****************************************************************************/

#include <cstdlib>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_io.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <outfile> < <infile>" << std::endl
              << "    where <infile> is a file formatted as follows with "
              << "0-indexing:" << std::endl << std::endl
              << "      <src> <dest>" << std::endl
              << "      <src> <dest>" << std::endl
              << "      ..... ......" << std::endl << std::endl
              << "    with no more than 100 edges." << std::endl;
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

  unsigned long i = 0;
  unsigned long n = 100;
  unsigned long m = 100;
  unsigned long* srcs = (unsigned long*) malloc(m * sizeof(unsigned long));
  unsigned long* dests = (unsigned long*) malloc(m * sizeof(unsigned long));

  while (std::cin >> srcs[i] >> dests[i]) ++i;

  m = i;

  unsigned long maxv = 0;
  for (unsigned long i = 0; i < m; ++i)
  {
    if (srcs[i] > maxv) maxv = srcs[i];
    if (dests[i] > maxv) maxv = dests[i];
  }

  n = maxv + 1;

  Graph g;
  init(n, m, srcs, dests, g);

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

  free(srcs);
  free(dests);

  return 0;
}
