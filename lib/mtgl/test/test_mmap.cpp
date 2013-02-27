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
/*! \file test_mmap.cpp

    \brief Tests that writing and reading a graph via an mmapped memory
           region works by intitializing a graph from an mmaped memory region
           and printing it.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 11/12/2010
*/
/****************************************************************************/

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/mtgl_io.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <mmap filename> <csru|csrd|csrb"
              << "|sting|vec|ht>" << std::endl << std::endl;
    exit(1);
  }

  if (strcmp(argv[2], "csru") == 0)
  {
    compressed_sparse_row_graph<undirectedS> ga;

    // Intialize with mmapped graph.
    if (!read_mmap(ga, argv[1]))
    {
      std::cerr << "Error reading file or mapping memory.  Exiting."
                << std::endl;
      exit(1);
    }

    std::cout << "Graph: (" << num_vertices(ga) << ", " << num_edges(ga)
              << ")" << std::endl;
    print(ga);
  }
  else if (strcmp(argv[2], "csrd") == 0)
  {
    compressed_sparse_row_graph<directedS> ga;

    // Intialize with mmapped graph.
    if (!read_mmap(ga, argv[1]))
    {
      std::cerr << "Error reading file or mapping memory.  Exiting."
                << std::endl;
      exit(1);
    }

    std::cout << "Graph: (" << num_vertices(ga) << ", " << num_edges(ga)
              << ")" << std::endl;
    print(ga);
  }
  else if (strcmp(argv[2], "csrb") == 0)
  {
    compressed_sparse_row_graph<bidirectionalS> ga;

    // Intialize with mmapped graph.
    if (!read_mmap(ga, argv[1]))
    {
      std::cerr << "Error reading file or mapping memory.  Exiting."
                << std::endl;
      exit(1);
    }

    std::cout << "Graph: (" << num_vertices(ga) << ", " << num_edges(ga)
              << ")" << std::endl;
    print(ga);
  }
  else if (strcmp(argv[2], "sting") == 0)
  {
    typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
    typedef stinger_graph_adapter<SGraph> Graph;

    SGraph sg(1);
    Graph ga(sg);

    // Intialize with mmapped graph.
    if (!read_mmap(ga, argv[1]))
    {
      std::cerr << "Error reading file or mapping memory.  Exiting."
                << std::endl;
      exit(1);
    }

    std::cout << "Graph: (" << num_vertices(ga) << ", " << num_edges(ga)
              << ")" << std::endl;
    print(ga);
  }
  else if (strcmp(argv[2], "vec") == 0)
  {
    long* vec;

    // Initialize with mmapped vector.
    unsigned long size = read_mmap(vec, argv[1]);

    if (size == 0)
    {
      std::cerr << "Error reading file or mapping memory.  Exiting."
                << std::endl;
      exit(1);
    }

    std::cout << "Vector: " << size << std::endl;
    for (unsigned long i = 0; i < size; ++i) std::cout << vec[i] << std::endl;
  }
  else if (strcmp(argv[2], "ht") == 0)
  {
    xmt_hash_table<int, int> ht;

    // Initialize with mmapped hash table.
    if (!read_mmap(ht, argv[1]))
    {
      std::cerr << "Error reading file or mapping memory.  Exiting."
                << std::endl;
      exit(1);
    }

    std::cout << "Hash Table:" << std::endl;
    ht.print();
  }
  else
  {
    std::cerr << "Invalid object type.  Must be csru, csrd, csrb, sting, "
              << "vec, ht." << std::endl;
  }

  return 0;
}
