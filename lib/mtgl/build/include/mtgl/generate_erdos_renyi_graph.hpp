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
/*! \file generate_erdos_renyi_graph.hpp

    \brief Generates an Erdos-Renyi graph.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 4/22/2011
*/
/****************************************************************************/

#ifndef MTGL_GENERATE_ERDOS_RENYI_GRAPH_HPP
#define MTGL_GENERATE_ERDOS_RENYI_GRAPH_HPP

#include <iostream>

#include <mtgl/util.hpp>
#include <mtgl/random.hpp>
#include <mtgl/mtgl_adapter.hpp>

namespace mtgl {

/****************************************************************************/
/*! \brief Generates an Erdos-Renyi graph with n vertices where p represents
           the probability an edge will exist between any 2 given vertices.
*/
/****************************************************************************/
template <typename Graph>
void
generate_erdos_renyi_graph(Graph& g,
                           typename graph_traits<Graph>::size_type n, double p)
{
  typedef typename graph_traits<Graph>::size_type size_type;

#ifdef DEBUG
  std::cout << "Generating ER graph." << std::endl;
#endif

  drand48_generator randVals(n * n);

  size_type max_e = (n * (n-1) / 2) + 1;
  assert(max_e < 1e10);

  size_type* srcs = (size_type *) malloc(sizeof(size_type) * max_e);
  size_type* dests = (size_type *) malloc(sizeof(size_type) * max_e);

  size_type m = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < n; ++i)
  {
    #pragma mta assert nodep
    for (size_type j = i + 1; j < n; ++j)
    {
      double my_p = randVals[i * n + j];

      if (my_p < p)
      {
        size_type next = mt_incr(m, 1);
        srcs[next] = i;
        dests[next] = j;
      }
    }
  }

#ifdef PERMUTE_NODES
  std::cout << "Permuting nodes." << std::endl;
  random_permutation(n, m, srcs, dests);
#endif

  init(n, m, srcs, dests, g);

  free(srcs);
  free(dests);
}

}

#endif
