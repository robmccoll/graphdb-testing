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
/*! \file graph_traits.hpp

    \brief This simple class implements the Boost Graph Library's graph_traits
           idiom.

    \author Jon Berry (jberry@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/4/2007
*/
/****************************************************************************/

#ifndef MTGL_GRAPH_TRAITS_HPP
#define MTGL_GRAPH_TRAITS_HPP

namespace mtgl {

template <typename Graph>
class graph_traits {
public:
  typedef typename Graph::size_type                    size_type;
  typedef typename Graph::vertex_descriptor            vertex_descriptor;
  typedef typename Graph::edge_descriptor              edge_descriptor;
  typedef typename Graph::vertex_iterator              vertex_iterator;
  typedef typename Graph::adjacency_iterator           adjacency_iterator;
  typedef typename Graph::in_adjacency_iterator        in_adjacency_iterator;
  typedef typename Graph::edge_iterator                edge_iterator;
  typedef typename Graph::out_edge_iterator            out_edge_iterator;
  typedef typename Graph::in_edge_iterator             in_edge_iterator;
  typedef typename Graph::thread_vertex_iterator       thread_vertex_iterator;
  typedef typename Graph::thread_adjacency_iterator
          thread_adjacency_iterator;
  typedef typename Graph::thread_in_adjacency_iterator
          thread_in_adjacency_iterator;
  typedef typename Graph::thread_edge_iterator         thread_edge_iterator;
  typedef typename Graph::thread_out_edge_iterator     thread_out_edge_iterator;
  typedef typename Graph::thread_in_edge_iterator      thread_in_edge_iterator;
  typedef typename Graph::directed_category            directed_category;
  typedef typename Graph::iterator_category            iterator_category;
};

struct undirectedS {
  static bool is_directed() { return false; }
  static bool is_bidirectional() { return false; }
  static const unsigned long direction_type = 0;
};

struct directedS {
  static bool is_directed() { return true; }
  static bool is_bidirectional() { return false; }
  static const unsigned long direction_type = 1;
};

struct bidirectionalS {
  static bool is_directed() { return true; }
  static bool is_bidirectional() { return true; }
  static const unsigned long direction_type = 2;
};

struct vector_iterators {
  static const bool has_vector_iterators = true;
  static const bool has_thread_iterators = false;
};

struct thread_iterators {
  static const bool has_vector_iterators = false;
  static const bool has_thread_iterators = true;
};

struct vector_thread_iterators {
  static const bool has_vector_iterators = true;
  static const bool has_thread_iterators = true;
};

}

#endif
