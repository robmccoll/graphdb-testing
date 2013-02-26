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
/*! \file partitioning.hpp

    \brief Functions for partitioning work into equal-sized chunks.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric L. Goodman (elgoodm@sandia.gov)

    \date 6/23/2011

    These partitioning functions are useful when you want to iterate over a
    space using explicit thread-based iteration.  In our experience the best
    way to do this is to create many more chunks of works than streams.  This
    has consistently performed better for us than creating a single chunk of
    work for each stream even if you evenly divide the work between streams.
    In addition, we have observed that in many circumstances a block size of
    1024 seems to be the right size for achieving the best performance.
    However, this size may need to be adjusted based on the amount of work
    per iteration.

    BLOCK PARTITIONING

    Block patitioning divides a 1D space into equal-size chunks.  It is
    performed by the following functions:
      begin_block_range()
      end_block_range()

    Block partitioning is useful when you want to perform thread iteration
    over a 1D space.  Here is an example of using block partitioning.

      // Set a block size and calculate the number of blocks.
      unsigned long BLOCK_SIZE = 1024;
      unsigned long num_blocks = (total_work + BLOCK_SIZE - 1) / BLOCK_SIZE;

      // Iterate over the blocks.
      for (unsigned long block_id = 0; block_id < num_blocks; ++block_id)
      {
        // Find this block's portion of the work.
        unsigned long start_pos = begin_block_range(total_work,
                                                    block_id, num_blocks);
        unsigned long end_pos = end_block_range(total_work,
                                                block_id, num_blocks);

        // Iterate over this blocks portion of the work.
        for (unsigned long i = start_pos; i != end_pos; ++i)
        {
          // Do work.
        }
      }

    MANHATTAN PARTITIONING

    Manhattan partitioning divides a 2D space, where the inner dimension has
    differing sizes, into equal sized chunks.  It is performed by the
    following functions:
      begin_manhattan_outer_range()
      end_manhattan_outer_range()
      begin_manhattan_inner_range()
      end_manhattan_inner_range()

    Manhattan partitioning is useful when you want to perform thread iteration
    over a 2D space where the inner dimension has differning sizes.  Here is
    an example of using manhattan partitioning.

      // Create the accumulation array of work for the inner dimensions.
      unsigned long* accum_work =
        (unsigned long*) malloc((num_outer + 1) * sizeof(unsigned long));
      accum_work[0] = 0;
      for (unsigned long i = 0; i < num_outer; ++i)
      {
        accum_work[i + 1] = accum_work[i] + work[i];
      }

      // Set a block size and calculate the number of blocks.
      unsigned long BLOCK_SIZE = 1024;
      unsigned long num_blocks = (accum_work[num_outer] + BLOCK_SIZE - 1) /
                                 BLOCK_SIZE;

      // Iterate over the blocks.
      for (unsigned long block_id = 0; block_id < num_blocks; ++block_id)
      {
        // Find this block's portion of the total work.
        unsigned long start_pos = begin_block_range(accum_work[num_outer],
                                                    block_id, num_blocks);
        unsigned long end_pos = end_block_range(accum_work[num_outer],
                                                block_id, num_blocks);

        // Find this block's start and past-the-end indices for the outer
        // dimension.
        size_type start_outer =
          begin_manhattan_outer_range(accum_work, num_outer, start_pos);
        size_type end_outer =
          end_manhattan_outer_range(accum_work, num_outer, end_pos);

        // Iterate over the outer dimension.
        for (unsigned long i = start_outer; i != end_outer; ++i)
        {
          // Find this block's start and past-the-end indices for the inner
          // dimension of the i'th outer dimension.
          size_type start_inner =
            begin_manhattan_inner_range(accum_work, start_pos, start_outer, i);
          size_type end_inner =
            end_manhattan_inner_range(accum_work, end_pos, end_outer, i);

          // Iterate over the inner dimension.
          for (size_type j = start_inner; j != end_inner; ++j)
          {
            // Do work.
          }
        }
      }

    In some cases you may need to add a mta no scalar expansion pragma on the
    loop iterating over blocks to prevent the compiler from generating scalar
    expansions of start_pos and end_pos.  In our experience, this makes the
    iteration perform worse.
*/
/****************************************************************************/

#ifndef MTGL_PARTITIONING_HPP
#define MTGL_PARTITIONING_HPP

#include <cstddef>

#include <mtgl/graph_traits.hpp>

namespace mtgl {

/*! \brief Returns the beginning index for the range of work assigned to a
           block.

    \param num_elements Total number of elements that need to be paritioned.
    \param block_id The id of the current block.
    \param num_blocks Total number of blocks.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric L. Goodman (elgoodm@sandia.gov)

    \date 7/9/2010
*/
template <typename T>
inline
T begin_block_range(T num_elements, T block_id, T num_blocks)
{
  return static_cast<T>((static_cast<double>(num_elements) / num_blocks) *
                        block_id);
}

/*! \brief Returns the past-the-end index for the range of work assigned to a
           block.

    \param num_elements Total number of elements that need to be paritioned.
    \param block_id The id of the current block.
    \param num_blocks Total number of blocks.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric L. Goodman (elgoodm@sandia.gov)

    \date 7/9/2010
*/
template <typename T>
inline
T end_block_range(T num_elements, T block_id, T num_blocks)
{
  return (block_id + 1 < num_blocks) ?
         static_cast<T>((static_cast<double>(num_elements) / num_blocks) *
                        (block_id + 1)) :
         num_elements;
}

/*! \brief Returns the beginning outer index corresponding with start_pos.

    \param accum An accumulation array of the number of inner iterations for
                 each outer iteration.  Entry i represents the total number
                 of inner iteration for outer iterations 0 through i-1.
    \param accum_size The size of the accum array.
    \param start_pos The starting position for a chunk of work in the range
                     [0, accum[accum_size]).

    \author Greg Mackey (gemacke@sandia.gov)

    \date 6/23/2011
*/
template <typename T>
inline
T begin_manhattan_outer_range(T* accum, T accum_size, T start_pos)
{
  // Binary search accum to find the position in the list that
  // contains the first value > start_pos.  Note that we need to search
  // the range [1, accum_size + 1).
  T first = 1;
  T count = accum_size;
  while (count > 0)
  {
    T step = count / 2;
    T mid = first + step;
    if (!(start_pos < accum[mid]))
    {
      first = mid + 1;
      count -= step + 1;
    }
    else
    {
      count = step;
    }
  }

  return start_pos == 0 ? 0 : first - 1;
}

/*! \brief Returns the past-the-end outer index corresponding with end_pos.

    \param accum An accumulation array of the number of inner iterations for
                 each outer iteration.  Entry i represents the total number
                 of inner iteration for outer iterations 0 through i-1.
    \param accum_size The size of the accum array.
    \param end_pos The ending position for a chunk of work in the range
                   [0, accum[accum_size]).

    \author Greg Mackey (gemacke@sandia.gov)

    \date 6/23/2011
*/
template <typename T>
inline
T end_manhattan_outer_range(T* accum, T accum_size, T end_pos)
{
  // Binary search accum to find the position in the list that
  // contains the first value > end_pos.  Note that we need to search
  // the range [1, accum_size + 1).
  T first = 1;
  T count = accum_size;
  while (count > 0)
  {
    T step = count / 2;
    T mid = first + step;
    if (!(end_pos < accum[mid]))
    {
      first = mid + 1;
      count -= step + 1;
    }
    else
    {
      count = step;
    }
  }

  return accum[first - 1] == end_pos ? first - 1 : first;
}

/*! \brief Returns the beginning inner index for outer iteration i in the
           range [start_outer, end_outer].

    \param accum An accumulation array of the number of inner iterations for
                 each outer iteration.  Entry i represents the total number
                 of inner iteration for outer iterations 0 through i-1.
    \param start_pos The starting position for a chunk of work in the range
                     [0, accum[accum_size]).
    \param start_outer The starting outer iteration corresponding to the chunk
                       of work defined by [start_pos, end_pos).
    \param i The current outer iteration.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 6/23/2011
*/
template <typename T>
inline
T begin_manhattan_inner_range(T* accum, T start_pos, T start_outer, T i)
{
  return i == start_outer ? start_pos - accum[i] : 0;
}

/*! \brief Returns the past-the-end inner iteration for outer iteration i in
           the range [start_outer, end_outer].

    \param accum An accumulation array of the number of inner iterations for
                 each outer iteration.  Entry i represents the total number
                 of inner iteration for outer iterations 0 through i-1.
    \param end_pos The ending position for a chunk of work in the range
                   [0, accum[accum_size]).
    \param end_outer The ending outer iteration corresponding to the chunk
                     of work defined by [end_pos, end_pos).
    \param i The current outer iteration.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 6/23/2011
*/
template <typename T>
inline
T end_manhattan_inner_range(T* accum, T end_pos, T end_outer, T i)
{
  return i + 1 == end_outer ? end_pos - accum[i] : accum[i + 1] - accum[i];
}

#ifdef USING_QTHREADS
namespace detail {

template <typename Graph, typename T>
class qt_accum_out_deg {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  qt_accum_out_deg(Graph& gg, T* ad) : g(gg), accum_deg(ad) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      accum_deg[start_pos + 1] = out_degree(*tverts, g);
    }
  }

private:
  Graph& g;
  T* accum_deg;
};

template <typename Graph, typename T>
class qt_accum_out_deg_list {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  qt_accum_out_deg_list(Graph& gg, T* ad, vertex_descriptor* vl) :
    g(gg), accum_deg(ad), vlist(vl) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i < stop; ++i)
    {
      accum_deg[i + 1] = out_degree(vlist[i], g);
    }
  }

private:
  Graph& g;
  T* accum_deg;
  vertex_descriptor* vlist;
};

template <typename Graph, typename T>
class qt_accum_in_deg {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  qt_accum_in_deg(Graph& gg, T* ad) : g(gg), accum_deg(ad) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type start_pos = start;
    size_type end_pos = stop;

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      accum_deg[start_pos + 1] = in_degree(*tverts, g);
    }
  }

private:
  Graph& g;
  T* accum_deg;
};

template <typename Graph, typename T>
class qt_accum_in_deg_list {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  qt_accum_in_deg_list(Graph& gg, T* ad, vertex_descriptor* vl) :
    g(gg), accum_deg(ad), vlist(vl) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i < stop; ++i)
    {
      accum_deg[i + 1] = in_degree(vlist[i], g);
    }
  }

private:
  Graph& g;
  T* accum_deg;
  vertex_descriptor* vlist;
};

}
#endif

/*! \brief Creates an accumulation array of the out degree of all the vertices
           in a graph suitable for the manhattan partitioning functions.

    \param accum_deg The array holding the accumulated out degrees.  Must be
                     of size order + 1.
    \param g The graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/20/2012
*/
template <typename Graph, typename T>
inline
void accumulate_out_degree(T* accum_deg, Graph& g)
{
  #pragma mta noalias *accum_deg
  #pragma mta noalias g

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  accum_deg[0] = 0;

#ifdef USING_QTHREADS
  detail::qt_accum_out_deg<Graph, T> qaod(g, accum_deg);
  qt_loop_balance(0, order, qaod);
#else
  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      accum_deg[start_pos + 1] = out_degree(*tverts, g);
    }
  }
#endif

  for (size_type i = 1; i <= order; ++i) accum_deg[i] += accum_deg[i - 1];
}

/*! \brief Creates an accumulation array of the out degree of the vertices
           in a list suitable for the manhattan partitioning functions.

    \param accum_deg The array holding the accumulated out degrees.  Must be
                     of size vlist_size + 1.
    \param vlist The list of vertices from g.
    \param vlist_size The size of vlist.
    \param g The graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/20/2012
*/
template <typename Graph, typename T, typename T2>
inline
void
accumulate_out_degree(T* accum_deg,
                      typename graph_traits<Graph>::vertex_descriptor* vlist,
                      T2 vlist_size, Graph& g)
{
  #pragma mta noalias *accum_deg
  #pragma mta noalias *vlist
  #pragma mta noalias g

  accum_deg[0] = 0;

#ifdef USING_QTHREADS
  detail::qt_accum_out_deg_list<Graph, T> qaodl(g, accum_deg, vlist);
  qt_loop_balance(0, vlist_size, qaodl);
#else
  for (T2 i = 0; i < vlist_size; ++i)
  {
    accum_deg[i + 1] = out_degree(vlist[i], g);
  }
#endif

  for (T2 i = 1; i <= vlist_size; ++i) accum_deg[i] += accum_deg[i - 1];
}

/*! \brief Creates an accumulation array of the in degree of all the vertices
           in a graph suitable for the manhattan partitioning functions.

    \param accum_deg The array holding the accumulated in degrees.  Must be
                     of size order + 1.
    \param g The graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/20/2012
*/
template <typename Graph, typename T>
inline
void accumulate_in_degree(T* accum_deg, Graph& g)
{
  #pragma mta noalias *accum_deg
  #pragma mta noalias g

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  accum_deg[0] = 0;

#ifdef USING_QTHREADS
  detail::qt_accum_in_deg<Graph, T> qaid(g, accum_deg);
  qt_loop_balance(0, order, qaid);
#else
  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      accum_deg[start_pos + 1] = in_degree(*tverts, g);
    }
  }
#endif

  for (size_type i = 1; i <= order; ++i) accum_deg[i] += accum_deg[i - 1];
}

/*! \brief Creates an accumulation array of the in degree of the vertices
           in a list suitable for the manhattan partitioning functions.

    \param accum_deg The array holding the accumulated in degrees.  Must be
                     of size vlist_size + 1.
    \param vlist The list of vertices from g.
    \param vlist_size The size of vlist.
    \param g The graph.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/20/2012
*/
template <typename Graph, typename T, typename T2>
inline
void
accumulate_in_degree(T* accum_deg,
                     typename graph_traits<Graph>::vertex_descriptor* vlist,
                     T2 vlist_size, Graph& g)
{
  #pragma mta noalias *accum_deg
  #pragma mta noalias *vlist
  #pragma mta noalias g

  accum_deg[0] = 0;

#ifdef USING_QTHREADS
  detail::qt_accum_in_deg_list<Graph, T> qaidl(g, accum_deg, vlist);
  qt_loop_balance(0, vlist_size, qaidl);
#else
  for (T2 i = 0; i < vlist_size; ++i)
  {
    accum_deg[i + 1] = in_degree(vlist[i], g);
  }
#endif

  for (T2 i = 1; i <= vlist_size; ++i) accum_deg[i] += accum_deg[i - 1];
}

}

#endif
