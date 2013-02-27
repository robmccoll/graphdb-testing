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
/*! \file test_partitioning_schemes.cpp

    \brief Compares different manhattan partitioning schemes including the
           compiler's automatic one and manual ones.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/30/2011

    The outer loop iterates over a list of vertices, while the inner loop
    iterates over all the adjacent vertices of each vertex.
*/
/****************************************************************************/

#include <iostream>
#include <iomanip>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/random.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/partitioning.hpp>

using namespace mtgl;

//#define DEBUG
//#define TEST_STINGER

#define MY_NUM_STREAMS 3
#define MY_BLOCK_SIZE 1024

/// Compiler-generated manhattan loop collapse.
template <typename Graph>
void
count_adjacencies_higher_id(Graph& g, int* indeg,
                            typename graph_traits<Graph>::size_type* vlist,
                            typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;

  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  for (size_type i = 0; i < vsize; ++i)
  {
    size_type uid = vlist[i];
    size_type begin = index[uid];
    size_type end = index[uid + 1];

    for (size_type j = begin; j < end; ++j)
    {
      size_type vid = end_points[j];
      if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
      std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                << std::setw(3) << vid << std::endl;
#endif
    }
  }
}

/// Manual manhattan loop collapse where each stream performs two binary
/// searches to find it's beginning and ending outer iteration.
template <typename Graph>
void
count_adjacencies_higher_id2(Graph& g, int* indeg,
                             typename graph_traits<Graph>::size_type* vlist,
                             typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;

  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    size_type vid = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + index[vid + 1] - index[vid];
  }

  size_type stream_id = 0;
  size_type num_streams = MY_NUM_STREAMS;

#ifdef __MTA__
  #pragma mta for all streams stream_id of num_streams
#else
  for ( ; stream_id < num_streams; ++stream_id)
#endif
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            stream_id, num_streams);
    size_type end_pos = end_block_range(accum_deg[vsize],
                                        stream_id, num_streams);

//    std::cout << "stream_id: " << stream_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > start_pos.  Note that we need to search
    // the range [1, vsize + 1).
    size_type first = 1;
    size_type count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(start_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "start_first: " << first << "    accum_deg: "
//              << accum_deg[first] << std::endl;

    size_type start_outer = start_pos == 0 ? 0 : first - 1;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > end_pos.  Note that we need to search
    // the range [1, vsize + 1).
    first = 1;
    count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(end_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "end_first: " << first << std::endl;

    size_type end_outer = accum_deg[first - 1] == end_pos ? first - 1 : first;

//    std::cout << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = vlist[i];
      size_type start_inner = (i == start_outer) ? start_pos - accum_deg[i] : 0;
      size_type end_inner = (i + 1 == end_outer) ? end_pos - accum_deg[i] :
                            accum_deg[i + 1] - accum_deg[i];

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type vid = end_points[j + index[uid]];
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

/// Manual manhattan loop collapse where two binary searches are performed for
/// each block to find it's beginning and ending outer iteration.  This is
/// different from *id2 because it creates many more blocks than streams.
template <typename Graph>
void
count_adjacencies_higher_id3(Graph& g, int* indeg,
                             typename graph_traits<Graph>::size_type* vlist,
                             typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;

  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    size_type vid = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + index[vid + 1] - index[vid];
  }

#ifdef __MTA__
  size_type num_blocks = (accum_deg[vsize] + MY_BLOCK_SIZE - 1) / MY_BLOCK_SIZE;
#else
  size_type num_blocks = MY_NUM_STREAMS;
#endif

  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            block_id, num_blocks);
    size_type end_pos = end_block_range(accum_deg[vsize], block_id, num_blocks);

//    std::cout << "block_id: " << block_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > start_pos.  Note that we need to search
    // the range [1, vsize + 1).
    size_type first = 1;
    size_type count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(start_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "start_first: " << first << "    accum_deg: "
//              << accum_deg[first] << std::endl;

    size_type start_outer = start_pos == 0 ? 0 : first - 1;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > end_pos.  Note that we need to search
    // the range [1, vsize + 1).
    first = 1;
    count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(end_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "end_first: " << first << std::endl;

    size_type end_outer = accum_deg[first - 1] == end_pos ? first - 1 : first;

//    std::cout << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = vlist[i];
      size_type start_inner = (i == start_outer) ? start_pos - accum_deg[i] : 0;
      size_type end_inner = (i + 1 == end_outer) ? end_pos - accum_deg[i] :
                            accum_deg[i + 1] - accum_deg[i];

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type vid = end_points[j + index[uid]];
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

/// Manual manhattan loop collapse where two binary searches are performed for
/// each block to find it's beginning and ending outer iteration.  This
/// is the same as *id3 except that it uses the functions defined in
/// partitioning.hpp to encapsulate some of the logic and simplify the code.
template <typename Graph>
void
count_adjacencies_higher_id4(Graph& g, int* indeg,
                             typename graph_traits<Graph>::size_type* vlist,
                             typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;

  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    size_type vid = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + index[vid + 1] - index[vid];
  }

#ifdef __MTA__
  size_type num_blocks = (accum_deg[vsize] + MY_BLOCK_SIZE - 1) / MY_BLOCK_SIZE;
#else
  size_type num_blocks = MY_NUM_STREAMS;
#endif

  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            block_id, num_blocks);
    size_type end_pos = end_block_range(accum_deg[vsize], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, vsize, start_pos);
    size_type end_outer = end_manhattan_outer_range(accum_deg, vsize, end_pos);

//    std::cout << "block_id: " << block_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl
//              << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = vlist[i];

      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type vid = end_points[j + index[uid]];
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

/// Manual manhattan loop collapse where each inner iteration performs a
/// binary search to find it's outer iteration.
template <typename Graph>
void
count_adjacencies_higher_id_slow(Graph& g, int* indeg,
                                 typename graph_traits<Graph>::size_type* vlist,
                                 typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;

  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    size_type vid = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + index[vid + 1] - index[vid];
  }

  for (size_type k = 0; k < accum_deg[vsize]; ++k)
  {
    size_type first = 1;
    size_type count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(k < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

    size_type i = first - 1;
    size_type j = k - accum_deg[i];
    size_type uid = vlist[i];
    size_type ep = j + index[uid];
    size_type vid = end_points[ep];

    if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
    std::cout << "i: " << std::setw(3) << i << "    j: " << std::setw(3) << ep
              << "    uid: " << std::setw(3) << uid << "    vid: "
              << std::setw(3) << vid << std::endl;
#endif
  }

  free(accum_deg);
}

/// Compiler-generated manhattan loop collapse using MTGL interface.
template <typename Graph>
void count_adjacencies_higher_id_mtgl(
    Graph& g, int* indeg,
    typename graph_traits<Graph>::vertex_descriptor* vlist,
    typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  for (size_type i = 0; i < vsize; ++i)
  {
    size_type uid = get(vid_map, vlist[i]);
    adjacency_iterator adjs = adjacent_vertices(vlist[i], g);
    size_type end = out_degree(vlist[i], g);

    for (size_type j = 0; j < end; ++j)
    {
      size_type vid = get(vid_map, adjs[j]);
      if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
    }
  }
}

/// Manual manhattan loop collapse where each stream performs two binary
/// searches to find it's beginning and ending outer iteration.  Uses MTGL
/// interface.
template <typename Graph>
void
count_adjacencies_higher_id2_mtgl(
    Graph& g, int* indeg,
    typename graph_traits<Graph>::vertex_descriptor* vlist,
    typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    vertex_descriptor v = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + out_degree(v, g);
  }

  size_type stream_id = 0;
  size_type num_streams = MY_NUM_STREAMS;

#ifdef __MTA__
  #pragma mta for all streams stream_id of num_streams
#else
  for ( ; stream_id < num_streams; ++stream_id)
#endif
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            stream_id, num_streams);
    size_type end_pos = end_block_range(accum_deg[vsize],
                                        stream_id, num_streams);

//    std::cout << "stream_id: " << stream_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > start_pos.  Note that we need to search
    // the range [1, vsize + 1).
    size_type first = 1;
    size_type count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(start_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "start_first: " << first << "    accum_deg: "
//              << accum_deg[first] << std::endl;

    size_type start_outer = start_pos == 0 ? 0 : first - 1;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > end_pos.  Note that we need to search
    // the range [1, vsize + 1).
    first = 1;
    count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(end_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "end_first: " << first << std::endl;

    size_type end_outer = accum_deg[first - 1] == end_pos ? first - 1 : first;

//    std::cout << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = get(vid_map, vlist[i]);
      size_type start_inner = (i == start_outer) ? start_pos - accum_deg[i] : 0;
      size_type end_inner = (i + 1 == end_outer) ? end_pos - accum_deg[i] :
                            accum_deg[i + 1] - accum_deg[i];

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      adjacency_iterator adjs = adjacent_vertices(vlist[i], g);

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type vid = get(vid_map, adjs[j]);
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

/// Manual manhattan loop collapse where two binary searches are performed for
/// each block to find it's beginning and ending outer iteration.  This is
/// different from *id2_mtgl because it creates many more blocks than streams.
/// Uses MTGL interface.
template <typename Graph>
void
count_adjacencies_higher_id3_mtgl(
    Graph& g, int* indeg,
    typename graph_traits<Graph>::vertex_descriptor* vlist,
    typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    vertex_descriptor v = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + out_degree(v, g);
  }

#ifdef __MTA__
  size_type num_blocks = (accum_deg[vsize] + MY_BLOCK_SIZE - 1) / MY_BLOCK_SIZE;
#else
  size_type num_blocks = MY_NUM_STREAMS;
#endif

  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            block_id, num_blocks);
    size_type end_pos = end_block_range(accum_deg[vsize], block_id, num_blocks);

//    std::cout << "block_id: " << block_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > start_pos.  Note that we need to search
    // the range [1, vsize + 1).
    size_type first = 1;
    size_type count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(start_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "start_first: " << first << "    accum_deg: "
//              << accum_deg[first] << std::endl;

    size_type start_outer = start_pos == 0 ? 0 : first - 1;

    // Binary search accum_deg to find the position in the list that
    // contains the first value > end_pos.  Note that we need to search
    // the range [1, vsize + 1).
    first = 1;
    count = vsize;
    while (count > 0)
    {
      size_type step = count / 2;
      size_type mid = first + step;
      if (!(end_pos < accum_deg[mid]))
      {
        first = mid + 1;
        count -= step + 1;
      }
      else
      {
        count = step;
      }
    }

//    std::cout << "end_first: " << first << std::endl;

    size_type end_outer = accum_deg[first - 1] == end_pos ? first - 1 : first;

//    std::cout << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = get(vid_map, vlist[i]);
      size_type start_inner = (i == start_outer) ? start_pos - accum_deg[i] : 0;
      size_type end_inner = (i + 1 == end_outer) ? end_pos - accum_deg[i] :
                            accum_deg[i + 1] - accum_deg[i];

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      adjacency_iterator adjs = adjacent_vertices(vlist[i], g);

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type vid = get(vid_map, adjs[j]);
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

/// Manual manhattan loop collapse where two binary searches are performed for
/// each block to find it's beginning and ending outer iteration.  This
/// is the same as *id3_mtgl except that it uses the functions defined in
/// partitioning.hpp to encapsulate some of the logic and simplify the code.
/// Uses MTGL interface.
template <typename Graph>
void
count_adjacencies_higher_id4_mtgl(
    Graph& g, int* indeg,
    typename graph_traits<Graph>::vertex_descriptor* vlist,
    typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  // Create the accumulation array of the inner dimension.
  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    vertex_descriptor v = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + out_degree(v, g);
  }

#ifdef __MTA__
  size_type num_blocks = (accum_deg[vsize] + MY_BLOCK_SIZE - 1) / MY_BLOCK_SIZE;
#else
  size_type num_blocks = MY_NUM_STREAMS;
#endif

  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            block_id, num_blocks);
    size_type end_pos = end_block_range(accum_deg[vsize], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, vsize, start_pos);
    size_type end_outer = end_manhattan_outer_range(accum_deg, vsize, end_pos);

//    std::cout << "block_id: " << block_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl
//              << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = get(vid_map, vlist[i]);

      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      adjacency_iterator adjs = adjacent_vertices(vlist[i], g);

      for (size_type j = start_inner; j != end_inner; ++j)
      {
        size_type vid = get(vid_map, adjs[j]);
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

/// Manual manhattan loop collapse where two binary searches are performed for
/// each block to find it's beginning and ending outer iteration.  This is
/// different from *id4_mtgl because it uses thread iterators.  Uses MTGL
/// interface.
template <typename Graph>
void
count_adjacencies_higher_id5_mtgl(
    Graph& g, int* indeg,
    typename graph_traits<Graph>::vertex_descriptor* vlist,
    typename graph_traits<Graph>::size_type vsize)
{
  #pragma mta noalias g
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *vlist

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  size_type* accum_deg = (size_type*) malloc((vsize + 1) * sizeof(size_type));
  accum_deg[0] = 0;
  for (size_type i = 0; i < vsize; ++i)
  {
    vertex_descriptor v = vlist[i];
    accum_deg[i + 1] = accum_deg[i] + out_degree(v, g);
  }

#ifdef __MTA__
  size_type num_blocks = (accum_deg[vsize] + MY_BLOCK_SIZE - 1) / MY_BLOCK_SIZE;
#else
  size_type num_blocks = MY_NUM_STREAMS;
#endif

  #pragma mta assert nodep
  #pragma mta no scalar expansion
  #pragma mta dynamic schedule
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos = begin_block_range(accum_deg[vsize],
                                            block_id, num_blocks);
    size_type end_pos = end_block_range(accum_deg[vsize], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, vsize, start_pos);
    size_type end_outer = end_manhattan_outer_range(accum_deg, vsize, end_pos);

//    std::cout << "block_id: " << block_id << "    start_pos: " << start_pos
//              << "    end_pos: " << end_pos << std::endl
//              << "start_outer: " << start_outer << "    end_outer: "
//              << end_outer << std::endl;

    for (size_type i = start_outer; i != end_outer; ++i)
    {
      size_type uid = get(vid_map, vlist[i]);

      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

//      std::cout << "start_inner: " << start_inner << "    end_inner: "
//                << end_inner << std::endl;

      thread_adjacency_iterator adjs =
        thread_adjacent_vertices(vlist[i], start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++adjs)
      {
        size_type vid = get(vid_map, *adjs);
        if (uid < vid) mt_incr(indeg[vid], 1);

#ifdef DEBUG
        std::cout << "uid: " << std::setw(3) << uid << "    vid: "
                  << std::setw(3) << vid << std::endl;
#endif
      }
    }
  }

  free(accum_deg);
}

template <typename T, typename T2>
void checkError(T* indeg, T* indeg2, T2 order)
{
  #pragma mta assert noalias *indeg
  #pragma mta assert noalias *indeg2

  T2 error = std::numeric_limits<T2>::max();

  for (T2 i = 0; i < order; ++i)
  {
    if (indeg[i] != indeg2[i]) error = i;
  }

  if (error != std::numeric_limits<T2>::max())
  {
    std::cout << "Error in computation: pos " << error << std::endl;
  }
}

int main(int argc, char* argv[])
{
#ifdef TEST_STINGER
  typedef stinger_graph<int64_t, int64_t, int64_t> SGraph;
  typedef stinger_graph_adapter<SGraph> Graph;
#else
  typedef compressed_sparse_row_graph<directedS> Graph;
//  typedef adjacency_list<directedS> Graph;
#endif
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::thread_vertex_iterator thread_vertex_iterator;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;

  mt_srand48(0);

  init_test(argc, argv);

#ifdef TEST_STINGER
  SGraph sg(1);
  Graph g(sg);
#else
  Graph g;
#endif

  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  std::cout << "order: " << order << ",  size: " << size << std::endl;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  int* indeg = new int[order];
  int* indeg2 = new int[order];

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  size_type vsize = order;
  vertex_descriptor* vlist = new vertex_descriptor[vsize];

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(vsize, stream_id, num_streams);
    size_type end_pos = end_block_range(vsize, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vlist[start_pos] = *verts;
    }
  }

  mt_timer timer;
  int issues, memrefs, concur, streams, traps;
#ifdef __MTA__
  int phantoms, ready;
#endif

  // Get the ground truth.
#ifdef TEST_STINGER
  count_adjacencies_higher_id5_mtgl(g, indeg, vlist, vsize);
#else
  count_adjacencies_higher_id_mtgl(g, indeg, vlist, vsize);
#endif

  int sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

#if 0
  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id2(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id2():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id3(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id3():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id4(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id4():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);
#endif

#if 0
  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id_slow(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id5():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id_mtgl(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id_mtgl():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id2_mtgl(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id2_mtgl():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id3_mtgl(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id3_mtgl():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id4_mtgl(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id4_mtgl():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);
#endif

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  #pragma mta fence
  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

  count_adjacencies_higher_id5_mtgl(g, indeg2, vlist, vsize);

  #pragma mta fence
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  std::cout << "count_adjacencies_higher_id5_mtgl():" << std::endl;
  print_mta_counters(timer, size, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  std::cout << "phantoms: " << phantoms << std::endl
            << "ready: " << ready << std::endl;
#endif

  sum = 0;
  for (size_type i = 0; i < order; ++i) sum += indeg[i];
  std::cout << "Sum: " << sum << std::endl << std::endl;

  checkError(indeg, indeg2, order);

  delete [] indeg;
  delete [] indeg2;
  delete [] vlist;

  return 0;
}
