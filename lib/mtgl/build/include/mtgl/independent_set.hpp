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
/*! \file independent_set.hpp

    \brief This file contains multithreaded independent set and maximal
           independent set algorithms.  These are currently set to favor
           high degree vertices in the independent set rather than the most
           vertices.  This affects the performance of the triangle finding
           algorithm, for example.

    \author Jon Berry (jberry@sandia.gov)
            Greg Mackey (gemacke@sandia.gov)

    \date 2007
*/
/****************************************************************************/

#ifndef MTGL_INDEPENDENT_SET_HPP
#define MTGL_INDEPENDENT_SET_HPP

#include <climits>
#include <iostream>
#include <iomanip>

#define IS_CHUNK 256

namespace mtgl {

template <typename Graph, typename VertexMap>
void validate_maximal_independent_set(Graph& g, VertexMap& ind_set)
{
  #pragma mta noalias g
  #pragma mta noalias ind_set

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type order = num_vertices(g);

  const int ERRMAX = 10;

  mtgl::pair<size_type, size_type>* conflicts =
    new mtgl::pair<size_type, size_type>[ERRMAX];

  vertex_property_map<Graph, unsigned long> isolated(g);

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
      isolated[*tverts] = 0;
    }
  }

  // Get the accumulation array of the vertex out degrees.
  size_type* accum_deg = new size_type[order + 1];

  accumulate_out_degree(accum_deg, g);

#ifdef USING_QTHREADS
  size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
  size_type num_blocks = (accum_deg[order] + IS_CHUNK - 1) / IS_CHUNK;
#else
  size_type num_blocks = 1;
#endif

  int num_conflicts = 0;

  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos =
      begin_block_range(accum_deg[order], block_id, num_blocks);
    size_type end_pos =
      end_block_range(accum_deg[order], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start_pos);
    size_type end_outer =
      end_manhattan_outer_range(accum_deg, order, end_pos);

    thread_vertex_iterator tverts = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

      vertex_descriptor u = *tverts;
      size_type uid = get(vid_map, u);
      thread_adjacency_iterator tadj_verts =
        thread_adjacent_vertices(u, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
      {
        vertex_descriptor v = *tadj_verts;
        size_type vid = get(vid_map, v);

        if (uid == vid) continue;

        if (ind_set[v]) mt_incr(isolated[u], 1);
        if (ind_set[u]) mt_incr(isolated[v], 1);

        if (ind_set[u] && ind_set[v] && num_conflicts < ERRMAX)
        {
          int pos = mt_incr(num_conflicts, 1);
          conflicts[pos] = mtgl::pair<size_type, size_type>(uid, vid);
        }
      }
    }
  }

  for (int i = 0; i < num_conflicts; ++i)
  {
    std::cout << " conflict: " << std::setw(8) << conflicts[i].first << "  "
              << std::setw(8) << conflicts[i].second << std::endl;
  }

  thread_vertex_iterator tverts = thread_vertices(0, g);
  for (size_type i = 0; i < order; ++i, ++tverts)
  {
    vertex_descriptor v = *tverts;

    if (!ind_set[v] && isolated[v] == 0 && out_degree(v, g) > 0)
    {
      std::cout << " isolated: " << std::setw(8) << i << std::endl;
    }
  }

  delete [] accum_deg;
  delete [] conflicts;
}

void print_independent_set(bool* active, int size)
{
  int num_in_line = 0;
  for (int i = 0; i < size; ++i)
  {
    if (active[i])
    {
      std::cout << " " << std::setw(7) << i;
      ++num_in_line;

      if (num_in_line == 9)
      {
        std::cout << std::endl;
        num_in_line = 0;
      }
    }
  }

  std::cout << std::endl;
  if (num_in_line != 0) std::cout << std::endl;
}

template <typename Graph, typename VertexMap>
typename graph_traits<Graph>::size_type
independent_set(Graph& g, VertexMap& ind_set,
                typename graph_traits<Graph>::size_type* accum_deg,
                typename graph_traits<Graph>::size_type num_blocks)
{
  #pragma mta noalias g
  #pragma mta noalias ind_set
  #pragma mta noalias *accum_deg

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  size_type order = num_vertices(g);
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos =
      begin_block_range(accum_deg[order], block_id, num_blocks);
    size_type end_pos =
      end_block_range(accum_deg[order], block_id, num_blocks);

    size_type start_outer =
      begin_manhattan_outer_range(accum_deg, order, start_pos);
    size_type end_outer =
      end_manhattan_outer_range(accum_deg, order, end_pos);

    thread_vertex_iterator tverts = thread_vertices(start_outer, g);

    for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
    {
      size_type start_inner =
        begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
      size_type end_inner =
        end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

      vertex_descriptor u = *tverts;
      size_type uid = get(vid_map, u);
      size_type u_out_deg = out_degree(u, g);
      thread_adjacency_iterator tadj_verts =
        thread_adjacent_vertices(u, start_inner, g);

      for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
      {
        vertex_descriptor v = *tadj_verts;
        size_type vid = get(vid_map, v);
        size_type v_out_deg = out_degree(v, g);

        // The winner is the vertex with the lowest (JWB) out degree with ties
        // being broken by the winner being the vertex with the smaller id.
        // The loser is knocked out of the independent set.
        if (ind_set[u] && ind_set[v])
        {
          if (u_out_deg < v_out_deg)
          {
            ind_set[v] = false;
          }
          else if (u_out_deg > v_out_deg)
          {
            ind_set[u] = false;
          }
          else if (uid < vid)
          {
            ind_set[v] = false;
          }
          else if (uid > vid)
          {
            ind_set[u] = false;
          }
        }
      }
    }
  }

  size_type issize = 0;

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type my_issize = 0;

    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;
      my_issize += ind_set[v];
    }

    mt_incr(issize, my_issize);
  }

  return issize;
}

template <typename Graph, typename VertexMap>
typename graph_traits<Graph>::size_type
independent_set(Graph& g, VertexMap& ind_set)
{
  #pragma mta noalias g
  #pragma mta noalias ind_set

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // Get the accumulation array of the vertex out degrees.
  size_type* accum_deg = new size_type[order + 1];

  accumulate_out_degree(accum_deg, g);

#ifdef USING_QTHREADS
  size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
  size_type num_blocks = (accum_deg[order] + IS_CHUNK - 1) / IS_CHUNK;
#else
  size_type num_blocks = 1;
#endif

  size_type retval = independent_set(g, ind_set, accum_deg, num_blocks);

  delete [] accum_deg;

  return retval;
}

template <typename Graph, typename VertexMap>
typename graph_traits<Graph>::size_type
maximal_independent_set(Graph& g, VertexMap& ind_set)
{
  #pragma mta noalias g
  #pragma mta noalias ind_set

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  size_type order = num_vertices(g);

  vertex_property_map<Graph, bool> active_v(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // The global independent set is initially empty.
  // The active graph initially contains all the vertices.
  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;

      ind_set[v] = false;
      active_v[v] = true;
    }
  }

  // Get the accumulation array of the vertex out degrees.
  size_type* accum_deg = new size_type[order + 1];

  accumulate_out_degree(accum_deg, g);

#ifdef USING_QTHREADS
  size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
  size_type num_blocks = (accum_deg[order] + IS_CHUNK - 1) / IS_CHUNK;
#else
  size_type num_blocks = 1;
#endif

  size_type issize = 0;

  while (1)
  {
    // Find an independent set in the active graph.
    size_type next_issize = independent_set(g, active_v, accum_deg, num_blocks);

#ifdef DEBUG
    std::cout << "New independent set: " << next_issize << std::endl;
    print_independent_set(active_v, order);
#endif

    // Set the global independent set to the union of itself and the new
    // independent set in the active graph.
    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;
        ind_set[v] = active_v[v] | ind_set[v];
      }
    }

    if (next_issize <= 0) break;

    issize += next_issize;

    // Set the active graph to the opposite of the global independent set.
    // The active graph vertices will be all those not in the global
    // independent set.
    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;
        active_v[v] = !ind_set[v];
      }
    }

    // Deactivate the neighbors of global independent set vertices in the
    // active graph.
    #pragma mta assert parallel
    #pragma mta block dynamic schedule
    #pragma mta assert nodep
    #pragma mta no scalar expansion
    for (size_type block_id = 0; block_id < num_blocks; ++block_id)
    {
      size_type start_pos =
        begin_block_range(accum_deg[order], block_id, num_blocks);
      size_type end_pos =
        end_block_range(accum_deg[order], block_id, num_blocks);

      size_type start_outer =
        begin_manhattan_outer_range(accum_deg, order, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_deg, order, end_pos);

      thread_vertex_iterator tverts = thread_vertices(start_outer, g);

      for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

        vertex_descriptor u = *tverts;
        thread_adjacency_iterator tadj_verts =
          thread_adjacent_vertices(u, start_inner, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor v = *tadj_verts;
          if (ind_set[u] && active_v[v]) active_v[v] = false;
          if (ind_set[v] && active_v[u]) active_v[u] = false;
        }
      }
    }
  }

  delete [] accum_deg;

  return issize;
}

// TODO: pull out common denominator
template <typename Graph, typename VertexMap>
typename graph_traits<Graph>::size_type
maximal_independent_set(Graph& g, VertexMap& ind_set,
                        VertexMap& global_active_v)
{
  #pragma mta noalias g
  #pragma mta noalias ind_set
  #pragma mta noalias global_active_v

  // perform the whole computation in the subgraph induced by global_active_v
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  size_type order = num_vertices(g);

  vertex_property_map<Graph, bool> active_v(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  // The global independent set is initially empty.
  // The active graph initially contains all the vertices.
  #pragma mta assert nodep
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      vertex_descriptor v = *tverts;

      ind_set[v] = false;
      active_v[v] = global_active_v[v];
    }
  }

  // Get the accumulation array of the vertex out degrees.
  size_type* accum_deg = new size_type[order + 1];

  accum_deg[0] = 0;

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

  for (size_type i = 1; i <= order; ++i) accum_deg[i] += accum_deg[i - 1];

#ifdef USING_QTHREADS
  size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
  size_type num_blocks = (accum_deg[order] + IS_CHUNK - 1) / IS_CHUNK;
#else
  size_type num_blocks = 1;
#endif

  size_type issize = 0;

  while (1)
  {
    // Find an independent set in the active graph.
    size_type next_issize = independent_set(g, active_v, accum_deg, num_blocks);

#ifdef DEBUG
    std::cout << "New independent set: " << next_issize << std::endl;
    print_independent_set(active_v, order);
#endif

    // Set the global independent set to the union of itself and the new
    // independent set in the active graph.
    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;
        // active_v holds the result of the most recent independent set
        // computation.
        ind_set[v] = global_active_v[v] && (active_v[v] | ind_set[v]);
      }
    }

    if (next_issize <= 0) break;

    issize += next_issize;

    // Set the active graph to the opposite of the global independent set,
    // with respect to global_active_v.
    // The active graph vertices will be all those not in the global
    // independent set, but in global_active_v.
    #pragma mta assert nodep
    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;
        active_v[v] = (global_active_v[v] && !ind_set[v]);
      }
    }

    // Deactivate the neighbors of global independent set vertices in the
    // active graph.
    #pragma mta assert parallel
    #pragma mta block dynamic schedule
    #pragma mta assert nodep
    #pragma mta no scalar expansion
    for (size_type block_id = 0; block_id < num_blocks; ++block_id)
    {
      size_type start_pos =
        begin_block_range(accum_deg[order], block_id, num_blocks);
      size_type end_pos =
        end_block_range(accum_deg[order], block_id, num_blocks);

      size_type start_outer =
        begin_manhattan_outer_range(accum_deg, order, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_deg, order, end_pos);

      thread_vertex_iterator tverts = thread_vertices(start_outer, g);

      for (size_type i = start_outer; i != end_outer; ++i, ++tverts)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

        vertex_descriptor u = *tverts;
        thread_adjacency_iterator tadj_verts =
          thread_adjacent_vertices(u, start_inner, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++tadj_verts)
        {
          vertex_descriptor v = *tadj_verts;
          if (ind_set[u] && active_v[v]) active_v[v] = false;
          if (ind_set[v] && active_v[u]) active_v[u] = false;
        }
      }
    }
  }

  delete [] accum_deg;

  return issize;
}

}

#endif
