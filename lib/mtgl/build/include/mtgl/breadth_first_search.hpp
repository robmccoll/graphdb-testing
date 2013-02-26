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
/*! \file breadth_first_search.hpp

    \brief  This is the version of breadth-first search that currently
            performs the best on graphs with power law degree distributions.
            It is also the version used in the 2009 MTAAP paper on Qthreads
            and the MTGL.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 9/15/2010
*/
/****************************************************************************/

#ifndef MTGL_BREADTH_FIRST_SEARCH_HPP
#define MTGL_BREADTH_FIRST_SEARCH_HPP

#include <limits>
#include <iostream>
#include <iomanip>

#include <mtgl/util.hpp>
#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/partitioning.hpp>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

#define BFS_CHUNK 256

namespace mtgl {

template <typename Graph>
class default_bfs_visitor {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  default_bfs_visitor() {}

  void discover_vertex(vertex_descriptor& u, Graph& g) {}
  void examine_vertex(vertex_descriptor& u, Graph& g) {}
  void examine_edge(edge_descriptor& e, Graph &g) {}
  void tree_edge(edge_descriptor& e, Graph &g) {}
  void non_tree_edge(edge_descriptor& e, Graph &g) {}
  void finish_vertex(vertex_descriptor& u, Graph& g) {}
  bool visit_test(bool isWhite, edge_descriptor& e, Graph &g)
  { return isWhite; }
};

template <typename Graph, typename PredMap>
class parents_bfs_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  parents_bfs_visitor(PredMap& p) : parents(p) {}

  void tree_edge(edge_descriptor& e, Graph &g)
  {
    vertex_descriptor u = source(e, g);
    vertex_descriptor v = target(e, g);
    put(parents, v, u);
  }

private:
  PredMap& parents;
};

#ifdef USING_QTHREADS
namespace detail {

template <typename Graph, typename BFSVisitor, typename ColorMap,
          typename Queue>
class eoe_loop {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_out_edge_iterator
          thread_out_edge_iterator;

  eoe_loop(Graph& gg, BFSVisitor& v, size_type nb, size_type bs,
           size_type ne, size_type& t, vertex_descriptor* b,
           vertex_descriptor* tv, Queue& q, ColorMap& c, size_type* ad) :
      g(gg), vis(v), num_blocks(nb), block_size(bs), num_elements(ne),
      tail(t), buffer(b), to_visit(tv), Q(q), color(c), accum_deg(ad) {}

  void operator()(size_t start_at, size_t stop_at)
  {
    for (size_t block_id = start_at; block_id < stop_at; ++block_id)
    {
      BFSVisitor my_vis = vis;

      size_type start_pos =
        begin_block_range(accum_deg[num_elements], block_id, num_blocks);
      size_type end_pos =
        end_block_range(accum_deg[num_elements], block_id, num_blocks);

      size_type start_outer =
        begin_manhattan_outer_range(accum_deg, num_elements, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_deg, num_elements, end_pos);

      size_type my_count = 0;
      size_type my_start = block_id * block_size;
      vertex_descriptor* my_buffer = buffer + my_start;

      for (size_type i = start_outer; i != end_outer; ++i)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

        vertex_descriptor u = to_visit[i];
        thread_out_edge_iterator oedges = thread_out_edges(u, start_inner, g);

        if (start_inner == 0) my_vis.examine_vertex(u, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++oedges)
        {
          edge_descriptor e = *oedges;
          vertex_descriptor v = target(e, g);
          size_type vcolor = mt_incr(color[v], 1);

          my_vis.examine_edge(e, g);

          if (my_vis.visit_test(vcolor == 0, e, g))
          {
            my_vis.tree_edge(e, g);
            my_vis.discover_vertex(v, g);
            my_buffer[my_count++] = v;
          }
          else
          {
            my_vis.non_tree_edge(e, g);
          }
        }

        if (end_inner == accum_deg[i+1] - accum_deg[i])
        {
          my_vis.finish_vertex(u, g);
        }
      }

      size_type my_q_ptr = mt_incr(tail, my_count);

      for (size_type i = 0; i < my_count; ++i) Q[my_q_ptr++] = my_buffer[i];
    }
  }

private:
  Graph& g;
  BFSVisitor& vis;
  size_type num_blocks;
  size_type block_size;
  size_type num_elements;
  size_type& tail;
  vertex_descriptor* buffer;
  vertex_descriptor* to_visit;
  Queue& Q;
  ColorMap& color;
  size_type* accum_deg;
};

}
#endif

template <typename Graph, typename BFSVisitor, typename ColorMap,
          typename Queue>
void
inline
expand_one_edge(Graph& g,
                Queue& Q,
                typename graph_traits<Graph>::size_type& head,
                typename graph_traits<Graph>::size_type& tail,
                BFSVisitor& vis, ColorMap& color,
                typename graph_traits<Graph>::vertex_descriptor*& buffer,
                typename graph_traits<Graph>::size_type& buf_size,
                typename graph_traits<Graph>::size_type*& accum_deg,
                typename graph_traits<Graph>::size_type& accum_deg_size)
{
  #pragma mta noalias g
  #pragma mta noalias Q
  #pragma mta noalias vis
  #pragma mta noalias color
  #pragma mta noalias *buffer
  #pragma mta noalias *accum_deg

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_out_edge_iterator
          thread_out_edge_iterator;

  size_type num_elements = tail - head;

#ifdef PHASE_DEBUG
  std::cout << "Expanding " << std::setw(10) << num_elements
            << " vertices (head: " << head << ", tail: " << tail << ")"
            << std::endl;
  mt_timer phase_timer;
  phase_timer.start();
#endif

  vertex_descriptor* to_visit = &Q[head];

  if (accum_deg_size < num_elements + 1)
  {
    if (accum_deg_size > 0) free(accum_deg);
    accum_deg_size = num_elements + 1;
    accum_deg = (size_type*) malloc(sizeof(size_type) * accum_deg_size);
  }

  accumulate_out_degree(accum_deg, to_visit, num_elements, g);

  head = tail;

  size_type work_this_phase = accum_deg[num_elements];

  if (work_this_phase > 0)
  {
#ifdef USING_QTHREADS
    size_type num_blocks = qthread_num_shepherds();
#elif defined(__MTA__)
    size_type num_blocks = (work_this_phase + BFS_CHUNK - 1) / BFS_CHUNK;
#else
    size_type num_blocks = 1;
#endif

    size_type block_size = (work_this_phase + num_blocks - 1) / num_blocks;

    if (buf_size < num_blocks * block_size)
    {
      if (buf_size > 0) free(buffer);
      buf_size = num_blocks * block_size;
      buffer =
        (vertex_descriptor*) malloc(sizeof(vertex_descriptor) * buf_size);
    }

#ifdef PHASE_DEBUG
    phase_timer.stop();

    double total_time = phase_timer.getElapsedSeconds();

    std::cout << "work_this_phase: " << work_this_phase << std::endl
              << "     num_blocks: " << num_blocks << std::endl
              << "Load Balance Time: " << std::setw(9) << std::setprecision(6)
              << std::fixed << total_time << std::endl;

    phase_timer.start();
#endif

#ifdef USING_QTHREADS
    detail::eoe_loop<Graph, BFSVisitor, ColorMap, Queue>
      eoel(g, vis, num_blocks, block_size, num_elements, tail,
           buffer, to_visit, Q, color, accum_deg);
    qt_loop_balance(0, num_blocks, eoel);
#else
    #pragma mta assert parallel
    #pragma mta assert nodep
    #pragma mta no scalar expansion
    for (size_type block_id = 0; block_id < num_blocks; ++block_id)
    {
      size_type start_pos =
        begin_block_range(work_this_phase, block_id, num_blocks);
      size_type end_pos =
        end_block_range(work_this_phase, block_id, num_blocks);

      size_type start_outer =
        begin_manhattan_outer_range(accum_deg, num_elements, start_pos);
      size_type end_outer =
        end_manhattan_outer_range(accum_deg, num_elements, end_pos);

      size_type my_count = 0;
      size_type my_start = block_id * block_size;
      vertex_descriptor* my_buffer = buffer + my_start;

      for (size_type i = start_outer; i != end_outer; ++i)
      {
        size_type start_inner =
          begin_manhattan_inner_range(accum_deg, start_pos, start_outer, i);
        size_type end_inner =
          end_manhattan_inner_range(accum_deg, end_pos, end_outer, i);

        vertex_descriptor u = to_visit[i];
        thread_out_edge_iterator oedges = thread_out_edges(u, start_inner, g);

        if (start_inner == 0) vis.examine_vertex(u, g);

        for (size_type j = start_inner; j != end_inner; ++j, ++oedges)
        {
          edge_descriptor e = *oedges;
          vertex_descriptor v = target(e, g);
          size_type vcolor = mt_incr(color[v], 1);

          vis.examine_edge(e, g);

          if (vis.visit_test(vcolor == 0, e, g))
          {
            vis.tree_edge(e, g);
            vis.discover_vertex(v, g);
            my_buffer[my_count++] = v;
          }
          else
          {
            vis.non_tree_edge(e, g);
          }
        }

        if (end_inner == accum_deg[i+1] - accum_deg[i])
        {
          vis.finish_vertex(u, g);
        }
      }

      size_type my_q_ptr = mt_incr(tail, my_count);

      for (size_type i = 0; i < my_count; ++i) Q[my_q_ptr++] = my_buffer[i];
    }
#endif

#ifdef PHASE_DEBUG
  phase_timer.stop();
  total_time += phase_timer.getElapsedSeconds();
  std::cout << "      Search Time: " << std::setw(9) << std::setprecision(6)
            << std::fixed << phase_timer.getElapsedSeconds() << std::endl
            << "       Phase Time: " << std::setw(9) << std::setprecision(6)
            << std::fixed << total_time << std::endl;
#endif
  }
}

template <typename Graph, typename BFSVisitor, typename ColorMap,
          typename Queue>
inline
void
expand_one_edge(Graph& g,
                Queue& Q,
                typename graph_traits<Graph>::size_type& head,
                typename graph_traits<Graph>::size_type& tail,
                BFSVisitor& vis, ColorMap& color)
{
  #pragma mta noalias g
  #pragma mta noalias Q
  #pragma mta noalias vis
  #pragma mta noalias color

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  vertex_descriptor* buffer = 0;
  size_type buf_size = 0;
  size_type* accum_deg = 0;
  size_type accum_deg_size = 0;

  expand_one_edge(g, Q, head, tail, vis, color,
                  buffer, buf_size, accum_deg, accum_deg_size);

  if (buf_size > 0) free(buffer);
  if (accum_deg_size > 0) free(accum_deg);
}

template <typename Graph, typename BFSVisitor, typename ColorMap,
          typename Queue>
typename graph_traits<Graph>::size_type
inline
breadth_first_search(Graph& g,
                     typename graph_traits<Graph>::vertex_descriptor root,
                     BFSVisitor& vis, ColorMap& color, Queue& Q,
                     typename graph_traits<Graph>::vertex_descriptor*& buffer,
                     typename graph_traits<Graph>::size_type& buf_size,
                     typename graph_traits<Graph>::size_type*& accum_deg,
                     typename graph_traits<Graph>::size_type& accum_deg_size)
{
  #pragma mta noalias g
  #pragma mta noalias vis
  #pragma mta noalias color
  #pragma mta noalias Q
  #pragma mta noalias *buffer
  #pragma mta noalias *accum_deg

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  size_type maxDist = 0;

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts) put(color, *verts, 0);
  }

  size_type head = 0;
  size_type tail = 1;

  // Put the root in the queue and set its color.
  vis.discover_vertex(root, g);
  Q[head] = root;
  put(color, root, 1);

  while (tail > head)
  {
    expand_one_edge(g, Q, head, tail, vis, color, buffer, buf_size,
                    accum_deg, accum_deg_size);

    ++maxDist;
  }

#ifdef DEBUG
  size_type visited_edges = 0;

  stream_id = 0;
  num_streams = 1;
  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      visited_edges += (get(color, v) > 0) * out_degree(v, g);
    }
  }

#ifdef PHASE_DEBUG
  std::cout << std::endl;
#endif

  std::cout << "Number of levels: " << maxDist << std::endl
            << "Visited vertices: " << tail << std::endl
            << "   Visited edges: " << visited_edges << std::endl;
#endif

  return tail;
}

template <typename Graph, typename BFSVisitor, typename ColorMap,
          typename Queue>
typename graph_traits<Graph>::size_type
breadth_first_search(Graph& g,
                     typename graph_traits<Graph>::vertex_descriptor root,
                     BFSVisitor& vis, ColorMap& color, Queue& Q)
{
  #pragma mta noalias g
  #pragma mta noalias vis
  #pragma mta noalias color
  #pragma mta noalias *Q

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  vertex_descriptor* buffer = 0;
  size_type buf_size = 0;
  size_type* accum_deg = 0;
  size_type accum_deg_size = 0;

  size_type retval = breadth_first_search(g, root, vis, color, Q, buffer,
                                          buf_size, accum_deg, accum_deg_size);

  if (buf_size > 0) free(buffer);
  free(accum_deg);

  return retval;
}

template <typename Graph, typename BFSVisitor, typename ColorMap>
typename graph_traits<Graph>::size_type
breadth_first_search(Graph& g,
                     typename graph_traits<Graph>::vertex_descriptor root,
                     BFSVisitor& vis, ColorMap& color)
{
  #pragma mta noalias g
  #pragma mta noalias vis
  #pragma mta noalias color

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  vertex_descriptor* Q =
    (vertex_descriptor*) malloc(num_vertices(g) * sizeof(vertex_descriptor));

  vertex_descriptor* buffer = 0;
  size_type buf_size = 0;
  size_type* accum_deg = 0;
  size_type accum_deg_size = 0;

  size_type retval = breadth_first_search(g, root, vis, color, Q, buffer,
                                          buf_size, accum_deg, accum_deg_size);

  free(Q);
  if (buf_size > 0) free(buffer);
  free(accum_deg);

  return retval;
}

template <typename Graph, typename BFSVisitor>
typename graph_traits<Graph>::size_type
breadth_first_search(Graph& g,
                     typename graph_traits<Graph>::vertex_descriptor root,
                     BFSVisitor& vis)
{
  #pragma mta noalias g
  #pragma mta noalias vis

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  vertex_property_map<Graph, size_type> color(g);

  vertex_descriptor* Q =
    (vertex_descriptor*) malloc(num_vertices(g) * sizeof(vertex_descriptor));

  vertex_descriptor* buffer = 0;
  size_type buf_size = 0;
  size_type* accum_deg = 0;
  size_type accum_deg_size = 0;

  size_type retval = breadth_first_search(g, root, vis, color, Q, buffer,
                                          buf_size, accum_deg, accum_deg_size);

  free(Q);
  if (buf_size > 0) free(buffer);
  free(accum_deg);

  return retval;
}

template <typename Graph, typename PredMap>
void
tree_check(Graph& g, typename graph_traits<Graph>::vertex_descriptor root,
           PredMap& bfs_parent)
{
  #pragma mta noalias g
  #pragma mta noalias bfs_parent

  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_adjacency_iterator
          thread_adjacency_iterator;

  std::cout << "Checking that BFS tree is valid." << std::endl;

  size_type order = num_vertices(g);
  vertex_descriptor null_vert = null_vertex(g);

  typedef vertex_property_map<Graph, size_type> VertexPropMap;
  VertexPropMap dist(g);

  size_type stream_id = 0;
  size_type num_streams = 1;

  #pragma mta for all streams stream_id of num_streams
  {
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      dist[*verts] = ULONG_MAX;
    }
  }
  dist[root] = 0;

  bool next_level = true;
  size_type l = 1;
  while (next_level)
  {
    next_level = false;

    #pragma mta for all streams stream_id of num_streams
    {
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++verts)
      {
        vertex_descriptor parent = bfs_parent[*verts];
        if (parent != null_vert && dist[*verts] == ULONG_MAX &&
            dist[parent] == l - 1)
        {
          dist[*verts] = l;
          next_level = true;
        }
      }
    }

    ++l;
  }

  unsigned long shortcuts = 0;
  unsigned long notreachable = 0;
  const unsigned long MAXSHORTCUT = 1000;

  struct Short {
    size_type u;
    size_type v;
  } shortcut[MAXSHORTCUT];

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  size_type num_blocks = (order + 128 - 1) / 128;

  #pragma mta assert parallel
  #pragma mta assert nodep
  #pragma mta no scalar expansion
  for (size_type block_id = 0; block_id < num_blocks; ++block_id)
  {
    size_type start_pos = begin_block_range(order, block_id, num_blocks);
    size_type end_pos = end_block_range(order, block_id, num_blocks);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);

    for (size_type i = start_pos; i != end_pos; ++i, ++verts)
    {
      vertex_descriptor u = *verts;
      size_type uid = get(vid_map, u);
      size_type uDist = dist[u];

      if (uDist != ULONG_MAX)
      {
        thread_adjacency_iterator adjs = thread_adjacent_vertices(u, 0, g);
        size_type odeg = out_degree(u, g);

        // If any of the edges is a shortcut, the tree is not
        // a shortest path tree.
        #pragma mta assert nodep
        for (size_type j = 0; j < odeg; ++j, ++adjs)
        {
          vertex_descriptor v = *adjs;
          size_type vid = get(vid_map, v);

          if (dist[v] == ULONG_MAX || uDist + 1 < dist[v])
          {
            unsigned long s = mt_incr(shortcuts, 1);
            if (s < MAXSHORTCUT)
            {
              shortcut[s].u = uid;
              shortcut[s].v = vid;
            }
          }
        }
      }
      else
      {
        if (bfs_parent[u] != null_vert) ++notreachable;
      }
    }
  }

  if (notreachable)
  {
    std::cout << notreachable << " nodes are marked, but not reachable from "
              << get(vid_map, root) << "." << std::endl;
  }

  if (shortcuts)
  {
    std::cout << "There are " << shortcuts << " shortcut edges." << std::endl;

    size_type minshort = (MAXSHORTCUT < shortcuts) ? MAXSHORTCUT : shortcuts;
    for (size_type i = 0; i < minshort; ++i)
    {
      size_type u = shortcut[i].u;
      size_type v = shortcut[i].v;
      std::cout << u << " " << v << " " << dist[u] << " " << dist[v] << " "
                << bfs_parent[u] << " " << bfs_parent[v] << std::endl;
    }
  }
}

}

#undef BFS_CHUNK

#endif
