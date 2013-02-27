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
/*! \file test_filtering.cpp

    \brief Tests filtering of graph vertices and edges.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 8/16/2010
*/
/****************************************************************************/

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/filter_adapter.hpp>
#include <mtgl/mtgl_test.hpp>

//#define BIDIR_GRAPH
//#define NO_FILTER
//#define ADJ_LIST

using namespace mtgl;

template <typename Graph>
class vertex_filter {
public:
  vertex_filter(const Graph& _g) : g(_g) {}

  bool operator()(typename graph_traits<Graph>::vertex_descriptor v) const
  { return get(get(_vertex_id_map, g), v) < 200; }

private:
  const Graph& g;
};

template <typename Graph>
class edge_filter {
public:
  edge_filter(const Graph& _g) : g(_g) {}

  bool operator()(typename graph_traits<Graph>::edge_descriptor e) const
  { return get(get(_edge_id_map, g), e) < 300; }

private:
  const Graph& g;
};

template <typename Graph>
void count_indegree(Graph& g,
                    typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type order = num_vertices(g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);

    for (size_type j = 0; j < end; ++j)
    {
      vertex_descriptor v = adjs[j];
      size_type vid = get(vid_map, v);

      mt_incr(indeg[vid], 1);
    }
  }
}

template <typename Graph, typename Filter>
void count_indegree_filtered(Graph& g,
                             typename graph_traits<Graph>::size_type* indeg,
                             const Filter& f)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type order = num_vertices(g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = verts[i];

    if (f(u))
    {
      adjacency_iterator adjs = adjacent_vertices(u, g);
      size_type end = out_degree(u, g);

      for (size_type j = 0; j < end; ++j)
      {
        vertex_descriptor v = adjs[j];

        if (f(v))
        {
          size_type vid = get(vid_map, v);
          mt_incr(indeg[vid], 1);
        }
      }
    }
  }
}

template <typename Graph>
inline
void count_indegree_filtered(Graph& g,
                             typename graph_traits<Graph>::size_type* indeg)
{
  count_indegree_filtered(g, indeg, keep_all());
}

template <typename Graph, typename Vertex_Filter, typename Edge_Filter>
void count_indegree_filtered(Graph& g,
                             typename graph_traits<Graph>::size_type* indeg,
                             const Vertex_Filter& vf, const Edge_Filter& ef)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type size = num_edges(g);
  edge_iterator edgs = edges(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];

    if (ef(e) && vf(source(e, g)) && vf(target(e, g)))
    {
      vertex_descriptor v = target(e, g);
      size_type vid = get(vid_map, v);
      mt_incr(indeg[vid], 1);
    }
  }
}

template <typename Graph>
void count_indegree_graph_filter(Graph& g,
                                 typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type order = num_vertices(g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    if (is_valid(verts, i, g))
    {
      vertex_descriptor u = verts[i];
      adjacency_iterator adjs = adjacent_vertices(u, g);
      size_type end = out_degree(u, g);

      for (size_type j = 0; j < end; ++j)
      {
        if (is_valid(adjs, j, g))
        {
          vertex_descriptor v = adjs[j];
          size_type vid = get(vid_map, v);

          mt_incr(indeg[vid], 1);
        }
      }
    }
  }
}

template <typename Graph>
void count_indegree_edge(Graph& g,
                         typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type size = num_edges(g);
  edge_iterator edgs = edges(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];
    vertex_descriptor v = target(e, g);
    size_type vid = get(vid_map, v);

    mt_incr(indeg[vid], 1);
  }
}

template <typename Graph, typename Filter>
void count_indegree_edge_filtered(
    Graph& g, typename graph_traits<Graph>::size_type* indeg, const Filter& f)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type size = num_edges(g);
  edge_iterator edgs = edges(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];

    if (f(e))
    {
      vertex_descriptor v = target(e, g);
      size_type vid = get(vid_map, v);
      mt_incr(indeg[vid], 1);
    }
  }
}

template <typename Graph>
inline
void count_indegree_edge_filtered(
    Graph& g, typename graph_traits<Graph>::size_type* indeg)
{
  count_indegree_edge_filtered(g, indeg, keep_all());
}

template <typename Graph>
void count_indegree_edge_graph_filter(
    Graph& g, typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type size = num_edges(g);
  edge_iterator edgs = edges(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < size; ++i)
  {
    if (is_valid(edgs, i, g))
    {
      edge_descriptor e = edgs[i];
      vertex_descriptor v = target(e, g);
      size_type vid = get(vid_map, v);
      mt_incr(indeg[vid], 1);
    }
  }
}

template <typename Graph>
void count_indegree_graph_outedge_filter(
    Graph& g, typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type order = num_vertices(g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    if (is_valid(verts, i, g))
    {
      vertex_descriptor u = verts[i];
      out_edge_iterator oedges = out_edges(u, g);
      size_type end = out_degree(u, g);

      for (size_type j = 0; j < end; ++j)
      {
        if (is_valid(oedges, j, g))
        {
          edge_descriptor e = oedges[j];
          size_type vid = get(vid_map, target(e, g));

          mt_incr(indeg[vid], 1);
        }
      }
    }
  }
}

#ifdef BIDIR_GRAPH
template <typename Graph>
void count_indegree_graph_inedge_filter(
    Graph& g, typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::in_edge_iterator in_edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type order = num_vertices(g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    if (is_valid(verts, i, g))
    {
      vertex_descriptor u = verts[i];
      size_type uid = get(vid_map, u);
      in_edge_iterator iedges = in_edges(u, g);
      size_type end = in_degree(u, g);

      for (size_type j = 0; j < end; ++j)
      {
        if (is_valid(iedges, j, g))
        {
          edge_descriptor e = iedges[j];

          mt_incr(indeg[uid], 1);
        }
      }
    }
  }
}
#endif

int main(int argc, char* argv[])
{
#ifdef BIDIR_GRAPH
#ifdef ADJ_LIST
  typedef adjacency_list<bidirectionalS> Graph;
#else
  typedef compressed_sparse_row_graph<bidirectionalS> Graph;
#endif
#else
#ifdef ADJ_LIST
  typedef adjacency_list<directedS> Graph;
#else
  typedef compressed_sparse_row_graph<directedS> Graph;
#endif
#endif
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;

  init_test(argc, argv);

  Graph g;
  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  printf("order: %lu, size: %lu\n", order, size);

  filter_adapter<Graph, keep_all, vertex_filter<Graph> >
    fg_v(g, keep_all(), vertex_filter<Graph>(g));

  filter_adapter<Graph, edge_filter<Graph>, keep_all>
    fg_e(g, edge_filter<Graph>(g), keep_all());

  filter_adapter<Graph, edge_filter<Graph>, vertex_filter<Graph> >
    fg_ve(g, edge_filter<Graph>(g), vertex_filter<Graph>(g));

  filter_adapter<Graph, keep_all, keep_all> fg(g, keep_all(), keep_all());

  size_type* indeg = new size_type[order];

  mt_timer timer;
  int issues, memrefs, concur, streams, traps;

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  int phantoms = mta_get_task_counter(RT_PHANTOM);
  int ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree(g, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  size_type indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("                     count_indegree():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_edge(g, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("                count_indegree_edge():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_filtered(g, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("            count_indegree_filtered():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_edge_filtered(g, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("       count_indegree_edge_filtered():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

#ifndef NO_FILTER
  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_filtered(g, indeg, vertex_filter<Graph>(g));

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("            count_indegree_filtered():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

#ifdef NO_FILTER
  count_indegree_graph_filter(g, indeg);
#else
  count_indegree_graph_filter(fg_v, indeg);
#endif

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("        count_indegree_graph_filter():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

#ifdef NO_FILTER
  count_indegree_graph_outedge_filter(g, indeg);
#else
  count_indegree_graph_outedge_filter(fg_v, indeg);
#endif

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("count_indegree_graph_outedge_filter():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

#ifdef BIDIR_GRAPH
  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

#ifdef NO_FILTER
  count_indegree_graph_inedge_filter(g, indeg);
#else
  count_indegree_graph_inedge_filter(fg_v, indeg);
#endif

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf(" count_indegree_graph_inedge_filter():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif
#endif

#ifndef NO_FILTER
  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_edge_filtered(g, indeg, edge_filter<Graph>(g));

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("       count_indegree_edge_filtered():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

#ifdef NO_FILTER
  count_indegree_edge_graph_filter(g, indeg);
#else
  count_indegree_edge_graph_filter(fg_e, indeg);
#endif

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("   count_indegree_edge_graph_filter():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

#ifndef NO_FILTER
  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_filtered(g, indeg, vertex_filter<Graph>(g),
                          edge_filter<Graph>(g));

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("            count_indegree_filtered():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_graph_outedge_filter(fg_ve, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf("count_indegree_graph_outedge_filter():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

#ifdef BIDIR_GRAPH
  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_indegree_graph_inedge_filter(fg_ve, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  indeg_sum = 0;
  for (size_type i = 0; i < order; ++i) indeg_sum += indeg[i];

  printf("------------------------------------------------\n");
  printf(" count_indegree_graph_inedge_filter():%8lu\n\n", indeg_sum);
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif
#endif
#endif

  printf("------------------------------------------------\n");

  delete [] indeg;

  return 0;
}
