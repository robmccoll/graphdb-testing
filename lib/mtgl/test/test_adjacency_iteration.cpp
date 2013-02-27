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
/*! \file test_adjacency_iteration.cpp

    \brief Compares visiting the adjacencies of a graph by accessing a CSR
           graph structure directly with using the MTGL interface with using
           the manually load-balanced version of visit_adj().

    \author Jon Berry (jberry@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    /date 8/15/2009

    The inlined function test_pred_func() below represents a typical predicate
    that users may wish to introduce as an algorithmic filter.  Unfortunately,
    the use of this predicate via a function pointer in the inner loop of the
    function count_adjacencies_higher_id_func_ptr() prevents the loop-merge
    necessary for good performance.  The result is a serial inner loop, which
    has severe consequences when processing power-law data.
     
    The predicate is easily introduced into the function object
    test_pred_func_obj.  Using this function object instead of the function
    pointer restores the loop merge as demonstrated in the function
    count_adjacencies_higher_id_func_obj().
*/
/****************************************************************************/

#include <cstdlib>
#include <climits>

#include <mtgl/visit_adj.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/random.hpp>

using namespace mtgl;

typedef compressed_sparse_row_graph<directedS> Graph;
typedef graph_traits<Graph>::size_type size_type;
typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<Graph>::vertex_iterator vertex_iterator;

/// Counts the number of adjacencies with a higher id than the source using a
/// pure C traveral.
template <typename Graph>
void
count_adjacencies_higher_id(Graph& g,
                            typename graph_traits<Graph>::size_type* indeg)
{
  typedef typename graph_traits<Graph>::size_type size_type;

  size_type order = g.get_order();
  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    size_type begin = index[i];
    size_type end = index[i + 1];

    for (size_type j = begin; j < end; ++j)
    {
      if (i < end_points[j]) mt_incr(indeg[end_points[j]], 1);
    }
  }
}

typedef bool (*pred_t)(size_type, size_type);

#pragma mta inline
inline bool test_pred_func(size_type i, size_type j)
{
  return (i < j);
}

/// Counts the number of adjacencies with a higher id than the source using a
/// pure C traveral with a function pointer.
template <typename Graph>
void
count_adjacencies_higher_id_func_ptr(
    Graph& g, typename graph_traits<Graph>::size_type* indeg, pred_t my_func)
{
  typedef typename graph_traits<Graph>::size_type size_type;

  size_type order = g.get_order();
  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  #pragma mta assert noalias *indeg
  #pragma mta assert parallel
  for (size_type i = 0; i < order; ++i)
  {
    size_type begin = index[i];
    size_type end = index[i + 1];

    #pragma mta assert parallel
    for (size_type j = begin; j < end; ++j)
    {
      if (my_func(i, end_points[j])) mt_incr(indeg[end_points[j]], 1);
    }
  }
}

template <typename Graph>
class test_pred_func_obj {
public:
  typedef typename graph_traits<Graph>::size_type size_type;

  inline bool operator()(size_type i, size_type j) { return (i < j); }
};

/// Counts the number of adjacencies with a higher id than the source using a
/// pure C traveral with a function object.
template <typename Graph, typename pred>
void
count_adjacencies_higher_id_func_obj(
    Graph& g, typename graph_traits<Graph>::size_type* indeg, pred my_func)
{
  typedef typename graph_traits<Graph>::size_type size_type;

  size_type order = g.get_order();
  size_type* index = g.get_index();
  size_type* end_points = g.get_end_points();

  #pragma mta assert noalias *indeg
  for (size_type i = 0; i < order; ++i)
  {
    size_type begin = index[i];
    size_type end = index[i + 1];

    for (size_type j = begin; j < end; ++j)
    {
      if (my_func(i, end_points[j])) mt_incr(indeg[end_points[j]], 1);
    }
  }
}

/// Counts the number of adjacencies with a higher id than the source using an
/// MTGL traveral.
template <typename Graph, typename predicate>
void count_adjacencies_higher_id_mtgl(
    Graph& g, typename graph_traits<Graph>::size_type* indeg, predicate f)
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
    size_type uid = get(vid_map, u);
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);

    for (size_type j = 0; j < end; ++j)
    {
      vertex_descriptor v = adjs[j];
      size_type vid = get(vid_map, v);

      if (f(uid, vid)) mt_incr(indeg[vid], 1);
    }
  }
}

/// A visitor to the load-balanced visit_adj function.
template <typename Graph, typename predicate>
class adjacencies_higher_id_computer {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  adjacencies_higher_id_computer(vertex_id_map<Graph> vm, size_type *ind,
                                 predicate ff) :
    vid_map(vm), in_degree(ind), f(ff) {}

  inline void operator()(vertex_descriptor src, vertex_descriptor dest)
  {
    size_type v1id = get(vid_map, src);
    size_type v2id = get(vid_map, dest);

    if (f(v1id, v2id)) mt_incr(in_degree[v2id], 1);
  }

  void post_visit() {}

private:
  vertex_id_map<Graph>& vid_map;
  size_type* in_degree;
  predicate f;
};

template <typename Graph, typename visitor>
void visit_adj_partial(Graph& g, visitor f)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  const size_type* index = g.get_index();
  const vertex_descriptor* end_points = g.get_end_points();
  const size_type order = num_vertices(g);

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert parallel
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = verts[i];
    size_type begin = index[get(vid_map, u)];
    size_type end = index[get(vid_map, u) + 1];
    #pragma mta assert parallel
    for (size_type j = begin; j < end; ++j)
    {
      vertex_descriptor v = verts[end_points[j]];
      f(u, v);
    }
  }
}

template <typename Graph, typename visitor>
void visit_adj_full(Graph& g, visitor f)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  const size_type *index = g.get_index();
  const vertex_descriptor *end_points = g.get_end_points();
  const size_type order = num_vertices(g);

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  vertex_iterator verts = vertices(g);

  #pragma mta assert parallel
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = verts[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type end = out_degree(u, g);
    #pragma mta assert parallel
    for (size_type j = 0; j < end; ++j)
    {
      vertex_descriptor v = adjs[j];
      f(u, v);
    }
  }
}

template <typename Graph, typename visitor>
void visit_adj_partial(Graph& g,
         typename graph_traits<Graph>::vertex_descriptor* to_visit,
         typename graph_traits<Graph>::size_type num_to_visit,
         visitor f)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  const size_type *index = g.get_index();
  const vertex_descriptor *end_points = g.get_end_points();

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  #pragma mta assert parallel
  for (size_type i = 0; i < num_to_visit; ++i)
  {
    vertex_descriptor u = to_visit[i];
    size_type begin = index[get(vid_map, u)];
    size_type end = index[get(vid_map, u) + 1];
    #pragma mta assert parallel
    for (size_type j = begin; j < end; ++j)
    {
      f(end_points[i], end_points[j]);
    }
  }
}

template <typename Graph, typename visitor>
void visit_adj(Graph& g,
         typename graph_traits<Graph>::vertex_descriptor* to_visit,
         typename graph_traits<Graph>::size_type num_to_visit,
         visitor f)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  const size_type *index = g.get_index();
  const vertex_descriptor *end_points = g.get_end_points();
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  #pragma mta assert parallel
  for (size_type i = 0; i < num_to_visit; ++i)
  {
    vertex_descriptor u = to_visit[i];
    adjacency_iterator adjs = adjacent_vertices(u, g);
    size_type deg = out_degree(u, g);

    #pragma mta assert parallel
    for (size_type j = 0; j < deg; ++j)
    {
      vertex_descriptor v = adjs[j];
      f(u, j);
    }
  }
}

void checkError(size_type* indeg, size_type* indeg2, size_type order)
{
  size_type error = (std::numeric_limits<size_type>::max)();

  for (size_type i = 0; i < order; ++i)
  {
    if (indeg[i] != indeg2[i]) error = i;
  }

  if (error != (std::numeric_limits<size_type>::max)())
  {
    std::cout << "Error in computation: pos " << error << std::endl;
  }
}

int main(int argc, char* argv[])
{
  mt_srand48(0);

  init_test(argc, argv);

  Graph g;
  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  size_type* indeg = new size_type[order];
  size_type* indeg2 = new size_type[order];

  test_pred_func_obj<Graph> tpc;

  printf("order: %lu, size: %lu\n", order, size);

  for (size_type i = 0; i < order; ++i) indeg[i] = 0;

  mt_timer timer;
  int issues, memrefs, concur, streams, traps;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  int phantoms = mta_get_task_counter(RT_PHANTOM);
  int ready = mta_get_task_counter(RT_READY);
#endif

  count_adjacencies_higher_id(g, indeg);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  printf("---------------------------------------------\n");
  printf("count_adjacencies_higher_id(): \n");
  printf("---------------------------------------------\n");
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);
#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_adjacencies_higher_id_func_ptr(g, indeg2, test_pred_func);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  printf("---------------------------------------------\n");
  printf("count_adjacencies_higher_id_func_ptr(): \n");
  printf("---------------------------------------------\n");
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_adjacencies_higher_id_func_obj(g, indeg2, tpc);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  printf("---------------------------------------------\n");
  printf("count_adjacencies_higher_id_func_obj(): \n");
  printf("---------------------------------------------\n");
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  count_adjacencies_higher_id_mtgl(g, indeg2, tpc);

  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  printf("---------------------------------------------\n");
  printf("count_adjacencies_higher_id_mtgl(): \n");
  printf("---------------------------------------------\n");
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  checkError(indeg, indeg2, order);

  for (size_type i = 0; i < order; ++i) indeg2[i] = 0;

  adjacencies_higher_id_computer<Graph, test_pred_func_obj<Graph> >
    idc(get(_vertex_id_map, g), indeg2, tpc);

  size_type* prefix_counts = 0;
  size_type* started_nodes = 0;
  size_type num_threads;

  init_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM);
  ready = mta_get_task_counter(RT_READY);
#endif

  visit_adj(g, idc, prefix_counts, started_nodes, num_threads);
  sample_mta_counters(timer, issues, memrefs, concur, streams, traps);

#ifdef __MTA__
  phantoms = mta_get_task_counter(RT_PHANTOM) - phantoms;
  ready = mta_get_task_counter(RT_READY) - ready;
#endif

  printf("---------------------------------------------\n");
  printf("visit_adj(): \n");
  printf("---------------------------------------------\n");
  print_mta_counters(timer, num_edges(g), issues, memrefs, concur, streams,
                     traps);

  free(prefix_counts);
  free(started_nodes);

#ifdef __MTA__
  printf("phantoms: %d\n", phantoms);
  printf("ready: %d\n", ready);
#endif

  checkError(indeg, indeg2, order);

  delete [] indeg;
  delete [] indeg2;

  return 0;
}
