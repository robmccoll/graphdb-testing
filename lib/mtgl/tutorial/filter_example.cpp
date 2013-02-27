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

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/filter_adapter.hpp>

using namespace mtgl;

template <typename Graph>
class vertex_filter {
public:
  vertex_filter(const Graph& _g) : g(_g) {}

  bool operator()(typename graph_traits<Graph>::vertex_descriptor v) const
  {
    return get(get(_vertex_id_map, g), v) > 0 &&
           get(get(_vertex_id_map, g), v) < 5;
  }

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

template <typename Graph>
void count_indegree_filter(Graph& g,
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
void count_indegree_filter_graph(Graph& g,
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

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  typedef filter_adapter<Graph, keep_all, vertex_filter<Graph> > FilterGraph;

  const size_type numVerts = 6;
  const size_type numEdges = 8;

  size_type sources[numEdges] = { 0, 0, 1, 1, 2, 3, 3, 4 };
  size_type targets[numEdges] = { 1, 2, 2, 3, 4, 4, 5, 5 };

  // Initialize the graph.
  Graph g;
  init(numVerts, numEdges, sources, targets, g);

  FilterGraph fg_v(g, keep_all(), vertex_filter<Graph>(g));

  size_type order = num_vertices(g);

  size_type* indeg = (size_type*) malloc(order * sizeof(size_type));
  for (int i = 0; i < order; ++i) indeg[i] = 0;

  count_indegree(g, indeg);

  size_type indeg_sum = 0;
  for (int i = 0; i < order; ++i) indeg_sum += indeg[i];
  printf("              indegree sum full graph:%8lu\n", indeg_sum);

  for (int i = 0; i < order; ++i) indeg[i] = 0;

  count_indegree_filter(g, indeg);

  indeg_sum = 0;
  for (int i = 0; i < order; ++i) indeg_sum += indeg[i];
  printf("  indegree sum full graph filter code:%8lu\n", indeg_sum);

  for (int i = 0; i < order; ++i) indeg[i] = 0;

  count_indegree_filter_graph(fg_v, indeg);

  indeg_sum = 0;
  for (int i = 0; i < order; ++i) indeg_sum += indeg[i];
  printf("indegree sum filter graph filter code:%8lu\n", indeg_sum);

  free(indeg);

  return 0;
}
