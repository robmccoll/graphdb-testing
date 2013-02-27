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

#include <mtgl/subgraph_isomorphism.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/filter_graph.hpp>
#include <mtgl/subgraph_adapter.hpp>

using namespace mtgl;

template <typename Graph, typename vertex_property_map,
          typename edge_property_map>
void print_my_graph(Graph& g, vertex_property_map& vTypes,
                    edge_property_map& eTypes)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> eid_map = get(_edge_id_map, g);

  edge_iterator edgs = edges(g);
  size_type size = num_edges(g);
  for (size_type i = 0; i < size; i++)
  {
    edge_descriptor e = edgs[i];
    size_type eid = get(eid_map, e);

    vertex_descriptor u = source(e, g);
    size_type uid = get(vid_map, u);

    vertex_descriptor v = target(e, g);
    size_type vid = get(vid_map, v);

    printf("%lu: (%lu, %lu)   [%d, %d, %d]\n", eid, uid, vid, vTypes[u],
           eTypes[e], vTypes[v]);
  }
}

template <typename Graph, typename Subgraph>
class match_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Subgraph>::size_type size_type_sg;
  typedef typename graph_traits<Subgraph>::edge_iterator edge_iterator_sg;

  match_visitor(Graph& gg, Subgraph& sgg) : g(gg), sg(sgg) {}

  void operator()(dynamic_array<size_type_sg>& path)
  {
    typedef typename dynamic_array<size_type_sg>::size_type d_size_type;

    edge_id_map<Graph> eid_map = get(_edge_id_map, g);
    edge_iterator_sg edgs_sg = edges(sg);

    printf("Path:");
    for (d_size_type i = 1; i < path.size(); i += 2)
    {
      size_type eid = get(eid_map, sg.local_to_global(edgs_sg[path[i]]));
      printf("%6lu", eid);
    }
    printf("\n");
  }

private:
  Graph& g;
  Subgraph& sg;
};

int main(int argc, char *argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  // Initialize the original graph.
  Graph g;

  const size_type numVerts = 6;
  const size_type numEdges = 7;

  size_type sources[numEdges] = { 1, 2, 1, 2, 3, 3, 4 };
  size_type targets[numEdges] = { 0, 0, 2, 3, 4, 5, 5 };

  init(numVerts, numEdges, sources, targets, g);

  // Create the big graph vertex and edge type property maps.
  typedef array_property_map<int, vertex_id_map<Graph> > vertex_property_map;
  typedef array_property_map<int, edge_id_map<Graph> > edge_property_map;

  int vTypes[numVerts] = { 1, 2, 3, 2, 3, 1 };
  int eTypes[numEdges] = { 1, 2, 3, 3, 3, 1, 2 };

  vertex_id_map<Graph> g_vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> g_eid_map = get(_edge_id_map, g);

  vertex_property_map vTypemap(vTypes, g_vid_map);
  edge_property_map eTypemap(eTypes, g_eid_map);

  // Print the original graph with its vertex and edge types.
  printf("Big graph (%lu, %lu)\n", num_vertices(g), num_edges(g));
  print_my_graph(g, vTypemap, eTypemap);

  // Initialize target graph.
  Graph target;

  const size_type numVertsTarget = 3;
  const size_type numEdgesTarget = 3;

  size_type sourcesTarget[numEdgesTarget] = { 1, 2, 1 };
  size_type targetsTarget[numEdgesTarget] = { 0, 0, 2 };

  init(numVertsTarget, numEdgesTarget, sourcesTarget, targetsTarget, target);

  // Create the target graph vertex and edge type property maps.
  int vTypesTarget[numVertsTarget] = { 1, 2, 3 };
  int eTypesTarget[numEdgesTarget]  = { 1, 2, 3 };

  vertex_id_map<Graph> trg_vid_map = get(_vertex_id_map, target);
  edge_id_map<Graph> trg_eid_map = get(_edge_id_map, target);

  vertex_property_map vTypemapTarget(vTypesTarget, trg_vid_map);
  edge_property_map eTypemapTarget(eTypesTarget, trg_eid_map);

  // Print the target graph with its vertex and edge types.
  printf("\nTarget graph (%lu, %lu)\n", numVertsTarget, numEdgesTarget);
  print_my_graph(target, vTypemapTarget, eTypemapTarget);
  printf("\n");

  // Set up walk through target graph.
  const size_type walk_length = 7;
  size_type walk_ids[walk_length] = { 0, 0, 1, 2, 2, 1, 0 };

  // Filter large graph to only edges matching edges in target graph.
  typedef subgraph_adapter<Graph> FiltGraph;
  typedef graph_traits<FiltGraph>::size_type size_type_filt;

  FiltGraph filteredG(g);

  f_default_comparator<Graph, vertex_property_map, edge_property_map>
    f_comp(vTypemap, eTypemap, vTypemapTarget, eTypemapTarget);

  filter_graph_by_edges(g, target, f_comp, filteredG);

  // Create the filtered graph vertex and edge type property maps.
  typedef array_property_map<int, vertex_id_map<FiltGraph> >
          vertex_property_map_filt;
  typedef array_property_map<int, edge_id_map<FiltGraph> >
          edge_property_map_filt;

  size_type_filt numVertsFiltered = num_vertices(filteredG);
  size_type_filt numEdgesFiltered = num_edges(filteredG);

  int* vTypesFiltered = new int[numVertsFiltered];
  int* eTypesFiltered = new int[numEdgesFiltered];

  vertex_id_map<FiltGraph> filteredG_vid_map = get(_vertex_id_map, filteredG);
  edge_id_map<FiltGraph> filteredG_eid_map = get(_edge_id_map, filteredG);

  vertex_property_map_filt vTypemapFiltered(vTypesFiltered, filteredG_vid_map);
  edge_property_map_filt eTypemapFiltered(eTypesFiltered, filteredG_eid_map);

  // Copy the vertex types for the filtered graph.
  graph_traits<FiltGraph>::vertex_iterator verts_filt = vertices(filteredG);
  for (size_type_filt i = 0; i < numVertsFiltered; ++i)
  {
    graph_traits<FiltGraph>::vertex_descriptor fg_v = verts_filt[i];

    graph_traits<Graph>::vertex_descriptor g_v =
      filteredG.local_to_global(fg_v);

    vTypemapFiltered[fg_v] = vTypemap[g_v];
  }

  // Copy the edge types for the filtered graph.
  graph_traits<FiltGraph>::edge_iterator edgs_filt = edges(filteredG);
  for (size_type_filt i = 0; i < numEdgesFiltered; ++i)
  {
    graph_traits<FiltGraph>::edge_descriptor fg_e = edgs_filt[i];

    graph_traits<Graph>::edge_descriptor g_e = filteredG.local_to_global(fg_e);

    eTypemapFiltered[fg_e] = eTypemap[g_e];
  }

  // Print the filtered graph with its vertex and edge types.
  printf("\nFiltered graph (%lu, %lu)\n", numVertsFiltered, numEdgesFiltered);
  print_my_graph(filteredG, vTypemapFiltered, eTypemapFiltered);
  printf("\n");

  // Declare the comparator to define the matching and the visitor to apply
  // to matches.
  si_default_comparator<FiltGraph, vertex_property_map_filt,
                        edge_property_map_filt,
                        Graph, vertex_property_map, edge_property_map>
    si_comp(vTypemapFiltered, eTypemapFiltered, vTypemapTarget, eTypemapTarget);

  match_visitor<Graph, FiltGraph> mvisit(g, filteredG);

  subgraph_isomorphism(filteredG, target, walk_ids, walk_length,
                       si_comp, mvisit);

  delete [] vTypesFiltered;
  delete [] eTypesFiltered;

  return 0;
}
