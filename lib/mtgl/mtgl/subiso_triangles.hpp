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
/*! \file subiso_triangles.hpp

    \brief Searches for a certain bipartite structure in a graph.  This
     provides an example for using subgraph_isomorphism.

    \author Jon Berry (jberry@sandia.gov)

    \date 9/2010
*/
/****************************************************************************/

#ifndef MTGL_FIND_BIPARTITE_HPP
#define MTGL_FIND_BIPARTITE_HPP

#include <mtgl/subgraph_isomorphism.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/dynamic_array.hpp>

namespace mtgl {

// Example "don't care" comparator to allow WILD CARD conditions.
// It is hard coded to allow the second vertex in the walk to be any type.
template <typename Graph, typename VertexPropertyMap,
          typename EdgePropertyMap, typename TargetGraph = Graph,
          typename TargetVertexPropertyMap = VertexPropertyMap,
          typename TargetEdgePropertyMap = EdgePropertyMap>
class subtri_dontcare_comparator {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex;
  typedef typename graph_traits<Graph>::edge_descriptor edge;
  typedef typename graph_traits<TargetGraph>::vertex_descriptor vertex_trg;
  typedef typename graph_traits<TargetGraph>::edge_descriptor edge_trg;

  subtri_dontcare_comparator(VertexPropertyMap& gvt, EdgePropertyMap& get,
                             TargetVertexPropertyMap& trgvt,
                             TargetEdgePropertyMap& trget) :
    g_vtype(gvt), g_etype(get), trg_vtype(trgvt), trg_etype(trget) {}

  bool operator()(const edge& g_e, const edge_trg& trg_e,
                  Graph& g, TargetGraph& trg, size_type walk_level)
  {
#ifdef DOUBLE_DEBUG
    printf("si_comp: level: %lu\n", walk_level);
#endif

    // Get vertices.
    vertex g_src = source(g_e, g);
    vertex g_dest = target(g_e, g);
    vertex_trg trg_src = source(trg_e, trg);
    vertex_trg trg_dest = target(trg_e, trg);

    // We know that the triangle walk goes 0->1->2->0, with all edges
    // pointed in the same direction.  In this case, the second vertex
    // (Vertex 1) will be the target of the walk edge at Level 0 and
    // the source of the walk edge at Level 1.  In those, cases, omit
    // the vertex type check (invoke the wild card).  In the other case,
    // match all types strictly as the default comparator would.
    if (walk_level == 0)
    {
      return (g_vtype[g_src] == trg_vtype[trg_src] &&
              g_etype[g_e] == trg_etype[trg_e] &&
              out_degree(g_src, g) >= out_degree(trg_src, trg) &&
              out_degree(g_dest, g) >= out_degree(trg_dest, trg));
    }
    else if (walk_level == 1)
    {
      return (g_etype[g_e] == trg_etype[trg_e] &&
              g_vtype[g_dest] == trg_vtype[trg_dest] &&
              out_degree(g_src, g) >= out_degree(trg_src, trg) &&
              out_degree(g_dest, g) >= out_degree(trg_dest, trg));
    }
    else
    {
      return (g_vtype[g_src] == trg_vtype[trg_src] &&
              g_etype[g_e] == trg_etype[trg_e] &&
              g_vtype[g_dest] == trg_vtype[trg_dest] &&
              out_degree(g_src, g) >= out_degree(trg_src, trg) &&
              out_degree(g_dest, g) >= out_degree(trg_dest, trg));
    }
  }

private:
  VertexPropertyMap& g_vtype;
  EdgePropertyMap& g_etype;
  TargetVertexPropertyMap& trg_vtype;
  TargetEdgePropertyMap& trg_etype;
};

template<typename Graph>
class default_path_printer {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  default_path_printer(Graph& gg) : g(gg) {}

  void operator()(dynamic_array<size_type>& path)
  {
    printf("default_path_visitor; path size: %lu\n", path.size());

    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> eid_map = get(_edge_id_map, g);
    edge_iterator ei = edges(g);
    vertex_iterator vi = vertices(g);

    for (size_type i = 0; i < path.size(); ++i)
    {
      if (i % 2 == 0)
      {
        size_type vid = path[i];
        printf("vertex [ %lu  ]\n", vid);
      }
      else
      {
        size_type eid = path[i];
        edge_descriptor e = ei[eid];
        size_type uid = get(vid_map, source(e, g));
        size_type vid = get(vid_map, target(e, g));

        printf("edge ( %lu , %lu )\n", uid, vid);
      }
    }
  }

private:
  Graph& g;
};

template <typename Graph, typename Visitor, typename VertexPropertyMap>
class tri_match_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  tri_match_visitor(Graph& gg, Visitor uv, VertexPropertyMap& vp) :
    g(gg), vid_map(get(_vertex_id_map, g)),
    eid_map(get(_edge_id_map, g)), user_vis(uv), vpm(vp) {}

  void operator()(dynamic_array<vertex_descriptor>& path_verts,
                  dynamic_array<edge_descriptor>& path_edges)
  {

    vertex_iterator vi = vertices(g);
    edge_iterator ei = edges(g);

    size_type v0 = get(vid_map, path_verts[0]);
    size_type v1 = get(vid_map, path_verts[1]);
    size_type v2 = get(vid_map, path_verts[2]);
    size_type v3 = get(vid_map, path_verts[3]);

    if (v0 != v1 && v1 != v2 && v0 != v2 && v3 == v0)
    {
      user_vis(path_verts, path_edges, "subiso");
    }
  }

  void operator()(size_type* path, size_type length)
  {
#ifdef DEBUG
    assert(length == 7);
#endif

    vertex_iterator vi = vertices(g);
    edge_iterator ei = edges(g);

    size_type v0 = path[0];
    size_type v1 = path[2];
    size_type v2 = path[4];
    size_type v3 = path[6];

    if (v0 != v1 && v1 != v2 && v0 != v2 && v3 == v0)
    {
      printf("calling user vis\n");

      user_vis(path, length, "subiso");
    }
  }

  void operator()(dynamic_array<size_type>& path)
  {
#ifdef DEBUG
    assert(path.size() == 7);
#endif

    vertex_iterator vi = vertices(g);
    edge_iterator ei = edges(g);

    size_type v0 = path[0];
    size_type v1 = path[2];
    size_type v2 = path[4];
    size_type v3 = path[6];

    if (v0 != v1 && v1 != v2 && v0 != v2 && v3 == v0)
    {
      user_vis(path, "subiso");
    }
  }

private:
  Graph& g;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
  Visitor user_vis;
  VertexPropertyMap& vpm;
};

template <typename Graph, typename VertexProperty, typename EdgeProperty,
          typename Visitor>
typename graph_traits<Graph>::size_type
subiso_triangles(Graph& g, vertex_property_map<Graph, VertexProperty>& vPropmap,
                 edge_property_map<Graph, EdgeProperty>& ePropmap,
                 Visitor& user_vis)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef vertex_property_map<Graph, VertexProperty> vertex_property_map_t;
  typedef edge_property_map<Graph, EdgeProperty> edge_property_map_t;
  typedef pair<size_type, size_type> pair_t;

  typedef Graph TargetGraph;
  typedef size_type size_type_trg;
  typedef array_property_map<size_type, vertex_id_map<TargetGraph> >
          vertex_property_map_trg;
  typedef array_property_map<size_type, edge_id_map<TargetGraph> >
          edge_property_map_trg;

  // ************************************************************
  // ** finding bipartite structure *****************************
  // ************************************************************
  TargetGraph target;
  const size_type_trg numVertsTarget = 3;
  const size_type_trg numEdgesTarget = 3;

  size_type_trg sourcesTarget[numEdgesTarget] = { 0, 1, 2 };
  size_type_trg targetsTarget[numEdgesTarget] = { 1, 2, 0 };
  size_type vTypesTarget[numVertsTarget] = { 0, 1, 2 };
  size_type eTypesTarget[numEdgesTarget] = { 0, 0, 0 };

  init(numVertsTarget, numEdgesTarget, sourcesTarget, targetsTarget, target);

#ifdef DEBUG
  printf("TARGET:\n");
  print(target);
#endif

  // ************************************************************

  vertex_id_map<TargetGraph> trg_vid_map = get(_vertex_id_map, target);
  edge_id_map<TargetGraph> trg_eid_map = get(_edge_id_map, target);
  vertex_property_map_trg vTypemapTarget(vTypesTarget, trg_vid_map);
  edge_property_map_trg eTypemapTarget(eTypesTarget, trg_eid_map);

  // *************************************************************
  // SET UP THE TARGET WALK
  // *************************************************************
  dynamic_array<vertex_descriptor> walk_verts;
  walk_verts.push_back(vertices(target)[0]);
  walk_verts.push_back(vertices(target)[1]);
  walk_verts.push_back(vertices(target)[2]);
  walk_verts.push_back(vertices(target)[0]);

  dynamic_array<edge_descriptor> walk_edges;
  walk_edges.push_back(edges(target)[0]);
  walk_edges.push_back(edges(target)[1]);
  walk_edges.push_back(edges(target)[2]);

  // *************************************************************

  typedef default_path_printer<Graph> path_printer;
  typedef tri_match_visitor<Graph, Visitor, vertex_property_map_t>
          path_evaluator;

  // **********************************************************************
  // WILD CARDS:  for an example of allowing "wild cards" (don't care
  //              conditions), comment in the dontcare_comparator to
  //              replace si_default_comparator.
  // **********************************************************************
//  typedef subtri_dontcare_comparator<Graph, vertex_property_map_t,
//                                     edge_property_map_t, TargetGraph,
//                                     vertex_property_map_trg,
//                                     edge_property_map_trg>
//          si_comparator;
  typedef si_default_comparator<Graph, vertex_property_map_t,
                                edge_property_map_t, TargetGraph,
                                vertex_property_map_trg,
                                edge_property_map_trg>
          si_comparator;

  si_comparator si_comp(vPropmap, ePropmap, vTypemapTarget, eTypemapTarget);

  path_printer si_print(g);
  path_evaluator si_visit(g, user_vis, vPropmap);

/*
  dynamic_array<pair_t> level_pairs;
  level_pairs.push_back(pair_t(0, 3));
  size_type num_candidates =
    subgraph_isomorphism(g, target, walk_verts, walk_edges, si_comp,
                         si_visit, level_pairs);
*/

  // The last argument forces vertex 0 in the walk to be the same as vertex 3
  // in the walk.  (Constraining the candidate matches to be closed walks in
  // this case.)

  size_type num_candidates = subgraph_isomorphism(g, target, si_comp, si_visit);

  return num_candidates;
}

}

#endif
