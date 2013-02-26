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
/*! \file subiso_5cycles.hpp

    \brief Searches for a certain bipartite structure in a graph.  This
     provides an example for using subgraph_isomorphism.

    \author Jon Berry (jberry@sandia.gov)

    \date 3/2011
*/
/****************************************************************************/

#ifndef MTGL_SUBISO_5CYCLES_HPP
#define MTGL_SUBISO_5CYCLES_HPP

#include <mtgl/subgraph_isomorphism.hpp>
#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>

#define absdiff(x, y)  ((x) > (y) ? (x) - (y) : (y) - (x))

namespace mtgl {

template <typename Graph>
class hash_edges {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef xmt_hash_table<size_type, size_type> hash_table_t;

  hash_edges(Graph& gg, hash_table_t& edgs) :
    g(gg), vid_map(get(_vertex_id_map, g)), eid_map(get(_edge_id_map, g)),
    the_edges(edgs), order(num_vertices(gg)) {}

  void operator()(const edge_descriptor& e)
  {
    vertex_descriptor v1 = source(e, g);
    vertex_descriptor v2 = target(e, g);
    int v1id = get(vid_map, v1);
    int v2id = get(vid_map, v2);
    int eid = get(eid_map, e);

    order_pair(v1id, v2id);

    int64_t key = v1id * order + v2id;
    the_edges.insert(key, eid);
  }

private:
  Graph& g;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
  hash_table_t& the_edges;
  uint64_t order;
};

template <typename Graph>
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

template <typename T>
bool correct_color_sequence(T cv1, T cv2)
{
  // return 1;          // JWB DEBUG
  bool way1 = (absdiff(cv1, cv2) == 1);
  bool way2 = (absdiff(cv1, cv2) == 4);

  return (way1 || way2);    // 5 types.
}

template <typename Graph, typename VertexPropertyMap, typename Visitor>
class si_5cycle_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef xmt_hash_table<size_type, size_type> hash_table_t;

  si_5cycle_visitor(Graph& gg, hash_table_t& edgs,
                    VertexPropertyMap& vp, Visitor uv) :
    g(gg), the_edges(edgs), vertex_prop(vp), user_vis(uv) {}

  void operator()(dynamic_array<vertex_descriptor>& path_verts,
                  dynamic_array<edge_descriptor>& path_edges)
  {
    vertex_id_map<Graph> g_vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> g_eid_map = get(_edge_id_map, g);
    size_type order = num_vertices(g);
    vertex_iterator vi = vertices(g);
    edge_iterator ei = edges(g);

    size_type path_vert_ids[6];
    path_vert_ids[0] = get(g_vid_map, path_verts[0]);
    path_vert_ids[1] = get(g_vid_map, path_verts[1]);
    path_vert_ids[2] = get(g_vid_map, path_verts[2]);
    path_vert_ids[3] = get(g_vid_map, path_verts[3]);
    path_vert_ids[4] = get(g_vid_map, path_verts[4]);
    path_vert_ids[5] = get(g_vid_map, path_verts[5]);

    xmt_hash_table<size_type, size_type> ht(4 * path_verts.size());

    // Quit because incorrect color sequence.
    if (!correct_color_sequence(vertex_prop[path_verts[0]],
                                vertex_prop[path_verts[1]]))
    {
#ifdef DEBUG
      printf("No match because of color sequence (%lu, %lu).\n",
             path_vert_ids[0], path_vert_ids[1]);
#endif
      return;
    }

    if (path_vert_ids[0] != path_vert_ids[5])
    {
#ifdef DEBUG
      printf("No match because start != end\n");
#endif
      return;
    }

    ht.insert(path_vert_ids[0], 0);
    ht.insert(path_vert_ids[1], 0);

#ifdef DEBUG
    printf("match: ht insert %d %d\n", path_vert_ids[0], path_vert_ids[1]);
    printf("match: path vertices size: %d\n", path_verts.size());
#endif

    for (size_type i = 2; i < path_verts.size(); ++i)
    {
      // Quit because multiple edge.
      if (path_verts[i] == path_verts[i - 2]) return;

      // Quit because incorrect color sequence.
      if (!correct_color_sequence(vertex_prop[path_verts[i - 1]],
                                  vertex_prop[path_verts[i]]))
      {
#ifdef DEBUG
        printf("No match because of color sequence (%lu, %lu).\n",
               path_vert_ids[i - 1], path_vert_ids[i]);
#endif
        return;
      }

      ht.insert(path_verts[i], 0);

#ifdef DEBUG
      printf("match: ht insert %d\n", path[i]);
#endif
    }

    if (ht.size() == 5)
    {
      // This is a 5-cycle with sequential types.  Quit if there are any
      // chords.
      for (size_type i = 0; i < path_verts.size() - 1; ++i)
      {
        for (size_type j = i + 2; j < path_verts.size() - 1 && j <= (i + 3);
             ++j)
        {
          vertex_descriptor u = path_verts[i];
          vertex_descriptor v = path_verts[j];
          size_type u_id = get(g_vid_map, u);
          size_type v_id = get(g_vid_map, v);

          order_pair(u_id, v_id);

          uint64_t uv_key = u_id * order + v_id;

#ifdef DEBUG
          printf("testing edge(%d,%d)\n", u_id, v_id);
#endif

          if (the_edges.member(uv_key))
          {
#ifdef DEBUG
            printf("No match because of chord(%d, %d)\n", u_id, v_id);

            for (size_type i = 0; i < path_verts.size(); ++i)
            {
              printf("bad path[%d]: %d\n", i, path_verts[i]);
            }
#endif

            return;
          }
        }
      }

#ifdef DEBUG
      for (size_type i = 0; i < path_verts.size(); ++i)
      {
        printf("good path[%d]: %d\n", i, path_vert_ids[i]);
      }
#endif

      user_vis(path_verts, path_edges);

#ifdef DOUBLE_DEBUG
      printf("MATCH FOUND\n");

      for (size_type i = 0; i < path_edges.size(); ++i)
      {
        edge_descriptor e = path_edges[i];
        size_type s_vid = get(g_vid_map, source(e, g));
        size_type t_vid = get(g_vid_map, target(e, g));
        edge_descriptor e_print = e;

        printf("%6lu (%6lu , %6lu)\n", get(g_eid_map, e_print) + 1,
               s_vid + 1, t_vid + 1);
      }
#endif
    }
#ifdef DEBUG
    else
    {
      printf("no match because of hash table size\n");
    }
#endif
  }

  void operator()(dynamic_array<size_type>& path)
  {
    vertex_id_map<Graph> g_vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> g_eid_map = get(_edge_id_map, g);
    size_type order = num_vertices(g);
    vertex_iterator vi = vertices(g);
    edge_iterator ei = edges(g);

    xmt_hash_table<size_type, size_type> ht(4 * path.size());

    // Quit because incorrect color sequence.
    if (!correct_color_sequence(vertex_prop[vi[path[0]]],
                                vertex_prop[vi[path[2]]]))
    {
#ifdef DEBUG
      printf("No match because of color sequence (%lu, %lu).\n",
             path[0], path[2]);
#endif

      return;
    }

    if (path[0] != path[10])
    {
#ifdef DEBUG
      printf("No match because start != end\n");
#endif

      return;
    }

    ht.insert(path[0], 0);
    ht.insert(path[2], 0);

#ifdef DEBUG
    printf("match: ht insert %d %d\n", path[0], path[2]);
    printf("match: path size: %d\n", path.size());
#endif

    for (size_type i = 4; i < path.size(); i += 2)
    {
      // Quit because multiple edge.
      if (path[i] == path[i - 4]) return;

      // Quit because incorrect color sequence.
      if (!correct_color_sequence(vertex_prop[vi[path[i - 2]]],
                                  vertex_prop[vi[path[i]]]))
      {
#ifdef DEBUG
        printf("No match because of color sequence (%lu, %lu).\n",
               path[i - 2], path[i]);
#endif

        return;
      }

      ht.insert(path[i], 0);

#ifdef DEBUG
      printf("match: ht insert %d\n", path[i]);
#endif
    }

    if (ht.size() == 5)
    {
      // This is a 5-cycle with sequential types.  Quit if there are any
      // chords.
      for (size_type i = 0; i < path.size() - 2; i += 2)
      {
        for (size_type j = i + 4; j < path.size() - 2 && j <= (i + 6); j += 2)
        {
          // for (size_type j = i + 4; j < path.size(); j+= 2)
          vertex_descriptor u = vi[path[i]];
          vertex_descriptor v = vi[path[j]];
          size_type u_id = get(g_vid_map, u);
          size_type v_id = get(g_vid_map, v);

          order_pair(u_id, v_id);

          uint64_t uv_key = u_id * order + v_id;

#ifdef DEBUG
          printf("testing edge(%d,%d)\n", u_id, v_id);
#endif

          if (the_edges.member(uv_key))
          {
#ifdef DEBUG
            printf("No match because of chord(%d, %d)\n", u_id, v_id);
            for (size_type i = 0; i < path.size(); ++i)
            {
              printf("bad path[%d]: %d\n", i, path[i]);
            }
#endif

            return;
          }
        }
      }

#ifdef DEBUG
      for (size_type i = 0; i < path.size(); ++i)
      {
        printf("good path[%d]: %d\n", i, path[i]);
      }
#endif

      user_vis(path);

#ifdef DOUBLE_DEBUG
      printf("MATCH FOUND\n");

      for (size_type i = 1; i < path.size(); i += 2)
      {
        edge_descriptor e = ei[path[i]];
        size_type s_vid = get(g_vid_map, source(e, g));
        size_type t_vid = get(g_vid_map, target(e, g));
        edge_descriptor e_print = e;

        printf("%6lu (%6lu , %6lu)\n", get(g_eid_map, e_print) + 1,
               s_vid + 1, t_vid + 1);
      }
#endif
    }
#ifdef DEBUG
    else
    {
      printf("no match because of hash table size\n");
    }
#endif
  }

private:
  Graph& g;
  hash_table_t& the_edges;
  VertexPropertyMap& vertex_prop;
  Visitor user_vis;
};

template <typename Graph, typename VertexPropertyMap, typename EdgePropertyMap>
class bowl_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;
  typedef xmt_hash_table<size_type, size_type> hash_table_t;

  bowl_visitor(Graph& gg, hash_table_t& edgs, VertexPropertyMap& vp,
               EdgePropertyMap& ep, size_type& nep) :
    g(gg), vid_map(get(_vertex_id_map, g)), eid_map(get(_edge_id_map, g)),
    the_edges(edgs), order(num_vertices(gg)), vertex_prop(vp), edge_prop(ep),
    num_edges_processed(nep)
  {
#ifdef DEBUG
    vertex_iterator verts = vertices(g);
    for (size_type i = 0; i < num_vertices(g); ++i)
    {
      vertex_descriptor v = verts[i];
      printf("bv:vprop[%lu]: %lu\n", get(vid_map, v), vertex_prop[v]);
    }
#endif
  }

  void operator()(const edge_descriptor& e)
  {
#ifdef DEBUG
    size_type cur_proc = mt_incr(num_edges_processed, 1);
    if (cur_proc % 1000 == 0) printf(".\n");
    fflush(stdout);
#endif

    vertex_descriptor v1 = source(e, g);
    vertex_descriptor v2 = target(e, g);

    // Check for correct sequence of types.
    size_type v1_color = vertex_prop[v1];
    size_type v2_color = vertex_prop[v2];
    if (!correct_color_sequence(v1_color, v2_color)) return;

    size_type v1_deg = out_degree(v1, g);
    size_type v2_deg = out_degree(v2, g);

    edge_prop[e] = 0;

    adjacency_iterator v1_neighs = adjacent_vertices(v1, g);
    adjacency_iterator v2_neighs = adjacent_vertices(v2, g);

    for (size_type i = 0; i < v1_deg; ++i)
    {
      vertex_descriptor v1_ngh = v1_neighs[i];

      if (v1_ngh != v2)
      {
        for (size_type j = 0; j < v2_deg; ++j)
        {
          vertex_descriptor v2_ngh = v2_neighs[j];

          if (v2_ngh != v1)
          {
            check_bowl(e, v1, v2, v1_ngh, v2_ngh, vertex_prop,
                       edge_prop, the_edges);
          }
        }
      }
    }
  }

  void check_bowl(const edge_descriptor& the_edge, const vertex_descriptor& v1,
                  const vertex_descriptor& v2, const vertex_descriptor& v1_ngh,
                  const vertex_descriptor& v2_ngh,
                  VertexPropertyMap& vertex_prop,
                  EdgePropertyMap& edge_prop, hash_table_t& the_edges)
  {
    vertex_descriptor a = v1;
    vertex_descriptor b = v2;
    vertex_descriptor c = v2_ngh;
    vertex_descriptor e = v1_ngh;

    // Check for correct sequence of types.
    size_type a_color = vertex_prop[a];
    size_type b_color = vertex_prop[b];
    size_type c_color = vertex_prop[c];
    size_type e_color = vertex_prop[e];

    if (!correct_color_sequence(b_color, c_color)) return;
    if (!correct_color_sequence(a_color, e_color)) return;

    // Check for absence of (b,e).
    size_type b_id = get(vid_map, b);
    size_type e_id = get(vid_map, e);
    order_pair(b_id, e_id);

    uint64_t be_key = b_id * order + e_id;

    if (the_edges.member(be_key)) return;

    // Check for absence of (a,c).
    size_type a_id = get(vid_map, a);
    size_type c_id = get(vid_map, c);
    order_pair(a_id, c_id);

    uint64_t ac_key = a_id * order + c_id;

    if (the_edges.member(ac_key)) return;

    // Check for absence of (c,e).
    size_type c_id2 = get(vid_map, c);
    size_type e_id2 = get(vid_map, e);
    order_pair(c_id2, e_id2);

    uint64_t ce_key = c_id2 * order + e_id2;

    if (the_edges.member(ce_key)) return;

    edge_prop[the_edge] = 1;

#ifdef DEBUG
    printf("epm(%lu, %lu) = 1\n", a_id, b_id);
#endif
  }

private:
  Graph& g;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
  hash_table_t& the_edges;
  uint64_t order;
  VertexPropertyMap& vertex_prop;
  EdgePropertyMap& edge_prop;
  size_type num_edges_processed;
};


template <typename Graph, typename vptype, typename eptype, typename Visitor>
typename graph_traits<Graph>::size_type
subiso_5cycles(Graph& g, vertex_property_map<Graph, vptype>& vPropmap,
               edge_property_map<Graph, eptype>& ePropmap,
               Visitor& user_vis)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef vertex_property_map<Graph, vptype> vertex_property_map_t;
  typedef edge_property_map<Graph, eptype> edge_property_map_t;
  typedef xmt_hash_table<size_type, size_type> hash_table_t;
  typedef pair<size_type, size_type> pair_t;
  typedef dynamic_array<pair_t> pair_array;

  typedef Graph TargetGraph;
  typedef size_type size_type_trg;
  typedef array_property_map<size_type, vertex_id_map<TargetGraph> >
          vertex_property_map_trg;
  typedef array_property_map<size_type, edge_id_map<TargetGraph> >
          edge_property_map_trg;

  // ************************************************************
  // ** Compute the edge types: 1 if chordless 5-cycle with
  // ** cyclic types is possible; 0 otherwise.
  // ************************************************************
  printf("setting edge types\n");
  fflush(stdout);

  size_type num_edges_processed = 0;

  printf("done setting edge types\n");
  fflush(stdout);

  // MAKE TARGET GRAPH: a 5-clique on (a,b,c,d,e).
  typedef Graph TargetGraph;
  typedef size_type size_type_trg;
  typedef array_property_map<size_type, vertex_id_map<TargetGraph> >
          vertex_property_map_trg;
  typedef array_property_map<size_type, edge_id_map<TargetGraph> >
          edge_property_map_trg;

  TargetGraph target;

  // ************************************************************
  // ** Finding a chordless, typed 5-cycle. *********************
  // ************************************************************
  const size_type_trg numVertsTarget = 5;
  const size_type_trg numEdgesTarget = 5;

  size_type vTypesTarget[numVertsTarget] = { 0, 1, 2, 3, 4 };
//  size_type vTypesTarget[numVertsTarget] = { 0, 0, 0, 0, 0 };
  size_type_trg sourcesTarget[numEdgesTarget] = { 0, 1, 2, 3, 4};
  size_type_trg targetsTarget[numEdgesTarget] = { 1, 2, 3, 4, 0};
  size_type eTypesTarget[numEdgesTarget] = { 1, 1, 1, 1, 1 };

  init(numVertsTarget, numEdgesTarget, sourcesTarget, targetsTarget, target);

  printf("\nTarget graph (%lu, %lu)\n", numVertsTarget, numEdgesTarget);

#ifdef DEBUG
  print(target);
  printf("\n");
#endif

  // ************************************************************

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
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
  walk_verts.push_back(vertices(target)[3]);
  walk_verts.push_back(vertices(target)[4]);
  walk_verts.push_back(vertices(target)[0]);

  dynamic_array<edge_descriptor> walk_edges;
  walk_edges.push_back(edges(target)[0]);
  walk_edges.push_back(edges(target)[1]);
  walk_edges.push_back(edges(target)[2]);
  walk_edges.push_back(edges(target)[3]);
  walk_edges.push_back(edges(target)[4]);

  typedef default_path_printer<Graph> path_printer;
  typedef si_5cycle_visitor<Graph, vertex_property_map_t, Visitor>
          path_evaluator;

  hash_table_t the_edges(2 * num_edges(g));
  hash_edges<Graph> hre(g, the_edges);
  visit_edges(g, hre);

#ifdef DEBUG
  size_type order = num_vertices(g);
  vertex_iterator verts = vertices(g);
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = verts[i];
    printf("vprop[%lu]: %lu\n", get(vid_map, v), vPropmap[v]);
  }
#endif

  bowl_visitor<Graph, vertex_property_map_t, edge_property_map_t>
    bv(g, the_edges, vPropmap, ePropmap, num_edges_processed);

  visit_edges(g, bv);

  typedef si_default_comparator<Graph, vertex_property_map_t,
                                edge_property_map_t, TargetGraph,
                                vertex_property_map_trg,
                                edge_property_map_trg>
          si_comparator;

  si_comparator si_comp(vPropmap, ePropmap, vTypemapTarget, eTypemapTarget);

  path_evaluator si_visit(g, the_edges, vPropmap, user_vis);
  path_printer si_print(g);

/*
  pair_t p(0, 5);
  pair_array pa;
  pa.push_back(p);
  size_type num_candidates =
       subgraph_isomorphism(g, target, walk_verts, walk_edges, si_comp,
                            si_visit, pa);
 */

  // The last argument forces the source vertex of edge 0 in the walk to
  // be the same as the target of edge 5 in the walk (constraining the
  // candidate matches to be closed walks).

  size_type num_candidates = subgraph_isomorphism(g, target, si_comp, si_visit);

  return num_candidates;
}

}

#endif
