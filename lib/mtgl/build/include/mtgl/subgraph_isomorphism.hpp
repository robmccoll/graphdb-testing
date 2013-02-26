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
/*! \file subgraph_isomorphism.hpp

    \brief Searches for a subgraph in the given graph that is isomorphic to
           a target graph.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 12/4/2007

    The big graph and target graphs are independent graph types that can be
    either a primary graph adapter or a wrapper adapter.

    The algorithm takes a function object that defines the comparison made
    for the isomorphic match.  A default comparator is provided
    (si_default_comparator), but a user can define their own.  The only
    restriction is that the () operator needs to take the same paramters as
    in si_default_comparator.  The algorithm can handle forward or reverse
    edges in the walk, and the user doesn't have to do anything special to
    make this work.
*/
/****************************************************************************/

#ifndef MTGL_SUBGRAPH_ISOMORPHISM_HPP
#define MTGL_SUBGRAPH_ISOMORPHISM_HPP

#include <limits>
#include <vector>
#include <cmath>
#include <algorithm>

#ifdef USING_QTHREADS
#include <qthread.h>
#include <qthread/qloop.h>
#include <qthread/qloop.hpp>
#include <qthread/futurelib.h>
#endif

#include <mtgl/partitioning.hpp>
#include <mtgl/psearch.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/subgraph_adapter.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/breadth_first_search.hpp>
#include <mtgl/random_walk.hpp>
#include <mtgl/random.hpp>
#include <mtgl/connected_components.hpp>
#include <mtgl/shiloach_vishkin.hpp>
#include <mtgl/euler_tour.hpp>

#define MAX_LEVELS 100    // Bound on walk length.
#define MAX_RECURSION 10  // Bound on depth of subproblem tree
// #define DEBUG

namespace mtgl {

/*! \brief Default comparison function object that determines if an edge
           from the big graph matches the given walk edge from the target
           graph.

    This function object checks if the type of the edge and the types of the
    source and target vertices from the big graph match those of the currently
    considered walk edge.  It also checkes if the out degrees of the source
    and target vertices in the big graph are at least as big as those from the
    currently considered walk edge from the target graph.  For directed and
    bidirectional graphs, the comparator is called once for each edge.  For
    undirected graphs the comparator is called twice for each edge where each
    endpoint gets a chance to be the "source".

    The comparator function object is passed as a parameter to
    subgraph_isomorphism().  Thus, the user can define their own comparator
    function object to customize the matching.  The only requirement is that
    the declaration of operator() must be identical to the one in the default
    comparator.
*/
template <typename Graph, typename VertexPropertyMap,
          typename EdgePropertyMap, typename GraphTrg = Graph,
          typename VertexPropertyMapTrg = VertexPropertyMap,
          typename EdgePropertyMapTrg = EdgePropertyMap>
class si_default_comparator {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex;
  typedef typename graph_traits<Graph>::edge_descriptor edge;
  typedef typename graph_traits<GraphTrg>::vertex_descriptor vertex_trg;
  typedef typename graph_traits<GraphTrg>::edge_descriptor edge_trg;

  si_default_comparator(VertexPropertyMap& gvt, EdgePropertyMap& get,
                        VertexPropertyMapTrg& trgvt,
                        EdgePropertyMapTrg& trget) :
    g_vtype(gvt), g_etype(get), trg_vtype(trgvt), trg_etype(trget) {}

  bool operator()(const edge& g_e, const edge_trg& trg_e,
                  Graph& g, GraphTrg& trg, size_type walk_level)
  {
#ifdef DOUBLE_DEBUG
    printf("si_comp: level: %lu\n", walk_level);
#endif
    // Get vertices.
    vertex g_src = source(g_e, g);
    vertex g_dest = target(g_e, g);
    vertex_trg trg_src = source(trg_e, trg);
    vertex_trg trg_dest = target(trg_e, trg);

    // The edge is considered a match to the walk edge if its type triple
    // (src vertex type, edge type, dest vertex type) matches the walk edge's
    // corresponding triple, and if the edge's vertices have out degree at
    // least as big as the out degree of the corresponding edge's vertices
    // from the target graph.
    return (g_vtype[g_src] == trg_vtype[trg_src] &&
            g_etype[g_e] == trg_etype[trg_e] &&
            g_vtype[g_dest] == trg_vtype[trg_dest] &&
            out_degree(g_src, g) >= out_degree(trg_src, trg) &&
            out_degree(g_dest, g) >= out_degree(trg_dest, trg));
  }

private:
  VertexPropertyMap& g_vtype;
  EdgePropertyMap& g_etype;
  VertexPropertyMapTrg& trg_vtype;
  EdgePropertyMapTrg& trg_etype;

  /*!
    \fn bool operator()(const edge& g_e, const edge_trg& trg_e,
                        Graph& g, GraphTrg& trg)

    \brief Checks if edge g_e matches edge trg_e where a match is when the
           type triples (source vertex type, edge type, destination vertex
           type) match and the vertices of g_e have out degree at least as
           big as the out degree of the corresponding vertices of trg_e.

    \param g_e The edge from the big graph.
    \param trg_e The edge from the target graph.
    \param g The big graph.
    \param trg The target graph.
  */
};


template <class Graph>
bool verify_graph_and_rev(Graph& g, Graph& rg, char* str, bool print = true)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  if (print) printf("verify_graph_and_rev: %x %x\n", &g, &rg);

  edge_iterator cedges = edges(g);
  edge_iterator credges = edges(rg);

  vertex_id_map<Graph> vid_map1 = get(_vertex_id_map, g);
  vertex_id_map<Graph> vid_map2 = get(_vertex_id_map, rg);
  edge_id_map<Graph> eid_map1 = get(_edge_id_map, g);
  edge_id_map<Graph> eid_map2 = get(_edge_id_map, rg);

  bool rev_correct = true;

  for (size_type i = 0; i < num_edges(g); ++i)
  {
    edge_descriptor c_edge = cedges[i];
    edge_descriptor cr_edge = credges[i];

    if (print)
    {
      printf("%s: {%lu}: (%lu,%lu), rev: {%lu} (%lu, %lu)\n", str,
             get(eid_map1, c_edge),
             get(vid_map1, source(c_edge, g)),
             get(vid_map1, target(c_edge, g)),
             get(eid_map2, cr_edge),
             get(vid_map2, source(cr_edge, rg)),
             get(vid_map2, target(cr_edge, rg)));
    }

    if ((get(vid_map1, source(c_edge, g)) !=
        get(vid_map2, target(cr_edge, rg))) ||
        get(vid_map1, target(c_edge, g)) !=
        get(vid_map2, source(cr_edge, rg)))
    {
      rev_correct = false;
      printf("OOOPS!\n");
    }
  }

  return rev_correct;
}

namespace detail {

/*! \brief This function, which is useful for sorting, returns the value u - v.

    \param u, v: ptrs to uint64_t values

    \author Jon Berry (jberry@sandia.gov)

    \date 8/2010
*/
static inline int uint_cmp(const void* u, const void* v)
{
  return static_cast<int>((*(uint64_t*) u) - (*(uint64_t*) v));
}

template <typename Graph>
class iso_bipartite_visitor : public default_psearch_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex;
  typedef typename graph_traits<Graph>::edge_descriptor edge;
  typedef typename graph_traits<Graph>::size_type size_type;

  iso_bipartite_visitor(Graph& bg, size_type* ec) :
    edge_color(ec), vid_map(get(_vertex_id_map, bg)),
    eid_map(get(_edge_id_map, bg)), bipartite_graph(bg) {}

  // When the other vertex of an edge is visited during a search, the visitor
  // for the other vertex is initialized using the copy constructor.
  iso_bipartite_visitor(const iso_bipartite_visitor& vis) :
    edge_color(vis.edge_color), vid_map(vis.vid_map), eid_map(vis.eid_map),
    bipartite_graph(vis.bipartite_graph) {}

  int visit_test(edge& e, vertex& src)
  {
    size_type eid = get(eid_map, e);
    edge_color[eid] = 1;

    return true;
  }

protected:
  size_type* edge_color;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
  Graph& bipartite_graph;
};

template <typename HT, typename SizeType>
class bipartite_vertex_visitor {
public:
  bipartite_vertex_visitor(typename HT::key_type o, typename HT::key_type* bvl,
                           SizeType* dtos, SizeType* bvid) :
    order(o), bipartite_vertex_level(bvl),
    dense_to_sparse(dtos), b_vid_map(bvid) {}

  void operator()(typename HT::value_type* i)
  {
    bipartite_vertex_level[i->second] = i->first / order;
    b_vid_map[i->second] = dense_to_sparse[i->first % order];

#ifdef DEBUG
    printf("sparse bipartite vertex id: %lu ---> original id: %lu\n", i->first,
           b_vid_map[i->second]);
#endif
  }

private:
  typename HT::key_type order;
  typename HT::key_type* bipartite_vertex_level;
  SizeType* dense_to_sparse;
  SizeType* b_vid_map;
};

template <typename HT, typename SizeType, typename SizeTypeB>
class colored_vertex_visitor {
public:
  colored_vertex_visitor(SizeType* bvid, SizeType* cvid,
                         SizeTypeB* b_v_level, SizeTypeB* c_v_level) :
    b_vid_map(bvid), c_vid_map(cvid), bipartite_vertex_level(b_v_level),
    colored_vertex_level(c_v_level) {}

  void operator()(typename HT::value_type* i)
  {
    c_vid_map[i->second] = b_vid_map[i->first];
    colored_vertex_level[i->second] = bipartite_vertex_level[i->first];

#ifdef DEBUG
    printf("colored bipartite vertex id: %lu ---> original id: %lu\n",
           i->second, c_vid_map[i->second]);
#endif
  }

private:
  SizeType* b_vid_map;
  SizeType* c_vid_map;
  SizeTypeB* bipartite_vertex_level;
  SizeTypeB* colored_vertex_level;
};

template <typename Graph>
class betweenness_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  betweenness_visitor(size_type* vb, size_type* eb) :
    vertex_betweenness(vb), edge_betweenness(eb) {}

  bool visit_test(bool isWhite, edge_descriptor& e, Graph& g)
  {
    vertex_descriptor u = source(e, g);
    vertex_descriptor v = target(e, g);

    size_type uid = get(get(_vertex_id_map, g), u);
    size_type vid = get(get(_vertex_id_map, g), v);
    size_type eid = get(get(_edge_id_map, g), e);

    size_type u_betweenness = vertex_betweenness[uid];
    edge_betweenness[eid] *= u_betweenness;
    size_type val = mt_incr(vertex_betweenness[vid], u_betweenness);

    return val == 0;
  }

protected:
  size_type* vertex_betweenness;
  size_type* edge_betweenness;
};

/*  This was to speed up the special case in which there's only one path
 *  through the colored graph.  It didn't speed things up, but remains for
 *  now as an example of searching the colored graph and referencing the
 *  big graph.
 */
template <typename BipartiteGraph, typename Graph>
class only_one_path_visitor : public default_bfs_visitor<BipartiteGraph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<BipartiteGraph>::size_type size_type_b;
  typedef typename graph_traits<BipartiteGraph>::vertex_descriptor
          vertex_descriptor_b;
  typedef typename graph_traits<BipartiteGraph>::edge_descriptor
          edge_descriptor_b;
  typedef xmt_hash_table<size_type, size_type> dense_map_t;

  only_one_path_visitor(BipartiteGraph& cg, Graph& bigG,
                        dynamic_array<vertex_descriptor>& va,
                        dynamic_array<edge_descriptor>& ea,
                        size_type* c_v,
                        dynamic_array<size_type>& c_e,
                        dense_map_t& bgdvi,
                        dense_map_t& bgdei) :
    colored_graph(cg), path_verts(va), path_edges(ea), c_vid_map(c_v),
    c_eid_map(c_e), big_graph_dense_vertex_id(bgdvi),
    big_graph_dense_edge_id(bgdei), verts(vertices(bigG)),
    edgs(edges(bigG)), vid_map(get(_vertex_id_map, colored_graph)),
    eid_map(get(_edge_id_map, colored_graph)) {}

  void discover_vertex(vertex_descriptor_b& v, BipartiteGraph& g)
  {
    size_type dense_id = (std::numeric_limits<size_type>::max)();
    big_graph_dense_vertex_id.lookup(c_vid_map[get(vid_map, v)], dense_id);
    path_verts.push_back(verts[dense_id]);
  }

  void tree_edge(edge_descriptor_b& e, BipartiteGraph& g)
  {
    size_type eid = c_eid_map[get(eid_map, e)];
    size_type dense_id = (std::numeric_limits<size_type>::max)();
    big_graph_dense_edge_id.lookup(eid, dense_id);
    edge_descriptor ge = edgs[dense_id];
    path_edges.push_back(ge);
  }

protected:
  BipartiteGraph& colored_graph;
  dynamic_array<vertex_descriptor>& path_verts;
  dynamic_array<edge_descriptor>& path_edges;
  vertex_descriptor* c_vid_map;
  dynamic_array<size_type>& c_eid_map;
  dense_map_t& big_graph_dense_vertex_id;
  dense_map_t& big_graph_dense_edge_id;
  vertex_iterator verts;
  edge_iterator edgs;
  vertex_id_map<BipartiteGraph> vid_map;
  edge_id_map<BipartiteGraph> eid_map;
};

template <typename Graph>
class subproblem_visitor : public default_bfs_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  subproblem_visitor(size_type* bc, size_type inc) :
    bread_crumbs(bc), incr(inc) {}

  void tree_edge(edge_descriptor& e, Graph& g)
  {
    vertex_descriptor v = target(e, g);
    size_type vid = get(get(_vertex_id_map, g), v);

    mt_incr(bread_crumbs[vid], incr);

#ifdef DEBUG
    vertex_descriptor u = source(e, g);
    size_type uid = get(get(_vertex_id_map, g), u);
    size_type eid = get(get(_edge_id_map, g), e);

    printf("after tree_edge(%lu, %lu)[%lu]: %lu (incr was %lu)\n",
           uid, vid, eid, bread_crumbs[vid], incr);
#endif
  }

protected:
  size_type* bread_crumbs;
  size_type incr;
};

template <typename GAdapter, typename BAdapter, typename TAdapter,
          typename SICompare>
void assign_sparse_bipartite_ids(
  GAdapter& bigG,
  BAdapter& bipartite_graph,
  TAdapter& targetG,
  bool** active_verts,
  dynamic_array<typename graph_traits<TAdapter>::vertex_descriptor>& walk_verts,
  dynamic_array<typename graph_traits<TAdapter>::edge_descriptor>& walk_edges,
  typename graph_traits<BAdapter>::size_type& b_size,
  dynamic_array<typename graph_traits<BAdapter>::size_type>& b_sources,
  dynamic_array<typename graph_traits<BAdapter>::size_type>& b_dests,
  dynamic_array<typename graph_traits<BAdapter>::size_type>& b_eid_map,
  xmt_hash_table<typename graph_traits<BAdapter>::size_type,
                 typename graph_traits<BAdapter>::size_type>&
    bipartite_vertex_id,
  xmt_hash_table<typename graph_traits<GAdapter>::size_type,
                 typename graph_traits<GAdapter>::size_type>&
    big_graph_dense_vertex_id,
  SICompare& si_compare)
{
  typedef typename graph_traits<GAdapter>::size_type size_type;
  typedef typename graph_traits<BAdapter>::size_type size_type_b;
  typedef typename graph_traits<TAdapter>::size_type size_type_trg;
  typedef typename graph_traits<GAdapter>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<GAdapter>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<TAdapter>::vertex_descriptor
          vertex_descriptor_trg;
  typedef typename graph_traits<TAdapter>::edge_descriptor
          edge_descriptor_trg;
  typedef typename graph_traits<GAdapter>::edge_iterator edge_iterator;
  typedef typename graph_traits<TAdapter>::vertex_iterator vertex_iterator_trg;
  typedef typename graph_traits<TAdapter>::edge_iterator edge_iterator_trg;

  // TODO: Can I modify this second pass somehow to get only the colored edges
  // without doing the extra filtering step?

  // 3. Keep track of the unique vertices and get the bipartite edges for the
  // first walk edge.

  vertex_id_map<GAdapter> vid_map = get(_vertex_id_map, bigG);
  edge_id_map<GAdapter> eid_map = get(_edge_id_map, bigG);
  vertex_id_map<TAdapter> vid_map_trg = get(_vertex_id_map, targetG);

  size_type bigG_order = num_vertices(bigG);
  size_type bigG_size = num_edges(bigG);

  b_size = 0;
  size_type active_level = 0;

  edge_iterator bg_edges = edges(bigG);
  edge_iterator_trg tg_edges = edges(targetG);
  vertex_iterator_trg tg_verts = vertices(targetG);

  vertex_descriptor_trg trg_walk_src = walk_verts[0];
  edge_descriptor_trg trg_e = walk_edges[0];

  // Determine if the target walk edge is in the forward or reversed
  // direction.  If the source vertex in the walk is the source of the
  // target edge, then the target edge is a forward edge.  Otherwise, the
  // target edge is a backward edge.
  size_type_trg trg_walk_src_id = get(vid_map_trg, trg_walk_src);
  size_type_trg trg_src_id = get(vid_map_trg, source(trg_e, targetG));
  bool trg_e_reversed = trg_walk_src_id != trg_src_id;

  // Add the unique vertices to the bipartite graph.
  #pragma mta assert parallel
  #pragma mta assert nodep
  for (size_type i = 0; i < bigG_size; ++i)
  {
    edge_descriptor e = bg_edges[i];

    // Get the "source" and "target" vertices for the big graph edge.  For
    // forward target edges, the "source" of the big graph edge is the edge's
    // actual source.  For backward target edges, the "source" of the big
    // graph edge is the edge's target.
    vertex_descriptor bg_src = trg_e_reversed ? target(e, bigG) :
                               source(e, bigG);
    vertex_descriptor bg_dest = bg_src == source(e, bigG) ? target(e, bigG) :
                                source(e, bigG);

    // Get the vertex and edge ids.
    size_type bg_src_id = get(vid_map, bg_src);
    size_type bg_dest_id = get(vid_map, bg_dest);
    size_type e_id = get(eid_map, e);

    // Get the dense vertex ids.
    size_type bg_dense_src_id = (size_type) 0, bg_dense_dest_id = (size_type) 0;
    big_graph_dense_vertex_id.lookup(bg_src_id, bg_dense_src_id);
    big_graph_dense_vertex_id.lookup(bg_dest_id, bg_dense_dest_id);

    // Check if the first walk edge matches this edge from the big graph.
    if (si_compare(e, trg_e, bigG, targetG, active_level))
    {
      size_type_b b_eid = mt_incr(b_size, 1);

      b_sources[b_eid] = bg_dense_src_id;
      b_dests[b_eid] = bigG_order + bg_dense_dest_id;
      b_eid_map[b_eid] = e_id;

      bipartite_vertex_id.insert(bg_dense_src_id, 0);
      bipartite_vertex_id.insert(bigG_order + bg_dense_dest_id, 0);

#ifdef DEBUG
      printf("(%lu, %lu): active level %d: added (%lu, %lu) "
             "(dense: (%lu, %lu))\n",
             bg_src_id, bg_dest_id, 0,
             b_sources[b_eid], b_dests[b_eid],
             bg_dense_src_id, bg_dense_dest_id);
      fflush(stdout);
#endif
    }

    // If the graph is undirected, also check the reverse of this edge.
    if (is_undirected(bigG))
    {
      edge_descriptor e_rev = edge_descriptor(target(e, bigG),
                                              source(e, bigG), get(eid_map, e));

      if (si_compare(e_rev, trg_e, bigG, targetG, active_level))
      {
        size_type_b b_eid = mt_incr(b_size, 1);

        b_sources[b_eid] = bg_dense_dest_id;
        b_dests[b_eid] = bigG_order + bg_dense_src_id;
        b_eid_map[b_eid] = e_id;

        bipartite_vertex_id.insert(bg_dense_dest_id, 0);
        bipartite_vertex_id.insert(bigG_order + bg_dense_src_id, 0);

#ifdef DEBUG
        printf("(%lu, %lu): active level %d: added (%lu, %lu) "
               "(dense: (%lu, %lu))\n", bg_dest_id, bg_src_id, 0,
               b_sources[b_eid], b_dests[b_eid],
               bg_dense_dest_id, bg_dense_src_id);
        fflush(stdout);
#endif
      }
    }
  }

  // 4. Keep track of the unique vertices and get the bipartite edges for the
  // remaining walk edges.

  active_level = 2;  // Level in bipartite graph.

  // Try to match the remainder of the walk through the target graph against
  // the filtered graph.
  for (size_type_trg j = 1; j < walk_edges.size(); ++j)
  {
    vertex_descriptor_trg trg_walk_src = walk_verts[j];
    edge_descriptor_trg trg_e = walk_edges[j];

    // Determine if the target walk edge is in the forward or reversed
    // direction.  If the source vertex in the walk is the source of the
    // target edge, then the target edge is a forward edge.  Otherwise, the
    // target edge is a backward edge.
    size_type_trg trg_walk_src_id = get(vid_map_trg, trg_walk_src);
    size_type_trg trg_src_id = get(vid_map_trg, source(trg_e, targetG));
    bool trg_e_reversed = trg_walk_src_id != trg_src_id;

    // Add the unique vertices to the bipartite graph.
    #pragma mta assert parallel
    #pragma mta assert nodep
    for (size_type i = 0; i < bigG_size; ++i)
    {
      edge_descriptor e = bg_edges[i];

      // Get the "source" and "target" vertices for the big graph edge.  For
      // forward target edges, the "source" of the big graph edge is the edge's
      // actual source.  For backward target edges, the "source" of the big
      // graph edge is the edge's target.
      vertex_descriptor bg_src = trg_e_reversed ? target(e, bigG) :
                                 source(e, bigG);
      vertex_descriptor bg_dest = bg_src == source(e, bigG) ? target(e, bigG) :
                                  source(e, bigG);

      // Get the vertex and edge ids.
      size_type bg_src_id = get(vid_map, bg_src);
      size_type bg_dest_id = get(vid_map, bg_dest);
      size_type e_id = get(eid_map, e);

      // Get the dense vertex ids.
      size_type bg_dense_src_id = 0;
      size_type bg_dense_dest_id = 0;
      big_graph_dense_vertex_id.lookup(bg_src_id, bg_dense_src_id);
      big_graph_dense_vertex_id.lookup(bg_dest_id, bg_dense_dest_id);

      // Check if the current walk edge matches this edge from the big graph.
      if (active_verts[active_level - 1][bg_dense_src_id] &&
          si_compare(e, trg_e, bigG, targetG, active_level - 1))
      {
        size_type_b b_eid = mt_incr(b_size, 1);

        b_sources[b_eid] = (active_level - 1) * bigG_order + bg_dense_src_id;
        b_dests[b_eid] = active_level * bigG_order + bg_dense_dest_id;
        b_eid_map[b_eid] = e_id;

        bipartite_vertex_id.insert(
          active_level * bigG_order + bg_dense_dest_id, 0);

#ifdef DEBUG
        printf("(%lu, %lu): active level %lu: added (%lu, %lu) "
               "(dense: (%lu, %lu))\n",
               bg_src_id, bg_dest_id, active_level,
               b_sources[b_eid], b_dests[b_eid],
               bg_dense_src_id, bg_dense_dest_id);
        fflush(stdout);
#endif
      }

      // If the graph is undirected, also check the reverse of this edge.
      if (is_undirected(bigG))
      {
        edge_descriptor e_rev =
          edge_descriptor(target(e, bigG), source(e, bigG), get(eid_map, e));

        if (active_verts[active_level - 1][bg_dense_dest_id] &&
            si_compare(e_rev, trg_e, bigG, targetG, active_level - 1))
        {
          size_type_b b_eid = mt_incr(b_size, 1);

          b_sources[b_eid] = (active_level - 1) * bigG_order + bg_dense_dest_id;
          b_dests[b_eid] = active_level * bigG_order + bg_dense_src_id;
          b_eid_map[b_eid] = e_id;

          bipartite_vertex_id.insert(
            active_level * bigG_order + bg_dense_src_id, 0);

#ifdef DEBUG
          printf("(%lu, %lu): active level %lu: added (%lu, %lu) "
                 "(dense: (%lu, %lu))\n",
                 bg_dest_id, bg_src_id, active_level,
                 b_sources[b_eid], b_dests[b_eid],
                 bg_dense_src_id, bg_dense_dest_id);
          fflush(stdout);
#endif
        }
      }
    }

    ++active_level;
  }
}

template <typename Graph>
bool
select_subproblem(Graph& colored_graph, Graph& colored_graph_rev,
                  typename graph_traits<Graph>::size_type c_order,
                  typename graph_traits<Graph>::size_type* bread_crumbs,
                  typename graph_traits<Graph>::vertex_descriptor u,
                  typename graph_traits<Graph>::vertex_descriptor v,
                  typename graph_traits<Graph>::size_type* colored_vertex_level)
{
  typedef typename graph_traits<Graph>::size_type size_type;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, colored_graph);
  size_type uid = get(vid_map, u);
  size_type vid = get(vid_map, v);

  // We will search downward from u and upward from v.  All edges that
  // are hit twice get bread crumbs. Then we search upward from u and
  // downward from v, distributing bread crumbs.

  assert(colored_vertex_level[u] < colored_vertex_level[v]);

  // Initialize the bread crumbs.
  for (size_type i = 0; i < c_order; ++i) bread_crumbs[i] = 0;
  bread_crumbs[uid] = 1;

  // Search down from u, laying one bread crumb apiece.
  detail::subproblem_visitor<Graph> bc1(bread_crumbs, 1);

  breadth_first_search(colored_graph, u, bc1);

#ifdef DEBUG
  printf("After forward search from %lu\n", uid);
  for (size_type i = 0; i < num_vertices(colored_graph); ++i)
  {
    printf("bread_crumbs[%lu]: %d\n", i, bread_crumbs[i]);
  }
#endif

  ++bread_crumbs[vid];

  // Search up from v, laying one bread crumb apiece.
  breadth_first_search(colored_graph_rev, v, bc1);

#ifdef DEBUG
  printf("After backward search from %lu\n", vid);
  for (size_type i = 0; i < num_vertices(colored_graph); ++i)
  {
    printf("bread_crumbs[%lu]: %d\n", i, bread_crumbs[i]);
  }
#endif

  if (!((bread_crumbs[uid] == 2) && (bread_crumbs[vid] == 2)))
  {
    // Empty subproblem.
    return false;
  }

  // Wipe clean extra bread crumbs.
  for (size_type i = 0; i < c_order; ++i)
  {
    if (bread_crumbs[i] < 2) bread_crumbs[i] = 0;
  }

  // Search up from u, laying two breads crumb apiece.
  detail::subproblem_visitor<Graph> bc2(bread_crumbs, 2);

  breadth_first_search(colored_graph_rev, u, bc2);

#ifdef DEBUG
  printf("After backward search from %lu\n", uid);
  for (size_type i = 0; i < num_vertices(colored_graph); ++i)
  {
    printf("bread_crumbs[%lu]: %d\n", i, bread_crumbs[i]);
  }
#endif

  // Search down from v, laying one bread crumb apiece.
  breadth_first_search(colored_graph, v, bc2);

  // The subproblem is induced by all vertices that have 2 bread crumbs.

#ifdef DEBUG
  printf("After forward search from %lu\n", vid);
  for (size_type i = 0; i < num_vertices(colored_graph); ++i)
  {
    printf("bread_crumbs[%lu]: %d\n", i, bread_crumbs[i]);
  }
#endif

  return true;
}

template <typename Graph>
void
special_betweenness(Graph& colored_graph, Graph& colored_graph_rev,
                    typename graph_traits<Graph>::size_type c_order,
                    typename graph_traits<Graph>::size_type c_size,
                    typename graph_traits<Graph>::size_type* edge_betweenness,
                    typename graph_traits<Graph>::size_type cv_src,
                    typename graph_traits<Graph>::size_type cv_trg)
{
  typedef typename graph_traits<Graph>::size_type size_type_b;

  // For each of the betweenness directions, we assign a betweenness of 1 to
  // the initial vertex.  We use the vertex betweenness as a temporary value
  // in the calculation of the edge betweenness.  It doesn't really represent
  // an actual vertex betweenness.  We assign the edge betweenness to be its
  // source vertex's temporary betweenness value.  Except for the initial
  // vertex, the temporary vertex betweenness value for a vertex is the sum
  // of the edge betweenness of its in edges.

  // Declare the vertex and edge betweenness arrays.
  size_type_b* vertex_betweenness =
    (size_type_b*) malloc(c_order * sizeof(size_type_b));

  // Initialize the betweenness for all vertices to 0.
  for (size_type_b i = 0; i < c_order; ++i) vertex_betweenness[i] = 0;

  // Set the vertex betweenness for the supersource to 1.
  vertex_betweenness[cv_src] = 1;

  // Initialize the betweenness for all edges to 1.
  #pragma mta assert nodep
  for (size_type_b i = 0; i < c_size; ++i) edge_betweenness[i] = 1;

  // Declare the betweennes visitor, and find the forward betweenness.
  detail::betweenness_visitor<Graph>
  b_vis(vertex_betweenness, edge_betweenness);

  breadth_first_search(colored_graph, vertices(colored_graph)[cv_src], b_vis);

  // Initialize the betweenness for all vertices to 0.
  for (size_type_b i = 0; i < c_order; ++i) vertex_betweenness[i] = 0;

  // Set the vertex betweenness for the supersink to 1.
  vertex_betweenness[cv_trg] = 1;

  // Find the reverse betweenness.
  breadth_first_search(colored_graph_rev,
                       vertices(colored_graph_rev)[cv_trg], b_vis);

  free(vertex_betweenness);
}

template <typename Graph, typename BipartiteAdapter, typename GlobalSizeType,
          typename Visitor, int DEPTH = 0>
class process_subgraph_of_colored {
public:
  typedef typename graph_traits<BipartiteAdapter>::size_type size_type;
#ifdef USING_QTHREADS
  typedef size_t qt_size_t;
#else
  typedef size_type qt_size_t;
#endif
  typedef typename graph_traits<BipartiteAdapter>::vertex_iterator
          vertex_iterator;
  typedef typename graph_traits<BipartiteAdapter>::vertex_descriptor
          vertex_descriptor;
  typedef typename graph_traits<BipartiteAdapter>::edge_descriptor
          edge_descriptor;
  typedef typename graph_traits<BipartiteAdapter>::edge_iterator edge_iterator;
  typedef typename graph_traits<BipartiteAdapter>::out_edge_iterator
          out_edge_iterator;
  typedef subgraph_adapter<BipartiteAdapter> subgraph_t;
  typedef pair<size_type, size_type> pair_s;
  typedef dynamic_array<pair<size_type, size_type> > pair_array;
  typedef pair<vertex_descriptor, vertex_descriptor> pair_v;
  typedef xmt_hash_table<typename graph_traits<BipartiteAdapter>::size_type,
                         typename graph_traits<BipartiteAdapter>::size_type>
          dense_map_t;

  process_subgraph_of_colored(Graph& g, BipartiteAdapter& cg,
                              BipartiteAdapter& cgr,
                              dense_map_t& dvm,
                              dense_map_t& dem,
                              Visitor& vs) :
    bigG(g), colored_graph(cg), colored_graph_rev(cgr),
    big_graph_dense_vertex_id(dvm), big_graph_dense_edge_id(dem), vis(vs) {}

  void operator()(size_type num_sub_verts, vertex_descriptor* sub_verts,
                  size_type levels_in_colored,
                  size_type c_order, size_type c_size,
                  size_type cv_src, size_type cv_trg,
                  GlobalSizeType* c_vid_map,
                  dynamic_array<GlobalSizeType>& c_eid_map,
                  size_type& num_candidates,
                  size_type* colored_vertex_level,
                  pair_array pa = pair_array())
  {
    subgraph_t sub_colored(colored_graph);
    subgraph_t sub_colored_rev(colored_graph_rev);

    init_vertices(sub_verts, num_sub_verts, sub_colored);
    init_vertices(sub_verts, num_sub_verts, sub_colored_rev);
//    verify_graph_and_rev(colored_graph, colored_graph_rev, "next cc (colored)");
//    verify_graph_and_rev(sub_colored, sub_colored_rev, "next cc (sub)");

    process(sub_colored, sub_colored_rev, levels_in_colored, c_order, c_size,
            cv_src, cv_trg, c_vid_map, c_eid_map, num_candidates,
            colored_vertex_level, pa);
  }

  void operator()(size_type num_sub_edges, edge_descriptor* sub_edges,
                  edge_descriptor* sub_edges_rev,
                  size_type levels_in_colored,
                  size_type c_order, size_type c_size,
                  size_type cv_src, size_type cv_trg,
                  GlobalSizeType* c_vid_map,
                  dynamic_array<GlobalSizeType>& c_eid_map,
                  size_type& num_candidates,
                  size_type* colored_vertex_level,
                  pair_array pa = pair_array())
  {
    subgraph_t sub_colored(colored_graph);
    subgraph_t sub_colored_rev(colored_graph_rev);

    init_edges(sub_edges, num_sub_edges, sub_colored);
    init_edges(sub_edges_rev, num_sub_edges, sub_colored_rev);
//    verify_graph_and_rev(colored_graph, colored_graph_rev, "next cc (colored)");
//    verify_graph_and_rev(sub_colored, sub_colored_rev, "next cc (sub)");

    process(sub_colored, sub_colored_rev, levels_in_colored, c_order, c_size,
            cv_src, cv_trg, c_vid_map, c_eid_map, num_candidates,
            colored_vertex_level, pa);
  }

  void operator()(edge_property_map<BipartiteAdapter, bool>& emask,
                  size_type levels_in_colored,
                  size_type c_order, size_type c_size,
                  size_type cv_src, size_type cv_trg,
                  GlobalSizeType* c_vid_map,
                  dynamic_array<GlobalSizeType>& c_eid_map,
                  size_type& num_candidates,
                  size_type* colored_vertex_level,
                  pair_array pa = pair_array())
  {
    subgraph_t sub_colored(colored_graph);
    subgraph_t sub_colored_rev(colored_graph_rev);

    init_edges(emask, sub_colored);
    init_edges(emask, sub_colored_rev);
//    verify_graph_and_rev(colored_graph, colored_graph_rev, "colored_graph");
//    verify_graph_and_rev(sub_colored, sub_colored_rev, "sub_colored");



#ifdef DEBUG
    for (size_type i = 0; i < num_edges(colored_graph); ++i)
    {
      printf("emask[%lu]: %lu\n", i, emask[edges(colored_graph)[i]]);
    }
#endif

    process(sub_colored, sub_colored_rev, levels_in_colored, c_order, c_size,
            cv_src, cv_trg, c_vid_map, c_eid_map, num_candidates,
            colored_vertex_level, pa);
  }

  void process(subgraph_t& sub_colored, subgraph_t& sub_colored_rev,
               size_type levels_in_colored, size_type c_order, size_type c_size,
               size_type cv_src, size_type cv_trg, GlobalSizeType* c_vid_map,
               dynamic_array<GlobalSizeType>& c_eid_map,
               size_type& num_candidates, size_type* colored_vertex_level,
               pair_array pa = pair_array())
  {
#ifdef DEBUG
    printf("colored graph: \n");
    print(colored_graph);
    printf("colored graph_rev: \n");
    print(colored_graph_rev);
    printf("subgraph of colored: \n");
    print(sub_colored);
    printf("subgraph of colored rev: \n");
    print(sub_colored_rev);
    printf("process_sub: %lu, %lu\n", num_vertices(sub_colored),
           num_edges(sub_colored));
    fflush(stdout);
#endif

    size_type sc_order = num_vertices(sub_colored);
    size_type sc_size = num_edges(sub_colored);

    GlobalSizeType* sc_vid_map =
      (GlobalSizeType*) malloc(sizeof(GlobalSizeType) * c_order);

    dynamic_array<GlobalSizeType> sc_eid_map(c_size);

    vertex_iterator cvrts = vertices(colored_graph);
    edge_iterator cedgs = edges(colored_graph);

    vertex_descriptor cv_src_desc = cvrts[cv_src];
    vertex_descriptor cv_trg_desc = cvrts[cv_trg];

    vertex_descriptor scv_src = global_to_local(cv_src_desc, sub_colored);
    vertex_descriptor scv_trg = global_to_local(cv_trg_desc, sub_colored);

    vertex_iterator scvrts = vertices(sub_colored);
    edge_iterator scedgs = edges(sub_colored);

    vertex_id_map<BipartiteAdapter> cvid = get(_vertex_id_map, colored_graph);
    edge_id_map<BipartiteAdapter> ceid = get(_edge_id_map, colored_graph);

    vertex_id_map<subgraph_t> scvid = get(_vertex_id_map, sub_colored);
    vertex_id_map<subgraph_t> csvid = get(_vertex_id_map, sub_colored);
    edge_id_map<subgraph_t> cseid = get(_edge_id_map, sub_colored);

    for (size_type i = 0; i < sc_order; ++i)
    {
      vertex_descriptor v = local_to_global(scvrts[i], sub_colored);
      sc_vid_map[get(csvid, scvrts[i])] = c_vid_map[get(cvid, v)];

#ifdef DEBUG
      printf("sc_vid: %lu, orig id: %lu\n",
             get(csvid, scvrts[i]), c_vid_map[get(cvid, v)]);
#endif
    }

    for (size_type i = 0; i < sc_size; ++i)
    {
      edge_descriptor e = local_to_global(scedgs[i], sub_colored);
      sc_eid_map[get(cseid, scedgs[i])] = c_eid_map[get(ceid, e)];
    }

    if (pa.size() == 0)
    {
      process_colored_graph(bigG, sub_colored, sub_colored_rev,
                            big_graph_dense_vertex_id, big_graph_dense_edge_id,
                            levels_in_colored,
                            sc_order, sc_size, scv_src, scv_trg,
                            sc_vid_map, sc_eid_map, num_candidates, vis);
    }
    else
    {
      pair_s p = pa[pa.size() - 1];
      assert(p.first != p.second);

      pair_array pa_cdr(pa);
      pa_cdr.pop_back();

      // We're going to enforce a constraint that the original graph
      // vertex at level p.first must be the same as the original graph
      // vertex at level p.second.
      //
      // Note that the number of levels in any subproblem will never
      // change.  Each subproblem can be thought of as a vertical slice
      // of the bipartite graph that includes the super-source at the top
      // and the super-sink at the bottom.
      size_type* sub_colored_level =
        (size_type*) malloc(sizeof(size_type) * sc_order);

      for (size_type i = 0; i < sc_order; ++i)
      {
        vertex_descriptor v = local_to_global(scvrts[i], sub_colored);

        sub_colored_level[get(csvid, scvrts[i])] =
          colored_vertex_level[get(cvid, v)];
      }

      std::vector<pair_v> pairs;
      find_like_pairs(p.first, p.second, sub_colored, sc_vid_map,
                      sub_colored_level, pairs);

#ifdef DEBUG
      printf("COLORED GRAPH: %lu vertices; %lu pairs\n",
             num_vertices(sub_colored), pairs.size());
#endif

      // For each pair, we need to search and mark the current portion of the
      // colored graph ("sub_colored"), then generate a recursive call to
      // this function with no pair argument.
      subproblem_loop_body splb(bigG, sub_colored, sub_colored_rev,
                                big_graph_dense_vertex_id,
                                big_graph_dense_edge_id,
                                sub_colored_level, levels_in_colored,
                                sc_order, sc_size, scv_src, scv_trg,
                                sc_vid_map, sc_eid_map,
                                num_candidates, vis, pairs, pa_cdr);

      size_type sz = pairs.size();

#ifdef USING_QTHREADS
      qt_loop_balance<subproblem_loop_body>(0, sz, splb);
#else
      for (size_type i = 0; i < sz; ++i) splb((qt_size_t) i, (qt_size_t) i + 1);
#endif

      free(sub_colored_level);
    }

    free(sc_vid_map);
  }

  class subproblem_loop_body {
  public:
    subproblem_loop_body(Graph& g, subgraph_t& sub_colord,
                         subgraph_t& sub_colord_rev,
                         dense_map_t& dvm,
                         dense_map_t& dem,
                         size_type* sub_colord_level,
                         size_type levels_in_colord,
                         size_type sc_ordr, size_type sc_siz,
                         size_type scv_sr, size_type scv_tr,
                         GlobalSizeType* sc_vid_mp,
                         dynamic_array<GlobalSizeType>& sc_eid_mp,
                         size_type& num_candidats, Visitor& vs,
                         std::vector<pair_v>& prs,
                         pair_array& remning_level_prs) :
      bigG(g), sub_colored(sub_colord), sub_colored_rev(sub_colord_rev),
      big_graph_dense_vertex_id(dvm), big_graph_dense_edge_id(dem),
      sub_colored_level(sub_colord_level), levels_in_colored(levels_in_colord),
      sc_order(sc_ordr), sc_size(sc_siz), scv_src(scv_sr), scv_trg(scv_tr),
      sc_vid_map(sc_vid_mp), sc_eid_map(sc_eid_mp),
      num_candidates(num_candidats), vis(vs), pairs(prs),
      remaining_level_pairs(remning_level_prs) {}

  private:
    Graph& bigG;
    subgraph_t& sub_colored;
    subgraph_t& sub_colored_rev;
    dense_map_t& big_graph_dense_vertex_id;
    dense_map_t& big_graph_dense_edge_id;
    size_type* sub_colored_level;
    size_type levels_in_colored;
    size_type sc_order;
    size_type sc_size;
    size_type scv_src;
    size_type scv_trg;
    GlobalSizeType* sc_vid_map;
    dynamic_array<GlobalSizeType>& sc_eid_map;
    size_type& num_candidates;
    Visitor& vis;
    std::vector<pair_v>& pairs;
    pair_array& remaining_level_pairs;

  public:
    void operator()(const qt_size_t start_at, const qt_size_t stop_at)
    {
      vertex_id_map<subgraph_t> csvid = get(_vertex_id_map, sub_colored);

      for (qt_size_t i = start_at; i < stop_at; ++i)
      {
        // Save loop for startat, stopat.
        pair_v p = pairs[i];
        vertex_descriptor u = p.first;
        vertex_descriptor v = p.second;

#ifdef DEBUG
        printf("COLORED pair ( %lu [%lu], %lu [%lu] )\n",
               get(csvid, u), sc_vid_map[get(csvid, u)],
               get(csvid, v), sc_vid_map[get(csvid, v)]);
#endif

        size_type* bread_crumbs =
          (size_type*) malloc(sizeof(size_type) * sc_order);

        if (!select_subproblem(sub_colored, sub_colored_rev, sc_order,
                               bread_crumbs, u, v, sub_colored_level))
        {
          free(bread_crumbs);
          return;
        }

        edge_property_map<subgraph_t, bool> emask(sub_colored);
        bool nonempty = false;

        edge_iterator edgs = edges(sub_colored);

        for (size_type i = 0; i < sc_size; ++i)
        {
          edge_descriptor e = edgs[i];
          vertex_descriptor u = source(e, sub_colored);
          vertex_descriptor v = target(e, sub_colored);

          emask[e] = ((bread_crumbs[get(csvid, u)] == 2) &&
                      (bread_crumbs[get(csvid, v)] == 2));

          if (emask[e]) nonempty = true;
        }

        if (nonempty)
        {
          process_subgraph_of_colored<Graph, subgraph_t, GlobalSizeType,
                                      Visitor, DEPTH + 1>
          psc(bigG, sub_colored, sub_colored_rev,
              big_graph_dense_vertex_id, big_graph_dense_edge_id, vis);

          psc(emask, levels_in_colored, sc_order, sc_size, scv_src, scv_trg,
              sc_vid_map, sc_eid_map, num_candidates, sub_colored_level,
              remaining_level_pairs);
        }

        free(bread_crumbs);
      }
    }
  };

private:
  Graph& bigG;
  BipartiteAdapter& colored_graph;
  BipartiteAdapter& colored_graph_rev;
  dense_map_t& big_graph_dense_vertex_id;
  dense_map_t& big_graph_dense_edge_id;
  Visitor& vis;
};

template <typename Graph, typename BipartiteAdapter, typename GlobalSizeType,
          typename Visitor>
class process_subgraph_of_colored<Graph, BipartiteAdapter, GlobalSizeType,
                                  Visitor, MAX_RECURSION> {
public:
  typedef typename graph_traits<BipartiteAdapter>::size_type size_type;
  typedef typename graph_traits<BipartiteAdapter>::vertex_iterator
          vertex_iterator;
  typedef typename graph_traits<BipartiteAdapter>::vertex_descriptor
          vertex_descriptor;
  typedef typename graph_traits<BipartiteAdapter>::edge_descriptor
          edge_descriptor;
  typedef typename graph_traits<BipartiteAdapter>::edge_iterator edge_iterator;
  typedef typename graph_traits<BipartiteAdapter>::out_edge_iterator
          out_edge_iterator;
  typedef subgraph_adapter<BipartiteAdapter> subgraph_t;
  typedef pair<size_type, size_type> pair_s;
  typedef dynamic_array<pair_s> pair_array;
  typedef xmt_hash_table<typename graph_traits<BipartiteAdapter>::size_type,
                         typename graph_traits<BipartiteAdapter>::size_type>
          dense_map_t;

  process_subgraph_of_colored(Graph& g, BipartiteAdapter& cg,
                              BipartiteAdapter& cgr,
                              dense_map_t& dvm,
                              dense_map_t& dem,
                              Visitor& vs) :
    bigG(g), colored_graph(cg), colored_graph_rev(cgr),
    big_graph_dense_vertex_id(dvm), big_graph_dense_edge_id(dem), vis(vs) {}

  void operator()(edge_property_map<BipartiteAdapter, bool>& emask,
                  size_type levels_in_colored,
                  size_type c_order, size_type c_size,
                  size_type cv_src, size_type cv_trg,
                  GlobalSizeType* c_vid_map,
                  dynamic_array<GlobalSizeType>& c_eid_map,
                  size_type& num_candidates, size_type* colored_vertex_level,
                  pair_array& pa = pair_array())
  {
    fprintf(stderr, "error: too many constraints led to template ");
    fprintf(stderr, "overflow (max = %d constraints. ", MAX_RECURSION);
    exit(1);
  }

  void operator()(size_type num_sub_verts, vertex_descriptor* sub_verts,
                  size_type levels_in_colored,
                  size_type c_order, size_type c_size,
                  size_type cv_src, size_type cv_trg,
                  size_type* c_vid_map, dynamic_array<size_type>& c_eid_map,
                  size_type& num_candidates, size_type* colored_vertex_level,
                  pair_array& pa = pair_array())
  {
    fprintf(stderr, "error: too many constraints led to template ");
    fprintf(stderr, "overflow (max = %d constraints. ", MAX_RECURSION);
    exit(1);
  }

private:
  Graph& bigG;
  BipartiteAdapter& colored_graph;
  BipartiteAdapter& colored_graph_rev;
  dense_map_t& big_graph_dense_vertex_id;
  dense_map_t& big_graph_dense_edge_id;
  Visitor& vis;
};

template <typename Graph, typename BipartiteAdapter, typename GlobalSizeType,
          typename Visitor>
void
process_colored_graph(
  Graph& bigG, BipartiteAdapter& colored_graph,
  BipartiteAdapter& colored_graph_rev,
  xmt_hash_table<typename graph_traits<BipartiteAdapter>::size_type,
                 typename graph_traits<BipartiteAdapter>::size_type>&
    big_graph_dense_vertex_id,
  xmt_hash_table<typename graph_traits<BipartiteAdapter>::size_type,
                 typename graph_traits<BipartiteAdapter>::size_type>&
    big_graph_dense_edge_id,
  typename graph_traits<BipartiteAdapter>::size_type levels_in_colored,
  typename graph_traits<BipartiteAdapter>::size_type c_order,
  typename graph_traits<BipartiteAdapter>::size_type c_size,
  typename graph_traits<BipartiteAdapter>::size_type cv_src,
  typename graph_traits<BipartiteAdapter>::size_type cv_trg,
  GlobalSizeType* c_vid_map,
  dynamic_array<GlobalSizeType>& c_eid_map,
  typename graph_traits<BipartiteAdapter>::size_type& num_candidates,
  Visitor& vis)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<BipartiteAdapter>::size_type size_type_b;
  typedef typename graph_traits<BipartiteAdapter>::vertex_iterator
          vertex_iterator_b;
  typedef typename graph_traits<BipartiteAdapter>::vertex_descriptor
          vertex_descriptor_b;
  typedef typename graph_traits<BipartiteAdapter>::edge_descriptor
          edge_descriptor_b;
  typedef typename graph_traits<BipartiteAdapter>::edge_iterator
          edge_iterator_b;
  typedef typename graph_traits<BipartiteAdapter>::out_edge_iterator
          out_edge_iterator_b;
  typedef array_property_map<size_type_b, edge_id_map<BipartiteAdapter> >
          edge_property_map_b;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, bigG);
  edge_id_map<Graph> eid_map = get(_edge_id_map, bigG);
  vertex_id_map<BipartiteAdapter> vid_map_c = get(_vertex_id_map,
                                                   colored_graph);
  edge_id_map<BipartiteAdapter> eid_map_c = get(_edge_id_map, colored_graph);

  // For each of the betweenness directions, we assign a betweenness of 1 to
  // the initial vertex.  We use the vertex betweenness as a temporary value
  // in the calculation of the edge betweenness.  It doesn't really represent
  // an actual vertex betweenness.  We assign the edge betweenness to be its
  // source vertex's temporary betweenness value.  Except for the initial
  // vertex, the temporary vertex betweenness value for a vertex is the sum
  // of the edge betweenness of its in edges.

  size_type_b* edge_betweenness =
    (size_type_b*) malloc(c_size * sizeof(size_type_b));

  special_betweenness(colored_graph, colored_graph_rev, c_order, c_size,
                      edge_betweenness, cv_src, cv_trg);

#ifdef DEBUG
  for (size_type_b i = 0; i < c_size; ++i)
  {
    printf("edge_betweenness of edge %lu: %lu\n", i, edge_betweenness[i]);
  }
#endif

  size_type_b* cumulative_edge_betweenness =
    (size_type_b*) malloc(c_size * sizeof(size_type_b));

  for (size_type_b i = 0; i < c_size; ++i) cumulative_edge_betweenness[i] = 1;

  num_candidates = 0;

// size_type_b cv_src;
// colored_vertex_id.lookup(b_order + 1, cv_src);
  out_edge_iterator_b src_out_edges = out_edges(cv_src, colored_graph);
  size_type_b deg_cv_src = out_degree(cv_src, colored_graph);

  // Set the number of candidates to be the sum of the edge_betweenness
  // values of the edges coming out of the source.
  #pragma mta assert nodep
  for (size_type_b i = 0; i < deg_cv_src; ++i)
  {
    edge_descriptor_b e = src_out_edges[i];

    // TODO: Had better make sure this reduction gets removed.
    num_candidates += edge_betweenness[get(eid_map_c, e)];

#ifdef DEBUG
    printf("incrementing num_candidates by %lu\n",
           edge_betweenness[get(eid_map_c, e)]);
#endif
  }

  edge_property_map_b eweight(edge_betweenness, eid_map_c);
  edge_property_map_b cumulative_eweight(cumulative_edge_betweenness,
                                         eid_map_c);

  vertex_iterator_b cverts = vertices(colored_graph);
  edge_iterator_b cedgs = edges(colored_graph);
  vertex_descriptor_b walk_src = cverts[cv_src];

  levels_in_colored = levels_in_colored * 2 + 1;

  // If the number of random walks is set to k * (T * log(T)), then the
  // probability of missing a candidate match is 1 / n^k.  Find k such that
  // this probability is at most one in a million.

  // ** First, get T log T (being careful in case T == 1).
  if (num_candidates > 1)
  {
    double logT = log(static_cast<double>(num_candidates)) / log(2.0);

    // Now, compute the minimum k we need.
    double sevenLogTen = 7 * (log(10.0) / log(2.0));
    int minK = (int) ceil(sevenLogTen / logT);
    num_candidates = static_cast<size_type_b>(minK * num_candidates * logT);
  }
  else
  {
/*
    dynamic_array<vertex_descriptor> path_verts;
    dynamic_array<edge_descriptor> path_edges;

    vertex_iterator verts = vertices(bigG);
    edge_iterator edgs = edges(bigG);

    detail::only_one_path_visitor<BipartiteAdapter, Graph>
      oopv(colored_graph, bigG, path_verts, path_edges, c_vid_map, c_eid_map,
           big_graph_dense_vertex_id, big_graph_dense_edge_id);

    breadth_first_search(colored_graph, walk_src, oopv);

    // Get rid of super-source, super-sink.
    dynamic_array<vertex_descriptor> path_verts2;
    dynamic_array<edge_descriptor> path_edges2;

    size_type sz = path_verts.size();
    for (size_type i = 1; i < sz - 1; ++i) path_verts2.push_back(path_verts[i]);

    sz = path_edges.size();
    for (size_type i = 1; i < sz - 1; ++i) path_edges2.push_back(path_edges[i]);

    vis(path_verts2, path_edges2);

    free(edge_betweenness);
    free(cumulative_edge_betweenness);

    return;
*/
  }

#ifdef DEBUG
  printf("num candidate matches: %lu\n", num_candidates);
  fflush(stdout);
#endif

  // Apply the visitor to each candidate match.  The candidate matches are
  // random walks through the colored graph where the choice of the next edge
  // in the walk is weighted based on edge betweenness.
  compute_weights_for_random_walk(colored_graph, eweight, cumulative_eweight);

#ifdef DEBUG
  edge_iterator_b edgs = edges(colored_graph);
  for (size_type_b i = 0; i < c_size; ++i)
  {
    printf("eweight[%lu]: %lu\n", get(eid_map_c, edgs[i]), eweight[edgs[i]]);
  }

  for (size_type_b i = 0; i < c_size; ++i)
  {
    printf("cum_ebetweenness[%lu]: %lu\n", get(eid_map_c, edgs[i]),
           cumulative_edge_betweenness[get(eid_map_c, edgs[i])]);
  }

  for (size_type_b i = 0; i < c_size; ++i)
  {
    printf("cum_eweight[%lu]: %lu\n", get(eid_map_c, edgs[i]),
           cumulative_eweight[edgs[i]]);
  }
#endif

#ifdef __MTA__
  size_type_b this_stream = 0;
#endif
  size_type_b num_streams = 1;

  #pragma mta use 75 streams
  #pragma mta for all streams this_stream of num_streams
  {
  }

  // Generate all the random numbers needed for the random walks through all
  // candidate paths.
  size_type_b max_num_candidates_per_stream =
    (num_candidates + num_streams - 1) / num_streams;
  size_type_b edges_in_walk = levels_in_colored / 2 + 1;
  size_type_b num_rand_per_stream = max_num_candidates_per_stream *
                                    edges_in_walk;
  size_type_b num_rand = num_rand_per_stream * num_streams;

  lrand48_generator rvals(num_rand);

  #pragma mta assert parallel
  for (size_type_b i = 0; i < num_streams; ++i)
  {
    size_type_b beg = begin_block_range(num_candidates, i, num_streams);
    size_type_b end = end_block_range(num_candidates, i, num_streams);

    dynamic_array<vertex_descriptor_b> path_verts;
    dynamic_array<edge_descriptor_b> path_edges;

    assert(levels_in_colored - 4 <= MAX_LEVELS);

    dynamic_array<GlobalSizeType> new_path(levels_in_colored - 4);
    dynamic_array<vertex_descriptor>
      path_vertex_descs((levels_in_colored - 4) / 2 + 1);
    dynamic_array<edge_descriptor> path_edge_descs((levels_in_colored - 4) / 2);

    vertex_iterator verts = vertices(bigG);
    edge_iterator edgs = edges(bigG);

//    printf("thread %lu about to do %lu paths\n", i, end - beg + 1);
//    fflush(stdout);

    size_type_b rval_start = i * num_rand_per_stream;

    for (size_type_b j = beg; j < end; ++j)
    {
      weighted_random_walk(colored_graph, walk_src, edges_in_walk, path_verts,
                           path_edges, cumulative_eweight, rvals, rval_start);

      // Create the match to send to the match visitor.
      GlobalSizeType walk_id = c_vid_map[get(vid_map_c, path_verts[1])];
      new_path[0] = walk_id;

      size_type next_vertex = 0;
      size_type next_edge = 0;

      size_type dense_id = (std::numeric_limits<size_type>::max)();
      big_graph_dense_vertex_id.lookup(walk_id, dense_id);

      path_vertex_descs[next_vertex++] = verts[dense_id];

      // Don't include the first and last edges.  They are the ones we
      // added to connect to the added supersource and supersink.
      for (size_type_b k = 3; k < levels_in_colored - 3; k += 2)
      {
        // Convert to the bigG ids.
        walk_id = c_eid_map[get(eid_map_c, path_edges[k / 2])];
        new_path[k - 2] = walk_id;
        big_graph_dense_edge_id.lookup(walk_id, dense_id);
        edge_descriptor e = edgs[dense_id];
        path_edge_descs[next_edge++] = e;

        walk_id = c_vid_map[get(vid_map_c, path_verts[k / 2 + 1])];
        new_path[k - 1] = walk_id;
        big_graph_dense_vertex_id.lookup(walk_id, dense_id);
        vertex_descriptor v = verts[dense_id];
        path_vertex_descs[next_vertex++] = v;
      }

#ifdef DEBUG
      printf("Path:\n");
      for (typename dynamic_array<size_type_b>::size_type k = 0;
           k < new_path.size(); k += 2)
      {
        printf("%lu\n", new_path[k]);
      }

      printf("Vertices and edges:\n");
      typename dynamic_array<edge_descriptor>::size_type k = 0;
      for ( ; k < path_edge_descs.size(); ++k)
      {
        printf("v: %lu, e: %lu\n", get(vid_map, path_vertex_descs[k]),
               get(eid_map, path_edge_descs[k]));
      }
      printf("v: %lu\n", get(vid_map, path_vertex_descs[k]));
#endif

      vis(path_vertex_descs, path_edge_descs);
    }

//    printf("thread %lu done\n", i);
//    fflush(stdout);
  }

  free(edge_betweenness);
  free(cumulative_edge_betweenness);
}

template <typename SizeType>
class cv_comp {
public:
  cv_comp(SizeType* cvm) : c_vid_map(cvm) {};

  bool operator()(SizeType x, SizeType y)
  {
    return c_vid_map[x] < c_vid_map[y];
  }

private:
  SizeType* c_vid_map;
};

template <typename Graph, typename GlobalSizeType>
void
find_like_pairs(
  typename graph_traits<Graph>::size_type first_level,
  typename graph_traits<Graph>::size_type second_level,
  Graph& colored_graph, GlobalSizeType* c_vid_map,
  typename graph_traits<Graph>::size_type* colored_vertex_level,
  std::vector<pair<typename graph_traits<Graph>::vertex_descriptor,
                   typename graph_traits<Graph>::vertex_descriptor> >& pairs)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  assert(first_level != second_level);

  // Need level counts for efficiency.
  // STL for now (this is not for the XMT!)

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, colored_graph);
  std::vector<vertex_descriptor> first_level_vertices, second_level_vertices;

  vertex_iterator verts = vertices(colored_graph);
  size_type order = num_vertices(colored_graph);

  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = verts[i];
    size_type c_id = get(vid_map, v);

    if (colored_vertex_level[c_id] == first_level)
    {
      first_level_vertices.push_back(v);
    }
    else if (colored_vertex_level[c_id] == second_level)
    {
      second_level_vertices.push_back(v);
    }
  }

  cv_comp<GlobalSizeType> cvc(c_vid_map);

// first_level_vertices.sort(cvc);
// second_level_vertices.sort(cvc);

  std::sort(first_level_vertices.begin(), first_level_vertices.end(), cvc);
  std::sort(second_level_vertices.begin(), second_level_vertices.end(), cvc);

  // Now, walk the lists picking out pairs that correspond to the same
  // vertex in the original graph.
  typename std::vector<vertex_descriptor>::iterator lev1_it =
    first_level_vertices.begin();
  typename std::vector<vertex_descriptor>::iterator lev2_it =
    second_level_vertices.begin();
  typename std::vector<vertex_descriptor>::iterator lev1_end =
    first_level_vertices.end();
  typename std::vector<vertex_descriptor>::iterator lev2_end =
    second_level_vertices.end();

#ifdef DEBUG
  printf("find_pairs: level1 vertices:\n");

  for ( ; lev1_it != lev1_end; ++lev1_it)
  {
    vertex_descriptor lev1_v = *lev1_it;
    printf("%lu\n", get(vid_map, lev1_v));
  }

  printf("find_pairs: level2 vertices:\n");

  for ( ; lev2_it != lev2_end; ++lev2_it)
  {
    vertex_descriptor lev2_v = *lev2_it;
    printf("%lu\n", get(vid_map, lev2_v));
  }

  lev1_it = first_level_vertices.begin();
  lev2_it = second_level_vertices.begin();
  lev1_end = first_level_vertices.end();
  lev2_end = second_level_vertices.end();
#endif

  while (lev1_it != lev1_end && lev2_it != lev2_end)
  {
    vertex_descriptor lev1_v = *lev1_it;
    vertex_descriptor lev2_v = *lev2_it;
    size_type l1v_id = get(vid_map, lev1_v);
    size_type l2v_id = get(vid_map, lev2_v);

#ifdef DEBUG
    printf("find_pairs: testing %lu [%lu] and %lu [%lu]\n",
           l1v_id, c_vid_map[l1v_id], l2v_id, c_vid_map[l2v_id]);
#endif

    if (c_vid_map[l1v_id] == c_vid_map[l2v_id])
    {
      pair<vertex_descriptor, vertex_descriptor> p(lev1_v, lev2_v);
      pairs.push_back(p);
      ++lev1_it;
      ++lev2_it;
    }
    else if (c_vid_map[l1v_id] < c_vid_map[l2v_id])
    {
      ++lev1_it;
    }
    else
    {
      ++lev2_it;
    }
  }
}

}

/*! \brief The walk-based subgraph isomorphism algorithm of Berry,
           Hendrickson, Kahan, and Konecny, Cray User Group 2006.

    \param bigG The big graph that is to be searched.
    \param targetG The target or pattern graph that is to be matched in the
                   big graph.
    \param walk_verts Ordered vertex descriptors of the walk through the
                      target graph.
    \param walk_edges Ordered edgs descriptors of the walk through the
                      target graph.
    \param si_compare The function object comparator used to customize the
                      matching.
    \param vis The visitor to apply to matches to the walk.
*/
template <typename Graph, typename GraphTrg,
          typename SIComparator, typename MatchVisitor>
typename graph_traits<Graph>::size_type
subgraph_isomorphism(
  Graph& bigG, GraphTrg& targetG,
  dynamic_array<typename graph_traits<GraphTrg>::vertex_descriptor>& walk_verts,
  dynamic_array<typename graph_traits<GraphTrg>::edge_descriptor>& walk_edges,
  SIComparator& si_compare, MatchVisitor& vis,
  dynamic_array<pair<typename graph_traits<Graph>::size_type,
                     typename graph_traits<Graph>::size_type> > lev_pairs =
    dynamic_array<pair<typename graph_traits<Graph>::size_type,
                       typename graph_traits<Graph>::size_type> >())
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  typedef typename graph_traits<GraphTrg>::vertex_descriptor vertex_trg;
  typedef typename graph_traits<GraphTrg>::edge_descriptor edge_trg;
  typedef typename graph_traits<GraphTrg>::size_type size_type_trg;
  typedef typename graph_traits<GraphTrg>::vertex_iterator vertex_iterator_trg;
  typedef typename graph_traits<GraphTrg>::edge_iterator edge_iterator_trg;

  typedef compressed_sparse_row_graph<directedS> bipartite_adapter;
  typedef typename graph_traits<bipartite_adapter>::size_type size_type_b;
  typedef typename graph_traits<bipartite_adapter>::vertex_descriptor
          vertex_descriptor_b;
  typedef typename graph_traits<bipartite_adapter>::edge_descriptor
          edge_descriptor_b;
  typedef typename graph_traits<bipartite_adapter>::vertex_iterator
          vertex_iterator_b;
  typedef typename graph_traits<bipartite_adapter>::edge_iterator
          edge_iterator_b;
  typedef typename graph_traits<bipartite_adapter>::out_edge_iterator
          out_edge_iterator_b;
  typedef array_property_map<size_type_b, vertex_id_map<bipartite_adapter> >
          vertex_property_map_b;
  typedef array_property_map<size_type_b, edge_id_map<bipartite_adapter> >
          edge_property_map_b;
  typedef vertex_property_map<bipartite_adapter, size_type_b> vpropmap_b;

#ifdef DEBUG
  printf("Subgraph_isomorphism called.\n");
  fflush(stdout);
#endif

  // Get the vertex and edge id maps.
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, bigG);
  vertex_id_map<GraphTrg> vid_map_trg = get(_vertex_id_map, targetG);
  edge_id_map<Graph> eid_map = get(_edge_id_map, bigG);
  edge_id_map<GraphTrg> eid_map_trg = get(_edge_id_map, targetG);

  size_type bigG_order = num_vertices(bigG);
  size_type bigG_size = num_edges(bigG);

  // **************************************************************************
  // The vertex id's might not be dense.  We need a flag like directedS to
  // indicate is.  In the absence of one, we'll assume that we need to create
  // dense vertex ids.
  // orig --> dense     :     use hash table big_graph_dense_vertex_id
  // dense --> original :     big_graph_original_id[dense_id]
  // **************************************************************************
  xmt_hash_table<size_type, size_type>
    big_graph_dense_vertex_id(static_cast<size_type_b>(2.2 * bigG_order));
  xmt_hash_table<size_type, size_type>
    big_graph_dense_edge_id(static_cast<size_type_b>(2.2 * bigG_size));

  vertex_iterator bgverts = vertices(bigG);

  #pragma mta assert parallel
  for (size_type i = 0; i < bigG_order; ++i)
  {
    vertex_descriptor v = bgverts[i];
    size_type vid = get(vid_map, v);
    big_graph_dense_vertex_id.insert(vid, i);
  }

  edge_iterator bgedges = edges(bigG);

  #pragma mta assert parallel
  for (size_type i = 0; i < bigG_size; ++i)
  {
    edge_descriptor e = bgedges[i];
    size_type eid = get(eid_map, e);
    big_graph_dense_edge_id.insert(eid, i);
  }

//  assign_contiguous_ids(big_graph_dense_vertex_id);

#ifdef DEBUG
  printf("done setting dense ids for big graph\n");
  fflush(stdout);
#endif

  size_type* big_graph_original_id =
    (size_type*) malloc(sizeof(size_type) * bigG_order);

  #pragma mta assert parallel
  for (size_type i = 0; i < bigG_order; ++i)
  {
    vertex_descriptor v = bgverts[i];
    size_type orig_id = get(vid_map, v);
    size_type dense_id = (std::numeric_limits<size_type>::max)();
    big_graph_dense_vertex_id.lookup(orig_id, dense_id);
    big_graph_original_id[dense_id] = orig_id;

#ifdef DEBUG
    printf("sparse->dense mapping of big graph vertex: %lu %lu\n",
           orig_id, dense_id);
    fflush(stdout);
#endif
  }
  // **************************************************************************

  // Get memory for and initialize active_verts which is a representation of
  // the vertices in the bipartite graph.
  size_type num_walk_verts = walk_verts.size();
  bool** active_verts = new bool*[num_walk_verts];
  #pragma mta assert nodep
  for (size_type_trg i = 0; i < num_walk_verts; ++i)
  {
    active_verts[i] = new bool[bigG_order];
  }

  for (size_type_trg i = 0; i < num_walk_verts; ++i)
  {
    for (size_type j = 0; j < bigG_order; ++j) active_verts[i][j] = false;
  }

  // 1. Initialize active_verts and count the bipartite edges for the first
  // walk edge.

  size_type_b b_size = 0;

  edge_iterator bg_edges = edges(bigG);
  edge_iterator_trg tg_edges = edges(targetG);
  vertex_iterator_trg tg_verts = vertices(targetG);

  vertex_trg trg_walk_src = walk_verts[0];
  edge_trg trg_e = walk_edges[0];

  // Determine if the target walk edge is in the forward or reversed
  // direction.  If the source vertex in the walk is the source of the
  // target edge, then the target edge is a forward edge.  Otherwise, the
  // target edge is a backward edge.
  size_type_trg trg_walk_src_id = get(vid_map_trg, trg_walk_src);
  size_type_trg trg_src_id = get(vid_map_trg, source(trg_e, targetG));
  bool trg_e_reversed = trg_walk_src_id != trg_src_id;

  // Add the unique vertices to the bipartite graph.
  #pragma mta assert parallel
  for (size_type i = 0; i < bigG_size; ++i)
  {
    edge_descriptor e = bg_edges[i];

    // Get the "source" and "target" vertices for the big graph edge.  For
    // forward target edges, the "source" of the big graph edge is the edge's
    // actual source.  For backward target edges, the "source" of the big
    // graph edge is the edge's target.
    vertex_descriptor bg_src = trg_e_reversed ? target(e, bigG) :
                               source(e, bigG);
    vertex_descriptor bg_dest = bg_src == source(e, bigG) ? target(e, bigG) :
                                source(e, bigG);

    // Get the vertex ids.
    size_type bg_src_id = get(vid_map, bg_src);
    size_type bg_dest_id = get(vid_map, bg_dest);

    // Get the dense vertex ids.
    size_type bg_dense_src_id = 0;
    size_type bg_dense_dest_id = 0;
    big_graph_dense_vertex_id.lookup(bg_src_id, bg_dense_src_id);
    big_graph_dense_vertex_id.lookup(bg_dest_id, bg_dense_dest_id);

    // Check if the first walk edge matches this edge from the big graph.
    if (si_compare(e, trg_e, bigG, targetG, 0))
    {
      mt_incr(b_size, 1);

      active_verts[0][bg_dense_src_id] = true;
      active_verts[1][bg_dense_dest_id] = true;
    }

    // If the graph is undirected, also check the reverse of this edge.
    if (is_undirected(bigG))
    {
      edge_descriptor e_rev =
        edge_descriptor(target(e, bigG), source(e, bigG), get(eid_map, e));

      if (si_compare(e_rev, trg_e, bigG, targetG, 0))
      {
        mt_incr(b_size, 1);

        active_verts[0][bg_dense_dest_id] = true;
        active_verts[1][bg_dense_src_id] = true;
      }
    }
  }

#ifdef DEBUG
  printf("subgraph_isomorphism: the active vertices at Level 0: \n");
  printf("------------------------------------------------------\n");
  fflush(stdout);

  for (size_type i = 0; i < bigG_size; ++i)
  {
    edge_descriptor e = bg_edges[i];
    vertex_descriptor bg_src = trg_e_reversed ? target(e, bigG) :
                               source(e, bigG);
    vertex_descriptor bg_dest = bg_src == source(e, bigG) ? target(e, bigG) :
                                source(e, bigG);

    size_type bg_src_id = get(vid_map, bg_src);
    size_type bg_dest_id = get(vid_map, bg_dest);
    size_type bg_dense_src_id, bg_dense_dest_id;

    big_graph_dense_vertex_id.lookup(bg_src_id, bg_dense_src_id);
    big_graph_dense_vertex_id.lookup(bg_dest_id, bg_dense_dest_id);

    if ((active_verts[0][bg_dense_src_id]) &&
        (active_verts[1][bg_dense_dest_id]))
    {
      printf("%lu  sparse (%lu, %lu) dense (%lu, %lu) bip_id (%lu, %lu) \n",
             bg_src_id, bg_src_id, bg_dest_id,
             bg_dense_src_id, bg_dense_dest_id,
             bg_dense_src_id, bigG_order * 1 + bg_dense_dest_id);
    }

    if ((active_verts[0][bg_dense_dest_id]) &&
        (active_verts[1][bg_dense_src_id]))
    {
      printf("%lu  sparse (%lu, %lu) dense (%lu, %lu) bip_id (%lu, %lu) \n",
             bg_dest_id, bg_dest_id, bg_src_id,
             bg_dense_dest_id, bg_dense_src_id,
             bg_dense_dest_id, bigG_order * 1 + bg_dense_src_id);
    }
  }
#endif

  // 2. Initialize active_verts and count the bipartite edges for the remaining
  // walk edges.

  size_type_b active_level = 2;  // Level in bipartite graph.

  // Try to match the remainder of the walk through the target graph against
  // the filtered graph.
  for (size_type_trg j = 1; j < walk_edges.size(); ++j)
  {
    vertex_trg trg_walk_src = walk_verts[j];
    edge_trg trg_e = walk_edges[j];

    // Determine if the target walk edge is in the forward or reversed
    // direction.  If the source vertex in the walk is the source of the
    // target edge, then the target edge is a forward edge.  Otherwise, the
    // target edge is a backward edge.
    size_type_trg trg_walk_src_id = get(vid_map_trg, trg_walk_src);
    size_type_trg trg_src_id = get(vid_map_trg, source(trg_e, targetG));
    bool trg_e_reversed = trg_walk_src_id != trg_src_id;

    // Add the unique vertices to the bipartite graph.
    #pragma mta assert parallel
    for (size_type i = 0; i < bigG_size; ++i)
    {
      edge_descriptor e = bg_edges[i];

      // Get the "source" and "target" vertices for the big graph edge.  For
      // forward target edges, the "source" of the big graph edge is the edge's
      // actual source.  For backward target edges, the "source" of the big
      // graph edge is the edge's target.
      vertex_descriptor bg_src = trg_e_reversed ? target(e, bigG) :
                                 source(e, bigG);
      vertex_descriptor bg_dest = bg_src == source(e, bigG) ? target(e, bigG) :
                                  source(e, bigG);

      // Get the vertex ids.
      size_type bg_src_id = get(vid_map, bg_src);
      size_type bg_dest_id = get(vid_map, bg_dest);

      // Get the dense vertex ids.
      size_type bg_dense_src_id = 0;
      size_type bg_dense_dest_id = 0;
      big_graph_dense_vertex_id.lookup(bg_src_id, bg_dense_src_id);
      big_graph_dense_vertex_id.lookup(bg_dest_id, bg_dense_dest_id);

      // Check if the current walk edge matches this edge from the big graph.
      if (active_verts[active_level - 1][bg_dense_src_id] &&
          si_compare(e, trg_e, bigG, targetG, active_level - 1))
      {
        mt_incr(b_size, 1);

        active_verts[active_level][bg_dense_dest_id] = true;
      }

      // If the graph is undirected, also check the reverse of this edge.
      if (is_undirected(bigG))
      {
        edge_descriptor e_rev =
          edge_descriptor(target(e, bigG), source(e, bigG), get(eid_map, e));

        if (active_verts[active_level - 1][bg_dense_dest_id] &&
            si_compare(e_rev, trg_e, bigG, targetG, active_level - 1))
        {
          mt_incr(b_size, 1);

          active_verts[active_level][bg_dense_src_id] = true;
        }
      }
    }

    ++active_level;
  }

#ifdef DEBUG
  for (size_type lev = 1; lev < active_level; ++lev)
  {
    printf("subgraph_isomorphism: the active vertices and edges at Level "
           "%lu: \n", lev);
    printf("------------------------------------------------------\n");
    fflush(stdout);

    for (size_type i = 0; i < bigG_size; ++i)
    {
      edge_descriptor e = bg_edges[i];
      vertex_descriptor bg_src = trg_e_reversed ? target(e, bigG) :
                                 source(e, bigG);
      vertex_descriptor bg_dest = bg_src == source(e, bigG) ? target(e, bigG) :
                                  source(e, bigG);

      size_type bg_src_id = get(vid_map, bg_src);
      size_type bg_dest_id = get(vid_map, bg_dest);
      size_type bg_dense_src_id, bg_dense_dest_id;

      big_graph_dense_vertex_id.lookup(bg_src_id, bg_dense_src_id);
      big_graph_dense_vertex_id.lookup(bg_dest_id, bg_dense_dest_id);

      if ((active_verts[lev - 1][bg_dense_src_id]) &&
          (active_verts[lev][bg_dense_dest_id]))
      {
        printf("%lu  sparse (%lu, %lu) dense (%lu, %lu) bip_id (%lu, %lu) \n",
               bg_src_id, bg_src_id, bg_dest_id,
               bg_dense_src_id, bg_dense_dest_id,
               bigG_order * (lev - 1) + bg_dense_src_id,
               bigG_order * lev + bg_dense_dest_id);
      }

      if ((active_verts[lev - 1][bg_dense_dest_id]) &&
          (active_verts[lev][bg_dense_src_id]))
      {
        printf("%lu  sparse (%lu, %lu) dense (%lu, %lu) bip_id (%lu, %lu) \n",
               bg_dest_id, bg_dest_id, bg_src_id,
               bg_dense_dest_id, bg_dense_src_id,
               bigG_order * (lev - 1) + bg_dense_dest_id,
               bigG_order * lev + bg_dense_src_id);
      }
    }
  }
#endif

  // Create the data structures that will hold the bipartite vertex and edge
  // information.  We know the exact number of edges and estimate the number
  // of vertices as 2.2 * # edges.  There can't be more than 2 * # edges, and
  // the .2 guarantees the hash table has some empty space.

  // The hash table bipartite_vertex_id stores a mapping from the spaced out
  // ids originally assigned to the bipartite graph to the compressed ids
  // later assigned to the bipartite graph.  The purpose of the hash table is
  // to create a set of the unique vertices in the bipartite graph and then
  // create compressed ids for those vertices.
  xmt_hash_table<size_type_b, size_type_b>
    bipartite_vertex_id(static_cast<size_type_b>(2.2 * b_size));

  // The arrays b_sources and b_dests will initially hold a mapping back to
  // the vertex ids in the original graph, but these will be replaced by the
  // compressed unique ids for the bipartite graph.  The array b_eid_map will
  // hold a mapping back to the edge ids in the original graph.  You get an
  // original edge id for edge i in the bipartite graph by b_eid_map[i].
  dynamic_array<size_type_b> b_sources(b_size);
  dynamic_array<size_type_b> b_dests(b_size);
  dynamic_array<size_type_b> b_eid_map(b_size);

  bipartite_adapter bipartite_graph;
  detail::assign_sparse_bipartite_ids(bigG, bipartite_graph, targetG,
                                      active_verts, walk_verts, walk_edges,
                                      b_size, b_sources, b_dests, b_eid_map,
                                      bipartite_vertex_id,
                                      big_graph_dense_vertex_id,
                                      si_compare);


  #pragma mta assert nodep
  for (size_type_trg i = 0; i < num_walk_verts; ++i) delete [] active_verts[i];
  delete [] active_verts;

  size_type_b b_order = bipartite_vertex_id.size();

  // The array b_vid_map will hold a mapping back to the vertex ids in the
  // original graph.  You get an original vertex id for vertex i in the
  // bipartite graph by b_vid_map[i].
  size_type* b_vid_map = (size_type*) malloc((b_order + 2) * sizeof(size_type));

  // The array bipartite_vertex_level is indexed by the compressed ids of the
  // bipartite graph.
  size_type_b* bipartite_vertex_level =
    (size_type_b*) malloc((b_order + 2) * sizeof(size_type_b));

  // Create the compressed unique ids for the bipartite vertices; set their
  // mappings back to the vertices in bigG; and set their levels.
  assign_contiguous_ids(bipartite_vertex_id);
  detail::bipartite_vertex_visitor<xmt_hash_table<size_type_b, size_type_b>,
                                   size_type>
  bv_vis(bigG_order, bipartite_vertex_level,
         big_graph_original_id, b_vid_map);

  bipartite_vertex_id.visit(bv_vis);
  free(big_graph_original_id);

  // Replace the spaced out ids of b_sources and b_dests with the compressed
  // ids.
  #pragma mta assert nodep
  #pragma mta assert parallel
  for (size_type_b i = 0; i < b_size; ++i)
  {
    size_type_b bsi = 0;
    size_type_b bdi = 0;

    bipartite_vertex_id.lookup(b_sources[i], bsi);
    bipartite_vertex_id.lookup(b_dests[i], bdi);

#ifdef DEBUG
    printf("%lu'th sparse bip edge: (%lu, %lu): dense bip: (%lu, %lu)\n",
           i, b_sources[i], b_dests[i], bsi, bdi);
#endif

    b_sources[i] = bsi;
    b_dests[i] = bdi;
  }

  // Add entries for the new supersource and supersink in the enhanced
  // bipartite graph.
  bipartite_vertex_level[b_order] = active_level;
  bipartite_vertex_level[b_order + 1] =
    (std::numeric_limits<size_type_b>::max)();
  b_vid_map[b_order] = (std::numeric_limits<size_type>::max)();
  b_vid_map[b_order + 1] = (std::numeric_limits<size_type>::max)();

  // Find the number of max level vertices.
  size_type_b num_max_level_verts = 0;
  #pragma mta assert parallel
  for (size_type_b i = 0; i < b_order; ++i)
  {
    mt_incr(num_max_level_verts,
            (bipartite_vertex_level[i] == active_level - 1));
  }

  if (num_max_level_verts == 0)
  {
    free(bipartite_vertex_level);
    free(b_vid_map);
    printf("subgraph_isomorphism: no candidates found\n");
    return 0;
  }

  size_type_b* max_level_vertices =
    (size_type_b*) malloc(num_max_level_verts * sizeof(size_type_b));

  // Fill the max level vertices array.
  num_max_level_verts = 0;

  #pragma mta assert nodep
  #pragma mta assert parallel
  for (size_type_b i = 0; i < b_order; ++i)
  {
    if (bipartite_vertex_level[i] == active_level - 1)
    {
      size_type_b pos = mt_incr(num_max_level_verts, 1);
      max_level_vertices[pos] = i;
    }
  }

  // Find the number of zero level vertices.
  size_type_b num_zero_level_verts = 0;

  #pragma mta assert parallel
  for (size_type_b i = 0; i < b_order; ++i)
  {
    mt_incr(num_zero_level_verts, (bipartite_vertex_level[i] == 0));
  }

  size_type_b* zero_level_vertices =
    (size_type_b*) malloc(num_zero_level_verts * sizeof(size_type_b));

  // Fill the zero level vertices array.
  num_zero_level_verts = 0;

  #pragma mta assert parallel
  for (size_type_b i = 0; i < b_order; ++i)
  {
    if (bipartite_vertex_level[i] == 0)
    {
      size_type_b pos = mt_incr(num_zero_level_verts, 1);
      zero_level_vertices[pos] = i;
    }
  }

  // Get the number of vertices and edges for the enhanced biparatite graph.
  size_type_b b_order2 = b_order + 2;
  size_type_b b_size2 = b_size + num_max_level_verts + num_zero_level_verts;

  b_sources.resize(b_size2);
  b_dests.resize(b_size2);
  b_eid_map.resize(b_size2);

  // Add the edges going from the max level vertices to the sink.
  #pragma mta assert parallel
  for (size_type_b i = 0; i < num_max_level_verts; ++i)
  {
    b_sources[b_size + i] = max_level_vertices[i];
    b_dests[b_size + i] = b_order;
    b_eid_map[b_size + i] = (std::numeric_limits<size_type>::max)();
  }

  // Add the edges going from the source to the zero level vertices.
  #pragma mta assert nodep
  for (size_type_b i = 0; i < num_zero_level_verts; ++i)
  {
    b_sources[b_size + num_max_level_verts + i] = b_order + 1;
    b_dests[b_size + num_max_level_verts + i] = zero_level_vertices[i];
    b_eid_map[b_size + num_max_level_verts + i] =
      (std::numeric_limits<size_type>::max)();
  }

  free(max_level_vertices);
  free(zero_level_vertices);

#ifdef DEBUG
  printf("initializing bipartite graph with %lu vertices and %lu edges \n",
         b_order2, b_size2);
  printf("-------------------------------------------------------------- \n");

  for (size_type_b i = 0; i < b_size2; ++i)
  {
    printf("%lu %lu\n", b_dests.get_data()[i], b_sources.get_data()[i]);
    fflush(stdout);
  }
#endif

  // Flip the sources and destinations because we want the directed edges
  // going from the end of the match to the beginning.  That way we don't have
  // to worry about levels.  Moving forward an edge in the graph moves forward
  // one level.
  init(b_order2, b_size2, b_dests.get_data(), b_sources.get_data(),
       bipartite_graph);

#ifdef DEBUG
  printf("bipartite graph: \n");
  printf("---------------- \n");
  fflush(stdout);
  print(bipartite_graph);
  fflush(stdout);

  vertex_iterator_b bverts = vertices(bipartite_graph);
  vertex_id_map<bipartite_adapter> b_vmap = get(_vertex_id_map,
                                                bipartite_graph);
  printf("bipartite graph vertex mappings \n");
  printf("------------------------------- \n");
  fflush(stdout);

  for (size_type_b i = 0; i < num_vertices(bipartite_graph); ++i)
  {
    size_type_b bv = get(b_vmap, bverts[i]);
    printf("%lu  --->  %lu\n", bv, b_vid_map[bv]);
  }

  edge_iterator_b bedgs = edges(bipartite_graph);
  edge_id_map<bipartite_adapter> b_emap = get(_edge_id_map, bipartite_graph);

  printf("bipartite graph edge mappings \n");
  printf("------------------------------- \n");
  fflush(stdout);

  for (size_type_b i = 0; i < num_edges(bipartite_graph); ++i)
  {
    size_type_b be = get(b_emap, bedgs[i]);
    printf("%lu  --->  %lu\n", be, b_eid_map[be]);
  }

  fflush(stdout);
#endif

  size_type_b* b_edge_color =
    (size_type_b*) calloc(b_size2, sizeof(size_type_b));

  int* search_color = (int*) calloc(b_order2, sizeof(int));

  detail::iso_bipartite_visitor<bipartite_adapter> ibv(bipartite_graph,
                                                       b_edge_color);

  psearch<bipartite_adapter, int*,
          detail::iso_bipartite_visitor<bipartite_adapter>,
          AND_FILTER, DIRECTED>
    psrch(bipartite_graph, search_color, ibv);

  psrch.run(vertices(bipartite_graph)[b_order]);

  free(search_color);

  // Count the number of colored edges.  Consider even the edges added for the
  // source and sink.
  #pragma mta assert parallel
  size_type_b num_colored_edges = 0;

  for (size_type_b i = 0; i < b_size2; ++i)
  {
    mt_incr(num_colored_edges, b_edge_color[i]);
  }

  // Create the data structures that will hold the colored vertex and edge
  // information.  We know the exact number of edges and estimate the number
  // of vertices as 2.2 * # edges.  There can't be more than 2 * # edges, and
  // the .2 guarantees the hash table has some empty space.

  // The hash table colored_vertex_id stores a mapping from the spaced out
  // ids originally assigned to the colored graph (ids from the bipartite
  // graph) to the compressed ids later assigned to the colored graph.  The
  // purpose of the hash table is to create a set of the unique vertices in
  // the colored graph and then create compressed ids for those vertices.
  xmt_hash_table<size_type_b, size_type_b>
    colored_vertex_id(static_cast<size_type_b>(2.2 * num_colored_edges));

  // The arrays c_sources and c_dests will initially hold a mapping back to
  // the vertex ids in the bipartite graph, but these will be replaced by the
  // compressed unique ids for the colored graph.  The array c_eid_map will
  // hold a mapping back to the edge ids in the original graph.  You get an
  // original edge id for edge i in the colored graph by c_eid_map[i].
  dynamic_array<size_type_b> c_sources(num_colored_edges);
  dynamic_array<size_type_b> c_dests(num_colored_edges);
  dynamic_array<size_type> c_eid_map(num_colored_edges);

  // Create list of colored edges.
  num_colored_edges = 0;

  #pragma mta assert parallel
  for (size_type_b i = 0; i < b_size2; ++i)
  {
    if (b_edge_color[i])
    {
      size_type_b c_eid = mt_incr(num_colored_edges, 1);

      c_sources[c_eid] = b_sources[i];
      c_dests[c_eid] = b_dests[i];
      c_eid_map[c_eid] = b_eid_map[i];

      colored_vertex_id.insert(b_sources[i], 0);
      colored_vertex_id.insert(b_dests[i], 0);
    }
  }

  free(b_edge_color);

  size_type_b c_order = colored_vertex_id.size();
  size_type_b c_size = num_colored_edges;

  // The array c_vid_map will hold a mapping back to the vertex ids in the
  // original graph.  You get an original vertex id for vertex i in the
  // colored graph by c_vid_map[i].
  size_type* c_vid_map = (size_type*) malloc(c_order * sizeof(size_type));

  // The array colored_vertex_level is indexed by the compressed ids of the
  // colored graph.
  size_type_b* colored_vertex_level =
    (size_type_b*) malloc(c_order * sizeof(size_type_b));

  // Create the compressed unique ids for the colored vertices, and set their
  // mappings back to the vertices in bigG.
  assign_contiguous_ids(colored_vertex_id);

  detail::colored_vertex_visitor<xmt_hash_table<size_type_b, size_type_b>,
                                 size_type, size_type_b>
    cv_vis(b_vid_map, c_vid_map, bipartite_vertex_level, colored_vertex_level);

  colored_vertex_id.visit(cv_vis);

  // Replace the spaced out ids of c_sources and c_dests with the compressed
  // ids.
  #pragma mta assert parallel
  for (size_type_b i = 0; i < c_size; ++i)
  {
    size_type_b csi = 0;
    size_type_b cdi = 0;

    colored_vertex_id.lookup(c_sources[i], csi);
    colored_vertex_id.lookup(c_dests[i], cdi);

    c_sources[i] = csi;
    c_dests[i] = cdi;
  }

  // Create the forward and reverse graphs.  This isn't more expensive in
  // space than creating a bipartite graph, and I'm not sure if there is a
  // BFS implementation that will work on reverse edges.
  size_type_b levels_in_colored = active_level + 1; // Add source and sink.
  bipartite_adapter colored_graph;
  init(c_order, c_size, c_sources.get_data(), c_dests.get_data(),
       colored_graph);

  vertex_id_map<bipartite_adapter> c_vmap = get(_vertex_id_map, colored_graph);
  edge_id_map<bipartite_adapter> c_emap = get(_edge_id_map, colored_graph);

#ifdef DEBUG
  printf("Colored graph:\n");
  print(colored_graph);

  vertex_iterator_b cverts = vertices(colored_graph);

  for (size_type_b i = 0; i < c_order; ++i)
  {
    size_type_b cid = get(c_vmap, cverts[i]);
    size_type orig_id = c_vid_map[cid];
    size_type c_lev = colored_vertex_level[cid];
    printf("colored_graph_id: %lu, orig_id: %lu, colored_level: %lu\n",
           cid, orig_id, c_lev);
  }
#endif

  bipartite_adapter colored_graph_rev;
  init(c_order, c_size, c_dests.get_data(), c_sources.get_data(),
       colored_graph_rev);

#ifdef DEBUG
  printf("Colored graph reversed:\n");
  print(colored_graph_rev);

  edge_iterator_b cedges = edges(colored_graph);
  edge_iterator_b credges = edges(colored_graph_rev);
  vertex_id_map<bipartite_adapter> vid_map1 =
    get(_vertex_id_map, colored_graph);
  vertex_id_map<bipartite_adapter> vid_map2 =
    get(_vertex_id_map, colored_graph_rev);
  edge_id_map<bipartite_adapter> eid_map1 = get(_edge_id_map, colored_graph);
  edge_id_map<bipartite_adapter> eid_map2 =
    get(_edge_id_map, colored_graph_rev);

  bool rev_correct = true;

  for (size_type i = 0; i < num_edges(colored_graph); ++i)
  {
    edge_descriptor_b c_edge = cedges[i];
    edge_descriptor_b cr_edge = credges[i];

    printf("colored: {%lu}: (%lu,%lu), rev: {%lu} (%lu, %lu)\n",
           get(eid_map1, c_edge),
           get(vid_map1, source(c_edge, colored_graph)),
           get(vid_map1, target(c_edge, colored_graph)),
           get(eid_map2, cr_edge),
           get(vid_map2, source(cr_edge, colored_graph_rev)),
           get(vid_map2, target(cr_edge, colored_graph_rev)));

    if ((get(vid_map1, source(c_edge, colored_graph)) !=
        get(vid_map2, target(cr_edge, colored_graph_rev))) ||
        get(vid_map1, target(c_edge, colored_graph)) !=
        get(vid_map2, source(cr_edge, colored_graph_rev)))
    {
      rev_correct = false;
      printf("OOOPS!\n");
    }
  }
#endif

  size_type_b cv_src = 0;
  size_type_b cv_trg = 0;
  colored_vertex_id.lookup(b_order + 1, cv_src);
  colored_vertex_id.lookup(b_order, cv_trg);

  // Create a filter to send to connected components.  We'll filter out the
  // edges adjacent to the super-source and super-sink.
  edge_iterator_b c_edges = edges(colored_graph);
  edge_iterator_b cr_edges = edges(colored_graph_rev);
  edge_property_map<bipartite_adapter, bool> c_internal_edge(colored_graph);

  for (size_type_b i = 0; i < c_size; ++i)
  {
    edge_descriptor_b e = c_edges[i];
    vertex_descriptor_b u = source(e, colored_graph);
    vertex_descriptor_b v = target(e, colored_graph);
    c_internal_edge[e] = !((get(c_vmap, u) == cv_src) ||
                           (get(c_vmap, v) == cv_trg));

#ifdef DEBUG
    printf("emask[(%lu, %lu)]: %d\n", get(c_vmap, u), get(c_vmap, v),
           c_internal_edge[e]);
#endif
  }

  size_type_b num_candidates = 0;

  vpropmap_b component(colored_graph);

  // That was kind of cheating, but we control the type of the colored graph.
  // So, we know that size_type is the same as vertex_descriptor_b.

  shiloach_vishkin(colored_graph, component, c_internal_edge);

  vertex_partition_iterator<bipartite_adapter, vpropmap_b>
  next_comp(colored_graph, component);

  size_type_b num_parts = next_comp.num_partitions();

  // Treat each connected component of the colored graph separately.
  for (size_type i = 0; i < num_parts; ++i)
  {
    pair<size_type_b, vertex_descriptor_b*> next = next_comp[i];

    // TODO: Get output into a nice format for subgraph formation,
    // where the each subgraph includes the super-source and super-sink,
    // yet uses array initialization rather than a O(n) mask.  Using
    // the latter leads to O(n^2) computation when there are O(n)
    // components.

    // It looks to be a hassle to avoid alloc's.  Doing them for now. JWB
    size_type_b comp_size = next.first;
    vertex_descriptor_b* comp_verts = next.second;

    // We filtered out super-source and super-sink for components
    // computation; they will appear alone and must not be considered
    // components on their own.
    if (!((comp_size == 1) &&
          ((comp_verts[0] == cv_src) || (comp_verts[0] == cv_trg))))
    {
      vertex_descriptor_b* sub_verts = (vertex_descriptor_b*)
        malloc(sizeof(vertex_descriptor_b) * (comp_size + 2));

      for (size_type i = 0; i < comp_size; ++i) sub_verts[i] = comp_verts[i];

      sub_verts[comp_size] = cv_src;
      sub_verts[comp_size + 1] = cv_trg;

      size_type s_size = comp_size + 2;

      // *****************************************************************
      // But since colored_graph and colored_graph_rev are separate,
      // we can't induce the subgraphs based on vertices and expect the
      // edge id's to be the same.  We have to construct arrays of edges
      // to pass to the subgraph_adapter initialization (init_edges).
      // Note that this would be unnecessary if we simply used an emask,
      // but we wish to avoid alloc'ing O(E) space.
      // *****************************************************************
      qsort(sub_verts, s_size, sizeof(vertex_descriptor_b), detail::uint_cmp);

      size_type num_s_edges = 0;

      for (size_type i = 0; i < s_size; ++i)
      {
        vertex_descriptor_b u = sub_verts[i];
        out_edge_iterator_b outedgs = out_edges(u, colored_graph);
        size_type deg = out_degree(u, colored_graph);

        for (size_type j = 0; j < deg; ++j)
        {
          edge_descriptor_b e = outedgs[j];
          vertex_descriptor_b v = target(e, colored_graph);
          vertex_descriptor_b* key = &v;

          vertex_descriptor_b* result =
            (vertex_descriptor_b*) bsearch(key, sub_verts, s_size,
                                           sizeof(vertex_descriptor_b),
                                           detail::uint_cmp);

          if (result) ++num_s_edges;
        }
      }

      edge_descriptor_b* sub_edges =
        (edge_descriptor_b*) malloc(sizeof(edge_descriptor_b) * num_s_edges);

      size_type next_e = 0;

      for (size_type i = 0; i < s_size; ++i)
      {
        vertex_descriptor_b u = sub_verts[i];
        out_edge_iterator_b outedgs = out_edges(u, colored_graph);
        size_type deg = out_degree(u, colored_graph);

        for (size_type j = 0; j < deg; ++j)
        {
          edge_descriptor_b e = outedgs[j];
          vertex_descriptor_b v = target(e, colored_graph);
          vertex_descriptor_b* key = &v;

          vertex_descriptor_b* result =
            (vertex_descriptor_b*) bsearch(key, sub_verts, s_size,
                                           sizeof(vertex_descriptor_b),
                                           detail::uint_cmp);

          if (result) sub_edges[next_e++] = e;
        }
      }

      // Now we "cheat" based on the knowledge that colored_graph and
      // colored_graph_rev are CSR graphs where size_type is the same as
      // vertex_descriptor.  We'll pass the same array of edge descriptors
      // to the subgraph initialization of both sub_colored and sub_colored_rev.
      // This wouldn't be correct in general.
      edge_descriptor_b* sub_edges_rev =
        (edge_descriptor_b*) malloc(sizeof(edge_descriptor_b) * num_s_edges);
      for (size_type i = 0; i < num_s_edges; ++i)
      {
        size_type eid = get(c_emap, sub_edges[i]);
        sub_edges_rev[i] = cr_edges[eid];
      }

      detail::process_subgraph_of_colored<Graph, bipartite_adapter, size_type,
                                          MatchVisitor, 0>
        psc(bigG, colored_graph, colored_graph_rev,
            big_graph_dense_vertex_id, big_graph_dense_edge_id, vis);

      psc(num_s_edges, sub_edges, sub_edges_rev, levels_in_colored,
          c_order, c_size, cv_src, cv_trg, c_vid_map, c_eid_map,
          num_candidates, colored_vertex_level, lev_pairs);

      free(sub_verts);
      free(sub_edges);
      free(sub_edges_rev);
    }
  }

/*
#ifdef DEBUG
  for (size_type_b i = 0; i < c_order; ++i)
  {
    printf("component[%lu]: %lu\n", i, component[i]);
  }
#endif

  // THIS WORKED BUT HAD O(n^2) WORST CASE BEHAVOIR.
  // if there are enough connected components,
  // partition by connected component and treat independently
  // naive impl.: worst case O(n^2); we may need to rewrite.  In that
  // case, need a method that will return an array of boolean masks,
  // or even a subgraph_adapter init method that takes a list rather
  // than a mask.  Haven't thought through the details. JWB
  dynamic_array<vertex_descriptor_b> leaders;
  component_leaders(colored_graph, component, leaders);
  size_type num_comp = leaders.size();

#ifdef DEBUG
  printf("number of components in colored graph %lu\n", num_comp);
#endif

  for (size_type i = 0; i < num_comp; ++i)
  {
    size_type next_leader = leaders[i];
    bool *emask = (bool*) calloc(sizeof(bool), c_size);
    size_type component_size = 0;

    for (size_type i = 0; i < c_order; ++i)
    {
      if (component[i] == next_leader) mt_incr(component_size, 1);
    }

    for (size_type i = 0; i < c_size; ++i)
    {
      edge_descriptor_b e = c_edges[i];
      vertex_descriptor_b u = source(e, colored_graph);
      vertex_descriptor_b v = target(e, colored_graph);
      size_type_b uid = get(c_vmap, u);
      size_type_b vid = get(c_vmap, v);

#ifdef DEBUG
      printf("c[%lu]: %lu, c[%lu]: %lu, cv_src: %lu, cv_trg: %lu\n",
             uid, component[uid], vid, component[vid],
             get(c_vmap, cv_src), get(c_vmap, cv_trg));
#endif

      if (((component[uid] == next_leader) &&
           (component[vid] == next_leader)) || (u == cv_src) || (v == cv_trg))
      {
        emask[i] = true;
      }
    }

#ifdef DEBUG
    printf("component leader %lu\n", next_leader);
    for (size_type i = 0; i < c_size; ++i)
    {
      printf("emask[%lu]: %d\n", i, emask[i]);
    }
#endif

    // At least one edge not adj. to super.
    if (component_size >= 2)
    {
      process_subgraph_of_colored<bipartite_adapter, MatchVisitor, 0>
        psc(colored_graph, colored_graph_rev, vis);

      psc(emask, levels_in_colored, c_order, c_size,
          cv_src, cv_trg, c_vid_map, c_eid_map, num_candidates,
          colored_vertex_level, lev_pair);

      // Need interface for specifying level constraints.
    }
  }
*/

  free(b_vid_map);
  free(c_vid_map);
  free(colored_vertex_level);
  free(bipartite_vertex_level);

  return num_candidates;
}

template <typename Graph, typename GraphTrg,
          typename SIComparator, typename MatchVisitor>
typename graph_traits<Graph>::size_type
subgraph_isomorphism(Graph& bigG, GraphTrg& targetG,
                     SIComparator& si_compare, MatchVisitor& vis)
{
  typedef typename graph_traits<GraphTrg>::size_type size_type_trg;
  typedef typename graph_traits<GraphTrg>::vertex_descriptor
          vertex_descriptor_trg;
  typedef typename graph_traits<GraphTrg>::edge_descriptor edge_descriptor_trg;

  dynamic_array<pair<size_type_trg, size_type_trg> > lev_pairs;
  dynamic_array<vertex_descriptor_trg> walk_verts;
  dynamic_array<edge_descriptor_trg> walk_edges;

  // Walk the graph.
  bool success = duplicate_euler_tour(targetG, vertices(targetG)[0],
                                      walk_verts, walk_edges);

  assert(success);

  // Find level pairs.
  vertex_id_map<GraphTrg> tvid_map = get(_vertex_id_map, targetG);
  edge_id_map<GraphTrg> teid_map = get(_edge_id_map, targetG);

  for (size_type_trg i = 0; i < walk_verts.size(); ++i)
  {
    vertex_descriptor_trg u = walk_verts[i];

    for (size_type_trg j = i + 1; j < walk_verts.size(); ++j)
    {
      vertex_descriptor_trg v = walk_verts[j];
      size_type_trg uid = get(tvid_map, u);
      size_type_trg vid = get(tvid_map, v);

      if (uid == vid)
      {
        lev_pairs.push_back(pair<size_type_trg, size_type_trg>(i, j));

        if (lev_pairs.size() >= (MAX_RECURSION - 1)) break;
      }
    }

    if (lev_pairs.size() >= (MAX_RECURSION - 1)) break;
  }

#ifdef DEBUG
  for (size_type_trg i = 0; i < walk_verts.size(); ++i)
  {
    printf("walk_verts[%lu]: %lu\n", i, get(tvid_map, walk_verts[i]));
  }

  for (size_type_trg i = 0; i < walk_edges.size(); ++i)
  {
    printf("walk_edges[%lu]: %lu\n", i, get(teid_map, walk_edges[i]));
  }

  for (size_type_trg i = 0; i < lev_pairs.size(); ++i)
  {
    pair<size_type_trg, size_type_trg> p = lev_pairs[i];

    printf("lev_pair: (%lu[%lu],%lu[%lu])\n",
           p.first, get(tvid_map, walk_verts[p.first]),
           p.second, get(tvid_map, walk_verts[p.second]));
  }
#endif

  return subgraph_isomorphism(bigG, targetG, walk_verts, walk_edges,
                              si_compare, vis, lev_pairs);
}

}

#undef MAX_LEVELS
#undef MAX_RECURSION

#endif
