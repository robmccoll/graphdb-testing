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
/*! \file subiso_triangles.cpp

    \brief This code uses the subgraph isomorphism algorithm to find
           triangles.

    \author Jon Berry (jberry@sandia.gov)

    \date 7/23/2010
*/
/****************************************************************************/

//#define DEBUG
//#define DOUBLE_DEBUG

#define NUM_SHEP 1

#include <set>

#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/subgraph_adapter.hpp>
#include <mtgl/subiso_triangles.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/triangles.hpp>
#include <mtgl/stats.hpp>
#include <mtgl/random.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace std;
using namespace mtgl;

//#define FILTERED
#define DIRECTION undirectedS
//#define TRIALS 30
//#define TRIALS 3
#define TRIALS 1

/***/

template <typename size_type>
struct result_triple {
  result_triple(size_type aa, size_type bb, size_type cc) :
    a(aa), b(bb), c(cc) {}

  size_type a, b, c;

  int operator<(const result_triple& other) const
  {
    if (a < other.a)
    {
      return true;
    }
    else if (a > other.a)
    {
      return false;
    }
    else if (b < other.b)
    {
      return true;
    }
    else if (b > other.b)
    {
      return false;
    }
    else if (c < other.c)
    {
      return true;
    }
    else
    {
      return (c < other.c);
    }
  }
};

static int cmp_uint64_t(const void* i1, const void* i2)
{
  return (int) ((*(uint64_t*) i1) - (*(uint64_t*) i2));
}

template <typename Graph>
class count_tricolor_triangles : public default_triangles_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  count_tricolor_triangles(Graph& gg, vertex_property_map<Graph, size_type>& vp,
                           set<result_triple<size_type> >& res,
                           size_type& ct, size_type& tcct) :
    g(gg), vpm(vp), count(ct), tricount(tcct), result(res), eiter(edges(g)) {}

  void operator()(size_type e1, size_type e2, size_type e3)
  {
    mt_incr(count, 1);

    edge_descriptor edg1 = eiter[e1];
    edge_descriptor edg2 = eiter[e2];
    edge_descriptor edg3 = eiter[e3];

    vertex_descriptor e1_src = source(edg1, g);
    vertex_descriptor e1_dest = target(edg1, g);
    vertex_descriptor e2_src = source(edg2, g);
    vertex_descriptor e2_dest = target(edg2, g);
    vertex_descriptor e3_src = source(edg3, g);
    vertex_descriptor e3_dest = target(edg3, g);

    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

    size_type v[6];
    v[0] = get(vid_map, e1_src);
    v[1] = get(vid_map, e1_dest);
    v[2] = get(vid_map, e2_src);
    v[3] = get(vid_map, e2_dest);
    v[4] = get(vid_map, e3_src);
    v[5] = get(vid_map, e3_dest);

    qsort(v, 6, sizeof(size_type), cmp_uint64_t);

    size_type eid[3];
    eid[0] = e1;
    eid[1] = e2;
    eid[2] = e3;

    qsort(eid, 3, sizeof(size_type), cmp_uint64_t);

    if ((vpm[e1_src] != vpm[e1_dest]) && (vpm[e2_src] != vpm[e2_dest]) &&
        (vpm[e3_src] != vpm[e3_dest]))
    {
      mt_incr(tricount, 1);

      // not yet threadsafe
      // size_type before_size = result.size();
      // result.insert(result_triple<size_type>(eid[0], eid[1], eid[2]));
      // size_type after_size = result.size();

#ifdef DEBUG
      printf("triangle: (%lu[%lu], %lu[%lu]) (%lu[%lu], %lu[%lu]) "
             "(%lu[%lu], %lu[%lu])\n",
             get(vid_map, e1_src), vpm[e1_src],
             get(vid_map, e1_dest), vpm[e1_dest],
             get(vid_map, e2_src), vpm[e2_src],
             get(vid_map, e2_dest), vpm[e2_dest],
             get(vid_map, e3_src), vpm[e3_src],
             get(vid_map, e3_dest), vpm[e3_dest]);
#endif
    }
  }

private:
  Graph& g;
  vertex_property_map<Graph, size_type>& vpm;
  size_type& count;
  size_type& tricount;
  set<result_triple<size_type> >& result;
  edge_iterator eiter;
};

template <typename Graph, typename Visitor>
void compare_to_purely_random_walks(
  Graph& g,
  vertex_property_map<Graph, typename graph_traits<Graph>::size_type>& color,
  double p,
  typename graph_traits<Graph>::size_type nc,
  int trial,
  double* purely_random_expected,
  double* purely_random_expected_semi_cheat,
  double* purely_random_expected_cheat,
  double exp_tricolor_tri,
  double computed_tricolor_tri,
  Visitor& vis)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  // nc *= 10;  // JWB DEBUG

  size_type n = num_vertices(g);
  size_type m = num_edges(g);

  double exp3paths4verts = p * p * p * n * (n - 1) * (n - 2) * (n - 3);
  double exp3paths3verts3edges = p * p * p * n * (n - 1) * (n - 2);
  double exp3paths3verts2edges = p * p * n * (n - 1) * (n - 2) * 2;
  double exp3paths2verts = p * n * (n - 1);

  double denom = exp3paths4verts + exp3paths3verts3edges +
                 exp3paths3verts2edges + exp3paths2verts;

  double exp_percentage_4verts = exp3paths4verts / denom;
  double exp_percentage_3verts3edges = exp3paths3verts3edges / denom;
  double exp_percentage_3verts2edges = exp3paths3verts2edges / denom;
  double exp_percentage_2verts = exp3paths2verts / denom;

  double pr_tricolor_triangle = (6 * exp_tricolor_tri) / denom;
  double exp_tricolor_hits = nc * pr_tricolor_triangle;
  double correction_factor = 1 - exp(-(exp_tricolor_hits / exp_tricolor_tri));

//  double correction_factor = 1 - pow(1 - 1 / exp_tricolor_tri,
//                                     exp_tricolor_hits);

#ifdef DEBUG
  printf("denom: %f\n", denom);
  printf("pr_tricolor_from_walk: %f\n", pr_tricolor_triangle);
  printf("exp_tricolor_hits (balls): %f\n", exp_tricolor_hits);
  printf("exp_tricolor_tri (bins): %f\n", exp_tricolor_tri);
  printf("correction_factor: %f\n", correction_factor);
  printf("exp_percentage_4verts: %f\n", exp_percentage_4verts);
  printf("exp_percentage_3verts3edges: %f\n", exp_percentage_3verts3edges);
  printf("exp_percentage_3verts2edges: %f\n", exp_percentage_3verts2edges);
  printf("exp_percentage_2verts: %f\n", exp_percentage_2verts);
#endif

  purely_random_expected[trial] = exp_tricolor_tri * correction_factor;

  printf("exp_full_bins: %f\n", purely_random_expected[trial]);

  vertex_iterator verts = vertices(g);
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> eid_map = get(_edge_id_map, g);

  size_type discovered = 0;
  size_type actual_tricolor_hits = 0;
  size_type actual_3paths_4vertices = 0;
  size_type actual_3paths_3vertices3edges = 0;
  size_type actual_3paths_3vertices2edges = 0;
  size_type actual_3paths_2vertices = 0;

  for (size_type i = 0; i < nc; ++i)
  {
    dynamic_array<vertex_descriptor> pverts;
    dynamic_array<edge_descriptor> pedges;

    size_type start = mt_lrand48() % n;

    vertex_descriptor src = verts[start];

    random_walk(g, src, 4, pverts, pedges);

    size_type v[4];
    v[0] = get(vid_map, pverts[0]);
    v[1] = get(vid_map, pverts[1]);
    v[2] = get(vid_map, pverts[2]);
    v[3] = get(vid_map, pverts[3]);

    if (v[0] != v[1] && v[1] != v[2] && v[0] != v[2] && v[3] == v[0] &&
        color[pverts[0]] != color[pverts[1]] &&
        color[pverts[1]] != color[pverts[2]] &&
        color[pverts[0]] != color[pverts[2]])
    {
      ++actual_tricolor_hits;

      dynamic_array<size_type> path(7);
      path[0] = v[0];
      path[1] = get(eid_map, pedges[0]);
      path[2] = v[1];
      path[3] = get(eid_map, pedges[1]);
      path[4] = v[2];
      path[5] = get(eid_map, pedges[2]);
      path[6] = v[3];

      vis(path, "pure");
    }

    // **************************************************************
    // ** Debug
    // **************************************************************
    xmt_hash_table<size_type, size_type> path_vertices(10);
    xmt_hash_table<size_type, size_type> path_edges(10);

    path_vertices.insert(v[0], 0);
    path_vertices.insert(v[1], 0);
    path_vertices.insert(v[2], 0);
    path_vertices.insert(v[3], 0);

    path_edges.insert(get(eid_map, pedges[0]), 0);
    path_edges.insert(get(eid_map, pedges[1]), 0);
    path_edges.insert(get(eid_map, pedges[2]), 0);

    if (path_vertices.size() == 4)
    {
      ++actual_3paths_4vertices;
    }
    else if (path_vertices.size() == 3 && path_edges.size() == 3)
    {
      ++actual_3paths_3vertices3edges;
    }
    else if (path_vertices.size() == 3 && path_edges.size() == 2)
    {
      ++actual_3paths_3vertices2edges;
    }
    else if (path_vertices.size() == 2)
    {
      ++actual_3paths_2vertices;
    }
    else
    {
      printf("OOOOPS! %lu %lu\n", path_vertices.size(), path_edges.size());
    }
    // **************************************************************
  }

  printf("actual_percentage_4verts: %f\n",
         actual_3paths_4vertices / (double) nc);
  printf("actual_percentage_3verts3edges: %f\n",
         actual_3paths_3vertices3edges / (double) nc);
  printf("actual_percentage_3verts2edges: %f\n",
         actual_3paths_3vertices2edges / (double) nc);
  printf("actual_percentage_2verts: %f\n",
         actual_3paths_2vertices / (double) nc);
  printf("actual_tricolor_hits(balls): %lu\n", actual_tricolor_hits);

  purely_random_expected_semi_cheat[trial] = computed_tricolor_tri *
                                             correction_factor;

  double correction_factor_cheat =
    1 - exp(-(actual_tricolor_hits / computed_tricolor_tri));

  printf("actual correction factor: %f\n", correction_factor_cheat);
  printf("actual_tricolor_tri (bins): %f\n", computed_tricolor_tri);
  printf("exp_full_bins based on actual balls and bins: %f\n",
         computed_tricolor_tri * correction_factor_cheat);

  purely_random_expected_cheat[trial] = computed_tricolor_tri *
                                        correction_factor_cheat;
}

// tri_match_visitor will be called once per match. subiso_triangles.hpp
// has its own routine ("si_bipartite_visitor") to determine which
// candidate matches will trigger calls to si_match_visitor.
template <typename Graph>
class tri_user_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  tri_user_visitor(Graph& gg,
                   set<result_triple<size_type> >& s) :
    g(gg), result_set(s), order(num_vertices(g)),
    vid_map(get(_vertex_id_map, g)), eid_map(get(_edge_id_map, g)){}

  void operator()(dynamic_array<vertex_descriptor>& path_verts,
                  dynamic_array<edge_descriptor>& path_edges,
                  const char* str = 0)
  {
    size_type vid[3];
    vid[0] = get(vid_map, path_verts[0]);
    vid[1] = get(vid_map, path_verts[1]);
    vid[2] = get(vid_map, path_verts[2]);

    qsort(vid, 3, sizeof(size_type), cmp_uint64_t);

    size_type eid[3];
    eid[0] = get(eid_map, path_edges[0]);
    eid[1] = get(eid_map, path_edges[1]);
    eid[2] = get(eid_map, path_edges[2]);

    qsort(eid, 3, sizeof(size_type), cmp_uint64_t);

    result_triple<size_type> t(eid[0], eid[1], eid[2]);
    result_set.insert(t);

#ifdef DEBUG
    if (str == 0)
    {
      printf("subiso: v( %lu %lu %lu )  e( %lu %lu %lu )\n",
             vid[0], vid[1], vid[2], eid[0], eid[1], eid[2]);
    }
    else
    {
      printf("%s subiso: v( %lu %lu %lu )  e( %lu %lu %lu )\n", str,
             vid[0], vid[1], vid[2], eid[0], eid[1], eid[2]);
    }
#endif
  }

private:
  Graph& g;
  set<result_triple<size_type> >& result_set;
  size_type order;
  const vertex_id_map<Graph>& vid_map;
  const edge_id_map<Graph>& eid_map;
};

#ifdef _WIN32
int subiso_main(int argc, char* argv[])
#else
int main(int argc, char* argv[])
#endif
{
  // typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;

  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;
  typedef graph_traits<Graph>::edge_descriptor edge;
  typedef graph_traits<Graph>::size_type size_type;
  typedef map_property_map<size_type, vertex_id_map<Graph> >
          vertex_property_map;
  typedef map_property_map<size_type, edge_id_map<Graph> > edge_property_map;

  typedef subgraph_adapter<Graph> SubGraph;
  typedef graph_traits<SubGraph>::vertex_descriptor vertex_descriptor_sg;
  typedef graph_traits<SubGraph>::edge_descriptor edge_sg;
  typedef graph_traits<SubGraph>::size_type size_type_sg;
  typedef array_property_map<int, vertex_id_map<SubGraph> >
          vertex_property_map_sg;
  typedef array_property_map<int, edge_id_map<SubGraph> > edge_property_map_sg;

  mt_srand48(0);

  init_test(argc, argv);

  Graph g;
  create_test_graph(g, argc, argv);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  double p = strcmp(argv[1], "er") == 0 ?
             atof(argv[3]) : size / (double) (order * (order - 1) / 2);

  printf("n: %lu,  m: %lu\n", order, size);

  const size_type numVerts = order;
  const size_type numEdges = size;

  // Compute expected number of triangles if ER random graph.
  double Etri = 0;
  double Etc_tri = 0;

  xmt_hash_table<size_type, size_type> vTypes(2 * numVerts);
  xmt_hash_table<size_type, size_type> eTypes(2 * numEdges);

  // Make type maps.
  vertex_id_map<Graph> g_vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> g_eid_map = get(_edge_id_map, g);

  mtgl::vertex_property_map<Graph, size_type> vpm(g);
  mtgl::edge_property_map<Graph, size_type> epm(g);

  vertex_iterator verts = vertices(g);
  edge_iterator edgs = edges(g);

  // Set vertex types (a, b, c, d, e).
  #pragma mta assert parallel
  for (size_type i = 0; i < order; ++i) vpm[verts[i]] = mt_lrand48() % 3;

#ifdef DEBUG
  for (size_type i = 0; i < order; ++i)
  {
    printf("vpm[%lu]: %lu\n", i, vpm[verts[i]]);
  }
#endif

  size_type num_of_color[3];
  num_of_color[0] = 0;
  num_of_color[1] = 0;
  num_of_color[2] = 0;

  for (size_type i = 0; i < order; ++i) ++num_of_color[vpm[verts[i]]];

  printf("observed number of color 0: %lu\n", num_of_color[0]);
  printf("observed number of color 1: %lu\n", num_of_color[1]);
  printf("observed number of color 2: %lu\n", num_of_color[2]);

  #pragma mta assert nodep
  for (size_type i = 0; i < size; ++i) epm[edgs[i]] = 0;

  size_type actual_tri = 0;
  size_type actual_tricolor_tri = 0;

  set<result_triple<size_type> > result;

  count_tricolor_triangles<Graph> ctv(g, vpm, result,
                                      actual_tri, actual_tricolor_tri);

  mt_timer ft_time;
  ft_time.start();

  find_triangles(g, ctv);

  ft_time.stop();
  printf("actual_tricolor_triangle count: %lu (%f seconds)\n",
         actual_tricolor_tri, ft_time.getElapsedSeconds());

  if (strcmp(argv[1], "er") == 0)
  {
    double n_choose_3 = (order * (order - 1) * (order - 2)) / 6.0;
    Etri = n_choose_3 * p * p * p;
    Etc_tri = (Etri / 27.0) * 6;

    printf("p: %f\n", p);
    printf("n_choose_3: %f\n", n_choose_3);
    printf("E[#triangles]: %f\n", Etri);
    printf("E[#tri-color triangles]: %f\n", Etc_tri);
  }

  size_type results[TRIALS];
  size_type nc_sizes[TRIALS];

#ifdef DOUBLE_DEBUG
  double purely_random_expected[TRIALS];
  double purely_random_expected_cheat[TRIALS];
  double purely_random_expected_semi_cheat[TRIALS];
  size_type purely_random_discovered[TRIALS];
#endif

  for (int i = 0; i < TRIALS; ++i)
  {
    set<result_triple<size_type> > result;

    tri_user_visitor<Graph> simv(g, result);

    mt_timer timer;
    timer.start();

    size_type nc = subiso_triangles(g, vpm, epm, simv);

    timer.stop();

    printf("time for %d'th trial: %f (nc: %lu)\n", i, timer.getElapsedSeconds(),
           nc);
    fflush(stdout);
    printf("result.size: %lu\n", result.size());

    results[i] = result.size();
    nc_sizes[i] = nc;

#ifdef DOUBLE_DEBUG
    xmt_hash_table<uint64_t, uint64_t> xt2(1000000);

    tri_user_visitor<Graph> purv(g, xt2);

    compare_to_purely_random_walks(g, vpm, p, nc, i,
                                   purely_random_expected,
                                   purely_random_expected_semi_cheat,
                                   purely_random_expected_cheat,
                                   Etc_tri, actual_tricolor_tri, purv);

    purely_random_discovered[i] = xt2.size();

    printf("expected: %f\n", purely_random_expected[i]);
    fflush(stdout);
    printf("expect_semi: %f\n", purely_random_expected_semi_cheat[i]);
    fflush(stdout);
    printf("discovered: %lu\n", purely_random_discovered[i]);
    fflush(stdout);
#endif
  }

#ifdef DOUBLE_DEBUG
  double e_rand_exp_tri = E(purely_random_expected, TRIALS);
  double e_rand_exp_tri_cheat = E(purely_random_expected_cheat, TRIALS);
  double e_rand_exp_tri_semi_cheat = E(purely_random_expected_semi_cheat,
                                       TRIALS);
  double e_rand_disc_tri = E(purely_random_discovered, TRIALS);
#endif

  printf("\n");

  double e_nc_size = E(nc_sizes, TRIALS);
  double e_tri = E(results, TRIALS);
  double sigma = std_dev(results, TRIALS);

#ifdef DOUBLE_DEBUG
  printf("%lu %lu %f %f %f %lu %lu %f %f %f %f %f %f\n", order, size, e_nc_size,
         Etri, Etc_tri, actual_tri, actual_tricolor_tri, e_tri, sigma,
         e_rand_exp_tri, e_rand_exp_tri_semi_cheat,
         e_rand_exp_tri_cheat, e_rand_disc_tri);
#else
  printf("%lu %lu %f %f %f %lu %lu %f %f\n", order, size, e_nc_size,
         Etri, Etc_tri, actual_tri, actual_tricolor_tri, e_tri, sigma);
#endif

  return 0;
}
