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
/*! \file subiso_5cycles.cpp

    \brief This code uses the subgraph isomorphism algorithm to find
           chordless 5-cycles.

    \author Jon Berry (jberry@sandia.gov)

    \date 7/23/2010
*/
/****************************************************************************/

//#define DEBUG

#include <mtgl/adjacency_list.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/subgraph_adapter.hpp>
#include <mtgl/write_graphlets.hpp>
#include <mtgl/subiso_5cycles.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/stats.hpp>
#include <mtgl/mtgl_test.hpp>

using namespace mtgl;

//#define FILTERED
#define DIRECTION undirectedS
//#define TRIALS 30
//#define TRIALS 3
#define TRIALS 1

template <typename Graph, typename Visitor>
void compare_to_purely_random_walks(Graph& g,
                                    vertex_property_map<Graph, uint8_t>& color,
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

static int cmp_uint64_t(const void* i1, const void* i2)
{
  return (int) ((*(uint64_t*) i1) - (*(uint64_t*) i2));
}

// cycle_user_visitor:  useful only for small graphs.
template <typename Graph>
class cycle_user_visitor {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  cycle_user_visitor(Graph& gg, xmt_hash_table<uint64_t, uint64_t>& ht) :
    g(gg), t_table(ht), order(num_vertices(g)),
    vid_map(get(_vertex_id_map, g)) {}

  void operator()(dynamic_array<vertex_descriptor>& path_verts,
                  dynamic_array<edge_descriptor>& path_edges, char* str = 0)
  {
    size_type v[5];
    v[0] = get(vid_map, path_verts[0]);
    v[1] = get(vid_map, path_verts[1]);
    v[2] = get(vid_map, path_verts[2]);
    v[3] = get(vid_map, path_verts[3]);
    v[4] = get(vid_map, path_verts[4]);

    qsort(v, 5, sizeof(size_type), cmp_uint64_t);

    // Won't work with big graphs!  Need to hash the string.
    uint64_t key = v[0] * order * order * order * order +
                   v[1] * order * order * order +
                   v[2] * order * order +
                   v[3] * order + v[4];

    if ((t_table.insert(key, key)).second == false) return;

#ifdef DEBUG
    if (str == 0)
    {
      printf("sub: 5cycle: %lu %lu %lu %lu %lu (key: %lu) \n",
             v[0], v[1], v[2], v[3], v[4], key);
    }
    else
    {
      printf("%s sub: 5cycle: %lu %lu %lu %lu %lu (key: %lu)\n",
             str, v[0], v[1], v[2], v[3], v[4], key);
    }
#endif

//    printf(".");
//    fflush(stdout);
  }

  void operator()(size_type* path, size_type path_length, char* str = 0)
  {
    size_type v[5];
    v[0] = path[0];
    v[1] = path[2];
    v[2] = path[4];
    v[3] = path[6];
    v[4] = path[8];

    qsort(v, 5, sizeof(size_type), cmp_uint64_t);

    uint64_t key = v[0] * order * order * order * order +
                   v[1] * order * order * order +
                   v[2] * order * order +
                   v[3] * order + v[4];

    if ((t_table.insert(key, key)).second == false) return;

#ifdef DEBUG
    if (str == 0)
    {
      printf("sub: 5cycle: %lu %lu %lu %lu %lu (key: %lu) \n",
             v[0], v[1], v[2], v[3], v[4], key);
    }
    else
    {
      printf("%s sub: 5cycle: %lu %lu %lu %lu %lu (key: %lu)\n",
             str, v[0], v[1], v[2], v[3], v[4], key);
    }
#endif

//    printf(".");
//    fflush(stdout);
  }

  void operator()(dynamic_array<size_type>& path, char* str = 0)
  {
    size_type v[5];
    v[0] = path[0];
    v[1] = path[2];
    v[2] = path[4];
    v[3] = path[6];
    v[4] = path[8];

    qsort(v, 5, sizeof(size_type), cmp_uint64_t);

    uint64_t key = v[0] * order * order * order * order +
                   v[1] * order * order * order +
                   v[2] * order * order +
                   v[3] * order + v[4];

    if ((t_table.insert(key, key)).second == false) return;

#ifdef DEBUG
    if (str == 0)
    {
      printf("sub: 5cycle: %lu %lu %lu %lu %lu (key: %lu) \n",
             v[0], v[1], v[2], v[3], v[4], key);
    }
    else
    {
      printf("%s sub: 5cycle: %lu %lu %lu %lu %lu (key: %lu)\n",
             str, v[0], v[1], v[2], v[3], v[4], key);
    }
#endif

//    printf(".");
//    fflush(stdout);
  }

private:
  Graph& g;
  xmt_hash_table<uint64_t, uint64_t>& t_table;
  size_type order;
  vertex_id_map<Graph> vid_map;
};

double n_choose_5(double n)
{
  return n * (n - 1) / 120 * (n - 2) * (n - 3) * (n - 4);
}

int main(int argc, char* argv[])
{
  // typedef adjacency_list<DIRECTION> Graph;
  typedef compressed_sparse_row_graph<DIRECTION> Graph;

  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::edge_iterator edge_iterator;
  typedef graph_traits<Graph>::edge_descriptor edge;
  typedef graph_traits<Graph>::size_type size_type;
  typedef xmt_hash_table<size_type, size_type> hash_table_t;
  typedef mtgl::vertex_property_map<Graph, uint8_t> vertex_property_map_t;
  typedef mtgl::edge_property_map<Graph, uint8_t> edge_property_map_t;

  typedef subgraph_adapter<Graph> SubGraph;
  typedef graph_traits<SubGraph>::vertex_descriptor vertex_descriptor_sg;
  typedef graph_traits<SubGraph>::edge_descriptor edge_sg;
  typedef graph_traits<SubGraph>::size_type size_type_sg;

  mt_srand48(0);

  init_test(argc, argv);

  Graph g;
  create_test_graph(g, argc, argv);

#if 0
  write_graphlets(g, "/zap/jberry/07-06-graphlets/test.mtx");
#endif

  size_type order = num_vertices(g);
  size_type size = num_edges(g);
  printf("n: %lu,  m: %lu\n", order, size);

  const size_type numVerts = order;
  const size_type numEdges = size;

  // compute expected number of triangles if ER random graph
  double Ech5c = 0;
  double Efc_ch5c = 0;
  double p = strcmp(argv[1], "er") == 0 ?
             atof(argv[3]) : size / (double) (order * (order - 1) / 2);

  if (strcmp(argv[1], "er") == 0)
  {
    double nc5 = n_choose_5((double) order);
    double five_fac = 5 * 4 * 3 * 2 * 1;
    double p_fifth = p * p * p * p * p;
    double oneminusp_fifth = (1 - p) * (1 - p) * (1 - p) * (1 - p) * (1 - p);

    // Pick 5 vertices; order them; each and its reverse 5-count
    // a set of 5 edges that could be a 5-cycle.  Thus, we divide
    // the number of permutations by 10 to obtain the number of
    // 5-edge sets to test for chordless 5-cycle-ness.
    // Force these edges to be there; force the other 5 to
    // be missing.
    Ech5c = nc5 * five_fac / 10 * p_fifth * oneminusp_fifth;
    Efc_ch5c = Ech5c / (5 * 5 * 5 * 5 * 5) * 10;

    printf("p: %f\n", p);
    printf("n_choose_5: %f\n", nc5);
    printf("p^5: %f\n", p_fifth);
    printf("(1-p)^5: %f\n", oneminusp_fifth);
    printf("E[#chordless 5cycles]: %f\n", Ech5c);
    printf("E[#5-color chordless 5cycles]: %f\n", Efc_ch5c);
  }

  xmt_hash_table<size_type, size_type> vTypes(2 * numVerts);
  xmt_hash_table<size_type, size_type> eTypes(2 * numEdges);

  // Make type maps.
  vertex_id_map<Graph> g_vid_map = get(_vertex_id_map, g);
  edge_id_map<Graph> g_eid_map = get(_edge_id_map, g);

  vertex_property_map_t vpm(g);
  edge_property_map_t epm(g);

  vertex_iterator verts = vertices(g);
  edge_iterator edgs = edges(g);

  // Set vertex types (a,b,c,d,e)
  #pragma mta assert parallel
#ifdef DEBUG
  for (size_type i = 0; i < order; ++i) vpm[verts[i]] = i % 5;
#else
  for (size_type i = 0; i < order; ++i) vpm[verts[i]] = lrand48() % 5;
#endif

  size_type num_of_color[5];
  num_of_color[0] = 0;
  num_of_color[1] = 0;
  num_of_color[2] = 0;
  num_of_color[3] = 0;
  num_of_color[4] = 0;

  for (size_type i = 0; i < order; ++i) ++num_of_color[vpm[verts[i]]];

  printf("observed number of color 0 : %lu\n", num_of_color[0]);
  printf("observed number of color 1 : %lu\n", num_of_color[1]);
  printf("observed number of color 2 : %lu\n", num_of_color[2]);
  printf("observed number of color 3 : %lu\n", num_of_color[3]);
  printf("observed number of color 4 : %lu\n", num_of_color[4]);

#ifdef DEBUG
  print(g);
#endif

  #pragma mta assert nodep
  for (size_type i = 0; i < size; ++i) epm[edgs[i]] = 0;

  size_type results[TRIALS];
  size_type nc_sizes[TRIALS];
//  double purely_random_expected[TRIALS];
//  double purely_random_expected_cheat[TRIALS];
//  double purely_random_expected_semi_cheat[TRIALS];
//  size_type purely_random_discovered[TRIALS];

  for (int i = 0; i < TRIALS; ++i)
  {
    xmt_hash_table<uint64_t, uint64_t> xt(1000000);

    cycle_user_visitor<Graph> mv(g, xt);

    mt_timer timer;
    timer.start();

    size_type nc = subiso_5cycles(g, vpm, epm, mv);

    timer.stop();
    printf("time for %d'th trial: %f\n", i, timer.getElapsedSeconds());
    fflush(stdout);

    results[i] = xt.size();
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

//  double e_rand_exp_tri = E(purely_random_expected, TRIALS);
//  double e_rand_exp_tri_cheat = E(purely_random_expected_cheat, TRIALS);
//  double e_rand_exp_tri_semi_cheat = E(purely_random_expected_semi_cheat,
//                                       TRIALS);
//  double e_rand_disc_tri = E(purely_random_discovered, TRIALS);

  printf("\n");

  double e_nc_size = E(nc_sizes, TRIALS);
  double e_tri = E(results, TRIALS);
  double sigma = std_dev(results, TRIALS);

#ifdef DOUBLE_DEBUG
  printf("%lu %lu %f %f %f %lu %lu %f %f %f %f %f %f\n", order, size, e_nc_size,
         Ech5c, Efc_ch5c,
         actual_tri, actual_tricolor_tri, e_tri, sigma, e_rand_exp_tri,
         e_rand_exp_tri_semi_cheat, e_rand_exp_tri_cheat, e_rand_disc_tri);
#else
  printf("%lu %lu %f %f %f %f %f\n", order, size, e_nc_size, Ech5c, Efc_ch5c,
         e_tri, sigma);
#endif

  return 0;
}
