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
/*! \file metrics.hpp

    \brief This code contains the classes used by the graph metrics such as
           assortativity, cluster coefficient, degree degree correlation,
           degree distribution, and modularity.

    \author Jon Berry (jberry@sandia.gov)
    \author Vitus Leung (vjleung@sandia.gov)

    \date 11/2006
*/
/****************************************************************************/

#ifndef MTGL_METRICS_HPP
#define MTGL_METRICS_HPP

#include <cstdio>
#include <climits>
#include <cassert>
#include <cmath>

#include <mtgl/util.hpp>
#include <mtgl/xmt_hash_set.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/psearch.hpp>
#include <mtgl/triangles.hpp>
#include <mtgl/dynamic_array.hpp>

namespace mtgl {

// We'll store the neighbor id's of high-degree vertices in hash sets.
// The following structure will serve as an input to the xmt_hash_set::visit()
// method, which will be called when two high-degree vertices are adjacent.
struct hvis {
  hvis(int sid, int tid, xmt_hash_set<int>* tnghs,
       dynamic_array<triple<int, int, int> >& res) :
    srcid(sid), destid(tid), tneighs(tnghs), result(res) {}

  void operator()(int neigh)
  {
    if (tneighs->member(neigh))
    {
      triple<int, int, int> t(srcid, destid, neigh);
      result.push_back(t);
    }
  }

  int srcid, destid;
  xmt_hash_set<int>* tneighs;
  dynamic_array<triple<int, int, int> >& result;
};

class htvis {
public:
  htvis(int types) : same_type_sum(0), ass_types(types) {}

  void operator()(const int& k, int& v)
  {
    if (k / ass_types == k % ass_types)
    {
      same_type_sum += v;
    }
  }

  int same_type_sum;

private:
  int ass_types;
};

class intvis {
public:
  intvis() : all_types_sum(0){}

  void operator()(int data) { all_types_sum += data; }

  int all_types_sum;
};

/// \brief Describes what happens during a graph search to find the type
///        mixing matrix used to compute assortativity and modularity.
template <typename Graph>
class type_mixing_matrix_visitor : public default_psearch_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  type_mixing_matrix_visitor(int* res, int* typs, int* map, Graph& _g,
                             int* subp = 0, int my_subp = 0) :
    result(res), types(typs), subproblem(subp), id2comm(map),
    my_subproblem(my_subp), vid_map(get(_vertex_id_map, _g)), g(_g) {}

  bool visit_test(edge_descriptor e, vertex_descriptor v)
  {
    return (!subproblem || subproblem[get(vid_map, v)] == my_subproblem);
  }

  void tree_edge(edge_descriptor e, vertex_descriptor src)
  {
    #pragma mta trace "te"
    int st = id2comm[types[get(vid_map, src)]];
    int dt = id2comm[types[get(vid_map, target(e, g))]];

    if (st == dt)
    {
      mt_incr(result[0], 1);
    }
    else
    {
      mt_incr(result[1], 1);
    }
  }

  void back_edge(edge_descriptor e, vertex_descriptor src)
  { tree_edge(e, src); }

protected:
  int* result;
  int* types;
  int* subproblem;
  int my_subproblem;
  int* id2comm;         // Given id of leader node, get commun. #
  vertex_id_map<Graph> vid_map;
  Graph& g;
};

/// \brief Describes what happens during a graph search to find the
///        degree degree correlation.
template <typename Graph>
class degree_degree_correlation_visitor : public default_psearch_visitor<Graph>
{
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  degree_degree_correlation_visitor(int** res, Graph& _g,
                                    int* subp = 0, int my_subp = 0) :
    result(res), subproblem(subp), my_subproblem(my_subp),
    vid_map(get(_vertex_id_map, g)), g(_g) {}

  bool visit_test(edge_descriptor e, vertex_descriptor v)
  {
    return (!subproblem || subproblem[get(vid_map, v)] == my_subproblem);
  }

  void tree_edge(edge_descriptor e, vertex_descriptor src)
  {
    if (source(e, g) == src)
    {
      #pragma mta trace "te"
      vertex_descriptor dest = target(e, g);

      int log2d = (int) ceil(log(out_degree(dest, g) + 1.0) / log(2.0)) + 1;
      int log2s = (int) ceil(log(out_degree(src, g) + 1.0) / log(2.0)) + 1;

      mt_incr(result[log2d][log2s], 1);
    }
  }

  void back_edge(edge_descriptor e, vertex_descriptor src)
  { tree_edge(e, src); }

protected:
  int** result;
  int* subproblem;
  int my_subproblem;
  vertex_id_map<Graph> vid_map;
  Graph& g;
};

/// \brief Describes what happens during a graph search to find the
///        degree distribution.
template <typename Graph>
class degree_distribution_visitor : public default_psearch_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  degree_distribution_visitor(Graph& gg, int* res,
                              int* subp = 0, int my_subp = 0) :
    g(gg), result(res), subproblem(subp), my_subproblem(my_subp),
    vid_map(get(_vertex_id_map, g)) {}

  void pre_visit(vertex_descriptor v)
  {
    #pragma mta trace "te"
    int log2d = (int) ceil(log(out_degree(v, g) + 1.0) / log(2.0)) + 1;
    mt_incr(result[log2d], 1);
  }

  bool visit_test(edge_descriptor e, vertex_descriptor v)
  {
    return (!subproblem || subproblem[get(vid_map, v)] == my_subproblem);
  }

protected:
  Graph& g;
  int* result;
  int* subproblem;
  int my_subproblem;
  vertex_id_map<Graph> vid_map;
};

/// \brief Describes what happens during a graph search to find the
///        connected triples used to find the cluster coefficient.
template <typename Graph>
class connected_triples_visitor : public default_psearch_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  connected_triples_visitor(Graph& gg, int* c, int* subp = 0, int my_subp = 0) :
    g(gg), count(c), subproblem(subp), my_subproblem(my_subp),
    vid_map(get(_vertex_id_map, g)) {}

  void pre_visit(vertex_descriptor v)
  {
    #pragma mta trace "te"
    int sd = out_degree(v, g);
    if (sd > 1)
    {
      int act = sd * (sd - 1) / 2;
      mt_incr(*count, act);
    }
  }

  bool visit_test(edge_descriptor e, vertex_descriptor v)
  {
    return (!subproblem || subproblem[get(vid_map, v)] == my_subproblem);
  }

protected:
  Graph& g;
  int* count;
  int* subproblem;
  int my_subproblem;
  vertex_id_map<Graph> vid_map;
};

/// \brief Given a set of community assignments (on nodes) and an edge
///        support, compute the "support variance" -- a measure of how well
///        the community assigment conforms to the support.
template <typename Graph>
class support_variance_visitor {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  support_variance_visitor(Graph& gg, int* comm, double* supp, double& res,
                           vertex_id_map<Graph>& vm, edge_id_map<Graph>& em) :
    g(gg), result(res), community(comm), support(supp),
    vid_map(vm), eid_map(em) {}

  void operator()(edge_descriptor e)
  {
    int eid = get(eid_map, e);

    vertex_descriptor src = source(e, g);
    vertex_descriptor trg = target(e, g);
    int sid = get(vid_map, src);
    int tid = get(vid_map, trg);

    int same = (community[sid] == community[tid]);
    double diff = same - support[eid];
    double res = mt_readfe(result);

    mt_write(result, res + diff * diff);
  }

protected:
  Graph& g;
  double& result;
  int* community;
  double* support;
  vertex_id_map<Graph>& vid_map;
  edge_id_map  <Graph>& eid_map;
};

/// \brief This is the algorithm that users will invoke to find the type
///        mixing matrix of a graph.  Used to find assortativity and
///        modularity.
template <typename Graph>
void type_mixing_matrix(Graph& g, int* result, int* types, int* id2comm,
                        int* subproblem = 0, int my_subproblem = 0)
{
  type_mixing_matrix_visitor<Graph> av(result, types, id2comm, g,
                                       subproblem, my_subproblem);
  psearch_high_low<Graph, type_mixing_matrix_visitor<Graph>,
                   AND_FILTER, UNDIRECTED, 10> psrch(g, av);
  psrch.run();
}

/// \brief This is the algorithm that users will invoke to find the
///        assortativity of a graph.
template <typename Graph>
double assortativity(Graph& g, int num_types)
{
  int ord = num_vertices(g);
  int* server = new int[ord];
  int* id_inversemap = new int[ord];

  int tmm[2] = { 0, 0 };

  #pragma mta assert nodep
  for (int i = 0; i < ord; ++i)
  {
    server[i] = i;
    id_inversemap[i] = i % num_types;
  }

  //mta_resume_event_logging();
  type_mixing_matrix(g, tmm, server, id_inversemap);
  //mta_suspend_event_logging();

  int same = tmm[0];
  int all = tmm[0] + tmm[1];

  delete [] server;
  delete [] id_inversemap;

  return (double) same / all;
}

/// \brief This is the algorithm object that users will invoke to find the
///        degree degree correlation of a graph.
template <typename Graph>
void degree_degree_correlation(Graph& g, double** result,
                               int* subproblem = 0, int my_subproblem = 0)
{
  unsigned int ord = num_vertices(g);
  int size = num_edges(g);
  int maxdegree = (int) ceil(log(2.0 * ord - 1) / log(2.0)) + 1;

  int** res = (int**) malloc(maxdegree * sizeof(int*));
  #pragma mta assert parallel
  for (int i = 0; i < maxdegree; ++i)
  {
    res[i] = (int*) calloc(maxdegree, sizeof(int));
  }

  degree_degree_correlation_visitor<Graph> ddcv(res, g, subproblem,
                                                my_subproblem);
  psearch_high_low<Graph, degree_degree_correlation_visitor<Graph>,
                   AND_FILTER, UNDIRECTED, 10> psrch(g, ddcv);
  psrch.run();

  #pragma mta assert nodep
  for (int i = 1; i < maxdegree; ++i)
  {
    for (int j = 1; j < maxdegree; ++j)
    {
      result[i][j] = (double) res[i][j] / size;
    }
  }

  for (int i = 0; i < maxdegree; ++i) free(res[i]);
  free(res);
}

/// \brief This is the algorithm object that users will invoke to find the
///        degree distribution of a graph.
template <typename Graph>
void
find_degree_distribution(Graph& g, double* result,
                         int* subproblem = 0, int my_subproblem = 0)
{
  unsigned int ord = num_vertices(g);
  unsigned int m = (int) ceil(log(2.0 * ord - 1) / log(2.0)) + 1;

  int* res = (int*) calloc(m, sizeof(int));

  degree_distribution_visitor<Graph> ddv(g, res, subproblem, my_subproblem);
  psearch_high_low<Graph, degree_distribution_visitor<Graph>,
                   AND_FILTER, UNDIRECTED, 10> psrch(g, ddv);
  psrch.run();

#ifdef DEBUG
  printf("degree distribution totals:\n");
  for (int i = 0; i < m; ++i) printf("%d %d\n", i, res[i]);
#endif

  #pragma mta assert nodep
  for (int i = 0; i < m; ++i) result[i] = (double) res[i] / ord;

  free(res);
}

/// \brief This is the algorithm object that users will invoke to find the
///        connected triples of a graph.  Used to find the cluster coefficient
///        of a graph.
template <typename Graph>
int find_connected_triples(Graph& g, int* subproblem = 0, int my_subproblem = 0)
{
  int count2 = 0;

  connected_triples_visitor<Graph> ctv(g, &count2, subproblem, my_subproblem);
  psearch_high_low<Graph, connected_triples_visitor<Graph>,
                   AND_FILTER, UNDIRECTED, 10> psrch(g, ctv);
  psrch.run();

  return count2;
}

template <typename Graph>
class count_triangles : public default_triangles_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;

  count_triangles(Graph& gg,
                  dynamic_array<triple<size_type, size_type, size_type> >& res,
                  size_type& ct) :
    g(gg), count(ct), result(res) {}

  void operator()(size_type e1, size_type e2, size_type e3)
  { mt_incr(count, 1); }

private:
  Graph& g;
  size_type& count;
  dynamic_array<triple<size_type, size_type, size_type> >& result;
};

/// \brief This is the algorithm object that users will invoke to find the
///        cluster coefficient of a graph.
template <typename Graph>
double cluster_coefficient(Graph& g)
{
  typedef typename count_triangles<Graph>::size_type ct_size_type;

  ct_size_type count = 0;

  dynamic_array<triple<ct_size_type, ct_size_type, ct_size_type> > result;

  count_triangles<Graph> ctv(g, result, count);

  find_triangles(g, ctv);

#ifdef DEBUG
  fprintf(stdout, "found %lu triangles\n", (long unsigned) count);
#endif

  int count2 = find_connected_triples(g);

#ifdef DEBUG
  fprintf(stdout, "found %d connected triples\n", count2);
#endif

  return (double) 3 * count / count2;
}

template <typename Graph>
class edge_type_matcher {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  edge_type_matcher(Graph& gg, int* clust, bool* same_clust,
                    vertex_id_map<Graph>& vmap,
                    edge_id_map<Graph>& emap, int& mism) :
    g(gg), cluster(clust), same_cluster(same_clust),
    vid_map(vmap), eid_map(emap), mismatches(mism) {}

  void operator()(edge_descriptor e)
  {
    int eid = get(eid_map, e);

    // in this context, 0 -> same cluster
    //                  1 -> diff cluster
    vertex_descriptor v1 = source(e, g);
    vertex_descriptor v2 = target(e, g);
    int v1id = get(vid_map, v1);
    int v2id = get(vid_map, v2);

    bool test_same_cluster = (cluster[v1id] == cluster[v2id]);

    if (test_same_cluster != !same_cluster[eid]) mt_incr(mismatches, 1);
  }

private:
  Graph& g;
  vertex_id_map<Graph>& vid_map;
  edge_id_map<Graph>& eid_map;
  int* cluster;
  bool* same_cluster;
  int& mismatches;
};

template <typename Graph>
class edge_hamming_distance {
public:
  edge_hamming_distance(Graph& gg, int* clust, bool* same_clust) :
    g(gg), cluster(clust), same_cluster(same_clust) {}

  int run()
  {
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> eid_map = get(_edge_id_map, g);

    int mismatches = 0;
    edge_type_matcher<Graph> etm(g, cluster, same_cluster,
                                 vid_map, eid_map, mismatches);
    visit_edges(g, etm);

    return mismatches;
  }

private:
  Graph& g;
  int* cluster;
  bool* same_cluster;
};

template <typename Graph>
class support_variance {
public:
  support_variance(Graph& gg, int* comm, double* supp) :
    g(gg), community(comm), support(supp) {}

  double run()
  {
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> eid_map =  get(_edge_id_map, g);

    double sv = 0.0;

    support_variance_visitor<Graph> svv(g, community, support,
                                        sv, vid_map, eid_map);
    visit_edges(g, svv);

    return sv;
  }

private:
  Graph& g;
  int* community;
  double* support;
};

template <typename Graph>
int* degree_distribution(Graph& g, int num_bins = 32)
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;

  unsigned int order = num_vertices(g);

  int* count  = new int[num_bins];
  int* degs  = new int[order];

  #pragma mta assert nodep
  for (int i = 0; i < num_bins; ++i) count[i] = 0;

  vertex_iterator verts = vertices(g);

  #pragma mta assert nodep
  for (unsigned int i = 0; i < order; ++i)
  {
    vertex_descriptor v = verts[i];
    degs[i] = out_degree(v, g);
  }

  for (unsigned i = 0; i < order; ++i)
  {
    int log2deg = (int) (log((double) degs[i]) / log(2.0));
    if (log2deg < 0) log2deg = 0;

    int index = (log2deg >= num_bins) ? num_bins - 1 : log2deg;
    ++count[index];
  }

  delete [] degs;

  return count;
}

}

#endif
